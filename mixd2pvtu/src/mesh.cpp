#include "mesh.hpp"
#include "utils.hpp"

#include <string>
#include <fstream>
#include <sstream>

#include <algorithm>
#include <experimental/filesystem>

#include <vtk/vtkCellArray.h>
#include <vtk/vtkDoubleArray.h>
#include <vtk/vtkPointData.h>
#include <vtk/vtkUnstructuredGrid.h>
#include <vtk/vtkXMLUnstructuredGridWriter.h>
#include <vtk/vtkXMLPUnstructuredGridWriter.h>
#include <vtk/vtkSmartPointer.h>

/*
* @brief: Mesh class
* @details: Reads mesh files and data, distributes them, and writes out.
* @param: [in] Config
*/
Mesh::Mesh(Config config)
{
    ndf        = config.getNdf();
    nrec       = config.getNrec();
    nrecoffset = config.getNrecoffset();
    nrecstride = config.getNrecstride();
    spacetime  = config.getSpacetime();
    dataFiles  = config.getDataFiles();
    title      = config.getTitle();
    outpath    = config.getOutpath();

    readminf(config.getMinfFile(), config.getSpacetime());
    distribute();

    elemType = config.getElemType();

    if (elemType == VTK_TETRA)
    {
        nen = 4;
        elem = new Tetrahedron[nec];
    }
    else if (elemType == VTK_TRIANGLE)
    {
        nen = 3;
        elem = new Triangle[nec];
    }

    xyz = new double [nnc*nsd];

    MPI_Barrier(MPI_COMM_WORLD);

    readmxyz(config.getMxyzFile());
    readmien(config.getMienFile());

    getTimesteps(config.getDtFile(), config.getDt(), config.getNrec());
    formLocalNodeList();
    localizeNodeCoordinates();

    MPI_Barrier(MPI_COMM_WORLD); //Need a barrier to delete xyz
    delete[] xyz;

    dataG = new double [nnc*ndf];
    dataL = new double [nnl*ndf];

}

/*
* @brief: read data per timestep and write to vtk
* @details: Read -> Localize -> VTK
*/
void Mesh::write()
{
    if (nrec == 0) vtkVisualization(-1, title, outpath, elemType);

    for (int irec = 0; irec < nrec; irec++) 
    {
        readDataRecord(irec, nrecoffset, nrecstride, spacetime, dataFiles);
        localizeData();
        vtkVisualization(irec, title, outpath, elemType);
    }
}

/*
* @readminf: Read the MIXD minf file
* @param: [in] string: filename
* @param: [in] bool: spacetime
*/
void Mesh::readminf(std::string filename, bool spacetime)
{
    std::ifstream file;
    std::string dummy;
    
    char* readStream;
    double dummyDouble;

    std::string lineString;
    std::string dummyString;

    int mype, npes;
    MPI_Comm_rank(MPI_COMM_WORLD, &mype);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);

    /// For parallel file input

    MPI_Status status;
    MPI_Offset offset;                      // offset from the beginning of file for parallel file read
    MPI_Offset size;
    MPI_Datatype mxyzftype,mienftype;       // mpi datatype used in parallel file read
    MPI_File fileptr;                       // file pointer for parallel file read


    file.open(filename.c_str(), std::ios::in);
    if (file.is_open()==false)
    {
        std::cout << "Unable to open minf file for pe: " << mype << "! Aborting... " << std::endl;
        MPI_Finalize();
        exit(0);
    }

    while (!file.eof())
    {
        // Get a line and store in lineString
        getline(file, lineString, '\n');
        // If the first character of the line is not a '#'
        if (lineString.c_str()[0] != '#')
        {
            std::istringstream iss(lineString);
            iss >> dummyString;
            if(dummyString == "nn")
                iss >> nn;
            else if(dummyString == "ne")
                iss >> ne;
        }
    }

    // We are only interested in the semi-discrete part of the mesh
    if (spacetime == 1) nn = nn/2;

    if (mype==0)
    {
        std::cout << "> Number of mesh elements : " << ne << std::endl;
        std::cout << "> Number of mesh nodes : " << nn << std::endl;
        std::cout << "> Number of mesh nodes (space): " << nnspace << std::endl;
        std::cout << "> File read complete: minf" << std::endl;
    }
    file.close();
}

void Mesh::distribute()
{
    int mype, npes;
    MPI_Comm_rank(MPI_COMM_WORLD, &mype);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);

    // Determine nec, mec
    nec = (ne-1)/npes + 1;
    mec = nec;
    if ((mype+1)*mec > ne)
        nec = ne - mype*mec;
    if (nec < 0)
        nec = 0;

    // Determine nnc, mnc
    nnc = (nn-1)/npes + 1;
    mnc = nnc;
    if ((mype+1)*mnc > nn)
        nnc = nn - mype*mnc;
    if (nnc < 0)
        nnc = 0;

}

void Mesh::readmxyz(std::string filename)
{
    std::ifstream file;
    char* readStream;

    int mype, npes;
    MPI_Comm_rank(MPI_COMM_WORLD, &mype);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);

    /// For parallel file input

    MPI_Status status;
    MPI_Offset offset;                      // offset from the beginning of file for parallel file read
    MPI_Offset size;
    MPI_Datatype mxyzftype;       // mpi datatype used in parallel file read
    MPI_File fileptr;                       // file pointer for parallel file read

    offset = mype*nsd*mnc*sizeof(double);
    MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fileptr);
    MPI_File_get_size(fileptr, &size);

    MPI_Type_contiguous(nnc*nsd, MPI_DOUBLE, &mxyzftype);
    MPI_Type_commit(&mxyzftype);

    if ((size != nn * nsd * sizeof(double)) && (size != 2 * nn * nsd * sizeof(double))) 
    {
        MPI_File_close(&fileptr);
        MPI_Barrier(MPI_COMM_WORLD);

        if(mype == 0) std::cout << "ERROR: MXYZ file size mismatch."  << std::endl;
        MPI_Finalize();
        exit(-1);
    }

    MPI_File_set_view(fileptr, offset, MPI_DOUBLE, mxyzftype, "native", MPI_INFO_NULL);
    readStream = new char [nsd*nnc*sizeof(double)];
    MPI_File_read(fileptr,readStream, nsd*nnc, MPI_DOUBLE, &status);
    swapBytes(readStream, nsd*nnc, sizeof(double));

    int xsd = 0;
    int ysd = 1;
    int zsd = 2;

    for(int i=0; i<nnc; i++)
    {
        /* node[i].setX(*((double*)readStream + nsd*i)); */
        /* node[i].setY(*((double*)readStream + nsd*i+1)); */
        /* node[i].setZ(*((double*)readStream + nsd*i+2)); */

        xyz[i*nsd+xsd] = *((double*)readStream + nsd*i);
        xyz[i*nsd+ysd] = *((double*)readStream + nsd*i+1);
        xyz[i*nsd+zsd] = *((double*)readStream + nsd*i+2);

    }

    if (mype==0) std::cout << "> File read complete: " << filename << std::endl;

    delete[] readStream;

    MPI_File_close(&fileptr);
    MPI_Type_free(&mxyzftype);
    MPI_Barrier(MPI_COMM_WORLD);

}

void Mesh::readmien(std::string filename)
{
    std::ifstream file;
    char* readStream;

    int mype, npes;
    MPI_Comm_rank(MPI_COMM_WORLD, &mype);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);

    MPI_Status status;
    MPI_Offset offset;                      // offset from the beginning of file for parallel file read
    MPI_Offset size;
    MPI_Datatype mienftype;                 // mpi datatype used in parallel file read
    MPI_File fileptr;                       // file pointer for parallel file read

    offset = mype*nen*mec*sizeof(int);

    MPI_Type_contiguous(nec*nen, MPI_INT, &mienftype);
    MPI_Type_commit(&mienftype);

    MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fileptr);
    MPI_File_get_size(fileptr, &size);

    std::cout << "size: " << size << std::endl;
    std::cout << "ne: " << ne << std::endl;
    std::cout << "nen: " << nen << std::endl;

    if (size != ne * nen * sizeof(int))
    {
        MPI_File_close(&fileptr);
        MPI_Barrier(MPI_COMM_WORLD);

        if(mype == 0) std::cout << "ERROR: MIEN file size mismatch."  << std::endl;
        MPI_Finalize();
        exit(-1);
    }

    MPI_File_set_view(fileptr, offset, MPI_INT, mienftype, "native", MPI_INFO_NULL);
    readStream = new char [nen*nec*sizeof(int)];
    MPI_File_read(fileptr, readStream, nen*nec, MPI_INT, &status);
    swapBytes(readStream, nen*nec, sizeof(int));
    int connValue;

    for(int i=0; i<nec; i++)
    {
        for(int j=0; j<nen; j++)
        {
            connValue = *((int*)readStream + nen*i+j);
            elem[i].setConn(j, connValue-1);
        }
    }

    if (mype==0) std::cout << "> File read complete: " << filename << std::endl;

    delete[] readStream;

    MPI_File_close(&fileptr);
    MPI_Type_free(&mienftype);
    MPI_Barrier(MPI_COMM_WORLD);
}


void Mesh::getTimesteps(std::string dtFile, double dt, int nrec)
{
    int mype;
    MPI_Comm_rank(MPI_COMM_WORLD, &mype);

    if (!dtFile.empty())
    {
        std::ifstream dtFileD;
        std::string lineString;
        char * pEnd;

        std::cout << "Using " << dtFile << "for timesteps" << std::endl;

        dtFileD.open(dtFile,std::ios::in);
        if (dtFileD.is_open()==false)
        {
            std::cout << "Unable to open input file for pe: " << mype << "! Aborting... " << std::endl;
            MPI_Finalize();
            exit(0);
        }

        while (!dtFileD.eof())
        {
            // Get a line and store in lineString
            getline(dtFileD, lineString, '\n');
            if (!lineString.empty())
                timesteps.push_back(std::strtod(lineString.c_str(), &pEnd));
        }

    }
    else
    {
        std::cout << "Assuming equal timesteps of size: " << dt << std::endl;
        for (auto it = 0; it != nrec; it++)
            timesteps.push_back(it * dt);
    }

}

/***************************************************************************************************
void tetMesh::formLocalConnectivity()
****************************************************************************************************
Every processor have only a portion of the elements and a portion of the connectivity.
Nodal info stored in a process is not necessarly belong to the elements stored in that p.
Therefore, we need to collect the necessary node information from other processes.
In order to do this, a local node list is produced with this function.

    7____8____9                                3____4____5
    |\   |\   |        e0 = 0,1,7              |\   |\   |        e0 = 0,1,3
    | \  | \  |        e1 = 1,8,7        ==>   | \  | \  |        e1 = 1,4,3
    |  \ |  \ |        e2 = 1,2,8              |  \ |  \ |        e2 = 1,2,4
    |___\|___\|        e3 = 2,9,8              |___\|___\|        e3 = 2,5,4
    0     1      2                             0    1    2

    Step 1:
    allLocalNodes = 0,1,7,1,8,7,1,2,8,2,9, 8
    index         = 0,1,2,3,4,5,6,7,8,9,10,11

    Step 2: sorts
    allLocalNodes = 0,1,1,1,2,2,7,7,8,8,8, 9
    index         = 0,1,3,6,7,9,2,5,4,8,11,10

    Step 3:
    nnl = 6
    rawLocalConn  = 0,1,1,1,2,2,3,3,4,4,4,5

    Step 4:
    nodeLToG = 0,1,2,7,8,9

    Step 5:
    allLocalNodes = 0,1,3,1,4,3,1,2,4,2,5,4

***************************************************************************************************/
void Mesh::formLocalNodeList()
{
    int mype, npes;                     // my processor rank and total number of processors
    MPI_Comm_rank(MPI_COMM_WORLD, &mype);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    int n_allLocalNodes = nec*nen;                      //number of all local nodes
    int * allLocalNodes = new int [n_allLocalNodes];    //array for all local node values
    int * index = new int [n_allLocalNodes];
    int * rawLocalConn = new int [n_allLocalNodes];

    /*
     * 1. Prepare allLocalNodes
     * All the node numbers of all the elements on processor are gathered in a single array.
     * Obviously, there are some duplicate node values in this array. We will get rid of the
     * duplicates by first sorting this array and then selecting the unique values.
    */
    for(int i=0; i<nec; i++)
        for(int j=0; j<nen; j++)
            allLocalNodes[i*nen+j] = elem[i].getConn(j);

    /*
     * 2. Sort allLocalNodes
     * Now we can sort our allLocalNodes array.
     * An index array accompanying the allLocalNodes array will also be sorted.
     * With the index array it is easier to form the local connectivity later.
     */
    for(int i=0; i<n_allLocalNodes; i++)
        index[i] = i;
    quickSort(allLocalNodes, index, 0, n_allLocalNodes-1);

    /*
     * 3. Determine nnl (number of unique local nodes) and form raw connectivity numbering.
     * This means that for each node in the sorted allLocalNodes array, we give a local node number.
     */
    nnl = 1;
    rawLocalConn[0] = 0;
    for(int i=1; i<n_allLocalNodes; i++)
    {
        if (allLocalNodes[i-1] != allLocalNodes[i])
            nnl++;
        rawLocalConn[i] = nnl-1;
    }

    /*
     * 4. Collect local to global node number converter: nodeLToG (ieng)
     * Now we know the number of unique local nodes(nnl). We can collect the unique node values.
     */
    nodeLToG = new int [nnl];            //allocate space for nodeLToG array
    nodeLToG[0] = allLocalNodes[0];        //first value is already known
    int i_unique = 0;                    //iterator on nodeLToG
    for(int i=1; i<n_allLocalNodes; i++)
        if (nodeLToG[i_unique] != allLocalNodes[i])
        {
            i_unique++;
            nodeLToG[i_unique] = allLocalNodes[i];
        }

    /*
     * 5. Form local connectivity
     * Using the index array, allLocalNodes is overwritten with rawLocalConn
     */
    for(int i=0; i<n_allLocalNodes; i++)
        allLocalNodes[index[i]] = rawLocalConn[i];

    //here a new local level tetNode structure is created.
    lNode = new Node [nnl];

    // finally element level lConn array is filled with the allLocalNodes
    for(int i=0; i<nec; i++)
        for(int j=0; j<nen; j++)
            elem[i].setLConn(j,allLocalNodes[i*nen+j]);

/*    if (mype==0)
    for(int i=0; i<nec; i++)
    {
        for(int j=0; j<nen; j++)
            cout << elem[i].getLConn(j) << '\t';
        cout << endl;
    }
*/

    delete[] allLocalNodes;
    delete[] index;
    delete[] rawLocalConn;

    return;
}

void Mesh::localizeNodeCoordinates()
{
    int n;      // counter
    int in;     // node number in global
    int ipe;    // pe number where 'in' lies
    int inc;    // offset of in on ipe
    int npes;

    MPI_Comm_size(MPI_COMM_WORLD, &npes);

    double * buffer;
    buffer = new double [mnc*nsd];
    int nncTarget;

    MPI_Win win;

    int xsd = 0;
    int ysd = 1;
    int zsd = 2;

    MPI_Win_create(xyz, nnc*nsd*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
    MPI_Win_fence(0, win);

    for(int ipes=0; ipes<npes; ipes++)
    {
        for (int i = 0; i < mnc*nsd; i++) buffer[i] = 0;
        if ((ipes+1)*mnc > nn) nncTarget = nn - ipes*mnc;
        else nncTarget = mnc;
        // cout << "ipes" << ipes << "npes" << npes << endl;
        // cout << "mype " << mype << " ipes " << ipes << " npes " << npes << " nncTarget " << nncTarget << endl;
        MPI_Get(&buffer[0], nncTarget*nsd, MPI_DOUBLE, ipes, 0, nncTarget*nsd, MPI_DOUBLE, win);
        MPI_Win_fence(0, win);

        for(int inl=0; inl<nnl; inl++)
        {
            in = nodeLToG[inl];
            ipe = in/mnc;
            inc = in%mnc;
            if (ipes == ipe)
            {
                lNode[inl].setX(buffer[inc*nsd+xsd]);
                lNode[inl].setY(buffer[inc*nsd+ysd]);
                lNode[inl].setZ(buffer[inc*nsd+zsd]);
            }
        }
    }

    MPI_Win_fence(0, win);
    MPI_Win_free(&win);

    delete[] buffer;

    return;
}

void Mesh::readDataRecord(
    int irec,
    int nrecoffset,
    int nrecstride,
    bool spacetime,
    std::vector<std::string> dataFiles
)
{
    std::ifstream    file;            // file name object for serial file read
    char*       readStream;      // temperory var used for strings read from files

    /// these are general mpi parameters
    int mype, npes;              // my processor rank and total number of processors
    MPI_Comm_rank(MPI_COMM_WORLD, &mype);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);

    /// these guys are used for parallel file input
    MPI_Status status;
    MPI_Offset offset;           // offset from the beginning of file for parallel file read
    MPI_Offset size;
    MPI_File fileptr;            // file pointer for parallel file read

    MPI_Datatype dataftype;
    MPI_Type_contiguous(nnc*ndf, MPI_DOUBLE, &dataftype);
    MPI_Type_commit(&dataftype);

    /* vector<string> dataFiles = settings->getDataFiles(); */
    std::string filename;

    if (spacetime) nrecstride *= 2;
    if (spacetime) nrecoffset *= 2;

    getFileAndOffset(irec, nrecstride, nrecoffset, ndf, dataFiles, filename, offset);
    /* cout << "Reading File: " << filename << " starting at offset: " << offset <<endl; */
    if (mype == 0) std::cout << "rec " << irec << ": reading..." << std::endl;


    MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fileptr);
    MPI_File_get_size(fileptr, &size);

    if (offset >= size)
    {
        MPI_File_close(&fileptr);
        MPI_Barrier(MPI_COMM_WORLD);

        if(mype == 0) std::cout << "ERROR: Overstepping file limits: Ciao:)"  << std::endl;
        MPI_Finalize();
        exit(-1);
    }

    MPI_File_set_view(fileptr, offset, MPI_DOUBLE, dataftype, "native", MPI_INFO_NULL);
    readStream = new char [ndf*nnc*sizeof(double)];
    MPI_File_read(fileptr, readStream, ndf*nnc, MPI_DOUBLE, &status);
    swapBytes(readStream, ndf*nnc, sizeof(double));

    for(int i=0; i<nnc*ndf; i++) { dataG[i] = 0;}; //TODO: convert all indices to long

    for(int inc=0; inc<nnc; inc++)
    {
        for(int idf=0; idf<ndf; idf++)
        {
            dataG[ndf*inc+idf]= *((double*)readStream + ndf*inc+idf);
        }
    }

    /* if (mype==0) cout << "> File read complete: " << filename << endl; */

    delete[] readStream;

    MPI_File_close(&fileptr);
    MPI_Type_free(&dataftype);
    MPI_Barrier(MPI_COMM_WORLD);

    return;

}

void Mesh::getFileAndOffset(
    int irec, int nrecstride, int nrecoffset, int ndf,
    std::vector<std::string> dataFiles,
    std::string& dataFile,
    MPI_Offset& totalOffset
)
{

    int mype, npes;              // my processor rank and total number of processors
    MPI_Comm_rank(MPI_COMM_WORLD, &mype);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);

    MPI_Offset processOffset, timestepOffset, recordOffset, fileOffset, rawOffset;

    /* rawOffset is the offset assuming we have one single large file that we're timestepping over
     * since we deal with Restart files, data is duplicated at the start of each data file: It has one extra timestep at the beginning.
     *
     * So we calculate the rawOffset and subtract a certain fileOffset to account for stepping over files.
     * Now to calculate the index of each file is a doozy. We have:
     *      - indexR (Restart) pointing to a file assuming we neglect the first timestep for every data file except the first.
     *
     *  There's probably an easier way to do this... like considering the LAST timestep of EVERY file as redundant, but this is what I went with.
     */
    processOffset = mype*ndf*mnc*sizeof(double);
    timestepOffset = (irec*nrecstride)*ndf*nn*sizeof(double);
    recordOffset = (nrecoffset)*ndf*nn*sizeof(double);
    rawOffset = processOffset + timestepOffset + recordOffset;

    //calculate fileOffset based on rawOffset and filesizes
    std::vector<size_t> dataFileSizes, dataFileSizesC, dataFileSizesR;
    size_t size = 0, cumulativeFileSize = 0;
    for(std::vector<std::string>::iterator it = dataFiles.begin(); it != dataFiles.end(); ++it)
    {
        int index = it - dataFiles.begin();
        size = std::experimental::filesystem::file_size(*it);

        // // Uncomment the following and delete the top line if the libc wasn't
        // // compiled with experimental features. i.e. if std::experimental or
        // // std::filesystem isn't available for the above line
        // // Also include <sys/stat.h>
        // struct stat sb{};
        // if (!stat((*it).c_str(), &sb)) {
        //     cout << sb.st_size << endl;
        // } else {
        //     perror("stat");
        // }
        // size = sb.st_size;


        cumulativeFileSize += size;
        dataFileSizes.push_back(size);

        // Cumulative file size array.
        // Necessary for finding fileOffset.
        dataFileSizesC.push_back(cumulativeFileSize);

        //Cumulative file size array without accounting for first timestep
        //Necessary for finding indexR
        dataFileSizesR.push_back( cumulativeFileSize - (index)*nn*ndf*sizeof(double));

        /* cout << "File: " << *it << " | size: " << size << " | cumulative: " << cumulativeFileSize << endl; */
    }

    // Find the file containing the current record.
    // using dataFileSizesR, we don't account for the first timestep in restart files, pinpointing the right file for us
    std::vector<size_t>::iterator iFileR = upper_bound(dataFileSizesR.begin(), dataFileSizesR.end(), rawOffset);
    int fileIndexR = iFileR - dataFileSizesR.begin();

    if (fileIndexR == 0)
        fileOffset = 0;
    else
    {
        // fileOffset = ACTUAL cumulative sum of all previous filesizes - one record per file except the first.
        // dataFileSizesC is used here (as opposed to R) because we essentially want to reset the reference to the start of
        // the file pointed to by fileIndexR. ALL the data in the previous files is irrelevant.
        // The second term pushes the reference forward by 1 step, avoiding the redundant initial timestep in restart files.
        fileOffset = dataFileSizesC.at(fileIndexR-1) - (fileIndexR)*nn*ndf*sizeof(double);
    }

    //totalOffset. Probably better called finalOffset or trueOffset
    totalOffset = rawOffset - fileOffset;
    dataFile = dataFiles.at(fileIndexR);

    /* cout << "irec: " << irec << " file: " << "[" << fileIndexR << "] " << dataFiles.at(fileIndexR) << " rawOffset: " << rawOffset << " fileOffset: " << fileOffset << " totalOffset: " << totalOffset << endl; */

}

// Read scalar data and localize it to each process so that the data can be split.
void Mesh::localizeData()
{

    int n;      // counter
    int in;     // node number in global
    int ipe;    // pe number where 'in' lies
    int inc;    // offset of in on ipe
    int npes;

    MPI_Comm_size(MPI_COMM_WORLD, &npes);

    double * buffer;
    buffer = new double [mnc*ndf];
    int nncTarget;

    MPI_Win win;

    MPI_Win_create(dataG, nnc*ndf*sizeof(double), sizeof(double), MPI_INFO_NULL, MPI_COMM_WORLD, &win);
    MPI_Win_fence(0, win);

    for(int i=0; i<nnl*ndf; i++) { dataL[i] = 0;};

    int mype;                // my processor rank and total number of processors
    MPI_Comm_rank(MPI_COMM_WORLD, &mype);

    for(int ipes=0; ipes<npes; ipes++)
    {
        for (int i = 0; i < mnc*ndf; i++) buffer[i] = 0;
        if ((ipes+1)*mnc > nn) nncTarget = nn - ipes*mnc;
        else nncTarget = mnc;
        // cout << "ipes" << ipes << "npes" << npes << endl;
        // cout << "mype " << mype << " ipes " << ipes << " npes " << npes << " nncTarget " << nncTarget << endl;
        MPI_Get(&buffer[0], nncTarget*ndf, MPI_DOUBLE, ipes, 0, nncTarget*ndf, MPI_DOUBLE, win);
        MPI_Win_fence(0, win);

        for(int inl=0; inl<nnl; inl++)
        {
            in = nodeLToG[inl];
            ipe = in/mnc;
            inc = in%mnc;

            if (ipes == ipe)
            {
                for(int idf=0; idf<ndf; idf++)
                {
                    dataL[inl*ndf+idf] = buffer[inc*ndf+idf];
                }
            }
        }
    }

    MPI_Win_fence(0, win);
    MPI_Win_free(&win);

    delete[] buffer;

    return;
}

void Mesh::vtkVisualization(int irec, std::string title, std::string outpath, int VTKElemType)
{

    // int nn  = mesh->getNn();
    // int nnc = mesh->getNnc();
    // int nnl = mesh->getNnl();
    // int mnc = mesh->getMnc();
    // int ne  = mesh->getNe();
    // int nec = mesh->getNec();
    // int mec = mesh->getMec();
    // int ndf = mesh->getNdf();


    int mype, npes;             // my processor rank and total number of processors

    MPI_Comm_rank(MPI_COMM_WORLD, &mype);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    std::ostringstream int2str;
    int2str << mype;
    std::string strMype = int2str.str();
    std::string dummy;


    // VTK Double Array
    vtkSmartPointer<vtkDoubleArray> pcoords = vtkSmartPointer<vtkDoubleArray>::New();
    pcoords->SetNumberOfComponents(nsd);
    pcoords->SetNumberOfTuples(nnl);


    // vtkDoubleArray type pcoords is filled with the data in meshPoints.
    // Typically, I deal with meshes embedded in 3D space
    for (int i=0; i<nnl; i++)
        pcoords->SetTuple3(i,getLNode(i)->getX(),getLNode(i)->getY(),getLNode(i)->getZ());

    // NOTE: This is only necessary for 2D meshes embedded in 2D space
    // It is just easier to add a z coordinate to your mesh, with 0 value
    // by modifying your mxyz file, making your space 3D and nsd==3.
    /* for (int i=0; i<nnl; i++) */
    /*     pcoords->SetTuple2(i,mesh->getLNode(i)->getX(),mesh->getLNode(i)->getY()); */

    //vtkPoints type outputPoints is filled with the data in pcoords.
    vtkSmartPointer<vtkPoints> outputPoints = vtkSmartPointer<vtkPoints>::New();
    outputPoints->SetData(pcoords);


    //Connectivity is written to vtkCellArray type outputCells
    vtkSmartPointer<vtkCellArray> connectivity = vtkSmartPointer<vtkCellArray>::New();
    for(int i=0; i<nec; i++)
    {
        connectivity->InsertNextCell(nen);
        for(int j=0; j<nen; j++)
            connectivity->InsertCellPoint(getElem(i)->getLConn(j));
    }


    vtkSmartPointer<vtkUnstructuredGrid> unsGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    unsGrid->SetPoints(outputPoints);
    unsGrid->SetCells(VTKElemType, connectivity);

    if (irec != -1)
    {

        vtkSmartPointer<vtkFieldData> TimeValue = vtkSmartPointer<vtkFieldData>::New();
        vtkSmartPointer<vtkDoubleArray> TimeValueArray = vtkSmartPointer<vtkDoubleArray>::New();
        TimeValueArray->SetName("TimeValue");
        TimeValueArray->SetNumberOfTuples(1);
        TimeValueArray->SetNumberOfValues(1);
        TimeValueArray->SetTuple1(0, timesteps[irec]);
        TimeValue->AddArray(TimeValueArray);
        unsGrid->SetFieldData(TimeValue);

        double * dataL = getDataL();
        std::vector<vtkSmartPointer<vtkDoubleArray>> scalars(ndf);
        for(int idf=0; idf < ndf ; idf++)
        {
            scalars[idf] = vtkSmartPointer<vtkDoubleArray>::New();
            dummy = "scalar_";
            dummy.append(std::to_string(idf));
            scalars[idf]->SetName(dummy.c_str());
            for(int inl=0; inl<nnl; inl++)
                scalars[idf]->InsertNextValue(dataL[inl*ndf + idf]);
            unsGrid->GetPointData()->AddArray(scalars[idf]);
            unsGrid->GetPointData()->SetActiveScalars(dummy.c_str());
        }
    }

    vtkSmartPointer<vtkXMLPUnstructuredGridWriter> pwriter = vtkSmartPointer<vtkXMLPUnstructuredGridWriter>::New();

    dummy = outpath + "/";

    dummy.append(title);
    dummy.append("_");
    dummy.append(std::to_string(irec));
    dummy.append(".pvtu");
    pwriter->SetFileName(dummy.c_str());
    pwriter->SetNumberOfPieces(npes);
    pwriter->SetStartPiece(mype);
    pwriter->SetEndPiece(mype);
    pwriter->SetInputData(unsGrid);
    /* pwriter->Write(); */
    pwriter->Update();

    return;
}
