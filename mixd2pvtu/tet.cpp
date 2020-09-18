#include "tet.h"
#include <vector>
#include <iterator>
#include <experimental/filesystem>
#include <algorithm>

/***************************************************************************************************
void preProcessor::prepareMesh()
****************************************************************************************************
Does everythin' realted to mesh generation and transfers between processors.
***************************************************************************************************/
void tetMesh::prepareMesh(inputSettings* settings)
{
    int mype;                // my processor rank and total number of processors
    MPI_Comm_rank(MPI_COMM_WORLD, &mype);
    if (mype==0) cout << endl << "====================== MESH ======================" << endl;

    readMeshFiles(settings);
    formLocalNodeList();
    localizeNodeCoordinates();

    MPI_Barrier(MPI_COMM_WORLD); //Need a barrier to delete xyz
    delete[] xyz;

    dataG = new double [nnc*ndf];
    dataL = new double [nnl*ndf];

    return;
}

void tetMesh::processData(inputSettings* settings, int irec)
{
    //given irec, find file and offset
    // first, list number of files = nDataFiles
    // find sizes of each files and put it in a vector: vDataFileSizes


    readDataFile(settings, irec);
    localizeData();

}

void tetMesh::getFileAndOffset(inputSettings* settings, int irec, string& dataFile, MPI_Offset& totalOffset)
{
    int nrecstride = settings->getNrecstride();
    if (settings->getSpacetime() == 1) nrecstride *= 2;
    int nrecoffset = settings->getNrecoffset();
    if (settings->getSpacetime() == 1) nrecoffset *= 2;

    int mype, npes;              // my processor rank and total number of processors
    MPI_Comm_rank(MPI_COMM_WORLD, &mype);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    int ndf = settings->getNdf();

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
    vector<string> dataFiles = settings->getDataFiles();
    vector<size_t> dataFileSizes, dataFileSizesC, dataFileSizesR;
    size_t size = 0, cumulativeFileSize = 0;
    for(vector<string>::iterator it = dataFiles.begin(); it != dataFiles.end(); ++it)
    {
        int index = it - dataFiles.begin();
        size = experimental::filesystem::file_size(*it);
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
    vector<size_t>::iterator iFileR = upper_bound(dataFileSizesR.begin(), dataFileSizesR.end(), rawOffset);
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

/***************************************************************************************************
void preProcessor::readMeshFiles()
****************************************************************************************************
File read procedure :
1- Name of the file to be opened is retrieved from the inputSetting obj.
2- File is opened in appropriate format, this is ascii format for minf or
   other text files and binary format for binary mesh files.
3- Read operation for minf file is straight forward. Binary files are read
   as size of a double or int and stored in readStream. Then swapbytes
   function is called to swap the bytes for the correct endianness.
4- Finally obtained data is deep-copied to the mesh data structure.
***************************************************************************************************/
void tetMesh::readMeshFiles(inputSettings* settings)
{
    ifstream    file;            // file name object for serial file read
    string      dummy;           // dummy string to hold file names etc.
    char*       readStream;      // temperory var used for strings read from files
    double      dummyDouble;     // temperory var used for double values read from files

    /// these are general mpi parameters
    int mype, npes;              // my processor rank and total number of processors
    MPI_Comm_rank(MPI_COMM_WORLD, &mype);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);

    /// these guys are used for parallel file input
    MPI_Status status;
    MPI_Offset offset;           // offset from the beginning of file for parallel file read
    MPI_Offset size;
    MPI_Datatype mxyzftype,mienftype;        // mpi datatype used in parallel file read
    MPI_File fileptr;            // file pointer for parallel file read



    /***********************************************************************************************
    * READ THE MINF FILE
    * This file should hold the number of elements, nodes, space dimensions, element nodes and
    * element faces.
    ***********************************************************************************************/
    dummy = settings->getMinfFile();
    file.open(dummy.c_str(), ios::in);
    if (file.is_open()==false)
    {
        cout << "Unable to open minf file for pe: " << mype << "! Aborting... " << endl;
        MPI_Finalize();
        exit(0);
    }

    string lineString;
    string dummyString;

    while (!file.eof())
    {
        // Get a line and store in lineString
        getline(file, lineString, '\n');
        // If the first character of the line is not a '#'
        if (lineString.c_str()[0] != '#')
        {
            istringstream iss(lineString);
            iss >> dummyString;
            if(dummyString == "nn")
                iss >> nn;
            else if(dummyString == "ne")
                iss >> ne;
        }
    }

    if (settings->getSpacetime() == 1) nn = nn/2;

    if (mype==0)
    {
        cout << "> Number of mesh elements : " << ne << endl;
        cout << "> Number of nodes : " << nn << endl;
        cout << "> File read complete: minf" << endl;
    }
    file.close();


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

    //cout << "mype:" << mype << " nnc:" << nnc << " mnc:" << mnc << " nec:" << nec << endl;

    ndf = settings->getNdf();

    //Allocation of memory for the mesh data structure
    /* cout << "Allocating " << nnc*nsd << " doubles for coordinate data: " << double(nnc)*nsd*sizeof(double) / (1024*1024) << "MB" << endl; */
    xyz = new double [nnc*nsd];

    /* cout << "Allocating " << nnc*ndf << " doubles for scalar data: " << double(nnc)*ndf*sizeof(double) / (1024*1024)<< "MB" << endl; */

    /* cout << "Allocating " << nnc << " tetNodes for coordinate data: " << double(nnc)*sizeof(tetNode) / (1024*1024) << "MB" << endl; */
    /* cout << "Allocating " << nec << " tetElements for coordinate data: " << double(nec)*sizeof(tetElement) / (1024*1024) << "MB" << endl; */
    /* node = new tetNode[nnc]; */
    elem = new tetElement[nec];

    /* if (mype==0) cout << "> Mesh data structure is created." << endl; */
    MPI_Barrier(MPI_COMM_WORLD);


    MPI_Type_contiguous(nnc*nsd, MPI_DOUBLE, &mxyzftype);
    MPI_Type_contiguous(nec*nen, MPI_INT, &mienftype);
    MPI_Type_commit(&mxyzftype);
    MPI_Type_commit(&mienftype);

    /***********************************************************************************************
    * READ THE MXYZ FILE
    * This file contains the node coordinates
    ***********************************************************************************************/
    dummy = settings->getMxyzFile();
    /* char * writable = new char[dummy.size() + 1]; */
    /* std::copy(dummy.begin(), dummy.end(), writable); */
    /* writable[dummy.size()] = '\0'; */

    offset = mype*nsd*mnc*sizeof(double);
    MPI_File_open(MPI_COMM_WORLD, dummy.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fileptr);
    MPI_File_get_size(fileptr, &size);

    if (size < nn * nsd * sizeof(double))
    {
        MPI_File_close(&fileptr);
        MPI_Barrier(MPI_COMM_WORLD);

        if(mype == 0) cout << "ERROR: MXYZ file is smaller than expected."  << endl;
        MPI_Finalize();
        exit(-1);
    }

    MPI_File_set_view(fileptr, offset, MPI_DOUBLE, mxyzftype, "native", MPI_INFO_NULL);
    readStream = new char [nsd*nnc*sizeof(double)];
    MPI_File_read(fileptr,readStream, nsd*nnc, MPI_DOUBLE, &status);
    swapBytes(readStream, nsd*nnc, sizeof(double));
    for(int i=0; i<nnc; i++)
    {
        /* node[i].setX(*((double*)readStream + nsd*i)); */
        /* node[i].setY(*((double*)readStream + nsd*i+1)); */
        /* node[i].setZ(*((double*)readStream + nsd*i+2)); */
        xyz[i*nsd+xsd] = *((double*)readStream + nsd*i);
        xyz[i*nsd+ysd] = *((double*)readStream + nsd*i+1);
        xyz[i*nsd+zsd] = *((double*)readStream + nsd*i+2);

    }

    if (mype==0) cout << "> File read complete: " << dummy << endl;

    delete[] readStream;

    MPI_File_close(&fileptr);
    MPI_Barrier(MPI_COMM_WORLD);

    /***********************************************************************************************
    * READ THE MIEN FILE
    * This file contains the element connectivity
    ***********************************************************************************************/
    dummy = settings->getMienFile();
    /* writable = new char[dummy.size() + 1]; */
    /* std::copy(dummy.begin(), dummy.end(), writable); */
    /* writable[dummy.size()] = '\0'; */

    offset = mype*nen*mec*sizeof(int);
    MPI_File_open(MPI_COMM_WORLD, dummy.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fileptr);
    MPI_File_get_size(fileptr, &size);

    if (size < ne * nen * sizeof(int))
    {
        MPI_File_close(&fileptr);
        MPI_Barrier(MPI_COMM_WORLD);

        if(mype == 0) cout << "ERROR: MIEN file is smaller than expected."  << endl;
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

    if (mype==0) cout << "> File read complete: " << dummy << endl;

    delete[] readStream;

    MPI_File_close(&fileptr);

    MPI_Type_free(&mxyzftype);
    MPI_Type_free(&mienftype);

    MPI_Barrier(MPI_COMM_WORLD);



    return;
}

void tetMesh::readDataFile(inputSettings* settings, int irec)
{
    ifstream    file;            // file name object for serial file read
    string      dummy;           // dummy string to hold file names etc.
    char*       readStream;      // temperory var used for strings read from files
    double      dummyDouble;     // temperory var used for double values read from files

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
    /***********************************************************************************************
    * READ THE DATA FILE
    * This file contains output scalar data
    ***********************************************************************************************/
    /* vector<string> dataFiles = settings->getDataFiles(); */
    string filename;

    int nrecstride = settings->getNrecstride();
    if (settings->getSpacetime() == 1) nrecstride *= 2;
    int nrecoffset = settings->getNrecoffset();
    if (settings->getSpacetime() == 1) nrecoffset *= 2;

    getFileAndOffset(settings, irec, filename, offset);
    /* cout << "Reading File: " << filename << " starting at offset: " << offset <<endl; */
    if (mype == 0) cout << "rec " << irec << ": reading...\r" << flush;


    MPI_File_open(MPI_COMM_WORLD, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fileptr);
    MPI_File_get_size(fileptr, &size);

    if (offset >= size)
    {
        MPI_File_close(&fileptr);
        MPI_Barrier(MPI_COMM_WORLD);

        if(mype == 0) cout << "ERROR: Overstepping file limits: Ciao:)"  << endl;
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

void tetMesh::swapBytes (char *array, int nelem, int elsize)
{
    int sizet, sizem, i, j;
    char *bytea, *byteb;
    sizet = elsize;
    sizem = sizet - 1;
    bytea = new char [sizet];
    byteb = new char [sizet];
    for (i = 0; i < nelem; i++)
    {
        memcpy((void *)bytea, (void *)(array+i*sizet), sizet);
        for (j = 0; j < sizet; j++)
            byteb[j] = bytea[sizem - j];
        memcpy((void *)(array+i*sizet), (void *)byteb, sizet);
    }
    delete[] bytea;
    delete[] byteb;

    return;
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
void tetMesh::formLocalNodeList()
{
    int mype, npes;                // my processor rank and total number of processors
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
    lNode = new tetNode [nnl];

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

void tetMesh::localizeNodeCoordinates()
{
    int n; //counter
    int in;     // node number in global
    int ipe;    // pe number where 'in' lies
    int inc;    // offset of in on ipe
    int npes;

    MPI_Comm_size(MPI_COMM_WORLD, &npes);

    double * buffer;
    buffer = new double [mnc*nsd];
    int nncTarget;

    MPI_Win win;

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


/***************************************************************************************************
void tetMesh::quickSort()
****************************************************************************************************
Necessary while forming of the local connectivity array.
***************************************************************************************************/
void tetMesh::quickSort(int* arr, int* index, int left, int right)
{
    int i = left, j = right;
    int tmp, tmpIndex;
    int pivot = arr[(left + right) / 2];

    /* partition */
    while (i <= j)
    {
        while (arr[i] < pivot)
            i++;
        while (arr[j] > pivot)
            j--;
        if (i <= j)
        {
            tmp = arr[i];       tmpIndex = index[i];
            arr[i] = arr[j];    index[i] = index[j];
            arr[j] = tmp;       index[j] = tmpIndex;
            i++;
            j--;
        }
    };

    /* recursion */
    if (left < j)
        quickSort(arr, index, left, j);
    if (i < right)
        quickSort(arr, index, i, right);

    return;
}

void tetMesh::localizeData()
{
    int n; //counter
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

