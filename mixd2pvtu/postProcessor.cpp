#include "postProcessor.h"

/***************************************************************************************************
  preProcessorControl
 **************************************************************************************************/

void postProcessor::postProcessorControl(inputSettings* argSettings, tetMesh* argMesh, int irec)
{
    int mype, npes;             // my processor rank and total number of processors
    MPI_Comm_rank(MPI_COMM_WORLD, &mype);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    mesh = argMesh;
    settings = argSettings;

    /* if (mype==0) cout << "Writing irec: " << irec << endl; */
    if (mype == 0) cout << "rec " << irec << ": writing...\r" << flush;
    vtkVisualization(irec);
    if (mype == 0) cout << "rec " << irec << ": done!     " << endl;;

    return;
}

/***************************************************************************************************
// Main visualization function
 ***************************************************************************************************/
void postProcessor::vtkVisualization(int irec)
{
    int nn  = mesh->getNn();
    int nnc = mesh->getNnc();
    int nnl = mesh->getNnl();
    int mnc = mesh->getMnc();
    int ne  = mesh->getNe();
    int nec = mesh->getNec();
    int mec = mesh->getMec();
    int ndf = mesh->getNdf();

    double * dataL = mesh->getDataL();
    double * dataG = mesh->getDataG();


    int * nodeLToG = mesh->getNodeLToG();

    int mype, npes;             // my processor rank and total number of processors

    MPI_Comm_rank(MPI_COMM_WORLD, &mype);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    ostringstream int2str;
    int2str << mype;
    string strMype = int2str.str();
    string dummy;


    // VTK Double Array
    vtkSmartPointer<vtkDoubleArray> pcoords = vtkSmartPointer<vtkDoubleArray>::New();
    pcoords->SetNumberOfComponents(nsd);
    pcoords->SetNumberOfTuples(nnl);


    //vtkDoubleArray type pcoords is filled with the data in meshPoints.
    for (int i=0; i<nnl; i++)
        pcoords->SetTuple3(i,mesh->getLNode(i)->getX(),mesh->getLNode(i)->getY(),mesh->getLNode(i)->getZ());

    //vtkPoints type outputPoints is filled with the data in pcoords.
    vtkSmartPointer<vtkPoints> outputPoints = vtkSmartPointer<vtkPoints>::New();
    outputPoints->SetData(pcoords);


    //Connectivity is written to vtkCellArray type outputCells
    vtkSmartPointer<vtkCellArray> connectivity = vtkSmartPointer<vtkCellArray>::New();
    for(int i=0; i<nec; i++)
    {
        connectivity->InsertNextCell(nen);
        for(int j=0; j<nen; j++)
            connectivity->InsertCellPoint(mesh->getElem(i)->getLConn(j));
    }

    vtkSmartPointer<vtkUnstructuredGrid> unsGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
    unsGrid->SetPoints(outputPoints);
    unsGrid->SetCells(10,connectivity);


    vector<vtkSmartPointer<vtkDoubleArray>> scalars(ndf);
    for(int idf=0; idf < ndf ; idf++)
    {
        scalars[idf] = vtkSmartPointer<vtkDoubleArray>::New();
        dummy = "scalar_";
        dummy.append(to_string(idf));
        scalars[idf]->SetName(dummy.c_str());
        for(int inl=0; inl<nnl; inl++)
            scalars[idf]->InsertNextValue(dataL[inl*ndf + idf]);
        unsGrid->GetPointData()->AddArray(scalars[idf]);
        unsGrid->GetPointData()->SetActiveScalars(dummy.c_str());
    }


    // Now we write the "Title.pvtu" file which contains the informtaiton about other files.
    // Can be used to add more prefixes to the piece filenames. Folders/subfolders and such.
    // Mostly unnecessary
    /* auto pwriter = vtkSmartPointer<LDEMVTKXMLPUnstructuredDataWriter>::New(); */

    vtkSmartPointer<vtkXMLPUnstructuredGridWriter> pwriter = vtkSmartPointer<vtkXMLPUnstructuredGridWriter>::New();


    dummy = "output/";
    dummy.append(settings->getTitle());
    dummy.append("_");
    dummy.append(to_string(irec));
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
