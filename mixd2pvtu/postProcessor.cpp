#include "postProcessor.h"

/*********************************************************************************************
 * Control function for postprocessor class.
**********************************************************************************************/
void postProcessor::postProcessorControl(inputSettings* argSettings, tetMesh* argMesh, int irec)
{
    int mype, npes;             // my processor rank and total number of processors
    MPI_Comm_rank(MPI_COMM_WORLD, &mype);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    mesh = argMesh;
    settings = argSettings;

    if (mype == 0) cout << "\rrec " << irec << ": writing..." << flush;
    vtkVisualization(irec);
    if (mype == 0) cout << "\rrec " << irec << ": done!     " << endl;;

    return;
}

/*********************************************
 * Main visualization function: writes pvtu files
*********************************************/
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


    if (irec != -1)
    {

        vtkSmartPointer<vtkFieldData> TimeValue = vtkSmartPointer<vtkFieldData>::New();
        vtkSmartPointer<vtkDoubleArray> TimeValueArray = vtkSmartPointer<vtkDoubleArray>::New();
        TimeValueArray->SetName("TimeValue");
        TimeValueArray->SetNumberOfTuples(1);
        TimeValueArray->SetNumberOfValues(1);
        TimeValueArray->SetTuple1(0, mesh->timesteps[irec]);
        TimeValue->AddArray(TimeValueArray);
        unsGrid->SetFieldData(TimeValue);

        double * dataL = mesh->getDataL();
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
    }

    vtkSmartPointer<vtkXMLPUnstructuredGridWriter> pwriter = vtkSmartPointer<vtkXMLPUnstructuredGridWriter>::New();

    /* dummy = "output/" */
    dummy = settings->getOutpath() + "/";

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
