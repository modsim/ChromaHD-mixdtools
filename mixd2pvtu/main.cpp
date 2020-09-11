/***************************************************************************************************
Name        : 2D_Unsteady_Diffusion.cpp
 ***************************************************************************************************/

#include "settings.h"
#include "tet.h"
#include "postProcessor.h"
#include "vtkMPIController.h"
#include <sys/stat.h>
#include <sys/types.h>

int main(int argc, char **argv) {
    /***************************************************************************************************
      MAIN PROGRAM FLOW
      1. Pre-Processing Stage
      1.1. Settings
      1.2. Mesh
      3. Post-Processing Stage
     ***************************************************************************************************/

    int mype, npes;
    double starttime;

    inputSettings*  settings    = new inputSettings(argc, argv);
    tetMesh*        mesh        = new tetMesh;
    postProcessor*  postP       = new postProcessor;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mype);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);

    vtkMPIController* controller = vtkMPIController::New();
    controller->Initialize(&argc, &argv, 1);
    controller->SetGlobalController(controller);

    if(mype == 0) cout << "MIXD to PVTU converter."  << endl;


    // Pre-Processing Stage
    starttime = MPI_Wtime(); //timer starts
    settings->prepareSettings();
    mesh->prepareMesh(settings);

    if (mype==0) cout << endl << "================ POST-PROCESSING =================" << endl;

    int check = mkdir("output", 0777);
    if (!check)
        printf("Output directory created\n");
    else {
        printf("Output directory already exists\n");
    }

    int nrec = settings->getNrec();
    for (int irec=0; irec<nrec; irec++)
    {
        mesh->processData(settings, irec);
        postP->postProcessorControl(settings, mesh, irec);
    }


    if (mype==0) cout << endl << "====================== END =======================" << endl;
    if(mype==0) cout << "Elapsed time is " << fixed << MPI_Wtime()-starttime << endl;

    // Cleanup
    delete settings;
    delete mesh;
    delete postP;
    if(mype == 0) cout << ": Ciao:)"  << endl;

    controller->Finalize(1);
    MPI_Finalize();

    return 0;
}
