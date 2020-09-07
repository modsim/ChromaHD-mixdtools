/***************************************************************************************************
Name        : 2D_Unsteady_Diffusion.cpp
***************************************************************************************************/

#include "settings.h"
#include "tet.h"
#include "postProcessor.h"

int main(int argc, char **argv) {
/***************************************************************************************************
MAIN PROGRAM FLOW
1. Pre-Processing Stage
    1.1. Settings
    1.2. Mesh
2. Solution Stage
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
    if(mype == 0) cout << "MIXD to PVTU converter."  << endl;

    // Pre-Processing Stage
    starttime = MPI_Wtime(); //timer starts
    settings->prepareSettings();
    mesh->prepareMesh(settings);

    int nrec = settings->getNrec();
    for (int irec=0; irec<nrec; irec++)
    {
        mesh->processData(settings, irec);
        postP->postProcessorControl(settings, mesh, irec);
    }


    if(mype==0) cout << "Elapsed time is " << fixed << MPI_Wtime()-starttime << endl;

    // Cleanup
    delete settings;
    delete mesh;
    delete postP;
    if(mype == 0) cout << ": Ciao:)"  << endl;
    MPI_Finalize();

    return 0;
}
