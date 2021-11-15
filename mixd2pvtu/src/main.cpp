#include "config.hpp"
#include "mesh.hpp"

#include <vtk/vtkMPIController.h>
#include <sys/stat.h>
#include <sys/types.h>

int main(int argc, char* argv[]) 
{
    int mype;
    int npes;
    double starttime;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mype);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    vtkMPIController* controller = vtkMPIController::New();
    controller->Initialize(&argc, &argv, 1);
    controller->SetGlobalController(controller);

    if(mype == 0) cout << "MIXD to PVTU converter."  << endl;
    starttime = MPI_Wtime(); //timer starts

    // NOTE: All processes read config file at once
    Config config(argc, argv);

    Mesh *mesh = new Mesh(config);

    mesh->write();

    if(mype==0) std::cout << "Elapsed time is " << std::fixed << MPI_Wtime()-starttime << std::endl;
    if(mype == 0) cout << ": Ciao:)"  << endl;

    // delete mesh;

    controller->Finalize(1);
    MPI_Finalize();

    return 0;
}
