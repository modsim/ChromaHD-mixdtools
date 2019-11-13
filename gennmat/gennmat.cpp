
#include "mixd.hpp"


void printUsage(char *binaryName)
{
    std::cout << "MIXD gennmat tool (operates on mesh in directory of call)" << std::endl;
    std::cout << "Usage: " << binaryName << " <nen> <st/sd>" << std::endl;
}


int main(int argc, char **argv)
{
    using namespace std;

    if(argc != 3)
    {
        printUsage(argv[0]);
        return 1;
    }

    const int nen = atoi(argv[1]);
    bool spacetime;

    if(strcmp(argv[2], "st") == 0) spacetime = true;
    else if(strcmp(argv[2], "sd") == 0) spacetime = false;
    else return 1;

    cout << "MIXD gennmat tool..." << endl << endl;

    try{


    long ne, nn;
    mixd::readminf("./minf", &nn, &ne);


    mixd::MixdFile<int> mien("./mien", ne, nen);
    mixd::MixdFile<int> mmat("./mmat", ne);
    mixd::MixdFile<int> nmat("./nmat", nn);

    mien.read();
    mmat.read();


    cout << "Finding nodal materials... " << flush;

    int nid;
    for(int i=0; i<ne; i++)
    {
        for(int n=0; n<nen; n++)
        {
            nid = mien(i,n) - 1;

            nmat(nid) = mmat(i);

            if(spacetime)
                nmat(nid + nn/2) = mmat(i);
        }
    }
    cout << "done!" << endl;

    cout << endl;

    nmat.write();

    } catch(mixd::MixdException e)
    { std::cout << e.msg() << std::endl; }



    return 0;
}
