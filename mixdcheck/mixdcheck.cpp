
#include "../mixdclass/mixd.hpp"
#include <set>
#include <assert.h>


void printUsage(char *binaryName)
{
    std::cout << "MIXD check tool (operates on mesh in directory of call)" << std::endl;
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

    cout << "MIXD check tool..." << endl << endl;

    try{


    long ne, nn;
    mixd::readminf("./minf", &nn, &ne);


    mixd::MixdFile<int> mien("./mien", ne, nen);

    mien.read();

    std::set<int> foundNodes(mien.data(), mien.data() + nen*ne);

    cout << "done!" << endl;

    cout << "Found " << foundNodes.size() << " nodes used in the given mien file." << endl;


    cout << endl;

    if(spacetime)
        assert(foundNodes.size() == nn/2);
    else
        assert(foundNodes.size() == nn);


    } catch(mixd::MixdException e)
    { std::cout << e.msg() << std::endl; }



    return 0;
}
