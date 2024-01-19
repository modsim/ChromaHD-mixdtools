#include "mixd.hpp"


void printUsage(char *binaryName)
{
    std::cout << "MIXD tool to map semi-discrete flow solution to complete column" << std::endl;
    std::cout << "Run the binary in the directory with full mesh minf" << std::endl;
    std::cout << "Usage: " << binaryName << " [-tri|-tet|-qua|-hex] <flow mesh dir> <flow solve dir> <output flowfield dir>" << std::endl;
    std::cout << "This tool assumes nmap file between full mesh and flow mesh exists in the full mesh directory, and is named nmap" << std::endl;
}


int main(int argc, char **argv)
{
    using namespace std;
    using namespace mixd;

    if(argc != 5)
    {
        printUsage(argv[0]);
        return 1;
    }

    string etype(argv[1]);
    string meshdir(argv[2]);
    string solvdir(argv[3]);
    string outdir(argv[4]);


    int nsd = 0;

    if(etype == "-tri")
    {
        nsd = 2;
    }
    else if(etype == "-tet")
    {
        nsd = 3;
    }
    else if(etype == "-qua")
    {
        nsd = 2;
    }
    else if(etype == "-hex")
    {
        nsd = 3;
    }
    else if(etype == "-tetP2")
    {
        nsd = 3;
    }
    else
    {
        printUsage(argv[0]);
        return 1;
    }


    cout << "MIXD mapflow tool..." << endl << endl;


    try{

    long ne_in, nn_in;
    readminf(meshdir+"/minf", &nn_in, &ne_in);

    MixdFile<double> data(solvdir+"/data.out", 2*nn_in, nsd+1);
    data.read();


    long ne, nn;
    readminf("minf", &nn, &ne);

    MixdFile<double> flow(outdir+"/flowfield", nn, nsd+1);
    flow.init();

    MixdFile<int> nmap("nmap", nn);
    nmap.read();


    for(long i=0; i<nn; i++)
    {
        if(nmap(i) > 0)
        {
            for(int j=0; j<nsd+1; j++)
                flow(i,j) = data(nn_in+nmap(i)-1,j);
        }
    }

    flow.write();


    } catch(mixd::MixdException e)
    { std::cout << e.msg() << std::endl; }



    return 0;
}
