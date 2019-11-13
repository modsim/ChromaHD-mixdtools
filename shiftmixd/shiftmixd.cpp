
#include "mixd.hpp"

#include <iostream>
#include <cstdlib>


void printUsage(char *binaryName)
{
    std::cout << "MIXD mesh shifting tool (operates on mxyz in directory of call)" << std::endl;
    std::cout << "Usage: " << binaryName << " [-tri|-tet|-qua|-hex] <dx> <dy> [<dz>]" << std::endl;
}


int main(int argc, char **argv)
{
    using namespace std;
    using namespace mixd;
    
    if(argc != 4 && argc != 5)
    {
        printUsage(argv[0]);
        return 1;
    }
    
    string etype(argv[1]);
    
    const double dx = atof(argv[2]);
    const double dy = atof(argv[3]);
    const double dz = (argc==5) ? atof(argv[4]) : 0.0;
    
    
    cout << "MIXD shifting tool..." << endl << endl;
    
    try{
    
    long ne, nn;
    readminf("./minf", &nn, &ne);
    
    if(nn<=0 || ne<=0) return 1;
    
    int nsd;
    parse_elem_type(etype, &nsd, NULL);
    
    
    MixdFile<double> mxyz("./mxyz", nn, nsd);
    mxyz.read();
    
    cout << endl;
    
    
    cout << "Computing new coordinates... " << flush;
    for(long i=0; i<nn; i++)
    {
        mxyz(i,0) += dx;
        mxyz(i,1) += dy;
        if(nsd>2) mxyz(i,2) += dz;
    }
    cout << "done!" << endl << endl;
    
    mxyz.write();
    
    } catch(MixdException e)
    { cout << e.msg() << endl; }    
    
    
    return 0;
}
