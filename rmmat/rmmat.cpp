#include "mixd.hpp"


void printUsage(char *binaryName)
{
    std::cout << "MIXD tool to remove material from mesh" << std::endl;
    std::cout << "Usage: " << binaryName << " [-tri|-tet|-qua|-hex] [-sd|-st] <output directory> <mat1 to remove> [<mat2 to remove>]" << std::endl;
    std::cout << "The -sd/-st parameter specifies whether the initial mesh is semi-discrete or space-time." << std::endl;
    std::cout << "The newly generated mesh (with material removed) will always be semi-discrete!" << std::endl;
}


int main(int argc, char **argv)
{
    using namespace std;
    using namespace mixd;

    if(argc != 5 && argc != 6)
    {
        printUsage(argv[0]);
        return 1;
    }

    string etype(argv[1]);
    string sdst(argv[2]);

    string outdir(argv[3]);

    const int rmmat1 = atoi(argv[4]);

    int rmmat2 = -999;
    if(argc == 6)
        rmmat2 = atoi(argv[5]);

    int nen = 0;
    int nef = 0;
    int nsd = 0;

    if(etype == "-tri")
    {
        nen = 3;
        nef = 3;
        nsd = 2;
    }
    else if(etype == "-tet")
    {
        nen = 4;
        nef = 4;
        nsd = 3;
    }
    else if(etype == "-qua")
    {
        nen = 4;
        nef = 4;
        nsd = 2;
    }
    else if(etype == "-hex")
    {
        nen = 8;
        nef = 6;
        nsd = 3;
    }
    else if (etype == "-tetP2")
    {
        nen = 10;
        nef = 4;
        nsd = 3;
    }
    else
    {
        printUsage(argv[0]);
        return 1;
    }

    bool spacetime = false;
    if(sdst == "-sd")
        spacetime = false;
    else if(sdst == "-st")
        spacetime = true;
    else
    {
        printUsage(argv[0]);
        return 1;
    }


    cout << "MIXD rmmat tool..." << endl << endl;


    try{

    long ne, nn, nnspace;
    readminf("./minf", &nn, &ne);

    if(spacetime)
    {
        if(nn % 2 != 0)
            throw MixdException("space-time mesh must have even number of nodes!");

        nnspace = nn/2;
    }
    else
        nnspace = nn;


    MixdFile<int> mien("./mien", ne, nen);
    mien.read();

    MixdFile<int> mmat("./mmat", ne);
    mmat.read();

    MixdFile<int> nmap("./nmap", nn);
    nmap.init();


    // count number of retained elements and nodes
    long ne_new = 0;
    long nn_new = 0;
    for(long i=0; i<ne; i++)
    {
        if(mmat(i) != rmmat1 && mmat(i) != rmmat2)
        {
            ne_new++;

            for(int j=0; j<nen; j++)
            {
                if(nmap(mien(i,j)-1) == 0)
                {
                    nn_new++;
                    nmap(mien(i,j)-1) = nn_new;
                    if(spacetime)
                        nmap(mien(i,j)-1 + nnspace) = nn_new;
                }
            }
        }
    }

    nmap.write();


    MixdFile<int> mrng("./mrng", ne, nef);
    mrng.read();

    MixdFile<double> mxyz("./mxyz", nn, nsd);
    mxyz.read();

    MixdFile<int> mien_new(outdir+"/mien", ne_new, nen);
    MixdFile<int> mmat_new(outdir+"/mmat", ne_new);
    MixdFile<int> mrng_new(outdir+"/mrng", ne_new, nef);

    long eid = 0;
    for(long i=0; i<ne; i++)
    {
        if(mmat(i) != rmmat1 && mmat(i) != rmmat2)
        {
            for(int j=0; j<nen; j++)
                mien_new(eid,j) = nmap(mien(i,j)-1);

            for(int j=0; j<nef; j++)
                mrng_new(eid,j) = mrng(i,j);

            mmat_new(eid) = mmat(i);

            eid++;
        }
    }

    mien_new.write();
    mrng_new.write();
    mmat_new.write();


    MixdFile<double> mxyz_new(outdir+"/mxyz", nn_new, nsd);

    for(long i=0; i<nnspace; i++)
    {
        if(nmap(i) > 0)
        {
            for(int j=0; j<nsd; j++)
                mxyz_new(nmap(i)-1, j) = mxyz(i,j);
        }
    }

    mxyz_new.write();



    writeminf(outdir+"/minf", nn_new, ne_new);





    } catch(mixd::MixdException e)
    { std::cout << e.msg() << std::endl; }



    return 0;
}
