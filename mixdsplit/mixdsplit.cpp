#include "mixd.hpp"
#include "SimpleProgress.hpp"
#include <cstdint>
#include <unistd.h>
#include <filesystem>

/* Use this program in conjunction with rmmat to split data from a full multi-domain simulation run.
 *
* rmmat -tet -st interstitial_mesh_dir 2 # Remove packed bed region
* mixdsplit -m <minffile> -N <nmapfile> -o interstitial_c.all
* rmmat -tet -st bed-mesh 1 # Remove interstitial region
* mixdsplit -m <minffile> -N <nmapfile> -o bed_c.all -i 0
* mixdsplit -m <minffile> -N <nmapfile> -o bed_q.all -i 1
*/

void printUsage(char *binaryName)
{
    std::cout << binaryName << " [-m {minffile}] [-N {nmapfile}] [-d {data.all}] [-o {splitdata.all}] [-i {idf}]" << std::endl;
    std::cout << "MIXD tool to separate a chromatography DATA into the interstitial and particle domains." << std::endl;
    std::cout << "Run in directory with spacetime data.all. Expects to find ../mesh/{minf,nmap} by default." << std::endl;
    std::cout << "Automatically uses the top timeslab data. New data is semidiscrete." << std::endl;
    std::cout << "DOES NOTE GENERATE A MIXD MESH! ONLY SPLITS DATA FOR THE NODES MAPPED IN NMAP!" << std::endl;
}


/*
 * @author : Jayghosh S. Rao
 * @date   : 15 Jul. 2021
 * @program: mixdsplit
 * @brief  :
 *           - The node relationship between FLOW and MASS meshes are given in the nmap file generated by rmmat.
 *           - This program maps the full spacetime 2-domain 2-ndf data.all to an interstitial 1-ndf data.all
 *           - TODO: split the mesh fully using mmat, split out the particle data as well
 */


int main(int argc, char **argv)
{
    using namespace std;
    using namespace mixd;

    /* printUsage(argv[0]); */

    bool spacetime = true;
    bool spacetimeupper = true;

    int ndf = 2;
    int idf = 0; // NOTE: Only extracting the first dof

    string datafile = "data.all";
    string splitdatafile = "splitdata.all";
    string minffile = "../mesh/minf";
    string nmapfile = "../mesh/nmap";

    char c;
    while ((c = getopt(argc, argv, "N:m:d:o:n:i:")) != -1)
    {
        switch(c)
        {
            case 'm':
                minffile = optarg;
                break;
            case 'd':
                datafile = optarg;
                break;
            case 'o':
                splitdatafile = optarg;
                break;
            case 'n':
                ndf = std::atoi(optarg);
                break;
            case 'N':
                nmapfile = optarg;
                break;
            case 'i':
                idf = std::atoi(optarg);
                break;
            default:
                printUsage(argv[0]);
                return 1;
        }
    }


    if(optind != argc)
    {
        printUsage(argv[0]);
        return 1;
    }

    std::cout << "minffile: " << minffile << std::endl;
    std::cout << "nmapfile: " << nmapfile << std::endl;
    std::cout << "datafile: " << datafile << std::endl;
    std::cout << "splitdatafile: " << splitdatafile << std::endl;
    std::cout << "ndf: " << ndf << std::endl;
    std::cout << "idf: " << idf << std::endl;

    std::cout << "NOTE: This program only extracts the data for ONE DOF specified with -i, default i = 0.\n" << std::endl;

    try{

        long ne, nn, nnspace;
        readminf(minffile, &nn, &ne);

        if(spacetime)
        {
            if(nn % 2 != 0)
                throw MixdException("space-time mesh must have even number of nodes!");

            nnspace = nn/2;
        }
        else
            nnspace = nn;

        std::cout << "nn: " << nn << std::endl;
        std::cout << "nnspace: " << nnspace << std::endl;

        // Read only half the file, because it just repeats
        MixdFile<int> nmap(nmapfile, nn/2, 1, false);
        nmap.read();

        int new_nn = nmap.max();
        std::cout << "nn split: " << new_nn << std::endl;

        mixd::MixdFile<double> data(datafile, nn, ndf, false);
        std::remove(splitdatafile.c_str());
        mixd::MixdFile<double> newdata(splitdatafile, new_nn, 1, false);

        //spacetime upper offset
        int stu_offset = 0;

        //calculate nts from filesize and nn
        std::uintmax_t filesize = std::filesystem::file_size(datafile);
        int nts = filesize / 8 / nn / ndf;

        std::cout << "Found " << nts << " timesteps!" << std::endl;

        if (spacetime && spacetimeupper)
            stu_offset = nnspace;


        SimpleProgress sp(0, nts, 20);
        for(int its=0; its<nts; its++)
        {
            data.read(its);
            for (long i=0; i<nnspace; i++)
            {
                if (nmap(i) > 0)
                {
                    newdata(nmap(i)-1, 0) = data(i + stu_offset, idf);
                }
            }
            newdata.append();
            sp.printIfHitNext(its);
        }

        std::cout << "done!" << std::endl;

    } catch(mixd::MixdException e)
    { std::cout << e.msg() << std::endl; }

    return 0;
}
