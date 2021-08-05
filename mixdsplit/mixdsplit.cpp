#include "mixd.hpp"
#include "SimpleProgress.hpp"
#include <cstdint>
#include <unistd.h>
#include <filesystem>

void printUsage(char *binaryName)
{
    std::cout << binaryName << " [ -m {spacetime mesh directory} ] [-d {data.all}] [-o {splitdata.all}]" << std::endl;
    std::cout << "MIXD tool to separate a chromatography mesh+data into the interstitial and particle domains" << std::endl;
    std::cout << "Run in directory with spacetime data.all. Expects to find ../mesh/{minf,nmap}" << std::endl;
    std::cout << "Automatically uses the top timeslab data. New data is semidiscrete." << std::endl;
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

    char c;


    string meshdir  = "../mesh";
    string datafile = "data.all";
    string splitdatafile = "splitdata.all";

    while ((c = getopt(argc, argv, "m:d:o:")) != -1)
    {
        switch(c)
        {
            case 'm':
                meshdir = optarg;
                break;
            case 'd':
                datafile = optarg;
                break;
            case 'o':
                splitdatafile = optarg;
                break;
            default:
                printUsage(argv[0]);
                return 1;
        }
    }

    string minffile = meshdir + "/minf";
    string nmapfile = meshdir + "/nmap";

    if(optind != argc)
    {
        printUsage(argv[0]);
        return 1;
    }

    std::cout << "meshdir: " << meshdir << std::endl;
    std::cout << "minffile: " << minffile << std::endl;
    std::cout << "nmapfile: " << nmapfile << std::endl;
    std::cout << "datafile: " << datafile << std::endl;
    std::cout << "splitdatafile: " << splitdatafile << std::endl;


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

        // TODO: calculate nts from filesize and nn
        int stu_offset = 0;

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
