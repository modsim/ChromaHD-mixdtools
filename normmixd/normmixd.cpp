#include "mixd.hpp"
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>
#include <unistd.h>
#include <vector>

void print_usage(char* binaryName)
{
    std::cout << "Norm MIXD: Compare two binary data files" << std::endl;
    std::cout << "Usage: " << binaryName << " <ref mixd data> <mixd data> -n <ndf>" << std::endl;
}

int main(int argc, char **argv)
{
    int option;
    int ndf = 0; // number of degrees of freedom == number of columns in mixd data
    int its = 0; // timestep to read
    std::string minffile = "../mesh/minf";

    bool spacetime = false;

    // Use getopt to parse command-line arguments
    while ((option = getopt(argc, argv, "t:m:n:s")) != -1) {
        switch (option) {
            case 't':
                its = std::stoi(optarg);
                break;
            case 'n':
                ndf = std::stoi(optarg);
                break;
            case 'm':
                minffile = optarg;
                break;
            case 's':
                spacetime=true;
                break;
            case '?':
                // Handling unknown options or missing arguments
                if (optopt == 'n') {
                    std::cerr << "Option -n requires an argument.\n";
                } else {
                    std::cerr << "Unknown option: -" << char(optopt) << "\n";
                }
                return 1;
            default:
                // Unexpected case
                std::cerr << "Unexpected case!" << std::endl;
                return 1;
        }
    }

    // Check if two file arguments are provided
    if (argc - optind != 2) {
        std::cerr << "Usage: " << argv[0] << " <file1> <file2> -n <int>\n";
        return 1;
    }

    // Access the file names after the options
    std::string file1 = argv[optind];
    std::string file2 = argv[optind + 1];

    long ne, nn, nnspace;
    mixd::readminf(minffile, &nn, &ne);

    if(spacetime)
    {
        if(nn % 2 != 0)
        {
            throw mixd::MixdException("space-time mesh must have even number of nodes!");
        }
        nnspace = nn/2;
    }
    else
    nnspace = nn;


    int stu_offset = 0;
    if (spacetime)
        stu_offset = nnspace;

    // Display parsed values
    std::cout << "mixd ref data file: " << file1 << std::endl;
    std::cout << "mixd cmp data file: " << file2 << std::endl;
    std::cout << "mixd minf file: " << minffile << std::endl;
    std::cout << "spacetime: " << spacetime << std::endl;
    std::cout << "ndf: " << ndf << std::endl;
    std::cout << "its: " << its << std::endl;
    std::cout << "nn: " << nn << std::endl;
    std::cout << "nnspace: " << nnspace << std::endl;
    std::cout << "stu_offset: " << stu_offset << std::endl;

    std::vector<double> vec_norms(ndf, 0.0);
    mixd::MixdFile<double> ref(file1, nn, ndf, false);
    mixd::MixdFile<double> cmp(file2, nn, ndf, false);

    ref.read(its);
    cmp.read(its);

    for(int i=stu_offset; i<nn; i++)
    {
        for(int j=0; j<ndf; j++) 
        {
            vec_norms[j] += pow((cmp(i,j) - ref(i,j)), 2);
        }
    }

    std::ofstream ofile("normmixd.csv");

    std::cout << "L2-norm, normalized-L2-norm, refmin, refmax, cmpmin, cmpmax" << std::endl;
    ofile << "L2-norm, normalized-L2-norm, refmin, refmax, cmpmin, cmpmax" << std::endl;

    for(int idf=0; idf<ndf; idf++)
    {
        vec_norms[idf] /= nnspace;
        vec_norms[idf] = sqrt(vec_norms[idf]);
        double relNorm = vec_norms[idf] / (ref.maxcol(idf) - ref.mincol(idf));
        std::cout << std::scientific;
        std::cout << vec_norms[idf] << ", " << relNorm << ", ";
        std::cout << ref.mincol(idf) << ", " << ref.maxcol(idf) << ", ";
        std::cout << cmp.mincol(idf) << ", " << cmp.maxcol(idf) << std::endl;
        ofile << vec_norms[idf] << ", " << relNorm << ", ";
        ofile << ref.mincol(idf) << ", " << ref.maxcol(idf) << ", ";
        ofile << cmp.mincol(idf) << ", " << cmp.maxcol(idf) << std::endl;
    }
    ofile.close();

    return 0;
}


