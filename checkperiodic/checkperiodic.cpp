#include "../mixdclass/mixd.hpp"
#include <getopt.h>
/* #include "SimpleProgress.hpp" */

int main(int argc, char * argv[])
{

    if (argc < 2)
    {
        std::cout << "checkperiodic: check if the solutions are truly periodic." <<std::endl;
        std::cout << "Run this program in the mesh directory." << std::endl;
        std::cout << "Usage: ./checkperiodic <data_file> -n/--ndf <ndf> -t/--nts <nts>" << std::endl;
        exit(-1);
    }

    long ne, nn;
    const int nsd = 3;  // 3D only!!!
    const int nen = 4;  // tetrahedra only!!!
    const int nef = 4;  // tetrahedra only!!!
    const int nnf = 3;  // number of nodes per face... tetra only!


    int c;
    int ndf = 0;
    int nts = 1;

    while (1)
    {
        static struct option long_options[] =
        {
            /* These options set a flag. */
            /* {"verbose", no_argument,       &verbose_flag, 1}, */
            /* {"brief",   no_argument,       &verbose_flag, 0}, */
            /* These options don’t set a flag.
               We distinguish them by their indices. */
            {"ndf",  required_argument, 0, 'n'},
            {"nts",  required_argument, 0, 't'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "n:t:",
                long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {
            case 0:
                /* If this option set a flag, do nothing else now. */
                if (long_options[option_index].flag != 0)
                    break;
                printf ("option %s", long_options[option_index].name);
                if (optarg)
                    printf (" with arg %s", optarg);
                printf ("\n");
                break;

            case 'n': ndf = std::atoi(optarg); break;
            case 't': nts = std::atoi(optarg); break;
            case '?':
                /* getopt_long already printed an error message. */
                break;
            default: abort ();
        }
    }

    if (optind < argc)
        std::cout << "data: " << argv[optind] << std::endl;
    else
    {
        std::cout << "No data file provided!" << std::endl;
        exit(-1);
    }

    std::cout << "ndf: " << ndf << std::endl;
    std::cout << "nts: " << nts << std::endl;

    mixd::readminf("minf", &nn, &ne);
    mixd::MixdFile<double> mxyz("mxyz", nn, nsd);
    mixd::MixdFile<int> mprd("mprd", nn);

    mxyz.read();
    mprd.read();

    std::cout << "Analyzing " << nts << " time steps..." << std::flush;
    for(int its=0; its<nts; its++)
    {
        /* std::cout << std::setw(5) << its << std::flush; */

        // only reads the first file; one file
        mixd::MixdFile<double> data(std::string(argv[optind]), nn, ndf, false);
        data.read(its);

        bool same=false;
        for (int in=0; in<nn; in++)
            if (mprd(in) != in + 1)
            {
                for (int idf=0; idf<ndf; idf++)
                {
                    same = (data(in,idf) == data(mprd(in)-1, idf));

                    /* std::cout << "(" << in+1 << ", " << mprd(in) << ")\r" << std::flush; */

                    if (!same)
                    {
                        std::cout << "ERROR: NOT SAME!" << std::endl;
                        exit(-1);
                    }

                }
            }

/* #pragma omp parallel for */
/*             for(int i=0; i<nSelBeads; i++) */
/*                 selBeads.at(i)->checkValues(mien, mxyz, data, its, 2); */

    }

    /* std::cout << std::endl; */
    std::cout << "done!"<< std::endl;
    std::cout << "Data is periodic."<< std::endl;

}

