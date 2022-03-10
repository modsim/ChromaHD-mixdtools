#include "config.hpp"
#include "utils.hpp"

#include <getopt.h>
#include <mpi.h>

#include <iostream>
#include <fstream>
#include <sstream>

Config::Config(int largc, char** largv)
{
    // Default values for the input parameters

    argc       = largc;                //Command line parameters are stored at private var
    argv       = largv;                //Command line parameters are stored at private var

    title      = "vis";
    outpath    = "output";

    elemType   = VTK_TETRA;
    minfFile   = "minf";
    mxyzFile   = "mxyz";
    mienFile   = "mien";
    dataFiles  = { "data.all" };

    spacetime  = false;
    ndf        = 1;
    nrec       = 0;
    nrecoffset = 0;
    nrecstride = 1;

    dt         = 1.0;
    dtFile     = "";

    if (argc < 2)
    {
        std::cout << "mixd2pvtu: convert from mixd to pvtu." <<std::endl;
        std::cout << "Usage: ./mixd2pvtu <settings.in file> <override arguments>" << std::endl;
        exit(-1);
    }
    else if ((argc >= 2) && (argv[1][0] != '-'))
    {
        read(argv[1]);
    }

    readCommandlineArguments();
    print();

}

Config::~Config(){};

// NOTE: Every process reads the input file at the same time
void Config::read(std::string filename)
{
    int mype, npes;                // my processor rank and total number of processors

    MPI_Comm_rank(MPI_COMM_WORLD, &mype);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);

    char        dummyChar;
    std::string lineString;
    std::string dummyString;
    std::string dummyString2;
    std::ifstream inputFile;

    inputFile.open(filename,std::ios::in);
    if (inputFile.is_open()==false)
    {
        std::cout << "Unable to open input file for pe: " << mype << "! Using defaults and commandline arguments... " << std::endl;
        return;
    }

    while (!inputFile.eof())
    {
        // Get a line and store in lineString
        getline(inputFile, lineString, '\n');
        // If the first character of the line is not a '#'
        if ((lineString.c_str()[0] != '#') && (!lineString.empty()))
        {
            std::istringstream iss(lineString);
            iss >> dummyString;
            if(dummyString == "title")
                iss >> title;
            else if(dummyString == "minf")
                iss >> minfFile;
            else if(dummyString == "mxyz")
                iss >> mxyzFile;
            else if(dummyString == "mien")
                iss >> mienFile;
            else if(dummyString == "elemtype")
            {
                getline(iss, dummyString2, ' ');
                getline(iss, dummyString2, ' ');
                elemType = processElementType(dummyString2);
            }
            else if(dummyString == "data")
            {
                dataFiles.clear();
                while (getline(iss, dummyString2, ' '))
                {
                    if (dummyString2 != " " && (!dummyString2.empty()))
                        dataFiles.push_back(dummyString2);
                }
            }
            else if(dummyString == "ndf")
                iss >> ndf;
            else if(dummyString == "nrec")
                iss >> nrec;
            else if(dummyString == "nrecstride")
                iss >> nrecstride;
            else if(dummyString == "nrecoffset")
                iss >> nrecoffset;
            else if(dummyString == "spacetime")
                spacetime = true;
            else if(dummyString == "dt")
                iss >> dt;
            else if(dummyString == "dtFile")
                iss >> dtFile;
            else if(dummyString == "outpath")
                iss >> outpath;
            else
            {
                std::cout << std::endl << "Unknown keyword in the settings file : " << dummyString;
                std::cout << std::endl << "Aborting...";
                inputFile.close();
                MPI_Finalize();
                exit(0);
            }
        }
    }
    inputFile.close();

    return;
}


void Config::print()
{
    int mype, npes;                // my processor rank and total number of processors
    MPI_Comm_rank(MPI_COMM_WORLD, &mype);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);

//    cout << fixed;
    std::cout.precision(4);
    std::cout << std::scientific;
    if (mype == 0)
    {
        std::cout << std::endl << "==================== SETTINGS ====================" << std::endl;

        std::cout << "Title of the simualation              : " << title   		      << std::endl;
        std::cout << "Output Directory                      : " << outpath            << std::endl;
        std::cout << "Spacetime                             : " << spacetime          << std::endl;
        std::cout << "Element Type                          : " << elemType           << std::endl;
        std::cout << "Name of the minf file                 : " << minfFile	          << std::endl;
        std::cout << "Name of the mxyz file                 : " << mxyzFile	          << std::endl;
        std::cout << "Name of the mien file                 : " << mienFile	          << std::endl;
        std::cout << "Number of variables in data           : " << ndf     		      << std::endl;
        std::cout << "Number of Timesteps                   : " << nrec    		      << std::endl;
        std::cout << "Number of Timesteps Stride            : " << nrecstride         << std::endl;
        std::cout << "Number of Timesteps Offset            : " << nrecoffset         << std::endl;
        std::cout << "Number of Data Files                  : " << dataFiles.size()   << std::endl;
        std::cout << "Data Files                            : " << std::flush;

        for (auto it = dataFiles.begin(); it != dataFiles.end(); it++)
            std::cout << *it << " " << std::flush;
        std::cout << std::endl;
        std::cout << "Time step size                        : " << dt		<< std::endl;
        std::cout << "Timesteps File                        : " << dtFile             << std::endl;
    }

    return;
}

void Config::readCommandlineArguments()
{

    int c;

    while (1)
    {
        static struct option long_options[] =
        {
            /* These options set a flag. */
            /* {"verbose", no_argument,       &verbose_flag, 1}, */
            /* {"brief",   no_argument,       &verbose_flag, 0}, */
            /* {"spacetime-upper",  no_argument, &spacetimeupper, 1}, */
            /* {"s",  no_argument, &spacetimeupper, 1}, */
            /* These options don’t set a flag.
               We distinguish them by their indices. */
            {"spacetime",  no_argument, 0, 's'},
            {"ndf",  required_argument, 0, 'n'},
            {"nrec",  required_argument, 0, 'r'},
            {"title",  required_argument, 0, 'T'},
            {"minf", required_argument, 0, 'i'},
            {"mxyz", required_argument, 0, 'x'},
            {"mien", required_argument, 0, 'e'},
            {"outpath", required_argument, 0, 'o'},
            {"elemtype", required_argument, 0, 'E'},
            {"data", required_argument, 0, 'd'},
            {"timesteps",  required_argument, 0, 't'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "n:r:t:i:x:e:o:E:s", long_options, &option_index);

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
            case 'r': nrec = std::atoi(optarg); break;
            case 't': dt = std::atoi(optarg); break;
            case 'T': title = optarg; break;
            case 'o': outpath = optarg; break;
            case 'i': minfFile = optarg; break;
            case 'x': mxyzFile = optarg; break;
            case 'e': mienFile = optarg; break;
            case 'd': dataFiles.push_back(optarg); break;
            case 'E': elemType = processElementType(optarg); break;
            case 's': spacetime = true; break;
            case '?':
                /* getopt_long already printed an error message. */
                break;
            default: abort ();
        }

    }

    // TODO: Check this
    if (optind >= argc)
    {
        std::cout << "No settings file provided!" << std::endl;
    }
}
