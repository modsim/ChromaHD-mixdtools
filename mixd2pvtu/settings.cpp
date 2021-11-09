#include "settings.h"
#include <getopt.h>

/***************************************************************************************************
void inputSettings::inputSettings()
private - Default constructor
***************************************************************************************************/
inputSettings::inputSettings(int largc, char** largv)
{
    // Default values for the input parameters

    argc       = largc;                //Command line parameters are stored at private var
    argv       = largv;                //Command line parameters are stored at private var
    title      = "vis";
    minfFile   = "minf";
    mxyzFile   = "mxyz";
    mienFile   = "mien";
    mrngFile   = "mrng";
    nrec       = 0;
    nrecstride = 1;
    nrecoffset = 0;
    ndf        = 1;
    dt         = 1.0;
    spacetime  = 0;
    outpath    = "output";
    dtFile     = "";
}

/***************************************************************************************************
void inputSettings::readSettingsFile()
public - main interface function, calls all other, private, functions of the class
***************************************************************************************************/
void inputSettings::prepareSettings()
{
    readSettingsFile();
    readCommandlineArguments();
    printSettings();

    return;
}


/***************************************************************************************************
void inputSettings::readSettingsFile()
private -reads settings.in file
***************************************************************************************************/
void inputSettings::readSettingsFile()
{
    int mype, npes;                // my processor rank and total number of processors
    MPI_Comm_rank(MPI_COMM_WORLD, &mype);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);
    char * fileName;
    string lineString;
    string dummyString;
    string dummyString2;
    char   dummyChar;

    ifstream inputFile;
    if (argc > 1)
    {
        fileName = argv[1];
    }
    else
    {
        fileName = (char*)"settings.in";
    }

    inputFile.open(fileName,ios::in);
    if (inputFile.is_open()==false)
    {
        cout << "Unable to open input file for pe: " << mype << "! Aborting... " << endl;
        MPI_Finalize();
        exit(0);
    }

    while (!inputFile.eof())
    {
        // Get a line and store in lineString
        getline(inputFile, lineString, '\n');
        // If the first character of the line is not a '#'
        if (lineString.c_str()[0] != '#')
        {
            istringstream iss(lineString);
            iss >> dummyString;
            if(dummyString == "title")
                iss >> title;
            else if(dummyString == "minf")
                iss >> minfFile;
            else if(dummyString == "mxyz")
                iss >> mxyzFile;
            else if(dummyString == "mien")
                iss >> mienFile;
            else if(dummyString == "mrng")
                iss >> mrngFile;
            else if(dummyString == "data")
            {
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
                iss >> spacetime;
            else if(dummyString == "dt")
                iss >> dt;
            else if(dummyString == "dtFile")
                iss >> dtFile;
            else if(dummyString == "outpath")
                iss >> outpath;
            else
            {
                cout << endl << "Unknown keyword in the settings file : " << dummyString;
                cout << endl << "Aborting...";
                inputFile.close();
                MPI_Finalize();
                exit(0);
            }
        }
    }
    inputFile.close();

    return;
}


/***************************************************************************************************
void inputSettings::printSettings()
***************************************************************************************************/
void inputSettings::printSettings()
{
    int mype, npes;                // my processor rank and total number of processors
    MPI_Comm_rank(MPI_COMM_WORLD, &mype);
    MPI_Comm_size(MPI_COMM_WORLD, &npes);

//    cout << fixed;
    cout.precision(4);
    cout << scientific;
    if (mype == 0)
    {
        cout << endl << "==================== SETTINGS ====================" << endl;

        cout << "Title of the simualation              : " << title		      << endl;
        cout << "Spacetime                             : " << spacetime       << endl;
        cout << "Name of the minf file                 : " << minfFile	      << endl;
        cout << "Name of the mxyz file                 : " << mxyzFile	      << endl;
        cout << "Name of the mien file                 : " << mienFile	      << endl;
        cout << "Name of the mrng file                 : " << mrngFile	      << endl;
        cout << "Number of variables in data           : " << ndf		      << endl;
        cout << "Number of Timesteps                   : " << nrec		      << endl;
        cout << "Number of Timesteps Stride            : " << nrecstride      << endl;
        cout << "Number of Timesteps Offset            : " << nrecoffset      << endl;
        cout << "Number of Data Files                  : " << dataFiles.size()<< endl;
        cout << "Output Directory                      : " << outpath         << endl;
        cout << "Timesteps File                        : " << dtFile          << endl;
        cout << "Data Files                            : " << flush;

        for(vector<string>::iterator it = dataFiles.begin(); it != dataFiles.end(); ++it)
            cout << *it << " " << flush;
        cout << endl;
        /* cout << "Time step size                        : " << dt		<< endl; */
    }

    return;
}

void inputSettings::readCommandlineArguments()
{
    if (argc < 2)
    {
        std::cout << "mixd2pvtu: convert from mixd to pvtu." <<std::endl;
        std::cout << "Usage: ./mixd2pvtu <settings.in file> <override arguments>" << std::endl;
        exit(-1);
    }

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
            {"ndf",  required_argument, 0, 'n'},
            {"nrec",  required_argument, 0, 'r'},
            {"spacetime",  required_argument, 0, 's'},
            {"title",  required_argument, 0, 't'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "n:r:st:", long_options, &option_index);

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
            case 's': spacetime = 1; break;
            case 't': title = optarg; break;
            case '?':
                /* getopt_long already printed an error message. */
                break;
            default: abort ();
        }
    }

    if (optind >= argc)
    {
        std::cout << "No settings file provided!" << std::endl;
        exit(-1);
    }

}
