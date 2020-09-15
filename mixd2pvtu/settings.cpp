#include "settings.h"

/***************************************************************************************************
void inputSettings::inputSettings()
private - Default constructor
***************************************************************************************************/
inputSettings::inputSettings(int largc, char** largv)
{
    // Default values for the input parameters

    argc       = largc;                //Command line parameters are stored at private var
    argv       = largv;                //Command line parameters are stored at private var
    title      = "Default Title";
    minfFile   = "minf";
    mxyzFile   = "mxyz";
    mienFile   = "mien";
    mrngFile   = "mrng";
    dataFile   = "data";
    nrec       = 1;
    nrecstride = 1;
    nrecoffset = 0;
    ndf        = 1;
    dt         = 1.0;
    spacetime  = 0;
}

/***************************************************************************************************
void inputSettings::readSettingsFile()
public - main interface function, calls all other, private, functions of the class
***************************************************************************************************/
void inputSettings::prepareSettings()
{
    readSettingsFile();
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
                iss >> dataFile;
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
        cout << "Title of the simualation              : " << title		<< endl;
        cout << "Name of the minf file                 : " << minfFile	<< endl;
        cout << "Name of the mxyz file                 : " << mxyzFile	<< endl;
        cout << "Name of the mien file                 : " << mienFile	<< endl;
        cout << "Name of the mrng file                 : " << mrngFile	<< endl;
        cout << "Name of the initial distribution file : " << dataFile	<< endl;
        cout << "Number of variables in data           : " << ndf		<< endl;
        cout << "Number of Timesteps                   : " << nrec		<< endl;
        cout << "Number of Timesteps Stride            : " << nrecstride<< endl;
        cout << "Number of Timesteps Offset            : " << nrecoffset<< endl;
        cout << "Spacetime                             : " << spacetime << endl;
        /* cout << "Time step size                        : " << dt		<< endl; */
    }

    return;
}
