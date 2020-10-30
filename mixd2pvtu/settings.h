#ifndef SETTINGS_H_
#define SETTINGS_H_

#include "constants.h"
#include <vector>

/***************************************************************************************************
Input parameters are stored in this class. readSettingsFile() method does the job: It opens the
settings.in file, reads the parameters and assigns them to the private variables.
Why triMEsh is a friend? : triMesh method readMeshFiles needs to set ne and nn. no other funtions
need to change a private variable of inputSettings. Thus, previliges are given to triMesh.
***************************************************************************************************/

class inputSettings
{
    private:
        // VARIABLES

        int       argc;       // command line argument
        char**    argv;        // command line argument
        string    title;       // title of the document
        string    minfFile;    // information file name
        string    mxyzFile;    // node coordinates file name
        string    mienFile;    // connctivity file name
        string    mrngFile;    // boundary info file name
        vector<string>    dataFiles;
        int       ndf; // number of scalar values
        int       nrec;
        int       nrecstride;  // get every nrecstride'th timestep record.
        int       nrecoffset;  // initial offset
        int       spacetime;    // if mesh is spacetime
        double    dt;          // time step size
        string    outpath;     // output directory

        // METHODS
        void readSettingsFile();
        void printSettings();

    protected:

    public:
        /// CONSTRUCTORS ///
        inputSettings();
        inputSettings(int, char**);
        /// GETTERS ///

        int            getArgc()                {return argc;};
        char**         getArgv()                {return argv;};
        string         getTitle()               {return title;};
        string         getMinfFile()            {return minfFile;};
        string         getMxyzFile()            {return mxyzFile;};
        string         getMienFile()            {return mienFile;};
        string         getMrngFile()            {return mrngFile;};
        vector<string> getDataFiles()           {return dataFiles;};
        int            getNdf()                 {return ndf;};
        int            getNrec()                {return nrec;};
        int            getNrecstride()          {return nrecstride;};
        int            getNrecoffset()          {return nrecoffset;};
        int            getSpacetime()           {return spacetime;};
        double         getDt()                  {return dt;};
        string         getOutpath()             {return outpath;};

        //INTERFACE METHOD
        void prepareSettings();

};

#endif /* SETTINGS_H_ */
