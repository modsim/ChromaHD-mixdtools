#ifndef CONFIG_H_
#define CONFIG_H_

#include <vector>
#include <string>
#include <vtk/vtkCellType.h>

class Config
{
    private:
        int       argc;
        char**    argv;

        std::string    title;
        std::string    outpath;

        VTKCellType elemType;

        std::string                 minfFile;
        std::string                 mxyzFile;
        std::string                 mienFile;
        std::vector<std::string>    dataFiles;
        bool      spacetime;

        int       ndf;
        int       nrec;
        int       nrecstride;
        int       nrecoffset;

        double         dt;
        std::string    dtFile;

    protected:

    public:
        /// CONSTRUCTORS ///
        Config(int, char**);
        ~Config();

        /// GETTERS ///
        int                      getArgc()                {return argc;};
        char**                   getArgv()                {return argv;};

        std::string              getTitle()               {return title;};
        std::string              getOutpath()             {return outpath;};

        std::string              getMinfFile()            {return minfFile;};
        std::string              getMxyzFile()            {return mxyzFile;};
        std::string              getMienFile()            {return mienFile;};
        std::vector<std::string> getDataFiles()           {return dataFiles;};

        int                      getNdf()                 {return ndf;};
        int                      getNrec()                {return nrec;};
        int                      getNrecstride()          {return nrecstride;};
        int                      getNrecoffset()          {return nrecoffset;};
        int                      getSpacetime()           {return spacetime;};

        double                   getDt()                  {return dt;};
        std::string              getDtFile()              {return dtFile;};

        VTKCellType getElemType()            {return elemType;};

        //INTERFACE METHOD
        void read(std::string filename);
        void print();
        void readCommandlineArguments();

};

#endif /* CONFIG_H_ */
