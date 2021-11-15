#ifndef MESH_H_
#define MESH_H_

#include "config.hpp"
#include "node.hpp"
#include "element.hpp"

#include <mpi.h>

class Mesh {

    private: 

    protected:

    public:
    Mesh(Config config);

    ~Mesh()
    {
        /* delete[] node;    // Recovers memory for node level data structure. */
        delete[] lNode;   // Recovers memory for local node level data structure.
        // delete[] elem;    // Recovers memory for element level data structure.
        delete[] dataL;
        delete[] dataG;
        delete[] nodeLToG;
    };


    long nn;
    long ne;
    long nnspace;

    long nec;
    long mec;
    long nnc;
    long nnl;
    long mnc;

    int ndf; // from config

    int nsd = 3;
    int nen;

    double * xyz;
    double * dataL; // size nnl * ndf
    double * dataG; // size nnc * ndf

    //TODO: Consider making this long
    int * nodeLToG; // array for local to global connectivity conversion 

    Node *node;
    Node *lNode;

    //NOTE: These bits are from the config class
    VTKCellType elemType;
    std::vector<std::string> dataFiles;
    int nrec;
    int nrecoffset;
    int nrecstride;
    bool spacetime;
    std::string title;
    std::string outpath;

    Element *elem;
    // std::vector<Element> elem;

    double* getXyz()   {return xyz;};
    double* getDataG() {return dataG;};
    double* getDataL() {return dataL;};

    Node*    getNode (int  index){return &node[index];};
    Node*    getLNode(int index) {return &lNode[index];};
    Element* getElem (int  index){return &elem[index];};


    std::vector<double> timesteps;

    void distribute();
    void readminf(std::string filename, bool spacetime);
    void readmxyz(std::string filename);
    void readmien(std::string filename);
    void getTimesteps(std::string dtFile, int dt, int nrec);

    void formLocalNodeList();
    void localizeNodeCoordinates();
    void localizeData();

    void readDataRecord(
        int irec,
        int nrecstride,
        int nrecoffset,
        bool spacetime,
        std::vector<std::string> dataFiles
    );

    void getFileAndOffset(
        int irec, int nrecstride, int nrecoffset, int ndf,
        std::vector<std::string> dataFiles,
        std::string& dataFile,
        MPI_Offset& totalOffset
    );

    void write();
    void vtkVisualization(int irec, std::string title, std::string outpath, int VTKElemType);


};
#endif // !MESH_H_
