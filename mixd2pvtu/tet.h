#ifndef TET_H_
#define TET_H_

#include "settings.h"

/***************************************************************************************************
NODE LEVEL DATA STRUCTURE
****************************************************************************************************
Each node has
    - coordinates
    - a variable
Comments:
    - For the current application, the only variable is temperature since we are solving transient
      diffusion equation.
***************************************************************************************************/
class tetNode
{
    private:
        // VARIABLES
        double x;    // x-coordinate
        double y;    // y-coordinate
        double z;    // z-coordinate
        int ndf;
        /* double * data; */
    protected:

    public:
        // DEFAULT CONSTRUCTOR
        tetNode(){x=0.0f; y=0.0f; z=0.0f;};
        ~tetNode(){}

        // SETTERS
        inline void setX         (double value)    {x = value;};
        inline void setY         (double value)    {y = value;};
        inline void setZ         (double value)    {z = value;};
        /* inline void setData      (int index, double value) {data[index] = value;}; */

        // GETTERS
        inline double    getX()       {return x;};
        inline double    getY()       {return y;};
        inline double    getZ()       {return z;};
        inline double    getCoord(int i) { if (i==0) return x;
                                            else if (i==1) return y;
                                            else return z; };
        /* inline double    getData(int index)  {return data[index];}; */
};

/***************************************************************************************************
ELEMENT LEVEL DATA STRUCTURE
****************************************************************************************************
An element has
    - connectivity [nen]
    - boundary condition [nef]
    - Jacobian determinant [nGQP]
    - dsdX [3xnGQP]
    - dSdY [3xnGQP]
***************************************************************************************************/
class tetElement
{
    private:
        // VARIABLES
        int conn[nen];
        int lConn[nen];
    protected:

    public:
        // DEFAULT CONSTRUCTOR
        // I must aboulutly write a proper constructor here :)
        // SETTERS
        void setConn    (int i, int value)              {conn[i] = value;};
        void setLConn   (int i, int value)              {lConn[i] = value;};
        // GETTERS
        int           getConn (int index)       {return conn[index];};
        int           getLConn(int index)       {return lConn[index];};
};

/***************************************************************************************************
MESH DATA STRUCTURE
****************************************************************************************************
This class is used to keep the pointers to node and element data structures. It is more
convinient to pass the pointer to this class rather than passin the pointers to both element and
node arrays during function calls. So it is just for simplification.
***************************************************************************************************/
class tetMesh
{
    private:
        // VARIABLES
        int ne;                     // total number of elements
        int nec;                    // number of elements per cpu
        int mec;                    // max. number of elements among all cpus
        int nn;                     // number of nodes
        int nnc;                    // number of nodes per core
        int nnl;                    // number of local nodes
        int mnc;                    // max. number of nodes among all cpus
        int ndf;
        int * nodeLToG;             // array for local to global connectivity conversion
        tetNode*            node;   // pointer to partition, node level data structure
        tetNode*            lNode;  // pointer to local, node level data structure
        tetElement*         elem;   // pointer to element level data structure

        double * xyz;
        double * dataL; // size nnl * ndf
        double * dataG; // size nnc * ndf

        //METHODS
        void readMeshFiles(inputSettings*);
        void readDataFile(inputSettings*, int irec);
        void swapBytes(char*, int, int);    // Used during file read operations
        void formLocalNodeList();
        void quickSort(int*, int*, int, int);
        void localizeNodeCoordinates();
        void localizeData();
        int isSemiDiscrete(const char* filename, int nn);
    protected:

    public:
        // DEFAULT CONSTRUCTOR
        tetMesh(){};
        // DESTRUCTOR
        ~tetMesh()
        {
            delete[] node;    // Recovers memory for node level data structure.
            delete[] lNode;   // Recovers memory for local node level data structure.
            delete[] elem;    // Recovers memory for element level data structure.
        };
        // GETTERS
        int        getNe()                 {return ne;};
        int        getNec()                {return nec;};
        int        getMec()                {return mec;};
        int        getNn()                 {return nn;};
        int        getNnc()                {return nnc;};
        int        getNnl()                {return nnl;};
        int        getMnc()                {return mnc;};
        int        getNdf()                {return ndf;};
        int*       getNodeLToG()           {return nodeLToG;};

        double*     getXyz()                {return xyz;};
        double*     getDataG()              {return dataG;};
        double*     getDataL()              {return dataL;};

        tetNode*           getNode     (int index)    {return &node[index];};
        tetNode*           getLNode    (int index)    {return &lNode[index];};
        tetElement*        getElem     (int index)    {return &elem[index];};

        // INTERFACE METHOD
        void prepareMesh(inputSettings*);
        void processData(inputSettings* settings, int irec);
};

#endif /* TET_H_ */
