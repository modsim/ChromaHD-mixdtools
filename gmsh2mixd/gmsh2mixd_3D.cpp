#include <iostream>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <set>
#include <map>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <unistd.h>
#include "SimpleProgress.hpp"






bool isBigEndian()
{
    short word = 0x4321;
    if ((* (char*) & word) != 0x21) return true;
    else      return false;
}

void swapbytes(char *array, int nelem, int elsize)
{
    register int sizet, sizem, i, j;
    char *bytea, *byteb;
    sizet = elsize;
    sizem = sizet - 1;
    bytea = (char*)malloc(sizet);
    byteb = (char*)malloc(sizet);

    for (i=0; i<nelem; i++)
    {
        memcpy((void *)bytea, (void *)(array+i*sizet), sizet);
        for (j = 0; j < sizet; j++) byteb[j] = bytea[sizem - j];
            memcpy((void *)(array+i*sizet), (void *)byteb, sizet);
    }

    free(bytea);
    free(byteb);
}



// Unique triangle keys generated from their node IDs
class TriKey
{
public:
    int nodeIDs[3];

    TriKey(const int *nids)
    {
        if(nids[0]<=nids[1] && nids[1]<=nids[2]){ nodeIDs[0]=nids[0]; nodeIDs[1]=nids[1]; nodeIDs[2]=nids[2]; }
        if(nids[0]<=nids[2] && nids[2]<=nids[1]){ nodeIDs[0]=nids[0]; nodeIDs[1]=nids[2]; nodeIDs[2]=nids[1]; }
        if(nids[1]<=nids[0] && nids[0]<=nids[2]){ nodeIDs[0]=nids[1]; nodeIDs[1]=nids[0]; nodeIDs[2]=nids[2]; }
        if(nids[1]<=nids[2] && nids[2]<=nids[0]){ nodeIDs[0]=nids[1]; nodeIDs[1]=nids[2]; nodeIDs[2]=nids[0]; }
        if(nids[2]<=nids[0] && nids[0]<=nids[1]){ nodeIDs[0]=nids[2]; nodeIDs[1]=nids[0]; nodeIDs[2]=nids[1]; }
        if(nids[2]<=nids[1] && nids[1]<=nids[0]){ nodeIDs[0]=nids[2]; nodeIDs[1]=nids[1]; nodeIDs[2]=nids[0]; }
    }
};


// functor class for comparison of TriKeys
class compareTriKey
{
public:
    bool operator()(const TriKey & key1, const TriKey & key2) const
    {
        if(key1.nodeIDs[0] < key2.nodeIDs[0]) return true;
        if(key1.nodeIDs[0] > key2.nodeIDs[0]) return false;

        if(key1.nodeIDs[1] < key2.nodeIDs[1]) return true;
        if(key1.nodeIDs[1] > key2.nodeIDs[1]) return false;

        if(key1.nodeIDs[2] < key2.nodeIDs[2]) return true;
        if(key1.nodeIDs[2] > key2.nodeIDs[2]) return false;

        return false;
    }
};



// a triangle element consists of an ID, its three nodes IDs, and a related boundary index
class Triangle
{
public:
    int elementID;
    int boundaryID;

    int nodes[3];

    int adjacentTetras[2];   // IDs of the two tetras between which this triangle is situated
    int nAdjoiningTetras;    // number of adjacent tetras

    Triangle(int eid_, int bid_, const int *nodes_) : elementID(eid_), boundaryID(bid_)
    {
        for(int i=0; i<3; i++)
            nodes[i] = nodes_[i];

        nAdjoiningTetras  = 0;
        adjacentTetras[0] = 0;
        adjacentTetras[1] = 0;
    }


    inline void setAdjacentTetra(int tetid)
    {
        if(nAdjoiningTetras>1)
        {
            std::cout << "FATAL ERROR in setAdjacentTetra!" << std::endl;

            std::cout << "Triangle nodes: " << nodes[0] << "    " << nodes[1] << "    " << nodes[2] << std::endl;
            std::cout << " boundary id: " << boundaryID << std::endl;

            std::cout << "adjacentTetras: " << adjacentTetras[0] << "     " << adjacentTetras[1] << std::endl;



            exit(1);
        }

        adjacentTetras[nAdjoiningTetras] = tetid;
        nAdjoiningTetras++;
    }


    inline int getNeighborTetra(int tetid) const
    {
        if(adjacentTetras[0] == tetid)
            return adjacentTetras[1];

        if(adjacentTetras[1] == tetid)
            return adjacentTetras[0];

        std::cout << "FATAL ERROR in getNeighborTetra!" << std::endl;
        exit(1);
    }

};



// a tetrahedron element consists of an ID, its four node IDs, and boundary indices related to its four faces
class Tetrahedron
{
public:
    int elementID;
    int boundaryID[4];
    int materialID;

    // face ordering (according to MIXD format):
    // face 0 consists of nodes 0 1 2
    // face 1 consists of nodes 0 1 3
    // face 2 consists of nodes 1 2 3
    // face 3 consists of nodes 0 2 3

    int nodes[4];

    const Triangle* tris[4];

    Tetrahedron(int eid_, const int *nodes_, int matid_) : elementID(eid_), materialID(matid_)
    {
        for(int i=0; i<4; i++)
        {
            boundaryID[i] = 0;
            nodes[i] = nodes_[i];
            tris[i] = NULL;
        }
    }



    // If triangle tri is identical to a face of this tetra, then the corresponding face gets same boundaryID as tri.
    // Nothing happens, if tri is not identical to any face of this tetra.
    inline void setBoundaryID(const Triangle & tri)
    {
        bool triIncludesNode[4];

        triIncludesNode[0] = (tri.nodes[0]==nodes[0]) || (tri.nodes[1]==nodes[0]) || (tri.nodes[2]==nodes[0]);
        triIncludesNode[1] = (tri.nodes[0]==nodes[1]) || (tri.nodes[1]==nodes[1]) || (tri.nodes[2]==nodes[1]);
        triIncludesNode[2] = (tri.nodes[0]==nodes[2]) || (tri.nodes[1]==nodes[2]) || (tri.nodes[2]==nodes[2]);
        triIncludesNode[3] = (tri.nodes[0]==nodes[3]) || (tri.nodes[1]==nodes[3]) || (tri.nodes[2]==nodes[3]);

             if(triIncludesNode[0] && triIncludesNode[1] && triIncludesNode[2]){ boundaryID[0] = tri.boundaryID;  tris[0] = &tri; }
        else if(triIncludesNode[0] && triIncludesNode[1] && triIncludesNode[3]){ boundaryID[1] = tri.boundaryID;  tris[1] = &tri; }
        else if(triIncludesNode[1] && triIncludesNode[2] && triIncludesNode[3]){ boundaryID[2] = tri.boundaryID;  tris[2] = &tri; }
        else if(triIncludesNode[0] && triIncludesNode[2] && triIncludesNode[3]){ boundaryID[3] = tri.boundaryID;  tris[3] = &tri; }
    }
};





/*
 *  Retrieves the next line with a maximum length
 *  of maxlen from stream fid and stores it in line
 */
int mygetline(FILE *fid, char *line, int maxlen)
{
    char c = fgetc(fid);
    int pos = 0;

    while(c!=EOF && pos < maxlen-1)
    {
        if(c=='\n')
            break;

        line[pos] = c;
        pos++;

        c = fgetc(fid);
    }

    line[pos] = '\0';
    return 0;
}





/*
 *  Checks if the given line begins
 *  with the specified keyword
 *  (ignores white spaces and tabs at the beginning)
 */
static int beginsWith(char *line, const char *keyword)
{
    int pos = 0;

    while (*line == ' ' || *line == '\t') line++;

    while (*line == keyword[pos] && *line != '\0' && keyword[pos] != '\0')
    {
        line++;
        pos++;
    }

    if (keyword[pos] == '\0' && (*line == ' ' || *line == '\t' || *line == '\0'))
        return 1;
    else
        return 0;
}





/*
 *  Removes white spaces, tabs and " at
 *  the beginning and at the end of s
 */
char *trim(char *s)
{
    if(s==NULL || *s=='\0')
        return NULL;

    while (*s == ' ' || *s == '\t' || *s == '\"')
        s++;

    char *begin = s;

    while (*s != '\"')
        s++;

    *s = '\0';

    return begin;
}



// reads the first integer from string s, writes it to *i, and returns the remaining portion of s
char *getFirstInt(char *s, int *i)
{
    while(*s == ' ' || *s == '\t') s++;

    *i = atoi(s);

    while(*s != ' ' && *s != '\t' && *s != '\0') s++;

    return s;
}


void printUsage(char *binaryName)
{
    std::cout << "Converter GMSH to MIXD for 3D meshes" << std::endl;
    std::cout << "Usage: " << binaryName << " [-m {st|sd} -d {entity number}] <GMSH file>" << std::endl;
}


int main(int argc, char **argv)
{
    bool generate_st = false;  // generate space-time mesh?
    int  entity_dbl  = -1;     // index of physical entity to be doubled for internal BC

    char c;

    while ((c = getopt(argc, argv, "m:d:")) != -1)
    {
        switch(c)
        {
            case 'm':
                if(strcmp(optarg, "st")==0 || strcmp(optarg, "spacetime")==0)
                    generate_st = true;
                else if(strcmp(optarg, "sd")==0 || strcmp(optarg, "semidiscrete")==0)
                    generate_st = false;
                else{ printUsage(argv[0]); return 1; }
                break;

            case 'd':
                entity_dbl = atoi(optarg);
                break;

            default:
                printUsage(argv[0]);
                return 1;
        }
    }


    if(optind != argc-1)
    {
        printUsage(argv[0]);
        return 1;
    }

    const char *fname = argv[optind];

    std::cout << "Reading file " << fname << std::endl;



    FILE *fid = fopen(fname, "r");

    if (fid == NULL)
    {
        printf("\nCould not open %s\n", fname);
        return 1;
    }


    int linelen = 1024;
    char line[linelen];



    //============================================================================
    // PHYSICAL GROUP NAMES
    //============================================================================
    mygetline(fid, line, linelen);

    // find the physical names section in the msh file.
    // CAUTION: this assumes that the section is located somewhere at the beginning
    // of the file!!! ... might also be placed somewhere at the end!!!
    while(!beginsWith(line, "$PhysicalNames"))
        mygetline(fid, line, linelen);


    mygetline(fid, line, linelen);

    int n_phys_grp = atoi(line);

    std::cout << "Found physical names: " << n_phys_grp << std::endl;


    int *phys_grp_dim = new int[n_phys_grp];
    std::string *phys_grp_names = new std::string[n_phys_grp];


    // mapping from the GMSH name group indices to MIXD material indices
    std::map<int,int> fromGmshToMixdMatIDs;

    int matid = 1;


    for(int i=0; i<n_phys_grp; i++)
    {
        int dim, ind;
        char name[512];

        fscanf(fid, "%d %d %s", &dim, &ind, name);

        char *tmp = trim(name);

        phys_grp_dim  [ind-1] = dim;
        phys_grp_names[ind-1] = std::string(tmp);

        if(dim == 3)
        {
            fromGmshToMixdMatIDs.insert(std::pair<int,int>(ind,matid));
            matid++;
        }

    }
    //============================================================================






    //============================================================================
    // NODES
    //============================================================================
    mygetline(fid, line, linelen);

    while(!beginsWith(line, "$Nodes"))
        mygetline(fid, line, linelen);


    mygetline(fid, line, linelen);

    int n_nodes = atoi(line);

    std::cout << "Found nodes: " << n_nodes << std::endl;


    double *nodes_x = new double[n_nodes];
    double *nodes_y = new double[n_nodes];
    double *nodes_z = new double[n_nodes];

    // mapping from the old GMSH indices to the new MIXD indices
    std::map<int,int> fromGmshToMixdNodeIDs;

    SimpleProgress spn(0, n_nodes, 20);

    for(int i=0; i<n_nodes; i++)
    {
        int ind;
        double x, y, z;

        fscanf(fid, "%d %lf %lf %lf", &ind, &x, &y, &z);

        nodes_x[i] = x;
        nodes_y[i] = y;
        nodes_z[i] = z;

        fromGmshToMixdNodeIDs.insert(std::pair<int,int>(ind, i+1));

        spn.printIfHitNext(i);
    }
    std::cout << std::endl;
    //============================================================================







    //============================================================================
    // ELEMENTS
    //============================================================================
    mygetline(fid, line, linelen);

    while(!beginsWith(line, "$Elements"))
        mygetline(fid, line, linelen);

    mygetline(fid, line, linelen);

    int n_elems = atoi(line);

    std::cout << "Found elements: " << n_elems << std::endl;


    // triangles are stored in a map for fast access O(log N)
    // they are sorted by their unique keys
    std::map<TriKey, Triangle, compareTriKey> tris;

    // a simple vector is sufficient for the tetras, since we need only sequential access
    std::vector<Tetrahedron> tets;

    // boundary nodes are sorted by their IDs ==> set gives logN access
    std::set<int> boundaryNodes;
    std::set<int> dblBoundaryNodes;


    int tet_id = 1;

    SimpleProgress spe(0, n_elems, 20);

    for(int i=0; i<n_elems; i++)
    {
        mygetline(fid, line, linelen);


        char *s = line;

        int id, type, ntags;

        s = getFirstInt(s, &id);
        s = getFirstInt(s, &type);
        s = getFirstInt(s, &ntags);

//         std::cout << "id = " << id << "   type = " << type << "   ntags = " << ntags << std::endl;

        int tags[ntags];

        for(int t=0; t<ntags; t++)
            s = getFirstInt(s, &(tags[t]));

        int n_elem_nodes = 0;

        switch(type)
        {
            case 2:    n_elem_nodes = 3;  break;  // triangle
//             case 3:    n_elem_nodes = 4;  break;  // quadrangle
            case 4:    n_elem_nodes = 4;  break;  // tetrahedron
//             case 5:    n_elem_nodes = 8;  break;  // hexahedron
//             case 6:    n_elem_nodes = 6;  break;  // prism
//             case 7:    n_elem_nodes = 5;  break;  // pyramid
            default:   std::cout << "ERROR: Found unknown element type: " << type << std::endl; return 1;
        }


        int elem_nodes[n_elem_nodes];

        int tmp;

        for(int n=0; n<n_elem_nodes; n++)
        {
            s = getFirstInt(s, &tmp);
            elem_nodes[n] = fromGmshToMixdNodeIDs.find(tmp)->second;
        }


        if(type == 2) // triangle
        {
            // tags[0] is the physical entity of this element!
            tris.insert(std::pair<TriKey,Triangle>(TriKey(elem_nodes), Triangle(id, tags[0], elem_nodes)));

//             if(tags[0] == 1)
//                 std::cout << "found inlet" << std::endl;

            boundaryNodes.insert(elem_nodes[0]);
            boundaryNodes.insert(elem_nodes[1]);
            boundaryNodes.insert(elem_nodes[2]);

            // if we are on a boundary where nodes are to be doubled,
            // then insert these face's nodes into the doubled set.
            if(tags[0] == entity_dbl)
            {
                dblBoundaryNodes.insert(elem_nodes[0]);
                dblBoundaryNodes.insert(elem_nodes[1]);
                dblBoundaryNodes.insert(elem_nodes[2]);
            }
        }

        else if(type == 4) // tetrahedron
        {
            // find material ID of this tetra element
            int matid = fromGmshToMixdMatIDs.find(tags[0])->second;
            tets.push_back(Tetrahedron(tet_id, elem_nodes, matid));
            tet_id++;
        }

        spe.printIfHitNext(i);

    }

    std::cout << std::endl;

    std::cout << "number of tris:  " << tris.size() << std::endl;
    std::cout << "number of tets: "  << tets.size() << std::endl << std::endl;
    std::cout << "boundary nodes:  "  << boundaryNodes.size() << std::endl << std::endl << std::endl;


    fclose(fid);


//     int n_bound_tets = 0;
//     int zero  = 0;
//     int one   = 0;
//     int two   = 0;
//     int three = 0;
//     int four  = 0;
//     int other = 0;


//     int n_found_tris = 0;


    std::cout << "now relating boundary conditions..." << std::endl;

    // This set stores the (MIXD) material IDs of all element domains
    // which touch the boundary with doubled nodes
    std::set<int> dblNodeDomains;

    // let's iterate over all tetras to find their boundary affiliation
    for (std::vector<Tetrahedron>::iterator it_tet=tets.begin(); it_tet!=tets.end(); ++it_tet)
    {
        //------------------------------------------------------------
        // how many of this tetra's nodes are on boundaries?
        // --> this is cheap: we search in a relatively small
        //     binary tree with cost O(log(N)), where N is max. 1e5
        int n_bound_nodes = 0;

        for(int n=0; n<4; n++)
        {
            n_bound_nodes += boundaryNodes.count(it_tet->nodes[n]);

            // is any of this tetra's nodes on the boundary of doubled nodes?
            if(dblBoundaryNodes.count(it_tet->nodes[n]))
                // if yes, remember this tetra's material as touching the boundary of double nodes
                dblNodeDomains.insert(it_tet->materialID);
        }
        //------------------------------------------------------------


        //----------------------------------------------------------------------------------------------
        // Tetra might have boundary faces, only if at least 3 of its nodes are on boundaries!
        // The vast majority of all tetras are not connected to the boundary, and thus n_bound_nodes<=2.
        // So we can just skip those
        if(n_bound_nodes > 2)
        {
            // get the nodes of all 4 faces of this tetra, i.e. the potential boundary triangles
            int faceNodes1[3] = { it_tet->nodes[0], it_tet->nodes[1], it_tet->nodes[2] };
            int faceNodes2[3] = { it_tet->nodes[0], it_tet->nodes[1], it_tet->nodes[3] };
            int faceNodes3[3] = { it_tet->nodes[0], it_tet->nodes[2], it_tet->nodes[3] };
            int faceNodes4[3] = { it_tet->nodes[1], it_tet->nodes[2], it_tet->nodes[3] };

            // generate the unique search keys of the four candidate triangles
            TriKey faces[4] = { TriKey(faceNodes1), TriKey(faceNodes2), TriKey(faceNodes3), TriKey(faceNodes4) };

            // loop over the four candidates...
            for(int f=0; f<4; f++)
            {
                // ... and check if the candidate is part of the boundary triangles list
                // (again a cheap bin tree search with O(log(N)) )
                std::map<TriKey,Triangle>::iterator it_tri = tris.find(faces[f]);

                // If we find the candidate in the list, then
                // we set the boundary condition for this tetra
                if(it_tri != tris.end())
                {
                    it_tri->second.setAdjacentTetra(it_tet->elementID);
                    it_tet->setBoundaryID(it_tri->second);
                }
            }
        }
        //----------------------------------------------------------------------------------------------
    }

    std::cout << "boundary conditions found!" << std::endl;

    // exactly two domains must be touching the boundary of double nodes
    if(entity_dbl!=-1 && dblNodeDomains.size() != 2)
    {
        std::cout << "ERROR: Boundary with nodes to be doubled is touched by more or less than 2 material domains!" << std::endl;
        std::cout << "Number of touching domains: " << dblNodeDomains.size() << std::endl;
        return 1;
    }


    std::map<int,int> fromSingleToDoubleNodeIDs;

    if(entity_dbl!=-1)
    {


    // set the material domain with the largest number as the one which uses the doubled nodes on the boundary
    int domainOfDoubledNodes = -1;
    for(std::set<int>::iterator it=dblNodeDomains.begin(); it!=dblNodeDomains.end(); ++it)
    {
        if(*it > domainOfDoubledNodes)
            domainOfDoubledNodes = *it;
    }





    int nid = n_nodes+1;

    for(std::set<int>::iterator it=dblBoundaryNodes.begin(); it!=dblBoundaryNodes.end(); ++it)
    {
        fromSingleToDoubleNodeIDs.insert(std::pair<int,int>(*it, nid));

        nid++;
    }



    std::cout << "now doubling nodes..." << std::endl;

    // iterate over the tetras to change their connectivity to doubled boundary nodes
    for (std::vector<Tetrahedron>::iterator it_tet=tets.begin(); it_tet!=tets.end(); ++it_tet)
    {
        // should this tetra be connected to doubled nodes?
        if(it_tet->materialID == domainOfDoubledNodes)
        {
            // iterate over tetras nodes
            for(int n=0; n<4; n++)
            {
                // is the current node on the boundary of doubled nodes?
                if(dblBoundaryNodes.count(it_tet->nodes[n]))
                {
                    // change node ID to doubled ID
                    it_tet->nodes[n] = fromSingleToDoubleNodeIDs.find(it_tet->nodes[n])->second;
                }
            }
        }
    }

    std::cout << "doubling of nodes done!" << std::endl;

    }












    // now let's start writing the MIXD files

    std::cout << "writing minf... ";

    std::ofstream of;

    of.open("minf");

    of << "ne" << std::setw(12) << tets.size() << std::endl;
    if(!generate_st)
        of << "nn" << std::setw(12) <<   (n_nodes+dblBoundaryNodes.size()) << std::endl;
    else
        of << "nn" << std::setw(12) << 2*(n_nodes+dblBoundaryNodes.size()) << std::endl;

    of.close();

    std::cout << "done!" << std::endl;


    // what machine type are we on?
    bool isBigEnd = isBigEndian();

    // buffers for swaping bytes and writing to binary file
//     char *id_buf = (char*) malloc(4);    // IDs are 4 byte ints

    char *en_buf = (char*) malloc(4*4);  // buffer for element node lists --> 4 ints of 4 bytes each
    char *ef_buf = (char*) malloc(4*4);  // buffer for element face lists --> 4 ints of 4 bytes each
    char *em_buf = (char*) malloc(4);    // buffer for element material ID--> 4 byte int

    char *nx_buf = (char*) malloc(8);    // node x coord buffer --> 8 byte double
    char *ny_buf = (char*) malloc(8);    // node y coord buffer --> 8 byte double
    char *nz_buf = (char*) malloc(8);    // node z coord buffer --> 8 byte double


    std::cout << "writing mien... ";

    fid = fopen("mien", "wb");

    if (fid == NULL)
    {
        printf("\nCould not open %s\n", "mien");
        return 1;
    }


//     of.open("mien_ascii");

    for (std::vector<Tetrahedron>::iterator it_tet=tets.begin(); it_tet!=tets.end(); ++it_tet)
    {
//         //=== ASCII =======================================
//         of << std::setw(12) << it_tet->elementID;
//
//         for(int n=0; n<4; n++)
//             of << std::setw(12) << it_tet->nodes[n];
//
//         of << std::endl;
//         //=================================================


        //=== BINARY ======================================
//         memcpy(id_buf, &(it_tet->elementID), 4);
        memcpy(en_buf, it_tet->nodes, 4*4);

        if(!isBigEnd)  // swap bytes if neccessary
        {
//             swapbytes(id_buf, 1, 4);
            swapbytes(en_buf, 4, 4);
        }

//         fwrite(id_buf, 4, 1, fid);
        fwrite(en_buf, 4, 4, fid);
        //=================================================
    }

//     of.close();

    fclose(fid);

    std::cout << "done!" << std::endl;




    std::cout << "writing mmat... ";

    fid = fopen("mmat", "wb");

    if (fid == NULL)
    {
        printf("\nCould not open %s\n", "mmat");
        return 1;
    }


//     of.open("mien_ascii");

    for (std::vector<Tetrahedron>::iterator it_tet=tets.begin(); it_tet!=tets.end(); ++it_tet)
    {
//         //=== ASCII =======================================
//         of << std::setw(12) << it_tet->elementID;
//
//         for(int n=0; n<4; n++)
//             of << std::setw(12) << it_tet->nodes[n];
//
//         of << std::endl;
//         //=================================================


        //=== BINARY ======================================
//         memcpy(id_buf, &(it_tet->elementID), 4);
        memcpy(em_buf, &(it_tet->materialID), 4);

        if(!isBigEnd)  // swap bytes if neccessary
        {
//             swapbytes(id_buf, 1, 4);
            swapbytes(em_buf, 1, 4);
        }

//         fwrite(id_buf, 4, 1, fid);
        fwrite(em_buf, 4, 1, fid);
        //=================================================
    }

//     of.close();

    fclose(fid);

    std::cout << "done!" << std::endl;



    std::cout << "writing mxyz... ";

    fid = fopen("mxyz", "wb");

    if (fid == NULL)
    {
        printf("\nCould not open %s\n", "mxyz");
        return 1;
    }

//     of.open("mxyz_ascii");

    int nloops = (generate_st)? 2 : 1;

    for(int loop=0; loop<nloops; loop++)
    {
    for(int n=0; n<n_nodes; n++)
    {
//         int id = n+1;

//         //=== ASCII =======================================
//         of << std::setw(12) << id;
//
//         of << std::setw(30) << std::scientific << std::setprecision(12) << nodes_x[n];
//         of << std::setw(30) << std::scientific << std::setprecision(12) << nodes_y[n];
//         of << std::setw(30) << std::scientific << std::setprecision(12) << nodes_z[n];
//
//         of << std::endl;
//         //=================================================


        //=== BINARY ======================================
//         memcpy(id_buf, &id, 4);
        memcpy(nx_buf, &(nodes_x[n]), 8);
        memcpy(ny_buf, &(nodes_y[n]), 8);
        memcpy(nz_buf, &(nodes_z[n]), 8);

        if(!isBigEnd)  // swap bytes if neccessary
        {
//             swapbytes(id_buf, 1, 4);
            swapbytes(nx_buf, 1, 8);
            swapbytes(ny_buf, 1, 8);
            swapbytes(nz_buf, 1, 8);
        }

//         fwrite(id_buf, 4, 1, fid);
        fwrite(nx_buf, 8, 1, fid);
        fwrite(ny_buf, 8, 1, fid);
        fwrite(nz_buf, 8, 1, fid);
        //=================================================
    }

    for(std::set<int>::iterator it=dblBoundaryNodes.begin(); it!=dblBoundaryNodes.end(); ++it)
    {

        memcpy(nx_buf, &(nodes_x[*it-1]), 8);
        memcpy(ny_buf, &(nodes_y[*it-1]), 8);
        memcpy(nz_buf, &(nodes_z[*it-1]), 8);


        if(!isBigEnd)  // swap bytes if neccessary
        {

            swapbytes(nx_buf, 1, 8);
            swapbytes(ny_buf, 1, 8);
            swapbytes(nz_buf, 1, 8);

        }

        fwrite(nx_buf, 8, 1, fid);
        fwrite(ny_buf, 8, 1, fid);
        fwrite(nz_buf, 8, 1, fid);


    }


    }

//     of.close();

    fclose(fid);

    std::cout << "done!" << std::endl;








    if(entity_dbl != -1)
    {

	    std::cout << "writing mtbl... ";
    fid = fopen("mtbl", "wb");

    if (fid == NULL)
    {
        printf("\nCould not open %s\n", "mtbl");
        return 1;
    }

//     of.open("mxyz_ascii");

    int nloops = (generate_st)? 2 : 1;

    int partner_node = 0;

    for(int loop=0; loop<nloops; loop++)
    {
    for(int n=0; n<n_nodes; n++)
    {
//         int id = n+1;

//         //=== ASCII =======================================
//         of << std::setw(12) << id;
//
//         of << std::setw(30) << std::scientific << std::setprecision(12) << nodes_x[n];
//         of << std::setw(30) << std::scientific << std::setprecision(12) << nodes_y[n];
//         of << std::setw(30) << std::scientific << std::setprecision(12) << nodes_z[n];
//
//         of << std::endl;
//         //=================================================


        if(dblBoundaryNodes.count(n+1))
        {
            partner_node = fromSingleToDoubleNodeIDs.find(n+1)->second;
        }
        else
            partner_node = n+1;

	  partner_node += loop*(n_nodes + dblBoundaryNodes.size());

        //=== BINARY ======================================
//         memcpy(id_buf, &id, 4);
        memcpy(em_buf, &(partner_node), 4);
//         memcpy(ny_buf, &(nodes_y[n]), 8);
//         memcpy(nz_buf, &(nodes_z[n]), 8);

        if(!isBigEnd)  // swap bytes if neccessary
        {
//             swapbytes(id_buf, 1, 4);
            swapbytes(em_buf, 1, 4);
//             swapbytes(ny_buf, 1, 8);
//             swapbytes(nz_buf, 1, 8);
        }

//         fwrite(id_buf, 4, 1, fid);
        fwrite(em_buf, 4, 1, fid);
//         fwrite(ny_buf, 8, 1, fid);
//         fwrite(nz_buf, 8, 1, fid);
        //=================================================
    }

    for(std::set<int>::iterator it=dblBoundaryNodes.begin(); it!=dblBoundaryNodes.end(); ++it)
    {

        partner_node = *it + loop*(n_nodes + dblBoundaryNodes.size());

        memcpy(em_buf, &(partner_node), 4);
//         memcpy(ny_buf, &(nodes_y[n]), 8);
//         memcpy(nz_buf, &(nodes_z[n]), 8);

        if(!isBigEnd)  // swap bytes if neccessary
        {
//             swapbytes(id_buf, 1, 4);
            swapbytes(em_buf, 1, 4);
//             swapbytes(ny_buf, 1, 8);
//             swapbytes(nz_buf, 1, 8);
        }

//         fwrite(id_buf, 4, 1, fid);
        fwrite(em_buf, 4, 1, fid);


    }


    }

//     of.close();

    fclose(fid);

    std::cout << "done!" << std::endl;





    std::cout << "writing mtbl.dual... ";

    fid = fopen("mtbl.dual", "wb");

    if (fid == NULL)
    {
        printf("\nCould not open %s\n", "mtbl.dual");
        return 1;
    }


//     of.open("mien_ascii");

    for (std::vector<Tetrahedron>::iterator it_tet=tets.begin(); it_tet!=tets.end(); ++it_tet)
    {
//         //=== ASCII =======================================
//         of << std::setw(12) << it_tet->elementID;
//
//         for(int n=0; n<4; n++)
//             of << std::setw(12) << it_tet->nodes[n];
//
//         of << std::endl;
//         //=================================================

        int nbtet[4];
        for(int iface=0; iface<4; iface++)
        {
            if(it_tet->tris[iface] != NULL)
                nbtet[iface] = - (it_tet->tris[iface]->getNeighborTetra(it_tet->elementID));
            else
                nbtet[iface] = 0;
        }

        //=== BINARY ======================================
//         memcpy(id_buf, &(it_tet->elementID), 4);
        memcpy(ef_buf, nbtet, 4*4);

        if(!isBigEnd)  // swap bytes if neccessary
        {
//             swapbytes(id_buf, 1, 4);
            swapbytes(ef_buf, 4, 4);
        }

//         fwrite(id_buf, 4, 1, fid);
        fwrite(ef_buf, 4, 4, fid);
        //=================================================
    }

//     of.close();

    fclose(fid);

    std::cout << "done!" << std::endl;






    }




    std::cout << "writing mrng... ";

    fid = fopen("mrng", "wb");

    if (fid == NULL)
    {
        printf("\nCould not open %s\n", "mrng");
        return 1;
    }


//     of.open("mrng_ascii");

    for (std::vector<Tetrahedron>::iterator it_tet=tets.begin(); it_tet!=tets.end(); ++it_tet)
    {
//         //=== ASCII =======================================
//         of << std::setw(12) << it_tet->elementID;
//
//         for(int f=0; f<4; f++)
//             of << std::setw(6) << it_tet->boundaryID[f];
//
//         of << std::endl;
//         //=================================================


        //=== BINARY ======================================
//         memcpy(id_buf, &(it_tet->elementID), 4);
        memcpy(ef_buf, it_tet->boundaryID, 4*4);

        if(!isBigEnd)  // swap bytes if neccessary
        {
//             swapbytes(id_buf, 1, 4);
            swapbytes(ef_buf, 4, 4);
        }

//         fwrite(id_buf, 4, 1, fid);
        fwrite(ef_buf, 4, 4, fid);
        //=================================================
    }

//     of.close();

    fclose(fid);


    std::cout << "done!" << std::endl;














//     free(id_buf);

    free(en_buf);
    free(ef_buf);
    free(em_buf);

    free(nx_buf);
    free(ny_buf);
    free(nz_buf);













    delete[] phys_grp_dim;
    delete[] phys_grp_names;
    delete[] nodes_x;
    delete[] nodes_y;
    delete[] nodes_z;



    return 0;
}


