#include <iostream>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <unistd.h>
#include <vector>

#include "physicalGroup.h"
#include "node.h"
#include "elements.h"
#include "writeMIXD.h"
#include "mesh.h"

// TODO: test 2-order meshes
// TODO: stop using endian library? be consistent with write functions?
// TODO: change from templating ORDER to using it as a function parameter or mesh attribute.
// NOTE: Spacetime features are not required because the relevant mtbl/mprd files just have to doubled. XNS handles the +nnspace offset within the code.

int readfile(std::string filename, Mesh &mesh)
{
    std::ifstream infile(filename);

    std::string line;
    std::string section;
    std::string section_size;

    //MESH FORMAT
    std::getline(infile, section);
    std::getline(infile, section);
    std::getline(infile, section);

    //Physical Names
    std::cout << "Reading Physical Groups... " << std::flush;
    std::getline(infile, section);
    std::getline(infile, section_size);
    char *end; 
    size_t num = (size_t) std::strtoull(section_size.c_str(), &end, 10);
    int matid = 0;
    for (size_t i=0; i<num; i++)
    {
        int dim, tag;
        std::string name;
        std::getline(infile, line);
        std::istringstream iss(line);
        /* std::cout << line << std::endl; */
        iss >> dim >> tag >> name;
        if (dim == 3)
            matid++;
        mesh.pg.push_back(new PhysicalGroup(dim, tag, name, matid));

    }

    std::getline(infile, section);
    std::cout << "done!" << std::endl;


    //Nodes
    std::cout << "Reading Nodes... " << std::flush;
    std::getline(infile, section);
    std::getline(infile, section_size);
    num = (size_t) std::strtoull(section_size.c_str(), &end, 10);

    for (size_t i=0; i<num; i++)
    {
        std::getline(infile,line);
        /* std::cout << line << std::endl; */
        mesh.nodes.push_back(new Node(line));
    }
    std::getline(infile, section);
    std::cout << "done!" << std::endl;


    //Elements
    std::cout << "Reading Elements... " << std::flush;
    std::getline(infile, section);
    std::getline(infile, section_size);
    num = (size_t) std::strtoull(section_size.c_str(), &end, 10);

    for (size_t i=0; i<num; i++)
    {
        std::getline(infile,line);
        /* std::cout << line << std::endl; */

        size_t eid, etype, ntags;
        std::istringstream iss(line);
        iss >> eid >> etype >> ntags;

        int * test;


        switch(etype)
        {
            case 2:
            /* case 9: mesh.tris.push_back(new Triangle(line)); break; */
            case 9: mesh.addToTriMap(line);break;
            case 4:
            case 11: mesh.tets.push_back(new Tetrahedron(line)); break;
            default: std::cout << "ERROR: Unknown element type: " << etype << std::endl; return 1;
        }

    }
    std::getline(infile, section);
    std::cout << "done!" << std::endl <<std::endl;;

    std::cout << "number of physical groups: " << mesh.pg.size() << std::endl;;
    std::cout << "number of nodes: " << mesh.nodes.size() << std::endl;;
    /* std::cout << "number of elements: " << mesh.tris.size() + mesh.tets.size() << std::endl; */
    std::cout << "number of elements: " << mesh.trimap.size() + mesh.tets.size() << std::endl;
    std::cout << "number of tris:  " << mesh.trimap.size() << std::endl;
    std::cout << "number of tets: "  << mesh.tets.size() << std::endl << std::endl;

    return 0;
}

int boundaryConditions(Mesh &mesh)
{

    //Check all triangles for boundary nodes and interface(doubled) nodes
    for(std::map<std::vector<size_t>, Triangle *>::iterator it_tri=mesh.trimap.begin(); it_tri!=mesh.trimap.end(); ++it_tri)
    {
        for(size_t i=0; i<it_tri->second->nen; i++)
        {
            mesh.boundaryNodes.insert((it_tri)->second->nodes[i]);
            if ((it_tri)->second->bID == mesh.entity_dbl)
                mesh.dblBoundaryNodes.insert((it_tri)->second->nodes[i]);
        }

    }

    std::cout << "number of boundary nodes: " << mesh.boundaryNodes.size() << std::endl;
    std::cout << "number of doubled nodes: " << mesh.dblBoundaryNodes.size() << std::endl;

    //Check all tets for boundary nodes
    std::cout << "matching tets and tris... " << std::flush;
    #pragma omp parallel for shared(mesh)
    for(std::vector<Tetrahedron *>::iterator it_tet=mesh.tets.begin(); it_tet!=mesh.tets.end(); ++it_tet)
    {
        size_t n_bound_nodes = 0;
        for(size_t i = 0; i< (*it_tet)->nen; i++)
        {
            n_bound_nodes += mesh.boundaryNodes.count((*it_tet)->nodes[i]);
            if(mesh.dblBoundaryNodes.count((*it_tet)->nodes[i]))
            {
                #pragma omp critical
                mesh.dblNodeDomains.insert((*it_tet)->matID);
            }

        }

        //faces on the boundary
        if (n_bound_nodes > 2)
        {
            /* std::cout << "n_bound_nodes = " << n_bound_nodes << std::endl; */
            //uses only vertex nodes for faces
            std::vector<std::vector<size_t>> faces;
            faces.push_back({(*it_tet)->nodes[0], (*it_tet)->nodes[1], (*it_tet)->nodes[2]});
            faces.push_back({(*it_tet)->nodes[0], (*it_tet)->nodes[1], (*it_tet)->nodes[3]});
            faces.push_back({(*it_tet)->nodes[0], (*it_tet)->nodes[2], (*it_tet)->nodes[3]});
            faces.push_back({(*it_tet)->nodes[1], (*it_tet)->nodes[2], (*it_tet)->nodes[3]});


            for(size_t f=0; f<4; f++)
            {
                //sort nodes of the faces and find matching tris.
                //uses only vertex nodes
                std::sort(faces[f].begin(), faces[f].end());

                /* auto it_tri = std::find_if(mesh.tris.begin(), mesh.tris.end(), [&faces,&f](Triangle * t){ return t->sortedNodes==faces[f]; }); */

                std::map<std::vector<size_t>,Triangle*>::iterator it_tri = mesh.trimap.find(faces[f]);

                //if tri is found, set boundary ID and adjacent tets
                if(it_tri != mesh.trimap.end())
                {
                    // std::cout << "Found match bID: " << it_tri->second->bID << std::endl;
                    (it_tri)->second->setAdjacentTetra(it_tet - mesh.tets.begin() + 1);
                    /* (*it_tet)->setBoundaryID(it_tri->second); */
                    // TODO: This is probably unnecessary since we explicitly find the matched nodes already. Simplify the function
                    (*it_tet)->setBoundaryID(*it_tri);

                }
            }

        }

    }
    std::cout << "done!" << std::endl;

    // exactly two domains must be touching the boundary of double nodes
    if(mesh.entity_dbl!=-1 && mesh.dblNodeDomains.size() != 2)
    {
        std::cout << "ERROR: Boundary with nodes to be doubled is touched by more or less than 2 material domains!" << std::endl;
        std::cout << "Number of touching domains: " << mesh.dblNodeDomains.size() << std::endl;
        return 1;
    }

    if(mesh.entity_dbl!=-1)
    {

        // set the material domain with the largest number as the one which uses the doubled nodes on the boundary
        // TODO: mesh.dblNodeDomains.max()
        int domainOfDoubledNodes = -1;
        for(std::set<size_t>::iterator it=mesh.dblNodeDomains.begin(); it!=mesh.dblNodeDomains.end(); ++it)
        {
            if(*it > domainOfDoubledNodes)
                domainOfDoubledNodes = *it;
        }

        size_t n_nodes = mesh.nodes.size();
        size_t nid = n_nodes+1;

        for(std::set<size_t>::iterator it=mesh.dblBoundaryNodes.begin(); it!=mesh.dblBoundaryNodes.end(); ++it)
        {
            mesh.fromSingleToDoubleNodeIDs.insert(std::pair<size_t,size_t>(*it, nid));
            nid++;
        }

        std::cout << "now doubling nodes..." << std::flush;

        // iterate over the tetras to change their connectivity to doubled boundary nodes
        #pragma omp parallel for shared(mesh)
        for (std::vector<Tetrahedron *>::iterator it_tet=mesh.tets.begin(); it_tet!=mesh.tets.end(); ++it_tet)
        {
            // should this tetra be connected to doubled nodes?
            if( (*it_tet)->matID == domainOfDoubledNodes)
            {
                // iterate over tetras nodes
                for(size_t n=0; n<(*it_tet)->nen; n++)
                {
                    // is the current node on the boundary of doubled nodes?
                    if(mesh.dblBoundaryNodes.count((*it_tet)->nodes[n]))
                    {
                        // change node ID to doubled ID
                        (*it_tet)->nodes[n] = mesh.fromSingleToDoubleNodeIDs.find((*it_tet)->nodes[n])->second;
                    }
                }
            }
        }

        std::cout << "done!" << std::endl;

    }

    if (mesh.order==2)
    {
        // iterate over the tetras to change the node order to match xns
        for (std::vector<Tetrahedron *>::iterator it_tet=mesh.tets.begin(); it_tet!=mesh.tets.end(); ++it_tet)
        {
            std::swap((*it_tet)->nodes[8], (*it_tet)->nodes[9]);
        }

    }

    return 0;
}

void printUsage(char *binaryName)
{
    std::cout << "Converter GMSH to MIXD for 3D meshes" << std::endl;
    std::cout << "Usage: " << binaryName << " [-m {st|sd} -d {entity number} -o {element order}] <GMSH file>" << std::endl;
}


int main(int argc, char * argv[])
{

    bool generate_st = false;  // generate space-time mesh?
    int  entity_dbl  = -1;     // index of physical entity to be doubled for internal BC
    int order = 1;

    char c;

    while ((c = getopt(argc, argv, "m:d:o:")) != -1)
    {
        switch(c)
        {
            case 'm':
                if(strcmp(optarg, "st")==0 || strcmp(optarg, "spacetime")==0)
                {
                    std::cout << "Spacetime meshes are not supported." << std::endl;
                    std::cout << "XNS handles the +nnspace offset internally." << std::endl;
                    std::cout << "Doubling the semidiscrete mxyz, mtbl, and mprd files are sufficient for conversion to spacetime." << std::endl;
                    exit(-1);
                    /* generate_st = true; */
                }
                else if(strcmp(optarg, "sd")==0 || strcmp(optarg, "semidiscrete")==0)
                    generate_st = false;
                else{ printUsage(argv[0]); return 1; }
                break;
            case 'd':
                entity_dbl = atoi(optarg);
                break;
            case 'o':
                order = atoi(optarg);
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

    std::string fname = argv[optind];

    Mesh mesh;
    mesh.order = order;
    mesh.entity_dbl = entity_dbl;

    readfile(fname, mesh);
    boundaryConditions(mesh);
    writeMIXD(mesh);

    return 0;
}
