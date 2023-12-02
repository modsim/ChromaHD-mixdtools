#ifndef MESH_H
#define MESH_H

#include<vector>
#include<set>
#include<map>

#include "elements.h"
#include "node.h"
/* #include "physicalGroup.h" */

class Mesh
{
    public:
        int order = 1;
        int entity_dbl = -1;
        std::vector<Tetrahedron *> tets;
        std::vector<Node *> nodes;
        std::vector<PhysicalGroup *> pg;

        std::set<size_t> boundaryNodes;        //All boundary nodes
        std::set<size_t> dblBoundaryNodes;     //doubled boundary nodes
        std::set<size_t> dblNodeDomains;       //mixd matID for all elements with doubled nodes
        std::map<size_t,size_t> fromSingleToDoubleNodeIDs;

        std::map<std::vector<size_t>, Triangle * > trimap;

        Mesh(){};
        ~Mesh(){};
        void addToTriMap(std::string line)
        {
            size_t eid, etype, ntags;
            std::istringstream iss(line);
            iss >> eid >> etype >> ntags;

            //only vertex nodes considered for sortednodes
            int nen = 3;

            size_t dummy;
            for (size_t i = 0; i < ntags; i++)
                iss >> dummy;

            std::vector<size_t> sortedNodes;

            for(size_t i=0; i<nen; i++)
            {
                iss >> dummy;
                sortedNodes.push_back(dummy);
            }
            std::sort(sortedNodes.begin(), sortedNodes.end());

            trimap.insert(std::pair<std::vector<size_t>, Triangle *>(sortedNodes, new Triangle(line)) );



        }
};

#endif
