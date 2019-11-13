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

        std::set<int> boundaryNodes;        //All boundary nodes
        std::set<int> dblBoundaryNodes;     //doubled boundary nodes
        std::set<int> dblNodeDomains;       //mixd matID for all elements with doubled nodes
        std::map<int,int> fromSingleToDoubleNodeIDs;

        std::map<std::vector<int>, Triangle * > trimap;

        Mesh(){};
        ~Mesh(){};
        void addToTriMap(std::string line)
        {
            int eid, etype, ntags;
            std::istringstream iss(line);
            iss >> eid >> etype >> ntags;

            //only vertex nodes considered for sortednodes
            int nen = 3;

            int dummy;
            for (int i = 0; i < ntags; i++)
                iss >> dummy;

            std::vector<int> sortedNodes;

            for(int i=0; i<nen; i++)
            {
                iss >> dummy;
                sortedNodes.push_back(dummy);
            }
            std::sort(sortedNodes.begin(), sortedNodes.end());

            trimap.insert(std::pair<std::vector<int>, Triangle *>(sortedNodes, new Triangle(line)) );



        }
};

#endif
