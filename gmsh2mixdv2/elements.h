#ifndef ELEMENTS_H
#define ELEMENTS_H

#include <memory>
#include <vector>
#include "physicalGroup.h"
#include <iostream>
#include <sstream>
#include <bits/stdc++.h>

class Triangle
{
    public:
        int eID;
        int bID;
        int nen;


        std::vector<int> nodes;
        std::vector<int> sortedNodes;
        /* std::vector<int> adjacentTetras; */
        int adjacentTetras[2];
        int nAdjoiningTetras;

        /* Triangle(int eID_, int _bid_, int * nodes_); */
        Triangle(std::string line);
        inline void setAdjacentTetra(int tetID_);
        inline int getNeighbourTetra(int tetID_);

        ~Triangle();

};


class Tetrahedron
{
    public:
        /* static std::vector<PhysicalGroup *>& pg; */
        int eID;
        int bID[4];
        int matID;
        int pgID;
        int nen;

        std::vector<int> nodes;
        std::vector<Triangle *> tris;
        std::vector<std::vector<int>> faces;

        /* Tetrahedron(int eID_, int * nodes_, int matID_); */
        Tetrahedron(std::string line);
        /* Tetrahedron(std::vector<PhysicalGroup *>& pg_); */
        ~Tetrahedron();
        /* inline void setBoundaryID(Triangle * tri); */
        inline void setBoundaryID(std::pair<std::vector<int>,Triangle *> tripair);
};



Triangle::Triangle(std::string line)
{
    int etype, ntags;
    std::istringstream iss(line);
    iss >> eID >> etype >> ntags;

    switch (etype)
    {
        case 2: nen = 3; break;
        case 9: nen = 6; break;
        default: std::cout<< "Triangle element order unsupported: " << etype << std::endl;
    }

    iss >> bID;

    // std::cout << "Found bID: " << bID << std::endl;

    int dummy;
    for (int i = 1; i < ntags; i++)
        iss >> dummy;

    for(int i=0; i<nen; i++)
    {
        iss >> dummy;
        nodes.push_back(dummy);
    }

    /* //Only uses vertex nodes */
    /* for(int i=0; i<3; i++) */
    /* { */
    /*     sortedNodes.push_back(nodes[i]); */
    /* } */
    /* std::sort(sortedNodes.begin(), sortedNodes.end()); */

    nAdjoiningTetras  = 0;
    adjacentTetras[0] = 0;
    adjacentTetras[1] = 0;

}


inline void Triangle::setAdjacentTetra(int tetID_)
{
    if(nAdjoiningTetras>1)
    {
        std::cout << "FATAL ERROR in setAdjacentTetra!" << std::endl;
        std::cout << "Triangle nodes: " << nodes[0] << "    " << nodes[1] << "    " << nodes[2] << std::endl;
        std::cout << "boundary id: " << bID << std::endl;
        std::cout << "adjacentTetras: " << adjacentTetras[0] << "     " << adjacentTetras[1] << std::endl;

        exit(1);
    }

    adjacentTetras[nAdjoiningTetras] = tetID_;
    nAdjoiningTetras++;
}

inline int Triangle::getNeighbourTetra(int tetid)
{
    if(adjacentTetras[0] == tetid)
        return adjacentTetras[1];

    if(adjacentTetras[1] == tetid)
        return adjacentTetras[0];

    std::cout << "FATAL ERROR in getNeighborTetra!" << std::endl;
    std::cout << "Triangle ID: " << eID << std::endl;
    /* std::cout << nodes[0] << "  " << */
    /*              nodes[1] << "  " << */
    /*              nodes[2] << std::endl; */
    std::cout << adjacentTetras[0] << "  " <<
                 adjacentTetras[1] << std::endl;
    exit(1);
}

Tetrahedron::Tetrahedron(std::string line)
{
    int etype, ntags;
    std::istringstream iss(line);
    iss >> eID >> etype >> ntags;

    for (int i=0; i<4; i++)
    {
        tris.push_back(NULL);
        bID[i] = 0;
    }

    switch (etype)
    {
        case 4: nen = 4; break;
        case 11: nen = 10; break;
        default: std::cout<< "Tetrahedron element order unsupported: " << etype << std::endl;
    }

    iss >> pgID;

    //TODO: Fix this hack
    matID = pgID - 4;

    int dummy;
    for (int i = 1; i < ntags; i++)
        iss >> dummy;

    for(int i=0; i<nen; i++)
    {
        iss >> dummy;
        nodes.push_back(dummy);
    }

    /* faces.push_back({nodes[0], nodes[1], nodes[2]}); */
    /* faces.push_back({nodes[0], nodes[1], nodes[3]}); */
    /* faces.push_back({nodes[0], nodes[2], nodes[3]}); */
    /* faces.push_back({nodes[1], nodes[2], nodes[3]}); */


}

/* inline void Tetrahedron::setBoundaryID(Triangle * tri) */
inline void Tetrahedron::setBoundaryID(std::pair<std::vector<int>,Triangle*> tripair)
{

    std::vector<std::vector<int>> faces;
    // faces following xns ordering
    faces.push_back({this->nodes[0], this->nodes[1], this->nodes[2]});
    faces.push_back({this->nodes[0], this->nodes[1], this->nodes[3]});
    faces.push_back({this->nodes[1], this->nodes[2], this->nodes[3]});
    faces.push_back({this->nodes[0], this->nodes[2], this->nodes[3]});

    for(int f=0; f<4; f++)
    {
        std::sort(faces[f].begin(), faces[f].end());
        if (faces[f] == tripair.first)
        {
            bID[f] = tripair.second->bID;
            // std::cout << "Found Tet Face:" << bID[f] << std::endl;
            tris[f] = tripair.second;
        }
    }

        /* bool triIncludesNode[4]; */

        /* triIncludesNode[0] = (tri->nodes[0]==nodes[0]) || (tri->nodes[1]==nodes[0]) || (tri->nodes[2]==nodes[0]); */
        /* triIncludesNode[1] = (tri->nodes[0]==nodes[1]) || (tri->nodes[1]==nodes[1]) || (tri->nodes[2]==nodes[1]); */
        /* triIncludesNode[2] = (tri->nodes[0]==nodes[2]) || (tri->nodes[1]==nodes[2]) || (tri->nodes[2]==nodes[2]); */
        /* triIncludesNode[3] = (tri->nodes[0]==nodes[3]) || (tri->nodes[1]==nodes[3]) || (tri->nodes[2]==nodes[3]); */

             /* if(triIncludesNode[0] && triIncludesNode[1] && triIncludesNode[2]){ bID[0] = tri->bID;  tris[0] = tri; } */
        /* else if(triIncludesNode[0] && triIncludesNode[1] && triIncludesNode[3]){ bID[1] = tri->bID;  tris[1] = tri; } */
        /* else if(triIncludesNode[1] && triIncludesNode[2] && triIncludesNode[3]){ bID[2] = tri->bID;  tris[2] = tri; } */
        /* else if(triIncludesNode[0] && triIncludesNode[2] && triIncludesNode[3]){ bID[3] = tri->bID;  tris[3] = tri; } */


}


Triangle::~Triangle()
{

}

Tetrahedron::~Tetrahedron()
{

}


#endif /* ELEMENTS_H */

