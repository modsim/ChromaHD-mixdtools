#ifndef WRITEMIXD_H
#define WRITEMIXD_H

#include <iostream>
#include <iomanip>
#include <fstream>

/* #include "endian/include/boost/endian/conversion.hpp" */
#include <boost/endian/conversion.hpp>
#include <boost/endian/buffers.hpp>
/* #include "endian/include/boost/endian/buffers.hpp" */

#include "elements.h"
#include "node.h"
#include "mesh.h"

//TODO: make sure endianness is ALWAYS big instead of reversing it.

bool isBigEndian()
{
    short word = 0x4321;
    if ((* (char*) & word) != 0x21) return true;
    else      return false;
}

void swapbytes(char *array, int nelem, int elsize)
{
    int sizet, sizem, i, j;
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

void writeMINF(Mesh &mesh)
{
    std::cout << "writing minf... ";
    std::ofstream of;
    of.open("minf");

    int ne = mesh.tets.size();
    int nn = mesh.nodes.size() + mesh.dblBoundaryNodes.size();

    of << "ne" << std::setw(12) << ne << std::endl;

    /* if(!generate_st) */
        /* of << "nn" << std::setw(12) <<   (n_nodes+dblBoundaryNodes.size()) << std::endl; */
        of << "nn" << std::setw(12) << nn << std::endl;
    /* else */
        /* of << "nn" << std::setw(12) << 2*(n_nodes+dblBoundaryNodes.size()) << std::endl; */

    of.close();

    std::cout << "done!" << std::endl;

}


void writeMIEN(Mesh &mesh)
{
    std::cout << "writing mien... " << std::flush;
    std::ofstream mien;
    mien.open("mien", std::ios::out | std::ios::binary);

    for (typename std::vector<Tetrahedron*>::iterator it_tet=mesh.tets.begin(); it_tet!=mesh.tets.end(); ++it_tet)
    {
        int nen = (*it_tet)->nen;

        /* char *buf = (char*) malloc(4); */
        for (int i=0; i<nen; i++)
        {
            int tmp = boost::endian::native_to_big((*it_tet)->nodes[i]);
            mien.write(reinterpret_cast<char*>(&tmp), sizeof(tmp));

            /* memcpy(buf, &((*it_tet)->nodes[i]), 4); */

            /* if(!isBigEndian())  // swap bytes if neccessary */
            /* { */
            /*     swapbytes(buf, 1, 4); */
            /* } */

            /* mien.write(buf, sizeof(buf)); */

        }
    }
    mien.close();

    std::cout << "done!" << std::endl;

}

void writeMMAT(Mesh &mesh)
{
    std::cout << "writing mmat... " << std::flush;
    std::fstream mmat;
    mmat.open("mmat", std::ios::out | std::ios::binary);


    for (typename std::vector<Tetrahedron*>::iterator it_tet=mesh.tets.begin(); it_tet!=mesh.tets.end(); ++it_tet)
    {
        int matID = (*it_tet)->matID;
        int tmp = boost::endian::native_to_big(matID);
        mmat.write(reinterpret_cast<char*>(&tmp), sizeof(tmp));
    }

    mmat.close();


    std::cout << "done!" << std::endl;

}

void writeMXYZ(Mesh &mesh)
{
    std::cout << "writing mxyz... " << std::flush;
    std::fstream mxyz;
    mxyz.open("mxyz", std::ios::out | std::ios::binary);

    char *nx_buf = (char*) malloc(8);    // node x coord buffer --> 8 byte double
    char *ny_buf = (char*) malloc(8);    // node y coord buffer --> 8 byte double
    char *nz_buf = (char*) malloc(8);    // node z coord buffer --> 8 byte double

    //TODO: if ST, loop twice
    for (std::vector<Node *>::iterator it_nodes = mesh.nodes.begin(); it_nodes != mesh.nodes.end(); it_nodes++)
    {

        memcpy(nx_buf, &((*it_nodes)->x), 8);
        memcpy(ny_buf, &((*it_nodes)->y), 8);
        memcpy(nz_buf, &((*it_nodes)->z), 8);

        if(!isBigEndian())  // swap bytes if neccessary
        {
            swapbytes(nx_buf, 1, 8);
            swapbytes(ny_buf, 1, 8);
            swapbytes(nz_buf, 1, 8);
        }

        mxyz.write(nx_buf, sizeof(nx_buf));
        mxyz.write(ny_buf, sizeof(ny_buf));
        mxyz.write(nz_buf, sizeof(nz_buf));

    }



    for(std::set<int>::iterator it=mesh.dblBoundaryNodes.begin(); it!=mesh.dblBoundaryNodes.end(); ++it)
    {
        //TODO: check how to deal with doubled nodes.
        memcpy(nx_buf, &(mesh.nodes[(*it)-1]->x), 8);
        memcpy(ny_buf, &(mesh.nodes[(*it)-1]->y), 8);
        memcpy(nz_buf, &(mesh.nodes[(*it)-1]->z), 8);

        if(!isBigEndian())  // swap bytes if neccessary
        {
            swapbytes(nx_buf, 1, 8);
            swapbytes(ny_buf, 1, 8);
            swapbytes(nz_buf, 1, 8);
        }

        mxyz.write(nx_buf, sizeof(nx_buf));
        mxyz.write(ny_buf, sizeof(ny_buf));
        mxyz.write(nz_buf, sizeof(nz_buf));
    }

    mxyz.close();
    std::cout << "done!" << std::endl;

}

void writeMTBL(Mesh &mesh)
{
    std::cout << "writing mtbl... " << std::flush;
    std::fstream mtbl;
    mtbl.open("mtbl", std::ios::out | std::ios::binary);

    //TODO:if(entity_dbl != -1)
    //int nloops = (generate_st)? 2 : 1;
    int loop = 0;

    int partner_node = 0;
    int n_nodes = mesh.nodes.size();

    for(int n=0; n<n_nodes; n++)
    {

        if(mesh.dblBoundaryNodes.count(n+1))
        {
            partner_node = mesh.fromSingleToDoubleNodeIDs.find(n+1)->second;
        }
        else
            partner_node = n+1;

        partner_node += loop*(n_nodes + mesh.dblBoundaryNodes.size());

        int tmp = boost::endian::native_to_big(partner_node);
        mtbl.write(reinterpret_cast<char*>(&tmp), sizeof(tmp));

    }

    for(std::set<int>::iterator it=mesh.dblBoundaryNodes.begin(); it!=mesh.dblBoundaryNodes.end(); ++it)
    {
        partner_node = *it + loop*(n_nodes + mesh.dblBoundaryNodes.size());
        int tmp = boost::endian::native_to_big(partner_node);
        mtbl.write(reinterpret_cast<char*>(&tmp), sizeof(tmp));
    }

    mtbl.close();
    std::cout << "done!" << std::endl;

}

void writeMTBLDUAL(Mesh &mesh)
{
    std::cout << "writing mtbldual... " << std::flush;
    std::fstream mtbldual;
    mtbldual.open("mtbl.dual", std::ios::out | std::ios::binary);

    for (typename std::vector<Tetrahedron*>::iterator it_tet=mesh.tets.begin(); it_tet!=mesh.tets.end(); ++it_tet)
    {

        /* int nbtet[4]; */
        /* char *ef_buf = (char*) malloc(4*4);  // buffer for element face lists --> 4 ints of 4 bytes each */
        int nbtet;
        for(int iface=0; iface<4; iface++)
        {
            if((*it_tet)->tris[iface] != NULL)
            {
                nbtet = - ((*it_tet)->tris[iface]->getNeighbourTetra((it_tet)- mesh.tets.begin() + 1));
            }
            else
                nbtet = 0;

            int tmp = boost::endian::native_to_big(nbtet);
            mtbldual.write(reinterpret_cast<char*>(&tmp), sizeof(tmp));
        }


        /* memcpy(ef_buf, nbtet, 4*4); */
        /* if(!isBigEndian())  // swap bytes if neccessary */
        /* { */
        /*     swapbytes(ef_buf, 4, 4); */
        /* } */
        /* mtbldual.write(ef_buf,sizeof(ef_buf)); */

    }

    mtbldual.close();

    std::cout << "done!" << std::endl;
}


void writeMRNG(Mesh &mesh)
{
    std::cout << "writing mrng... " << std::flush;
    std::fstream mrng;
    mrng.open("mrng", std::ios::out | std::ios::binary);

    for (typename std::vector<Tetrahedron *>::iterator it_tet=mesh.tets.begin(); it_tet!=mesh.tets.end(); ++it_tet)
    {
        for (int i=0; i<4 ; i++)
        {
            int bID = (*it_tet)->bID[i];
            int tmp = boost::endian::native_to_big(bID);
            mrng.write(reinterpret_cast<char*>(&tmp), sizeof(tmp));
        }
    }

    mrng.close();


    std::cout << "done!" << std::endl;

}

void writeMIXD(Mesh &mesh)
{
    writeMINF(mesh);
    writeMIEN(mesh);
    writeMMAT(mesh);
    writeMXYZ(mesh);
    writeMTBL(mesh);
    writeMTBLDUAL(mesh);
    writeMRNG(mesh);
}

#endif
