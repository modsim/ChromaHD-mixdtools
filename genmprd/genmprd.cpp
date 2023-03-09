
#include "mixd.hpp"
#include "SimpleProgress.hpp"

#include <cmath>
#include <set>


class Node
{
private:
    int _id;

public:
    double _coords[3];

    static constexpr double eps = 1e-10;

    Node(int id, double *coords) : _id(id)
    {
        _coords[0] = coords[0];
        _coords[1] = coords[1];
        _coords[2] = coords[2];
    }

    inline int id() const
    {
        return _id;
    }

    inline void printCoords() const
    {
        std::cout <<
            this->_coords[0] << ", " <<
            this->_coords[1] << ", " <<
            this->_coords[2] << std::endl;
    }

    inline bool isCloseTo(const Node & n, int neglect_dir=-1) const
    {
        double delta = 0.0;

        for(int dir=0; dir<3; ++dir)
        {
            if(dir!=neglect_dir)
            {
                delta +=   (_coords[dir] - n._coords[dir])
                         * (_coords[dir] - n._coords[dir]);
            }
        }

        return ( sqrt(delta) < eps );
    }

};



bool operator< (const Node & n1, const Node & n2)
{
    return ( n1.id() < n2.id() );
}




static const int facenodes[4][3] = { {0, 1, 2},
                                     {0, 1, 3},
                                     {1, 2, 3},
                                     {0, 2, 3}  };


void printUsage(char *binaryName)
{
    std::cout << "MIXD genmprd tool (operates on mesh in directory of call)" << std::endl;
    std::cout << "Usage: " << binaryName << " <RNG A> <RNG B> <direction of periodicity: x/y/z> [-readmprd] [-readmtbl]" << std::endl;

    std::cout << "If (RNG_A==0 && RNG_B==0), nodes at min/max coords in direction of periodicity will matched." << std::endl;

    std::cout << "Use -readmprd to extend an existing mprd file with additional perdiodicity information." << std::endl;
    std::cout << "Use -readmtbl to correctly handle doubled interface nodes in the given mesh." << std::endl;
    std::cout << "Ensure that the input mesh is semi-discrete!!" << std::endl;
    std::cout << "Run the program on the semi-discrete mesh generated after gmsh2mixdv2 -d 4 -o 1 ..." << std::endl;
    std::cout << "Double the mprd file for space-time meshes. XNS should handle the +nnspace offset." << std::endl;

}


int main(int argc, char **argv)
{
    using namespace std;

    if(argc < 4)
    {
        printUsage(argv[0]);
        return 1;
    }

    const int rng_per_a = atoi(argv[1]);
    const int rng_per_b = atoi(argv[2]);
    const char per_dir_char = argv[3][0];

    int per_dir = -1;
    switch(per_dir_char)
    {
        case 'x': per_dir =  0; break;
        case 'y': per_dir =  1; break;
        case 'z': per_dir =  2; break;
        default:  per_dir = -1; break;
    }

    bool readmprd = false;
    bool readmtbl = false;


    cout << "MIXD genmprd tool..." << endl << endl;

    for (int i=4; i<argc; i++)
    {
        if(string(argv[i]) == "-readmprd")
        {
            readmprd = true;
            cout << "extending existing mprd file" << endl;
        }
        else if(string(argv[i]) == "-readmtbl")
        {
            readmtbl = true;
            cout << "using mtbl file" << endl;
        }
    }



    try{

    long ne, nn;
    mixd::readminf("./minf", &nn, &ne);



    const int nsd = 3;  // 3D only!!!
    const int nen = 4;  // tetrahedra only!!!
    const int nef = 4;  // tetrahedra only!!!
    const int nnf = 3;  // number of nodes per face... tetra only!


    cout << "Finding nodes of periodic surfaces... " << flush;
    set<Node> nodes_rng_a;
    set<Node> nodes_rng_b;


    if(rng_per_a==0 && rng_per_b==0)
    {
        mixd::MixdFile<double> mxyz("./mxyz", nn, nsd);
        mxyz.read();

        double min = mxyz(0,per_dir);
        double max = mxyz(0,per_dir);

        // find min/max extents
        for(int i=1; i<nn; i++)
        {
            if(mxyz(i,per_dir) < min) min = mxyz(i,per_dir);
            if(mxyz(i,per_dir) > max) max = mxyz(i,per_dir);
        }

        cout << "Model extents in direction " << per_dir << " are: [" << min << " | " << max << "]" << endl;

        for(int i=0; i<nn; i++)
        {
            if( abs(mxyz(i,per_dir)-min) < Node::eps )
                nodes_rng_a.insert( Node(i+1, &(mxyz(i))) );

            else if( abs(mxyz(i,per_dir)-max) < Node::eps )
                nodes_rng_b.insert( Node(i+1, &(mxyz(i))) );
        }
    }
    else
    {
        mixd::MixdFile<int>    mien("./mien", ne, nen);
        mixd::MixdFile<int>    mrng("./mrng", ne, nef);
        mixd::MixdFile<double> mxyz("./mxyz", nn, nsd);

        mien.read();
        mrng.read();
        mxyz.read();

        int nodeid;

        for(int ie=0; ie<ne; ie++)
        {
            for(int iface=0; iface<nef; iface++)
            {
                if(mrng(ie,iface)==rng_per_a)
                {
                    for(int fn=0; fn<nnf; fn++)
                    {
                        nodeid = mien(ie, facenodes[iface][fn]);

                        nodes_rng_a.insert( Node(nodeid, &(mxyz(nodeid-1))) );
                    }
                }

                else if(mrng(ie,iface)==rng_per_b)
                {
                    for(int fn=0; fn<nnf; fn++)
                    {
                        nodeid = mien(ie, facenodes[iface][fn]);

                        nodes_rng_b.insert( Node(nodeid, &(mxyz(nodeid-1))) );
                    }
                }
            }
        }
    }
    cout << "done!" << endl;

    /* for(set<Node>::const_iterator ita=nodes_rng_a.begin(); ita!=nodes_rng_a.end(); ++ita) */
    /*     std::cout << ita->id() << std::endl; */

    /* std::cout << "RNG B" << std::endl; */
    /* for(set<Node>::const_iterator itb=nodes_rng_b.begin(); itb!=nodes_rng_b.end(); ++itb) */
    /*     std::cout << itb->id() << std::endl; */
    /*     std::cout << itb->id() << ": "<< itb->_coords[0]  << ", " << itb->_coords[1] << ", " << itb->_coords[2] << std::endl; */


    cout << "RNG A contains " << nodes_rng_a.size() << " nodes." << endl;
    cout << "RNG B contains " << nodes_rng_b.size() << " nodes." << endl;

    cout << endl;

    if(nodes_rng_a.size() != nodes_rng_b.size())
    {
        cout << "Differing number of nodes cannot be matched! Exit!" << endl << endl;
        return 1;
    }


    mixd::MixdFile<int> mprd("./mprd", nn);

    if(readmprd)
        mprd.read();
    else
    {
        // initialize mprd
        for(int i=0; i<nn; i++)
            mprd(i) = i+1;  // init: self connectivity
    }

    cout << "Matching nodes of periodic RNGs (brute force)... " << endl;

    SimpleProgress sp(0, nodes_rng_a.size(), 10);

    //DONE: Read mtbl
    //DONE: Calculate original number of nodes before interface doubling
    //      nn_original = mtbl.size() - n_changed / 2
    mixd::MixdFile<int> mtbl("./mtbl", nn);

    if(readmtbl)
    {
        mtbl.read();
    }
    else
    {
        // initialize mtbl
        for(int i=0; i<nn; i++)
            mtbl(i) = i+1;  // init: self connectivity
    }

    int n_doubled = 0;
    for (int i=0; i<nn; i++)
    {
        if (i+1 != mtbl(i))
            n_doubled++;
    }
    n_doubled /= 2;

    std::cout << "Number of doubled nodes: " << n_doubled << std::endl;

    int nn_orig = nn - n_doubled;
    std::cout << "Original number of nodes (before doubling bead surface): " << nn_orig << std::endl;

    int nodes_matched = 0;
    for(set<Node>::const_iterator ita=nodes_rng_a.begin(); ita!=nodes_rng_a.end(); ++ita)
    {
        if(mprd(ita->id() -1) == ita->id())
        {
            // if not doubled (in mtbl)
            if (mtbl(ita->id() - 1) == ita->id())
            {
                for(set<Node>::const_iterator itb=nodes_rng_b.begin(); itb!=nodes_rng_b.end(); ++itb)
                {
                    if(ita->isCloseTo(*itb, per_dir))
                    {
                        nodes_matched++;

                        mprd(ita->id() -1) = itb->id();
                        mprd(itb->id() -1) = ita->id();

                        break;
                    }
                }
            }
            else // if node on rng_a is doubled
            {
                for(set<Node>::const_iterator itb=nodes_rng_b.begin(); itb!=nodes_rng_b.end(); ++itb)
                {
                    if(ita->isCloseTo(*itb, per_dir))
                    {
                        if ( ((ita->id() <= nn_orig) && (itb->id() <= nn_orig)) || ((ita->id() > nn_orig) && (itb->id() > nn_orig)) )
                        {
                            nodes_matched++;

                            // Normal
                            mprd(ita->id() -1) = itb->id();
                            mprd(itb->id() -1) = ita->id();

                            // TODO: Figure this out properly.
                            // Because we only modify mprds if mprd(id-1) == id
                            // we don't need to increment the counter for the duplicates.
                            // we can calculate the ID from MTBL and silently update it.
                            // The counter update happens in the 'else' branch anyway

                            // Doubled
                            mprd(mtbl(ita->id() -1) -1) = mtbl(itb->id() -1);
                            mprd(mtbl(itb->id() -1) -1) = mtbl(ita->id() -1);

                            break;

                        }

                    }
                }

            }
        }
        else
            nodes_matched++;

        sp.printIfHitNext(nodes_matched);
    }
    cout << "done!" << endl;


    cout << "Nodes matched: " << nodes_matched << endl;

    if((unsigned)nodes_matched != nodes_rng_a.size())
    {
        cout << "Unknown error: Not all nodes could be matched! Exit!" << endl << endl;
        return 1;
    }

    cout << endl;



    mprd.write();

    } catch(mixd::MixdException e)
    { std::cout << e.msg() << std::endl; }


    return 0;
}
