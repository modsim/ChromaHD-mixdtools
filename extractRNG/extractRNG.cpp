#include "mixd.hpp"
#include <vector>
#include <cmath>
#include <set>
#include <getopt.h>
#include "SimpleProgress.hpp"

class Triangle
{
    public:
        int nen=3;

        double _nodes[3];

        Triangle(int n1, int n2, int n3)
        {
            _nodes[0] = n1;
            _nodes[1] = n2;
            _nodes[2] = n3;
        }

};

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

void getArgs(int argc, char * argv[], int& ndf, int& nts, int& rng, int& spacetimeupper, int& jump, std::string& meshdir)
{
    if (argc < 2)
    {
        std::cout << "extractRNG: Extract RNG mesh and data from an existing MIXD mesh + data." <<std::endl;
        std::cout << "Run this program in the solution directory." << std::endl;
        std::cout << "Usage: ./extractRNG <data_file> -n/--ndf <ndf> -t/--nts <nts> -r/--rng <rng> -s/--spacetime-upper" << std::endl;
        std::cout << "Writes data from lower slab if spacetime-upper (-s) is not specified." << std::endl;
        std::cout << "Copied from stitchperiodic" << std::endl;
        exit(-1);
    }

    int c;

    while (1)
    {
        static struct option long_options[] =
        {
            /* These options set a flag. */
            /* {"verbose", no_argument,       &verbose_flag, 1}, */
            /* {"brief",   no_argument,       &verbose_flag, 0}, */
            /* {"spacetime-upper",  no_argument, &spacetimeupper, 1}, */
            /* {"s",  no_argument, &spacetimeupper, 1}, */
            /* These options don’t set a flag.
               We distinguish them by their indices. */
            {"ndf",  required_argument, 0, 'n'},
            {"nts",  required_argument, 0, 't'},
            {"rng",  required_argument, 0, 'r'},
            {"meshdir",  required_argument, 0, 'm'},
            {"spacetime-upper",  no_argument, 0, 's'},
            {"jump",  no_argument, 0, 'j'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "n:t:r:sj:m:", long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;

        switch (c)
        {
            case 0:
                /* If this option set a flag, do nothing else now. */
                if (long_options[option_index].flag != 0)
                    break;
                printf ("option %s", long_options[option_index].name);
                if (optarg)
                    printf (" with arg %s", optarg);
                printf ("\n");
                break;

            case 'n': ndf = std::atoi(optarg); break;
            case 't': nts = std::atoi(optarg); break;
            case 'r': rng = std::atoi(optarg); break;
            case 'j': jump= std::atoi(optarg); break;
            case 's': spacetimeupper = 1; break;
            case 'm': meshdir=optarg; break;
            case '?':
                /* getopt_long already printed an error message. */
                break;
            default: abort ();
        }
    }

    /* if (optind >= argc) */
    /* { */
    /*     std::cout << "No data file provided!" << std::endl; */
    /*     exit(-1); */
    /* } */

}


int main(int argc, char * argv[])
{
    int ndf=0;
    int nts=1;
    int rng=0;
    int spacetimeupper=0;
    int jump=0; // To jump over timesteps from data.all/data.in. Helps with flow bdf2 output, set j=1
    std::string meshdir="../mesh";

    getArgs(argc, argv, ndf, nts, rng, spacetimeupper, jump, meshdir);

    std::cout << "ndf: " << ndf << std::endl;
    std::cout << "nts: " << nts << std::endl;
    std::cout << "rng: " << rng << std::endl;
    std::cout << "spacetime upper: " << spacetimeupper << std::endl;
    std::cout << "meshdir: " << meshdir << std::endl;
    /* std::cout << "data: " << argv[optind] << std::endl; */

    long ne, nn;
    mixd::readminf(meshdir + "/minf", &nn, &ne);

    if ((spacetimeupper == 1) && (nn % 2 != 0)) exit(-1); //exit if not properly spacetime

    const int nsd = 3;  // 3D only!!!
    const int nen = 4;  // tetrahedra only!!!
    const int nef = 4;  // tetrahedra only!!!
    const int nnf = 3;  // number of nodes per face... tetra only!

    // face ordering (according to MIXD format):
    // face 0 consists of nodes 0 1 2
    // face 1 consists of nodes 0 1 3
    // face 2 consists of nodes 1 2 3
    // face 3 consists of nodes 0 2 3
    const std::vector<std::vector<int>> facemap = {
        {0,1,2},
        {0,1,3},
        {1,2,3},
        {0,2,3}
    };

    mixd::MixdFile<int>    mien(meshdir + "/mien", ne, nen);
    mixd::MixdFile<int>    mrng(meshdir + "/mrng", ne, nef);
    mixd::MixdFile<double> mxyz(meshdir + "/mxyz", nn, nsd);

    mien.read();
    mrng.read();
    mxyz.read();

    std::set<Node> nodes_rng;
    std::vector<Triangle> tris;
    int n1 = -1; int n2 = -1; int n3 = -1;
    int newid;

    // Iterate over elements::faces
    // and save node and triangle information
    for(int ie=0; ie<ne; ie++)
    {
        for(int iface=0; iface<nef; iface++)
        {
            // Since mrng files are not changed for spacetime meshes,
            // We always get node information for the lower slab
            // Updates in data are handled later with an offset
            if(mrng(ie,iface)==rng)
            {
                std::vector<int> nodeids;
                for(int fn=0; fn<nnf; fn++)
                {
                    int nodeid = mien(ie, facemap[iface][fn]);
                    nodes_rng.insert( Node(nodeid, &(mxyz(nodeid-1))) );
                    nodeids.push_back(nodeid);
                }
                tris.push_back(Triangle(nodeids[0], nodeids[1], nodeids[2]));
            }
        }
    }

    // for every triangle(i)::node(j)
    // renumber nodes by position in vector
    std::vector<Node> vNodes(nodes_rng.begin(), nodes_rng.end());
    for (int i=0; i< tris.size(); i++)
        for(int j=0; j<3; j++)
        {
           int dummy = tris[i]._nodes[j];
           auto ret = std::find_if(vNodes.begin(), vNodes.end(), [dummy](const Node& n){return n.id() == dummy;});
           tris[i]._nodes[j] = ret - vNodes.begin() + 1;
        }

    std::cout << "Number of RNG elements: " << tris.size() << std::endl;
    std::cout << "Number of RNG nodes: " << vNodes.size() << std::endl;
    std::cout << "Number of RNG nodes: " << nodes_rng.size() << std::endl;

    mixd::MixdFile<int> rngmien("rng.mien", tris.size(), 3, false);
    mixd::MixdFile<double> rngmxyz("rng.mxyz", vNodes.size(), 3, false);



    // Prepare mien file
    for (int i=0; i<tris.size(); i++)
        for(int j=0; j<3; j++)
            rngmien(i,j) = tris.at(i)._nodes[j];

    std::cout << "rngmien ready..." << std::endl;

    // Prepare mxyz file
    for (int i=0; i<vNodes.size(); i++)
        for(int j=0;j<3; j++)
            rngmxyz(i,j) = vNodes.at(i)._coords[j];

    std::cout << "rngmxyz ready..." << std::endl;

    rngmien.write();
    rngmxyz.write();

    int ne_rng = tris.size();
    int nn_rng = vNodes.size();

    std::cout << "writing minf... ";
    std::ofstream of;
    of.open("rng.minf");
    of << "ne" << std::setw(12) << ne_rng << std::endl;
    of << "nn" << std::setw(12) << nn_rng << std::endl;
    of.close();

    mixd::MixdFile<double> data(argv[optind], nn, ndf,false);

    // remove previous rng.data file
    // Used because I append data instead of writing for now
    std::remove("rng.data");
    mixd::MixdFile<double> rngdata("rng.data", nodes_rng.size(), ndf, false);

    long offset=0;
    if (spacetimeupper == 1)
        offset=nn/2;

    SimpleProgress sp(0, nts, 20);
    for(int i=0; i<nts; i++)
    {
        data.read(jump+i);
        int nodes_count = 0;
        for(auto it:nodes_rng)
        {
            for(int idf=0; idf<ndf; idf++)
            {
                // if spacetimeupper, get data from upper slab
                // works because nodes_rng only captures lower slab nodes
                rngdata(nodes_count, idf) = data(offset + it.id()-1, idf);
            }
            nodes_count++;
        }
        rngdata.append();
        /* datarng.write(i); // NOTE: write(i) doesn't write properly... zeroes previous data */
        sp.printIfHitNext(i);

    }
    std::cout << "rngdata ready..." << std::endl;


    std::cout << "done!" << std::endl;
}
