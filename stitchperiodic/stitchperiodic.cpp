#include "mixd.hpp"
#include <vector>
#include <cmath>
#include <set>
#include <getopt.h>
#include "SimpleProgress.hpp"

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

void getArgs(int argc, char * argv[], int& ndf, int& nts, int& rng, int& spacetimeupper, int& jump, bool& dryRun)
{
    if (argc < 2)
    {
        std::cout << "stitchperiodic: Stitch periodically linked simulations together." <<std::endl;
        std::cout << "Run this program in the solution directory." << std::endl;
        std::cout << "Usage: ./stitchperiodic <data_file> -n/--ndf <ndf> -t/--nts <nts> -r/--rng <rng> -s/--spacetime-upper" << std::endl;
        std::cout << "      > FLOW: ./stitchperiodic data.out -r 2 -n 4 -t 2" << std::endl;
        std::cout << "      > MASS: ./stitchperiodic data.all -r 2 -n 2 -t <nts>" << std::endl;
        std::cout << "Generates rng.xyz and rng.data files that must be copied to the  new solution directory." << std::endl;
        std::cout << "Writes data from lower slab if spacetime-upper (-s) is not specified." << std::endl;
        std::cout << "Ensure that xns.in is updated accordingly." << std::endl;
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
            {"spacetime-upper",  no_argument, 0, 's'},
            {"jump",  required_argument, 0, 'j'},
            {"dry-run", no_argument, 0, 'd'},
            {0, 0, 0, 0}
        };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "n:t:r:sj:d", long_options, &option_index);

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
            case 'd': dryRun = true; break;
            case '?':
                /* getopt_long already printed an error message. */
                break;
            default: abort ();
        }
    }

    if ((optind >= argc) && !dryRun)
    {
        std::cout << "No data file provided!" << std::endl;
        exit(-1);
    }

}


int main(int argc, char * argv[])
{
    int ndf=0;
    int nts=1;
    int rng=0;
    int spacetimeupper=0;
    int jump=0; // To jump over timesteps from data.all/data.in. Helps with flow bdf2 output, set j=1
    bool dryRun=false;

    getArgs(argc, argv, ndf, nts, rng, spacetimeupper, jump, dryRun);

    std::cout << "ndf: " << ndf << std::endl;
    std::cout << "nts: " << nts << std::endl;
    std::cout << "rng: " << rng << std::endl;
    std::cout << "spacetime upper: " << spacetimeupper << std::endl;
    std::cout << "dry run: " << dryRun << std::endl;
    /* std::cout << "data: " << argv[optind] << std::endl; */

    long ne, nn;
    mixd::readminf("../mesh/minf", &nn, &ne);

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

    mixd::MixdFile<int>    mien("../mesh/mien", ne, nen);
    mixd::MixdFile<int>    mrng("../mesh/mrng", ne, nef);
    mixd::MixdFile<double> mxyz("../mesh/mxyz", nn, nsd);

    mien.read();
    mrng.read();
    mxyz.read();

    std::set<Node> nodes_rng;

    for(int ie=0; ie<ne; ie++)
    {
        for(int iface=0; iface<nef; iface++)
        {
            if(mrng(ie,iface)==rng)
            {
                for(int fn=0; fn<nnf; fn++)
                {
                    int nodeid = mien(ie, facemap[iface][fn]);
                    nodes_rng.insert( Node(nodeid, &(mxyz(nodeid-1))) );
                }
            }
        }
    }

    std::cout << "Number of RNG nodes (nn_rng): " << nodes_rng.size() << std::endl;

    if (dryRun) exit(0);

    /* mixd::MixdFile<int>    rngnodes("rng.nodeids", nodes_rng.size(), 1, false); */
    /* int index=0; */
    /* for(auto it:nodes_rng) */
    /* { */
    /*     rngnodes(index,1) = it.id(); */
    /*     index++; */
    /* } */
    /* rngnodes.write(); */
    /* exit(-1); */

    // remove previous rng.data file
    // Used because I append data instead of writing for now
    std::remove("rng.data");
    mixd::MixdFile<double> datarng("rng.data", nodes_rng.size(), ndf, false);
    mixd::MixdFile<double> rngxyz("rng.xyz", nodes_rng.size(), nsd, false);
    mixd::MixdFile<double> data(argv[optind], nn, ndf,false);

    int count =0;
    for (auto it:nodes_rng)
    {
        for(int isd=0; isd<nsd; isd++)
            rngxyz(count,isd) = it._coords[isd];
        count++;
    }
    rngxyz.write();

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
                datarng(nodes_count, idf) = data(offset + it.id()-1, idf);
            }
            nodes_count++;
        }
        datarng.append();
        /* datarng.write(i); // TODO: write(i) doesn't write properly... zeroes previous data */
        sp.printIfHitNext(i);

    }

    std::cout << std::endl;
    std::cout << "Number of RNG nodes (nn_rng): " << nodes_rng.size() << std::endl;

}
