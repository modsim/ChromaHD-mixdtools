#include "mixd.hpp"
#include "dunavant.hpp"
#include <cmath>

// TODO: make integrate() operate on only a given dof. We don't need to integrate every dof
// TODO: Allow specifying mesh directory instead of assuming current dir

static const int facemap_tri[3][2] = { {0,1}, {1,2}, {2,0} };
static const int facemap_tet[4][3] = { {0,1,2}, {0,3,1}, {1,3,2}, {0,2,3} };
static const int facemap_qua[4][2] = { {0,1}, {1,2}, {2,3}, {3,0} };
static const int facemap_hex[6][4] = { {0,1,2,3}, {0,1,5,4}, {1,2,6,5}, {2,3,7,6}, {0,4,7,3}, {4,5,6,7} };

template <int NEN, int NEF, int NFN>
int node(int f, int fn);

template <> int node<3,3,2>(int f, int fn)
{
    return facemap_tri[f][fn];
}

template <> int node<4,4,3>(int f, int fn)
{
    return facemap_tet[f][fn];
}

template <> int node<4,4,2>(int f, int fn)
{
    return facemap_qua[f][fn];
}

template <> int node<8,6,4>(int f, int fn)
{
    return facemap_hex[f][fn];
}




double getTime(const std::string & timefile, int tstep)
{
    if(timefile != "")
    {
        using namespace std;
        using namespace mixd;
        
        ifstream file(timefile.c_str(), ifstream::in);
        
        if(!file.is_open())
            throw MixdException("could not open file " + timefile);
        
        string line;
        
        for(int i=0; i<tstep; i++)
            getline(file, line);
        
        getline(file, line);
        
        file.close();
        
        return atof(line.c_str());
    }
    else
        return (double)tstep;
}



void printUsage(char *binaryName)
{
    std::cout << "MIXD tool for computing chromatograms" << std::endl;
    std::cout << "Usage: " << binaryName << " <data file> <flowfield file> <time file> <chromatogram rng> <ndf> <nsteps>" << std::endl;
}


void integrate
    (const mixd::MixdFile<int> & mien, const mixd::MixdFile<double> & mxyz, const mixd::MixdFile<double> & data, double * intvals);
    


int main(int argc, char **argv)
{
    if(argc != 7)
    {
        printUsage(argv[0]);
        return 1;
    }
    
    try{
    
    std::cout << "MIXD chromatogram tool ONLY FOR TETRA MESHES..." << std::endl << std::endl;
    
    using namespace std;
    using namespace mixd;
    
    const string datafile(argv[1]);
    const string flowfile(argv[2]);
    const string timefile(argv[3]);
    const int chromrng = atoi(argv[4]);
    const int ndf = atoi(argv[5]);
    const int nsteps = atoi(argv[6]);
    
    cout << "chromatogram rng: " << chromrng << endl << endl;
    
    const int nen = 4;
    const int nef = 4;
    const int nsd = 3;
    const int nfn = 3;
    
    long ne, nn;
    readminf("./minf", &nn, &ne);
    
    
    MixdFile<int> * mien = new MixdFile<int>("./mien", ne, nen);
    MixdFile<int> * mrng = new MixdFile<int>("./mrng", ne, nef);
    
    mien->read();
    mrng->read();
    
    
    
    //===============================================================
    // determine number of surface elements
    int nelems = 0;
    
    for(int i=0; i<ne; i++)
        for(int j=0; j<nef; j++)
            if((*mrng)(i,j) == chromrng)
            {
                // DOES NOT WORK FOR INTERIOR BOUNDARIES!!!!
                nelems++;
            }
    
    cout << "Number of surface elements: " << nelems << endl<< endl;
    MixdFile<int> mien_surf("./surf/mien", nelems, nfn);
    //===============================================================
    
    
    
    int * newnodeid = new int[nn];
    for(int i=0; i<nn; i++)
        newnodeid[i] = 0;
    
    int nnodes = 0;
    
    int surfelem = 0;
    
    for(int i=0; i<ne; i++)
    {
        for(int j=0; j<nef; j++)
        {
            if((*mrng)(i,j) == chromrng)
            {
                for(int k=0; k<nfn; k++)
                {
                    // global node id
                    // CAUTION!!!! here we apply immediately the space-time shift of nn/2
                    // because we want to deal only with the upper level!!!!
                    // THIS IS WORSE THAN QUICK-AND-DIRTY!!!!!!!!!
                    int globalnode = (*mien)(i, node<4,4,3>(j, k)) - 1 + nn/2;
                    
                    if(newnodeid[globalnode] == 0)
                    {
                        nnodes++;
                        newnodeid[globalnode] = nnodes;
                    }
                    
                    mien_surf(surfelem, k) = newnodeid[globalnode];
                }
                
                surfelem++;
            }
        }
    }
    
    cout << "Number of surface nodes: " << nnodes << endl<< endl;
    
    delete mien;
    delete mrng;
    
    
    MixdFile<double> * mxyz = new MixdFile<double>("./mxyz", nn, nsd);
    MixdFile<double> * flow = new MixdFile<double>(flowfile, nn, nsd+1);
    
    mxyz->read();
    flow->read();
    
    MixdFile<double> mxyz_surf("./surf/mxyz", nnodes, nsd);
    MixdFile<double> flow_surf("./surf/flow", nnodes, nsd);
    
    for(int i=0; i<nn; i++)
    {
        if(newnodeid[i] > 0)
        {
            for(int j=0; j<nsd; j++)
            {
//                 if(j<nsd-1)
                mxyz_surf(newnodeid[i]-1, j) = (*mxyz)(i,j);
                flow_surf(newnodeid[i]-1, j) = (*flow)(i,j);
            }
        }
    }
    
    
    delete mxyz;
    delete flow;
    
//     writeminf("./surf/minf", nnodes, nelems);
//     mien_surf.write();
//     mxyz_surf.write();
//     flow_surf.write();
    
    
//     // CAUTION!!!
//     int ndf = 1;
    
    MixdFile<double> data_surf("./surf/data", nnodes, ndf);
    
    
    
    double * chromatogram = new double[ndf*nsteps];

    double * intflow = new double[4];
    integrate(mien_surf, mxyz_surf, flow_surf, intflow);

    cout << "\nFlowrates: ";
    for (int i=0; i<4; i++)
    {
        cout << intflow[i] << "\t";
    }
    cout << endl << endl;;
    
    
//     // compute max breakthrough flux value
//     for(int i=0; i<nnodes; i++)
//     {
//         for(int j=0; j<ndf; j++)
//             data_surf(i, j) = maxconc * flow_surf(i, 2);
//     }
//     
//     double maxflux;
//     integrate(mien_surf, mxyz_surf, data_surf, maxflux);
    
    
    cout << "Analyzing " << nsteps << " time steps..." << flush << endl << endl;
    
    for(int step=0; step<nsteps; step++)
    {
        cout << setw(5) << step << flush;
        
        MixdFile<double> data(datafile, nn, ndf, false);
        data.read(step);
        
        // copy data to condensed surface node data array
        // maybe iterate only over upper space-time level???
        for(int i=0; i<nn; i++)
        {
            if(newnodeid[i] > 0)
            {
                for(int j=0; j<ndf; j++)
                    // ALSO VERY DANGEROUS!!!!
                    // ASSUMES THAT THE OUTLET IS A Z-PLANE!!!!!!!!!
                    // SHOULD USE A REAL INNER PRODUCT OF VELOCITY VECTOR AND OUTER NORMAL!!!
                    data_surf(newnodeid[i]-1, j) = data(i,j) * flow_surf(newnodeid[i]-1, 2);
            }
        }
        
        // NOTE: chromatogram here stores integral(c * u * dA) = mass flow rate
        integrate(mien_surf, mxyz_surf, data_surf, &(chromatogram[step*ndf]));
    }
    
    cout << endl << endl;
    
    
    // print to files...
    cout << "Generating output text file... " << flush;
    
    ofstream cfile("chromatogram_rng" + str(chromrng) + ".csv", ifstream::out);
    
    if(!cfile.is_open())
        throw MixdException("could not open file");
    
    for(int step=0; step<nsteps; step++)
    {
        double time = getTime(timefile, step);
        
        cfile << scientific << setprecision(10) << time ;
        
        // // NOTE: We write out mass-flowrate / average-volume-flowrate to get back the concentration at the exit
        // for(int idf=0; idf<ndf; idf++)
        //     cfile << ',' << chromatogram[step*ndf+idf] / intflow[2];

        // only the first dof is sufficient
        cfile << ',' << chromatogram[step*ndf+0] / intflow[2];
        
        cfile << endl;
    }
    
    cfile.close();
    
    cout << "done!" << endl;
    
    delete[] chromatogram;
    delete[] newnodeid;
    
    
    
    
    } catch(mixd::MixdException e)
    { std::cout << e.msg() << std::endl; }
    
    
    return 0;
}





    
    // Given is an arbitrary triangle with nodes (x1,y1), (x2,y2), (x3,y3).
    // These are stored in the array x in the following order: (x1,y1, x2,y2, x3,y3)
    // On each node we have given nv data values.
    // These are stored in the array v in the following order (for example nv=2): (v11,v12, v21,v22, v31,v32)
    // Then this function computes for each value the unique coefficients a,b,c for linear interpolation on this triangle.
    // v = a*x + b*y + c
    // They are stored in the array c in the following order: (a1,b1,c1, a2,b2,c2)
    inline void linCoeffs(const double *x, const double *v, int nv, double *c)
    {
        double denom = (x[0]-x[4]) * (x[3]-x[5]) - (x[2]-x[4]) * (x[1]-x[5]);
        
        for(int i=0; i<nv; i++)
        {
            c[i*3 + 0] = ((v[0*nv + i]-v[2*nv + i]) * (x[3]-x[5])-(v[1*nv + i]-v[2*nv + i]) * (x[1]-x[5])) / denom;
            c[i*3 + 1] = ((v[1*nv + i]-v[2*nv + i]) * (x[0]-x[4])-(v[0*nv + i]-v[2*nv + i]) * (x[2]-x[4])) / denom;
            c[i*3 + 2] = ( v[2*nv + i] - x[4]*c[i*3+0] - x[5]*c[i*3+1] );
        }
    }





// compute and store average concentrations for given time step
    // parameter spacetime may be:
    // 0 for semi-discrete mesh, 1 for lower space-time level, 2 for upper space-time level
    // CAUTION: this function works ONLY for tetrahedra meshes!!!!
    void integrate
    (const mixd::MixdFile<int> & mien, const mixd::MixdFile<double> & mxyz, const mixd::MixdFile<double> & data, double * intvals)
    {
        const int nen = mien.cols();  // always 4
        const int nsd = 2; //mxyz.cols();  // always 3
        const int ndf = data.cols();
        
//         // if we are using the upper space-time level,
//         // always skip the first nn/2 nodes
//         int nnskip = 0;
//         if(spacetime == 2)
//         {
//             int nn = mxyz.rows();
//             if(nn % 2 != 0)
//                 throw mixd::MixdException("space-time mesh must have even number of nodes!");
//             nnskip = nn / 2;
//         }
        
        // quadrature rule type
        // (first order is sufficient since we use only linear shape functions)
        int rule = 1;
        
        // quadrature rule order
        int order = dunavant_order_num(rule);
        
        // quad point locations and weights
        double *xqr = new double[2*order]; // reference coordinates
        double *xqp = new double[2*order]; // physical coordinates
        double *wq  = new double[  order];
        
        // initialize (reference) quad point locations and weights
        dunavant_rule(rule, order, xqr, wq);
        
        // element local values
        double *x = new double[nen * nsd];  // nodal coordinates
        double *v = new double[nen * ndf];  // nodal DOFs
        double *c = new double[  3 * ndf];  // coefficients of linear distribution
        
        // accumulated concentrations over bead
//         double *conc = new double[ndf];
        for(int i=0; i<ndf; i++) intvals[i] = 0.0;
        
        int node;
        
//         double vol = 0.0; // total volume of bead
        double det;       // local transformation determinant == volume of current tetrahedron
        
        // iterate over all elements that belong to this bead
        for(size_t i=0; i<mien.rows(); i++)
        {
            // localize element values
            for(int j=0; j<nen; j++)
            {
                node = mien(i, j) - 1;
                
                // nodal coordinates
                x[j*2 + 0] = mxyz(node,0);
                x[j*2 + 1] = mxyz(node,1);
                
                // nodal DOFs
                for(int k=0; k<ndf; k++)
                    v[j*ndf + k] = data(node,k);
            }
            
            // compute linear coefficients
            linCoeffs(x, v, ndf, c);
            
            // transform reference quad point locations to physical locations
            reference_to_physical_t3(x, order, xqr, xqp);
            
            // tetra volume by keast library
            det = triangle_area(x);
            
            // add to total bead volume
//             vol += fabs(det);
            
            // compute quad point contributions to integral
            for(int iq=0; iq<order; iq++)
                // for each DOF evaluate at current quad point
                for(int idf=0; idf<ndf; idf++)
                    intvals[idf] += wq[iq] * fabs(det) * (c[idf*3]*xqp[iq*2] + c[idf*3+1]*xqp[iq*2+1] + c[idf*3+2]);
            
        }
        
//         // adapt to multicomponent systems!
//         integral = conc[0];
        
//         // copy to large array
//         for(int i=0; i<ndf; i++)
//             std::cout << std::scientific << std::setw(20) << std::setprecision(8) << conc[i] << std::endl;
        
//         delete[] conc;
        
        delete[] xqr;
        delete[] xqp;
        delete[] wq;
        
        delete[] x;
        delete[] v;
        delete[] c;
    }
    
    
