
#include "mixd.hpp"

#include <iostream>
#include <iomanip>
#include <cmath>


int main(int argc, char **argv)
{
    using namespace std;
    using namespace mixd;
    
    cout << "MIXD mesh volume tool..." << endl << endl;
    cout << "ONLY FOR 3D TETRA MESHES!!!!!!" << endl << endl;
    
    
    double zmin = -99999999999.9;
    double zmax =  99999999999.9;
    
    if(argc == 3)
    {
        zmin = atof(argv[1]);
        zmax = atof(argv[2]);
    }
    else if(argc != 1)
    {
        cout << "invalid number of arguments" << endl;
        exit(1);
    }
    
    
    long ne, nn;
    readminf("minf", &nn, &ne);
    
    if(nn<=0 || ne<=0) return 1;
    
    
    
    
    MixdFile<double> mxyz("mxyz", nn, 3);
    mxyz.read();
    
    MixdFile<int> mien("mien", ne, 4);
    mien.read();
    
    MixdFile<int> mmat("mmat", ne);
    mmat.init(1);
    
    try{
    mmat.read();
    }catch(MixdException e){}
    
    
    double a[3];
    double b[3];
    double c[3];
    double d[3];
    
    double bxc[3];
    double dot;
    
    double vol = 0.0;
    double vol1 = 0.0;
    double vol2 = 0.0;
    
    
    for(int i=0; i<ne; i++)
    {
        int n1 = mien(i,0) -1;
        int n2 = mien(i,1) -1;
        int n3 = mien(i,2) -1;
        int n4 = mien(i,3) -1;
        
        a[0] = mxyz(n1, 0);
        a[1] = mxyz(n1, 1);
        a[2] = mxyz(n1, 2);
        
        b[0] = mxyz(n2, 0);
        b[1] = mxyz(n2, 1);
        b[2] = mxyz(n2, 2);
        
        c[0] = mxyz(n3, 0);
        c[1] = mxyz(n3, 1);
        c[2] = mxyz(n3, 2);
        
        d[0] = mxyz(n4, 0);
        d[1] = mxyz(n4, 1);
        d[2] = mxyz(n4, 2);
        
        if(  a[2]>zmin && b[2]>zmin && c[2]>zmin && d[2]>zmin
          && a[2]<zmax && b[2]<zmax && c[2]<zmax && d[2]<zmax)
        {
            
        
        for(int j=0; j<3; j++)
        {
            a[j] -= d[j];
            b[j] -= d[j];
            c[j] -= d[j];
        }
        
        
        
        bxc[0] = b[1]*c[2] - b[2]*c[1];
        bxc[1] = b[2]*c[0] - b[0]*c[2];
        bxc[2] = b[0]*c[1] - b[1]*c[0];
        
        
        dot = 0.0;
        for(int j=0; j<3; j++)
            dot += a[j] * bxc[j];
        
        
        vol += abs( dot ) / 6.0;
        
        if(mmat(i) == 1)
            vol1 += abs( dot ) / 6.0;
        
        if(mmat(i) == 2)
            vol2 += abs( dot ) / 6.0;
        }
        
        
    }
    
    cout << endl;
    cout << "mesh volume total = " << scientific << setprecision(20) << vol << endl;
    cout << "mesh volume mat1  = " << scientific << setprecision(20) << vol1 << endl;
    cout << "mesh volume mat2  = " << scientific << setprecision(20) << vol2 << endl;
    cout << endl;
    
    return 0;
}
