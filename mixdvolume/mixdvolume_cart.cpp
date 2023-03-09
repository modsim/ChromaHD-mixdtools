
#include "mixd.hpp"

#include <iostream>
#include <iomanip>
#include <cmath>


double norm(double x, double y, double z)
{
    return sqrt(x*x + y*y + z*z);
}

int main(int argc, char **argv)
{
    using namespace std;
    using namespace mixd;
    
    cout << "MIXD mesh volume tool..." << endl << endl;
    cout << "ONLY FOR 3D CARTESIAN MESHES!!!!!!" << endl << endl;
    
    
    long ne, nn;
    readminf("minf", &nn, &ne);
    
    if(nn<=0 || ne<=0) return 1;
    
    
    
    
    MixdFile<double> mxyz("mxyz", nn, 3);
    mxyz.read();
    
    MixdFile<int> mien("mien", ne, 8);
    mien.read();
    
    MixdFile<int> mmat("mmat", ne);
    mmat.init(1);
    
    try{
    mmat.read();
    }catch(MixdException e){}
    
    
    int n1 = mien(0,0) - 1;
    int n2 = mien(0,1) - 1;
    int n4 = mien(0,3) - 1;
    int n5 = mien(0,4) - 1;
    
    double a = norm( mxyz(n2,0)-mxyz(n1,0), mxyz(n2,1)-mxyz(n1,1), mxyz(n2,2)-mxyz(n1,2) );
    double b = norm( mxyz(n4,0)-mxyz(n1,0), mxyz(n4,1)-mxyz(n1,1), mxyz(n4,2)-mxyz(n1,2) );
    double c = norm( mxyz(n5,0)-mxyz(n1,0), mxyz(n5,1)-mxyz(n1,1), mxyz(n5,2)-mxyz(n1,2) );
    
    double elem_vol = a*b*c;
    
    cout << "edge length a:  " << scientific << setprecision(10) << a << endl;
    cout << "edge length b:  " << scientific << setprecision(10) << b << endl;
    cout << "edge length c:  " << scientific << setprecision(10) << c << endl;
    cout << "element volume: " << scientific << setprecision(10) << elem_vol << endl;
    
    
    
    double vol = 0.0;
    double vol0 = 0.0;
    double vol1 = 0.0;
    double vol2 = 0.0;
    
    
    for(int i=0; i<ne; i++)
    {
        vol += elem_vol;
        
        if(mmat(i) == 0)
            vol0 += elem_vol;
        
        if(mmat(i) == 1)
            vol1 += elem_vol;
        
        if(mmat(i) == 2)
            vol2 += elem_vol;
    }
    
    cout << endl;
    cout << "mesh volume total = " << scientific << setprecision(10) << vol << endl;
    cout << "mesh volume mat0  = " << scientific << setprecision(10) << vol0 << endl;
    cout << "mesh volume mat1  = " << scientific << setprecision(10) << vol1 << endl;
    cout << "mesh volume mat2  = " << scientific << setprecision(10) << vol2 << endl;
    cout << endl;
    
    return 0;
}
