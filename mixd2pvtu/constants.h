#ifndef CONSTANTS_H_
#define CONSTANTS_H_

#include <mpi.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <cstring>
#include <limits>
#include <time.h>
#include <vtkCellType.h>

using namespace std;
//using namespace MPI;

const double PI = 3.14159265; // Pi constant

const int VTKElemType = VTK_TETRA;       // Element type
const int nen = 4;                       // number of element nodes
const int nef = 4;                       // number of element faces

/* const int VTKElemType = VTK_TRIANGLE;       // Element type */
/* const int nen = 3;                          // number of element nodes */
/* const int nef = 2;                          // number of element faces */

// nsd remains 3 for VTK_TRIANGLE because of convenience. This way extractRNG
// and stitchperiodic can output mxyz with 3 dofs (x,y,z) and we only need to
// recompile mixd2pvtu.
const int nsd = 3;                          // Number of space dimensions.

const int xsd = 0;
const int ysd = 1;
const int zsd = 2;

/* const int edgeNodes[3][2] = {{0,1},{1,2},{2,0}}; */

#endif /* CONSTANTS_H_ */
