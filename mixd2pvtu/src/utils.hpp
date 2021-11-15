#ifndef UTILS_H_
#define UTILS_H_

#include <string>
#include <vtk/vtkCellType.h>

void swapBytes (char *array, int nelem, int elsize);
void quickSort(int* arr, int* index, int left, int right);
VTKCellType processElementType(std::string typestring);

#endif // !UTILS_H_
