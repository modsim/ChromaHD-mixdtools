#include "utils.hpp"
#include <cstring>

#include <iostream>

void swapBytes (char *array, int nelem, int elsize)
{
    int sizet, sizem, i, j;
    char *bytea, *byteb;
    sizet = elsize;
    sizem = sizet - 1;
    bytea = new char [sizet];
    byteb = new char [sizet];
    for (i = 0; i < nelem; i++)
    {
        std::memcpy((void *)bytea, (void *)(array+i*sizet), sizet);
        for (j = 0; j < sizet; j++)
            byteb[j] = bytea[sizem - j];
        std::memcpy((void *)(array+i*sizet), (void *)byteb, sizet);
    }
    delete[] bytea;
    delete[] byteb;

    return;
}

void quickSort(int* arr, int* index, int left, int right)
{
    int i = left, j = right;
    int tmp, tmpIndex;
    int pivot = arr[(left + right) / 2];

    /* partition */
    while (i <= j)
    {
        while (arr[i] < pivot)
            i++;
        while (arr[j] > pivot)
            j--;
        if (i <= j)
        {
            tmp = arr[i];       tmpIndex = index[i];
            arr[i] = arr[j];    index[i] = index[j];
            arr[j] = tmp;       index[j] = tmpIndex;
            i++;
            j--;
        }
    };

    /* recursion */
    if (left < j)
        quickSort(arr, index, left, j);
    if (i < right)
        quickSort(arr, index, i, right);

    return;
}

VTKCellType processElementType(std::string typestring)
{
    if (typestring == "tet")
        return VTK_TETRA;
    else if (typestring == "tri")
        return VTK_TRIANGLE;
    else
    {
        std::cerr << "Invalid element type specified: " << typestring << std::endl;
        std::exit(-1);
    }

}
