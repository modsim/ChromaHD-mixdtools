#include "endian.h"

#include <memory>
#include <string.h>

//TODO: use <algorithm>'s std::copy() instead of memcpy

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

void endianHandler(int invalue, char * buf)
{
    //how to ensure that buf has size of 4?

    memcpy(buf, &(invalue), 4);

    if(!isBigEndian())
        swapbytes(buf, 1, 4);

}
