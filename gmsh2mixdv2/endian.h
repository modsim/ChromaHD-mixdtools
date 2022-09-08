#ifndef ENDIAN_H
#define ENDIAN_H

#include <stdio.h>

bool isBigEndian();
void swapbytes(char *array, size_t nelem, int elsize);
void endianHandler(int invalue, char * buf);


#endif
