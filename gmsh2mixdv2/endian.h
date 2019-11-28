#ifndef ENDIAN_H
#define ENDIAN_H

bool isBigEndian();
void swapbytes(char *array, int nelem, int elsize);
void endianHandler(int invalue, char * buf);


#endif
