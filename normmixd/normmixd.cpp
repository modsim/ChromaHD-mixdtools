#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>





bool isBigEndian()
{
    short word = 0x4321;
    if ((* (char*) & word) != 0x21) return true;
    else      return false;
}

void swapbytes(char *array, int nelem, int elsize)
{
    register int sizet, sizem, i, j;
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



/* a tool to dump contents of binary MIXD files to the console */

int main(int argc, char **argv)
{

    if(argc != 3)
    {
        std::cout << "Dump MIXD" << std::endl;
        std::cout << "Usage: " << argv[0] << " <number of items per row> <MIXD file>" << std::endl;
        return 1;
    }

    /* const char type = argv[1][0]; */

	int itsize = 8; //double

	/* switch(type) */
	/* { */
	/* 	case 'i':   itsize = 4; break; */
	/* 	case 'f':   itsize = 4; break; */
	/* 	case 'd':   itsize = 8; break; */

		/* default: */
		/* 	std::cout << "Error: illegal data type: " << argv[1] << std::endl; */
		/* 	return 1; */
	/* } */

	const int nitems = atoi(argv[1]);

	const char *fname = argv[2];

    std::cout << "Reading file " << fname << std::endl;



    FILE *fid = fopen(fname, "rb");

    if (fid == NULL)
    {
        printf("\nCould not open %s\n", fname);
        return 1;
    }



    char *buf = (char*)malloc(nitems * itsize);


	bool isBigEnd = isBigEndian();

    double l2norm[nitems], infnorm[nitems];

    for(int i=0; i< nitems; i++)
    {
        l2norm[i] = 0;
        infnorm[i] = 0;

    }


	std::cout << std::endl;

	int cnt = 0;


    while(true)
	{

		int itemsread = fread((void*)buf, itsize, nitems, fid);

		if(itemsread != nitems) break;

		if(!isBigEnd)
			swapbytes(buf, nitems, itsize);


		cnt++;
		/* printf("%4d", cnt); */


		for(char *it=buf; it<buf+nitems*itsize; it+=itsize)
		{
		    /* double val = *((double*)it);   printf("%30.16f", val); */
            l2norm[(it-buf)/itsize] += *((double*)it)**((double*)it);
            if(*((double*)it) > infnorm[(it-buf)/itsize])
                infnorm[(it-buf)/itsize] = *((double*)it);
		}

		/* std::cout << std::endl; */



	}

    printf("\tL2\t\tL-inf\n");
    for(int i=0; i< nitems; i++)
	{
        l2norm[i] = sqrt(l2norm[i]);
		printf("[%d]: %e\t%e\n", i, l2norm[i], infnorm[i] );
		/* std::cout << l2norm[i] << std::endl; */
	}

    free(buf);

	fclose(fid);

    return 0;
}


