#include "rand48.h"
#include "stream48.h"

int main(void)
{
    int nSt = 3, i;

    Stream48Init(nSt);

    for(i=0; i<nSt; i++)
    {
	Stream48(0);
	/*Stream48Randomize();*/
    }

    while(1)
    {
	for(i=0; i<nSt; i++)
	{
	    Stream48(i);
	    printf("%.16g, ", drand48());
	}
	printf("\n");
    }
    return 0;
}
