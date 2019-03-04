#include <stdio.h>
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include "mt19937.h"

int main(){
    int i;
    MT19937 *a = Mt19937Alloc(time(0)+getpid());
    for(i=0;i<1000000;i++)
	printf("%.16lf\n", Mt19937NextDouble(a));
    Mt19937Free(a);
return 0;
}
