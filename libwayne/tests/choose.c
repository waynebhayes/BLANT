#include "misc.h"
#include "combin.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
    int n = atoi(argv[1]), m=atoi(argv[2]);
    printf("%llu\n", CombinChoose(n,m));
    return 0;
}
