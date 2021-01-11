#include "misc.h"
#include "htree.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#define DEPTH 3
int main(int argc, char *argv[])
{
    assert(argc == 2);
{
    FILE *fp = fopen(argv[1], "r");
    HTREE *h = HTreeAlloc(DEPTH);
    char buf[DEPTH][BUFSIZ];
    char *keys[DEPTH];
    foint data;
    int i;
    for(i=0;i<DEPTH;i++) keys[i] = buf[i];

    while(fscanf(fp, "%s%s%s%d", keys[0], keys[1], keys[2], &data.i) == 4) {
	printf("Inserting %s %s %s = %d\n", keys[0], keys[1], keys[2], data.i);
	HTreeInsert(h, keys, (foint)data);
	foint found;
	assert(HTreeLookup(h, keys, &found));
	assert(found.i == data.i);
    }
    fclose(fp);

    while(scanf("%s %s %s", keys[0], keys[1], keys[2]) == 3)
    {
	int sizes[DEPTH]={-1,-1,-1};
	int depth = HTreeSizes(h, keys, sizes);
	printf("filled %d sizes [%d %d %d]\n",depth, sizes[0],sizes[1],sizes[2]);
	if(HTreeLookup(h, keys, &data))
	    printf("%d\n", data.i);
	else
	    puts("not found");
    }

    HTreeFree(h);
    return 0;
}
}
