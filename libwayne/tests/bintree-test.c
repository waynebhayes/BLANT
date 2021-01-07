#include "misc.h"
#include "bintree.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>


int main(int argc, char *argv[])
{
    assert(argc == 2);
{
    FILE *fp = fopen(argv[1], "r");
    BINTREE *tree = BinTreeAlloc(strcmp, strdup, free, NULL, NULL);
    char buf[100];
    foint key, data;

    while(fscanf(fp, "%s %d", buf, &data.i) == 2)
	BinTreeInsert(tree, (foint)(key.s=buf), data);
    fclose(fp);

    while(scanf("%s", buf) == 1)
    {
	if(BinTreeLookup(tree, (foint)(key.s=buf), &data))
	    printf("%d\n", data.i);
	else
	    puts("nope");
    }

    BinTreeFree(tree);
    return 0;
}
}
