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
    char key[100];
    int data;

    while(fscanf(fp, "%s %d", key,&data) == 2)
	BinTreeInsert(tree, (foint)key, (foint)data);
    fclose(fp);

    while(scanf("%s", key) == 1)
    {
	if(BinTreeLookup(tree, (foint) key, (foint*)&data))
	    printf("%d\n", data);
	else
	    puts("nope");
    }

    BinTreeFree(tree);
    return 0;
}
}
