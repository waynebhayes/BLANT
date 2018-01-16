#include "misc.h"
#include "bintree.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

int Strcmp(foint a, foint b) { return strcmp((char*)a.v, (char*)b.v); }
foint Strdup(foint a) { return (void*) strdup((char*)a.v);}
void Vfree(foint a) { free(a.v);}

main(int argc, char *argv[])
{
    assert(argc == 2);
{
    FILE *fp = fopen(argv[1], "r");
    BINTREE *tree = BinTreeAlloc(unbalanced, Strcmp, Strdup, Vfree, Strdup, Vfree);
    char key[100], info[100];

    while(fscanf(fp, "%s%s", info, key) == 2)
	BinTreeInsert(tree, (void*)key, (void*)info);
    fclose(fp);

    while(scanf("%s", key) == 1)
    {
	foint name = BinTreeLookup(tree, key);
	printf("%s\n", name.i == ABSTRACT_ERROR.i ? "bozo" : (char*)name.v);
    }

    BinTreeFree(tree);
    return 0;
}
}
