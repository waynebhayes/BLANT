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
    BINTREE *tree = BinTreeAlloc((pCmpFcn)strcmp, (pFointCopyFcn)strdup, (pFointFreeFcn)free, NULL, NULL);
    char buf[BUFSIZ], bufs[BUFSIZ][BUFSIZ];
    foint key, data;
    int lines=0;

    while(fscanf(fp, "%s %d", &bufs[lines], &data.i) == 2)
	BinTreeInsert(tree, (foint)(key.s=bufs[lines++]), data);
    fclose(fp);

    while(scanf("%s", buf) == 1)
    {
	if(BinTreeLookup(tree, (foint)(key.s=buf), &data))
	    printf("%d\n", data.i);
	else
	    puts("nope");
    }
    assert(false == BinTreeLookup(tree, (foint)(key.s="foo"), &data));
    assert(false == BinTreeLookup(tree, (foint)(key.s="bar"), &data));
    assert(false == BinTreeLookup(tree, (foint)(key.s="foobar"), &data));

    while(lines>0) assert(BinTreeDelete(tree, (foint)(key.s=bufs[--lines])));
    assert(tree->n == 0);
    BinTreeFree(tree);
    return 0;
}
}
