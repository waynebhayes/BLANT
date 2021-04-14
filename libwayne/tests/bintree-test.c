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

    while(fscanf(fp, "%s %d", bufs[lines], &data.i) == 2) {
	BinTreeInsert(tree, (foint)(key.s=bufs[lines++]), data);
	BinTreeSanityCheck(tree);
    }
    fclose(fp);
    BinTreeSanityCheck(tree);
    BinTreeRebalance(tree);
    BinTreeSanityCheck(tree);

    while(scanf("%s", buf) == 1)
    {
	BinTreeSanityCheck(tree);
	if(BinTreeLookup(tree, (foint)(key.s=buf), &data))
	    printf("%d\n", data.i);
	else
	    puts("nope");
	++lines;
	BinTreeSanityCheck(tree);
    }
    assert(false == BinTreeLookup(tree, (foint)(key.s="foo"), &data)); BinTreeSanityCheck(tree);
    assert(false == BinTreeLookup(tree, (foint)(key.s="bar"), &data)); BinTreeSanityCheck(tree);
    assert(false == BinTreeLookup(tree, (foint)(key.s="foobar"), &data)); BinTreeSanityCheck(tree);

    while(lines>0) {
	BinTreeSanityCheck(tree);
	Boolean exists = BinTreeLookup(tree, (foint)(key.s=bufs[--lines]), &data);
	assert(exists == BinTreeDelete(tree, (foint)(key.s=bufs[lines])));
	printf("Deleted <%s>; n=%d, phys=%d\n", bufs[lines], tree->n, tree->physical_n);
	BinTreeSanityCheck(tree);
    }
    assert(lines == 0 && tree->n == 0 && tree->physical_n >= 0);
    BinTreeSanityCheck(tree);
    BinTreeFree(tree);
    return 0;
}
}
