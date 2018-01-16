/* Version 0.0
** From "Wayne's Little DSA Library" (DSA == Data Structures and
** Algorithms) Feel free to change, modify, or burn these sources, but if
** you modify them please don't distribute your changes without a clear
** indication that you've done so.  If you think your change is spiffy,
** send it to me and maybe I'll include it in the next release.
** 
** Wayne Hayes, wayne@cs.utoronto.ca (preffered), or wayne@cs.toronto.edu
*/ 

/* Binary tree algorithms from Lewis & Denenberg
*/
#include <string.h>

#include "bintree.h"

static foint CopyInt(foint i)
{
    return i;
}

static void FreeInt(foint i) {}

static int CmpInt(foint i, foint j) { return i.i - j.i; }


BINTREE *BinTreeAlloc(enum _treeType type, pCmpFcn cmpKey,
    pFointCopyFcn copyKey, pFointFreeFcn freeKey,
    pFointCopyFcn copyInfo, pFointFreeFcn freeInfo)
{
    BINTREE *tree = Malloc(sizeof(BINTREE));
    tree->type = type;
    tree->cmpKey = cmpKey ? cmpKey : CmpInt;
    tree->root = NULL;
    tree->copyKey = copyKey ? copyKey : CopyInt;
    tree->freeKey = freeKey ? freeKey : FreeInt;
    tree->copyInfo = copyInfo ? copyInfo : CopyInt;
    tree->freeInfo = freeInfo ? freeInfo : FreeInt;
    return tree;
}

void BinTreeInsert(BINTREE *tree, foint key, foint info)
{
    BINTREENODE *p = tree->root, **locative = &(tree->root);

    while(p)
    {
	int cmp = tree->cmpKey(key, p->key);
	if(cmp == 0)
	{
	    tree->freeInfo(p->info);
	    p->info = tree->copyInfo(info);
	    return;
	}
	else if(cmp < 0)
	{
	    locative = &(p->left);
	    p = p->left;
	}
	else
	{
	    locative = &(p->right);
	    p = p->right;
	}
    }

    p = (BINTREENODE*) Calloc(1,sizeof(BINTREENODE));
    p->key = tree->copyKey(key);
    p->info = tree->copyInfo(info);
    p->left = p->right = NULL;
    *locative = p;
}


foint BinTreeLookup(BINTREE *tree, foint key)
{
    BINTREENODE *p = tree->root;
    while(p)
    {
	int cmp = tree->cmpKey(key, p->key);
	if(cmp == 0)
	    return p->info;
	if(cmp < 0)
	    p = p->left;
	else
	    p = p->right;
    }
    return ABSTRACT_ERROR;
}


#if 0
/* LookupKey: the tree is ordered on key, not info.  So we have to do a
** full traversal, O(size of tree), not O(log(size of tree)).  This is
** only used when errors occur.
*/
foint BinTreeLookupKey(BINTREE *tree, foint info)
{
    BINTREENODE *p = tree->root;
    if(p)
    {
	if(tree->cmpKey(p->info, info) == 0)
	    return p->key;
	else
	{
	    foint key = BinTreeLookupKey(p->left, info);
	    if(key)
		return key;
	    else
		return BinTreeLookupKey(p->right, info);
	}
    }
    return NULL;
}
#endif


static void BinTreeNodeFree(BINTREE *tree, BINTREENODE *t)
{
    if(t)
    {
	BinTreeNodeFree(tree, t->left);
	BinTreeNodeFree(tree, t->right);
	tree->freeKey(t->key);
	tree->freeInfo(t->info);
	free(t);
    }
}

void BinTreeFree(BINTREE *tree)
{
    BinTreeNodeFree(tree, tree->root);
    free(tree);
}
