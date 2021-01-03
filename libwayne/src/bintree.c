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
#include <math.h> // for logarithm

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
    tree->root = NULL;
    tree->type = type;
    tree->cmpKey = cmpKey ? cmpKey : CmpInt;
    tree->copyKey = copyKey ? copyKey : CopyInt;
    tree->freeKey = freeKey ? freeKey : FreeInt;
    tree->copyInfo = copyInfo ? copyInfo : CopyInt;
    tree->freeInfo = freeInfo ? freeInfo : FreeInt;
    tree->n = tree->depth = 0;
    return tree;
}


static void BinTreeRebalance(BINTREE *tree);
void BinTreeInsert(BINTREE *tree, foint key, foint info)
{
    int depth = 0;
    BINTREENODE *p = tree->root, **locative = &(tree->root);
    while(p)
    {
	if(++depth > tree->depth) tree->depth = depth;
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
    tree->n++;

    if(tree->n > 4 && tree->depth > 4*log(tree->n)) BinTreeRebalance(tree);
}


Boolean BinTreeLookup(BINTREE *tree, foint key, foint *pInfo)
{
    int depth = 0;
    BINTREENODE *p = tree->root;
    while(p)
    {
	if(++depth > tree->depth) tree->depth = depth;
	int cmp = tree->cmpKey(key, p->key);
	if(cmp == 0)
	{
	    *pInfo = p->info;
	    return true;
	}
	if(cmp < 0)
	    p = p->left;
	else
	    p = p->right;
    }
    return false;
}

static Boolean BinTreeTraverseHelper ( BINTREENODE *p, pFointTraverseFcn f)
{
    Boolean cont = true;
    if(p) {
	if(p->left) cont = cont && BinTreeTraverseHelper(p->left, f);
	if(cont) cont = cont && f(p->key, p->info);
	if(cont && p->right) cont = cont && BinTreeTraverseHelper(p->right, f);
    }
    return cont;
}

Boolean BinTreeTraverse ( BINTREE *tree, pFointTraverseFcn f)
{
    return BinTreeTraverseHelper(tree->root, f);
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


//////////////////// REBALANCING CODE
static foint *keyArray, *dataArray;
static int arraySize, currentItem;
static Boolean inRebalance;

// Squirrel away all the item *in sorted order*
static Boolean TraverseTreeToArray(foint key, foint data) {
    assert(currentItem < arraySize);
    keyArray[currentItem] = key;
    dataArray[currentItem] = data;
    ++currentItem;
    return true;
}

static void BinTreeInsertMiddleElementOfArray(BINTREE *tree, int low, int high) // low to high inclusive
{
    if(low <= high) { // we're using low though high *inclusive* so = is a valid case.
	int mid = (low+high)/2;
	BinTreeInsert(tree, keyArray[mid], dataArray[mid]);
	BinTreeInsertMiddleElementOfArray(tree, low, mid-1);
	BinTreeInsertMiddleElementOfArray(tree, mid+1,high);
    }
}

static void BinTreeRebalance(BINTREE *tree)
{
    assert(!inRebalance);
    inRebalance = true;
    //fprintf(stderr,"BinTreeRebalance: old n %d depth %d\n", tree->n, tree->depth);
    if(tree->n > arraySize){
	arraySize = tree->n;
	keyArray = Realloc(keyArray, arraySize*sizeof(foint));
	dataArray = Realloc(dataArray, arraySize*sizeof(foint));
    }
    currentItem = 0;
    BinTreeTraverse (tree, TraverseTreeToArray);
    assert(currentItem == tree->n);
    
    BINTREE *newTree = BinTreeAlloc(tree->type, tree->cmpKey , tree->copyKey , tree->freeKey , tree->copyInfo , tree->freeInfo);
    // Now re-insert the items in *perfectly balanced* order.
    BinTreeInsertMiddleElementOfArray(newTree, 0, tree->n - 1);
    assert(tree->n == newTree->n);
    assert(newTree->depth < tree->depth);
    assert(newTree->depth <= 1+log(newTree->n)/log(2));
    // Swap the roots, record new depth.
    BINTREENODE *tmp = tree->root; tree->root = newTree->root; newTree->root = tmp;
    int itmp = tree->depth; tree->depth = newTree->depth; newTree->depth = itmp;
    BinTreeFree(newTree);
    //fprintf(stderr,"BinTreeRebalance: new n %d depth %d\n", tree->n, tree->depth);
    inRebalance = false;
}
