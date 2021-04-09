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


BINTREE *BinTreeAlloc(pCmpFcn cmpKey,
    pFointCopyFcn copyKey, pFointFreeFcn freeKey,
    pFointCopyFcn copyInfo, pFointFreeFcn freeInfo)
{
    BINTREE *tree = Malloc(sizeof(BINTREE));
    tree->root = NULL;
    tree->cmpKey = cmpKey ? cmpKey : CmpInt;
    tree->copyKey = copyKey ? copyKey : CopyInt;
    tree->freeKey = freeKey ? freeKey : FreeInt;
    tree->copyInfo = copyInfo ? copyInfo : CopyInt;
    tree->freeInfo = freeInfo ? freeInfo : FreeInt;
    tree->physical_n = tree->n = tree->depthSum = tree->depthSamples = 0;
    return tree;
}


void BinTreeRebalance(BINTREE *tree);
void BinTreeInsert(BINTREE *tree, foint key, foint info)
{
    int depth = 0;
    BINTREENODE *p = tree->root, **locative = &(tree->root);
    while(p)
    {
	++depth;
	int cmp = tree->cmpKey(key, p->key);
	if(cmp == 0)
	{
	    tree->freeInfo(p->info);
	    p->info = tree->copyInfo(info);
	    if(p->deleted) { p->deleted = false; tree->n++; }
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
    tree->n++; tree->physical_n++;

    tree->depthSum += depth; ++tree->depthSamples;
    double meanDepth = tree->depthSum/(double)tree->depthSamples;
    if(tree->n > 5 && tree->depthSamples > 20 && meanDepth > 3*log(tree->n)) BinTreeRebalance(tree);
}


Boolean BinTreeDelete(BINTREE *tree, foint key)
{
    int depth = 0;
    BINTREENODE *p = tree->root, **locative = &(tree->root);
    while(p)
    {
	++depth;
	int cmp = tree->cmpKey(key, p->key);
	if(cmp == 0)
	    break;
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
    if(!p) return false;

    // At this point, p points to the node we want to delete. If either child is NULL, then the other child moves up.

    if(p->left && p->right) // can't properly delete, so just mark as deleted and it'll go away when rebalance happens
	p->deleted = true;
    else if(p->left)  *locative = p->left;
    else if(p->right) *locative = p->right;
    else *locative = NULL;

    if(!p->deleted) { // if it was marked as deleted, everything needs to remain; otherwise we can nuke everything.
	tree->physical_n--;
	tree->freeKey(p->key);
	tree->freeInfo(p->info);
	Free(p);
    }

    tree->n--;
    return true;
}


Boolean BinTreeLookup(BINTREE *tree, foint key, foint *pInfo)
{
    int depth=0;
    BINTREENODE *p = tree->root;
    while(p)
    {
	++depth;
	int cmp = tree->cmpKey(key, p->key);
	if(cmp == 0)
	{
	    if(p->deleted) return false;
	    if(pInfo) *pInfo = p->info;
	    return true;
	}
	if(cmp < 0)
	    p = p->left;
	else
	    p = p->right;
    }
    tree->depthSum += depth; ++tree->depthSamples;
    double meanDepth = tree->depthSum/(double)tree->depthSamples;
    if(tree->n > 5 && tree->depthSamples > 20 && meanDepth > 3*log(tree->n)) BinTreeRebalance(tree);
    return false;
}

static Boolean BinTreeTraverseHelper ( BINTREENODE *p, pFointTraverseFcn f)
{
    Boolean cont = true;
    if(p) {
	if(p->left) cont = BinTreeTraverseHelper(p->left, f);
	if(cont && !p->deleted) cont = f(p->key, p->info);
	if(cont && p->right) cont = BinTreeTraverseHelper(p->right, f);
    }
    return cont;
}

Boolean BinTreeTraverse ( BINTREE *tree, pFointTraverseFcn f)
{
    return BinTreeTraverseHelper(tree->root, f);
}

static int _binTreeSanityNodeCount, _binTreeSanityPhysicalNodeCount;
static Boolean BinTreeSanityHelper ( BINTREENODE *p, pCmpFcn cmpKey )
{
    if(p) {
	++_binTreeSanityPhysicalNodeCount;
	if(p->left) { assert(cmpKey(p->left->key, p->key)<0); BinTreeSanityHelper(p->left, cmpKey); }
	if(!p->deleted) ++_binTreeSanityNodeCount;
	if(p->right) { assert(cmpKey(p->key, p->right->key)<0); BinTreeSanityHelper(p->right, cmpKey); }
    }
    return true;
}

Boolean BinTreeSanityCheck ( BINTREE *tree)
{
    _binTreeSanityNodeCount = 0;
    _binTreeSanityPhysicalNodeCount = 0;
    assert(0 <= tree->n && tree->n <= tree->physical_n && tree->physical_n >= 0);
    BinTreeSanityHelper(tree->root, tree->cmpKey);
    assert(_binTreeSanityNodeCount == tree->n);
    assert(_binTreeSanityPhysicalNodeCount == tree->physical_n);
    return true;
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


static void BinTreeFreeHelper(BINTREE *tree, BINTREENODE *t)
{
    if(t)
    {
	BinTreeFreeHelper(tree, t->left);
	BinTreeFreeHelper(tree, t->right);
	tree->freeKey(t->key);
	tree->freeInfo(t->info);
	if(!t->deleted) tree->n--;
	free(t);
	tree->physical_n--;
    }
}

void BinTreeFree(BINTREE *tree)
{
    BinTreeFreeHelper(tree, tree->root);
    assert(tree->n == 0 && tree->physical_n == 0);
    free(tree);
}


//////////////////// REBALANCING CODE
static foint *keyArray, *dataArray;
static int arraySize, currentItem;

// Squirrel away all the items *in sorted order*
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

void BinTreeRebalance(BINTREE *tree)
{
    static Boolean inRebalance; // only allow one tree to be rebalanced concurrently
    if(inRebalance) return; // even without threading, the BinTreeTraverse below may trigger the rebalance of another tree
    //Warning("inRebalance tree %x size %d mean depth %g", tree, tree->n, tree->depthSum/(double)tree->depthSamples);
    inRebalance = true;
    if(tree->n > arraySize){
	arraySize = tree->n;
	keyArray = Realloc(keyArray, arraySize*sizeof(foint));
	dataArray = Realloc(dataArray, arraySize*sizeof(foint));
    }
    currentItem = 0;
    BinTreeTraverse (tree, TraverseTreeToArray);
    assert(currentItem == tree->n);
    
    BINTREE *newTree = BinTreeAlloc(tree->cmpKey , tree->copyKey , tree->freeKey , tree->copyInfo , tree->freeInfo);
    // Now re-insert the items in *perfectly balanced* order.
    BinTreeInsertMiddleElementOfArray(newTree, 0, tree->n - 1);
    assert(tree->n == newTree->n);
    assert(newTree->n == newTree->physical_n);
    // Swap the roots, record new depth.
    BINTREENODE *tmp = tree->root; tree->root = newTree->root; newTree->root = tmp;
    BinTreeFree(newTree);
    inRebalance = false;
}
