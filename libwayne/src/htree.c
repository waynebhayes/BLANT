#include "htree.h" // bintree is included there

/*-------------------  Types  ------------------*/

HTREE *HTreeAlloc(int depth) {
    assert(depth>0);
    HTREE *h = Calloc(1, sizeof(HTREE));
    h->depth = depth;
    // Be nice and create the top-level tree.
    h->tree = BinTreeAlloc((pCmpFcn) strcmp, (pFointCopyFcn) strdup, (pFointFreeFcn) free, NULL,NULL);
    return h;
}


static void HTreeInsertHelper(HTREE *h, int currentDepth, BINTREE *tree, foint keys[], foint data)
{
    assert(tree && 0 <= currentDepth && currentDepth < h->depth);
    if(currentDepth == h->depth-1) // we're hit the lowest level tree; its data elements are the final elements.
	BinTreeInsert(tree, keys[currentDepth], data);
    else {
	// Otherwise, we are NOT at the lowest level tree; the data members of these nodes are themselves other trees,
	// so to find the next tree we use the key at this level to *look up* the binary tree at the next level down
	foint nextLevel;
	BINTREE *nextTree;
	if(!BinTreeLookup(tree, keys[currentDepth], &nextLevel)) {
	    nextTree = BinTreeAlloc((pCmpFcn) strcmp, (pFointCopyFcn) strdup, (pFointFreeFcn) free, NULL,NULL);
	    nextLevel.v = nextTree;
	    BinTreeInsert(tree, keys[currentDepth], nextLevel);
	}
	else
	    nextTree = nextLevel.v;
	assert(nextTree);
	HTreeInsertHelper(h, currentDepth+1, nextTree, keys, data);
    }
}

// key is an array with exactly "depth" elements, data is what you want to put at the lowest level.
void HTreeInsert(HTREE *h, foint keys[], foint data)
{
    HTreeInsertHelper(h, 0, h->tree, keys, data);
}

static Boolean HTreeLookupHelper(HTREE *h, int currentDepth, BINTREE *tree, foint keys[], foint *pData)
{
    assert(tree && 0 <= currentDepth && currentDepth < h->depth);
    if(currentDepth == h->depth-1) // we're hit the lowest level tree; its data elements are the final elements.
	return BinTreeLookup(tree, keys[currentDepth], pData);
    else {
	foint nextLevel;
	BINTREE *nextTree;
	if(!BinTreeLookup(tree, keys[currentDepth], &nextLevel))
	    return false;
	else
	    nextTree = nextLevel.v;
	assert(nextTree);
	return HTreeLookupHelper(h, currentDepth+1, nextTree, keys, pData);
    }
}

Boolean HTreeLookup(HTREE *h, foint keys[], foint *pData)
{
    return HTreeLookupHelper(h, 0, h->tree, keys, pData);
}

static int HTreeSizesHelper(HTREE *h, int currentDepth, BINTREE *tree, foint keys[], int sizes[])
{
    assert(tree && 0 <= currentDepth && currentDepth < h->depth);
    sizes[currentDepth] = tree->n;
    if(currentDepth == h->depth-1) // we're hit the lowest level tree; its data elements are the final elements.
	return 1;
    else {
	foint nextLevel;
	BINTREE *nextTree;
	if(!BinTreeLookup(tree, keys[currentDepth], &nextLevel))
	    return 1;
	else
	    nextTree = nextLevel.v;
	assert(nextTree);
	return 1 + HTreeSizesHelper(h, currentDepth+1, nextTree, keys, sizes);
    }
}

int HTreeSizes(HTREE *h, foint keys[], int sizes[])
{
    return HTreeSizesHelper(h, 0, h->tree, keys, sizes);
}


static void HTreeFreeHelper(HTREE *h, int currentDepth, BINTREE *tree);
static HTREE *_TraverseH;
static int _TraverseDepth;
static void TraverseFree(foint key, foint data) {
    assert(_TraverseDepth < _TraverseH->depth);
    BINTREE *t = data.v;
    int depth = _TraverseDepth;
    HTreeFreeHelper(_TraverseH, _TraverseDepth+1, t);
    _TraverseDepth = depth;
}

static void HTreeFreeHelper(HTREE *h, int currentDepth, BINTREE *tree)
{
    assert(tree && 0 <= currentDepth && currentDepth < h->depth);
    if(currentDepth == h->depth-1) // we're hit the lowest level tree; its data elements are the final elements.
	BinTreeFree(tree);
    else {
	_TraverseH = h; _TraverseDepth = currentDepth;
	BinTreeTraverse(tree, (pFointTraverseFcn) TraverseFree);
    }
}

void HTreeFree(HTREE *h)
{
    HTreeFreeHelper(h, 0, h->tree);
}
