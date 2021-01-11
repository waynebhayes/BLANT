#ifndef _HTREE_H
#define _HTREE_H

// A hierarchical binary tree is a binary tree in which the *data* member of each node is itself another binary tree.
// It is best used when each key K at one level has associated with it many sub-keys of a different sort; these sub-keys
// are associated only with K. And so it goes down the hierarchy.
// The hierarchy has a fixed depth, and when you insert or lookup an element, you specify a depth you wish to go,
// and an array of keys with length equal to the depth you're searching.
// For now, to keep things simple, we assume the key is *always* a string, and the data is always a simple foint.
// That way all key operations use strcmp, strdup, and free.

#include "bintree.h"

/*-------------------  Types  ------------------*/

typedef struct _hTree
{
    BINTREE *tree;
    unsigned char depth; // we're going to assume it's less than 255 layers deep, OK?
    int n; // total number of elements across all sub-trees.
} HTREE;

/*-----------   Function Prototypes  -----------*/

HTREE *HTreeAlloc(int depth); // Allocate an HTREE of specified depth.

// key is an array with exactly "depth" elements, info is what you want to put at the lowest level.
void HTreeInsert(HTREE *, foint keys[], foint info);

Boolean HTreeLookup(HTREE *, foint keys[], foint *pInfo); // as above, and the lowest level info is returned if it exists.

// number of elements in trees down the hierarchy along key path; returns number of sizes[] we managed to fill.
int HTreeSizes(HTREE *, foint keys[], int sizes[]);

void HTreeFree(HTREE *);

#endif  /* _HTREE_H */
