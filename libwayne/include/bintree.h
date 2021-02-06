#ifndef _BINTREE_H
#define _BINTREE_H

#include <malloc.h>
#include <stdio.h>
#include "misc.h"   /* for foint */


/*-------------------  Types  ------------------*/

typedef struct _binTreeNode
{
    foint key, info;
    struct _binTreeNode *left, *right;
} BINTREENODE;

typedef struct _binTree
{
    int n, depthSum, depthSamples; // number of entries, and tree depth statistics which can be averaged anytime.
    BINTREENODE *root;
    pCmpFcn cmpKey;
    pFointCopyFcn copyKey, copyInfo;
    pFointFreeFcn freeKey, freeInfo;
} BINTREE;

/*-----------   Function Prototypes  -----------*/

BINTREE *BinTreeAlloc(pCmpFcn cmpKey, pFointCopyFcn copyKey, pFointFreeFcn freeKey,
	    pFointCopyFcn copyInfo, pFointFreeFcn freeInfo);

void BinTreeInsert(BINTREE *, foint key, foint info); // replaces info if the key already exists

/* O(log n); returns false if failure, and true if found and then assigns to *pInfo if *pInfo is non-NULL */
Boolean BinTreeLookup(BINTREE *, foint key, foint *pInfo);

// O(log n), returns true if deleted successfully, false if item not found.
Boolean BinTreeDelete(BINTREE *, foint key);

/*
** BinTreeTraverse: Traverse a binary tree, calling your function pointer (pFointTraversalFcn) on each
** element, in order.
*/
Boolean BinTreeTraverse ( BINTREE *bt, pFointTraverseFcn);

void BinTreeFree(BINTREE *);

#endif  /* _BINTREE_H */
