#ifndef _BINTREE_H
#define _BINTREE_H

#include <malloc.h>
#include <stdio.h>
#include "misc.h"   /* for foint */


/*-------------------  Types  ------------------*/

typedef struct _binTreeNode
{
    foint key, info;
    struct _binTreeNode *left, *right, *parent;
} BINTREENODE;

enum _treeType { unbalanced = 0 };   /* AVL next? */

typedef struct _binTree
{
    enum _treeType type;
    BINTREENODE *root;
    pCmpFcn cmpKey;
    pFointCopyFcn copyKey, copyInfo;
    pFointFreeFcn freeKey, freeInfo;
} BINTREE;

/*-----------   Function Prototypes  -----------*/

BINTREE *BinTreeAlloc(enum _treeType type, pCmpFcn cmpKey,
    pFointCopyFcn copyKey, pFointFreeFcn freeKey,
    pFointCopyFcn copyInfo, pFointFreeFcn freeInfo);

void BinTreeInsert(BINTREE *, foint key, foint info);

/* O(log n); returns false if failure, and true if found and then sets returnInfo to the pointer to foint; */
Boolean BinTreeLookup(BINTREE *, foint key, foint *pInfo);

/*
** BinTreeTraverse: Traverse a binary tree, calling your function pointer (pFointTraversalFcn) on each
** element, in order.
*/
Boolean BinTreeTraverse ( BINTREE *bt, pFointTraverseFcn);

void BinTreeFree(BINTREE *);

#endif  /* _BINTREE_H */
