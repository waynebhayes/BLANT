#ifndef _MULTI_SETS_H
#define _MULTI_SETS_H
/*
    Array (not hash table) representation of multiset.
    Uses unsigned char to hold frequency of each element.
    Assumes all elements are between 0 and sizeof(set)-1.
    Designed to quickly return and maintain cardinality of multiset.
*/

typedef unsigned char FREQTYPE; //Using unsigned char limits multiset to 255 frequency
#define MAX_MULTISET_NUM 255 //The largest number that can be stored in freqtype
typedef struct _multisetType {
    unsigned int n; //is sizeof the array
    unsigned int cardinality; //number of distinct elements
    FREQTYPE* array;
} MULTISET;

MULTISET *MultisetAlloc(unsigned int n);
MULTISET *MultisetResize(MULTISET *mset, unsigned int new_n);
void MultisetFree(MULTISET *mset);
#define MultisetReset MultisetEmpty
MULTISET *MultisetEmpty(MULTISET *mset);
unsigned MultisetCardinality(MULTISET *mset);

MULTISET *MultisetAdd(MULTISET *mset, unsigned element);    /* add single element to multiset */
MULTISET *MultisetDelete(MULTISET *mset, unsigned element); /* delete a single element */


#endif