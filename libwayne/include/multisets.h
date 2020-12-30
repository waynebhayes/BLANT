#ifndef _MULTI_SETS_H
#define _MULTI_SETS_H
#include "sets.h"
/*
    Array (not hash table) representation of multiset.
    Uses unsigned char to hold multiplicity(frequency) of each element.
    Assumes all elements are between 0 and sizeof(set)-1.
    Designed to quickly return and maintain the support(number of distinct elements) of the multiset.
*/

typedef unsigned short int FREQTYPE; //Using unsigned char limits multiset to a multiplicty of 255
#define MAX_MULTISET_FREQ 65535 //The largest number that can be stored in freqtype
typedef struct _multisetType {
    unsigned int n; //is sizeof the array
    unsigned int support; //number of distinct elements
    FREQTYPE* array;
    SET *set; // The set version, without keeping count.
} MULTISET;

MULTISET *MultisetAlloc(unsigned int n);
MULTISET *MultisetResize(MULTISET *mset, unsigned int new_n);  /* resizes and copies elements that still fits */
void MultisetFree(MULTISET *mset);
#define MultisetReset MultisetEmpty
MULTISET *MultisetEmpty(MULTISET *mset); /* empties the multiset */
unsigned MultisetSupport(MULTISET *mset); /* returns the support of the multiset */
FREQTYPE MultisetMultiplicity(MULTISET *mset, unsigned element); /* returns the multiplicity of an element in the multiset */

MULTISET *MultisetAdd(MULTISET *mset, unsigned element);    /* add single element to multiset */
MULTISET *MultisetComplement(MULTISET *mset, FREQTYPE N);   /* array[e] = N-array[e] for every element e */
MULTISET *MultisetAddSet(MULTISET *mset, SET *s);    /* add a set of elements to a multiset */
SET *MultisetToSet(SET *s, MULTISET *mset);    /* SET will have same list of non-zero elements as MULTISET */
MULTISET *MultisetDelete(MULTISET *mset, unsigned element); /* delete a single element */
MULTISET *MultisetDeleteSet(MULTISET *mset, SET *s); /* remove a set of elements from a multiset */
MULTISET *MultisetSum(MULTISET *C, MULTISET *A, MULTISET *B);  /* C = sum of A and B (sum is disjoint union)*/
MULTISET *MultisetSubtract(MULTISET *C, MULTISET *A, MULTISET *B); /* C = difference of A and B. Negative result -> error */

#endif
