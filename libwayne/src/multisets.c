#include "multisets.h"
#include "misc.h"
#include <assert.h>
/* Allocates multiset with arraysize of n. Multiplicity allowed up to 1^sizeof(FREQTYPE) -1 (255)
*/
MULTISET *MultisetAlloc(unsigned int n) {
    MULTISET *mset = (MULTISET*) Calloc(1,sizeof(MULTISET));
    mset->n = n;
    mset->support = 0;
    mset->array = (FREQTYPE*) Calloc(sizeof(FREQTYPE), n);
    return mset;
}

/* Resizes a multiset array to a new size;
   Recalculates the support in case the size is smaller and elements are removed.
*/
MULTISET *MultisetResize(MULTISET *mset, unsigned int new_n) {
    unsigned char* newArray = (FREQTYPE*) Calloc(sizeof(FREQTYPE), new_n);
    int loop = MIN(mset->n, new_n);
    int i;
    mset->support = 0;
    for (i = 0; i < loop; i++) {
        newArray[i] = mset->array[i];
        if (newArray[i] > 0) mset->support++;
    }
    unsigned char* temp = mset->array;
    mset->array = newArray;
    free(temp);
}

void MultisetFree(MULTISET *mset) {
    if(mset)
    {
	if(mset->array)
	    free(mset->array);
	free(mset);
    }
}

/*
** erase all members from a multiset, but don't free it's memory.
*/
MULTISET *MultisetEmpty(MULTISET *mset) {
#if 1
    memset(mset->array, 0, mset->n * sizeof(mset->array[0]));
#else
    int i;
    for(i=0; i<arrayElem; i++)
	set->array[i] = 0;
#endif
    mset->support = 0;
    return mset;
}

/* Returns the number of distinct elements in a multiset. Its cardinality.
*/
unsigned MultisetSupport(MULTISET *mset) {
    return mset->support;
}

/*  Returns the multiplicity of a single element
*/
unsigned char MultisetMultiplicity(MULTISET *mset, unsigned element) {
    assert(element < mset->n);
    return mset->array[element];
}

/* Add an element to a multiset.  Returns the same set handle.
   Adds one to the multiplicity of the element and updates the cardinality if necessary.
*/
MULTISET *MultisetAdd(MULTISET *mset, unsigned element) {
    assert(element < mset->n);
    if (mset->array[element] == 0) mset->support++; //If there weren't any before cardinality goes up
    if (mset->array[element] == MAX_MULTISET_NUM) //Check for unsigned overflow before incrementing
        Fatal("Multiset attempted to incremement past limit: %d", MAX_MULTISET_NUM);
    else mset->array[element]++;
    return mset;
}

/* Remove an element from a multiset. If the multiplicity of that element reaches 0, cardinality is decremented.
   First checks for underflow
*/ 
MULTISET *MultisetDelete(MULTISET *mset, unsigned element) {
    assert(element < mset->n);
    if (mset->array[element] == 0) //Check for unsigned underflow before decrementing
        Fatal("Multiset attempted to decrement below 0");
    mset->array[element]--;
    if (mset->array[element] == 0) mset->support--;
    return mset;
}

/* The sum of two multisets is the disjoint union of the two sets. Comparable to SetUnion.
   This function assumes the multiset C allocated and clears previous values.
   A B and C should have the same size.
*/
MULTISET *MultisetSum(MULTISET *C, MULTISET *A, MULTISET *B) {
    int i, freq;
    int loop = SIZE(C->n);
    assert(A->n == B->n && B->n == C->n);
    C->support = 0;
    for(i=0; i < loop; i++) {
        freq = A->array[i] + B->array[i]; //An int can safely store two unsigned chars added
        if (freq > MAX_MULTISET_NUM)
            FATAL("MultisetSum overflow. Limit is %d", MAX_MULTISET_NUM);
        C->array[i] = (unsigned char) freq;
        if (freq > 0) C->support++;
    }
    return C;
}

/* The difference of two multisets is the difference of their multiplicities.
   Negative differences become 0. This function assumes the multiset C allocated and clears previous values.
   A B and C should have the same size.
*/
MULTISET *MultisetSubtract(MULTISET *C, MULTISET *A, MULTISET *B) {
    int i, freq;
    int loop = SIZE(C->n);
    assert(A->n == B->n && B->n == C->n);
    C->support = 0;
    for(i=0; i < loop; i++) {
        freq = A->array[i] - B->array[i]; //An int can safely store two unsigned chars subtracted
        if (freq < 0) freq = 0;
        C->array[i] = (unsigned char) freq;
        if (freq > 0) C->support++;
    }
    return C;
}
