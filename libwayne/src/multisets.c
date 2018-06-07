#include "multisets.h"
#include "misc.h"
#include <assert.h>
/* Allocates multiset with arraysize of n. Frequency allowed up to 1^sizeof(FREQTYPE) -1
*/
MULTISET *MultisetAlloc(unsigned int n) {
    MULTISET *mset = (MULTISET*) Calloc(1,sizeof(MULTISET));
    mset->n = n;
    mset->cardinality = 0;
    mset->array = (FREQTYPE*) Calloc(sizeof(FREQTYPE), n);
    return mset;
}

//MULTISET *MultisetResize(MULTISET *mset, unsigned int new_n); not implemented
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
    mset->cardinality = 0;
    return mset;
}

/* Returns the number of distinct elements in a multiset. Its cardinality.
*/
unsigned MultisetCardinality(MULTISET *mset) {
    return mset->cardinality;
}

/* Add an element to a multiset.  Returns the same set handle.
*/
MULTISET *MultisetAdd(MULTISET *mset, unsigned element) {
    assert(element < mset->n);
    if (mset->array[element] == 0) mset->cardinality++; //If there weren't any before cardinality goes up
    if (mset->array[element] == MAX_MULTISET_NUM) //Check for unsigned overflow before incrementing
        Fatal("Multiset attempted to incremement past limit: %d", MAX_MULTISET_NUM);
    else mset->array[element]++;
    return mset;
}

/* Remove an element from a multiset. If the frequency of that element reaches 0, cardinality is decremented.
   First checks for underflow
*/ 
MULTISET *MultisetDelete(MULTISET *mset, unsigned element) {
    assert(element < mset->n);
    if (mset->array[element] == 0) //Check for unsigned underflow before decrementing
        Fatal("Multiset attempted to decrement below 0");
    mset->array[element]--;
    if (mset->array[element] == 0) mset->cardinality--;
    return mset;
}
