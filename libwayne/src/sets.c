/* Version 0.0
** From "Wayne's Little DSA Library" (DSA == Data Structures and
** Algorithms) Feel free to change, modify, or burn these sources, but if
** you modify them please don't distribute your changes without a clear
** indication that you've done so.  If you think your change is spiffy,
** send it to me and maybe I'll include it in the next release.
**
** Wayne Hayes, wayne@cs.toronto.edu
*/

/* reasonably opaque implementation of bit fields, or sets of integers, by
** Wayne Hayes, wayne@cs.toronto.edu.
*/

#include "sets.h"
#include <stdio.h>
#include <stdarg.h>

static unsigned setBits = sizeof(SETTYPE)*8;

static int SIZE(int n)
{ return (n+setBits-1)/setBits; }   /* number of array elements needed to store n bits */


unsigned int lookupBitCount[LOOKUP_SIZE];
/* count the number of 1 bits in a long
*/
static unsigned DumbCountBits(unsigned long i)
{
    unsigned n = 0;
    while(i)
    {
	if(i&1) ++n;
	i >>= 1;
    }
    return n;
}


int SetStartup(void)
{
    assert(sizeof(SETTYPE) == 4);	// we assume 32-bit ints :-(
    if(lookupBitCount[1])
	return 0;
    else
    {
	unsigned long i;
	for(i=0; i<LOOKUP_SIZE; i++)
	    lookupBitCount[i] = DumbCountBits(i);
	return 1;
    }
}


/*
** SetAlloc: create a new empty set of max size n elements,
** and return its handle.
*/
SET *SetAlloc(unsigned n)
{
    SET *set = (SET*) Calloc(1,sizeof(SET));
    set->n = n;
    set->array = (SETTYPE*) Calloc(sizeof(SETTYPE), SIZE(n));
    return set;
}


/*
** SetResize: re-size a set.
*/
SET *SetResize(SET *set, unsigned new_n)
{
    int i, old_n = set->n;
    set->array = (SETTYPE*) Realloc(set->array, sizeof(SETTYPE) * SIZE(new_n));
    set->n = new_n;
    for(i=old_n; i < new_n; i++)
	SetDelete(set, i);
    return set;
}


/*
** erase all members from a set, but don't free it's memory.
*/
SET *SetEmpty(SET *set)
{
    int i, arrayElem=SIZE(set->n);
#if 1
    memset(set->array, 0, arrayElem * sizeof(set->array[0]));
#else
    for(i=0; i<arrayElem; i++)
	set->array[i] = 0;
#endif
    return set;
}


/* free all space occupied by a set
*/
void SetFree(SET *set)
{
    if(set)
    {
	if(set->array)
	    free(set->array);
	free(set);
    }
}


/* Copy a set.  If the destination is NULL, it will be alloc'd
*/
SET *SetCopy(SET *dst, SET *src)
{
    int i, numSrc = SIZE(src->n);

    if(!dst)
	dst = SetAlloc(src->n);
    assert(dst->n == src->n);

    for(i=0; i < numSrc; i++)
	dst->array[i] = src->array[i];
    return dst;
}


/* Add an element to a set.  Returns the same set handle.
*/
SET *SetAdd(SET *set, unsigned element)
{
    assert(element < set->n);
    set->array[element/setBits] |= (1UL << (element % setBits));
    return set;
}


/* Add a bunch of elements to a set.  End the list with (-1).
*/
SET *SetAddList(SET *set, ...)
{
#ifdef sgi
    Apology("stdarg doesn't work on the sgi");
    return NULL;
#else
    int e;
    va_list argptr;
    va_start(argptr, set);

    while((e = va_arg(argptr, int)) != -1)
	SetAdd(set, (unsigned)e);

    va_end(argptr);
    return set;
#endif
}


/* Delete an element from a set.  Returns the same set handle.
*/
SET *SetDelete(SET *set, unsigned element)
{
    assert(element < set->n);
    set->array[element/setBits] &= ~(1 << (element % setBits));
    return set;
}


/* query if an element is in a set; return 0 or non-zero.
*/
Boolean SetIn(SET *set, unsigned element)
{
    assert(element < set->n);
    if(set->array[element/setBits] & (1UL << (element % setBits)))
	return true;
    else
	return false;
}


/* See if A and B are the same set.
*/
Boolean SetEq(SET *A, SET *B)
{
    int i;
    int loop = SIZE(A->n);
    assert(A->n == B->n);
    for(i=0; i < loop; i++)
	if(A->array[i] != B->array[i])
	    return false;
    return true;
}


/* See if A is a (non-proper) subset of B.
*/
Boolean SetSubsetEq(SET *A, SET *B)
{
    int i;
    int loop = SIZE(A->n);
    assert(A->n == B->n);
    for(i=0; i < loop; i++)
	if((A->array[i] & B->array[i]) != A->array[i])
	    return false;
    return true;
}

Boolean SetSubsetProper(SET *A, SET *B)
{
    return SetSubsetEq(A,B) && !SetEq(A,B);
}


/* Union A and B into C.  Any or all may be the same pointer.
*/
SET *SetUnion(SET *C, SET *A, SET *B)
{
    int i;
    int loop = SIZE(C->n);
    assert(A->n == B->n && B->n == C->n);
    for(i=0; i < loop; i++)
	C->array[i] = A->array[i] | B->array[i];
    return C;
}


/* Intersection A and B into C.  Any or all may be the same pointer.
*/
SET *SetIntersect(SET *C, SET *A, SET *B)
{
    int i;
    int loop = SIZE(C->n);
    assert(A->n == B->n && B->n == C->n);
    for(i=0; i < loop; i++)
	C->array[i] = A->array[i] & B->array[i];
    return C;
}


/* Complement of A.  Both may be the same pointer.
*/
SET *SetComplement(SET *B, SET *A)
{
    int i;
    int loop = SIZE(B->n);
    assert(A->n == B->n);
    for(i=0; i < loop; i++)
	B->array[i] = ~A->array[i];
    return B;
}


unsigned SetCardinality(SET *A)
{
    unsigned n = 0, i, loop = SIZE(A->n);
    for(i=0; i < loop; i++)
	if(A->array[i]) n += SetCountBits(A->array[i]);
    return n;
}

/* populate the given array with the list of members currently present
** in the set.  The array is assumed to have enough space.
*/
unsigned SetToArray(unsigned *array, SET *set)
{
    int pos = 0;
    int i;
    for(i=0; i < set->n; i++)
	if(SetIn(set,i))
	    array[pos++] = i;

    assert(pos == SetCardinality(set));
    return pos;
}

unsigned SSetToArray(unsigned *array, SSET set)
{
    int pos = 0;
    int i;
    for(i=0; i < MAX_SSET; i++)
	if(SSetIn(set,i))
	    array[pos++] = i;

    assert(pos == SSetCountBits(set));
    return pos;
}


unsigned TSetToArray(unsigned *array, TSET set)
{
    int pos = 0;
    int i;
    for(i=0; i < MAX_TSET; i++)
	if(TSetIn(set,i))
	    array[pos++] = i;

    assert(pos == TSetCountBits(set));
    return pos;
}


/* Add the elements listed in the array to the set.
*/
SET *SetFromArray(SET *set, int n, unsigned *array)
{
    while(n > 0)
	SetAdd(set, array[--n]);
    return set;
}

/* Add the elements listed in the array to the small set.
*/
SSET SSetFromArray(int n, unsigned *array)
{
    SSET set;
    SSetEmpty(set);
    while(n > 0)
	SSetAdd(set, array[--n]);
    return set;
}

char *SSetToString(int len, char s[], SSET set)
{
    int i;
    for(i=0; i<len-1; i++)
	s[i] = '0' + !!SSetIn(set, i);
    s[len-1] = '\0';
    return s;
}

TSET TSetFromArray(int n, unsigned *array)
{
    TSET set;
    TSetEmpty(set);
    while(n > 0)
	TSetAdd(set, array[--n]);
    return set;
}

char *TSetToString(int len, char s[], TSET set)
{
    int i;
    for(i=0; i<len-1; i++)
	s[i] = '0' + !!TSetIn(set, i);
    s[len-1] = '\0';
    return s;
}

char *SetToString(int len, char s[], SET *set)
{
    int i;
    assert(len > set->n); /* need space for trailing '\0' */
    for(i=0; i<MIN(len, set->n); i++)
	s[i] = '0' + !!SetIn(set, i);
    s[i] = '\0';
    return s;
}

/* Use a seive to get all primes between 0 and n */
SET *SetPrimes(long n)
{
    SET *primes = SetAlloc(n+1);
    int i, loop=SIZE(n+1), p;

    for(i=0; i<loop; i++)
	--primes->array[i];     /* turn on all the bits */
    SetDelete(primes, 0);
    SetDelete(primes, 1);

    p=2;
    while(p <= n)
    {
	for(i = p+p; i <= n; i+=p) /* delete all multiples of p */
	    SetDelete(primes, i);
	/* find next prime */
	do  ++p;
	while(p <= n && !SetIn(primes, p));
    }
    return primes;
}


/*
** SSETDICT: a set of sets.  Idea is to be able to store sets and quickly
** query if a set is in the dictionary yet.  Used mostly by the circulant
** graph generation routines which need to quickly see if a circulant has
** already been generated yet.  Each set is stored only once -- adding
** an already existing set does nothing.  At allocation time, you need to
** know about how many items will eventually be stored.
** Implementation is an array of hash tables all of the same size.
** Collisions are resolved by putting it in the next hash table "down",
** and if you run out of hash tables then you move to the first hash
** table's next position, etc.
*/

/*
** Number of "stacked" hash tables.  Experiments show that as long as
** the initial allocation is a good estimate of the actual number of
** things eventually to be inserted, then there's no gain to have more
** than two.
*/
#define NUM_HASHES 2

struct _ssetDict {
    SSET *array[NUM_HASHES];
    int nCols;	/* number of columns - always prime */
    int nElem;	/* numbef of elements currently stored, not including NULLSET */
    Boolean containsNull; /* never explicitly store the empty Set */
};


SSETDICT *SSetDictAlloc(int n)
{
    int i;
    SSETDICT *ssd = Calloc(sizeof(SSETDICT), 1);
    assert(n>=1);
    ssd->nCols = n;
    ssd->nElem = 0;
    for(i=0; i< NUM_HASHES; i++)
	ssd->array[i] = Calloc(sizeof(SSET), ssd->nCols);
    return ssd;
}

/*
** Finds a place for ss; if ss is in the dictionary, it returns a pointer
** to ss's position; otherwise it returns a pointer to a blank position
** where ss can be put.
** Returns NULL if the dictionary is full.
*/
static SSET *SSetDict_find_place(SSETDICT *ssd, SSET ss)
{
    int position = ss % ssd->nCols, col;
    assert(ss != SSET_NULLSET);

    for(col=position; col != (position+ssd->nCols-1) % ssd->nCols; col = (col+1) % ssd->nCols)
    {
	int itable;
	for(itable=0; itable<NUM_HASHES; itable++)
	{
	    SSET *place = &(ssd->array[itable][col]);
	    if(*place == ss || *place == SSET_NULLSET)
		return place;
	}
    }

    return NULL;
}


/* Double the size of the dictionary */
static SSETDICT *SSetDictDouble(SSETDICT *ssd)
{
    int col, itable, nAdded = 0;
    SSETDICT *newDict = SSetDictAlloc(2*ssd->nCols);

    for(itable=0; itable<NUM_HASHES; itable++)
    {
	for(col=0; col < ssd->nCols; col++)
	    if(ssd->array[itable][col])
	    {
		SSetDictAdd(newDict, ssd->array[itable][col]);
		++nAdded;
	    }
	Free(ssd->array[itable]);
	ssd->array[itable] = newDict->array[itable];
    }
    assert(nAdded == ssd->nElem);
    
    ssd->nCols = newDict->nCols;
    Free(newDict);

    return ssd;
}


SSETDICT *SSetDictAdd(SSETDICT *ssd, SSET ss)
{
    SSET *place;

    if(ss == SSET_NULLSET)
    {
	ssd->containsNull = true;
	return ssd;
    }

    if(ssd->nElem == ssd->nCols)
	SSetDictDouble(ssd);

    place = SSetDict_find_place(ssd, ss);

    assert(place);

    if(*place != ss)
    {
	assert(*place == SSET_NULLSET);
	ssd->nElem++;
	*place = ss;	/* it's either already ss or NULLSET; assign it either way */
    }
		
    return ssd;
}

Boolean SSetDictIn(SSETDICT *ssd, SSET ss)
{
    SSET *place;
    if(ss == SSET_NULLSET)
	return ssd->containsNull;
    place = SSetDict_find_place(ssd, ss);
    if(place && *place == ss)
	return true;
    else
	return false;
}

void SSetDictFree(SSETDICT *ssd)
{
    int i;
    for(i=0; i<NUM_HASHES; i++)
	Free(ssd->array[i]);
    Free(ssd);
}

