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
#include <math.h> // for sqrt(n) in SPARSE_SETs

unsigned setBits = sizeof(SETTYPE)*8, setBits_1;
static int SIZE(int n) { return (n+setBits-1)/setBits; }   /* number of array elements needed to store n bits */

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


/* Currently this just initializes lookupBitCount[].
** SetStartup doesn't perform startup more than once, so it's safe
** (and costs little) to call it again if you're not sure.  It returns
** 1 if it did the initialization, else 0.
*/
int SetStartup(void)
{
    assert(sizeof(SETTYPE) == 4);	// we assume 32-bit ints
    if(setBits_1)
	return 0;
    else
    {
	setBits_1 = setBits-1;
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
    if(!setBits_1) SetStartup();
    SET *set = (SET*) Calloc(1,sizeof(SET));
    set->n = n;
    set->smallestElement = n; // ie., invalid
    set->array = (SETTYPE*) Calloc(sizeof(SETTYPE), SIZE(n));
    return set;
}


SPARSE_SET *SparseSetAlloc(unsigned long n)
{
    SPARSE_SET *set = (SPARSE_SET*) Calloc(1,sizeof(SPARSE_SET));
    set->n = n;
    set->sqrt_n = ceil(sqrt(n));
    set->sets = (SET**) Calloc(set->sqrt_n,sizeof(SET*));
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
    int arrayElem=SIZE(set->n);
    set->smallestElement = set->n;
#if 1
    memset(set->array, 0, arrayElem * sizeof(set->array[0]));
#else
    int i;
    for(i=0; i<arrayElem; i++)
	set->array[i] = 0;
#endif
    return set;
}


/*
** erase all members from a set, but don't free it's memory.
*/
SPARSE_SET *SparseSetEmpty(SPARSE_SET *set)
{
    int i;
    for(i=0; i < set->sqrt_n; i++)
	if(set->sets[i])
	{
	    SetFree(set->sets[i]);
	    set->sets[i] = NULL;
	}
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


/* free all space occupied by a set
*/
void SparseSetFree(SPARSE_SET *set)
{
    int i;
    for(i=0; i < set->sqrt_n; i++)
	if(set->sets[i])
	{
	    SetFree(set->sets[i]);
	    set->sets[i] = NULL;
	}
    free(set);
}


/* Copy a set.  If the destination is NULL, it will be alloc'd
*/
SET *SetCopy(SET *dst, SET *src)
{
    int i, numSrc = SIZE(src->n);

    if(!dst)
	dst = SetAlloc(src->n);
    assert(dst->n == src->n);
    dst->smallestElement = src->smallestElement;

    for(i=0; i < numSrc; i++)
	dst->array[i] = src->array[i];
    return dst;
}

SPARSE_SET *SparseSetCopy(SPARSE_SET *dst, SPARSE_SET *src)
{
    int i;

    if(!dst)
	dst = SparseSetAlloc(src->n);
    assert(dst->n == src->n);

    for(i=0; i < dst->sqrt_n; i++)
	SetCopy(dst->sets[i], src->sets[i]);
    return dst;
}


/* Add an element to a set.  Returns the same set handle.
*/
SET *SetAdd(SET *set, unsigned element)
{
    assert(element < set->n);
    set->array[element/setBits] |= SET_BIT(element);
    if(element < set->smallestElement) set->smallestElement = element;
    return set;
}

SPARSE_SET *SparseSetAdd(SPARSE_SET *set, unsigned long element)
{
    assert(element < set->n);
    int which = element / set->sqrt_n;
    if(!set->sets[which])
	set->sets[which] = SetAlloc(set->sqrt_n);
    SetAdd(set->sets[which], element - which*set->sqrt_n);
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


unsigned int SetAssignSmallestElement1(SET *set)
{
    int i, old=set->smallestElement;
    assert(old == set->n || !SetIn(set, old)); // it should not be in there!

    for(i=0; i<set->n; i++)
        if(SetIn(set, i)) // the next smallest element is here
	       break;
    set->smallestElement = i; // note this works even if there was no new smallest element, so it's now set->n
    if(i == set->n)
	assert(SetCardinality(set) == 0);
    return i;
}

/* Delete an element from a set.  Returns the same set handle.
*/
SET *SetDelete(SET *set, unsigned element)
{
    assert(element < set->n);
    set->array[element/setBits] &= ~SET_BIT(element);
    if(element == set->smallestElement)
    {
	SetAssignSmallestElement1(set);
	assert(set->smallestElement > element);
    }
    return set;
}


SPARSE_SET *SparseSetDelete(SPARSE_SET *set, unsigned long element)
{
    assert(element < set->n);
    int which = element / set->sqrt_n;
    if(set->sets[which])
    {
	SetDelete(set->sets[which], element - which*set->sqrt_n);
	if(SetCardinality(set->sets[which]) == 0)
	{
	    SetFree(set->sets[which]);
	    set->sets[which] = NULL;
	}
    }
    return set;
}


/* query if an element is in a set; return 0 or non-zero.
*/
#define SET_BIT_SAFE(e) (1UL<<((e)%setBits))
Boolean SetInSafe(SET *set, unsigned element)
{
    assert(element < set->n);
    return (set->array[element/setBits] & SET_BIT_SAFE(element));
}

Boolean SparseSetIn(SPARSE_SET *set, unsigned long element)
{
    assert(element < set->n);
    int which = element / set->sqrt_n;
    return set->sets[which] && SetIn(set->sets[which], element - which*set->sqrt_n);
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
Boolean SparseSetEq(SPARSE_SET *A, SPARSE_SET *B)
{
    int i;
    assert(A->n == B->n);
    for(i=0; i < A->sqrt_n; i++)
	if(!SetEq(A->sets[i], B->sets[i]))
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
    C->smallestElement = MIN(A->smallestElement, B->smallestElement);
    return C;
}
SPARSE_SET *SparseSetUnion(SPARSE_SET *C, SPARSE_SET *A, SPARSE_SET *B)
{
    int i;
    assert(A->n == B->n && B->n == C->n);
    for(i=0; i < A->sqrt_n; i++)
	SetUnion(C->sets[i], A->sets[i], B->sets[i]);
    return C;
}


unsigned int SetAssignSmallestElement3(SET *C,SET *A,SET *B)
{
    if(A->smallestElement == B->smallestElement)
	C->smallestElement = A->smallestElement;
    else if(SetIn(A, B->smallestElement))
    {
	assert(!SetIn(B, A->smallestElement));
	assert(A->smallestElement < B->smallestElement);
	C->smallestElement = A->smallestElement;
    }
    else if(SetIn(B, A->smallestElement))
    {
	assert(!SetIn(A, B->smallestElement));
	assert(B->smallestElement < A->smallestElement);
	C->smallestElement = B->smallestElement;
    }
    else
    {
	SetAssignSmallestElement1(C);
	assert(C->smallestElement > A->smallestElement);
	assert(C->smallestElement > B->smallestElement);
    }
    return C->smallestElement;
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
    SetAssignSmallestElement3(C,A,B);
    return C;
}

SPARSE_SET *SparseSetIntersect(SPARSE_SET *C, SPARSE_SET *A, SPARSE_SET *B)
{
    int i;
    assert(A->n == B->n && B->n == C->n);
    for(i=0; i < C->sqrt_n; i++)
	SetIntersect(C->sets[i], A->sets[i], B->sets[i]);
    return C;
}


/* XOR A and B into C.  Any or all may be the same pointer.
*/
SET *SetXOR(SET *C, SET *A, SET *B)
{
    int i;
    int loop = SIZE(C->n);
    if(!A) return SetCopy(C, B);
    if(!B) return SetCopy(C, A);
    assert(A->n == B->n && B->n == C->n);
    for(i=0; i < loop; i++)
	C->array[i] = A->array[i] ^ B->array[i];
    SetAssignSmallestElement3(C,A,B);
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
    SetAssignSmallestElement1(B);
    return B;
}


unsigned SetCardinality(SET *A)
{
    unsigned n = 0, i, loop = SIZE(A->n);
    for(i=0; i < loop; i++)
	if(A->array[i]) n += SetCountBits(A->array[i]);
    return n;
}

unsigned long SparseSetCardinality(SPARSE_SET *set)
{
    unsigned int i;
    unsigned long sum=0;
    for(i=0; i < set->sqrt_n; i++)
	if(set->sets[i])
	    sum += SetCardinality(set->sets[i]);
    return sum;
}

/* populate the given array with the list of members currently present
** in the set.  The array is assumed to have enough space.
*/
unsigned SetToArray(unsigned int *array, SET *set)
{
    int pos = 0;
    int i;
    for(i=0; i < set->n; i++)
	if(SetIn(set,i))
	    array[pos++] = i;

    assert(pos == SetCardinality(set));
    return pos;
}

unsigned SSetToArray(unsigned int *array, SSET set)
{
    int pos = 0;
    int i;
    for(i=0; i < MAX_SSET; i++)
	if(SSetIn(set,i))
	    array[pos++] = i;

    assert(pos == SSetCountBits(set));
    return pos;
}


unsigned TSetToArray(unsigned int *array, TSET set)
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
SET *SetFromArray(SET *set, int n, unsigned int *array)
{
    while(n > 0)
	SetAdd(set, array[--n]);
    return set;
}

/* Add the elements listed in the array to the small set.
*/
SSET SSetFromArray(int n, unsigned int *array)
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

TSET TSetFromArray(int n, unsigned int *array)
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

void SetPrint(SET *A)
{
    int i;
    for(i=0;i<A->n;i++) if(SetIn(A,i)) printf("%d ", i);
    printf("\n");
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

