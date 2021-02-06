#ifndef _SETS_H
#define _SETS_H
/* Blame this on Wayne Hayes, wayne@csri.utoronto.ca.
**
** simple implementation of sets.  You allocate a SET by specifying the max
** number N of things that could be in the set; each possible element is
** represented by an unsigned from 0..N-1, and it's presence in the set is
** represented by that bit in the array being a 1, else 0.
**
** You of course should not make your programs dependent on this implementation;
** you should act upon sets using *only* the defined functions.
**
** It is very space efficient, and reasonably time efficient, except for
** SetToArray in sets.c, which can be very slow for sparse sets inside a
** large SET structure.
**
** Any operations on more than one set *must* have both operands of the
** exact same size and type of set (this restriction will be laxed in
** later implementations).
*/

#include <stdlib.h>
#include <assert.h>
#include "misc.h"

typedef unsigned SETTYPE;
extern unsigned setBits, setBits_1;
//static unsigned SET_BIT(unsigned e) {assert((e%setBits) == (e&setBits_1)); return (1U<<((e)%setBits));}
//#define SET_BIT(e) (1U<<((e)%setBits))
#define SET_BIT(e) (1U<<((e)&setBits_1))

typedef struct _setType {
    unsigned n; /* in bits */
    unsigned smallestElement;
    SETTYPE* array;
} SET;

#define LOOKUP_NBITS 16
#define LOOKUP_SIZE (1 << LOOKUP_NBITS)
#define LOOKUP_MASK (LOOKUP_SIZE - 1)
extern unsigned lookupBitCount[LOOKUP_SIZE];
#define SetCountBits(i) \
    (lookupBitCount[((SETTYPE)(i)) & LOOKUP_MASK] + \
	lookupBitCount[(((SETTYPE)(i)) >> LOOKUP_NBITS) & LOOKUP_MASK])

Boolean SetStartup(void); // always succeeds, but returns whether it did anything or not.

/* allocate & return empty set capable of storing integers 0..n-1 inclusive */
SET *SetAlloc(unsigned n);
SET *SetResize(SET *s, unsigned new_n);
void SetFree(SET *set); /* free all memory used by a set */
SET *SetEmpty(SET *set);    /* make the set empty (set must be allocated )*/
#define SetReset SetEmpty
SET *SetCopy(SET *dst, SET *src);  /* if dst is NULL, it will be alloc'd */
SET *SetAdd(SET *set, unsigned element);    /* add single element to set */
SET *SetAddList(SET *set, ...); /* end list with (-1); uses varargs/stdarg */
SET *SetDelete(SET *set, unsigned element); /* delete a single element */
SET *SetUnion(SET *C, SET *A, SET *B);  /* C = union of A and B */
SET *SetIntersect(SET *C, SET *A, SET *B);  /* C = intersection of A and B */
SET *SetXOR(SET *C, SET *A, SET *B);  /* C = XOR of A and B */
SET *SetComplement(SET *B, SET *A);  /* B = complement of A */
unsigned SetCardinality(SET *A);    /* returns non-negative integer */
Boolean SetInSafe(SET *set, unsigned element); /* boolean: 0 or 1 */
#define SetSmallestElement(S) (S->smallestElement)
#if 1 //NDEBUG || !PARANOID_ASSERTS
// Note we do not check here if e is < set->n, which is dangerous
#define SetIn(set,e) ((set)->array[(e)/setBits] & SET_BIT(e))
//#define SetIn(set,e) ((e)>=0 && (e)<(set)->n && ((set)->array[(e)/setBits] & SET_BIT(e)))
#else
#define SetIn SetInSafe
#endif
Boolean SetEq(SET *set1, SET *set2);
Boolean SetSubsetEq(SET *sub, SET *super); /* is sub <= super? */
#define SetSupersetEq(spr,sb) SetSubsetEq((sb),(spr))
Boolean SetSubsetProper(SET *sub, SET *super);	/* proper subset */
#define SetSupersetProper(spr,sub) SetSubsetProper((sub),(spr))

/*
** You allocate an array big enough to hold the number of elements,
** and this function will populate it with the actual integers that are
** members of the set.  So if the set contains, say {0,17,324}, array
** will have array[0]=0, array[1]=17, array[2]=324, array[3..*] unchanged,
** and the function returns the cardinality (ie, number of array elements
** populated)
*/
unsigned SetToArray(unsigned *array, SET *set);
SET *SetFromArray(SET *s, int n, unsigned *array);
char *SetToString(int len, char s[], SET *set);

SET *SetPrimes(long n); /* return the set of all primes between 0 and n */
void SetPrint(SET *A); /* print elements of the set */

/*
*********  SPARSE_SET  ********
*/
typedef struct _sparseSetType {
    unsigned long n;
    unsigned sqrt_n;
    SET **sets;
} SPARSE_SET;
SPARSE_SET *SparseSetAlloc(unsigned long n);
void SparseSetFree(SPARSE_SET *set); /* free all memory used by a set */
SPARSE_SET *SparseSetEmpty(SPARSE_SET *set);    /* make the set empty (set must be allocated )*/
#define SetReset SetEmpty
SPARSE_SET *SparseSetCopy(SPARSE_SET *dst, SPARSE_SET *src);  /* if dst is NULL, it will be alloc'd */
SPARSE_SET *SparseSetAdd(SPARSE_SET *set, unsigned long element);    /* add single element to set */
SPARSE_SET *SparseSetAddList(SPARSE_SET *set, ...); /* end list with (-1); uses varargs/stdarg */
SPARSE_SET *SparseSetDelete(SPARSE_SET *set, unsigned long element); /* delete a single element */
SPARSE_SET *SparseSetUnion(SPARSE_SET *C, SPARSE_SET *A, SPARSE_SET *B);  /* C = union of A and B */
SPARSE_SET *SparseSetIntersect(SPARSE_SET *C, SPARSE_SET *A, SPARSE_SET *B);  /* C = intersection of A and B */
SPARSE_SET *SparseSetXOR(SPARSE_SET *C, SPARSE_SET *A, SPARSE_SET *B);  /* C = XOR of A and B */
SPARSE_SET *SparseSetComplement(SPARSE_SET *B, SPARSE_SET *A);  /* B = complement of A */
unsigned long SparseSetCardinality(SPARSE_SET *A);    /* returns non-negative integer */
Boolean SparseSetIn(SPARSE_SET *set, unsigned long element); /* boolean: 0 or 1 */
Boolean SparseSetEq(SPARSE_SET *set1, SPARSE_SET *set2);
Boolean SparseSetSubsetEq(SPARSE_SET *sub, SPARSE_SET *super); /* is sub <= super? */
#define SparseSetSupersetEq(spr,sb) SparseSetSubsetEq((sb),(spr))
Boolean SparseSetSubsetProper(SPARSE_SET *sub, SPARSE_SET *super);	/* proper subset */
#define SparseSetSupersetProper(spr,sub) SetSubsetProper((sub),(spr))


#define SMALL_SET_SIZE 64

#if SMALL_SET_SIZE == 128
    typedef __int128 SSET;
#elif SMALL_SET_SIZE == 64
    typedef unsigned long long SSET;    /* Small set */
#endif
#define SSET1 1ULL
#define SSET_NULLSET 0ULL
#define MAX_SSET (8*sizeof(SSET))

#define SSetEmpty(s) s = 0
#define SSetReset SSetEmpty
#define SSetAdd(s,e) (s |= (SSET1 << (e)))
#define SSetDelete(s,e) (s &= ~(SSET1 <<(e)))
#define SSetIn(s,e) ((s) & (SSET1 << (e)))
#define SSetEq(s1,s2) ((s1)==(s2))
#define SSetSubsetEq(sub,super) (((super)&(sub))==(sub))
#define SSetSupersetEq(a,b) SSetSubsetEq((b),(a))
#define SSetSubsetProper(sb,spr) (SSetSubsetEq((sb),(spr))&&!SSetEq((sb),(spr)))
#define SSetSupersetProper(spr,sb) SSetSubsetProper(sb,spr)
#define SSetUnion(a,b) ((a) | (b))
#define SSetIntersect(a,b) ((a) & (b))
#define SSetCountBits(s) (SetCountBits((s) & 0xffffffff) + SetCountBits((s) >> 32))
#define SSetCardinality SSetCountBits
SSET SSetFromArray(int n, unsigned *array);
unsigned SSetToArray(unsigned *array, SSET set);
char *SSetToString(int len, char s[], SSET set);

/* SSET dictionary - a set of sets */
typedef struct _ssetDict SSETDICT;
SSETDICT *SSetDictAlloc(int init_size);
SSETDICT *SSetDictAdd(SSETDICT*, SSET);
Boolean SSetDictIn(SSETDICT*, SSET);
void SSetDictFree(SSETDICT*);

#define TINY_SET_SIZE 16

#if TINY_SET_SIZE >= 64
    typedef SSET TSET;
    #define TSET1 SSET1
    #define TSET_NULLSET SSET_NULLSET
    #define MAX_TSET MAX_SSET

    #define TSetEmpty(s) SSetEmpty(s)
    #define TSetReset SSetEmpty
    #define TSetAdd(s,e) SSetAdd(s,e)
    #define TSetDelete(s,e) SSetDelete(s,e)
    #define TSetIn(s,e) SSetIn(s,e)
    #define TSetEq(s1,s2) SSetEq(s1,s2)
    #define TSetSubsetEq(sub,super) SSetSubsetEq(sub,super)
    #define TSetSupersetEq(a,b) SSetSupersetEq(a,b)
    #define TSetSubsetProper(sb,spr) SSetSubsetProper(sb,spr)
    #define TSetSupersetProper(spr,sb) SSetSupersetProper(spr,sb) 
    #define TSetUnion(a,b) SSetUnion(a,b)
    #define TSetIntersect(a,b) SSetIntersect(a,b)
    #define TSetCountBits(i) SSetCountBits(i) 
    #define TSetCardinality SSetCardinality 
    #define TSetFromArray SSetFromArray
    #define TSetToArray SSetToArray
    #define TSetToString SSetToString

#else
    #if TINY_SET_SIZE == 32
	typedef unsigned int TSET;
    #elif TINY_SET_SIZE == 16
	typedef unsigned short TSET;
    #elif TINY_SET_SIZE == 8
	typedef unsigned char TSET;
    #else
	#error unknown TINY_SET_SIZE
    #endif

    #define TSET1 ((TSET)1)
    #define TSET_NULLSET ((TSET)0)
    #define MAX_TSET (8*sizeof(TSET))

    #define TSetEmpty(s) s = 0
    #define TSetReset TSetEmpty
    #define TSetAdd(s,e) (s |= (TSET1 << (e)))
    #define TSetDelete(s,e) (s &= ~(TSET1 <<(e)))
    #define TSetIn(s,e) ((s) & (TSET1 << (e)))
    #define TSetEq(s1,s2) ((s1)==(s2))
    #define TSetSubsetEq(sub,super) (((super)&(sub))==(sub))
    #define TSetSupersetEq(a,b) TSetSubsetEq((b),(a))
    #define TSetSubsetProper(sb,spr) (TSetSubsetEq((sb),(spr))&&!TSetEq((sb),(spr)))
    #define TSetSupersetProper(spr,sb) TSetSubsetProper(sb,spr)
    #define TSetUnion(a,b) ((a) | (b))
    #define TSetIntersect(a,b) ((a) & (b))
    #define TSetCountBits(i) lookupBitCount[i]
    #define TSetCardinality TSetCountBits
    TSET TSetFromArray(int n, unsigned int *array);
    unsigned TSetToArray(unsigned int *array, TSET set);
    char *TSetToString(int len, char s[], TSET set);

#endif

#endif /* _SETS_H */
