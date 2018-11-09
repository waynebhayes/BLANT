#ifndef _SORTS_H
#define _SORTS_H

#include <assert.h>
#include <alloca.h>
#include <stdlib.h> /* to ensure our definition of SortFcn matches theirs */

typedef int (*pfnCmpFcn)(const void*, const void*);

/*
** The sort function type returns the number times memcpy was called to move
** a data item during sorting. (Except QuickSort, which just calls the
** system qsort, so we don't know how many copies were performed.)
** Otherwise, the interface is exactly like the system qsort.
*/
typedef int SortFcn(void *, size_t nel, size_t width, pfnCmpFcn);

/* some sort funcitons. MergeSort is the only nlogn stable one. */
extern SortFcn
    QuickSort,	    /* nlogn, in-place.  I just call the stdlib qsort */
    CombSort,	    /* simple, in-place nlogn sort; see BYTE Apr 1991 */
    HeapSort,	    /* nlogn, in-place. */
    MergeSort,	    /* nlogn, but requires n*w extra space. */
    InsertionSort,  /* n^2 insertion sort, in-place. */
    PileSort,       /* Wayne's sorting method... HayesSort? */
    FredSort;       /* Fredrickson's O(n^2 log n) space*time algorithm */
extern void InsertionSortInt(int A[],size_t k);
/*
** PointerSort is useful for sorting large data items.  In that case, it's
** more efficient to create an array of pointers, sort the pointers, and
** then copy the elements back in one O(n) operation.  This is what
** PointerSort does, in a black-box fashion.  It uses 4*n bytes extra
** space for pointers.  You have to tell it which SortFcn above to use.
** (One of the above.)
*/
#if 1
#define PointerSort(Sort,a,n,w,compare) Sort(a,n,w,compare)
#else
extern int PointerSort(SortFcn, void *, size_t nel, size_t width, pfnCmpFcn);
#endif

#endif  /* _SORTS_H */
