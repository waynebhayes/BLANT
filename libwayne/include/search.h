#ifndef _SEARCH_H
#define _SEARCH_H

#include <assert.h>
#include <alloca.h>
#include <stdlib.h> /* to ensure our definition of SortFcn matches theirs */

typedef int (*pfnCmpFcn)(const void*, const void*);

typedef void SearchFcn(const void *key, const void *a, size_t nel,
    size_t width, pfnCmpFcn);

/* some sort funcitons. MergeSort is the only nlogn stable one. */
extern SearchFcn
    bsearch,	/* the ANSI standard one */
    BinarySearch;	/* My version */

#endif  /* _SEARCH_H */
