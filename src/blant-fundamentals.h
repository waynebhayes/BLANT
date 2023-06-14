#ifndef _BLANT_FUNDAMENTALS_H
#define _BLANT_FUNDAMENTALS_H

#define SELF_LOOPS 1	// do we allow self-loops? MUST be 0 or 1, no other values allowed.

// This is the maximum graphlet size that BLANT supports when using a fixed lookup table. (Cannot be bigger than 8.)
// Currently only used to determine the amount of static memory to allocate.
#define MAX_K (8-SELF_LOOPS)
#define maxBk (1 << (MAX_K*(MAX_K-1)/2)) // maximum number of entries in the canon_map

#endif // _BLANT_FUNDAMENTALS_H
