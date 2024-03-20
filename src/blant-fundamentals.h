// This file defines the fundamental compile-time constants that all blant-related programs must know.

#ifndef _BLANT_FUNDAMENTALS_H
#define _BLANT_FUNDAMENTALS_H

#define SELF_LOOPS 0	// do we allow self-loops? MUST be 0 or 1, no other values allowed.

// MAX_K is the maximum number of nodes in a graphlet that is supported by BLANT when using a fixed lookup table (as
// opposed to one that uses associaive arrays).  Maximum value is 7 with self-loops, 8 without.
#define MAX_K (8-SELF_LOOPS)

// maximum number of entries in the canon_map (lookup table)
#define maxBk (1 << (MAX_K*(MAX_K-1)/2 + MAX_K*SELF_LOOPS))

// BLANT represents a graphlet using one-half of the adjacency matrix (since we are assuming symmetric, undirected graphs)
// We have a choice of using the upper or lower triangle. We prefer the lower triangle because that's what Jesse uses
// (the graphlet / orbit generation code of Ine Melckenbeeck and friends Ghent university), and they published first.
// NOTE THAT IF YOU CHOOSE UPPER TRIANGLE THEN THE TESTS IN THE MAKEFILE WILL FAIL.
#define LOWER_TRIANGLE	1

// Once we find which canonical graphlet corresponds to a sampled graphlet, we want to know the permutation between the
// two.  We default to the permutation from the canonical to the sampled non-canonical; thus, when we list the nodes
// in the graphlet on BLANT's output, the ordering means that the columns of identical graphlets correspond to an exact
// local alignment between the nodes listed.  Listing them in the other direction doesn't seem to have much use.
#define PERMS_CAN2NON	1

// This compile-time constant defines whether or not we perform dynamic on-the-fly construction of the canon_map lookup
// table (stored in _K) rather than reading in canon_map/* files. The default (for now) is 0, meaning read in the files.
// It would be nice to get this working with the value 1 rather than 0.
#define DYNAMIC_CANON_MAP 0 // it kinda does work now but let's keep it off to be safe

#define DEFAULT_DIGITS 3 // 3 digits of precision by defalut

#endif // _BLANT_FUNDAMENTALS_H
