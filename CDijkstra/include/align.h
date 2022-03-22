#ifndef CDIJKSTRA_ALIGN_HEADER
#define CDIJKSTRA_ALIGN_HEADER

#include "hashmap.h"
#include "sets.h"

#include <stdio.h>

// Implementation of the align combinations support data structure: a mapping of k-value to ordered set of aligned node pairs. See bst.h/.c for details on the custom BST implementation used in the implementation.

struct align_combinations {
    hashmap_t buckets;
};

// Pairs are stored in the set via a "pair-packing" technique that reduces the capacity of the set from O(n1 * n2) to O(n1 + n2). Tldr the two node indices in the pair are concatenated on top of each other in an 8-byte unsigned int; see align.c for details.

struct aligned_pairs {
    SPARSE_SET* pairs;
    unsigned int n1;
    unsigned int n2;
};

struct aligned_pair {
    unsigned int n1;
    unsigned int n2;
};

void align_combinations_init(struct align_combinations* align);
void align_combinations_free(struct align_combinations* align);
void align_combinations_insert(struct align_combinations* align, int k, struct aligned_pairs alignment);

struct aligned_pairs aligned_pairs_init();
void aligned_pairs_free(struct aligned_pairs alignment);
unsigned char aligned_pairs_contains(struct aligned_pairs alignment, unsigned int n1, unsigned int n2);
void aligned_pairs_insert(struct aligned_pairs alignment, unsigned int n1, unsigned int n2);
unsigned long aligned_pairs_to_array(struct aligned_pairs alignment, struct aligned_pair** buf);

// align combinations is an array of buckets
    // each bucket is a mapping from k -> bst
        // each bst contains nodes
            // each node contains a set of aligned pairs

#endif
