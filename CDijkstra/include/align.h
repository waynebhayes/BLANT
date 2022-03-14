#ifndef CDIJKSTRA_ALIGN_HEADER
#define CDIJKSTRA_ALIGN_HEADER

#include "hashmap.h"
#include "sets.h"

#include <stdio.h>

struct align_combinations {
    hashmap_t buckets;
};

struct aligned_pairs {
    SPARSE_SET* pairs;
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
