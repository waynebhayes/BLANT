#include "align.h"

#include "bst.h"
#include "graph.h"
#include "hashmap.h"
#include "sets.h"

#include <stdlib.h>

void align_combinations_init(struct align_combinations* align) {
    align->buckets = hashmap_new();
}

int align_combinations_bucket_free(int k, bst_t* bst) {
    bst_free(bst);
    return MAP_OK;
}

void align_combinations_free(struct align_combinations* align) {
    hashmap_iterate(align->buckets, (PFany)align_combinations_bucket_free);
    hashmap_free(align->buckets);
}

void align_combinations_insert(struct align_combinations* align, int k, struct aligned_pairs alignment) {
    bst_t* bst;

    if (hashmap_get(align->buckets, k, (void**)&bst) != MAP_OK) {
        // insert new bucket into buckets
        bst = malloc(sizeof(bst_t));
        bst_init(bst);
        hashmap_put(align->buckets, k, bst);

        // TODO: what if it fails
    }
    
    bst_insert(bst, alignment);
}

unsigned long aligned_pair_index(unsigned int n1, unsigned int n2) {
    unsigned long n1_component = (unsigned long)n1 << sizeof(unsigned int);
    unsigned long n2_component = (unsigned long)n2;

    return n1_component | n2_component;
}

struct aligned_pair aligned_pair_from_index(unsigned long i) {
    unsigned int n1 = i >> sizeof(unsigned int);
    unsigned int n2 = i;

    struct aligned_pair pair;
    pair.n1 = n1;
    pair.n2 = n2;

    return pair;
}

struct aligned_pairs aligned_pairs_init(GRAPH* g1, GRAPH* g2) {
    unsigned long last_pair = aligned_pair_index(g1->n - 1, g2->n - 1);

    struct aligned_pairs alignment;
    alignment.pairs = SparseSetAlloc(last_pair);

    return alignment;
}

void aligned_pairs_free(struct aligned_pairs alignment) {
    SparseSetFree(alignment.pairs);
}

unsigned char aligned_pairs_contains(struct aligned_pairs alignment, unsigned int n1, unsigned int n2) {
    return SparseSetIn(alignment.pairs, aligned_pair_index(n1, n2));
}

void aligned_pairs_insert(struct aligned_pairs alignment, unsigned int n1, unsigned int n2) {
    unsigned long pair_index = aligned_pair_index(n1, n2);
    SparseSetAdd(alignment.pairs, pair_index);
}

unsigned long aligned_pairs_to_array(struct aligned_pairs alignment, struct aligned_pair** buf) {
    unsigned long n_pairs = SparseSetCardinality(alignment.pairs);
    struct aligned_pair* array = malloc(sizeof(struct aligned_pair) * n_pairs);

    unsigned long i;
    unsigned long p = 0;
    for (i = 0; i < alignment.pairs->n; i++) {
        if (SparseSetIn(alignment.pairs, i)) {
            array[p++] = aligned_pair_from_index(i);
        }
    }

    // because the pair indices are implicitly in sorted order, this array is guaranteed to be sorted
    *buf = array;
    return n_pairs;
}
