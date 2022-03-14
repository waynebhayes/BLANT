#ifndef CDIJKSTRA_BST_HEADER
#define CDIJKSTRA_BST_HEADER

#include "align.h"

struct bst_node {
    struct aligned_pair* pairs;
    unsigned long n_pairs;

    struct bst_node* left;
    struct bst_node* right;
};

struct bst_node* bst_node_init(struct aligned_pairs alignment);
void bst_node_free(struct bst_node* node);
int bst_node_cmp(struct bst_node* n1, struct bst_node* n2);

typedef struct {
    struct bst_node* root;
} bst_t;

void bst_init(bst_t* bst);
void bst_free(bst_t* bst);
int bst_insert(bst_t* bst, struct aligned_pairs alignment);

#endif
