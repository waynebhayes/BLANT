#include "align.h"
#include "bst.h"

#include <stdlib.h>

struct bst_node* bst_node_init(struct aligned_pairs alignment) {
    struct bst_node* node = malloc(sizeof(struct bst_node));
    node->n_pairs = aligned_pairs_to_array(alignment, &node->pairs);
    node->left = NULL;
    node->right = NULL;

    return node;
}

void bst_node_free(struct bst_node* node) {
    bst_node_free(node->left);
    bst_node_free(node->right);

    free(node->pairs);
    free(node);
}

int bst_node_cmp(struct bst_node* n1, struct bst_node* n2) {
    if (n1 == NULL || n2 == NULL) {
        // this is an artifact from the original Python version:
        // if either node is null, it's equivalent to equality
        return 0; 
    }

    int i;
    for (i = 0; i < n1->n_pairs; i++) {
        struct aligned_pair p1 = n1->pairs[i];
        struct aligned_pair p2 = n2->pairs[i];

        if (p1.n1 < p2.n1) {
            return -1;
        } else if (p1.n1 > p2.n1) {
            return 1;
        }
    }

    for (i = 0; i < n1->n_pairs; i++) {
        struct aligned_pair p1 = n1->pairs[i];
        struct aligned_pair p2 = n2->pairs[i];

        if (p1.n2 < p2.n2) {
            return -1;
        } else if (p1.n2 > p2.n2) {
            return 1;
        }
    }

    return 0;
}

void bst_init(bst_t* bst) {
    bst->root = NULL;
}

void bst_free(bst_t* bst) {
    bst_node_free(bst->root);
}

int bst_insert(bst_t* bst, struct aligned_pairs alignment) {
    struct bst_node* new_node = bst_node_init(alignment);

    if (bst->root == NULL) {
        bst->root = new_node;
        return 1;
    } else {
        struct bst_node* current = bst->root;

        while (1) {
            switch (bst_node_cmp(current, new_node)) {
            case 1:
                if (current->right != NULL) {
                    current = current->right;
                } else {
                    current->right = new_node;
                    return 1;
                }

                break;

            case -1:
                if (current->left != NULL) {
                    current = current->left;
                } else {
                    current->left = new_node;
                    return 1;
                }

                break;

            case 0: // make it explicit, but still exhaustive
            default:
                return 0;
            }
        }
    }
}
