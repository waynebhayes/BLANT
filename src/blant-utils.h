#ifndef BLANT_UTILS_H
#define BLANT_UTILS_H

#include "blant.h"

typedef struct nwhn { // this struct is needed for sorting
    int node;
    int heur;
    char *name;
} node_whn;

struct orbitpair_bits {
// NOTE that we store k-3 for k (thus k3), thus k=3 maps to 0, 8 maps to 5, and 10 maps to 7 (the highest value in 3 bits).
#if MAX_K >= 9
    unsigned long long k3:4, ordinal: 38; // Kimia: I changed this from 54 to 38 since even at k=12 38 is enough
    unsigned int o:4, p:4;  
#else 
    unsigned int k3:3, ordinal: 14; // Kimia: I changed this from 23 to 14, which is enough to store unsigned ints < 12346
    unsigned int o:3, p:3; // k is stored as k-3, o and p correctly go from 0 through 7 max
#endif
};

int orbitpair_cmp(long int a, long int b);
long int orbitpair_copy(long int src);

Boolean arrayIn(unsigned *arr, int size, int item);
void printIntArray(int* arr, int n, char* name);
void ExtractPerm(unsigned char perm[_k], int i);
void InvertPerm(unsigned char inverse[_k], unsigned const char perm[_k]);
TINY_GRAPH *TinyGraphInducedFromGraph(TINY_GRAPH *Gv, GRAPH *G, unsigned *Varray);
int asccompFunc(const foint i, const foint j);
int descompFunc(const void *a, const void *b);
int nwhn_des_alph_comp_func(const void *a, const void *b);
int nwhn_des_rev_comp_func(const void *a, const void *b);
int nwhn_asc_alph_comp_func(const void *a, const void *b);
int nwhn_asc_rev_comp_func(const void *a, const void *b);
void SetGlobalCanonMaps(void);
void LoadMagicTable(void);

#endif
