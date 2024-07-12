#ifndef BLANT_UTILS_H
#define BLANT_UTILS_H

#include "blant.h"

typedef struct nwhn { // this struct is needed for sorting
    int node;
    int heur;
    char *name;
} node_whn;

struct orbitpair_bits {
    #if K_GE_9  
    unsigned int k : 4;
	unsigned int ordinal: 52;
    unsigned int o : 4;
    unsigned int p : 4;  
    #else 
    unsigned int k : 3;
	unsigned int ordinal: 23;
    unsigned int o : 3;
    unsigned int p : 3;
    #endif
};

int orbitpair_cmp(long int a, long int b);
long int orbitpair_copy(long int src);

Boolean arrayIn(unsigned *arr, int size, int item);
void printIntArray(int* arr, int n, char* name);
void ExtractPerm(unsigned char perm[_k], int i);
void InvertPerm(unsigned char inverse[_k], unsigned const char perm[_k]);
TINY_GRAPH *TinyGraphInducedFromGraph(TINY_GRAPH *Gv, GRAPH *G, unsigned *Varray);
int getMaximumIntNumber(int K);
int asccompFunc(const foint i, const foint j);
int descompFunc(const void *a, const void *b);
int nwhn_des_alph_comp_func(const void *a, const void *b);
int nwhn_des_rev_comp_func(const void *a, const void *b);
int nwhn_asc_alph_comp_func(const void *a, const void *b);
int nwhn_asc_rev_comp_func(const void *a, const void *b);
void SetGlobalCanonMaps(void);
void LoadMagicTable(void);

#endif
