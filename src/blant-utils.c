#include <sys/file.h>
#include <sys/mman.h>
#include <math.h>
#include "blant.h"
#include "blant-utils.h"

// The following is the most compact way to store the permutation between a non-canonical and its canonical representative,
// when k=8: there are 8 entries, and each entry is a integer from 0 to 7, which requires 3 bits. 8*3=24 bits total.
// For simplicity we use the same 3 bits per entry, and assume 8 entries, even for k<8.  It wastes memory for k<4, but
// makes the coding much simpler.
typedef unsigned char kperm[3]; // The 24 bits are stored in 3 unsigned chars.

// Here's where we're lazy on saving memory, and we could do better.  We're going to allocate a static array
// that is big enough for the 256 million permutations from non-canonicals to canonicals for k=8, even if k<8.
// So we're allocating 256MBx3=768MB even if we need much less.  I figure anything less than 1GB isn't a big deal
// these days. It needs to be aligned to a page boundary since we're going to mmap the binary file into this array.
//static kperm Permutations[maxBk] __attribute__ ((aligned (8192)));
kperm *Permutations = NULL; // Allocating memory dynamically

static int _magicTable[MAX_CANONICALS][12]; //Number of canonicals for k=8 by number of columns in magic table

unsigned int L_K_Func(Gint_type Gint) {Apology("L_K_Func() not yet implemented");}

// Assuming the global variable _k is set properly, go read in and/or mmap the big global
// arrays related to canonical mappings and permutations.
void SetGlobalCanonMaps(void)
{
    int i;
    char BUF[BUFSIZ];
    assert(3 <= _k && _k <= 8);
    _Bk = (1 <<(_k*(_k-1)/2));
    _connectedCanonicals = canonListPopulate(BUF, _canonList, _k);
    _numCanon = _connectedCanonicals->n;
    _numConnectedCanon = SetCardinality(_connectedCanonicals);
    _numOrbits = orbitListPopulate(BUF, _orbitList, _orbitCanonMapping, _orbitCanonNodeMapping, _numCanon, _k);
    _K = (short*) mapCanonMap(BUF, _K, _k);

    sprintf(BUF, "%s/%s/perm_map%d.bin", _BLANT_DIR, CANON_DIR, _k);
    int pfd = open(BUF, 0*O_RDONLY);
    Permutations = (kperm*) mmap(Permutations, sizeof(kperm)*_Bk, PROT_READ, MAP_PRIVATE, pfd, 0);
    assert(Permutations != MAP_FAILED);
    _numConnectedOrbits = 0;
    for (i=0; i < _numOrbits; i++)
	if (SetIn(_connectedCanonicals, _orbitCanonMapping[i]))
	    _connectedOrbits[_numConnectedOrbits++] = i;
    if ((_outputMode == outputODV && _k == 4) || _k == 5)
	orcaOrbitMappingPopulate(BUF, _orca_orbit_mapping, _k);
}

void LoadMagicTable()
{
    int i,j;
    char BUF[BUFSIZ];
    sprintf(BUF, "%s/orca_jesse_blant_table/UpperToLower%d.txt", _BLANT_DIR, _k);
    FILE *fp_ord=fopen(BUF, "r");
    if(!fp_ord) Fatal("cannot find %s\n", BUF);
    for(i=0; i<_numCanon; i++) {
        for(j=0; j<12 ;j++) {
            assert(1==fscanf(fp_ord, "%d", &_magicTable[i][j]));
        }
    }
    fclose(fp_ord);
    if(_displayMode == orca || _displayMode == jesse) {
        // Load mapping from lower_ordinal to ORCA/Jesse ID into table for fast lookup
	for (i=0; i < _numCanon; i++) _outputMapping[_magicTable[i][4]] = _magicTable[i][7];
    }
}

// You provide a permutation array, we fill it with the permutation extracted from the compressed Permutation mapping.
// There is the inverse transformation, called "EncodePerm", in createBinData.c.
void ExtractPerm(char perm[_k], int i)
{
    int j, i32 = 0;
    for(j=0;j<3;j++) i32 |= (Permutations[i][j] << j*8);
    for(j=0;j<_k;j++)
	perm[j] = (i32 >> 3*j) & 7;
}

void InvertPerm(char inv[_k], const char perm[_k])
{
    int j;
    for(j=0; j<_k; j++)
	inv[(int)perm[j]]=j;
}

Boolean arrayIn(int* arr, int size, int item) {
    int i;
    for(i=0; i<size; i++) {
        if (arr[i] == item) {
            return true;
        }
    }
    return false;
}

void printIntArray(int* arr, int n, char* name) {
    fprintf(stderr, "%s:");
    int i;

    for (i = 0; i < n; i++) {
        fprintf(stderr, " %d", arr[i]);
    }

    fprintf(stderr, "\n");
}

int asccompFunc(const foint i, const foint j) {return i.i > j.i ? 1 : i.i == j.i ? 0 : -1;} // ascending (smallest at top)

int descompFunc(const void *a, const void *b)
{
    int *i = (int*)a, *j = (int*)b;
    return (*j)-(*i);
}

int nwhn_des_alph_comp_func(const void *a, const void *b) {
    node_whn *nwhn_a = (node_whn*)a;
    node_whn *nwhn_b = (node_whn*)b;
    int i = nwhn_a->heur, j = nwhn_b->heur;

    if (j != i) {
        return j - i;
    } else {
        return strcmp(nwhn_a->name, nwhn_b->name);
    }
}

int nwhn_des_rev_comp_func(const void *a, const void *b) {
    node_whn *nwhn_a = (node_whn*)a;
    node_whn *nwhn_b = (node_whn*)b;
    int i = nwhn_a->heur, j = nwhn_b->heur;

    if (j != i) {
        return j - i;
    } else {
        return -1 * strcmp(nwhn_a->name, nwhn_b->name);
    }
}

int nwhn_asc_alph_comp_func(const void *a, const void *b) {
    node_whn *nwhn_a = (node_whn*)a;
    node_whn *nwhn_b = (node_whn*)b;
    int i = nwhn_a->heur, j = nwhn_b->heur;

    if (j != i) {
        return i - j;
    } else {
        return strcmp(nwhn_a->name, nwhn_b->name);
    }
}

int nwhn_asc_rev_comp_func(const void *a, const void *b) {
    node_whn *nwhn_a = (node_whn*)a;
    node_whn *nwhn_b = (node_whn*)b;
    int i = nwhn_a->heur, j = nwhn_b->heur;

    if (j != i) {
        return i - j;
    } else {
        return -1 * strcmp(nwhn_a->name, nwhn_b->name);
    }
}

// Given the big graph G and a set of nodes in V, return the TINY_GRAPH created from the induced subgraph of V on G.
TINY_GRAPH *TinyGraphInducedFromGraph(TINY_GRAPH *Gv, GRAPH *G, int *Varray)
{
    unsigned i, j;
    TinyGraphEdgesAllDelete(Gv);
    for(i=0; i < Gv->n; i++) for(j=i+1; j < Gv->n; j++)
        if(GraphAreConnected(G, Varray[i], Varray[j]))
            TinyGraphConnect(Gv, i, j);
    return Gv;
}

int getMaximumIntNumber(int K)
{
    assert(K >= 3 && K <= 8);
    int num_of_bits = K * (K-1) / 2;
    return pow(2, num_of_bits);
}
