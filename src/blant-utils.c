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

int asccompFunc(const foint i, const foint j) {return i.i > j.i ? 1 : i.i == j.i ? 0 : -1;} // ascending (smallest at top)

int descompFunc(const void *a, const void *b)
{
    int *i = (int*)a, *j = (int*)b;
    return (*j)-(*i);
}

int nwd_descompFunc(const void *a, const void *b) {
    int i = ((node_wdegree*)a)->degree_val, j = ((node_wdegree*)b)->degree_val;
    return j - i;
}

int nwd_asccompFunc(const void *a, const void *b) {
    int i = ((node_wdegree*)a)->degree_val, j = ((node_wdegree*)b)->degree_val;
    return i - j;
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

int* enumerateDegreeOrder(GRAPH *G) {
    node_wdegree orderArray[G->n];
    int i;

    for (i = 0; i < G->n; ++i) {
        orderArray[i].node = i;
    }

    enumerateDegreeOrderHelper(G, orderArray, 0, G->n, 0);
    int* degreeOrder = malloc(sizeof(int) * G->n);

    for (i = 0; i < G->n; ++i) {
        degreeOrder[orderArray[i].node] = i;
    }

    return degreeOrder;
}

void enumerateDegreeOrderHelper(GRAPH *G, node_wdegree* orderArray, int start, int end, int layer) {
    // this is a temporary solution that will stop this function if it's called on two nodes with identical neighbors. this is pretty common since there are a lot of nodes of degree 1 and if both of them have the same neighbor, the algorithm will go through the entire graph before declaring them a tie. this just does that earlier while barely missing out on any real ties. I'm leaving this temporary solution here because the whole thing is temporary; I'm gonna test out many different methods of sorting the array
    if (layer >= 3) {
        return;
    }

    // fill in degree vals in orderArray based on the layer, from start to end
    SET *old_nodes = SetAlloc(G->n);
    SET *new_nodes = SetAlloc(G->n);
    SET *all_nodes = SetAlloc(G->n);
    int curr_i, curr_node, curr_layer, old_i, old_node, neigh_i, neigh_node, layer_deg_sum, i;

    for (curr_i = start; curr_i < end; ++curr_i) {
        curr_node = orderArray[curr_i].node;
        SetAdd(old_nodes, curr_node);
        SetAdd(all_nodes, curr_node);

        for (curr_layer = 0; curr_layer < layer; ++curr_layer) {
            int old_nodes_count = SetCardinality(old_nodes);
            int old_nodes_arr[old_nodes_count];
            assert(SetToArray(old_nodes_arr, old_nodes) == old_nodes_count);

            for (old_i = 0; old_i < old_nodes_count; ++old_i) {
                old_node = old_nodes_arr[old_i];

                for (neigh_i = 0; neigh_i < G->degree[old_node]; ++neigh_i) {
                    neigh_node = G->neighbor[old_node][neigh_i];

                    if (!SetIn(all_nodes, neigh_node)) {
                        SetAdd(new_nodes, neigh_node);
                        SetAdd(all_nodes, neigh_node);
                    }
                }
            }

            SetFree(old_nodes); // doing this is faster than copying new_nodes to old_nodes
            old_nodes = new_nodes;
            new_nodes = SetAlloc(G->n);
        }

        int old_nodes_count = SetCardinality(old_nodes);
        int old_nodes_arr[old_nodes_count];
        assert(SetToArray(old_nodes_arr, old_nodes) == old_nodes_count);
        layer_deg_sum = 0;

        for (old_i = 0; old_i < old_nodes_count; ++old_i) {
            old_node = old_nodes_arr[old_i];
            layer_deg_sum += G->degree[old_node];
        }

        orderArray[curr_i].degree_val = layer_deg_sum;
        SetEmpty(old_nodes);
        SetEmpty(all_nodes);
    }

    SetFree(old_nodes);
    SetFree(new_nodes);

    // sort start through end by degree val
    qsort((void*)(orderArray + start), end - start, sizeof(node_wdegree), nwd_descompFunc);

    // sort next layer for all tied sections between start and end
    int new_start;
    int new_end;

    for (i = start; i < end; ++i) {
        if (i == start || orderArray[i].degree_val != orderArray[i - 1].degree_val) {
            new_start = i;
        }

        if (i == end || orderArray[i].degree_val != orderArray[i + 1].degree_val) {
            new_end = i + 1;

            if (new_end - new_start >= 2) { // make sure it's a tie with more than one number
                if (orderArray[i].degree_val != 0) { // if the degree is 0, it's an edge case because none of the nodes' degree_vals will change, causing an infinite loop
                    enumerateDegreeOrderHelper(G, orderArray, new_start, new_end, layer + 1);
                }
            }
        }
    }
}

void antidupFillNextStepSet(SET **next_step_pointer, SET **deg_set_pointer, GRAPH *G, int *prev_nodes_array, int prev_nodes_count, int *degreeOrder) {
    assert(_useAntidup); // make sure this setting is set
    assert(degreeOrder != NULL); // assert this because degreeOrder is NULL when useAntidup is true
    int i, j, neigh, rolling_max[prev_nodes_count];

    int last_degree_order = degreeOrder[prev_nodes_array[prev_nodes_count-1]];
    rolling_max[prev_nodes_count-1] = last_degree_order;
    for(i=prev_nodes_count-2; i>=0; i--) {
        int curr_node = prev_nodes_array[i];
        int curr_degree_order = degreeOrder[curr_node];
        if (curr_degree_order > rolling_max[i + 1]) {
            rolling_max[i] = curr_degree_order;
        } else {
            rolling_max[i] = rolling_max[i + 1];
        }
    }
    // Then, go forwards in the array and find valid neighbors
    SET *next_step = SetAlloc(G->n);
    SET *deg_set = SetAlloc(G->n);
    *next_step_pointer = next_step;
    *deg_set_pointer = deg_set;
    SET *all_neighbors = SetAlloc(G->n);
    int first_degree_order = degreeOrder[prev_nodes_array[0]];
    for(i=0; i<prev_nodes_count; i++) {
        for(j=0; j<G->degree[prev_nodes_array[i]]; j++) { // loop through all neighbors for the current node
            neigh = G->neighbor[prev_nodes_array[i]][j];
            int neigh_degree_order = degreeOrder[neigh];
            if (!SetIn(all_neighbors, neigh) && !arrayIn(prev_nodes_array, prev_nodes_count, neigh)) { // we're only processing it if it's not an old neighbor and it's not a node in the graphlet we're building
                Boolean neigh_is_valid = true;
                if (i == prev_nodes_count - 1) {
                    if (neigh_degree_order > first_degree_order) {
                        neigh_is_valid = true;
                    } else {
                        neigh_is_valid = false;
                    }
                } else {
                    if (neigh_degree_order > rolling_max[i + 1]) {
                        neigh_is_valid = true;
                    } else {
                        neigh_is_valid = false;
                    }
                }
                if (neigh_is_valid) {
                    SetAdd(next_step, neigh);
                    SetAdd(deg_set, G->degree[neigh]);
                }
                SetAdd(all_neighbors, neigh);
            }
        }
    }
    SetFree(all_neighbors);
}
