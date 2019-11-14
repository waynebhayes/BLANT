#include "tinygraph.h"
#include "sets.h"

// This is the maximum graphlet size that BLANT supports.  Cannot be bigger than 8.
// Currently only used to determine the amount of static memory to allocate.
#define maxK 8

#define maxBk (1 << (maxK*(maxK-1)/2)) // maximum number of entries in the canon_map

#define mcmc_d 2 // arbitrary d graphlet size < k for MCMC algorithm. Should always be 2 or k-1

#define MAX_CANONICALS	12346	// This is the number of canonical graphettes for k=8
#define MAX_ORBITS	79264	// This is the number of orbits for k=8

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

#define CANON_DIR "canon_maps"
//#define CANON_DIR "/var/preserve/Graphette/canon_maps" // if you happen to put it there...

#define DEFAULT_BLANT_DIR "."
extern char* _BLANT_DIR;

void BuildGraph(TINY_GRAPH* G, int Gint);
int TinyGraph2Int(TINY_GRAPH *g, int numNodes);
short int* mapCanonMap(char* BUF, short int *K, int k);
SET *canonListPopulate(char *BUF, int *canon_list, int k); // returns a SET containing list of connected ordinals
int orbitListPopulate(char *BUF, int orbit_list[MAX_CANONICALS][maxK],  int orbit_canon_mapping[MAX_ORBITS], int numCanon, int k);
void orcaOrbitMappingPopulate(char *BUF, int orca_orbit_mapping[58], int k);
char** convertToEL(char* file); // from convert.cpp
