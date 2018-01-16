// This is the maximum graphlet size that BLANT supports.  Cannot be bigger than 8.
// Currently only used to determine the amount of static memory to allocate.
#define maxK 8

#define MAX_CANONICALS 12346	// This is the number of canonical graphettes for k=8

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

// Try to use mmap? As long as mmap is supported, even if it doesn't work blant will still work 
// because it'll revert to simply loading the entire binary mapping if the mmap fails.
#define MMAP 1

#define CANON_DIR "canon_maps"
//#define CANON_DIR "/var/preserve/Graphette/data" // if you happen to put it there...
