#include <sys/file.h>
#include <sys/mman.h>
#include "blant.h"

char* _BLANT_DIR = DEFAULT_BLANT_DIR;

// Given a TINY_GRAPH and k, return the integer ID created from one triangle (upper or lower) of the adjacency matrix.
int TinyGraph2Int(TINY_GRAPH *g, int k)
{
    int i, j, bitPos=0, Gint = 0, bit;
    
#if LOWER_TRIANGLE	// Prefer lower triangle to be compatible with Ine Melckenbeeck's Jesse code.
    for(i=k-1;i>0;i--)
    {
        for(j=i-1;j>=0;j--)
#else   // UPPER_TRIANGLE // this is what we used in the original faye code and paper with Adib Hasan and Po-Chien Chung.
    for(i=k-2;i>=0;i--)
    {
        for(j=k-1;j>i;j--)
#endif
        {
	    if(TinyGraphAreConnected(g,i,j))
	    {
		bit = (1 << bitPos);
		Gint |= bit;
	    }
            bitPos++;
        }
    }
    return Gint;
}

/*
** Given an integer, build the graph into the TINY_GRAPH *G, which has already been allocated.
** Handles either upper or lower triangle representation depending upon compile-time option below.
*/
void BuildGraph(TINY_GRAPH* G, int Gint)
{
    int i, j, bitPos=0, k = G->n;
    int Gint2 = Gint;  // Gint2 has bits nuked as they're used, so when it's zero we can stop.
    TinyGraphEdgesAllDelete(G);
#if LOWER_TRIANGLE
    for(i=k-1;i>0;i--)
    {
	for(j=i-1;j>=0;j--)
#else	// UPPER_TRIANGLE
    for(i=k-2;i>=0;i--)
    {
	for(j=k-1;j>i;j--)
#endif
	{
	    if(!Gint2) break;
	    int bit = (1 << bitPos);
	    if(Gint & bit)
		TinyGraphConnect(G,i,j);
	    Gint2 &= ~bit;
	    bitPos++;
	}
	if(!Gint2) break;
    }
}

/*
** Given a pre-allocated filename buffer, a 256MB aligned array K, num nodes k
** Mmap the canon_map binary file to the aligned array. 
*/
short int* mapCanonMap(char* BUF, short int *K, int k) {
    int Bk = (1 <<(k*(k-1)/2));
    sprintf(BUF, "%s/%s/canon_map%d.bin", _BLANT_DIR, CANON_DIR, k);
    int Kfd = open(BUF, 0*O_RDONLY);
    assert(Kfd > 0);
    //short int *Kf = Mmap(K, Bk*sizeof(short int), Kfd); // Using Mmap will cause error due to MAP_FIXED flag
    short int *Kf = (short int*) mmap(K, sizeof(short int)*Bk, PROT_READ, MAP_PRIVATE, Kfd, 0);
    assert(Kf != MAP_FAILED);
    return Kf;
}

int canonListPopulate(char *BUF, int *canon_list, SET *connectedCanonicals, int k) {
    sprintf(BUF, "%s/%s/canon_list%d.txt", _BLANT_DIR, CANON_DIR, k);
    FILE *fp_ord=fopen(BUF, "r");
    if(!fp_ord) Fatal("cannot find %s\n", BUF);
    int numCanon, i, connected, numEdges;
    fscanf(fp_ord, "%d",&numCanon);
    for(i=0; i<numCanon; i++) {
	fscanf(fp_ord, "%d %d %d", &canon_list[i], &connected, &numEdges);
	if(connectedCanonicals && connected) SetAdd(connectedCanonicals, i);
    }
    fclose(fp_ord);
    return numCanon;
}
	
int orbitListPopulate(char *BUF, int orbit_list[MAX_CANONICALS][maxK], int orbit_canon_mapping[MAX_ORBITS], int numCanon, int k) {
    sprintf(BUF, "%s/%s/orbit_map%d.txt", _BLANT_DIR, CANON_DIR, k);
    FILE *fp_ord=fopen(BUF, "r");
    if(!fp_ord) Fatal("cannot find %s\n", BUF);
    int numOrbit, i, j;
    fscanf(fp_ord, "%d",&numOrbit);
    for(i=0; i<numCanon; i++) for(j=0; j<k; j++) { fscanf(fp_ord, "%d", &orbit_list[i][j]); orbit_canon_mapping[orbit_list[i][j]] = i; }
    fclose(fp_ord);
    return numOrbit;
}

void orcaOrbitMappingPopulate(char *BUF, int orca_orbit_mapping[58], int k) {
    sprintf(BUF, "%s/%s/orca_orbit_mapping%d.txt", _BLANT_DIR, "orca_jesse_blant_table", k);
    FILE *fp_ord=fopen(BUF, "r");
    if(!fp_ord) Fatal("cannot find %s\n", BUF);
    int numOrbit, i;
    fscanf(fp_ord, "%d",&numOrbit);
    for (i=0; i<numOrbit; i++) { fscanf(fp_ord, "%d", &orca_orbit_mapping[i]); }
}
