#include <sys/file.h>
#include "blant.h"

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
void mapCanonMap(char* BUF, short int *K, int k) {
    int Bk = (1 <<(k*(k-1)/2));
    sprintf(BUF, CANON_DIR "/canon_map%d.bin", k);
    int Kfd = open(BUF, 0*O_RDONLY);
    assert(Kfd > 0);
    short int *Kf = Mmap(K, Bk*sizeof(K[0]), Kfd);
    assert(Kf == K);
}

int canonListPopulate(char *BUF, int *canon_list, int k) {
    sprintf(BUF, CANON_DIR "/canon_list%d.txt", k);
    FILE *fp_ord=fopen(BUF, "r");
    if(!fp_ord) Fatal("cannot find %s/canon_list%d.txt\n", CANON_DIR, k);
    int numCanon, i, connected, numEdges;
    fscanf(fp_ord, "%d",&numCanon);
    for(i=0; i<numCanon; i++) fscanf(fp_ord, "%d %d %d", &canon_list[i], &connected, &numEdges);
    fclose(fp_ord);
    return numCanon;
}
	
int orbitListPopulate(char *BUF, int orbit_list[MAX_CANONICALS][maxK], int k) {
    sprintf(BUF, CANON_DIR "/orbit_map%d.txt", k);
    FILE *fp_ord=fopen(BUF, "r");
    if(!fp_ord) Fatal("cannot find %s/orbit_map%d.txt\n", CANON_DIR, k);
    int numOrbit, i, j;
    fscanf(fp_ord, "%d",&numOrbit);
    for(i=0; i<numOrbit; i++) for(j=0; j<k; j++)fscanf(fp_ord, "%d", &orbit_list[i][j]);
    fclose(fp_ord);
    return numOrbit;
}
