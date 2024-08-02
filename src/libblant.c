#include <sys/file.h>
#include <sys/mman.h>
#include <assert.h>
#include "blant.h"
#include "tinygraph.h"

const char* _BLANT_DIR = DEFAULT_BLANT_DIR;

// Given a TINY_GRAPH and k, return the integer ID created from one triangle (upper or lower) of the adjacency matrix.
Gint_type TinyGraph2Int(TINY_GRAPH *g, int k)
{
    int i, j, bitPos=0;
    Gint_type Gint = 0, bit;

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
		bit = (((Gint_type)1) << bitPos);
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
void Int2TinyGraph(TINY_GRAPH* G, Gint_type Gint)
{
    int i, j, bitPos=0, k = G->n;
    Gint_type Gint2 = Gint;  // Gint2 has bits nuked as they're used, so when it's zero we can stop.
    TinyGraphEdgesAllDelete(G);
#if LOWER_TRIANGLE
    for(i=k-1;i>(0-SELF_LOOPS);i--) // note that when SELF_LOOPS, the loop condition i (-1), so i must be signed.
    {
	for(j=i-!SELF_LOOPS;j>=0;j--)
#else	// UPPER_TRIANGLE
    assert(!SELF_LOOPS);
    for(i=k-2;i>=0;i--)
    {
	for(j=k-1;j>i;j--)
#endif
	{
	    if(!Gint2) break;
	    Gint_type bit = ((Gint_type)1 << bitPos);
	    if(Gint & bit)
		TinyGraphConnect(G,i,j);
	    Gint2 &= ~bit;
	    bitPos++;
	    assert(bitPos < 8*sizeof(Gint_type)); // technically they could be equal... change when that happens
	}
	if(!Gint2) break;
    }
}

/*
** Given a pre-allocated filename buffer, a 256MB aligned array K, num nodes k
** Mmap the canon_map binary file to the aligned array.
*/
Gordinal_type* mapCanonMap(char* BUF, Gordinal_type *K, int k) {
#if SELF_LOOPS
    if (k > 7) Fatal("cannot have k>7 when SELF_LOOPS");
    int Bk = (1 <<(k*(k+1)/2));
#else
    int Bk = (1 <<(k*(k-1)/2));
#endif
    sprintf(BUF, "%s/%s/canon_map%d.bin", _BLANT_DIR, CANON_DIR, k);
    int Kfd = open(BUF, 0*O_RDONLY);
    assert(Kfd > 0);
    //short int *Kf = Mmap(K, Bk*sizeof(short int), Kfd); // Using Mmap will cause error due to MAP_FIXED flag
    Gordinal_type *Kf = (Gordinal_type*) mmap(K, sizeof(Gordinal_type)*Bk, PROT_READ, MAP_PRIVATE, Kfd, 0);
    assert(Kf != MAP_FAILED);
    return Kf;
}

SET *canonListPopulate(char *BUF, Gint_type *canon_list, int k, char *canon_num_edges) {
    sprintf(BUF, "%s/%s/canon_list%d.txt", _BLANT_DIR, CANON_DIR, k);
    FILE *fp_ord=fopen(BUF, "r");
    if(!fp_ord) Fatal("cannot find %s\n", BUF);
    Gordinal_type numCanon, i;
    int connected;
    assert(1==fscanf(fp_ord, GORDINAL_FMT "\n",&numCanon));
    SET *connectedCanonicals = SetAlloc(numCanon);
    for(i=0; i<numCanon; i++) {
	char buf[BUFSIZ], *tmp;
	tmp = fgets(buf, sizeof(buf), fp_ord);
	assert(tmp == buf);
	assert(3==sscanf(buf, GINT_FMT "\t%d %hhd", &canon_list[i], &connected, &canon_num_edges[i]));
	if(connected) SetAdd(connectedCanonicals, i);
    }
    fclose(fp_ord);
    return connectedCanonicals;
}

Gint_type orbitListPopulate(char *BUF,
	Gint_type orbit_list[MAX_CANONICALS][MAX_K],
	Gordinal_type orbit_canon_mapping[MAX_ORBITS],
	char orbit_canon_node_mapping[MAX_ORBITS],
	Gordinal_type numCanon, int k) {
    sprintf(BUF, "%s/%s/orbit_map%d.txt", _BLANT_DIR, CANON_DIR, k);
    FILE *fp_ord=fopen(BUF, "r");
    if(!fp_ord) Fatal("cannot find %s\n", BUF);
    Gint_type o, numOrbit;
    assert(1==fscanf(fp_ord, GINT_FMT, &numOrbit));
    for(o=0;o<numOrbit;o++)
	orbit_canon_node_mapping[o] = -1;
    Gordinal_type c; int j;
    for(c=0; c<numCanon; c++) for(j=0; j<k; j++) {
	assert(1==fscanf(fp_ord, GINT_FMT, &orbit_list[c][j]));
	orbit_canon_mapping[orbit_list[c][j]] = c;
	if(orbit_canon_node_mapping[orbit_list[c][j]] < 0)
	    orbit_canon_node_mapping[orbit_list[c][j]] = j;
    }
    fclose(fp_ord);
    return numOrbit;
}

void orcaOrbitMappingPopulate(char *BUF, int orca_orbit_mapping[58], int k) {
    assert(k<=5);
    sprintf(BUF, "%s/%s/orca_orbit_mapping%d.txt", _BLANT_DIR, "orca_jesse_blant_table", k);
    FILE *fp_ord=fopen(BUF, "r");
    if(!fp_ord) Fatal("cannot find %s\n", BUF);
    int numOrbit, i;
    assert(1==fscanf(fp_ord, "%d", &numOrbit));
    for (i=0; i<numOrbit; i++) { assert(1==fscanf(fp_ord, "%d", &orca_orbit_mapping[i])); }
}
