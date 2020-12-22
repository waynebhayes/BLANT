#include "blant.h"
#include "blant-output.h"
#include "blant-kovacs.h"
#include "blant-utils.h"
#include <math.h>

int _kovacsOrdinal = -1;
int _kovacsOrbit1 = -1, _kovacsOrbit2=-1; // the columns (ie orbits) you want
float **_KovacsScore, **_KovacsNorm; // will be allocated lower triangle of n by n matrix of node-pair counts
int _KovacsMotifCount[MAX_CANONICALS][maxK][maxK]; //pre-computed matrices of 64-bit ints take only about 6MB
float _KovacsMotifNorm[MAX_CANONICALS][maxK][maxK]; 
unsigned _kovacsOrbitPairSeen[MAX_CANONICALS][maxK][maxK];
unsigned _kovacsOrbitPairEdge[MAX_CANONICALS][maxK][maxK];

// Recursively count the specified orbits in the specified canonical motifs of g.
// The topOrdinal is the canonical graphlet who's count we're computing as we descend the recursion to
// enumerate all the counts of the One True Motif's node pairs the user has given us (eg, the L3).
// The perm array tells us how to map nodes in g all the way back to nodes in topCanon.
void PreComputeKovacs(TINY_GRAPH *g, int topOrdinal, TINY_GRAPH *gTop, char perm[maxK])
{
    static Boolean initDone;
    static SET *seen; // size 2^B(k), *not* canonical but a specific set of nodes and edges in the *top* graphlet
    if(!initDone) {
	assert(_Bk>0);
	seen = SetAlloc(_Bk);
	assert(_k>= 3 && _k <= 8);
	gTop = TinyGraphAlloc(_k);
	initDone = true;
    }
    static int depth;
#if PARANOID_ASSERTS
    assert(g->n == _k);
    assert(TinyGraphDFSConnected(g, 0));
#endif
    int i,j, Gint = TinyGraph2Int(g,_k), GintOrdinal=_K[Gint];
    char gPerm[maxK], permComposed[maxK];
    memset(gPerm, 0, _k);
    ExtractPerm(gPerm, Gint);
    if(depth == 0) {
#if PARANOID_ASSERTS
	assert(GintOrdinal == topOrdinal && Gint == _canonList[topOrdinal]); // we should only be called on canonicals
	for(i=0;i<_k;i++) assert(gPerm[i]==i && perm[i] == i);
#endif
    }
    // OK, think carefully here: gPerm maps the nodes in this g back to the graphlet one level up,
    // which contains the extra edge we've had deleted. Our gPerm maps the nodes in *this* canonical
    // back up to nodes one level up; that one in turn maps up another level. We want to increment the
    // values of the node pair all the way up in topCanon. So, we need to compose our gPerm with the
    // one we were passed:
    assert(PERMS_CAN2NON);
    for(i=0;i<_k;i++) permComposed[i] = perm[(int)gPerm[i]];
    TinyGraphEdgesAllDelete(gTop); // build this motif at the top level locations
    for(i=0; i<_k-1; i++)for(j=i+1;j<_k;j++)
	if(TinyGraphAreConnected(g,i,j)) TinyGraphConnect(gTop,permComposed[i],permComposed[j]);
    int gIntTop = TinyGraph2Int(gTop,_k);
    assert(gIntTop < _Bk);

    // If we have seen this exact non-canonical motif in the top-level set of nodes, we don't need to see them again.
    if(SetIn(seen, gIntTop)) return;
    SetAdd(seen,gIntTop);
   
    // Check to see if the whole thing is the cononical of interest before we start removing edges.
    if(GintOrdinal == _kovacsOrdinal) {
	assert(PERMS_CAN2NON);
        for(i=0;i<_k-1;i++) for(j=i+1;j<_k;j++) {
	    int uTop = permComposed[i], vTop = permComposed[j];
            if(_orbitList[GintOrdinal][i] == _kovacsOrbit1 && _orbitList[GintOrdinal][j] == _kovacsOrbit2){
		    if(!TinyGraphAreConnected(g,i,j)) // increment the count if there's no edge already here in the *MOTIF*
			++_KovacsMotifCount[topOrdinal][uTop][vTop];
	    } else {
		if(TinyGraphAreConnected(g,i,j))
		    _KovacsMotifNorm[topOrdinal][uTop][vTop] += sqrt(gTop->degree[uTop]*gTop->degree[vTop]);
	    }
	}
	// Stop the recursion since removing more edges can't possibly get the _kovacsOrdinal again.
	return;
    }
    // Now go about deleting edges recursively.
    for(i=0; i<_k-1; i++)for(j=i+1;j<_k;j++)
    {
	if(TinyGraphAreConnected(g,i,j)) // if it's an edge, then delete it, recurse, and re-attach upon return.
	{
	    TinyGraphDisconnect(g,i,j);
	    if(TinyGraphDFSConnected(g,0)) {
		++depth;
		PreComputeKovacs(g,topOrdinal,gTop,permComposed);
		--depth;
	    }
	    TinyGraphConnect(g,i,j);
	}
    }
}

// Recursively count the specified orbits in the specified canonical motifs of g.
void ProcessKovacsNorm(TINY_GRAPH *g, GRAPH *G, unsigned Varray[])
{
    static int depth;
    static Boolean initDone;
    static SET *seen; // size 2^B(k), *not* canonical but a specific set of nodes and edges in the *top* graphlet
    if(!initDone) {
        assert(_Bk>0);
        seen = SetAlloc(_Bk);
        assert(_k>= 3 && _k <= 8);
        initDone = true;
    }

    if(depth==0) SetReset(seen);

#if PARANOID_ASSERTS
    assert(g->n == _k);
    assert(TinyGraphDFSConnected(g, 0));
#endif

    int i,j, Gint = TinyGraph2Int(g,_k), GintOrdinal=_K[Gint];
    if(SetIn(seen,Gint)) return;
    SetAdd(seen,Gint);

    // Check to see if the whole thing is the cononical of interest before we start removing edges.
    if(GintOrdinal == _kovacsOrdinal) {
        char perm[maxK];
        memset(perm, 0, _k);
        ExtractPerm(perm, Gint);
        for(i=0;i<_k-1;i++) for(j=i+1;j<_k;j++) // if(TinyGraphAreConnected(g,i,j))
        {
            int u = Varray[(int)perm[i]], v= Varray[(int)perm[j]];
            if(_orbitList[GintOrdinal][i] == _kovacsOrbit1 && _orbitList[GintOrdinal][j] == _kovacsOrbit2)
            {
            int l;
            double norm = 1;
            for(l=0;l<_k;l++) if(l!=i&&l!=j) norm *= G->degree[Varray[(int)perm[l]]];
            norm = pow(norm, 1./(_k-2)); // other values, such as _k-1, sometimes work better.
            _KovacsScore[MAX(u,v)][MIN(u,v)] += 1./norm;
            }
        }
    #if 0
        printf("\tM %d ",depth);
        PrintCanonical(GintOrdinal);
        for(l=0;l<_k;l++) PrintNode(' ',Varray[(int)perm[l]]);
        putchar('\n');
    #endif
        // Stop the recursion since removing more edges can't possibly get the canonical again.
        return;
    }

    // Now go about deleting edges recursively.
    for(i=0; i<_k-1; i++)for(j=i+1;j<_k;j++)
    {
        if(TinyGraphAreConnected(g,i,j)) // if it's an edge, delete it.
        {
            TinyGraphDisconnect(g,i,j);
            if(TinyGraphDFSConnected(g,0)) {
            ++depth;
            ProcessKovacsNorm(g,G,Varray);
            --depth;
            }
            TinyGraphConnect(g,i,j);
        }
    }
}

void ProcessKovacsPreComputed(TINY_GRAPH *g, unsigned Varray[])
{
    static int depth;
#if PARANOID_ASSERTS
    assert(PERMS_CAN2NON);
    assert(g->n == _k);
    assert(TinyGraphDFSConnected(g, 0));
#endif

    int i,j, Gint = TinyGraph2Int(g,_k), GintOrdinal=_K[Gint];
    char perm[maxK];
    memset(perm, 0, _k);
    ExtractPerm(perm, Gint);
    for(i=0;i<_k-1;i++)
	for(j=i+1;j<_k;j++)
	{
	    int u = Varray[(int)perm[i]], v = Varray[(int)perm[j]];
	    assert(u!=v);
	    // Order of nodes doesn't matter, and _KovacsScore is only the lower triangle of node pairs.
	    float norm = 1;
	    if(_KovacsMotifNorm[GintOrdinal][i][j]) norm = _KovacsMotifNorm[GintOrdinal][i][j];
	    _KovacsScore[MAX(u,v)][MIN(u,v)] += _KovacsMotifCount[GintOrdinal][i][j]/norm;
	}
}

void KovacsProcessIndexEntry(GRAPH *G, int Gint, int GintOrdinal, unsigned Varray[], char perm[], TINY_GRAPH *g, int k)
{
    int i, j;
    assert(0<=GintOrdinal && GintOrdinal<_numCanon);

    int baseOrbit = _orbitList[GintOrdinal][0]; // the 0th orbit of this ordinal.
    for(i=0;i<k-1;i++) for(j=i+1;j<k;j++) // loop through all pairs of nodes in the canonical
    {
        // We are trying to determine the frequency that a pair of nodes in G have an edge based on
        // their being located at a pair of orbits in a motif. However, it only makes sense to do
        // this if the edge in the motif does *not* exist; otherwise the edge *always* exists in G.
        if(!TinyGraphAreConnected(g,perm[i],perm[j])) {
            int u=Varray[(int)perm[i]];
            int v=Varray[(int)perm[j]];
            int edge = GraphAreConnected(G,u,v);

            int orbit0=_orbitList[GintOrdinal][i];
            int orbit1=_orbitList[GintOrdinal][j];
            assert(SetIn(_connectedCanonicals,_orbitCanonMapping[orbit0]));
            assert(SetIn(_connectedCanonicals,_orbitCanonMapping[orbit1]));
    #define OUTPUT_MOTIF_ORBIT_PAIR 1
    #if OUTPUT_MOTIF_ORBIT_PAIR
            // preserve order so we can associate the node with the orbit.
            PrintNode(0, u); PrintNode(' ', v);
	        printf(" %d\t%d %d\n", edge, orbit0, orbit1);
    #else
            int oMin=MIN(orbit0,orbit1), oMax=MAX(orbit0,orbit1);
            // when recording the orbit pair without the nodes, order doesn't matter
            ++_kovacsOrbitPairSeen[GintOrdinal][oMin-baseOrbit][oMax-baseOrbit];
            if(edge) ++_kovacsOrbitPairEdge[GintOrdinal][oMin][oMax];
    #endif
        }
    }
}

// Recursively gather stats on orbit pairs of motifs "under" this graphlet
void ProcessKovacsAllOrbits(TINY_GRAPH *g, GRAPH *G, unsigned Varray[])
{
    static int depth;
    static Boolean initDone;
    static SET *seen; // size 2^B(k), *not* canonical but a specific set of nodes and edges in the *top* graphlet
    if(!initDone) {
		assert(_Bk>0);
		seen = SetAlloc(_Bk);
		assert(_k>= 3 && _k <= 8);
		initDone = true;
    }

    if(depth==0) {
	SetReset(seen);
    }

#if PARANOID_ASSERTS
    assert(g->n == _k);
    assert(TinyGraphDFSConnected(g, 0));
#endif

    int Gint = TinyGraph2Int(g,_k);
    if(SetIn(seen,Gint)) return;
    SetAdd(seen,Gint);
    int i,j, GintOrdinal;
    GintOrdinal=_K[Gint];

    char perm[maxK];
    memset(perm, 0, _k);
    ExtractPerm(perm, Gint);
    KovacsProcessIndexEntry(G, Gint, GintOrdinal, Varray, perm, g, _k);

    // Now go about deleting edges recursively.
    for(i=0; i<_k-1; i++)for(j=i+1;j<_k;j++)
    {
	if(TinyGraphAreConnected(g,i,j)) // if it's an edge, delete it.
	{
	    TinyGraphDisconnect(g,i,j);
	    if(TinyGraphDFSConnected(g,0)) {
		++depth;
		ProcessKovacsAllOrbits(g,G,Varray);
		--depth;
	    }
	    TinyGraphConnect(g,i,j);
	}
    }
}
