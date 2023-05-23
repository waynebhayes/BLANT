#include "blant.h"
#include "blant-output.h"
#include "blant-predict.h"
#include "blant-utils.h"
#include "blant-sampling.h"
#include "sorts.h"

#define SORT_INDEX_MODE 0 // Note this destroys the columns-are-identical property, don't use by default.

char *PrintNode(char c, int v) {
    static char buf[BUFSIZ];
    char *s=buf;
    if(c) *s++ = c;
    if(_supportNodeNames)
	sprintf(s,"%s", _nodeNames[v]);
    else
	sprintf(s, "%d", v);
    return buf;
}

char *PrintNodePairSorted(int u, char c, int v) {
    static char buf[BUFSIZ];
    if(_supportNodeNames) {
	char *s1=_nodeNames[u], *s2=_nodeNames[v];
	if(strcmp(s1,s2)<0) { char *tmp=s1;s1=s2;s2=tmp; }
	sprintf(buf, "%s%c%s", s1,c,s2);
    }
    else
	sprintf(buf, "%d%c%d", MAX(u,v),c,MIN(u,v));
    return buf;
}

static int IntCmp(const void *a, const void *b)
{
    const int *i = (const int*)a, *j = (const int*)b;
    return (*i)-(*j);
}

void VarraySort(unsigned *Varray, int k)
{
#if USE_INSERTION_SORT
    //InsertionSortInt(Varray,k);
    InsertionSort((void*)Varray,k,sizeof(Varray[0]),IntCmp);
#else
    qsort((void*)Varray, k, sizeof(Varray[0]), IntCmp);
#endif
}

char *PrintCanonical(int GintOrdinal)
{
    static char buf[BUFSIZ];
    int j, GintNumBits = _k*(_k-1)/2;
    char GintBinary[GintNumBits+1]; // Only used in -db output mode for indexing
    switch (_displayMode) {
    case undefined:
    case ordinal:
	sprintf(buf, "%d", GintOrdinal);
	break;
    case decimal: // Prints the decimal integer form of the canonical
#if TINY_SET_SIZE <= 32
	sprintf(buf, "%u", _canonList[GintOrdinal]);
#elif TINY_SET_SIZE >= 64
	sprintf(buf, "%llu", _canonList[GintOrdinal]);
#else
#error "unknown TINY_SET_SIZE"
#endif
	break;
    case binary: // Prints the bit representation of the canonical
	for (j=0;j<GintNumBits;j++)
	    {GintBinary[GintNumBits-j-1]=((_canonList[GintOrdinal] >> j) & 1 ? '1' : '0');}
	GintBinary[GintNumBits] = '\0';
	strcpy(buf, GintBinary);
	break;
    case orca: // Prints the ORCA ID of the canonical. Jesse uses same number.
    case jesse:
	sprintf(buf, "%d", _outputMapping[GintOrdinal]);
	break;
    }
    return buf;
}

// Below is code to help reduce (mostly eliminite if we're lucky) MCMC's duplicate output, which is copious
// Empirically I've found that neither of these need to be very big since most repetition happens with high locality
//#define MCMC_CIRC_BUF 999983 // prime, not sure it needs to be but why not... 1M * 4b = 4MB RAM.
#define MCMC_CIRC_BUF 134217689 // 2^27-39, prime according to https://www.dcode.fr/closest-prime-number (512MB)

// But apparently the above are too big, not sure why, but they cause seg faults.
//#define MCMC_MAX_HASH 1610612741U // about 1.6 billion, taken from https://planetmath.org/goodhashtableprimes
//#define MCMC_MAX_HASH  402653189U // halfway from 2^28 to 2^29 (works)
//#define MCMC_MAX_HASH  805306457U // halfway from 2^29 to 2^30 (works)
//#define MCMC_MAX_HASH 1610612741U // halfway from 2^30 to 2^31 (works)

// Below is simply the largest number in the list at:
// https://github.com/leventov/Koloboke/lib/impl/src/main/java/net/openhft/koloboke/collect/impl/hash/DHashCapacities.java
#define MCMC_MAX_HASH 2136745621U // about 1% smaller than 2^31.
// The following is > 2^32, and would requires a SET implementation allowing members with value > 32 bits.
//#define MCMC_MAX_HASH 8589934591UL // 2^33-9, according to https://www.dcode.fr/closest-prime-number; about 1GB

// NOTE WE DO NOT CHECK EDGES. So if you call it with the same node set but as a motif, it'll (incorrectly) return TRUE
Boolean NodeSetSeenRecently(GRAPH *G, unsigned Varray[], int k) {
    static unsigned circBuf[MCMC_CIRC_BUF], bufPos;
    static BITVEC *seen;
    static unsigned Vcopy[MAX_K];
    unsigned i;
    if(!seen) seen=BitvecAlloc(MCMC_MAX_HASH);
    memcpy(Vcopy, Varray, k*sizeof(*Varray));
    VarraySort(Vcopy, k);

    if (_outputMode == indexGraphletsRNO) {
        // move the first node in Varray (the root node) to the start of Vcopy
        // this is because we now consider identical sets of nodes different if they were created in a different order (specifically, if the root node was different)
        unsigned base_node = Varray[0];

        if (Vcopy[0] != base_node) {
            unsigned stored_node = Vcopy[0];

            for (i = 1; i < k; i++) {
                unsigned tmp = stored_node;
                stored_node = Vcopy[i];
                Vcopy[i] = tmp;

                if (stored_node == base_node) {
                    Vcopy[0] = stored_node;
                    break;
                }
            }
        }
    }

    unsigned hash=Vcopy[0];
    for(i=1;i<k;i++) hash = hash*G->n + Vcopy[i]; // Yes this will likely overflow. Shouldn't matter.
    hash = hash % MCMC_MAX_HASH;
    if(BitvecInSafe(seen, hash)) return true; // of course false positives are possible but we hope they are rare.
    //for(i=0;i<k;i++) printf("%d ", Vcopy[i]); printf("\thash %d\n",hash); // checking for rareness.
    BitvecDelete(seen, circBuf[bufPos]); // this set hasn't been seen in at least the last MCMC_CIRC_BUF samples
    circBuf[bufPos] = hash;
    BitvecAdd(seen, hash);
    if(++bufPos >= MCMC_CIRC_BUF) bufPos=0;
    return false;
}

char *PrintIndexEntry(Gint_type Gint, int GintOrdinal, unsigned Varray[], TINY_GRAPH *g, int k)
{
    int j;
    unsigned char perm[MAX_K];
#if SORT_INDEX_MODE
    VarraySort(Varray, k);
    for(j=0;j<k;j++) perm[j]=j;
#else
    memset(perm, 0, k);
    ExtractPerm(perm, Gint);
    assert(PERMS_CAN2NON);
#endif
    static char buf[2][BUFSIZ];
    int which=0; // which should ALWAYS point to the one that HAS the data
    strcpy(buf[which], PrintCanonical(GintOrdinal));

    // IMPORTANT NOTE: this code prints the perm, not the orbit (ambiguous graphlets have repeating orbits but don't have repeating perms). If all graphlets are unambiguous, doing this is fine (since perm will be a bijection with orbit). However, if you want to extract ambiguous graphlets, you'll have to change the code here (and code in a lot of other places)
    if (_outputMode == indexGraphletsRNO) {
        sprintf(buf[1-which], "%s+o%d", buf[which], perm[0]);
	which=1-which;
    }

    for(j=0;j<k;j++) {
	sprintf(buf[1-which], "%s%s", buf[which], PrintNode(' ', Varray[(int)perm[j]]));
	which=1-which;
    }
    return buf[which];
}

char *PrintIndexOrbitsEntry(Gint_type Gint, int GintOrdinal, unsigned Varray[], TINY_GRAPH *g, int k) {
    assert(TinyGraphDFSConnected(g,0));
    int j;
    static SET* printed;
    if(!printed) printed = SetAlloc(k);
    SetEmpty(printed);
    unsigned char perm[MAX_K+1];
#if SORT_INDEX_MODE
    VarraySort(Varray, k);
    for(j=0;j<k;j++) perm[j]=j;
#else
    memset(perm, 0, k);
    ExtractPerm(perm, Gint);
    assert(PERMS_CAN2NON); // Apology("Um, don't we need to check PERMS_CAN2NON? See outputODV for correct example");
#endif
    static char buf[2][BUFSIZ];
    int which=0;
    strcpy(buf[which], PrintCanonical(GintOrdinal));
    for(j=0;j<k;j++) if(!SetIn(printed,j))
    {
	which=1-which; sprintf(buf[which], "%s%s", buf[1-which], PrintNode(' ', Varray[(int)perm[j]]));
	SetAdd(printed, j);
	int j1;
	for(j1=j+1;j1<k;j1++) if(_orbitList[GintOrdinal][j1] == _orbitList[GintOrdinal][j])
	{
	    assert(!SetIn(printed, j1));
	    which=1-which; sprintf(buf[which], "%s%s", buf[1-which], PrintNode(':', Varray[(int)perm[j1]]));
	    SetAdd(printed, j1);
	}
    }
    return buf[which];
}

void ProcessNodeOrbitNeighbors(GRAPH *G, Gint_type Gint, int GintOrdinal, unsigned Varray[], TINY_GRAPH *g, int k) {
    assert(TinyGraphDFSConnected(g,0));
    int c,d; // canonical nodes
    unsigned char perm[MAX_K+1];
#if SORT_INDEX_MODE
    VarraySort(Varray, k);
    for(c=0;c<k;c++) perm[c]=c;
#else
    assert(PERMS_CAN2NON); // Apology("Um, don't we need to check PERMS_CAN2NON? See outputODV for correct example");
    memset(perm, 0, k);
    ExtractPerm(perm, Gint);
#endif
    for(c=0;c<k;c++)
    {
	// FIXME: This could be made more memory efficient by indexing ONLY the CONNECTED canonicals, but...
	int u=Varray[(int)perm[c]], u_orbit=_orbitList[GintOrdinal][c];
	++ODV(u, u_orbit);
	for(d=c+1;d<k;d++) if(TinyGraphAreConnected(g,c,d)) {
	    int v=Varray[(int)perm[d]]; // v_orbit=_orbitList[GintOrdinal][d];
	    if(!_communityNeighbors[u]) _communityNeighbors[u] = (SET**) Calloc(_numOrbits, sizeof(SET*));
	    if(!_communityNeighbors[u][u_orbit]) _communityNeighbors[u][u_orbit] = SetAlloc(G->n);
	    SetAdd(_communityNeighbors[u][u_orbit], v);
	}
    }
}

#if 0	// FIXME: related to the above.... this was how I attempted to have ONLY the CONNECTED canonicals,
	// but it seemed buggy and I can't figure out what's wrong. - WH
	int u=Varray[(int)perm[c]], u_orbit=_orbitList[GintOrdinal][c], u_connectedOrb=_orbit2connectedIndex[u_orbit];
	assert(0 <= u_connectedOrb && u_connectedOrb < _numConnectedOrbits);
	if(_communityMode=='o') ++ODV(u, u_orbit); else ++GDV(u, GintOrdinal);
	int item = (_communityMode=='o' ? u_connectedOrb : _canon2connectedIndex[GintOrdinal]);
#endif

void ProcessNodeGraphletNeighbors(GRAPH *G, Gint_type Gint, int GintOrdinal, unsigned Varray[], TINY_GRAPH *g, int k) {
    assert(TinyGraphDFSConnected(g,0));
    int c,d; // canonical nodes
    unsigned char perm[MAX_K+1];
#if SORT_INDEX_MODE
    VarraySort(Varray, k);
    for(c=0;c<k;c++) perm[c]=c;
#else
    assert(PERMS_CAN2NON); // Apology("Um, don't we need to check PERMS_CAN2NON? See outputODV for correct example");
    memset(perm, 0, k);
    ExtractPerm(perm, Gint);
#endif
    for(c=0;c<k;c++)
    {
	// FIXME: This could be made more memory efficient by indexing ONLY the CONNECTED canonicals, but...
	int u=Varray[(int)perm[c]];
	++GDV(u, GintOrdinal);
	for(d=c+1;d<k;d++) if(TinyGraphAreConnected(g,c,d)) {
	    int v=Varray[(int)perm[d]];
	    if(!_communityNeighbors[u]) _communityNeighbors[u] = (SET**) Calloc(_numCanon, sizeof(SET*));
	    if(!_communityNeighbors[u][GintOrdinal]) _communityNeighbors[u][GintOrdinal] = SetAlloc(G->n);
	    SetAdd(_communityNeighbors[u][GintOrdinal], v);
	}
    }
}

Boolean ProcessGraphlet(GRAPH *G, SET *V, unsigned Varray[], const int k, TINY_GRAPH *g)
{
    Boolean processed = true;
    TinyGraphInducedFromGraph(g, G, Varray);
    Gint_type Gint = TinyGraph2Int(g,k), GintOrdinal=L_K(Gint), j;

#if PARANOID_ASSERTS
    assert(0 <= GintOrdinal && GintOrdinal < _numCanon);
#endif
    switch(_outputMode)
    {
	unsigned char perm[MAX_K];
    case graphletFrequency:
	++_graphletCount[GintOrdinal];
	break;
    case indexGraphlets: case indexGraphletsRNO:
	if(NodeSetSeenRecently(G, Varray,k) ||
	    (_sampleMethod == SAMPLE_INDEX && !SetIn(_windowRep_allowed_ambig_set, GintOrdinal)) ||
	    _canonNumEdges[GintOrdinal] < _min_edge_count) processed=false;
	else puts(PrintIndexEntry(Gint, GintOrdinal, Varray, g, k));
	break;
    case predict:
	if(NodeSetSeenRecently(G,Varray,k)) processed=false;
	else Predict_AccumulateMotifs(G,Varray,g,Gint,GintOrdinal);
	break;
    case indexOrbits:
	if(NodeSetSeenRecently(G,Varray,k) ||
	    (_sampleMethod == SAMPLE_INDEX && !SetIn(_windowRep_allowed_ambig_set, GintOrdinal))) processed=false;
	else puts(PrintIndexOrbitsEntry(Gint, GintOrdinal, Varray, g, k));
	break;
    case communityDetection:
	if(_canonNumEdges[GintOrdinal] < _min_edge_count) processed=false;
	else if(_communityMode == 'o') ProcessNodeOrbitNeighbors(G, Gint, GintOrdinal, Varray, g, k);
	else if(_communityMode == 'g') ProcessNodeGraphletNeighbors(G, Gint, GintOrdinal, Varray, g, k);
	else Fatal("unkwown _communityMode %c", _communityMode);
	break;
    case outputGDV:
	for(j=0;j<k;j++) ++GDV(Varray[j], GintOrdinal);
	break;
    case outputODV:
	memset(perm, 0, _k);
	ExtractPerm(perm, Gint);
#if PERMS_CAN2NON
	for(j=0;j<k;j++) ++ODV(Varray[(int)perm[j]], _orbitList[GintOrdinal][          j ]);
#else
	for(j=0;j<k;j++) ++ODV(Varray[          j ], _orbitList[GintOrdinal][(int)perm[j]]);
#endif
	break;

    default: Abort("ProcessGraphlet: unknown or un-implemented outputMode %d", _outputMode);
	break;
    }
    return processed;
}
