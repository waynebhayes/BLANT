#include "blant.h"
#include "blant-output.h"
#if PREDICT
#include "EdgePredict/blant-predict.h"
#endif
#include "blant-utils.h"
#include "blant-sampling.h"

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
	if(strcmp(s1,s2)>0) { char *tmp=s1;s1=s2;s2=tmp; }
	sprintf(buf, "%s%c%s", s1,c,s2);
    }
    else
	sprintf(buf, "%d%c%d", MIN(u,v),c,MAX(u,v));
    return buf;
}

static int IntCmp(const void *a, const void *b)
{
    int *i = (int*)a, *j = (int*)b;
    return (*i)-(*j);
}

void VarraySort(int *Varray, int k)
{
#if USE_INSERTION_SORT
    InsertionSortInt(Varray,k);
    //InsertionSort((void*)Varray,k,sizeof(Varray[0]),IntCmp);
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
	sprintf(buf, "%d", _canonList[GintOrdinal]);
	break;
    case binary: // Prints the bit representation of the canonical
	for (j=0;j<GintNumBits;j++)
	    {GintBinary[GintNumBits-j-1]=(((unsigned)_canonList[GintOrdinal] >> j) & 1 ? '1' : '0');}
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
#define MCMC_PREDICT_CIRC_BUF 999983 // prime, not sure it needs to be but why not... 1M * 4b = 4MB RAM.
#define MCMC_PREDICT_MAX_HASH 2147483647 // this value should be prime; at 2^31/8 it's 256MB of RAM.
// NOTE WE DO NOT CHECK EDGES. So if you call it with the same node set but as a motif, it'll (incorrectly) return TRUE
Boolean NodeSetSeenRecently(GRAPH *G, unsigned Varray[], int k) {
    static unsigned circBuf[MCMC_PREDICT_CIRC_BUF], bufPos;
    static SET *seen;
    static unsigned Vcopy[MAX_K];
    unsigned hash=0, i;
    if(!seen) seen=SetAlloc(MCMC_PREDICT_MAX_HASH);
    for(i=0;i<k;i++) Vcopy[i]=Varray[i];
    VarraySort(Vcopy, k);
    for(i=0;i<k;i++) hash = hash*G->n + Vcopy[i]; // Yes this will overflow. Shouldn't matter.
    hash = hash % MCMC_PREDICT_MAX_HASH;
    if(SetIn(seen, hash)) return true; // of course false positives are possible but we hope they are rare.
    //for(i=0;i<k;i++) printf("%d ", Vcopy[i]); printf("\thash %d\n",hash); // checking for rareness.
    SetDelete(seen, circBuf[bufPos]);
    SetAdd(seen, hash);
    circBuf[bufPos++] = hash;
    if(bufPos >= MCMC_PREDICT_CIRC_BUF) bufPos=0;
    return false;
}

char *PrintIndexEntry(int Gint, int GintOrdinal, unsigned Varray[], TINY_GRAPH *g, int k)
{
    int j;
    char perm[MAX_K];
    memset(perm, 0, k);
    ExtractPerm(perm, Gint);
    assert(PERMS_CAN2NON);
    static char buf[2][BUFSIZ];
    int which=0;
    strcpy(buf[which], PrintCanonical(GintOrdinal));
    for(j=0;j<k;j++) {
	which=1-which; sprintf(buf[which], "%s%s", buf[1-which], PrintNode(' ', Varray[(int)perm[j]]));
    }
    return buf[which];
}

char *PrintIndexOrbitsEntry(int Gint, int GintOrdinal, unsigned Varray[], TINY_GRAPH *g, int k) {
    assert(TinyGraphDFSConnected(g,0));
    int j;
    static SET* printed;
    if(!printed) printed = SetAlloc(k);
    SetEmpty(printed);
    char perm[MAX_K+1];
    memset(perm, 0, k);
    ExtractPerm(perm, Gint);
    static char buf[2][BUFSIZ];
    int which=0;
    strcpy(buf[which], PrintCanonical(GintOrdinal));
    assert(PERMS_CAN2NON); // Apology("Um, don't we need to check PERMS_CAN2NON? See outputODV for correct example");
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

// Recursively print all the motifs under this graphlet. Note that this one actually outputs directly, it does
// ** not ** print into a buffer and return a char*.
void PrintAllMotifs(TINY_GRAPH *g, int Gint, int GintOrdinal, GRAPH *G, unsigned Varray[])
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

    if(SetIn(seen,Gint)) return;
    SetAdd(seen,Gint);

    if(_outputMode == indexMotifOrbits)
	puts(PrintIndexOrbitsEntry(Gint, GintOrdinal, Varray, g, _k));
    else {
	assert(_outputMode == indexMotifs);
	puts(PrintIndexEntry(Gint, GintOrdinal, Varray, g, _k));
    }

    // Now go about deleting edges recursively.
    int i,j;
    for(i=0; i<_k-1; i++)for(j=i+1;j<_k;j++)
    {
	if(TinyGraphAreConnected(g,i,j)) // if it's an edge, delete it.
	{
	    TinyGraphDisconnect(g,i,j);
	    Gint = TinyGraph2Int(g,_k);
	    GintOrdinal = _K[Gint];
	    if(SetIn(_connectedCanonicals, GintOrdinal)) {
		++depth;
		PrintAllMotifs(g,Gint,GintOrdinal, G,Varray);
		--depth;
	    }
	    TinyGraphConnect(g,i,j);
	}
    }
}

Boolean ProcessGraphlet(GRAPH *G, SET *V, unsigned Varray[], const int k, TINY_GRAPH *g)
{
    Boolean processed = true;
    TinyGraphInducedFromGraph(g, G, Varray);
    int Gint = TinyGraph2Int(g,k), GintOrdinal=_K[Gint], j;

#if PARANOID_ASSERTS
    assert(0 <= GintOrdinal && GintOrdinal < _numCanon);
#endif
    switch(_outputMode)
    {
	char perm[MAX_K];
    case graphletFrequency:
	++_graphletCount[GintOrdinal];
	break;
    case indexGraphlets:
#if SORT_INDEX_MODE // Note this destroys the columns-are-identical property, don't use by default.
	VarraySort(Varray, k);
#endif
	if(NodeSetSeenRecently(G, Varray,k) ||
	    (_sampleMethod == SAMPLE_INDEX && !SetIn(_windowRep_allowed_ambig_set, GintOrdinal))) processed=false;
	else puts(PrintIndexEntry(Gint, GintOrdinal, Varray, g, k));
	break;
    case indexMotifs: case indexMotifOrbits:
	if(NodeSetSeenRecently(G,Varray,k)) processed=false;
	else PrintAllMotifs(g,Gint,GintOrdinal, G,Varray);
	break;
#if PREDICT
    case predict:
	if(NodeSetSeenRecently(G,Varray,k)) processed=false;
	else AccumulateGraphletParticipationCounts(G,Varray,g,Gint,GintOrdinal);
	break;
#endif
    case indexOrbits:
#if SORT_INDEX_MODE // Note this destroys the columns-are-identical property, don't use by default.
	VarraySort(Varray, k);
#endif
	if(NodeSetSeenRecently(G,Varray,k) ||
	    (_sampleMethod == SAMPLE_INDEX && !SetIn(_windowRep_allowed_ambig_set, GintOrdinal))) processed=false;
	else puts(PrintIndexOrbitsEntry(Gint, GintOrdinal, Varray, g, k));
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
