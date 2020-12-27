#include "blant.h"
#include "blant-output.h"
#include "blant-predict.h"
#include "blant-utils.h"
#include "blant-sampling.h"

void PrintNode(char c, int v) {
    if(c) putchar(c);
#if SHAWN_AND_ZICAN
    printf("%s", _nodeNames[v]);
#else
    if(_supportNodeNames)
	printf("%s", _nodeNames[v]);
    else
	printf("%d", v);
#endif
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

void PrintCanonical(int GintOrdinal)
{
    int j, GintNumBits;
    char GintBinary[GintNumBits+1]; // Only used in -db output mode for indexing
    switch (_displayMode) {
    case undefined:
    case ordinal:
	printf("%d", GintOrdinal);
	break;
    case decimal: // Prints the decimal integer form of the canonical
	printf("%d", _canonList[GintOrdinal]);
	break;
    case binary: // Prints the bit representation of the canonical
	GintNumBits = _k*(_k-1)/2;
	for (j=0;j<GintNumBits;j++)
	    {GintBinary[GintNumBits-j-1]=(((unsigned)_canonList[GintOrdinal] >> j) & 1 ? '1' : '0');}
	GintBinary[GintNumBits] = '\0';
	printf("%s", GintBinary);
	break;
    case orca: // Prints the ORCA ID of the canonical. Jesse uses same number.
    case jesse:
	printf("%d", _outputMapping[GintOrdinal]);
	break;
    }
}

// Below is code to help reduce (mostly eliminite if we're lucky) MCMC's duplicate output, which is copious
// Empirically I've found that neither of these need to be very big since most repetition happens with high locality
#define MCMC_PREDICT_CIRC_BUF 999983 // prime, not sure it needs to be but why not... 1M * 4b = 4MB RAM.
#define MCMC_PREDICT_MAX_HASH 2147483647 // this value should be prime; at 2^31/8 it's 256MB of RAM.
// NOTE WE DO NOT CHECK EDGES. So if you call it with the same node set but as a motif, it'll (incorrectly) return TRUE
Boolean NodeSetSeenRecently(GRAPH *G, unsigned Varray[], int k) {
    static unsigned circBuf[MCMC_PREDICT_CIRC_BUF], bufPos;
    static SET *seen;
    static unsigned Vcopy[maxK];
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

void PrintIndexEntry(int Gint, int GintOrdinal, unsigned Varray[], char perm[], TINY_GRAPH *g, int k)
{
    int j;
    memset(perm, 0, k);
    ExtractPerm(perm, Gint);
    PrintCanonical(GintOrdinal);
    assert(PERMS_CAN2NON);
    for(j=0;j<k;j++) {
	PrintNode(' ', Varray[(int)perm[j]]);
    }
    puts("");
}

void PrintIndexOrbitsEntry(int Gint, int GintOrdinal, unsigned Varray[], char perm[], TINY_GRAPH *g, int k) {
    assert(TinyGraphDFSConnected(g,0));
    int j;
    static SET* printed;
    if(!printed) printed = SetAlloc(k);
    SetEmpty(printed);
    memset(perm, 0, k);
    ExtractPerm(perm, Gint);
    PrintCanonical(GintOrdinal);
    assert(PERMS_CAN2NON); // Apology("Um, don't we need to check PERMS_CAN2NON? See outputODV for correct example");
    for(j=0;j<k;j++) if(!SetIn(printed,j))
    {
	PrintNode(' ', Varray[(int)perm[j]]);
	SetAdd(printed, j);
	int j1;
	for(j1=j+1;j1<k;j1++) if(_orbitList[GintOrdinal][j1] == _orbitList[GintOrdinal][j])
	{
	    assert(!SetIn(printed, j1));
	    PrintNode(':', Varray[(int)perm[j1]]);
	    SetAdd(printed, j1);
	}
    }
    putchar('\n');
}

// Recursively print all the motifs under this graphlet
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

    char perm[maxK];
    memset(perm, 0, _k);
    if(_outputMode == indexMotifOrbits)
		PrintIndexOrbitsEntry(Gint, GintOrdinal, Varray, perm, g, _k);
    else {
		assert(_outputMode == indexMotifs);
		PrintIndexEntry(Gint, GintOrdinal, Varray, perm, g, _k);
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

Boolean ProcessGraphlet(GRAPH *G, SET *V, unsigned Varray[], const int k, char perm[], TINY_GRAPH *g)
{
    Boolean processed = true;
    TinyGraphInducedFromGraph(g, G, Varray);
    int Gint = TinyGraph2Int(g,k), j, GintOrdinal=_K[Gint];

#if PARANOID_ASSERTS
    assert(0 <= GintOrdinal && GintOrdinal < _numCanon);
#endif
    switch(_outputMode)
    {
    case graphletFrequency:
	++_graphletCount[GintOrdinal];
	break;
    case indexGraphlets:
#if SORT_INDEX_MODE // Note this destroys the columns-are-identical property, don't use by default.
	VarraySort(Varray, k);
#endif
	if(NodeSetSeenRecently(G, Varray,k) || _sampleMethod == SAMPLE_INDEX && !SetIn(_windowRep_allowed_ambig_set, GintOrdinal)) processed=false;
	else PrintIndexEntry(Gint, GintOrdinal, Varray, perm, g, k);
	break;
    case indexMotifs: case indexMotifOrbits:
	if(NodeSetSeenRecently(G,Varray,k)) processed=false;
	else PrintAllMotifs(g,Gint,GintOrdinal, G,Varray);
	break;
    case predict:
	if(NodeSetSeenRecently(G,Varray,k)) processed=false;
	else CountSubmotifOrbitPairs(g,G,Varray);
	break;
	break;
    case indexOrbits:
#if SORT_INDEX_MODE // Note this destroys the columns-are-identical property, don't use by default.
	VarraySort(Varray, k);
#endif
	if(NodeSetSeenRecently(G,Varray,k) || _sampleMethod == SAMPLE_INDEX && !SetIn(_windowRep_allowed_ambig_set, GintOrdinal)) processed=false;
	else PrintIndexOrbitsEntry(Gint, GintOrdinal, Varray, perm, g, k);
	break;
    case outputGDV:
	for(j=0;j<k;j++) ++GDV(Varray[j], GintOrdinal);
	break;
    case outputODV:
	memset(perm, 0, k);
	ExtractPerm(perm, Gint);
#if PERMS_CAN2NON            
	for(j=0;j<k;j++) ++ODV(Varray[(int)perm[j]], _orbitList[GintOrdinal][          j ]);
#else
	for(j=0;j<k;j++) ++ODV(Varray[          j ], _orbitList[GintOrdinal][(int)perm[j]]);
#endif
	break;
	    
    default: Abort("ProcessGraphlet: unknown or un-implemented outputMode");
	break;
    }
    return processed;
}
