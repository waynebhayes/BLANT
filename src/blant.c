#include <unistd.h>
#include <sys/wait.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <limits.h>
#include "misc.h"
#include "tinygraph.h"
#include "graph.h"
#include "heap.h"
#include "blant.h"
#include "queue.h"
#include "multisets.h"
#include "sorts.h"
#include "blant-window.h"
#include "blant-output.h"
#include "blant-utils.h"
#include "blant-sampling.h"
#include "blant-synth-graph.h"
#include "rand48.h"
Boolean _earlyAbort; // Can be set true by anybody anywhere, and they're responsible for producing a warning as to why
#include "blant-predict.h"
#include "importance.h"
#include "odv.h"

static int *_pairs, _numNodes, _numEdges, _maxEdges=1024, _seed = -1; // -1 means "not initialized"
char **_nodeNames, _supportNodeNames = true;
Boolean _child; // are we a child process?

char * _sampleFileName;

// _k is the global variable storing k; _Bk=actual number of entries in the canon_map for given k.
unsigned int _k;
unsigned int _Bk, _k_small;

int _alphaList[MAX_CANONICALS];
int _numCanon, _numSamples;
Gint_type _canonList[MAX_CANONICALS]; // map ordinals to integer representation of the canonical
SET *_connectedCanonicals; // the SET of canonicals that are connected.
int _numConnectedCanon;
int _numConnectedComponents;
int *_componentSize;

int _numOrbits, _orbitList[MAX_CANONICALS][MAX_K]; // map from [ordinal][canonicalNode] to orbit ID.
int _orbitCanonMapping[MAX_ORBITS]; // Maps orbits to canonical (including disconnected)
int _orbitCanonNodeMapping[MAX_ORBITS]; // Maps orbits to canonical (including disconnected)
int *_whichComponent;

// char* _BLANT_DIR;

enum OutputMode _outputMode = undef;
unsigned long int _graphletCount[MAX_CANONICALS];
int **_graphletDistributionTable;
double _g_overcount, _graphletConcentration[MAX_CANONICALS];

enum CanonicalDisplayMode _displayMode = undefined;
enum FrequencyDisplayMode _freqDisplayMode = freq_display_mode_undef;

int _outputMapping[MAX_CANONICALS];

int _orca_orbit_mapping[58];
int _connectedOrbits[MAX_ORBITS];
int _numConnectedOrbits;

// A bit counter-intuitive: we need to allocate this many vectors each of length [_numNodes],
// and then the degree for node v, graphlet/orbit g is _degreeVector[g][v], NOT [v][g].
// We do this simply because we know the length of MAX_CANONICALS so we pre-know the length of
// the first dimension, otherwise we'd need to get more funky with the pointer allocation.
// Only one of these actually get allocated, depending upon outputMode.
unsigned long int *_graphletDegreeVector[MAX_CANONICALS];
unsigned long int    *_orbitDegreeVector[MAX_ORBITS];
double *_doubleOrbitDegreeVector[MAX_ORBITS];

double *_cumulativeProb;

// number of parallel threads required, and the maximum allowed at one time.
int _JOBS, _MAX_THREADS;

// Here's the actual mapping from non-canonical to canonical, same argument as above wasting memory, and also mmap'd.
// So here we are allocating 256MB x sizeof(short int) = 512MB.
// Grand total statically allocated memory is exactly 1.25GB.
//static short int _K[maxBk] __attribute__ ((aligned (8192)));
short *_K = NULL; // Allocating memory dynamically

/* AND NOW THE CODE */

// return how many nodes found. If you call it with startingNode == 0 then we automatically clear the visited array
static TSET _visited;
static int NumReachableNodes(TINY_GRAPH *g, int startingNode)
{
    if(startingNode == 0) TSetEmpty(_visited);
    TSetAdd(_visited,startingNode);
    unsigned int j, Varray[MAX_K], numVisited = 0;
    int numNeighbors = TSetToArray(Varray, g->A[startingNode]);
    assert(numNeighbors == g->degree[startingNode]);
    for(j=0; j<numNeighbors; j++)if(!TSetIn(_visited,Varray[j])) numVisited += NumReachableNodes(g,Varray[j]);
    return 1+numVisited;
}

static int **_componentList; // list of lists of components, largest to smallest.
static double _totalCombinations, *_combinations, *_probOfComponent;
SET **_componentSet;

void SetBlantDir() {
    char* temp = getenv("BLANT_DIR");
    if (temp != NULL)
	_BLANT_DIR = strdup(temp); // can't assume the string returned by getetv never changes, so copy it.
}

static int InitializeConnectedComponents(GRAPH *G)
{
    static unsigned int v, *Varray, j, i;
    assert(!Varray); // we only can be called once.
    assert(_numConnectedComponents == 0);
    SET *visited = SetAlloc(G->n);
    Varray = Calloc(G->n, sizeof(int));
    _whichComponent = Calloc(G->n, sizeof(int));
    _componentSize = Calloc(G->n, sizeof(int)); // probably bigger than it needs to be but...
    _componentList = Calloc(G->n, sizeof(int*)); // probably bigger...
    _combinations = Calloc(G->n, sizeof(double*)); // probably bigger...
    _probOfComponent = Calloc(G->n, sizeof(double*)); // probably bigger...
    _cumulativeProb = Calloc(G->n, sizeof(double*)); // probably bigger...
    _componentSet = Calloc(G->n, sizeof(SET*));

    int nextStart = 0;
    _componentList[0] = Varray;
    for(v=0; v < G->n; v++) if(!SetIn(visited, v))
    {
	_componentSet[_numConnectedComponents] = SetAlloc(G->n);
	_componentList[_numConnectedComponents] = Varray + nextStart;
	GraphVisitCC(G, v, _componentSet[_numConnectedComponents], Varray + nextStart, _componentSize + _numConnectedComponents);
	SetUnion(visited, visited, _componentSet[_numConnectedComponents]);
	for(j=0; j < _componentSize[_numConnectedComponents]; j++)
	{
	    assert(_whichComponent[Varray[nextStart + j]] == 0);
	    _whichComponent[Varray[nextStart + j]] = _numConnectedComponents;
	}
	nextStart += _componentSize[_numConnectedComponents];
	++_numConnectedComponents;
    }
    assert(nextStart == G->n);

    _totalCombinations = 0.0;
    for(i=0; i< _numConnectedComponents; i++)
    {
	//find the biggest one
	int biggest = i;
	for(j=i+1; j<_numConnectedComponents;j++)
	    if(_componentSize[j] > _componentSize[biggest])
		biggest = j;
	// Now swap the biggest one into position i;
	for(j=0; j < _componentSize[biggest]; j++)
	    _whichComponent[_componentList[biggest][j]] = i;
	int itmp, *pitmp;
	SET * stmp;
	itmp = _componentSize[i];
	_componentSize[i] = _componentSize[biggest];
	_componentSize[biggest] = itmp;
	pitmp = _componentList[i];
	_componentList[i] = _componentList[biggest];
	_componentList[biggest] = pitmp;
    stmp = _componentSet[i];
    _componentSet[i] = _componentSet[biggest];
    _componentSet[biggest] = stmp;
	_combinations[i] = CombinChooseDouble(_componentSize[i], _k);
	_totalCombinations += _combinations[i];
    }

    double cumulativeProb = 0.0;
    for(i=0; i< _numConnectedComponents; i++)
    {
	_probOfComponent[i] =  _combinations[i] / _totalCombinations;
	_cumulativeProb[i] = cumulativeProb + _probOfComponent[i];
	cumulativeProb = _cumulativeProb[i];
	if(cumulativeProb > 1)
	{
	    assert(cumulativeProb - 1.0 < _numConnectedComponents * 1e-15);  // allow some roundoff error
	    cumulativeProb = _cumulativeProb[i] = 1.0;
	}
	//printf("Component %d has %d nodes and probability %lf, cumulative prob %lf\n", i, _componentSize[i], _probOfComponent[i], _cumulativeProb[i]);
    }
    SetFree(visited);
    return _numConnectedComponents;
}

int alphaListPopulate(char *BUF, int *alpha_list, int k) {
	sprintf(BUF, "%s/%s/alpha_list_mcmc%d.txt", _BLANT_DIR, CANON_DIR, k);
    FILE *fp_ord=fopen(BUF, "r");
    if(!fp_ord) Fatal("cannot find %s\n", BUF);
    int numAlphas, i;
    assert(1==fscanf(fp_ord, "%d",&numAlphas));
	assert(numAlphas == _numCanon);
    for(i=0; i<numAlphas; i++) assert(1==fscanf(fp_ord, "%d", &alpha_list[i]));
    fclose(fp_ord);
    return numAlphas;
}

// Loads alpha values(overcounting ratios) for MCMC sampling from files
// The alpha value represents the number of ways to walk over that graphlet
// Global variable _MCMC_L is needed by many functions
// _MCMC_L represents the length of the sliding window in d graphlets for sampling
// Global variable _numSamples needed for the algorithm to reseed halfway through
// Concentrations are initialized to 0
void initializeMCMC(GRAPH* G, int k, int numSamples) {
	_MCMC_L = k - mcmc_d  + 1;
	// Count the number of valid edges to start from
	int i, validEdgeCount = 0;
	for (i = 0; i < G->numEdges; i++)
	    if (_componentSize[_whichComponent[G->edgeList[2*i]]] >= k)
		validEdgeCount++;
	_samplesPerEdge = (numSamples + (validEdgeCount / 2)) / validEdgeCount; // Division rounding up for samples per edge
	if(_sampleSubmethod == SAMPLE_MCMC_EC) {
	    //_samplesPerEdge = numSamples;
	    _EDGE_COVER_G = GraphCopy(NULL, G);
	}

	char BUF[BUFSIZ];
	_numSamples = numSamples;
	if(!_window)
	{
		alphaListPopulate(BUF, _alphaList, k);
		for (i = 0; i < _numCanon; i++)
		{
			_graphletConcentration[i] = 0.0;
		}
	}
}

// Convert the graphlet frequencies to concentrations
void finalizeMCMC() {
	double totalConcentration = 0;
	int i;
	for (i = 0; i < _numCanon; i++) {
#if PARANOID_ASSERTS
	    if(_graphletConcentration[i] < 0.0) {
		Warning("_graphletConcentration[%d] is %g\n",i, _graphletConcentration[i]);
		assert(false);
	    }
#endif
	    totalConcentration += _graphletConcentration[i];
	}
	for (i = 0; i < _numCanon; i++) {
	    _graphletConcentration[i] /= totalConcentration;
	}
}

// Compute the degree of the state in the state graph (see Lu&Bressen)
// Given the big graph G, and a set of nodes S (|S|==k), compute the
// degree of the *state* represented by these nodes, which is:
//     Degree(S) = k*|NN2| + (k-1)*|NN1|
// where NN1 is the set of nodes one step outside S that have only 1 connection back into S,
// and   NN2 is the set of nodes one step outside S that have more than 1 connection back into S.
static int StateDegree(GRAPH *G, SET *S)
{
#if PARANOID_ASSERTS
    assert(SetCardinality(S) == _k);
#endif
    int Varray[_k]; // array of elements in S
    SetToArray(Varray, S);

    SET *outSet=SetAlloc(G->n); // the set of nodes in G that are one step outside S, from *anybody* in S
    int connectCount[G->n]; // for each node in the outset, count the number of times it's connected into S
    memset(connectCount, 0, G->n * sizeof(int)); // set them all to zero

    int i, j, numOut = 0;
    for(i=0;i<_k;i++) for(j=0; j<G->degree[Varray[i]]; j++)
    {
	int neighbor = G->neighbor[Varray[i]][j];
	if(!SetIn(S, neighbor))
	{
	    SetAdd(outSet, neighbor);
	    if(connectCount[neighbor] == 0) ++numOut;
	    ++connectCount[neighbor];
	}
    }
#if PARANOID_ASSERTS
    assert(SetCardinality(outSet) == numOut);
#endif
    int outArray[numOut], Sdeg = 0;
    i = SetToArray(outArray, outSet);
    assert(i == numOut);
    for(i=0; i < numOut; i++)
    {
	switch(connectCount[outArray[i]])
	{
	case 0: break;
	case 1: Sdeg += _k-1; break;
	default:  Sdeg += _k; break;
	}
    }
    SetFree(outSet);
    return Sdeg;
}


// This converts graphlet frequencies to concentrations or integers based on the sampling algorithm and command line arguments
void convertFrequencies(int numSamples)
{
    int i;
    if (_sampleMethod == SAMPLE_MCMC) {
	if (_freqDisplayMode == count) {
	    for (i = 0; i < _numCanon; i++) {
		_graphletCount[i] = _graphletConcentration[i] * numSamples;
	    }
	}
    }
    else {
	if (_freqDisplayMode == concentration) {
	    for (i = 0; i < _numCanon; i++) {
		_graphletConcentration[i] = _graphletCount[i] / (double)numSamples;
	    }
	}
    }
}


// This is the single-threaded BLANT function. YOU PROBABLY SHOULD NOT CALL THIS.
// Call RunBlantInThreads instead, it's the top-level entry point to call once the
// graph is finished being input---all the ways of reading input call RunBlantInThreads.
// Note it does stuff even if numSamples == 0, because we may be the parent of many
// threads that finished and we have nothing to do except output their accumulated results.
int RunBlantFromGraph(int k, int numSamples, GRAPH *G)
{
    int i,j, windowRepInt, D;
    char perm[MAX_K+1];
    assert(k <= G->n);
    SET *V = SetAlloc(G->n);
    SET *prev_node_set = SetAlloc(G->n);
    SET *intersect_node = SetAlloc(G->n);
    TINY_GRAPH *empty_g = TinyGraphAlloc(k); // allocate it here once, so functions below here don't need to do it repeatedly
    int varraySize = _windowSize > 0 ? _windowSize : MAX_K + 1;
    unsigned Varray[varraySize];
    InitializeConnectedComponents(G);
    if (_sampleMethod == SAMPLE_MCMC)
	_window? initializeMCMC(G, _windowSize, numSamples) : initializeMCMC(G, k, numSamples);
    if (_outputMode == graphletDistribution) {
        SampleGraphlet(G, V, Varray, k, G->n);
        SetCopy(prev_node_set, V);
        TinyGraphInducedFromGraph(empty_g, G, Varray);
    }

    if ((_sampleMethod == SAMPLE_INDEX || _sampleSubmethod == SAMPLE_MCMC_EC) &&
	(_outputMode != indexGraphlets && _outputMode != indexGraphletsRNO && _outputMode != indexOrbits))
	    Fatal("currently only -mi and -mj output modes are supported for INDEX and EDGE_COVER sampling methods");

    if (_sampleMethod == SAMPLE_INDEX) {
        int i, count = 0;
        int prev_nodes_array[_k];

        // Get heuristic values based on orbit number, if ODV file provided
        double heuristicValues[G->n];

        if (_orbitNumber != -1) {
            getOdvValues(heuristicValues, _orbitNumber, _nodeNames, G->n);
        } else {
            getDoubleDegreeArr(heuristicValues, G); // since heuristic values are doubles, we need to convert degree values to doubles
        }

        int percentToPrint = 1;
        node_whn nwhn_arr[G->n]; // nodes sorted first by the heuristic function and then either alphabetically or reverse alphabetically

        // fill node order array with base values
        for (i = 0; i < G->n; i++) {
            nwhn_arr[i].node = i;
            nwhn_arr[i].heur = heuristicValues[i];
            nwhn_arr[i].name = _nodeNames[i];
        }

        // sort array
        int (*comp_func)(const void*, const void*);

        if (_alphabeticTieBreaking) {
            comp_func = nwhn_des_alph_comp_func;
        } else {
            comp_func = nwhn_des_rev_comp_func;
        }

        qsort((void*)nwhn_arr, G->n, sizeof(node_whn), comp_func);

        for(i=0; i<G->n; i++) {
            prev_nodes_array[0] = nwhn_arr[i].node;

            SampleGraphletIndexAndPrint(G, prev_nodes_array, 1, heuristicValues);
            count = 0;

            if (i * 100 / G->n >= percentToPrint) {
                fprintf(stderr, "%d%% done\n", percentToPrint);
                ++percentToPrint;
            }
	    }
    }
    else if (_sampleMethod == SAMPLE_MCMC_EC) {
	Fatal("should not get here--EDGE_COVER is a submethod of MCMC");
#if 0
	int e;
	for(e=0; e<G->numEdges; e++) if(SetIn(needEdge, e)) {
	    int whichCC = -(e+1); // encoding edge number in whichCC
	    SampleGraphlet(G, V, Varray, k, whichCC);
	    ProcessGraphlet(G, V, Varray, k, empty_g);

	    // Now remove all the edges in the graphlet from "needEdge" (SLOW AND DUMB but it's not really a problem)
	    int i,j,f;
	    for(i=0;i<k;i++) for(j=i+1;j<k;j++) {
		int u=Varray[i], v=Varray[j];
		if(GraphAreConnected(G,u,v)) { // find (u,v) in the edgeList and remove it from needEdge
		    Boolean found=false;
		    for(f=0;f<G->numEdges;f++) {
			if((G->edgeList[2*f]==u && G->edgeList[2*f+1]==v) || (G->edgeList[2*f]==v && G->edgeList[2*f+1]==u)) {
			    found=true;
			    SetDelete(needEdge, f);
			    break;
			}
		    }
		    assert(found);
		}
	    }
	}
#endif
    }
    else // sample numSamples graphlets for the entire graph
    {
        for(i=0; (i<numSamples || (_sampleFile && !_sampleFileEOF)) && !_earlyAbort; i++)
        {
            if(_window) {
                SampleGraphlet(G, V, Varray, _windowSize, G->n);
                _numWindowRep = 0;
                if (_windowSampleMethod == WINDOW_SAMPLE_MIN || _windowSampleMethod == WINDOW_SAMPLE_MIN_D ||
			_windowSampleMethod == WINDOW_SAMPLE_LEAST_FREQ_MIN)
                    windowRepInt = getMaximumIntNumber(_k);
                if (_windowSampleMethod == WINDOW_SAMPLE_MAX || _windowSampleMethod == WINDOW_SAMPLE_MAX_D ||
			_windowSampleMethod == WINDOW_SAMPLE_LEAST_FREQ_MAX)
                    windowRepInt = -1;
                D = _k * (_k - 1) / 2;
                FindWindowRepInWindow(G, V, &windowRepInt, &D, perm);
                if(_numWindowRep > 0)
                    ProcessWindowRep(G, Varray, windowRepInt);
            }
            else if (_outputMode == graphletDistribution)
                ProcessWindowDistribution(G, V, Varray, k, empty_g, prev_node_set, intersect_node);
            else {
		static int stuck;
		// HACK: make the graphlet overcount global; it should really be PASSED into ProcessGraphlet
                _g_overcount = SampleGraphlet(G, V, Varray, k, G->n); // weight will be 1.0 in most cases but if sample method is MCMC and it's not windowed it will be the count of the graphlet
                if(ProcessGraphlet(G, V, Varray, k, empty_g)) stuck = 0;
		else {
		    --i; // negate the sample count of duplicate graphlets
		    ++stuck;
		    if(stuck > numSamples) {
			Warning("Sampling aborted: no new graphlets discovered after %d attempts", stuck);
			_earlyAbort = true;
		    }
		}
            }
        }
	if(i<numSamples) Warning("only took %d samples out of %d", i, numSamples);
    }

    // Sampling done. Now generate output for output modes that require it.

    if(_window) {
        for(i=0; i<_numWindowRepArrSize; i++)
            free(_windowReps[i]);
        free(_windowReps);
        if(_windowRep_limit_method) HeapFree(_windowRep_limit_heap);
    }
    if (_sampleMethod == SAMPLE_MCMC && !_window)
	finalizeMCMC();
    if (_outputMode == graphletFrequency && !_window)
	convertFrequencies(numSamples);

    switch(_outputMode)
    {
	int canon;
	int orbit_index;
    case indexGraphlets: case indexGraphletsRNO: case indexOrbits: case indexMotifs: case indexMotifOrbits:
	break; // already printed on-the-fly in the Sample/Process loop above
    case graphletFrequency:
	for(canon=0; canon<_numCanon; canon++) {
	if (_freqDisplayMode == concentration) {
	    if (SetIn(_connectedCanonicals, canon)) {
		printf("%lf ", _graphletConcentration[canon]);
		puts(PrintCanonical(canon));
	    }
	}
	else {
	    if (SetIn(_connectedCanonicals, canon)) {
		printf("%lu ", _graphletCount[canon]);
		puts(PrintCanonical(canon));
	    }
	}
	}
	break;
    case predict_merge:
	assert(false); break; // shouldn't get here
    case predict:
	Predict_FlushMotifs(G);
	break;
    case outputGDV:
	for(i=0; i < G->n; i++)
	{
	    printf("%s", PrintNode(0,i));
	    for(canon=0; canon < _numCanon; canon++)
		printf(" %lu", GDV(i,canon));
	    puts("");
	}
	break;
    case outputODV:
        for(i=0; i<G->n; i++) {
	    printf("%s", PrintNode(0,i));
	    for(j=0; j<_numConnectedOrbits; j++) {
		if (k == 4 || k == 5) orbit_index = _connectedOrbits[_orca_orbit_mapping[j]];
		else orbit_index = _connectedOrbits[j];
		if (!_MCMC_EVERY_EDGE || _sampleMethod != SAMPLE_MCMC) printf(" %lu", ODV(i,orbit_index));
		else printf(" %.12f", _doubleOrbitDegreeVector[orbit_index][i]);
	    }
	    printf("\n");
	}
        break;
    case graphletDistribution:
        for(i=0; i<_numCanon; i++) {
            for(j=0; j<_numCanon; j++)
                printf("%d ", _graphletDistributionTable[i][j]);
            printf("\n");
        }
        break;
    default: Abort("RunBlantFromGraph: unknown or un-implemented outputMode");
	break;
	}

#if PARANOID_ASSERTS // no point in freeing this stuff since we're about to exit; it can take significant time for large graphs.
    if(_outputMode == outputGDV) for(i=0;i<_numCanon;i++)
	Free(_graphletDegreeVector[i]);
    if(_outputMode == outputODV) for(i=0;i<_numOrbits;i++) Free(_orbitDegreeVector[i]);
	if(_outputMode == outputODV && _MCMC_EVERY_EDGE) for(i=0;i<_numOrbits;i++) Free(_doubleOrbitDegreeVector[i]);
    TinyGraphFree(empty_g);
#endif
    if (_sampleMethod == SAMPLE_ACCEPT_REJECT)
    	fprintf(stderr,"Average number of tries per sample is %g\n", _acceptRejectTotalTries/(double)numSamples);
    SetFree(V);
    SetFree(prev_node_set);
    SetFree(intersect_node);
    return 0;
}

/*
** Fork a BLANT process and return a FILE pointer where it'll be sending stuff.
** Caller is responsible for reading all the stuff from the returned FILE pointer,
** detecting EOF on it, and fclose'ing it.
*/
FILE *ForkBlant(int k, int numSamples, GRAPH *G)
{
    int fds[2];
    assert(pipe(fds) >= 0);
    int threadSeed = INT_MAX*RandomUniform(); // this must happen BEFORE the fork for each thread to get a different seed
    int pid = fork();
    if(pid > 0) // we are the parent
    {
	close(fds[1]); // we will not be writing to the pipe, so close it.
	return fdopen(fds[0],"r");
    }
    else if(pid == 0) // we are the child
    {
	_child = true;
	_seed = threadSeed;
	RandomSeed(_seed);
	(void)close(fds[0]); // we will not be reading from the pipe, so close it.
	(void)close(1); // close our usual stdout
	assert(dup(fds[1])>=0); // copy the write end of the pipe to fd 1.
	(void)close(fds[1]); // close the original write end of the pipe since it's been moved to fd 1.

	// For any "counting" mode, use internal numbering when communicating through pipes to the parent
	if(_outputMode != indexGraphlets && _outputMode != indexGraphletsRNO && _outputMode != indexOrbits) _supportNodeNames = false;

	RunBlantFromGraph(k, numSamples, G);
	exit(0);
	_exit(0);
	Abort("Both exit() and _exit failed???");
    }
    else
	Abort("fork failed");
    assert(false);
    return NULL; //  should never get here, but quell compiler warning.
}


static FILE *fpThreads[MAX_POSSIBLE_THREADS]; // these will be the pipes reading output of the parallel blants

// This is the primary entry point into BLANT, even if THREADS=1.  We assume you've already
// read the graph into G, and will do whatever is necessary to run blant with the number of
// threads specified.  Also does some sanity checking.
int RunBlantInThreads(int k, int numSamples, GRAPH *G)
{
    int i,j;
    assert(k == _k);
    assert(G->n >= k); // should really ensure at least one connected component has >=k nodes. TODO
    if(_outputMode == outputGDV) for(i=0;i<_numCanon;i++)
	_graphletDegreeVector[i] = Calloc(G->n, sizeof(**_graphletDegreeVector));
    if(_outputMode == outputODV) for(i=0;i<_numOrbits;i++){
	_orbitDegreeVector[i] = Calloc(G->n, sizeof(**_orbitDegreeVector));
	for(j=0;j<G->n;j++) _orbitDegreeVector[i][j]=0;
    }
    if (_outputMode == outputODV) for(i=0;i<_numOrbits;i++){
	_doubleOrbitDegreeVector[i] = Calloc(G->n, sizeof(**_doubleOrbitDegreeVector));
	for(j=0;j<G->n;j++) _doubleOrbitDegreeVector[i][j]=0.0;
    }
    if(_outputMode == predict) Predict_Init(G);
    if (_outputMode == graphletDistribution) {
        _graphletDistributionTable = Calloc(_numCanon, sizeof(int*));
        for(i=0; i<_numCanon; i++) _graphletDistributionTable[i] = Calloc(_numCanon, sizeof(int));
        for(i=0; i<_numCanon; i++) for(j=0; j<_numCanon; j++) _graphletDistributionTable[i][j] = 0;
    }

    if(_JOBS == 1)
	return RunBlantFromGraph(k, numSamples, G);

    if (_sampleMethod == SAMPLE_INDEX || _sampleSubmethod == SAMPLE_MCMC_EC)
        Fatal("The sampling methods INDEX and EDGE_COVER do not yet support multithreading (feel free to add it!)");

    // At this point, _JOBS must be greater than 1.
    assert(_JOBS>1);
    int totalSamples = numSamples;
    double meanSamplesPerJob = totalSamples/(double)_JOBS;
    Warning("Parent %d starting about %d jobs of about %d samples each", getpid(), _JOBS, (int)meanSamplesPerJob);

    int threadsRunning = 0, jobsDone = 0;
    int thread, lineNum = 0, job=0;
    for(i=0; numSamples > 0 && i<_MAX_THREADS;i++) {
	int samples = meanSamplesPerJob;
	assert(samples>0);
	if(samples > numSamples) samples = numSamples;
	numSamples -= samples;
	fpThreads[i] = ForkBlant(_k, samples, G);
	Warning("Started job %d requesting %d samples; %d threads running, %d samples remaining to take",
	    job++, samples, ++threadsRunning, numSamples);
    }
    do
    {
	char line[MAX_ORBITS * BUFSIZ];
	for(thread=0;thread<_MAX_THREADS;thread++)	// process one line from each thread
	{
	    if(!fpThreads[thread]) continue;
	    char *tmp = fgets(line, sizeof(line), fpThreads[thread]);
	    assert(tmp>=0);
	    if(feof(fpThreads[thread]))
	    {
		wait(NULL); // reap the child
		clearerr(fpThreads[thread]);
		fclose(fpThreads[thread]);
		fpThreads[thread] = NULL;
		++jobsDone; --threadsRunning;
		Warning("Thead %d finished; jobsDone %d, threadsRunning %d", thread, jobsDone, threadsRunning);
		if(numSamples == 0) fpThreads[thread] = NULL; // signify this pointer is finished.
		else {
		    int samples = meanSamplesPerJob;
		    if(samples > numSamples) samples = numSamples;
		    numSamples -= samples;
		    fpThreads[thread] = ForkBlant(_k, samples, G);
		    assert(fpThreads[thread]);
		    ++threadsRunning;
		    Warning("Started job %d (thread %d) requesting %d samples, %d threads running, %d samples remaining to take",
			job++, thread, samples, threadsRunning, numSamples);
		}
		continue; // we'll ask for output next time around.
	    }
	    char *nextChar = line, *pch;
	    unsigned long int count;
	    int canon, orbit, numRead, nodeId, value;
	    float fValue;
	    //fprintf(stderr, "Parent received the following line from the child: <%s>\n", line);
	    switch(_outputMode)
	    {
	    case graphletFrequency:
		numRead = sscanf(line, "%lu%d", &count, &canon);
		assert(numRead == 2);
		_graphletCount[canon] += count;
		break;
	    case graphletDistribution:
		for(i=0; i<_numCanon; i++) {
		    pch = strtok(line, " ");
		    for(j=0; j<_numCanon; j++){
			_graphletDistributionTable[i][j] += atoi(pch);
			pch = strtok(NULL, " ");
		    }
		    char *OK = fgets(line, sizeof(line), fpThreads[thread]);
		    assert(OK || (i==_numCanon-1 && j == _numCanon));
		}
		break;
	    case outputGDV:
		assert(isdigit(*nextChar));
		numRead = sscanf(nextChar, "%d", &nodeId);
		assert(numRead == 1 && nodeId == lineNum);
		while(isdigit(*nextChar)) nextChar++; // read past current integer
		assert(*nextChar == ' ' || (canon == _numCanon-1 && *nextChar == '\n'));
		nextChar++;
		for(canon=0; canon < _numCanon; canon++)
		{
		    assert(isdigit(*nextChar));
		    numRead = sscanf(nextChar, "%lu", &count);
		    assert(numRead == 1);
		    GDV(lineNum,canon) += count;
		    while(isdigit(*nextChar)) nextChar++; // read past current integer
		    assert(*nextChar == ' ' || (canon == _numCanon-1 && *nextChar == '\n'));
		    nextChar++;
		}
		assert(*nextChar == '\0');
		break;
	    case outputODV:
		assert(isdigit(*nextChar));
		numRead = sscanf(nextChar, "%d", &nodeId);
		assert(numRead == 1 && nodeId == lineNum);
		while(isdigit(*nextChar)) nextChar++; // read past current integer
		assert(*nextChar == ' ' || (orbit == _numOrbits-1 && *nextChar == '\n'));
		nextChar++;
		for(orbit=0; orbit < _numOrbits; orbit++)
		{
		    assert(isdigit(*nextChar));
		    numRead = sscanf(nextChar, "%lu", &count);
		    assert(numRead == 1);
		    ODV(lineNum,orbit) += count;
		    while(isdigit(*nextChar)) nextChar++; // read past current integer
		    assert(*nextChar == ' ' || (orbit == _numOrbits-1 && *nextChar == '\n'));
		    nextChar++;
		}
		assert(*nextChar == '\0');
		break;
	    case indexGraphlets: case indexGraphletsRNO: case indexOrbits: case indexMotifs: case indexMotifOrbits:
		fputs(line, stdout);
		if(_window)
		    while(fgets(line, sizeof(line), fpThreads[thread]))
			fputs(line, stdout);
		break;
	    case predict_merge: assert(false); break; // should not be here
	    case predict:
		Predict_ProcessLine(G, line);
		break;
	    default:
		Abort("oops... unknown or unsupported _outputMode in RunBlantInThreads while reading child process");
		break;
	    }
	}
	lineNum++;
    } while(threadsRunning > 0);

    // if numSamples is not a multiple of _THREADS, finish the leftover samples
    int leftovers = numSamples % _JOBS;
    return RunBlantFromGraph(_k, leftovers, G);
}

void BlantAddEdge(int v1, int v2)
{
    if(!_pairs) _pairs = Malloc(2*_maxEdges*sizeof(_pairs[0]));
    assert(_numEdges <= _maxEdges);
    if(_numEdges >= _maxEdges)
    {
	_maxEdges *=2;
	_pairs = Realloc(_pairs, 2*_maxEdges*sizeof(int));
    }
    _numNodes = MAX(_numNodes, v1+1); // add one since, for example, if we see a node numbered 100, numNodes is 101.
    _numNodes = MAX(_numNodes, v2+1);
    _pairs[2*_numEdges] = v1;
    _pairs[2*_numEdges+1] = v2;
    if(_pairs[2*_numEdges] == _pairs[2*_numEdges+1])
	Fatal("BlantAddEdge: edge %d (%d,%d) has equal nodes; cannot have self-loops\n", _numEdges, v1, v2);
    if(_pairs[2*_numEdges] > _pairs[2*_numEdges+1])
    {
	int tmp = _pairs[2*_numEdges];
	_pairs[2*_numEdges] = _pairs[2*_numEdges+1];
	_pairs[2*_numEdges+1] = tmp;
    }
    assert(_pairs[2*_numEdges] < _pairs[2*_numEdges+1]);
    _numEdges++;
}

int RunBlantEdgesFinished(int k, int numSamples, int numNodes, char **nodeNames)
{
    GRAPH *G = GraphFromEdgeList(_numNodes, _numEdges, _pairs, SPARSE);
    Free(_pairs);
    _nodeNames = nodeNames;
    return RunBlantInThreads(k, numSamples, G);
}

// Initialize the graph G from an edgelist; the user must allocate the pairs array
// to have 2*numEdges elements (all integers), and each entry must be between 0 and
// numNodes-1. The pairs array MUST be allocated using malloc or calloc, because
// we are going to free it right after creating G (ie., before returning to the caller.)
int RunBlantFromEdgeList(int k, int numSamples, int numNodes, int numEdges, int *pairs)
{
    assert(numNodes >= k);
    GRAPH *G = GraphFromEdgeList(numNodes, numEdges, pairs, SPARSE);
    Free(pairs);
    return RunBlantInThreads(k, numSamples, G);
}

const char * const USAGE_SHORT =
"BLANT (Basic Local Alignment of Network Topology): sample graphlets of up to 8 nodes from a graph.\n"\
"USAGE: blant [OPTIONS] -k numNodes -n numSamples graphInputFile\n"\
" Common options: (use -h for longer help)\n"\
"    -s samplingMethod (default MCMC; NBE, EBE, RES, AR, INDEX, EDGE_COVER)\n"\
"    -m{outputMode} (default o=ODV; g=GDV, f=frequency, i=index, r=root(used only for INDEX), d=distribution of neighbors\n"\
"    -d{displayModeForCanonicalIDs} (o=ORCA, j=Jesse, b=binaryAdjMatrix, d=decimal, i=integerOrdinal)\n"\
"    -r seed (integer)\n\n"\
"    -t N[:M]: (CURRENTLY BROKEN): use threading (parallelism); break the task up into N jobs (default 1) allowing\n"\
"        at most M to run at one time; M can be anything from 1 to a compile-time-specified maximum possible value\n"\
"        (MAX_POSSIBLE_THREADS in blant.h), but defaults to 4 to be conservative.";

const char * const USAGE_LONG =
"BLANT: Basic Local Alignment of Network Topology (work in progress)\n"\
"PURPOSE: sample graphlets of up to 8 nodes from a graph. Default output is similar to ORCA, though via stochastic sampling\n"\
"    rather than exaustive enumeration. Our APPROXIMATE results come MUCH faster than ORCA on large or dense networks.\n"\
"USAGE: blant [OPTIONS] -k numNodes -n numSamples graphInputFile\n"\
"where the following are REQUIRED:\n"\
"    numNodes is an integer 3 through 8 inclusive, specifying the size (in nodes) of graphlets to sample;\n"\
"    numSamples is the number of graphlet samples to take (large samples are recommended), except in INDEX sampling mode,\n"\
"	where it specifies the maximum number of samples to take from each node in the graph.\n"\
"	(note: -c {confidence} option is mutually exclusive to -n but is pending implementation)\n"\
"    samplingMethod is:\n"\
"	MCMC (default): Markov Chain Monte Carlo: Build the first set S of k nodes using NBE; then randomly remove and add\n"\
"         nodes to S using an MCMC graph walking algorithm with restarts; gives asymptotically correct relative frequencies\n"\
"         when using purely counting modes like -m{o|g|f}, but biased counts in indexing modes like -m{i|j} since we remove\n"\
"           duplicates in indexing modes.)\n"\
"	NBE (Node-Based Expansion): pick a node at random and add it to the node set S; add new nodes by choosing uniformly\n"\
"         at random from all nodes one step outside S. (Fast on sparse networks, slightly biased counts)\n"\
"	EBE (Edge-Based Expansion): pick an edge at random and add its two nodes to S; add nodes to S by picking an edge\n"\
"         uniformly at random from those emanating from S. (faster than NBE on dense networks, but more biased)\n"\
"	RES (Lu Bressan's REServoir sampling): also asymptotically correct but much slower than MCMC.\n"\
"	AR (Accept-Reject): EXTREMELY SLOW but asymptotically correct: pick k nodes entirely at random, reject if\n"\
"	  resulting graphlet is disconnected (vast majority of such grpahlets are disconnected, thus VERY SLOW)\n"\
"	INDEX: deterministic: for each node v in the graph, build a topologically deterministic set of k-graphlets to\n"\
"         be used as indices for seed-and-extend local alignments (using, eg., our onw Dijkstra-inspired local aligner--\n"\
"         see Dijkstra diretory). When using INDEX sampling, the -n command-line option specifies the maximum number\n"\
"         of index entries per starting node v.\n"\
"	EDGE_COVER: starting with E=all edges in the input graph, pick one edge in E and build a k-graphlet using EBE;\n"\
"         then subtract ALL its edges from E. Continue until E is empty. The goal is to output a short list of graphlets\n"\
"         that, in comglomerate, cover each edge in E at least once.\n"\
"    graphInputFile: graph must be in one of the following formats with its extension name:\n"\
"	Edgelist (.el), LEDA(.leda), GML (.gml), GraphML (.xml), LGF(.lgf), CSV(.csv)\n"\
"	(extensions .gz and .xz are automatically decompressed using gunzip and unxz, respectively)\n"\
"	Duplicate edges (either direction) and self-loops should be removed!\n"\
"COMMON OPTIONS:\n"\
"    -m{outputMode}, where {outputMode} is a single character, one of:\n"\
"	o = the default, which is ODV (Orbit Degree Vector), identical to ORCA (commonly though incorrectly called a GDV)\n"\
"	g = GDV (Graphlet Degree Vector) Note this is NOT what is commonly called a GDV, which is actually an ODV (above).\n"\
"	NOTE: the difference is that an ODV counts the number of nodes that touch all possible *orbits*, while a GDV lists\n"\
"		only the smaller vector of how many nodes touch each possible *graphlet* (independent of orbit).\n"\
"	f = graphlet {f}requency, similar to Relative Graphlet Frequency, produces a raw count across our random samples.\n"\
"	    sub-option -mf{freqDispMode} can be i(integer or count) or d(decimal or concentration)\n"\
"	i = {i}ndex: each line is a graphlet with columns: canonical ID, then k nodes in canonical order; useful since\n"\
"	    two lines with the same first column constitutes a PERFECT k-node local alignment between the two graphlets.\n"\
"	r = index with {r}oot node orbit: each line is a canonical ID + the orbit of the root node, then k nodes in canonical order; produces better seeds when the index is queried by the alignment algorithm\n"\
"	d = graphlet neighbor {D}istribution\n"\
"    -d{displayMode} [no default--MANDATORY for indexing modes]: single character controls how canonical IDs are displayed:\n"\
"	o = ORCA numbering\n"\
"	j = JESSE numbering\n"\
"	b = explicit binary representation of the half-adjacency matrix of the canonical graphlet\n"\
"	d = decimal (base-10) integer representation of the above binary\n"\
"	i = integer ordinal = sorting the above integers and numbering them 0, 1, 2, 3, etc.\n"\
"Less Common OPTIONS:\n"\
"    -t N[:M]: use threading (parallelism); break the task up into N jobs (default 1) allowing at most M to run at one time.\n"\
"       M can be anything from 1 to a compile-time-specified maximum possible value (MAX_POSSIBLE_THREADS in blant.h),\n"\
"       but defaults to 4 to be conservative.\n"\
"    -r seed: pick your own random seed\n"\
"    -w windowSize: DEPRECATED. (use '-h' option for more)\n"\
"	-p windowRepSamplingMethod: (DEPRECATED) one of the below\n"\
"	    MIN (Minimizer); MAX (Maximizer); DMIN (Minimizer With Distance); DMAX (Maximizer with Distance);\n"\
"	    LFMIN (Least Frequent Minimizer); LFMAX (Least Frequent Maximizer)\n"\
"	-P windowRepIterationMethods is one of: COMB (Combination) or DFS\n" \
"	-l windowRepLimitMethod is one of: [suffix N: limit to Top N satisfied graphlets]\n"\
"	    DEG (graphlet Total Degree); EDGE (1-step away numEdges)\n"\
"   -M = multiplicity, meaning max allowed number of ambiguous permutations in found graphlets (M=0 is a special case and means no max)\n" \
"   -T = top percent to expand to in -sINDEX sampling method (default 0)\n" \
"   -o = the orbit to use for the heuristic function\n" \
"   -f = the .orca4 file for the network\n" \
"   -a = sets whether ties are broken alphabetically or reverse alphabetically"\
;

// The main program, which handles multiple threads if requested.  We simply fire off a bunch of parallel
// blant *processes* (not threads, but full processes), and simply merge all their outputs together here
// in the parent.
int main(int argc, char *argv[])
{
    int i, j, opt, numSamples=0, multiplicity=1;
    confidence = 0;
    double windowRep_edge_density = 0.0;
    int exitStatus = 0;

    if(argc == 1)
    {
	puts(USAGE_SHORT);
	exit(1);
    }

    _JOBS = 1;
    _MAX_THREADS = 4;

    _k = 0; _k_small = 0;

    int odv_fname_len = 0;

    while((opt = getopt(argc, argv, "hm:d:t:r:s:c:k:K:o:f:e:g:w:p:P:l:n:M:T:a:")) != -1)
    {
	switch(opt)
	{
	case 'h':
	    printf("%s\n", USAGE_LONG);
	    printf("Note: current TSET size is %ld bits\n", 8*sizeof(TSET));
	    exit(1); break;
	case 'm':
	    if(_outputMode != undef) Fatal("tried to define output mode twice");
	    switch(*optarg)
	    {
	    case 'm': _outputMode = indexMotifs; break;
	    case 'M': _outputMode = indexMotifOrbits; break;
	    case 'i': _outputMode = indexGraphlets; break;
	    case 'r': _outputMode = indexGraphletsRNO; break;
	    case 'j': _outputMode = indexOrbits; break;
	    case 'f': _outputMode = graphletFrequency;
		switch (*(optarg + 1))
		{
		    case 'i': _freqDisplayMode = count; break;
		    case 'd': _freqDisplayMode = concentration; break;
		    case '\0': _freqDisplayMode = freq_display_mode_undef; break;
		    default: Fatal("-mf%c: unknown frequency display mode;\n"
		    "\tmodes are i=integer(count), d=decimal(concentration)", *(optarg + 1));
		    break;
		}
	    break;
	    case 'g': _outputMode = outputGDV; break;
	    case 'o': _outputMode = outputODV; break;
	    case 'd': _outputMode = graphletDistribution; break;
	    case 'p': _outputMode = predict; break;
	    case 'q': _outputMode = predict_merge; break;
	    default: Fatal("-m%c: unknown output mode \"%c\"", *optarg,*optarg);
	    break;
	    }
	    break;
	case 'd':
	    if (_displayMode != undefined) Fatal("tried to define canonical display mode twice");
	    switch(*optarg)
	    {
	    case 'b': _displayMode = binary; break;
	    case 'd': _displayMode = decimal; break;
	    case 'i': _displayMode = ordinal; break;
	    case 'j': _displayMode = jesse; break;
	    case 'o': _displayMode = orca; break;
	    default: Fatal("-d%c: unknown canonical display mode:n"
		    "\tmodes are i=integer ordinal, d=decimal, b=binary, o=orca, j=jesse", *optarg);
	    break;
	    }
	    break;
	case 't': assert(1==sscanf(optarg, "%d", &_JOBS));
	    _MAX_THREADS = _JOBS;
	    assert(1 <= _JOBS && _MAX_THREADS <= MAX_POSSIBLE_THREADS);
	    break;
	case 'r': _seed = atoi(optarg); if(_seed==-1)Apology("seed -1 ('-r -1' is reserved to mean 'uninitialized'");
	    break;
	case 's':
	    if (_sampleMethod != -1) Fatal("Tried to define sampling method twice");
	    else if (strncmp(optarg, "NBE", 3) == 0)
		_sampleMethod = SAMPLE_NODE_EXPANSION;
	    else if (strncmp(optarg, "FAYE", 4) == 0)
		_sampleMethod = SAMPLE_FAYE;
	    else if (strncmp(optarg, "EBE", 3) == 0)
		_sampleMethod = SAMPLE_EDGE_EXPANSION;
	    else if (strncmp(optarg, "MCMC",4) == 0) {
		_sampleMethod = SAMPLE_MCMC;
		if (strchr(optarg, 'u') || strchr(optarg, 'U'))
		    _MCMC_EVERY_EDGE=true;
	    }
	    else if (strncmp(optarg, "EDGE_COVER", 10) == 0) {
		_sampleMethod = SAMPLE_MCMC;
		_sampleSubmethod = SAMPLE_MCMC_EC;
	    }
	    else if (strncmp(optarg, "RES", 3) == 0)
		_sampleMethod = SAMPLE_RESERVOIR;
	    else if (strncmp(optarg, "AR", 2) == 0)
		_sampleMethod = SAMPLE_ACCEPT_REJECT;
	    else if (strncmp(optarg, "INDEX", 5) == 0)
		_sampleMethod = SAMPLE_INDEX;
	    else
	    {
		_sampleFileName = optarg;
		if(strcmp(optarg,"STDIN") == 0) _sampleFile = stdin;
		else _sampleFile = fopen(_sampleFileName, "r");
		if(!_sampleFile)
		    Fatal("Unrecognized sampling method specified: '%s'. Options are: {NBE|EBE|MCMC|RES|FAYE|AR|{filename}}\n"
			"If unrecognized, we try opening a file by the name '%s', but no such file exists",
			_sampleFileName, _sampleFileName);
		_sampleMethod = SAMPLE_FROM_FILE;
	    }
	    break;
	case 'c': confidence = atof(optarg);
        if (confidence < 0 || confidence > 1) Fatal("Confidence level must be between 0 and 1");
	    // Apology("confidence intervals not implemented yet");
	    break;
	case 'k': _k = atoi(optarg);
		if (_GRAPH_GEN && _k >= 33) {
			_k_small = _k % 10;
			if (!(3 <= _k_small && _k_small <= 8))
			    Fatal("%s\nERROR: k [%d] must be between 3 and 8\n%s", USAGE_SHORT, _k_small);
			_k /= 10;
			assert(_k_small <= _k);
		} // First k indicates stamping size, second k indicates KS test size.
	    if (!(3 <= _k && _k <= 8)) Fatal("%s\nERROR: k [%d] must be between 3 and 8\n%s", USAGE_SHORT, _k);
	    break;
	case 'w': _window = true; _windowSize = atoi(optarg); break;
	case 'p':
	    if (_windowSampleMethod != -1) Fatal("Tried to define window sampling method twice");
	    else if (strncmp(optarg, "DMIN", 4) == 0)
		_windowSampleMethod = WINDOW_SAMPLE_MIN_D;
	    else if (strncmp(optarg, "DMAX", 4) == 0)
		_windowSampleMethod = WINDOW_SAMPLE_MAX_D;
	    else if (strncmp(optarg, "MIN", 3) == 0)
		_windowSampleMethod = WINDOW_SAMPLE_MIN;
	    else if (strncmp(optarg, "MAX", 3) == 0)
		_windowSampleMethod = WINDOW_SAMPLE_MAX;
	    else if (strncmp(optarg, "LFMIN", 5) == 0)
		_windowSampleMethod = WINDOW_SAMPLE_LEAST_FREQ_MIN;
	    else if (strncmp(optarg, "LFMAX", 5) == 0)
		_windowSampleMethod = WINDOW_SAMPLE_LEAST_FREQ_MAX;
	    else if (strncmp(optarg, "DEGMAX", 6) == 0)
		_windowSampleMethod = WINDOW_SAMPLE_DEG_MAX;
	    else
		Fatal("Unrecognized window searching method specified. Options are: -p[u|U]{MIN|MAX|DMIN|DMAX|LFMIN|LFMAX|DEGMAX}\n");
	    break;
	case 'P':
		if (strncmp(optarg, "COMB", 4) == 0)
			_windowIterationMethod = WINDOW_ITER_COMB;
		else if (strncmp(optarg, "DFS", 3) == 0)
			_windowIterationMethod = WINDOW_ITER_DFS;
		else
			Fatal("Unrecognized window Iteration method specified. Options are: -P{COMB|DFS}\n");
		break;
	case 'l':
		if (_windowRep_limit_method != WINDOW_LIMIT_UNDEF) Fatal("Tried to define window limiting method twice");
		if (strncmp(optarg, "n", 1) == 0 || strncmp(optarg, "N", 1) == 0) {
		    _windowRep_limit_neglect_trivial = true; optarg += 1;
		}
		if (strncmp(optarg, "DEG", 3) == 0) {
			_windowRep_limit_method = WINDOW_LIMIT_DEGREE; optarg += 3;
		}
		else if (strncmp(optarg, "EDGE", 4) == 0) {
			_windowRep_limit_method = WINDOW_LIMIT_EDGES; optarg += 4;
		}
		else
			Fatal("Unrecognized window limiting method specified. Options are: -l{DEG}{EDGE}{limit_num}\n");
		_numWindowRepLimit = atoi(optarg);
		if (!_numWindowRepLimit) {_numWindowRepLimit = 10; _numWindowRepArrSize = _numWindowRepLimit;}
		_windowRep_limit_heap = HeapAlloc(_numWindowRepLimit, asccompFunc, NULL);
		break;
	case 'n': numSamples = atoi(optarg);
		char lastChar = optarg[strlen(optarg)-1];
		if(!isdigit(lastChar))
		    switch(lastChar) {
		    case 'b': case 'B': case 'g': case 'G': numSamples *= 1024; // do NOT break, fall through
		    case 'm': case 'M': numSamples *= 1024;
		    case 'k': case 'K': numSamples *= 1024; break;
		    default: Fatal("%s\nERROR: numSamples can be appended by k, m, b, or g but not %c\n%s", USAGE_SHORT, lastChar);
		    break;
		    }
		if(numSamples < 0) Fatal("%s\nFatal Error: numSamples [%d] must be non-negative", USAGE_SHORT, numSamples);
		//fprintf(stderr, "numSamples set to %d\n", numSamples);
		break;
	case 'K': _KS_NUMSAMPLES = atoi(optarg);
	    break;
	case 'e':
	    _GRAPH_GEN_EDGES = atoi(optarg);
	    windowRep_edge_density = atof(optarg);
	    break;
	case 'g':
	    if (!GEN_SYN_GRAPH) Fatal("Turn on Global Variable GEN_SYN_GRAPH");
	    _GRAPH_GEN = true;
	    if (_genGraphMethod != -1) Fatal("Tried to define synthetic graph generating method twice");
	    else if (strncmp(optarg, "NBE", 3) == 0) _genGraphMethod = GEN_NODE_EXPANSION;
	    else if (strncmp(optarg, "MCMC", 4) == 0) Apology("MCMC for Graph Syn is not ready");  // _genGraphMethod = GEN_MCMC;
	    else Fatal("Unrecognized synthetic graph generating method specified. Options are: -g{NBE|MCMC}\n");
	    break;
	case 'M': multiplicity = atoi(optarg);
	    if(multiplicity < 0) Fatal("%s\nERROR: multiplicity [%d] must be non-negative\n", USAGE_SHORT, multiplicity);
    case 'T': _topThousandth = atoi(optarg);
	    break;
	case 'o':
        _orbitNumber = atoi(optarg);
        break;
    case 'f':
        odv_fname_len = strlen(optarg);
        _odvFile = malloc(sizeof(char) * odv_fname_len);
        strncpy(_odvFile, optarg, odv_fname_len);
        break;
    case 'a':
        _alphabeticTieBreaking = atoi(optarg) != 0;
        break;
	default: Fatal("%s\nERROR: unknown option %c", USAGE_SHORT, opt);
    }
    }

    if (_sampleMethod == -1) _sampleMethod = SAMPLE_MCMC; // the default

    if (_orbitNumber != -1) {
        if (_odvFile != NULL) {
            parseOdvFromFile(_odvFile);
        } else {
            Fatal("an ODV orbit number was provided, but no ODV file path was supplied");
        }
    }

    if (_sampleMethod == SAMPLE_INDEX && _k <= 5) Fatal("k is %d but must be at least 6 for INDEX sampling method because there are no unambiguous graphlets for k<=5",_k);

    if(_seed == -1) _seed = GetFancySeed(false);
    // This only seeds the main thread; sub-threads, if they exist, are seeded later by "stealing"
    // exactly _THREADS-1 values from near the beginning of this main random stream.
    RandomSeed(_seed);

    if(_outputMode == undef) _outputMode = outputODV; // default to the same thing ORCA and Jesse us
	if (_freqDisplayMode == freq_display_mode_undef) // Default to integer(count)
		_freqDisplayMode = count;

    if(numSamples!=0 && confidence>0 && !_GRAPH_GEN)
	Fatal("cannot specify both -n (sample size) and -c (confidence)");

    FILE *fpGraph;
    int piped = 0;
    if(!argv[optind])
    {
	fpGraph = stdin;
	if(isatty(0)) Warning("reading graph input file from terminal, press ^D to finish");
    }
    else {
	char *graphFileName = argv[optind];
	fpGraph = readFile(graphFileName, &piped);
	if(!fpGraph) Fatal("cannot open graph input file '%s'\n", argv[optind]);
	optind++;
    }
    assert(optind == argc || _GRAPH_GEN || _windowSampleMethod == WINDOW_SAMPLE_DEG_MAX);

    SetBlantDir(); // Needs to be done before reading any files in BLANT directory
    SetGlobalCanonMaps(); // needs _k to be set
    LoadMagicTable(); // needs _k to be set

    if (_window && _windowSize >= 3) {
        if (_windowSampleMethod == -1) Fatal("Haven't specified window searching method. Options are: -p{MIN|MAX|DMIN|DMAX|LFMIN|LFMAX}\n");
        if(_windowSize < _k) Fatal("windowSize must be at least size k\n");
        _MAXnumWindowRep = CombinChooseDouble(_windowSize, _k);
        _numWindowRepArrSize = _MAXnumWindowRep > 0 ? MIN(_numWindowRepArrSize, _MAXnumWindowRep) : _numWindowRepArrSize;
        _windowReps = Calloc(_numWindowRepArrSize, sizeof(int*));
        for(i=0; i<_numWindowRepArrSize; i++) _windowReps[i] = Calloc(_k+1, sizeof(int));
        if (windowRep_edge_density < 0) windowRep_edge_density = 0;
		if (windowRep_edge_density > 1) windowRep_edge_density = 1;
		_windowRep_min_num_edge = (int) CombinChooseDouble(_k, 2) * windowRep_edge_density;
		if (_windowRep_min_num_edge < 0) Fatal("WindowRep minimum number of edges must be larger than 0. Check edge density\n");
    }

    // This section of code computes the info necessary to implement "multiplicity" mode ('M' above), which controls
    // whether to output a sampled graphlet at all during INDEXING based on how much "ambiguity" there is in its
    // local alignment. The more permutations of the graphlet there are, the less useful it is for seeding local
    // alignments and therefore the less useful as a database index entry.
    _windowRep_allowed_ambig_set = SetAlloc(_numCanon);
    SET *orbit_temp = SetAlloc(_numOrbits);
    for(i=0; i<_numCanon; i++) {
        if (SetIn(_connectedCanonicals, i)) {
            // calculate number of permutations for the given canonical graphlet
	    // (loop through all unique orbits and count how many times they appear)
            // the formula is for every unique orbit, multiply the number of permutations by
	    // the factorial of how many appearances that unique orbit has
            // if there is one orbit with three nodes and a second orbit with 2 nodes,
	    // the number of permutations would be (3!)(2!)
            // NOTE: this is not the most efficient algorithm since it doesn't use hash tables.
	    // I didn't want to overcomplicate it because it only happens once per run.
            // however, if speed is important (this currently takes about 5 seconds on k=8) this can be sped up
            for(j=0; j<_k; j++) SetAdd(orbit_temp, _orbitList[i][j]);
            unsigned uniq_orbits[_k];
            unsigned num_uniq_orbits = SetToArray(uniq_orbits, orbit_temp);
            unsigned uniq_orbit_i;
            unsigned total_orbit_perms = 1;

            for (uniq_orbit_i=0; uniq_orbit_i<num_uniq_orbits; uniq_orbit_i++) {
                unsigned orbit_appearances = 0;
                unsigned orbit_i = 0;

                for (orbit_i=0; orbit_i<_k; orbit_i++) {
                    if (_orbitList[i][orbit_i] == uniq_orbits[uniq_orbit_i]) {
                        orbit_appearances++;
                        total_orbit_perms *= orbit_appearances;
                    }
                }
            }

            // I know it's inefficient to put multiplicity here instead of around the whole orbit perm calculation code but
	    // it increases readability, at least until orbit perm calculation is put into a function
            if(multiplicity == 0 || total_orbit_perms <= multiplicity) { // multiplicity = 0 means any ambiguity is allowed
                SetAdd(_windowRep_allowed_ambig_set, i);
            }
            SetEmpty(orbit_temp);
        }
    }
    SetFree(orbit_temp);

    // Read network using native Graph routine.
    GRAPH *G = GraphReadEdgeList(fpGraph, SPARSE, _supportNodeNames);
    if(_supportNodeNames)
    {
	assert(G->name);
	_nodeNames = G->name;
    }
    if(fpGraph != stdin) closeFile(fpGraph, &piped);

    if (_windowSampleMethod == WINDOW_SAMPLE_DEG_MAX)
    {
        FILE *fp;
        _graphNodeImportance = Calloc(G->n, sizeof(float));
        if((optind + 1) == argc) {
            _supportNodeImportance = true;
            fp = fopen(argv[optind++], "r");
            if (fp == NULL) Fatal("cannot open graph Node Importance File.");
            char line[BUFSIZ], nodeName[BUFSIZ];
            foint nodeNum;
            float importance;

            while(fgets(line, sizeof(line), fp))
            {
                if(sscanf(line, "%s%f ", nodeName, &importance) != 2)
                    Fatal("GraphNodeImportance: Error while reading\n");
                if(!BinTreeLookup(G->nameDict, (foint)nodeName, &nodeNum))
                    Fatal("Node Importance Error: %s is not in the Graph file\n", nodeName);
                _graphNodeImportance[nodeNum.i] = importance;
            }
        }
    }

#if GEN_SYN_GRAPH
    FILE *fpSynGraph = NULL;
    if (_GRAPH_GEN) {
    	if (_k_small == 0) _k_small = _k;
        if (numSamples == 0) Fatal("Haven't specified sample size (-n sampled_size)");
        if (confidence == 0) confidence = 0.05;
        if (_KS_NUMSAMPLES == 0) _KS_NUMSAMPLES = 1000;
        if (_GRAPH_GEN_EDGES == 0) _GRAPH_GEN_EDGES = G->numEdges;
        if((optind + 1) == argc) {
            fpSynGraph = fopen(argv[optind++], "w");
            if (fpSynGraph == NULL) Fatal("cannot open synthetic graph outputfile.");
        }
        exitStatus = GenSynGraph(_k, _k_small, numSamples, G, fpSynGraph);
    }
#endif

    if(_outputMode == predict_merge)
	exitStatus = Predict_Merge(G);
    else
	exitStatus = RunBlantInThreads(_k, numSamples, G);
    GraphFree(G);
    return exitStatus;
}
