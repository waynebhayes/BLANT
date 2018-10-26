#include <sys/file.h>
#include <unistd.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "misc.h"
#include "tinygraph.h"
#include "graph.h"
#include "blant.h"
#include "sets.h"

char USAGE[] = "synthetic -k maxk -s STAGNATED Gtarget.el Gsynth.el [--blant.Gt.index files--] [--blant.Gs.index files--]\n - output is new Gsynth.el\n";

Boolean _supportNodeNames = false;

// The following is the most compact way to store the permutation between a non-canonical and its canonical representative,
// when k=8: there are 8 entries, and each entry is a integer from 0 to 7, which requires 3 bits. 8*3=24 bits total.
// For simplicity we use the same 3 bits per entry, and assume 8 entries, even for k<8.  It wastes memory for k<4, but
// makes the coding much simpler.
typedef unsigned char kperm[3]; // The 24 bits are stored in 3 unsigned chars.


static int _k[maxK]; // stores what values of k have to be considered e.g. [3,4,5,6] or [3,4,5] or [2,7,8]. Will be followed by -1s
static int _numCanon[maxK];  // canonicals for particular value of k. So for k=5, _numCanon[5-1] stores ~32
static SET *_connectedCanonicals[maxK];
static int _maxNumCanon = -1;  // max number of canonicals
static int _numSamples = -1;  // same number of samples in each blant index file
static int _canonList[maxK][MAX_CANONICALS];
static int _stagnated = 1000;

// Here's where we're lazy on saving memory, and we could do better.  We're going to allocate a static array
// that is big enough for the 256 million permutations from non-canonicals to canonicals for k=8, even if k<8.
// So we're allocating 256MBx3=768MB even if we need much less.  I figure anything less than 1GB isn't a big deal
// these days. It needs to be aligned to a page boundary since we're going to mmap the binary file into this array.
static kperm Permutations[maxBk] __attribute__ ((aligned (8192)));
// Here's the actual mapping from non-canonical to canonical, same argument as above wasting memory, and also mmap'd.
// So here we are allocating 256MB x sizeof(short int) = 512MB.
// Grand total statically allocated memory is exactly 1.25GB.
static short int* _K[maxK];


// Assuming the global variable _k[] is set properly, go read in and/or mmap the big global
// arrays related to canonical mappings and permutations.
void SetGlobalCanonMaps(void){
	unsigned int _Bk;
	int i;
	for(i=0; i<maxK; i++){  // for all values of 'k'
		if (_k[i] == -1) 
			break;
		assert(3 <= _k[i] && _k[i] <= 8);
		_Bk = (1 <<(_k[i]*(_k[i]-1)/2));
		if ((_Bk * sizeof(short int)) < 8192)
			_Bk = 8192 / sizeof(short int);  // aligned_alloc constraint
		char BUF[BUFSIZ];
		_connectedCanonicals[_k[i]-1] = SetAlloc(_Bk);
		_numCanon[_k[i]-1] = canonListPopulate(BUF, _canonList[_k[i]-1], _connectedCanonicals[_k[i]-1], _k[i]);
		// Pushkar: need to add connectedCanonicals in above line
		if (_numCanon[_k[i]-1] > _maxNumCanon)  // set max number of canonicals for a k
			_maxNumCanon = _numCanon[_k[i]-1];
		_K[_k[i]-1] = (short int*) aligned_alloc(8192, _Bk * sizeof(short int));
		assert(_K[_k[i]-1] != NULL);
		mapCanonMap(BUF, _K[_k[i]-1], _k[i]);
		sprintf(BUF, CANON_DIR "/perm_map%d.bin", _k[i]);
		int pfd = open(BUF, 0*O_RDONLY);
		_Bk = (1 <<(_k[i]*(_k[i]-1)/2));  // reset
		kperm *Pf = Mmap(Permutations, _Bk*sizeof(Permutations[0]), pfd);
    	assert(Pf == Permutations);
	}
}

// Given the big graph G and a set of nodes in V, return the TINY_GRAPH created from the induced subgraph of V on G.
static TINY_GRAPH *TinyGraphInducedFromGraph(TINY_GRAPH *Gv, GRAPH *G, int *Varray){
    unsigned i, j;
    TinyGraphEdgesAllDelete(Gv);
    for(i=0; i < Gv->n; i++) for(j=i+1; j < Gv->n; j++)
        if(GraphAreConnected(G, Varray[i], Varray[j]))
            TinyGraphConnect(Gv, i, j);
    return Gv;
}

double FastObjective(double oldcost, double olddelta, double change){
    // olddelta is the difference in graphlet counts of D[0] and D[1] before the surgery
    // change is the increment/decrement in count of a graphlet, which will be +1/-1
    double unchanged = SQR(oldcost) - SQR(olddelta);
    double newcost_sq = unchanged + SQR(olddelta+change);
    return sqrt(newcost_sq);
}

double ReBLANT(int _maxNumCanon, int D[2][maxK][_maxNumCanon], GRAPH *G, SET ***samples, int ***Varrays, int ***BLANT, int v1, int v2, double oldcost){
    int i, j, line, s;
    static TINY_GRAPH *g[maxK];

    double olddelta, change;
    double newcost = oldcost;

    for (i=0; i<maxK; i++){
    	if (_k[i] == -1)
    		break;
    	
    	// allocate a tiny graph 
    	if (!g[_k[i]-1]) 
    		g[_k[i]-1] = TinyGraphAlloc(_k[i]);

        /*
    	{
    	int testCount = 0;
    	for(j=0; j<_numCanon[_k[i]-1]; j++)
    		testCount += D[1][_k[i]-1][j];
    	assert(testCount == _numSamples);
    	}*/

    	for (s=1; line=Varrays[_k[i]-1][v1][s], s<=Varrays[_k[i]-1][v1][0]; s++)    		
    		if(SetIn(samples[_k[i]-1][v2], line)){
    			
                // decrement a graphlet
                olddelta = D[1][_k[i]-1][BLANT[_k[i]-1][line][0]] - D[0][_k[i]-1][BLANT[_k[i]-1][line][0]];
                --D[1][_k[i]-1][BLANT[_k[i]-1][line][0]];
                change = -1;
                newcost = FastObjective(newcost, olddelta, change);

    		TinyGraphInducedFromGraph(g[_k[i]-1], G, BLANT[_k[i]-1][line]+1);
    		BLANT[_k[i]-1][line][0] = _K[_k[i]-1][TinyGraph2Int(g[_k[i]-1], _k[i])];
    			
                // increment a graphlet
                olddelta = D[1][_k[i]-1][BLANT[_k[i]-1][line][0]] - D[0][_k[i]-1][BLANT[_k[i]-1][line][0]];
                ++D[1][_k[i]-1][BLANT[_k[i]-1][line][0]];
                change = 1;
                newcost = FastObjective(newcost, olddelta, change);
    		}
    	    	
    	{
    	int testCount = 0;
    	for(j=0; j<_numCanon[_k[i]-1]; j++)
    		testCount += D[1][_k[i]-1][j];
    	assert(testCount == _numSamples);
    	}
    }

    return newcost;
}

double PoissonDistribution(double l, int k){
    //->   (e^(-l) x l^k) / k!
    double r = exp(-l);
    int i;
    for(i=k; i>0; i--) // divide by k!
	r *= l/i;
    return r;
}


double Objective(int _maxNumCanon, int D[2][maxK][_maxNumCanon]){
    int i,j;
    double logP = 0, sum2 = 0;

    for (i=0; i<maxK; i++){
    	if (_k[i] == -1) break;
    	for (j=0; j<_numCanon[_k[i]-1]; j++){
	    double pd = PoissonDistribution(D[0][_k[i]-1][j], D[1][_k[i]-1][j]);
	    if(pd > 1)
		Fatal("umm.... PoissonDistribution returned a number greater than 1");
	    if(pd>0) 
		logP += log(pd); // if we're close, use probability
	    // use this one when we're so far away the probability is zero
	    double term = SQR((double)D[0][_k[i]-1][j] - D[1][_k[i]-1][j]);
	    sum2 += term;
    	}
    }
    double returnVal = sqrt(sum2);
    assert(returnVal == returnVal);
    return returnVal; //exp(logP);
}

int main(int argc, char *argv[])
{
    srand48(time(0)+getpid());
    int i, opt, j, line;

    if(argc == 1)
    {
	fprintf(stderr, "%s\n", USAGE);
	exit(1);
    }

    // initialize _k[]
    for(i=0; i<maxK; i++)
    	_k[i] = -1; 

    /*
	Read max k and stagnation
    */
    int kmax;
    while((opt = getopt(argc, argv, "k:s:")) != -1){
		switch(opt){
			case 'k': 
				kmax = atoi(optarg);
			    if(!(3 <= kmax && kmax <= 8)) Fatal("k must be between 3 and 8\n%s", USAGE);

			    // _k[] = [3,4,5,6,k]
			    for(i=3; i<=kmax; i++)
			    	_k[i-3] = i;
			    break;
			case 's': _stagnated = atoi(optarg);
			    if(!(_stagnated>=10)) Fatal("STAGNATED must be > 10\n%s", USAGE);
			    break;
			default: Fatal("unknown option %c\n%s", opt, USAGE);
		}
    }

    SetGlobalCanonMaps(); // needs _k[] to be set

    GRAPH *G[2]; // G[0] is the target, G[1] is the synthetic
    for(i=0; i<2; i++){
		if(!argv[optind]) Fatal("no input graph file specified\n%s", USAGE);
		FILE *fpGraph = fopen(argv[optind], "r");
		if(!fpGraph) Fatal("cannot open graph input file '%s'\n", argv[optind]);
		// Read it in using native Graph routine.
		G[i] = GraphReadEdgeList(fpGraph, true, _supportNodeNames); // sparse=true
		if(_supportNodeNames)
		    assert(G[i]->name);
		fclose(fpGraph);
		optind++;
    }

    // The distribution of graphlets (squiggly plot vectors)
    assert(_maxNumCanon != -1);  // this should be set >0 by calling SetGlobalCanonMaps() first
    int D[2][maxK][_maxNumCanon];
    for(i=0; i<maxK;i++)
    	for (j=0; j<_maxNumCanon; j++) 
    		D[0][i][j] = D[1][i][j] = 0;

    // expect 2 blant files (target & synthetic for every _k value)
    // assume all blant files have same number of samples = _numSamples
    int **BLANT[2][maxK];
    for(i=0;i<2;i++){
    	for(j=0; j<maxK; j++){
    		if(_k[j] == -1)
    			break;

    		char cmd[BUFSIZ];
    		sprintf(cmd, "wc -l < %s", argv[optind]);
    		FILE *fp = popen(cmd, "r");
    		assert(fp);
    		if (_numSamples == -1){
    			fscanf(fp, "%d", &_numSamples);
    			assert(_numSamples > 0);
    		}
    		else{
    			int tempsamples;
    			fscanf(fp, "%d", &tempsamples);
    			assert(tempsamples == _numSamples);
    		}
    		pclose(fp);
			fp=fopen(argv[optind], "r");
			assert(fp);

			BLANT[i][_k[j]-1] = (int**) Malloc(_numSamples * sizeof(int*));
			for (line=0; line<_numSamples; line++){
				BLANT[i][_k[j]-1][line] = (int*) Malloc((maxK+1) * sizeof(int));
				int l;
				for (l=0; l<=_k[j]; l++)
					assert(1 == fscanf(fp, "%d", &(BLANT[i][_k[j]-1][line][l])));
                assert(BLANT[i][_k[j]-1][line][0] < _maxNumCanon);
				++D[i][_k[j]-1][BLANT[i][_k[j]-1][line][0]]; // squiggly plot
			}
			assert(fscanf(fp, "%d", &line) < 1); // ensure there's nothing left to read.
			fclose(fp);
			optind++;
    	}
    
    }

    // check (squiggly vectors count == numsamples)
    for(i=0; i<2; i++){
    	for (j=0; j<maxK; j++){
    		if (_k[j] == -1)
    			break;
    		int testCount = 0;
    		int l;
    		for (l=0; l<_numCanon[_k[j]-1]; l++)
    			testCount += D[i][_k[j]-1][l];
    		assert(testCount == _numSamples);
    	}
    }

    // store what samples is a particular node-part of
    SET **samples[maxK];
    for(i=0; i<maxK; i++){
    	if (_k[i] == -1) break;
    	samples[_k[i]-1] = (SET**) Malloc(G[1]->n * sizeof(SET*));
    }

    for(i=0; i<maxK; i++){
    	if (_k[i] == -1) break;

    	for (j=0; j<G[1]->n; j++)
    		samples[_k[i]-1][j] = SetAlloc(_numSamples);

    	for (line=0; line<_numSamples; line++)
    		for (j=1; j<= _k[i]; j++)
    			SetAdd(samples[_k[i]-1][BLANT[1][_k[i]-1][line][j]], line);	   // SetAdd(samples[k][nodenum], line);
    }

    int **Varrays[maxK];
    for (i=0; i<maxK; i++){
    	if (_k[i] == -1) break;
    	Varrays[_k[i]-1] = (int**) Malloc(G[1]->n * sizeof(int*));
    	
    	for (j=0; j < G[1]->n; j++){
    		Varrays[_k[i]-1][j] =  (int*) Malloc((1+SetCardinality(samples[_k[i]-1][j]))* sizeof(int));
    		Varrays[_k[i]-1][j][0] = SetToArray(Varrays[_k[i]-1][j]+1, samples[_k[i]-1][j]);
    	}
    }


    // while(not done---either some number of iterations, or objective function says we're too far away)
    double cost = Objective(_maxNumCanon, D), startCost = cost, newCost, maxCost = cost;  // evaluate Objective() once, at the start. 
    assert(cost == cost);
    long int sa_iter = 0;
    double temperature, pBad, unif_random;
    
    while(fabs(cost - startCost)/startCost < 0.5) // that's enough progress otherwise we're over-optimizing at this sample size
    {
	int edge = drand48() * G[1]->numEdges;
	int u1, u2, v1 = G[1]->edgeList[2*edge], v2 = G[1]->edgeList[2*edge+1];
	do {
	    u1 = drand48()*G[1]->n;
	    do u2 = drand48()*G[1]->n; while(u1==u2);
	} while(GraphAreConnected(G[1], u1, u2)); // find a non-edge  u1,u2
	assert(GraphAreConnected(G[1], v1, v2));  // find edge v1,v2

	GraphDisconnect(G[1], v1, v2); // remove edge e from Gs
	newCost = ReBLANT(_maxNumCanon, D, G[1], samples, Varrays, BLANT[1], v1, v2, cost);
        //printf("compare fast-objective & objective, cost=%g, OCOST=%g\n", newCost, Objective(_maxNumCanon, D));
	
	GraphConnect(G[1], u1, u2); // add an edge to Gs
	newCost = ReBLANT(_maxNumCanon, D, G[1], samples, Varrays, BLANT[1], u1, u2, newCost);
        //printf("compare fast-objective & objective, cost=%g, OCOST=%g\n", newCost, Objective(_maxNumCanon, D));

	//newCost = Objective(_maxNumCanon, D);
	maxCost = MAX(maxCost, newCost);
	assert(newCost == newCost);
	static int same;
#if 1 // HILLCLIMB
	if(newCost < cost)
#else
#define SA_START_TEMP 1.0
#define SA_DECAY 1e-4
	temperature = SA_START_TEMP * exp(-SA_DECAY*sa_iter);
	unif_random = drand48();
	pBad = exp((-(newCost - cost)/maxCost)/temperature);
	// assert (newCost < cost <=> pBad > 1)
	assert((newCost < cost) == (pBad > 1));
	if(newCost < cost || unif_random < pBad)
#endif
	{
	    fprintf(stderr,"A");fflush(stderr);
	    static double printVal, printInterval;
	    if(//fabs(newCost - printVal)/printVal >= 0.02 ||
		++printInterval > 100)
	    {
		fprintf(stderr, "\ntemp %g cost %g newCost %g maxCost %g pBad %g", temperature, cost, newCost, maxCost, pBad);
		//fprintf(stderr, "%g ", newCost);
		printVal = newCost;
		printInterval = 0;
	    }
	    cost = newCost;
	    same = 0;
	}
	else // revert
	{
	    fprintf(stderr,"R");fflush(stderr);
	    static double printVal, printInterval;
	    if(//fabs(newCost - printVal)/printVal >= 0.02 ||
		++printInterval > 100)
	    {
		fprintf(stderr, "\ntemp %g cost %g newCost %g maxCost %g pBad %g", temperature, cost, newCost, maxCost, pBad);
		//fprintf(stderr, "%g ", newCost);
		printVal = newCost;
		printInterval = 0;
	    }
	    ++same;
	    GraphDisconnect(G[1], u1, u2);
	    ReBLANT(_maxNumCanon, D, G[1], samples, Varrays, BLANT[1], u1, u2, newCost);  // ignore the returned newcost
	    GraphConnect(G[1], v1, v2);
	    ReBLANT(_maxNumCanon, D, G[1], samples, Varrays, BLANT[1], v1, v2, newCost);  // ignore the returned newcost
	}
	if(same > _stagnated) break;
	++sa_iter;
    }
    fprintf(stderr,"\n");
    for(i=0; i < G[1]->numEdges; i++)
	printf("%d %d\n", G[1]->edgeList[2*i], G[1]->edgeList[2*i+1]);
}
