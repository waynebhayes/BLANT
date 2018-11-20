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
#include "syntheticDS.h" // some data structures (Revert stacks, and GDV hash tables)

char USAGE[] = "synthetic -k maxk -s STAGNATED Gtarget.el Gsynth.el [--blant.Gt.index files--] [--blant.Gs.index files--]\n - output is new Gsynth.el\n";

Boolean _supportNodeNames = false;

// The following is the most compact way to store the permutation between a non-canonical and its canonical representative,
// when k=8: there are 8 entries, and each entry is a integer from 0 to 7, which requires 3 bits. 8*3=24 bits total.
// For simplicity we use the same 3 bits per entry, and assume 8 entries, even for k<8.  It wastes memory for k<4, but
// makes the coding much simpler.
typedef unsigned char kperm[3]; // The 24 bits are stored in 3 unsigned chars.


static int _k[maxK]; // stores what values of k have to be considered e.g. [3,4,5,6] or [3,4,5] or [2,7,8]. Will be followed by -1s
static int _numCanon[maxK];  // canonicals for particular value of k. So for k=5, _numCanon[5-1] stores ~32
static int totalCanons;
static SET *_connectedCanonicals[maxK];
static int _maxNumCanon = -1;  // max number of canonicals
static int _numSamples = -1;  // same number of samples in each blant index file
static int _maxNodes = -1;
static int _canonList[maxK][MAX_CANONICALS];
static int _stagnated = 1000, _numDisconnectedGraphlets;

#define IGNORE_DISCONNECTED_GRAPHLETS 1  // Note: this param is implemented only for the GDV objective as of now
#define USING_GDV_OBJECTIVE 1
#define PRINT_INTERVAL 10000
#define NUMPROPS 2 // degree-distribution is 0th, graphlets is 1st

// NORMALIZATION
static double weights[NUMPROPS] = {0, 1};
static double max_abscosts[NUMPROPS];

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
        char BUF[BUFSIZ];
        _connectedCanonicals[_k[i]-1] = SetAlloc(_Bk);
        _numCanon[_k[i]-1] = canonListPopulate(BUF, _canonList[_k[i]-1], _connectedCanonicals[_k[i]-1], _k[i]);
        _maxNumCanon = MAX(_maxNumCanon, _numCanon[_k[i]-1]);  // set max number of canonicals for a k
        _K[_k[i]-1] = (short int*) aligned_alloc(8192, MAX(_Bk * sizeof(short int), 8192));
        assert(_K[_k[i]-1] != NULL);
        mapCanonMap(BUF, _K[_k[i]-1], _k[i]);
        sprintf(BUF, CANON_DIR "/perm_map%d.bin", _k[i]);
        int pfd = open(BUF, 0*O_RDONLY);
        kperm *Pf = Mmap(Permutations, _Bk*sizeof(Permutations[0]), pfd);
        assert(Pf == Permutations);
    }
}

int getBinSize(int n, int GDVcolumn[n], int* scratchspace){
    // sort
    // getIQR
    // return  (2IQR)/(n^1/3)
    assert(n>=5);

    int* sorted = scratchspace;
    memcpy(sorted, GDVcolumn, n * sizeof(int));
    qsort(sorted, n, sizeof(int), compare_ints);

    int q1,q3;

    if (n%2 == 0){
        q1 = getMedian(sorted, 0, (n/2)-1);
        q3 = getMedian(sorted, n/2, n-1);
    }else{
        q1 = getMedian(sorted, 0, n/2);
        q3 = getMedian(sorted, n/2, n-1);
    }

    assert(q3 >= q1);
    double obs = ((double) (2 * (q3-q1)) / cbrt((double) n));
    int returnVal = (int) ceil(obs);
    returnVal = MAX(returnVal, 1);
    return returnVal;
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

double FastEuclideanObjective(double oldcost, double olddelta, double change){
    double unchanged = SQR(oldcost) - SQR(olddelta);
    double newcost_sq = unchanged + SQR(olddelta+change);
    return sqrt(newcost_sq);
}

double FastDiffObjective(double oldcost, int k, int canonNum, int D0, int D1, int change){
    double oldratio, newratio;
    assert(abs(change) == 1);
    int diff;

    if (!SetIn(_connectedCanonicals[k-1], canonNum))  // for disconnected graphlets
        return oldcost;

    // compute old ratio
    diff = abs(D1-change-D0);
    if (D0 == 0){
        if (diff == 0)
            oldratio = (double) 1;  // 0/0
        else
            oldratio = (double) diff; // diff/0
    }else{
        oldratio = (double) diff/D0;  // diff/target
    }

    // compute new ratio
    diff = abs(D1-D0);
    if (D0 == 0){
        if (diff == 0)
            newratio = (double) 1;
        else
            newratio = (double) diff;
    }else{
        newratio = (double) diff/D0;
    }

    double returnVal = oldcost - oldratio + newratio;
    return returnVal;
}

double FastSGKObjective(double oldcost, int D0, int D1, int change){
    // D0 and D1 are the graphlet counts in D[0] and D[1] for a particular k and canonNum
    double oldratio, newratio;
    //int newTotalCanons = totalCanons;
    int m,M;
    assert(abs(change) == 1);

    // compute old ratio
    m = MIN(D0, D1 - change);  // D1-change is the older value of D1 (D1 is the count after decrementing/incrementing in ReBLANT)
    M = MAX(D0, D1 - change);
    if (M==0){
        //oldratio = 0;
        //if (change == 1) newTotalCanons++;
        oldratio = 1;
    }else{
        oldratio = (double) m/M;
    }

    // compute new ratio
    m = MIN(D0, D1);
    M = MAX(D0, D1);
    if (M==0){
        //newratio = 0;
        //if (change == -1) newTotalCanons--;
        newratio = 1;
    }else{
        newratio = (double) m/M;
    }
    
    oldcost = 1 - oldcost;
    double unchanged = totalCanons * oldcost;
    //assert(newTotalCanons > 0);
    //double returnVal = (double) (unchanged - oldratio + newratio)/newTotalCanons;
    double returnVal = (unchanged - oldratio + newratio)/totalCanons;
    //totalCanons = newTotalCanons;
    returnVal = 1 - returnVal;

    assert((returnVal >= 0) && (returnVal <= 1));
    return returnVal;
}

void Revert(int ***BLANT, int D[2][maxK][_maxNumCanon], Dictionary histograms[2][maxK][_maxNumCanon], int binsize[2][maxK][_maxNumCanon], int GDV[2][maxK][_maxNumCanon][_maxNodes], RevertStack* rvStack){
    // restore the BLANT line
    // restore D vectorS AND _numDisconnectedGraphlets
    // restore GDV[1] matrix and the histograms
    Change change;
    int n, node, key, value, b;

    while (pop(rvStack, &change) == 0){

        if (USING_GDV_OBJECTIVE){
            // GDV updates
            for(n=1; n<=change.k; n++){
                node = BLANT[change.k-1][change.linenum][n];            
                
                b = binsize[1][change.k-1][change.new];
                key = GDV[1][change.k-1][change.new][node];
                key = (int) ((int) key/b) * b;
                value = dictionary_get(&histograms[1][change.k-1][change.new], key, 0);
                assert(value > 0);
                dictionary_set(&histograms[1][change.k-1][change.new], key, value-1);
                GDV[1][change.k-1][change.new][node] -= 1; // main
                key = GDV[1][change.k-1][change.new][node];
                key = (int) ((int) key/b) * b;
                value = dictionary_get(&histograms[1][change.k-1][change.new], key, 0);
                dictionary_set(&histograms[1][change.k-1][change.new], key, value+1);                


                b = binsize[1][change.k-1][change.original];
                key = GDV[1][change.k-1][change.original][node];
                key = (int) ((int) key/b) * b;
                value = dictionary_get(&histograms[1][change.k-1][change.original], key, 0);
                assert(value >= 0);
                dictionary_set(&histograms[1][change.k-1][change.original], key, value-1);
                GDV[1][change.k-1][change.original][node] += 1; // main 
                key = GDV[1][change.k-1][change.original][node];
                key = (int) ((int) key/b) * b;
                value = dictionary_get(&histograms[1][change.k-1][change.original], key, 0);
                dictionary_set(&histograms[1][change.k-1][change.original], key, value+1);

            }
        }

        Boolean wasConnected = SetIn(_connectedCanonicals[change.k-1], BLANT[change.k-1][change.linenum][0]);
        Boolean  isConnected = SetIn(_connectedCanonicals[change.k-1], change.original);
        BLANT[change.k-1][change.linenum][0] = change.original;
        if(wasConnected && !isConnected) 
            ++_numDisconnectedGraphlets;
        if(!wasConnected && isConnected) 
            --_numDisconnectedGraphlets;
        --D[1][change.k-1][change.new];
        ++D[1][change.k-1][change.original];
    }
}

double AdjustDegree(int x, int y, int connected, GRAPH* G, int maxdegree, int Degree[2][maxdegree+1], double oldcost){
    // returns the new absolute cost of degree-distribution
    // Should be called AFTER calling GraphConnect or GraphDisconnect

    assert(abs(connected) == 1);  // connected = +1 if nodes were connected, else -1
    assert(oldcost >= 0);
    double newcost, olddelta, change;

    // for node x
    int old_deg_x = G->degree[x] - connected;
    olddelta = Degree[1][old_deg_x] - Degree[0][old_deg_x];
    newcost = FastEuclideanObjective(oldcost, olddelta, -1);
    olddelta = Degree[1][G->degree[x]] - Degree[0][G->degree[x]];
    newcost = FastEuclideanObjective(newcost, olddelta, 1);
    Degree[1][old_deg_x] -= 1;
    Degree[1][G->degree[x]] += 1;

    // for node y
    int old_deg_y = G->degree[y] - connected;
    olddelta = Degree[1][old_deg_y] - Degree[0][old_deg_y];
    newcost = FastEuclideanObjective(newcost, olddelta, -1);
    olddelta = Degree[1][G->degree[y]] - Degree[0][G->degree[y]];
    newcost = FastEuclideanObjective(newcost, olddelta, 1);
    Degree[1][old_deg_y] -= 1;
    Degree[1][G->degree[y]] += 1;

    // sanity check
    // number of elements
    int i, elts=0;
    for (i=0; i<=maxdegree; i++)
        elts += Degree[1][i];
    assert(elts == (G->n));

    assert(newcost >= 0);
    return newcost;
}

double ReBLANT(int D[2][maxK][_maxNumCanon], Dictionary histograms[2][maxK][_maxNumCanon], int binsize[2][maxK][_maxNumCanon], int GDV[2][maxK][_maxNumCanon][_maxNodes], GRAPH *G, SET ***samples, int ***Varrays, int ***BLANT, int v1, int v2, double oldcost, RevertStack* rvStack){

    int i, j, line, s, n;
    int canon, key, b, value, node;
    static TINY_GRAPH *g[maxK];

    double olddelta, change;
    double newcost = oldcost;

    for (i=0; i<maxK; i++){
        if (_k[i] == -1)
            break;
        int k = _k[i];
        
        // allocate a tiny graph 
        if (!g[k-1]) 
            g[k-1] = TinyGraphAlloc(k);

        for (s=1; line=Varrays[k-1][v1][s], s<=Varrays[k-1][v1][0]; s++)            
            if(SetIn(samples[k-1][v2], line)){
            
                // decrement a graphlet
                canon = BLANT[k-1][line][0];
                b = binsize[1][k-1][canon];
                if (b<=0){
                    fprintf(stderr, "invalid bin size in reblant =%d\n", b);
                    assert(b>0);
                }
                Boolean wasConnected = SetIn(_connectedCanonicals[k-1], BLANT[k-1][line][0]);
                --D[1][k-1][canon];

                if (USING_GDV_OBJECTIVE){
                    for (n=1; n<=k; n++){
                        node = BLANT[k-1][line][n];

                        // existing GDV key gets decremented
                        change = -1;
                        key = GDV[1][k-1][canon][node];
                        key = (int) ((int) key/b) * b;
                        value = dictionary_get(&(histograms[1][k-1][canon]), key, 0);
                        assert(value > 0);
                        olddelta = value - dictionary_get(&(histograms[0][k-1][canon]), key, 0);  // difference in the histograms
                        dictionary_set(&(histograms[1][k-1][canon]), key, value-1);
                        if ((!IGNORE_DISCONNECTED_GRAPHLETS) || wasConnected) 
                            newcost = FastEuclideanObjective(newcost, olddelta, change); // update cost

                        GDV[1][k-1][canon][node] -= 1;
                        assert(GDV[1][k-1][canon][node] >= 0);

                        // new GDV key gets incremented
                        change = 1;
                        key = GDV[1][k-1][canon][node];
                        key = (int) ((int) key/b) * b;
                        value = dictionary_get(&(histograms[1][k-1][canon]), key, 0);
                        olddelta = value - dictionary_get(&(histograms[0][k-1][canon]), key, 0);
                        dictionary_set(&(histograms[1][k-1][canon]), key, value+1);
                        if ((!IGNORE_DISCONNECTED_GRAPHLETS) || wasConnected) 
                            newcost = FastEuclideanObjective(newcost, olddelta, change);  // update cost
                    }

                }else{
                    // olddelta is for the Incremental Objective; this will change depending on the Objective
                    change = -1;
                    olddelta = (D[1][k-1][canon] + 1) - D[0][k-1][canon];

                    //newcost = FastDiffObjective(newcost, k, canon, D[0][k-1][canon], D[1][k-1][canon], change);
                    //newcost = FastSGKObjective(newcost, D[0][k-1][canon], D[1][k-1][canon], change);
                    //newcost = FastEuclideanObjective(newcost, olddelta, change);
                }

                // Change object (to be pushed on the stack)
                Change newchange;
                newchange.k = k; 
                newchange.linenum = line; 
                newchange.original = (int) canon;

                TinyGraphInducedFromGraph(g[k-1], G, &(BLANT[k-1][line][1])); // address of the Varray without element 0
                //Boolean wasConnected = SetIn(_connectedCanonicals[k-1], BLANT[k-1][line][0]);
                BLANT[k-1][line][0] = _K[k-1][TinyGraph2Int(g[k-1], k)];
                Boolean  isConnected = SetIn(_connectedCanonicals[k-1], BLANT[k-1][line][0]);
                if(wasConnected && !isConnected) 
                    ++_numDisconnectedGraphlets;
                if(!wasConnected && isConnected) 
                    --_numDisconnectedGraphlets;

                // increment a graphlet
                canon = BLANT[k-1][line][0];
                b = binsize[1][k-1][canon];
                assert(b>0);
                ++D[1][k-1][canon];

                if (USING_GDV_OBJECTIVE){
                    for (n=1; n<=k; n++){
                        node = BLANT[k-1][line][n];

                        // existing GDV key gets decremented
                        change = -1;
                        key = GDV[1][k-1][canon][node];
                        key = (int) ((int) key/b) * b;
                        value = dictionary_get(&(histograms[1][k-1][canon]), key, 0); 
                        assert(value > 0);
                        olddelta = value - dictionary_get(&(histograms[0][k-1][canon]), key, 0);  // difference in the histograms
                        dictionary_set(&(histograms[1][k-1][canon]), key, value-1);
                        if ((!IGNORE_DISCONNECTED_GRAPHLETS) || isConnected)
                            newcost = FastEuclideanObjective(newcost, olddelta, change);  // update cost              

                        GDV[1][k-1][canon][node] += 1;

                        // new GDV key gets incremented
                        change = 1;
                        key = GDV[1][k-1][canon][node];
                        key = (int) ((int) key/b) * b;
                        value = dictionary_get(&(histograms[1][k-1][canon]), key, 0);
                        olddelta = value - dictionary_get(&(histograms[0][k-1][canon]), key, 0);  // difference in the histograms
                        dictionary_set(&(histograms[1][k-1][canon]), key, value+1);
                        if ((!IGNORE_DISCONNECTED_GRAPHLETS) || isConnected)
                            newcost = FastEuclideanObjective(newcost, olddelta, change);  // update cost
                    }

                }else{
                    change = 1;
                    olddelta = (D[1][k-1][canon] - 1) - D[0][k-1][canon];
                    //newcost = FastDiffObjective(newcost, k, canon, D[0][k-1][canon], D[1][k-1][canon], change);
                    //newcost = FastSGKObjective(newcost, D[0][k-1][canon], D[1][k-1][canon], change);
                    //newcost = FastEuclideanObjective(newcost, olddelta, change);
                }

                // change object
                newchange.new = (int) canon;
                assert(push(rvStack, newchange) == 0);
            }          
             
        /*{
        int testCount = 0;
        for(j=0; j<_numCanon[k-1]; j++)
            testCount += D[1][k-1][j];
        assert(testCount == _numSamples);
        }*/
    }

    assert(newcost >= 0);
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

double Objective(double abscosts[NUMPROPS]){
    // returns WEIGHTED and NORMALIZED total costs, given current absolute costs
    double cost = 0;
    int i = 0;
    
    for(i=0; i<NUMPROPS; i++){
        assert(abscosts[i] >= 0);
        assert(max_abscosts[i] > 0);
        cost += ((double) weights[i] * (abscosts[i] / max_abscosts[i]));
    }
    
    assert(cost >= 0);
    return cost;
}

double GraphletEuclideanObjective(int D[2][maxK][_maxNumCanon]){
    // returns ABSOLUTE cost
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
    assert(returnVal >= 0);
    return returnVal; //exp(logP);
}

double DiffObjective(int D[2][maxK][_maxNumCanon]){
    // returns ABSOLUTE cost
    int i,j;
    double sum = 0;

    for(i=0; i<maxK; i++){
        int k = _k[i];
        if (k == -1) 
            break;        
        for (j=0; j<_numCanon[k-1]; j++){
            int diff = abs(D[0][k-1][j] - D[1][k-1][j]);
            int target = D[0][k-1][j];
            if (SetIn(_connectedCanonicals[k-1], j)){ // only count if canonical is connected
                if (target == 0){
                    if (diff == 0)
                        sum += (double) 1;
                    else
                        sum += (double) diff;
                }else{
                    sum += (double) diff/target;
                }
            }
        }
    }

    return sum;
}

double SGKObjective(int D[2][maxK][_maxNumCanon]){
    // returns ABSOLUTE cost
    int i,j;
    double sum = 0;
    totalCanons = 0;

    for(i=0; i<maxK; i++){
        int k = _k[i];
        if (k == -1) 
            break;
        
        assert(_numCanon[k-1] >= 0);
        totalCanons += _numCanon[k-1];

        for (j=0; j<_numCanon[k-1]; j++){
            int m,M;
            m = MIN(D[0][k-1][j], D[1][k-1][j]);
            M = MAX(D[0][k-1][j], D[1][k-1][j]);
            if (M==0){
                m=M=1;
                sum += 1;
                //totalCanons -= 1;  // ignore denominator, when M==0
            }
            else
                sum += ((double) m/M);
            
        }
    }

    assert(totalCanons > 0);
    double returnVal = sum/totalCanons;  // totalCanons is a global variable
    returnVal = 1-returnVal;
    assert((returnVal>=0) && (returnVal<=1));
    return returnVal;
}

double DegreeDistObjective(int maxdegree, int Degree[2][maxdegree+1]){
    // returns ABSOLUTE cost
    double sum = 0;
    int i;
    for(i=0; i<=maxdegree; i++)
        sum += SQR((double) Degree[0][i] - Degree[1][i]);

    double returnVal = sqrt(sum);
    assert(returnVal == returnVal);
    assert(returnVal >= 0);
    return returnVal;
}

double GDVObjective(Dictionary histograms[2][maxK][_maxNumCanon]){
    // returns the cost b/w every graphlet's GDV histogram
    double sum = 0;
    int j,k,canon;
    int key_tar, val_tar, key_syn, val_syn;

    Dictionary hist_tar, hist_syn;
    KeyValue *iter_tar, *iter_syn;

    for(j=0; j<maxK; j++){
        k = _k[j];
        if (k == -1) break;
        
        for(canon=0; canon < _numCanon[k-1]; canon++){

            // skip this canon if it is disconnected
            if ((IGNORE_DISCONNECTED_GRAPHLETS) && (!SetIn(_connectedCanonicals[k-1], canon))) 
                continue;

            // the 2 histograms (target & synthetic)
            hist_tar = histograms[0][k-1][canon];
            iter_tar = getIterator(&hist_tar);

            hist_syn = histograms[1][k-1][canon];
            iter_syn = getIterator(&hist_syn);

            double newsum = sum;
            // find 'difference' b/w the 2 histograms
            // iterate over all keys in target (these could be in synthetic or not)
            while(getNext(&iter_tar, &key_tar, &val_tar) == 0){
                assert((key_tar >= 0) && (val_tar >= 0));
                key_syn = key_tar;
                val_syn = dictionary_get(&hist_syn, key_syn, 0);  // default value is 0
                newsum += (double) SQR(val_tar - val_syn);
            }

            // iterate over all keys *ONLY* in synthetic
            while(getNext(&iter_syn, &key_syn, &val_syn) == 0){
                assert((key_syn >= 0) && (val_syn >= 0));
                key_tar = key_syn;
                val_tar = dictionary_get(&hist_tar, key_tar, -1);  // -1 will only be returned if the key doesn't exist in target
                if (val_tar == -1)
                    newsum += (double) SQR(val_syn);
            }

            assert(newsum >= sum);  // sum cannot decrease
            sum = newsum;
            assert(sum >= 0);
        }
    }

    double returnval = (double) sqrt(sum);
    assert(returnval == returnval);
    assert(returnval >= 0);
    return returnval;
}


int main(int argc, char *argv[]){
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

    // GDV matrices
    _maxNodes = MAX(G[0]->n, G[1]->n);
    int GDV[2][maxK][_maxNumCanon][_maxNodes];  // 4 dimensional
    for(i=0; i++; i<2){
        for(j=0; j<maxK; j++){
            if (_k[j] == -1) break;
            int l;
            for (l=0; l<_numCanon[_k[j]-1]; l++){
                int m;
                for (m=0; m < G[i]->n; m++)
                    GDV[i][_k[j]-1][l][m] = 0;  // initialize to 0
            }
        }
    }

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
                for (l=0; l<=_k[j]; l++){
                    assert(1 == fscanf(fp, "%d", &(BLANT[i][_k[j]-1][line][l])));
                    if (l>0){
                        GDV[i][_k[j]-1][BLANT[i][_k[j]-1][line][0]][BLANT[i][_k[j]-1][line][l]] += 1; // update the GDV 
                    }
                }
                assert(BLANT[i][_k[j]-1][line][0] < _maxNumCanon);
                ++D[i][_k[j]-1][BLANT[i][_k[j]-1][line][0]]; // squiggly plot
            }
            assert(fscanf(fp, "%d", &line) < 1); // ensure there's nothing left to read.
            fclose(fp);
            optind++;
        }
    
    }

    // sanity check - squiggly vectors
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

    // sanity check - GDV vectors
    for(i=0; i<2; i++){
        int l,m;
        for(j=0; j<maxK; j++){
            int matrixsum = 0;
            if (_k[j] == -1) break;
            for (l=0; l<_numCanon[_k[j]-1]; l++){
                int columnsum = 0;
                for (m=0; m < G[i]->n; m++){
                    matrixsum += GDV[i][_k[j]-1][l][m];
                    columnsum += GDV[i][_k[j]-1][l][m];
                }
                assert(columnsum == (D[i][_k[j]-1][l] * _k[j]));
            }
            assert(matrixsum == (_numSamples * _k[j]));
        }
    }

    
    Dictionary histograms[2][maxK][_maxNumCanon];  // Histograms, derived from GDV
    int binsize[2][maxK][_maxNumCanon];  // // Bin-size for these histograms

    int* scratchspace = (int*) malloc(_maxNodes * sizeof(int));  // used for sorting the GDV column
    for(i=0; i<2; i++){
        for(j=0; j<maxK; j++){
            if (_k[j] == -1) break;
            
            int l, b;
            for(l=0; l<_numCanon[_k[j]-1]; l++){ // for every graphlet

                if (i == 1)
                    binsize[1][_k[j]-1][l] = binsize[0][_k[j]-1][l];  // synthetic gets the same binsize as the corresponding target
                else
                    if (!SetIn(_connectedCanonicals[_k[j]-1], l))
                        binsize[0][_k[j]-1][l] = 1;  // disconnected graphlet. All GDV counts will be 0
                    else
                        binsize[0][_k[j]-1][l] = getBinSize(G[0]->n, GDV[0][_k[j]-1][l], scratchspace);

                b = binsize[i][_k[j]-1][l];
                assert(b>0);

                Dictionary* this = &(histograms[i][_k[j]-1][l]);
                dictionary_create(this);
                int n, key, prev;
                
                for (n=0; n < G[i]->n; n++){  // traverse the nodes involved in a particular graphlet
                    key = GDV[i][_k[j]-1][l][n];  // actual key value
                    key = (int) ((int) key/b) * b;  // binned key value
                    prev = dictionary_get(this, key, 0);
                    dictionary_set(this, key, prev+1);
                }
            }
        }
    }
    free(scratchspace);

    // sanity check bin size
    for (i=0; i<2; i++){
        for(j=0; j<maxK; j++){
            if (_k[j] == -1)
                break;
            int l=0;
            for(l=0; l<_numCanon[_k[j]-1]; l++){
                assert(binsize[i][_k[j]-1][l] > 0);
            }
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
                SetAdd(samples[_k[i]-1][BLANT[1][_k[i]-1][line][j]], line);    // SetAdd(samples[k][nodenum], line);
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

    // degree distribution vectors
    int maxdegree = -1;
    for (i=0; i<2; i++)
        for (j=0; j < G[i]->n; j++)
            if (G[i]->degree[j] > maxdegree)
                maxdegree = G[i]->degree[j];
    maxdegree = MAX(maxdegree, G[0]->n);
    maxdegree = MAX(maxdegree, G[1]->n);
    
    // initialize space
    int Degree[2][maxdegree+1];   // indexing: nodes with degree=5, are at index 5
    for (i=0; i<2; i++)
        for (j=0; j <= maxdegree; j++)
            Degree[i][j] = 0;
    for (i=0; i<2; i++)
        for (j=0; j < G[i]->n; j++)
            Degree[i][G[i]->degree[j]] += 1;

    // degree dist sanity check
    for (i=0; i<2; i++){
        int elts = 0;
        for(j=0; j<=maxdegree; j++)
            elts += Degree[i][j];
        assert(elts == (G[i]->n));
    }

    RevertStack uv, xy;  // uv gets disconnected and xy gets connected
    create_stack(&uv, maxK * _numSamples);
    create_stack(&xy, maxK * _numSamples);

    max_abscosts[0] = DegreeDistObjective(maxdegree, Degree);
    max_abscosts[1] = GDVObjective(histograms);
    //max_abscosts[1] = GraphletEuclideanObjective(D);
    //max_abscosts[1] = SGKObjective(D);
    //max_abscosts[1] = DiffObjective(D);
 
    double abscosts[NUMPROPS];
    abscosts[0] = max_abscosts[0];
    abscosts[1] = max_abscosts[1];

    fprintf(stderr, "Starting ABSOLUTE costs, DegreeDist: %g, Graphlets: %g\n", abscosts[0], abscosts[1]);

    // while(not done---either some number of iterations, or objective function says we're too far away)
    double cost = Objective(abscosts), startCost = cost, newCost, maxCost = cost;  // evaluate Objective() once, at the start. 
    assert(cost == cost);
    long int sa_iter = 0;
    double pBad, unif_random;
    float temperature;  // it's okay to overflow and become 0
    
    while((startCost - cost)/startCost < 0.5) // that's enough progress otherwise we're over-optimizing at this sample size
    {
    int edge = drand48() * G[1]->numEdges;
    int u1, u2, v1 = G[1]->edgeList[2*edge], v2 = G[1]->edgeList[2*edge+1];
    do {
        u1 = drand48()*G[1]->n;
        do u2 = drand48()*G[1]->n; while(u1==u2);
    } while(GraphAreConnected(G[1], u1, u2)); // find a non-edge  u1,u2
    assert(GraphAreConnected(G[1], v1, v2));  // find edge v1,v2

    // initialize new stacks
    assert(init_stack(&uv) == 0);
    assert(init_stack(&xy) == 0);
    double newcosts[NUMPROPS];

    GraphDisconnect(G[1], v1, v2); // remove edge e from Gs
    newcosts[0] = AdjustDegree(v1, v2, -1, G[1], maxdegree, Degree, abscosts[0]);
    newcosts[1] = ReBLANT(D, histograms, binsize, GDV, G[1], samples, Varrays, BLANT[1], v1, v2, abscosts[1], &uv);
    
    GraphConnect(G[1], u1, u2); // add an edge to Gs
    newcosts[0] = AdjustDegree(u1, u2, 1, G[1], maxdegree, Degree, newcosts[0]);
    newcosts[1] = ReBLANT(D, histograms, binsize, GDV, G[1], samples, Varrays, BLANT[1], u1, u2, newcosts[1], &xy);

    //fprintf(stderr, "\nthese 2 numbers should be the same - %g %g", newcosts[1], GDVObjective(histograms));

    newCost = Objective(newcosts);
    maxCost = MAX(maxCost, newCost);
    assert(newCost == newCost);
    static int same;

#if 1 // HILLCLIMB
    if(newCost < cost)
#else
#define SA_START_TEMP 1.0
#define SA_DECAY 1e-4
    temperature = SA_START_TEMP * exp(-SA_DECAY*500*sa_iter);
    unif_random = drand48();
    pBad = exp((-(newCost - cost)/maxCost)/temperature);
    // assert (newCost < cost <=> pBad > 1)
    assert((newCost < cost) == (pBad > 1));
    if(newCost < cost || unif_random < pBad)
#endif
    {
        //fprintf(stderr,"A");fflush(stderr);
        static double printVal, printInterval;
        if(fabs(newCost - printVal)/printVal >= 0.02 || ++printInterval > PRINT_INTERVAL)
        {
        fprintf(stderr, "\ntemp %g cost %g newCost %g maxCost %g pBad %g numDis %d", temperature, cost, newCost, maxCost, pBad, _numDisconnectedGraphlets);
        //fprintf(stderr, "%g ", newCost);
        printVal = newCost;
        printInterval = 0;
        }
        cost = newCost;
        memcpy(abscosts, newcosts, NUMPROPS * (sizeof(double)));
        same = 0;
    }
    else // revert
    {
        //fprintf(stderr,"R");fflush(stderr);
        static double printVal, printInterval;
        if(fabs(newCost - printVal)/printVal >= 0.02 || ++printInterval > PRINT_INTERVAL)
        {
        fprintf(stderr, "\ntemp %g cost %g newCost %g maxCost %g pBad %g numDis %d", temperature, cost, newCost, maxCost, pBad, _numDisconnectedGraphlets);
        //fprintf(stderr, "%g ", newCost);
        printVal = newCost;
        printInterval = 0;
        }
        ++same;

        GraphDisconnect(G[1], u1, u2);
        AdjustDegree(u1, u2, -1, G[1], maxdegree, Degree, newcosts[0]);

        GraphConnect(G[1], v1, v2);
        AdjustDegree(v1, v2, 1, G[1], maxdegree, Degree, newcosts[0]);

        // revert changes to blant file and D vectors
        Revert(BLANT[1], D, histograms, binsize, GDV, &xy);
        Revert(BLANT[1], D, histograms, binsize, GDV, &uv);
    }

    if(same > _stagnated || _numDisconnectedGraphlets >= _numSamples*10){
        fprintf(stderr, "stagnated!\n");
        break;
    }
    ++sa_iter;
    }
    fprintf(stderr,"\n");
    for(i=0; i < G[1]->numEdges; i++)
        printf("%d %d\n", G[1]->edgeList[2*i], G[1]->edgeList[2*i+1]);
}
