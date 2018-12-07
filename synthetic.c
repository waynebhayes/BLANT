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
#include "syntheticDS.h"

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
static int maxdegree = -1;
static int _canonList[maxK][MAX_CANONICALS];
static int _stagnated = 1000, _numDisconnectedGraphlets;

#define IGNORE_DISCONNECTED_GRAPHLETS 1
#define PRINT_INTERVAL 10000
#define KHOP_INTERVAL 5000

#define NUMPROPS 6
#define GraphletEuclidean 0
#define SGK 1
#define SGKDiff 2
#define GraphletGDV 3
#define DegreeDist 4
#define ClustCoff 5

static double weights[NUMPROPS] =
// weights: 0 GraphletEuclidean; 1 SGK; 2 SGKDiff; 3 GDV;  4 DegreeDist; 5 ClustCoff
            {1,                0,     0,      0,   0,         0};  // SGK should be avoided for now
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
    double unchanged = (oldcost * oldcost) - (olddelta * olddelta);
    unchanged = MAX(0, unchanged);
    double newcost_sq = unchanged + SQR(olddelta+change);
    newcost_sq = MAX(0, newcost_sq);
    return sqrt(newcost_sq);
}

double FastGDVObjective(double oldcost, int olddelta, double change){
    double unchanged = (oldcost * oldcost) - (olddelta * olddelta);
    unchanged = MAX(0, unchanged);
    double newcost_sq = unchanged + (double) SQR(olddelta+change);
    newcost_sq = MAX(0, newcost_sq);
    return sqrt(newcost_sq);
}

double FastSGKDiffObjective(double oldcost, int k, int canonNum, int D0, int D1, int change){
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

double AdjustGDV(int k, int canon, int change, int line[k+1], Dictionary GDVhistograms[2][maxK][_maxNumCanon], int GDVbinsize[2][maxK][_maxNumCanon], int GDV[2][maxK][_maxNumCanon][_maxNodes], double oldcost){

    // k and canon identify a graphlet
    // change = +1 if this graphlet is the "NEW" graphlet. -1 is this is the "OLD" graphlet (w.r.t ReBLANT code)
    // line is the BLANT sample line (contains which nodes are involved)

    assert(abs(change) == 1);
    assert(line[0] == canon);

    int n, b, c, node, key, value, olddelta;
    double newcost = oldcost;

    b = GDVbinsize[1][k-1][canon];

    for (n=1; n<=k; n++){
        node = line[n];

        // existing GDV key gets decremented
        c = -1;
        key = GDV[1][k-1][canon][node];
        key = (int) ((int) key/b) * b;
        value = dictionary_get(&(GDVhistograms[1][k-1][canon]), key, 0);
        assert(value > 0);
        olddelta = value - dictionary_get(&(GDVhistograms[0][k-1][canon]), key, 0);  // difference in the GDVhistograms
        dictionary_set(&(GDVhistograms[1][k-1][canon]), key, value-1);
        newcost = FastGDVObjective(newcost, olddelta, c); // update cost

        GDV[1][k-1][canon][node] += change;
        if (change == -1)
            assert(GDV[1][k-1][canon][node] >= 0);

        // new GDV key gets incremented
        c = 1;
        key = GDV[1][k-1][canon][node];
        key = (int) ((int) key/b) * b;
        value = dictionary_get(&(GDVhistograms[1][k-1][canon]), key, 0);
        olddelta = value - dictionary_get(&(GDVhistograms[0][k-1][canon]), key, 0);
        dictionary_set(&(GDVhistograms[1][k-1][canon]), key, value+1);
        newcost = FastGDVObjective(newcost, olddelta, c);  // update cost
    }

    assert(newcost >= 0);
    return newcost;
}

double AdjustDegree(const int x, const int y, const int connected, GRAPH* G, int Degree[2][maxdegree+1], double oldcost){
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

double AdjustClustCoff(const int x, const int y, const int connected, GRAPH* G, int localConnections[2][_maxNodes], int* CCHistograms[2], int CCHistograms_size[2], double CCbinsize[2], double oldcost){

    assert(abs(connected) == 1);
    assert(CCHistograms_size[0] == CCHistograms_size[1]);
    int i, j, node, nodecount, c, nc2, oldkey, newkey;

    double oldcc, newcc;
    double newcost = oldcost;
    assert((newcost == newcost) && (newcost >= 0));

    SET* xn = SetAlloc(G->n);
    for(i=0; i < G->degree[x]; i++)
        if (G->neighbor[x][i] != y)
            SetAdd(xn, G->neighbor[x][i]);
    
    nodecount = 0;  // nodes which are connected both to x & y
    for(i=0; i < G->degree[y]; i++){
        node = G->neighbor[y][i];
        if ((node!=x) && (SetIn(xn, node))){
            nodecount += 1;
            assert(G->degree[node] >= 2);
            nc2 = (G->degree[node] * (G->degree[node] - 1))/2;
            c = localConnections[1][node];  // the original num of connections
            localConnections[1][node] = c + connected;  // the new num of connections
            if (nc2 > 0){
                oldcc = (double) c / (double) nc2;
                oldkey = (int) (oldcc / CCbinsize[1]);
                assert(oldkey < CCHistograms_size[1]);
                newcc = (double) (c+connected) / (double) nc2;
                newkey = (int) (newcc / CCbinsize[1]);
                assert(newkey < CCHistograms_size[1]);

                newcost = FastEuclideanObjective(newcost, (double) (CCHistograms[1][oldkey] - CCHistograms[0][oldkey]), (double) -1.0);
                assert((newcost == newcost) && (newcost >= 0));
                assert(CCHistograms[1][oldkey] > 0);
                CCHistograms[1][oldkey] -= 1;
                newcost = FastEuclideanObjective(newcost, (double) (CCHistograms[1][newkey] - CCHistograms[0][newkey]), (double) 1.0);
                assert((newcost == newcost) && (newcost >= 0));
                CCHistograms[1][newkey] += 1;
            }
            else{
                oldcc = 0.0;
                oldkey = 0;
                newcc = 0.0;
                newkey = 0;
            }
        }
    }

    SetFree(xn);


    // for x
    if ((G->degree[x] - connected) < 2){ // original degree
        c = localConnections[1][x];  // the old connections
        assert(c == 0);
        oldcc = 0;
        oldkey = 0;
    }else{
        nc2 = ((G->degree[x] - connected) * ((G->degree[x] - connected) - 1))/2;
        c = localConnections[1][x];  // the old connections
        oldcc = ((double) c) / ((double) nc2);
        oldkey = (int) (oldcc / CCbinsize[1]);
    }

    if(G->degree[x] < 2){  // current degree
        c += (connected * nodecount);
        assert(c == 0);
        newcc = 0;
        newkey = 0;
    }else{
        nc2 = (G->degree[x] * (G->degree[x] - 1))/2;
        c += (connected * nodecount);
        newcc = ((double) c) / ((double) nc2);
        newkey = (int) (newcc / CCbinsize[1]);
    }

    localConnections[1][x] = c;
    if(oldkey != newkey){
        newcost = FastEuclideanObjective(newcost, (double) (CCHistograms[1][oldkey] - CCHistograms[0][oldkey]), (double) -1.0);
        assert((newcost == newcost) && (newcost >= 0));
        assert(CCHistograms[1][oldkey] > 0);
        CCHistograms[1][oldkey] -= 1;
        newcost = FastEuclideanObjective(newcost, (double) (CCHistograms[1][newkey] - CCHistograms[0][newkey]), (double) 1.0);
        assert((newcost == newcost) && (newcost >= 0));
        CCHistograms[1][newkey] += 1;
    }


    // for y
    if ((G->degree[y] - connected) < 2){ // original degree
        c = localConnections[1][y];  // the old connections
        assert(c == 0);
        oldcc = 0;
        oldkey = 0;
    }else{
        nc2 = ((G->degree[y] - connected) * ((G->degree[y] - connected) - 1))/2;
        c = localConnections[1][y];  // the old connections
        oldcc = ((double) c) / ((double) nc2);
        oldkey = (int) (oldcc / CCbinsize[1]);
    }

    if(G->degree[y] < 2){  // current degree
        c += (connected * nodecount);
        assert(c == 0);
        newcc = 0;
        newkey = 0;
    }else{
        nc2 = (G->degree[y] * (G->degree[y] - 1))/2;
        c += (connected * nodecount);
        newcc = ((double) c) / ((double) nc2);
        newkey = (int) (newcc / CCbinsize[1]);
    }

    localConnections[1][y] = c;
    if(oldkey != newkey){
        newcost = FastEuclideanObjective(newcost, (double) (CCHistograms[1][oldkey] - CCHistograms[0][oldkey]), (double) -1.0);
        assert((newcost == newcost) && (newcost >= 0));
        assert(CCHistograms[1][oldkey] > 0);
        CCHistograms[1][oldkey] -= 1;
        newcost = FastEuclideanObjective(newcost, (double) (CCHistograms[1][newkey] - CCHistograms[0][newkey]), (double) 1.0);
        assert((newcost == newcost) && (newcost >= 0));
        CCHistograms[1][newkey] += 1;
    }

    return newcost;
}

void ReBLANT(int D[2][maxK][_maxNumCanon], Dictionary GDVhistograms[2][maxK][_maxNumCanon], int GDVbinsize[2][maxK][_maxNumCanon], int GDV[2][maxK][_maxNumCanon][_maxNodes], GRAPH *G, SET ***samples, int ***Varrays, int ***BLANT, int v1, int v2, double oldcost[4], RevertStack* rvStack){

    // updates cost ONLY for GraphletEuclidean, SGK, Diff, GDV
    int i, j, line, s, change;
    int canon;
    static TINY_GRAPH *g[maxK];
    double oldcanondiff, temp_gdv_newcost;

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
                Boolean wasConnected = SetIn(_connectedCanonicals[k-1], canon);
                oldcanondiff = D[1][k-1][canon] - D[0][k-1][canon];
                --D[1][k-1][canon];
                change = -1;

                // recompute cost
                temp_gdv_newcost = AdjustGDV(k, canon, change, BLANT[k-1][line], GDVhistograms, GDVbinsize, GDV, oldcost[3]);  // fix GDV
                if ((!IGNORE_DISCONNECTED_GRAPHLETS) || wasConnected){
                    oldcost[0] = FastEuclideanObjective(oldcost[0], oldcanondiff, change);
                    oldcost[1] = FastSGKObjective(oldcost[1], D[0][k-1][canon], D[1][k-1][canon], change);
                    oldcost[2] = FastSGKDiffObjective(oldcost[2], k, canon, D[0][k-1][canon], D[1][k-1][canon], change);
                    oldcost[3] = temp_gdv_newcost;
                }

                // Change object (to be pushed on the stack)
                Change newchange;
                newchange.k = k; 
                newchange.linenum = line; 
                newchange.original = (int) canon;


                TinyGraphInducedFromGraph(g[k-1], G, &(BLANT[k-1][line][1])); // address of the Varray without element 0
                BLANT[k-1][line][0] = _K[k-1][TinyGraph2Int(g[k-1], k)];
                Boolean  isConnected = SetIn(_connectedCanonicals[k-1], BLANT[k-1][line][0]);
                if(wasConnected && !isConnected) 
                    ++_numDisconnectedGraphlets;
                if(!wasConnected && isConnected) 
                    --_numDisconnectedGraphlets;


                // increment a graphlet
                canon = BLANT[k-1][line][0];
                oldcanondiff = D[1][k-1][canon] - D[0][k-1][canon];
                ++D[1][k-1][canon];
                change = 1;

                // recompute cost
                temp_gdv_newcost = AdjustGDV(k, canon, change, BLANT[k-1][line], GDVhistograms, GDVbinsize, GDV, oldcost[3]);  // fix GDV
                if ((!IGNORE_DISCONNECTED_GRAPHLETS) || wasConnected){
                    oldcost[0] = FastEuclideanObjective(oldcost[0], oldcanondiff, change);
                    oldcost[1] = FastSGKObjective(oldcost[1], D[0][k-1][canon], D[1][k-1][canon], change);
                    oldcost[2] = FastSGKDiffObjective(oldcost[2], k, canon, D[0][k-1][canon], D[1][k-1][canon], change);
                    oldcost[3] = temp_gdv_newcost;
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

    for(i=0; i<4; i++)
        assert(oldcost[i] >= 0);
}

void Revert(int ***BLANT, int D[2][maxK][_maxNumCanon], Dictionary GDVhistograms[2][maxK][_maxNumCanon], int GDVbinsize[2][maxK][_maxNumCanon], int GDV[2][maxK][_maxNumCanon][_maxNodes], RevertStack* rvStack){
    // restore the BLANT line
    // restore D vectorS AND _numDisconnectedGraphlets
    // restore GDV[1] matrix and the GDVhistograms
    Change change;
    int n, node, key, value, b;
    int* line;

    while (pop(rvStack, &change) == 0){
        line = BLANT[change.k-1][change.linenum];
        AdjustGDV(change.k, change.new, -1, line, GDVhistograms, GDVbinsize, GDV, 1000);  // the 1000 has no meaning
        Boolean wasConnected = SetIn(_connectedCanonicals[change.k-1], BLANT[change.k-1][change.linenum][0]);
        Boolean  isConnected = SetIn(_connectedCanonicals[change.k-1], change.original);
        BLANT[change.k-1][change.linenum][0] = change.original;
        AdjustGDV(change.k, change.original, 1, line, GDVhistograms, GDVbinsize, GDV, 1000);
        if(wasConnected && !isConnected) 
            ++_numDisconnectedGraphlets;
        if(!wasConnected && isConnected) 
            --_numDisconnectedGraphlets;
        --D[1][change.k-1][change.new];
        ++D[1][change.k-1][change.original];
    }
}

double Objective(double abscosts[NUMPROPS]){
    // returns WEIGHTED and NORMALIZED total costs, given current absolute costs
    double cost = 0;
    int i = 0;
    
    for(i=0; i<NUMPROPS; i++){
        assert(abscosts[i] >= 0);
        assert(max_abscosts[i] >= 0);
        if (max_abscosts[i] == 0)
            cost += ((double) weights[i] * (abscosts[i] / 1));
        else
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

double SGKDiffObjective(int D[2][maxK][_maxNumCanon]){
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

double GDVObjective(Dictionary GDVhistograms[2][maxK][_maxNumCanon]){
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

            // the 2 GDVhistograms (target & synthetic)
            hist_tar = GDVhistograms[0][k-1][canon];
            iter_tar = getIterator(&hist_tar);

            hist_syn = GDVhistograms[1][k-1][canon];
            iter_syn = getIterator(&hist_syn);

            double newsum = sum;
            // find 'difference' b/w the 2 GDVhistograms
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

double DegreeDistObjective(int Degree[2][maxdegree+1]){
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

void GetConnections(GRAPH *G, int localConnections[G->n]){
    // computes local clustering coefficient for every node
    int n, x, node, degree;
    int i, j, edges;
    
    // set for marking neighbors
    int scratch[G->n];
    for(i=0; i<G->n; i++){
        localConnections[i] = 0;
        scratch[i] = -1;
    }

    for(n=0; n < G->n; n++){
        edges = 0;
        degree = G->degree[n];

        if (degree < 2){
            localConnections[n] = 0;
            continue;
        }
        
        for(i=0; i<degree; i++)
            scratch[(G->neighbor[n])[i]] = n;
        
        for(i=0; i<degree; i++){
            node = (G->neighbor[n])[i];
            assert(node != n);
            for(j=0; j < G->degree[node]; j++){
                x = (G->neighbor[node])[j];
                if ((x!=n) && (scratch[x] == n))
                    edges += 1;
            }
        }

        edges = edges/2;
        localConnections[n] = edges;
    }
}

double ClustCoffObjective(int* CCHistograms[2], int CCHistograms_size[2]){
    // new implementation : euclidean distance b/w the 2 histograms (arrays)
    double cost = 0;
    double temp;
    int j;
    
    assert(CCHistograms_size[0] == CCHistograms_size[1]);

    for(j=0; j < CCHistograms_size[0]; j++){
        temp = (double) SQR(CCHistograms[0][j] - CCHistograms[1][j]);
        assert(temp >= 0);
        cost += temp;
    }

    double returnVal = sqrt(cost);
    assert(returnVal == returnVal);
    assert(returnVal >= 0);
    return returnVal;
}

void FindNodes(GRAPH* G, const Hops hops, int* u1, int* u2, int* v1, int* v2){
    // v1 v2 should be disconnected
    // u1 u2 should be connected
    int edge, x, y, d;

    // find existing edge
    do{
        edge = drand48() * G->numEdges;
        x = G->edgeList[2*edge];
        y = G->edgeList[2*edge+1]; 
    } while((G->degree[x] == 1) || (G->degree[y] == 1));  // don't disconnect nodes with degree=1
    assert(GraphAreConnected(G, x, y));
    memcpy(v1, &x, sizeof(int));
    memcpy(v2, &y, sizeof(int));

    // find existing non-edge
    if(hops.valid == 1){
        d = hops.lower + (int)((hops.upper - hops.lower + 1) * drand48());
        do {
            x = drand48()*G->n;
            y = getRandomNode(G, x, d);
        }while(GraphAreConnected(G, x, y));
    }else{
        do {
            x = drand48()*G->n;
            do y = drand48()*G->n; while(x==y);
        }while(GraphAreConnected(G, x, y));
    }
    assert(!GraphAreConnected(G, x, y));
    memcpy(u1, &x, sizeof(int));
    memcpy(u2, &y, sizeof(int));
}

void OptimumHops(Dictionary khop[2], Hops* hops){
    KeyValue* iter;
    int count, k, v, i, l, mi;

    int medians[2];
    int maxKeys[2] = {-1, -1};

    for(i=0; i<2; i++){
        count = 0;

        iter = getIterator(&(khop[i]));
        while((getNext(&iter, &k, &v)) == 0){
            maxKeys[i] = MAX(maxKeys[i], k);
            count += v;
        }

        mi = (int) count/2;
        count = 0;
        for(l=0; l<=maxKeys[i]; l++){
            count += dictionary_get(&(khop[i]), l, 0);
            if (count >= mi){
                medians[i] = l;
                break;
            }
        }
    }

    if(medians[0] < medians[1]){
        // make synthetic MORE small-world
        hops->valid = 1;
        hops->lower = medians[1];
        hops->upper = maxKeys[1];
    }else if(medians[1] < medians[0]){
        // make synthetic LESS small-world
        hops->valid = 1;
        hops->lower = 2;
        hops->upper = medians[1];
    }else{
        // all hop values are equally likely
        hops->valid = 0;
    }
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
    for(i=0; i<2; i++){
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
            long matrixsum = 0;
            if (_k[j] == -1) break;
            for (l=0; l<_numCanon[_k[j]-1]; l++){
                int columnsum = 0;
                for (m=0; m < G[i]->n; m++){
                    matrixsum += ((long) GDV[i][_k[j]-1][l][m]);
                    columnsum += GDV[i][_k[j]-1][l][m];
                }
                assert(columnsum == (D[i][_k[j]-1][l] * _k[j]));
            }
            assert(matrixsum == (((long)_numSamples) * _k[j]));
        }
    }

    Dictionary GDVhistograms[2][maxK][_maxNumCanon];  // GDVhistograms, derived from GDV
    int GDVbinsize[2][maxK][_maxNumCanon];  // // Bin-size for these GDVhistograms

    int* scratchspace = (int*) malloc(_maxNodes * sizeof(int));  // used for sorting the GDV column
    for(i=0; i<2; i++){
        for(j=0; j<maxK; j++){
            if (_k[j] == -1) break;
            
            int l, b;
            for(l=0; l<_numCanon[_k[j]-1]; l++){ // for every graphlet

                if (i == 1)
                    GDVbinsize[1][_k[j]-1][l] = GDVbinsize[0][_k[j]-1][l];  // synthetic gets the same GDVbinsize as the corresponding target
                else
                    if (!SetIn(_connectedCanonicals[_k[j]-1], l))
                        GDVbinsize[0][_k[j]-1][l] = 1;  // disconnected graphlet. All GDV counts will be 0
                    else
                        GDVbinsize[0][_k[j]-1][l] = getIntegerBinSize(G[0]->n, GDV[0][_k[j]-1][l], scratchspace);

                b = GDVbinsize[i][_k[j]-1][l];
                assert(b>0);

                Dictionary* this = &(GDVhistograms[i][_k[j]-1][l]);
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
                assert(GDVbinsize[i][_k[j]-1][l] > 0);
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
    maxdegree = MAX(G[0]->n, G[1]->n);
    for (i=0; i<2; i++)
        for (j=0; j < G[i]->n; j++)
            assert((G[i]->degree[j]) < (G[i]->n));
    
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

    // CLUSTERING COEFFICIENT
    int localConnections[2][_maxNodes];  // for every node, number of connections in it's neighborhood
    for (i=0; i<2; i++)
        GetConnections(G[i], localConnections[i]);   // populates the array with every nodes' local CONNECTIONS (not clustCoff)

    int* CCHistograms[2];  // stores num of elts in a bin
    int CCHistograms_size[2];  // stores num of bins
    double CCbinsize[2];  // stores bin size/width
    
    for(i=0; i<2; i++){
        int degree, nc2, key, value;
        double local_cc, b, localClustCoff[G[i]->n], scratchspace[G[i]->n];
        
        for(j=0; j < G[i]->n; j++){  // get local_CCoff from local_connections
            degree = G[i]->degree[j];
            nc2 = (degree * (degree-1))/2;
            if (nc2 <= 0) local_cc = 0;
            else local_cc = (double) localConnections[i][j] / (double) nc2;
            localClustCoff[j] = local_cc;
            assert(localClustCoff[j] >= 0);
        }

        if(i == 0)
            CCbinsize[i] = getDoubleBinSize(G[i]->n, localClustCoff, scratchspace);
        else
            CCbinsize[i] = CCbinsize[0];

        b = CCbinsize[i];
        assert(b > 0);

        CCHistograms_size[i] = ((int) (1/b)) + 1;
        CCHistograms[i] = (int*) Malloc(CCHistograms_size[i] * sizeof(int));
        for(j=0; j < CCHistograms_size[i]; j++)
            CCHistograms[i][j] = 0;

        for(j=0; j < G[i]->n; j++){
            key = (int) (localClustCoff[j]/b);
            assert(key < CCHistograms_size[i]);
            CCHistograms[i][key] += 1;
        }
    }

    Dictionary khop[2];
    Hops hops;

    RevertStack uv, xy;  // uv gets disconnected and xy gets connected
    create_stack(&uv, maxK * _numSamples);
    create_stack(&xy, maxK * _numSamples);

    max_abscosts[GraphletEuclidean] = GraphletEuclideanObjective(D);
    max_abscosts[SGK] = SGKObjective(D);
    max_abscosts[SGKDiff] = SGKDiffObjective(D);
    max_abscosts[GraphletGDV] = GDVObjective(GDVhistograms);
    max_abscosts[DegreeDist] = DegreeDistObjective(Degree);
    max_abscosts[ClustCoff] = ClustCoffObjective(CCHistograms, CCHistograms_size);
 
    double abscosts[NUMPROPS];
    memcpy(abscosts, max_abscosts, NUMPROPS * sizeof(double));

    fprintf(stderr, "Starting ABSOLUTE costs: GraphletEuclidean: %g, SGK: %g, GraphDiff: %g, GDV: %g, DegreeDist: %g, ClustCoff: %g\n", abscosts[0], abscosts[1], abscosts[2], abscosts[3], abscosts[4], abscosts[5]);

    // while(not done---either some number of iterations, or objective function says we're too far away)
    double cost = Objective(abscosts), startCost = cost, newCost, maxCost = cost;  // evaluate Objective() once, at the start. 
    assert(cost == cost);
    long int sa_iter = 0;
    double pBad, unif_random;
    float temperature;  // it's okay to overflow and become 0
    
    double newcosts[NUMPROPS];

    while((startCost - cost)/startCost < 0.5) // that's enough progress otherwise we're over-optimizing at this sample size
    {

    int u1, u2, v1, v2;
    if((sa_iter % KHOP_INTERVAL) == 0){  // compute k-hop
        sampleKHop(G[0], &(khop[0]), 0.8);
        sampleKHop(G[1], &(khop[1]), 0.8);
        OptimumHops(khop, &hops);
    }

    FindNodes(G[1], hops, &u1, &u2, &v1, &v2);

    // initialize new stacks
    assert(init_stack(&uv) == 0);
    assert(init_stack(&xy) == 0);
    memcpy(newcosts, abscosts, NUMPROPS * sizeof(double));

    GraphDisconnect(G[1], v1, v2); // remove edge e from Gs
    newcosts[DegreeDist] = AdjustDegree(v1, v2, -1, G[1], Degree, newcosts[4]);
    newcosts[ClustCoff] = AdjustClustCoff(v1, v2, -1, G[1], localConnections, CCHistograms, CCHistograms_size, CCbinsize, newcosts[ClustCoff]);
    ReBLANT(D, GDVhistograms, GDVbinsize, GDV, G[1], samples, Varrays, BLANT[1], v1, v2, newcosts, &uv);

    GraphConnect(G[1], u1, u2); // add an edge to Gs
    newcosts[DegreeDist] = AdjustDegree(u1, u2, 1, G[1], Degree, newcosts[4]);
    newcosts[ClustCoff] = AdjustClustCoff(u1, u2, 1, G[1], localConnections, CCHistograms, CCHistograms_size, CCbinsize, newcosts[ClustCoff]);
    ReBLANT(D, GDVhistograms, GDVbinsize, GDV, G[1], samples, Varrays, BLANT[1], u1, u2, newcosts, &xy);

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
        static double printVal=1e30, printInterval;
        if(fabs(cost - printVal)/printVal >= 0.02 || ++printInterval > PRINT_INTERVAL)
        {
            fprintf(stderr, "\ntemp %g cost %g newCost %g maxCost %g pBad %g numDis %d", temperature, cost, newCost, maxCost, pBad, _numDisconnectedGraphlets);
            //fprintf(stderr, "%g ", newCost);
            printVal = cost;
            printInterval = 0;
        }
        cost = newCost;
        memcpy(abscosts, newcosts, NUMPROPS * sizeof(double));
        // update max costs -- only if move is accepted
        for(i=0; i<NUMPROPS; i++)
            max_abscosts[i] = MAX(max_abscosts[i], abscosts[i]);

        same = 0;
    }
    else // revert
    {
        //fprintf(stderr,"R");fflush(stderr);
        static double printVal=1e30, printInterval;
        if(fabs(cost - printVal)/printVal >= 0.02 || ++printInterval > PRINT_INTERVAL)
        {
            fprintf(stderr, "\ntemp %g cost %g newCost %g maxCost %g pBad %g numDis %d", temperature, cost, newCost, maxCost, pBad, _numDisconnectedGraphlets);
            //fprintf(stderr, "%g ", newCost);
            printVal = cost;
            printInterval = 0;
        }
        ++same;

        GraphDisconnect(G[1], u1, u2);
        newcosts[DegreeDist] = AdjustDegree(u1, u2, -1, G[1], Degree, newcosts[DegreeDist]);
        newcosts[ClustCoff] = AdjustClustCoff(u1, u2, -1, G[1], localConnections, CCHistograms, CCHistograms_size, CCbinsize, newcosts[ClustCoff]);

        GraphConnect(G[1], v1, v2);
        newcosts[DegreeDist] = AdjustDegree(v1, v2, 1, G[1], Degree, newcosts[DegreeDist]);
        newcosts[ClustCoff] = AdjustClustCoff(v1, v2, 1, G[1], localConnections, CCHistograms, CCHistograms_size, CCbinsize, newcosts[ClustCoff]);
        
        // revert changes to blant file and D vectors
        Revert(BLANT[1], D, GDVhistograms, GDVbinsize, GDV, &xy);
        Revert(BLANT[1], D, GDVhistograms, GDVbinsize, GDV, &uv);
    }

    if(same > _stagnated || _numDisconnectedGraphlets >= _numSamples*10){
        fprintf(stderr, "stagnated!, iterations=%d\n", sa_iter);
        break;
    }
    ++sa_iter;
    }
    fprintf(stderr,"\n");
    for(i=0; i < G[1]->numEdges; i++)
        printf("%d %d\n", G[1]->edgeList[2*i], G[1]->edgeList[2*i+1]);
}
