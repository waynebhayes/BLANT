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

static int _k[MAX_K]; // stores what values of k have to be considered e.g. [3,4,5,6] or [3,5,7] or [2,7,8]. Will be followed by '-1's
static int _numCanon[MAX_K];  // canonicals for particular value of k. So for k=5, _numCanon[5-1] stores ~32
static SET *_connectedCanonicals[MAX_K];
static int _maxNumCanon = -1;  // max number of canonicals
static int _numSamples = -1;  // same number of samples in each blant index file
static int _numNodes = -1;  // number of nodes in the target/synthetic network
static int maxdegree = -1;  // equals _numNodes (a node can be connected to every other node)
static Gint_type _canonList[MAX_K][MAX_CANONICALS];
static int _stagnated = 1000, _numDisconnectedGraphlets;

#define PRINT_INTERVAL 100000
#define ISZERO(w) (fabs(w)<0.00001)
#define NC2(n) ((n* (n-1)) / 2)

/* khop distribution is computed (bfs from each node) and used to make
the synthetic more/less like a small-world network (or keep things as they are).
So, compute khop of target and synthetic & then compare the MEDIANS of the two distributions.
Now, do node-selection to make the network more/less small-worldy (details in GetNodes function)
This has the effect of tuning the centrality measures - node eccentricity; node & edge betweenness; KHOP; diameter; ....
*/
#define DEFAULT_KHOP_INTERVAL 400  // accepted-iterations after which khop distribution is computed (if things are to be kept the same)
#define FIX_KHOP_INTERVAL 200  // accepted-iterations after which khop distribution is computed (if network needs to made more/less small-worldy)
#define KHOP_QUALITY 0.8  // pick these fraction of nodes, do a BFS on each of these nodes (for computing KHOP)

// node-selection
#define NODE_SEL_ALWAYS_RANDOM 0
#define NODE_SEL_SHORT_PATH 1  // join & disconnect nodes which have more/less Shortest Paths going through them
#define NODE_SEL_BY_HOPS 2  // join & disconnect nodes which are a specific BFS hops apart (SLOW)
static int node_selection = NODE_SEL_SHORT_PATH;
// node-selection-strategy can be set using an env variable
// USAGE: export SYNTHETIC_NODE_SELECTION = 0   # 0 for random, 1 for shortestpaths, 2 for hops

// objective functions
#define IGNORE_DISCONNECTED_GRAPHLETS 1
#define NUMPROPS 7
#define GraphletEuclidean 0
#define GraphletKernel 1  // let u & v be two graphlet vectors (all k). Then, GK(u,v) = u dot v / ||u|| ||v||
#define SGKDiff 2   // sum of abs(observed-ideal)/ideal ; for all graphlets, all k
#define GraphletGDV 3
#define EdgeHammingDistance 4
#define DegreeDist 5  // euclidean distance between the two degree distribution vectors
#define ClustCoff 6  // binned histogram difference between LOCAL clustering coefficient distributions

static double weights[NUMPROPS] =
// weights: 0 GraphletEuclidean; 1 GraphletKernel; 2 SGKDiff; 3 GDV;  4 EdgeHammingDistance, 5 DegreeDist; 6 ClustCoff
           {0.1,                 0,                0,         0,      0.35,                     0.05,            0.50};
// weights can be set using an env variable
// USAGE: export SYNTHETIC_GRAPHLET_WEIGHTS = 'a b c d e f g'   #a+b+c+d+e+f+g == 1

static double max_abscost[NUMPROPS];

// Here's where we're lazy on saving memory, and we could do better.  We're going to allocate a static array
// that is big enough for the 256 million permutations from non-canonicals to canonicals for k=8, even if k<8.
// So we're allocating 256MBx3=768MB even if we need much less.  I figure anything less than 1GB isn't a big deal
// these days. It needs to be aligned to a page boundary since we're going to mmap the binary file into this array.
static kperm Permutations[maxBk] __attribute__ ((aligned (8192)));
// Here's the actual mapping from non-canonical to canonical, same argument as above wasting memory, and also mmap'd.
// So here we are allocating 256MB x sizeof(short int) = 512MB.
// Grand total statically allocated memory is exactly 1.25GB.
static short int* _K[MAX_K];

// Assuming the global variable _k[] is set properly, go read in and/or mmap the big global
// arrays related to canonical mappings and permutations.
void SetGlobalCanonMaps(void){
    unsigned int _Bk;
    int i;
    for(i=0; i<MAX_K; i++){  // for all values of 'k'
        if (_k[i] == -1)
            break;
        assert(3 <= _k[i] && _k[i] <= 8);
        _Bk = (1 <<(_k[i]*(_k[i]-1)/2));
        char BUF[BUFSIZ];
        _connectedCanonicals[_k[i]-1] = canonListPopulate(BUF, _canonList[_k[i]-1], _k[i]);
        _numCanon[_k[i]-1] = _connectedCanonicals[_k[i]-1]->n;
        _maxNumCanon = MAX(_maxNumCanon, _numCanon[_k[i]-1]);  // set max number of canonicals for a k
        _K[_k[i]-1] = (short int*) aligned_alloc(8192, MAX(_Bk * sizeof(short int), 8192));
        assert(_K[_k[i]-1] != NULL);
        mapCanonMap(BUF, _K[_k[i]-1], _k[i]);
        sprintf(BUF, "%s/%s/perm_map%d.bin", _BLANT_DIR, CANON_DIR, _k[i]);
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

// Incrementally returns new GraphletEuclidean objective cost
// O(1) time
double FastEuclideanObjective(double oldcost, double olddelta, double change){
    // oldcost is the original cost
    // olddelta is the original difference b/w 2 graphlet counts D[1][...canon...] - D[0][...canon...]
    // change = +1 (-1) if count of a graphlet increased (decreased) in the caller code
    double unchanged = (oldcost * oldcost) - (olddelta * olddelta);
    unchanged = MAX(0, unchanged);
    double newcost_sq = unchanged + SQR(olddelta+change);
    newcost_sq = MAX(0, newcost_sq);
    return sqrt(newcost_sq);
}

// Incrementally returns new GraphletGDV objective cost
// O(1) time
double FastGDVObjective(double oldcost, int olddelta, double change){
    // oldcost is the original cost
    // olddelta is the original difference b/w 2 graphlet counts D[1][...canon...] - D[0][...canon...]
    // change = +1 (-1) if count of a graphlet increased (decreased) in the caller code
    double unchanged = (oldcost * oldcost) - (olddelta * olddelta);
    unchanged = MAX(0, unchanged);
    double newcost_sq = unchanged + (double) SQR(olddelta+change);
    newcost_sq = MAX(0, newcost_sq);
    return sqrt(newcost_sq);
}

// Incrementally returns new SGKDiff objective cost
// O(1) time
double FastSGKDiffObjective(double oldcost, int k, int canonNum, int D0, int D1, int change){
    // oldcost is the original cost
    // k and canonNum specify a graphlet
    // D0 and D1 are the graphlet counts - D[0] and D[1]; for a particular k and canonNum
    // change = +1 (-1) if count of a graphlet increased (decreased) in the caller code
    double oldratio, newratio;
    assert(abs(change) == 1);
    int diff;

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

// Updates the GraphletKernel object (GKState), and incrementally returns new GraphletKernel cost
// O(1) time
double AdjustGraphletKernel(int D0, int D1, int change, GKState* gkstate){
    // D0 and D1 are the graphlet counts - D[0] and D[1]; for a particular k and canonNum
    // change = +1 (-1) if count of a graphlet increased (decreased) in the caller code
    // GKState stores 3 integers - (1) u dot v (2) ||u||^2 (3) ||v||^2
    assert(abs(change) == 1);

    // update u dot v
    gkstate->udotv += ((long) D0 * (long) change);  // after simplification :    += [(d1*d0) - ((d1-change)*d0)]

    // update ||v|| ** 2
    gkstate->sq_length_v += (long) ((2*D1*change) - 1);  // after simplification :    += [-(d1-change)**2 + (d1)**2]

    assert(sqrt(gkstate->sq_length_u) >= 0);
    assert(sqrt(gkstate->sq_length_v) >= 0);

    double gk = ((double) gkstate->udotv) / (((double) sqrt(gkstate->sq_length_u)) * ((double) sqrt(gkstate->sq_length_v)));
    double cost = 1-gk;
    assert(ISZERO(cost) || ISZERO(1-cost) || ((cost>0) && (cost<1)));
    return cost;
}

// Updates GDV matrices, and incrementally returns new GraphletGDV cost
double AdjustGDV(int k, int canon, int change, int line[k+1], Dictionary GDVhistograms[2][MAX_K][_maxNumCanon], int GDVbinsize[2][MAX_K][_maxNumCanon], int GDV[2][MAX_K][_maxNumCanon][_numNodes], double oldcost){
    // k and canon identify a graphlet
    // change = +1 (-1) if count of a graphlet increased (decreased) in the caller code
    // line is the BLANT sample line
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

// Updates Degree Distribution arrays, and incrementally returns new DegreeDist cost
// Should be called AFTER calling GraphConnect or GraphDisconnect
// O(1) time
double AdjustDegree(const int x, const int y, const int connected, GRAPH* G, int Degree[2][maxdegree+1], double oldcost){
    // x and y are the 2 nodes which were connected or disconnected
    // connected = +1 (-1) if x and y were connected (disconnected)

    assert(abs(connected) == 1);
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

// Updates Local Clustering CONNECTIONS arrays, and incrementally returns new ClustCoff objective cost
// Should be called AFTER calling GraphConnect or GraphDisconnect
// O(d) time, where d is the number of nodes connected to both x and y
double AdjustClustCoff(const int x, const int y, const int connected, GRAPH* G, int localConnections[2][_numNodes], ClustCoffState* ccstate, double oldcost){
    // x and y are the 2 nodes which were connected or disconnected
    // connected = +1 (-1) if x and y were connected (disconnected)
    // localConnections is an array where lc[i] stores number of edges b/w the neighbors of node `i`
    // ccstate->histograms are binned histograms storing how many nodes lie in a particular interval of local-clustering-coefficient

    /*
    Local clustering coefficient calculation is more precise if number of connections in the neighboorhood
    of a node (integers) are updated, instead of the exact local clustering coefficent (which is a float)
    */

    assert(abs(connected) == 1);
    assert(ccstate->histograms_size[0] == ccstate->histograms_size[1]);
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
                oldkey = (int) (oldcc / ccstate->histograms_bin_size[1]);
                assert(oldkey < ccstate->histograms_size[1]);
                newcc = (double) (c+connected) / (double) nc2;
                newkey = (int) (newcc / ccstate->histograms_bin_size[1]);
                assert(newkey < ccstate->histograms_size[1]);

                newcost = FastEuclideanObjective(newcost, (double) (ccstate->histograms[1][oldkey] - ccstate->histograms[0][oldkey]), (double) -1.0);
                assert((newcost == newcost) && (newcost >= 0));
                assert(ccstate->histograms[1][oldkey] > 0);
                ccstate->histograms[1][oldkey] -= 1;
                newcost = FastEuclideanObjective(newcost, (double) (ccstate->histograms[1][newkey] - ccstate->histograms[0][newkey]), (double) 1.0);
                assert((newcost == newcost) && (newcost >= 0));
                ccstate->histograms[1][newkey] += 1;
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
        oldkey = (int) (oldcc / ccstate->histograms_bin_size[1]);
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
        newkey = (int) (newcc / ccstate->histograms_bin_size[1]);
    }

    localConnections[1][x] = c;
    if(oldkey != newkey){
        newcost = FastEuclideanObjective(newcost, (double) (ccstate->histograms[1][oldkey] - ccstate->histograms[0][oldkey]), (double) -1.0);
        assert((newcost == newcost) && (newcost >= 0));
        assert(ccstate->histograms[1][oldkey] > 0);
        ccstate->histograms[1][oldkey] -= 1;
        newcost = FastEuclideanObjective(newcost, (double) (ccstate->histograms[1][newkey] - ccstate->histograms[0][newkey]), (double) 1.0);
        assert((newcost == newcost) && (newcost >= 0));
        ccstate->histograms[1][newkey] += 1;
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
        oldkey = (int) (oldcc / ccstate->histograms_bin_size[1]);
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
        newkey = (int) (newcc / ccstate->histograms_bin_size[1]);
    }

    localConnections[1][y] = c;
    if(oldkey != newkey){
        newcost = FastEuclideanObjective(newcost, (double) (ccstate->histograms[1][oldkey] - ccstate->histograms[0][oldkey]), (double) -1.0);
        assert((newcost == newcost) && (newcost >= 0));
        assert(ccstate->histograms[1][oldkey] > 0);
        ccstate->histograms[1][oldkey] -= 1;
        newcost = FastEuclideanObjective(newcost, (double) (ccstate->histograms[1][newkey] - ccstate->histograms[0][newkey]), (double) 1.0);
        assert((newcost == newcost) && (newcost >= 0));
        ccstate->histograms[1][newkey] += 1;
    }

    return newcost;
}

// EHD is a matrix which tells the ehd b/w 2 canonicals
// EHDaway is a matrix which tells which canonicals are 'x' distance away from 'y'th canonical
// D[ ] is the raw graphlet count matrix
double EHDObjective(int D[2][MAX_K][_maxNumCanon], int CanonicalEdges[MAX_K][_maxNumCanon], int EHD[MAX_K][_maxNumCanon][_maxNumCanon], int EHDaway[MAX_K][_maxNumCanon][NC2(MAX_K)+1][1 + _maxNumCanon]){
    int i,j,k,l,m,x,diff,temp;
    double sum = 0;

    for(i=0; i<MAX_K; i++){
        k = _k[i];
        if (k == -1)
            break;

        // loop through all canonicals for this k value
        for(j=0; j<_numCanon[k-1]; j++){
            diff = CanonicalEdges[k-1][j] * abs(D[0][k-1][j] - D[1][k-1][j]);
            for(l=1; l<=NC2(k); l++){
                temp = 0;
                // loop through all the canonicals which are 'l' EHD away from jth canonical
                for(m=0; m<EHDaway[k-1][j][l][0]; m++){  // EHDaway[k-1][j][l][0] gives number of canonicals which are 'l' away from 'j'
                    x = EHDaway[k-1][j][l][m+1]; // a canonical which is 'l' EHD edges away from j
                    assert(EHD[k-1][j][x] == l);
                    temp += abs(D[0][k-1][j] - D[1][k-1][x]);
                }
                temp = temp * l;
                assert(temp >= 0);
                diff = MIN(diff, temp);
            }

            // diff for 'j'th canonical is fixed now
            assert(diff >= 0);
            sum += diff;
            assert(sum >= 0);
        }
    }

    return (double) sum;
}

// updates the BLANT sample, and updates the cost for the GRAPHLET BASED objective functions : GraphletEuclidean, GraphletKernel, SGKDiff, GDV
void ReBLANT(int D[2][MAX_K][_maxNumCanon], GKState* gkstate, Dictionary GDVhistograms[2][MAX_K][_maxNumCanon], int GDVbinsize[2][MAX_K][_maxNumCanon], int GDV[2][MAX_K][_maxNumCanon][_numNodes], int CanonicalEdges[MAX_K][_maxNumCanon], int EHD[MAX_K][_maxNumCanon][_maxNumCanon], int EHDaway[MAX_K][_maxNumCanon][NC2(MAX_K)+1][1 + _maxNumCanon], GRAPH *G, SET ***samples, int ***Varrays, int ***BLANT, int v1, int v2, double oldcost[NUMPROPS], RevertStack* rvStack){
    // D stores the squiggly plot vectors
    // gkstate maintains variables for GraphletKernel objective
    // GDV histograms and matrices for GraphletGDV objective
    // SET*** samples stores which lines (samples) is a particular node present in
    // Varrays is the same thing as SET*** samples; used because there is no set iterator
    // BLANT stores the samples [canonNum n1 n2 n3 .... nk]
    // v1 and v2 are the nodes which WERE connected or disconnected, by the caller
    // oldcost contains the current cost of the 4 GRAPHLET BASED objectives
    // rvStack is a stack used to revert changes to the BLANT samples, if a move is rejected (optimization to avoid calling ReBLANT)

    int i, j, line, s, change;
    int canon;
    static TINY_GRAPH *g[MAX_K];
    double oldcanondiff, temp_gdv_newcost;

    for (i=0; i<MAX_K; i++){
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
                if ((!IGNORE_DISCONNECTED_GRAPHLETS) || wasConnected){
                    if (!ISZERO(weights[GraphletEuclidean]))
                        oldcost[0] = FastEuclideanObjective(oldcost[0], oldcanondiff, change);
                    if (!ISZERO(weights[GraphletKernel]))
                        oldcost[1] = AdjustGraphletKernel(D[0][k-1][canon], D[1][k-1][canon], change, gkstate);
                    if (!ISZERO(weights[SGKDiff]))
                        oldcost[2] = FastSGKDiffObjective(oldcost[2], k, canon, D[0][k-1][canon], D[1][k-1][canon], change);
                    if (!ISZERO(weights[GraphletGDV]))
                        oldcost[3] = AdjustGDV(k, canon, change, BLANT[k-1][line], GDVhistograms, GDVbinsize, GDV, oldcost[3]);  // GDV matrices might be out of sync from ideal if GDV weight = 0
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
                if ((!IGNORE_DISCONNECTED_GRAPHLETS) || isConnected){
                    if (!ISZERO(weights[GraphletEuclidean]))
                        oldcost[0] = FastEuclideanObjective(oldcost[0], oldcanondiff, change);
                    if (!ISZERO(weights[GraphletKernel]))
                        oldcost[1] = AdjustGraphletKernel(D[0][k-1][canon], D[1][k-1][canon], change, gkstate);
                    if (!ISZERO(weights[SGKDiff]))
                        oldcost[2] = FastSGKDiffObjective(oldcost[2], k, canon, D[0][k-1][canon], D[1][k-1][canon], change);
                    if (!ISZERO(weights[GraphletGDV]))
                        oldcost[3] = AdjustGDV(k, canon, change, BLANT[k-1][line], GDVhistograms, GDVbinsize, GDV, oldcost[3]);
                }


                // change object
                newchange.new = (int) canon;
                assert(push(rvStack, newchange) == 0);
            }

    }

    // recomputed from scratch, can be outside 'iterate over k' loop
    if (!ISZERO(weights[EdgeHammingDistance]))
        oldcost[4] = EHDObjective(D, CanonicalEdges, EHD, EHDaway);



    for(i=0; i<5; i++)
        assert(oldcost[i] >= 0);
}

// restores the BLANT samples, given a stack (rvStack)
void Revert(int ***BLANT, int D[2][MAX_K][_maxNumCanon], Dictionary GDVhistograms[2][MAX_K][_maxNumCanon], int GDVbinsize[2][MAX_K][_maxNumCanon], int GDV[2][MAX_K][_maxNumCanon][_numNodes], RevertStack* rvStack){
    // restore the BLANT line
    // restore D vectors AND _numDisconnectedGraphlets
    // restore GDV[1] matrix and the GDVhistograms

    /*
    a stack entry contains the BLANT line-number, k, original and new graphlet number
    so, count of 'new' gets decremented and 'original' gets incremented (in D)
    */

    Change change;
    int n, node, key, value, b;
    int* line;

    while (pop(rvStack, &change) == 0){
        line = BLANT[change.k-1][change.linenum];

        Boolean wasConnected = SetIn(_connectedCanonicals[change.k-1], BLANT[change.k-1][change.linenum][0]);
        if ((!IGNORE_DISCONNECTED_GRAPHLETS) || wasConnected){
            if (!ISZERO(weights[GraphletGDV]))
                AdjustGDV(change.k, change.new, -1, line, GDVhistograms, GDVbinsize, GDV, 1000);  // the 1000 has no meaning
        }

        Boolean  isConnected = SetIn(_connectedCanonicals[change.k-1], change.original);
        BLANT[change.k-1][change.linenum][0] = change.original;
        if ((!IGNORE_DISCONNECTED_GRAPHLETS) || isConnected){
            if (!ISZERO(weights[GraphletGDV]))
                AdjustGDV(change.k, change.original, 1, line, GDVhistograms, GDVbinsize, GDV, 1000);
        }

        if(wasConnected && !isConnected)
            ++_numDisconnectedGraphlets;
        if(!wasConnected && isConnected)
            --_numDisconnectedGraphlets;
        --D[1][change.k-1][change.new];
        ++D[1][change.k-1][change.original];
    }
}

// global objective function
// returns WEIGHTED and NORMALIZED total costs, given current absolute costs
double Objective(double abscost[NUMPROPS]){
    double cost = 0;
    int i = 0;

    for(i=0; i<NUMPROPS; i++){
        assert(abscost[i] >= 0);
        assert(max_abscost[i] >= 0);
        if (max_abscost[i] == 0)
            cost += ((double) weights[i] * (abscost[i] / 1));
        else
            cost += ((double) weights[i] * (abscost[i] / max_abscost[i]));
    }

    assert(cost >= 0);
    return cost;
}

// graphlet euclidean cost
// slow, used once to initialize. O(n) where n is num of canonicals for all k
double GraphletEuclideanObjective(int D[2][MAX_K][_maxNumCanon]){
    // returns ABSOLUTE cost
    int i,j;
    double logP = 0, sum2 = 0;

    for (i=0; i<MAX_K; i++){
        if (_k[i] == -1) break;
        for (j=0; j<_numCanon[_k[i]-1]; j++){
        /*
        double pd = PoissonDistribution(D[0][_k[i]-1][j], D[1][_k[i]-1][j]);
        if(pd > 1)
        Fatal("umm.... PoissonDistribution returned a number greater than 1");
        if(pd>0)
        logP += log(pd); // if we're close, use probability
        */
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

// graphlet kernel cost; GK(u, v) = (u dot v) / ||u|| ||v||; where u and v are two graphlet count vectors (single big vector for all k)
// As GK is a similarity measure, 1-GK is used as the objective function COST
// slow, used once to initialize. O(n) where n is num of canonicals for all k
double GraphletKernelObjective(const int D[2][MAX_K][_maxNumCanon], GKState* gkstate){
    // GKState stores 3 integers - (1) u dot v (2) ||u||^2 (3) ||v||^2
    // calculations are precise if we maintain 3 integers, rather than the exact GK value (which is a float)
    int i,j;
    gkstate->udotv = (long) 0;
    gkstate->sq_length_u = (long) 0;
    gkstate->sq_length_v = (long) 0;

    for(i=0; i<MAX_K; i++){
        int k = _k[i];
        if (k == -1)
            break;

        for(j=0; j<_numCanon[k-1]; j++){
            if ((!IGNORE_DISCONNECTED_GRAPHLETS) || (SetIn(_connectedCanonicals[k-1], j))){
                gkstate->sq_length_u += SQR((long) D[0][k-1][j]);
                assert(gkstate->sq_length_u >= 0);
                gkstate->sq_length_v += SQR((long) D[1][k-1][j]);
                assert(gkstate->sq_length_v >= (long) 0);
                gkstate->udotv += ((long) D[0][k-1][j] * (long) D[1][k-1][j]);
                assert(gkstate->udotv >= 0);
            }
        }
    }

    assert(sqrt(gkstate->sq_length_u) > 0);
    assert(sqrt(gkstate->sq_length_v) > 0);
    double gk = ((double) gkstate->udotv) / (((double) sqrt(gkstate->sq_length_u)) * ((double) sqrt(gkstate->sq_length_v)));
    double cost = 1-gk;
    assert(ISZERO(cost) || ISZERO(1-cost) || ((cost>0) && (cost<1)));
    return cost;
}

// Senatorial Graphlet Difference; SGK = sum of abs(observed-real)/real for all canonicals (all k)
// slow, used once to initialize. O(n) where n is num of canonicals for all k
double SGKDiffObjective(int D[2][MAX_K][_maxNumCanon]){
    int i,j;
    double sum = 0;

    for(i=0; i<MAX_K; i++){
        int k = _k[i];
        if (k == -1)
            break;
        for (j=0; j<_numCanon[k-1]; j++){
            int diff = abs(D[0][k-1][j] - D[1][k-1][j]);
            int target = D[0][k-1][j];
            if ((!IGNORE_DISCONNECTED_GRAPHLETS) || (SetIn(_connectedCanonicals[k-1], j))){
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

// computes the difference between the GDV distribution histograms (traget vs synthetic)
// these histograms are initialized/populated in the main driver code
// slow, used once to initialize
double GDVObjective(Dictionary GDVhistograms[2][MAX_K][_maxNumCanon]){
    double sum = 0;
    int j,k,canon;
    int key_tar, val_tar, key_syn, val_syn;

    Dictionary hist_tar, hist_syn;
    KeyValue *iter_tar, *iter_syn;

    for(j=0; j<MAX_K; j++){
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

// euclidean distance between two degree distribution vectors
// slow, used once to initialize
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

// computes the number of connections in the neighborhood (used for local-clustering coefficient)
// slow, used once to initialize
void GetConnections(GRAPH *G, int localConnections[G->n]){
    // localConnections[i] will contain the number of edges b/w the neighbors of node `i`
    // used for local-clustering-coeffienct computation later on
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

// computes euclidean distance between the local-clustering-coefficient histograms
double ClustCoffObjective(const ClustCoffState* ccstate){
    double cost = 0;
    double temp;
    int j;

    assert(ccstate->histograms_size[0] == ccstate->histograms_size[1]);

    for(j=0; j < ccstate->histograms_size[0]; j++){
        temp = (double) SQR(ccstate->histograms[0][j] - ccstate->histograms[1][j]);
        assert(temp >= 0);
        cost += temp;
    }

    double returnVal = sqrt(cost);
    assert(returnVal == returnVal);
    assert(returnVal >= 0);
    return returnVal;
}

// node selection (to connect or disconnect): as per `static int node_selection`
void GetNodes(GRAPH* G, const SmallWorld sw, int nodesBySp[G->n], int index_in_nodesBySp[G->n], int* u1, int* u2, int* v1, int* v2){
    // sw specifies if the network should be made more/less small-worldy (or nothing) --> read the definition in syntheticDS.h
    // nodesBySp is an array of nodes, where nodes at lower indexes have smaller number of shortest paths going through them, as compared to nodes at higher indexes
    // index_in_nodesBySp is a lookup table for nodesBySp. So, index_in_nodesBySp[i] provides the index of node `i` in nodesBySp
    // u1 u2 v1 v2 are the RETURN values from this function (nodes v1,v2 should be disconnected) (nodes u1,u2 should be connected)

    /*
    Case 1: Random selection explicitly set OR k-hop distributions match (medians)
    $>: randomly select edge and non-edge

    Case 2: k-hop distributions don't match (medians)
    todo: make the synthetic more/less like a small-world network (this tunes - betweenness, eccentricity, diameter, and khop)
    How to do this?
        Method 1: By BFS hops (select nodes to connect & disconnect that are at a distance 'd' away from eachother)
        Method 2: Use nodesBySp array to connect/disconnect nodes

    Method 1 -
    If you need to make the synthetic MORE small-world; connect nodes which are further away, disconnect any random edge
    If you need to make the synthetic LESS small-world; connect nodes which are closer, disconnect any random edge
    (these conditions are set on line 1358)

    $>:
        get a random node1
        get a random number of hops from a range (ex [2,5]), ex. 3
        do a bfs starting from the node, and get a random node2 at that distance from node1
        connect node1-node2
        disconnect any random edge

    Method 2 -
    If you need to make the synthetic MORE small-world; connect nodes with lesser num. of SPs through them, disconnect nodes with higher num. of SPs
    If you need to make the synthetic LESS small-world; connect nodes with higher num. of SPs through them, disconnect nodes with lower num. of SPs
    (these conditions are set on line 1358)

    $>:
        take random n1 & n2 in a particular interval (lower 35% or higher 35% indexes in nodesBySp)
        connect n1 & n2
        do:
            take random n3 in a particular interval (lower 35% or higher 35% indexes in nodesBySp)
            get random n4 connected to n3
        while:
            n4 not in same interval as n3
        disconnect n3 & n4

    */


    int x, y;

    if ((node_selection==NODE_SEL_ALWAYS_RANDOM) || (sw.make == 0)){
        // random selection

        // existing edge
        do {
            x = drand48()*G->n;
            y = getRandomConnectedNode(G, x);
        } while((G->degree[x] <= 1) || (G->degree[y] <= 1));
        assert(GraphAreConnected(G, x, y));
        assert(x!=y);
        memcpy(v1, &x, sizeof(int));
        memcpy(v2, &y, sizeof(int));

        // existing non-edge
        do {
            x = drand48()*G->n;
            y = drand48()*G->n;
        }while((x==y) || (GraphAreConnected(G, x, y)));
        assert(!GraphAreConnected(G, x, y));
        assert(x!=y);
        memcpy(u1, &x, sizeof(int));
        memcpy(u2, &y, sizeof(int));
        return;
    }

    // for node selection by BFS hops
    if(node_selection == NODE_SEL_BY_HOPS){
        int d = sw.lowerHops + (int)((sw.upperHops - sw.lowerHops + 1) * drand48());
        assert(d>1);

        // existing edge
        do{
            x = drand48()*G->n;
            y = getRandomConnectedNode(G, x);
        } while((G->degree[x] <= 1) || (G->degree[y] <= 1));
        assert(GraphAreConnected(G, x, y));
        assert(x!=y);
        memcpy(v1, &x, sizeof(int));
        memcpy(v2, &y, sizeof(int));

        // existing non-edge
        do {
            x = drand48()*G->n;
            y = getRandomNodeAtHops(G, x, d);
        }while((x==y) || (GraphAreConnected(G, x, y)));
        assert(x!=y);
        assert(!GraphAreConnected(G, x, y));
        memcpy(u1, &x, sizeof(int));
        memcpy(u2, &y, sizeof(int));
        return;
    }

    // for node selection by node ShortestPaths
    if(node_selection == NODE_SEL_SHORT_PATH){
        double lower = 0.0;
        double pd_interval = 0.05;
        double interval = 0.35;
        int random;
        int l, u, yi;

        // existing edge
        do{
            random = drand48() * interval * G->n;
            if (sw.make == -1){  // make network LESS small-world
                // disconnect a node with HIGH SP
                x = nodesBySp[(G->n -1) - (int)(MAX(lower, pd_interval) * G->n) - random];
                l = (G->n -1) - (int)(MAX(lower, pd_interval) * G->n) - (interval * G->n);
                u = G->n - 1;
            }else if(sw.make == 1){  // make network MORE small-world
                // disconnect a node with LOW SP
                x = nodesBySp[0 + (int)(lower*G->n) + random];
                l = 0;
                u = (int)(lower*G->n) + (interval * G->n);
            }
            y = getRandomConnectedNode(G, x);
            yi = index_in_nodesBySp[y];
            assert((yi>=0) && (yi<G->n));
        }while((G->degree[x] <= 1) || (G->degree[y] <= 1) || (yi<l) || (yi>u));
        assert(GraphAreConnected(G, x, y));
        assert(x!=y);
        memcpy(v1, &x, sizeof(int));
        memcpy(v2, &y, sizeof(int));

        // existing non-edge
        do{
            if(sw.make == -1){  // make network LESS small world
                // connect nodes with LOW SP
                random = drand48() * interval * G->n;
                x = nodesBySp[0 + (int)(lower*G->n) + random];
                random = drand48() * interval * G->n;
                y = nodesBySp[0 + (int)(lower*G->n) + random];
            }else if(sw.make == 1){  // make network MORE small world
                // connect nodes with HIGH SP
                random = drand48() * interval * G->n;
                x = nodesBySp[(G->n -1) - (int)(lower*G->n) - random];
                random = drand48() * interval * G->n;
                y = nodesBySp[(G->n -1) - (int)(lower*G->n) - random];
            }
        }while((x==y) || (GraphAreConnected(G, x, y)));
        assert(!GraphAreConnected(G, x, y));
        assert(x!=y);
        memcpy(u1, &x, sizeof(int));
        memcpy(u2, &y, sizeof(int));
        return;
    }
}

int main(int argc, char *argv[]){
    srand48(time(0)+getpid());
    int i, opt, j, line;

    if(argc==1){
    fprintf(stderr, "%s\n", USAGE);
    exit(1);
    }

    // initialize _k[]
    for(i=0; i<MAX_K; i++)
        _k[i] = -1;

    // Read objective function weights
    char* weightString = getenv("SYNTHETIC_GRAPHLET_WEIGHTS");
    if(weightString){
        i=0;
        int n = strlen(weightString);
        fprintf(stderr, "reading weights: %s\n", weightString);
        int l;
        for(j=0; j < n; j++){
            if(weightString[j] == ' ')
                continue;
            i += 1;
            l = j;
            while((l<n) && (weightString[l] != ' '))
                l += 1;
            if (weightString[l] == ' '){
                weightString[l] = '\0';
                weights[i-1] = atof(weightString + j);
                weightString[l] = ' ';
            }else
                weights[i-1] = atof(weightString + j);
            j = l;
        }
        assert(i == NUMPROPS);
    }
    {
        double wsum = 0;  // sum of weights should be 1.0
    fprintf(stderr, "weights:");
        for(i=0; i<NUMPROPS; i++){
        fprintf(stderr, " %g",weights[i]);
        wsum += weights[i];
    }
    fprintf(stderr,"\n");
        assert(fabs(wsum-1) < 0.0001);
    }


    //Read node selection strategy, 0 for random | 1 for bfs-hops | 2 for nodes by ShortestPaths
    char* nselect = getenv("SYNTHETIC_NODE_SELECTION");
    if(nselect)
        node_selection = atoi(nselect);
    assert((node_selection>=0) && (node_selection<=2));


    //Read max value of k and stagnation
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

    // read edge-lists
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
    assert(G[0]->n == G[1]->n);
    _numNodes = G[0]->n;

    // The distribution of graphlets (squiggly plot vectors) (initialization)
    assert(_maxNumCanon != -1);  // this should be set >0 by calling SetGlobalCanonMaps() first
    int D[2][MAX_K][_maxNumCanon];
    for(i=0; i<MAX_K;i++)
        for (j=0; j<_maxNumCanon; j++)
            D[0][i][j] = D[1][i][j] = 0;

    // GraphletKernel (initialization)
    GKState gkstate, newGkstate;
    gkstate.udotv = gkstate.sq_length_u = gkstate.sq_length_v = 0;

    // GDV matrices (initialization)
    /* The GDV is a matrix which has n rows and l columns, n = number of nodes in the graph
    and l = number of canonical graphlets for the size of k-graphlets you are working with.
    I can consider each column j of the GDV as the graphlet degree distribution of graphlet j in the graph.
    To get the degree distribution, you'd have to "summarize" the column as follows:
    Compute S(t)_j = The number of nodes appearing in a graphlet of type j, t times.
    In other words, in that column, how many times did t appear for t = 0, 1, 2, .... infinity
    (it won't actually be infinity, but it can get very large).*/
    int GDV[2][MAX_K][_maxNumCanon][_numNodes];  // 4 dimensional
    for(i=0; i<2; i++){
        for(j=0; j<MAX_K; j++){
            if (_k[j] == -1) break;
            int l;
            for (l=0; l<_numCanon[_k[j]-1]; l++){
                int m;
                for (m=0; m < G[i]->n; m++)
                    GDV[i][_k[j]-1][l][m] = 0;  // initialize to 0
            }
        }
    }

    // READ blant into squilly plot vectors and GDV matrices
    // expect 2 blant files (target & synthetic for every _k value)
    // assume all blant files have same number of samples = _numSamples
    int **BLANT[2][MAX_K];
    for(i=0;i<2;i++){
        for(j=0; j<MAX_K; j++){
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
                BLANT[i][_k[j]-1][line] = (int*) Malloc((MAX_K+1) * sizeof(int));
                int l;
                for (l=0; l<=_k[j]; l++){
                    assert(1 == fscanf(fp, "%d", &(BLANT[i][_k[j]-1][line][l])));
                    if (l>0){
                        GDV[i][_k[j]-1][BLANT[i][_k[j]-1][line][0]][BLANT[i][_k[j]-1][line][l]] += 1; // update GDV matrix
                    }
                }
                assert(BLANT[i][_k[j]-1][line][0] < _maxNumCanon);
                ++D[i][_k[j]-1][BLANT[i][_k[j]-1][line][0]]; // update squiggly plot vector
            }
            assert(fscanf(fp, "%d", &line) < 1); // ensure there's nothing left to read.
            fclose(fp);
            optind++;
        }

    }

    // sanity check - squiggly vectors
    for(i=0; i<2; i++){
        for (j=0; j<MAX_K; j++){
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
        for(j=0; j<MAX_K; j++){
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

    // Create GDV histograms from GDV matrices
    // Bin size : getIntegerBinSize in syntheticDS.c (Freedman-Diaconis’s Rule)
    Dictionary GDVhistograms[2][MAX_K][_maxNumCanon];  // GDVhistograms, derived from GDV
    int GDVbinsize[2][MAX_K][_maxNumCanon];  // // Bin-size for these GDVhistograms
    int* scratchspace = (int*) malloc(_numNodes * sizeof(int));  // used for sorting the GDV column
    for(i=0; i<2; i++){
        for(j=0; j<MAX_K; j++){
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

    // sanity check GDV bin size
    for (i=0; i<2; i++){
        for(j=0; j<MAX_K; j++){
            if (_k[j] == -1)
                break;
            int l=0;
            for(l=0; l<_numCanon[_k[j]-1]; l++){
                assert(GDVbinsize[i][_k[j]-1][l] > 0);
            }
        }
    }

    // store what samples is a particular node part of
    SET **samples[MAX_K];
    for(i=0; i<MAX_K; i++){
        if (_k[i] == -1) break;
        samples[_k[i]-1] = (SET**) Malloc(G[1]->n * sizeof(SET*));
    }

    for(i=0; i<MAX_K; i++){
        if (_k[i] == -1) break;

        for (j=0; j<G[1]->n; j++)
            samples[_k[i]-1][j] = SetAlloc(_numSamples);

        for (line=0; line<_numSamples; line++)
            for (j=1; j<= _k[i]; j++)
                SetAdd(samples[_k[i]-1][BLANT[1][_k[i]-1][line][j]], line);    // SetAdd(samples[k][nodenum], line);
    }

    // Varrays is the same as SET samples. It is used as an iterator
    int **Varrays[MAX_K];
    for (i=0; i<MAX_K; i++){
        if (_k[i] == -1) break;
        Varrays[_k[i]-1] = (int**) Malloc(G[1]->n * sizeof(int*));

        for (j=0; j < G[1]->n; j++){
            Varrays[_k[i]-1][j] =  (int*) Malloc((1+SetCardinality(samples[_k[i]-1][j]))* sizeof(int));
            Varrays[_k[i]-1][j][0] = SetToArray(Varrays[_k[i]-1][j]+1, samples[_k[i]-1][j]);
        }
    }

    // Edge-Hamming-Distances
    int CanonicalEdges[MAX_K][_maxNumCanon];  // given a canonical, get num edges
    int EHD[MAX_K][_maxNumCanon][_maxNumCanon];  // get EHD b/w two canonicals
    int EHDaway[MAX_K][_maxNumCanon][NC2(MAX_K)+1][1 + _maxNumCanon];  // given a canonical, get all canonicals, 'x' EHD away

    // 0. Read canon_map.txt files to get canonical vs edges
    // 1. Read the ehdk.txt files
    // 2. Populate EHD[MAX_K][__maxNumCanon][_maxNumCanon]
    // 3. Populate EHDaway[MAX_K][_maxNumCanon][NC2(MAX_K)+1][1 + _maxNumCanon]

    for(i=0; i<MAX_K; i++){
        int k = _k[i];
        if(k == -1) break;
        char FILENAME[100];
        sprintf(FILENAME, "%s/%s/canon_list%d.txt", _BLANT_DIR, CANON_DIR, k);
        FILE* fp = fopen(FILENAME, "r");
        assert(fp);

        int j, c, e;

        // 3 columns
        fscanf(fp, "%d", &e);
        assert(e == _numCanon[k-1]);

        for(c=0; c<_numCanon[k-1]; c++){
            for(j=0; j<3; j++){
                fscanf(fp, "%d", &e);
                if (j == 2){
                    CanonicalEdges[k-1][c] = e;
                }
            }
        }

        // check nothing left to read
        assert(fscanf(fp, "%d", &e) < 1);
        fclose(fp);
    }

    for (i=0; i<MAX_K; i++){
        if(_k[i] == -1) break;
        int k = _k[i];
        char FILENAME[100];
        sprintf(FILENAME, "%s/%s/EdgeHammingDistance%d.txt", _BLANT_DIR, CANON_DIR, k);
        FILE* fp = fopen(FILENAME, "r");
        if(!fp) Fatal("cannot open file %s\n", FILENAME);
        int c1,c2, x,y,z;

        // 3 columns
        for(c1=0; c1<_numCanon[k-1]; c1++){
            for(c2=0; c2<_numCanon[k-1]; c2++){
                fscanf(fp, "%d", &x);
                fscanf(fp, "%d", &y);
                fscanf(fp, "%d", &z);
                assert(x == c1);
                assert(y == c2);
                assert(z >= 0);
                if (x==y) assert(z == 0);
                EHD[k-1][c1][c2] = z;
            }
        }

        // check nothing left to read
        assert(fscanf(fp, "%d", &x) < 1);
        fclose(fp);
    }

    for (i=0; i<MAX_K; i++){
        if(_k[i] == -1) break;
        int k = _k[i];
        int c1,c2,d,index;

        for(c1=0; c1<_numCanon[k-1]; c1++){
            for(j=0; j<NC2(k); j++)
                EHDaway[k-1][c1][j][0] = 0;  // initialize the count to 0
            for(c2=0; c2<_numCanon[k-1]; c2++){
                d = EHD[k-1][c1][c2];
                EHDaway[k-1][c1][d][0] += 1;
                index = EHDaway[k-1][c1][d][0];
                EHDaway[k-1][c1][d][index] = c2;
            }
        }
    }

    // degree distribution vectors
    maxdegree = MAX(G[0]->n, G[1]->n);
    for (i=0; i<2; i++)
        for (j=0; j < G[i]->n; j++)
            assert((G[i]->degree[j]) < (G[i]->n));

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
    int localConnections[2][_numNodes];  // for every node, number of connections in it's neighborhood
    for (i=0; i<2; i++)
        GetConnections(G[i], localConnections[i]);   // populates the array with every nodes' local CONNECTIONS (not clustCoff)

    ClustCoffState ccstate;  // CC Histograms

    // popluate the local-clust-coff histograms
    for(i=0; i<2; i++){
        int degree, nc2, key, value;
        double local_cc, b, localClustCoff[G[i]->n], scratchspace[G[i]->n];

        for(j=0; j < G[i]->n; j++){
            degree = G[i]->degree[j];
            nc2 = (degree * (degree-1))/2;
            if (nc2 <= 0) local_cc = 0;
            else local_cc = (double) localConnections[i][j] / (double) nc2;
            localClustCoff[j] = local_cc;
            assert((localClustCoff[j] >= 0) && (localClustCoff[j] <= 1));
        }

        if(i == 0)
            ccstate.histograms_bin_size[i] = getDoubleBinSize(G[i]->n, localClustCoff, scratchspace);
        else
            ccstate.histograms_bin_size[i] = ccstate.histograms_bin_size[0];

        b = ccstate.histograms_bin_size[i];
        assert(b > 0);

        ccstate.histograms_size[i] = ((int) (1/b)) + 1;
        ccstate.histograms[i] = (int*) Malloc(ccstate.histograms_size[i] * sizeof(int));
        for(j=0; j < ccstate.histograms_size[i]; j++)
            ccstate.histograms[i][j] = 0;

        for(j=0; j < G[i]->n; j++){
            key = (int) (localClustCoff[j]/b);
            assert(key < ccstate.histograms_size[i]);
            ccstate.histograms[i][key] += 1;
        }
    }

    // initialize KHOP distribution
    Dictionary khop[2];  // khop distributions
    int nodesBySp[2][_numNodes];  // an array which contains nodes which have less SPs going through them at lower indexes (vice-versa)
    int index_in_nodesBySp[2][_numNodes];  // index `i` gives the position of vertex `i` in nodesBySp
    sampleKHop(G[0], &(khop[0]), KHOP_QUALITY, nodesBySp[0]);
    for(i=0; i<G[0]->n; i++)
        index_in_nodesBySp[0][nodesBySp[0][i]] = i;

    // initialize node-selection strategy (more/less small-worldy)
    SmallWorld sw;
    sw.khop_interval = DEFAULT_KHOP_INTERVAL;
    sw.make = 0;
    sw.lowerHops = sw.upperHops = -1;

    RevertStack rvStack;
    create_stack(&rvStack, 2 * MAX_K * _numSamples);

    max_abscost[GraphletEuclidean] = GraphletEuclideanObjective(D);
    max_abscost[GraphletKernel] = GraphletKernelObjective(D, &gkstate);
    max_abscost[SGKDiff] = SGKDiffObjective(D);
    max_abscost[GraphletGDV] = GDVObjective(GDVhistograms);
    max_abscost[EdgeHammingDistance] = EHDObjective(D, CanonicalEdges, EHD, EHDaway);
    max_abscost[DegreeDist] = DegreeDistObjective(Degree);
    max_abscost[ClustCoff] = ClustCoffObjective(&ccstate);

    double abscost[NUMPROPS];
    memcpy(abscost, max_abscost, NUMPROPS * sizeof(double));

    fprintf(stderr, "STAG=%d\n", _stagnated);
    fprintf(stderr, "BLANT samples=%d\n", _numSamples);
    fprintf(stderr, "Starting ABSOLUTE costs: GraphletEuclidean: %g, GraphletKernel: %g, GraphDiff: %g, GDV: %g, EHD: %g, DegreeDist: %g, ClustCoff: %g\n", abscost[0], abscost[1], abscost[2], abscost[3], abscost[4], abscost[5], abscost[6]);

    double cost = Objective(abscost), startCost = cost, newCost, maxCost = cost;  // evaluate Objective() once, at the start.
    assert(cost == cost);

    long int sa_iter = 0;
    long int a_iter = 0;  // accepted moves
    double pBad, unif_random;
    float temperature;  // it's okay to overflow and become 0

    double newcost[NUMPROPS];

    /* MAIN LOOP */
    // while(not done---either some number of iterations, or objective function says we're too far away)
    while(cost/startCost > 0.45) // that's enough progress otherwise we're over-optimizing at this sample size
    {

    int u1, u2, v1, v2;

    // compute k-hop -> make synthetic more/less like a small-world
    if((node_selection != NODE_SEL_ALWAYS_RANDOM) && ((a_iter % sw.khop_interval) == 0)){
        sampleKHop(G[1], &(khop[1]), KHOP_QUALITY, nodesBySp[1]);

        // reverse lookup
        for(i=0; i<G[1]->n; i++)
        index_in_nodesBySp[1][i] = -1;

        for(i=0; i<G[1]->n; i++){
        assert(index_in_nodesBySp[1][nodesBySp[1][i]] == -1);
        index_in_nodesBySp[1][nodesBySp[1][i]] = i;
        }

        int medians[2];
        int MAX_Keys[2];

        int m = compareKHopByMedian(khop, medians, MAX_Keys);
        assert(abs(m) <= 1);
        if(m==0){
        // make synthetic LESS like a small-world
        sw.make = -1;
        sw.khop_interval = FIX_KHOP_INTERVAL;
        sw.lowerHops = 2;
        sw.upperHops = medians[1];
        }else if(m==1){
        // make synthetic MORE like a small-world
        sw.make = 1;
        sw.khop_interval = FIX_KHOP_INTERVAL;
        sw.lowerHops = medians[1];
        sw.upperHops = MAX(MAX_Keys[0], MAX_Keys[1]);
        }else if(m==-1){
        // random node selection
        sw.make = 0;
        sw.khop_interval = DEFAULT_KHOP_INTERVAL;
        }
    }

    // get nodes to connect or disconnect (acoording to selection strategy)
    GetNodes(G[1], sw, nodesBySp[1], index_in_nodesBySp[1], &u1, &u2, &v1, &v2);

    // initialize the stack
    assert(init_stack(&rvStack) == 0);

    // copy current cost into newcost, then accept/reject move based on newcost
    memcpy(newcost, abscost, NUMPROPS * sizeof(double));

    // graphletKernelObjective
    memcpy(&newGkstate, &gkstate, sizeof(GKState));

    // Remove an edge
    GraphDisconnect(G[1], v1, v2);
    if (!ISZERO(weights[DegreeDist]))
        newcost[DegreeDist] = AdjustDegree(v1, v2, -1, G[1], Degree, newcost[DegreeDist]);
    if (!ISZERO(weights[ClustCoff]))
        newcost[ClustCoff] = AdjustClustCoff(v1, v2, -1, G[1], localConnections, &ccstate, newcost[ClustCoff]);
    ReBLANT(D, &newGkstate, GDVhistograms, GDVbinsize, GDV, CanonicalEdges, EHD, EHDaway, G[1], samples, Varrays, BLANT[1], v1, v2, newcost, &rvStack);

    // Add an edge
    GraphConnect(G[1], u1, u2);
    if (!ISZERO(weights[DegreeDist]))
        newcost[DegreeDist] = AdjustDegree(u1, u2, 1, G[1], Degree, newcost[DegreeDist]);
    if (!ISZERO(weights[ClustCoff]))
        newcost[ClustCoff] = AdjustClustCoff(u1, u2, 1, G[1], localConnections, &ccstate, newcost[ClustCoff]);
    ReBLANT(D, &newGkstate, GDVhistograms, GDVbinsize, GDV, CanonicalEdges, EHD, EHDaway, G[1], samples, Varrays, BLANT[1], u1, u2, newcost, &rvStack);

    newCost = Objective(newcost);
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
        memcpy(abscost, newcost, NUMPROPS * sizeof(double));

        // update max costs
        for(i=0; i<NUMPROPS; i++)
        max_abscost[i] = MAX(max_abscost[i], abscost[i]);

        // graphletKernel
        memcpy(&gkstate, &newGkstate, sizeof(GKState));

        same = 0;
        ++a_iter;
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
        if (!ISZERO(weights[DegreeDist]))
        newcost[DegreeDist] = AdjustDegree(u1, u2, -1, G[1], Degree, newcost[DegreeDist]);
        if (!ISZERO(weights[ClustCoff]))
        newcost[ClustCoff] = AdjustClustCoff(u1, u2, -1, G[1], localConnections, &ccstate, newcost[ClustCoff]);

        GraphConnect(G[1], v1, v2);
        if (!ISZERO(weights[DegreeDist]))
        newcost[DegreeDist] = AdjustDegree(v1, v2, 1, G[1], Degree, newcost[DegreeDist]);
        if (!ISZERO(weights[ClustCoff]))
        newcost[ClustCoff] = AdjustClustCoff(v1, v2, 1, G[1], localConnections, &ccstate, newcost[ClustCoff]);

        // revert changes
        Revert(BLANT[1], D, GDVhistograms, GDVbinsize, GDV, &rvStack);
    }

    if(same > _stagnated || _numDisconnectedGraphlets >= _numSamples*10){
        fprintf(stderr, "stagnated!, total-iterations=%d, accepted-iterations=%d\n", sa_iter, a_iter);
        break;
    }
    ++sa_iter;
    }
    fprintf(stderr,"\n");
    for(i=0; i < G[1]->numEdges; i++)
        printf("%d %d\n", G[1]->edgeList[2*i], G[1]->edgeList[2*i+1]);
}
