//#define NDEBUG 1
//#define PARANOID_ASSERTS 0
#include "misc.h"
#include "graph.h"
#include "sets.h"
#include "rand48.h"
#include <math.h>
#include <signal.h>

typedef struct _clustering {
    GRAPH *G;
    SET **clusters;
    unsigned nC, *clusSize;
} CLUSTERING;

// Allocate a clustering and return it with the zero'th cluster assigned to all nodes in G.
CLUSTERING *ClusteringAlloc(GRAPH *G) {
    CLUSTERING *C = Calloc(sizeof(CLUSTERING), G->n);
    C->G = G;
    C->nC=0;
    C->clusters = Calloc(G->n, sizeof(SET*));
    C->clusters[0] = SetAlloc(G->n);
    for(int i=0;i < G->n; i++) SetAdd(C->clusters[0], i);
    C->clusSize = Calloc(G->n, sizeof(unsigned));
    C->clusSize[0] = G->n;
    ++C->nC;
    return C;
}

double ScoreOneCluster(GRAPH *G, SET *c) {
    unsigned nodes[G->n], n = SetToArray(nodes, c), m=0;
    assert(n>0);
    if(n==1) return 0;
    if(G->sparse) {
	for(int i=0; i<n; i++) {
	    int u = nodes[i];
	    for(int j=0; j < G->degree[u]; j++) if(SetIn(c, G->neighbor[u][j])) ++m;
	}
    } else {
	for(int i=0; i<n; i++) for(int j=i+1;j<n;j++) if(GraphAreConnected(G,nodes[i],nodes[j])) ++m;
    }
    double ED= m/(n*(n-1)/2.0);
    return m*ED;
}

static int _largestCluster;

double ScoreClustering(CLUSTERING *C) {
    double score=0;
    _largestCluster=0;
    for(int c=0; C->clusters[c]; c++) {
	assert(C->clusSize[c] > 0);
	assert(C->clusSize[c] < C->G->n);
	if(C->clusSize[c] > _largestCluster) _largestCluster = C->clusSize[c];
	double CS = ScoreOneCluster(C->G, C->clusters[c]); // printf(" %g", CS);
	score += CS;
    }
    return score;
}


// Remove some random set of elements from C[who], put them in C->clusters[nC], and increment C->nC
void SplitCluster(CLUSTERING *C, unsigned who) {
    GRAPH *G=C->G;
    assert(who < G->n);
    unsigned nodes[G->n], n = SetToArray(nodes, C->clusters[who]);
    assert(n>1);
    assert(C->clusters[C->nC]==NULL);
    C->clusters[C->nC]=SetAlloc(G->n);
    int which, num2move = 1+drand48()*(n-1)/2; // remove at least 1, and at most half FIXME: might this introduce a bias?
    // int which, num2move = 1+drand48()*(n-1);   // remove at least 1, and at most (n-1) FIXME: might this introduce a bias?
    for(int i=0; i<num2move; i++) {
	until(SetIn(C->clusters[who], nodes[(which = drand48()*n)])) ;
	SetDelete(C->clusters[who], nodes[which]) ;
	SetAdd(C->clusters[C->nC], nodes[which]);
    }
    C->clusSize[C->nC]= num2move;
    C->clusSize[who] -= num2move;
    C->nC++;
}

Boolean TrySplit(CLUSTERING *C){
    GRAPH *G = C->G;
    assert(0 < C->nC && C->nC < G->n-1);
    assert(C->clusters[C->nC-1]);
    assert(C->clusters[C->nC]==NULL);

    int oldNC = C->nC, who=0, oldSize=0, tries = 0;
    // prefer to split larger clusters with the drand48()
    while(oldSize < 2 || 1.0*C->clusSize[who]/_largestCluster < drand48()) {
	who = oldNC * drand48();
	oldSize = C->clusSize[who];
	assert(++tries < G->n);
    }
    double oldScore = ScoreOneCluster(C->G, C->clusters[who]);
    SplitCluster(C, who);
    assert(C->nC == oldNC+1);
    double newScore = ScoreOneCluster(C->G, C->clusters[who]) + ScoreOneCluster(C->G, C->clusters[oldNC]);
    if(newScore > oldScore) return true;
    else {
	SetUnion(C->clusters[who], C->clusters[who], C->clusters[oldNC]);
	SetFree(C->clusters[oldNC]);
	C->clusters[oldNC] = NULL;
	C->clusSize[oldNC]= 0;
	C->clusSize[who]  = oldSize;
	C->nC = oldNC;
	return false;
    }
}

Boolean TryMerge(CLUSTERING *C){
    GRAPH *G = C->G;
    assert(1 < C->nC && C->nC <= G->n);
    assert(C->clusters[C->nC-1]);

    int who[2]={-1,-1}, oldSize[2]={G->n,G->n};
    for(int i=0;i<2;i++) {
	int tries = 0;
	while((who[i] = C->nC*drand48()) == who[1-i] ||
	    // prefer to merge smaller clusters together instead of bigger ones
	    C->clusSize[who[i]] > G->n*drand48()) assert(++tries < G->n);
    }
    assert(who[0] != who[1]);
    assert(0 <= who[0] && who[0] < C->nC);
    assert(0 <= who[1] && who[1] < C->nC);
    double oldScore = ScoreOneCluster(C->G, C->clusters[who[0]]) + ScoreOneCluster(C->G, C->clusters[who[1]]);
    static SET *merged;
    if(merged) SetReset(merged); else merged = SetAlloc(C->G->n);
    SetUnion(merged, C->clusters[who[0]], C->clusters[who[1]]);
    double newScore = ScoreOneCluster(C->G, merged);
    if(newScore < oldScore) return false;
    else {
	SetFree(C->clusters[who[0]]);
	SetFree(C->clusters[who[1]]);
	int keep;
	if(who[0] < who[1]) keep=0; else keep=1;
	C->clusters[who[keep]] = merged;
	C->clusSize[who[keep]] = SetCardinality(merged);
	merged = NULL;
	C->nC--;
	C->clusters[who[1-keep]] = C->clusters[C->nC];
	C->clusSize[who[1-keep]] = C->clusSize[C->nC];
	C->clusters[C->nC] = NULL;
	C->clusSize[C->nC] = 0;
	return true;
    }
}

Boolean TryMove(CLUSTERING *C){
    assert(0 < C->nC && C->nC < C->G->n);
    assert(C->clusters[C->nC-1]);
    assert(C->clusters[C->nC]==NULL);

    // We want to bias towards a merge if there's "too many" clusters, and bias towards split if there's too few,
    // but but we need at least 2 before we can merge.
    double mergeProb = (C->nC-1.0)/C->G->n;
    
    if(drand48() < mergeProb) return TryMerge(C);
    else return TrySplit(C);
}

static CLUSTERING *_C;

void OutputClustering(int sig) {
    double score=0;
    for(int c=0; _C->clusters[c]; c++) {
	double CS = ScoreOneCluster(_C->G, _C->clusters[c]);
	printf("C %d n %d score %g:",c, _C->clusSize[c], CS);
	unsigned nodes[_C->G->n], n = SetToArray(nodes, _C->clusters[c]);
	for(int i=0; i<=n; i++) printf(" %d", nodes[i]);
	puts("");
    }
    exit(0);
}

int main(int argc, char *argv[])
{
    Boolean sparse=true, supportNodeNames=false, weighted=false;
    FILE *fp=Fopen(argv[1], "r");
    GRAPH *G = GraphReadEdgeList(fp, sparse, supportNodeNames, weighted);

    CLUSTERING *C = ClusteringAlloc(G);
    _C = C;

    // At this point we have one big cluster consisting of the entire graph... split it once and get going....
    SplitCluster(C, 0);

    signal(SIGINT, OutputClustering);
    double fullScore=0;
    while(fullScore>=0) {
	fullScore = ScoreClustering(C);
	if(TryMove(C)) printf("Status: %u clusters, largest %d, full score %g\n", C->nC, _largestCluster, fullScore);
    }
    return 0;
}

