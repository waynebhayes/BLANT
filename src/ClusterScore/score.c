#include "score.h"
#include "anneal.h"

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

CLUSTERING *_C;

unsigned ClusterEdgeCount(GRAPH *G, SET *c, unsigned num) {
    unsigned nodes[G->n], n = SetToArray(nodes, c), m=0;
    assert(n>0);
    assert(n==num);
    if(n==1) return 0;
    if(G->sparse) {
	for(int i=0; i<n; i++) {
	    int u = nodes[i];
	    for(int j=0; j < G->degree[u]; j++){
		int v = G->neighbor[u][j];
		if(v>u && SetIn(c, v)) ++m; // yes, I've checked that the v>u part is correct.
	    }
	}
    } else {
	for(int i=0; i<n; i++) for(int j=i+1;j<n;j++) if(GraphAreConnected(G,nodes[i],nodes[j])) ++m;
    }
    return m;
}

double ScoreOneCluster(GRAPH *G, SET *c, unsigned n) {
    unsigned m=ClusterEdgeCount(G,c,n);
    double ED = m/(n*(n-1)/2.0);
    return m*ED;
}

int _largestCluster;

double ScoreClustering(CLUSTERING *C) {
    double score=0;
    _largestCluster=0;
    for(int c=0; C->clusters[c]; c++) {
	assert(C->clusSize[c] > 0);
	assert(C->clusSize[c] < C->G->n);
	if(C->clusSize[c] > _largestCluster) _largestCluster = C->clusSize[c];
	double CS = ScoreOneCluster(C->G, C->clusters[c], C->clusSize[c]); // printf(" %g", CS);
	score += CS;
    }
    return score;
}


#define MIN_SIZE 9

// Remove some random set of elements from C[who] (consistent with MIN_SIZE); put them in C->clusters[nC], and increment C->nC
void SplitCluster(CLUSTERING *C, unsigned who) {
    GRAPH *G=C->G;
    assert(who < G->n);
    assert(C->clusSize[who] >= 2*MIN_SIZE);
    unsigned nodes[G->n], n = SetToArray(nodes, C->clusters[who]);
    assert(n>1);
    assert(n == C->clusSize[who]);
    assert(C->clusters[C->nC]==NULL);
    C->clusters[C->nC]=SetAlloc(G->n);
    int num2move = n/2 - (n/2 - MIN_SIZE)*drand48();
    assert(num2move >= MIN_SIZE && n-num2move >= MIN_SIZE);
    for(int i=0; i<num2move; i++) {
	int which=drand48()*n;
	until(SetIn(C->clusters[who], nodes[which])) which = drand48()*n;
	SetDelete(C->clusters[who], nodes[which]) ;
	SetAdd(C->clusters[C->nC], nodes[which]);
    }
    C->clusSize[C->nC]= num2move;
    C->clusSize[who] -= num2move;
    assert(C->clusSize[who] >= MIN_SIZE);
    assert(C->clusSize[C->nC] >= MIN_SIZE);
    C->nC++;
}

Boolean TrySplit(CLUSTERING *C){
    GRAPH *G = C->G;
    assert(0 < C->nC && C->nC < G->n-1);
    assert(C->clusters[C->nC-1]);
    assert(C->clusters[C->nC]==NULL);

    int oldNC = C->nC, who=0, oldSize=0, tries = 0;
    // prefer to split larger clusters with the drand48()
    while(oldSize < 2*MIN_SIZE || 1.0*C->clusSize[who]/_largestCluster < drand48()) {
	who = oldNC * drand48();
	oldSize = C->clusSize[who];
	if(++tries > 100*G->n) return false; // Fatal("too many tries in TrySplit");
    }
    double oldScore = ScoreOneCluster(C->G, C->clusters[who], C->clusSize[who]);
    SplitCluster(C, who);
    assert(C->nC == oldNC+1);
    double newScore = ScoreOneCluster(C->G, C->clusters[who], C->clusSize[who]) +
	ScoreOneCluster(C->G, C->clusters[oldNC], C->clusSize[oldNC]);
    double p = AnnealAcceptProb(oldScore, newScore), r = drand48();
    if(r < p) // if(newScore > oldScore)
	return true;
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
	// prefer to merge smaller clusters together instead of bigger ones
	while((who[i] = C->nC*drand48()) == who[1-i] || C->clusSize[who[i]] > G->n*drand48())
	    assert(++tries < G->n);
    }
    assert(who[0] != who[1]);
    assert(0 <= who[0] && who[0] < C->nC);
    assert(0 <= who[1] && who[1] < C->nC);
    double oldScore = ScoreOneCluster(C->G, C->clusters[who[0]], C->clusSize[who[0]]) +
	ScoreOneCluster(C->G, C->clusters[who[1]], C->clusSize[who[1]]);
    static SET *merged;
    if(merged) SetReset(merged); else merged = SetAlloc(C->G->n);
    SetUnion(merged, C->clusters[who[0]], C->clusters[who[1]]);
    double newScore = ScoreOneCluster(C->G, merged, SetCardinality(merged));
    double p = AnnealAcceptProb(oldScore, newScore), r = drand48();
    if(r < p) { // accept the merge
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
    else return false;
}

Boolean TryMove(CLUSTERING *C){
    assert(0 < C->nC && C->nC < C->G->n);
    assert(C->clusters[C->nC-1]);
    assert(C->clusters[C->nC]==NULL);

    // We want to bias towards a merge if there's "too many" clusters, and bias towards split if there's too few,
    // but but we need at least 2 before we can merge.
    double mergeProb = 1.0*C->nC/C->G->n;
    assert(0 < mergeProb && mergeProb < 1);
    mergeProb = sqrt(mergeProb); // bias it more towards 1 since splitting is occurring too much.
    
    if(drand48() < mergeProb) return TryMerge(C);
    else return TrySplit(C);
}

void OutputClustering(int sig) {
    double score=0;
    for(int c=0; _C->clusters[c]; c++) {
	double CS = ScoreOneCluster(_C->G, _C->clusters[c], _C->clusSize[c]);
	printf("C %d n %d m %d score %g:", c, _C->clusSize[c], ClusterEdgeCount(_C->G, _C->clusters[c], _C->clusSize[c]), CS);
	unsigned nodes[_C->G->n], n = SetToArray(nodes, _C->clusters[c]);
	for(int i=0; i<n; i++) printf(" %d", nodes[i]);
	puts("");
    }
    exit(0);
}

// returns whether the iteration was successful (ie., move accepted)
Boolean ScoreIteration(double pBad, Boolean running) {
    static unsigned fails, delay;
    Boolean success = TryMove(_C);
    if(running && ++delay > 10000) {
	double fullScore = ScoreClustering(_C);
	delay=0;
	printf("Status: pBad %g, %u clusters, largest %d, full score %g\n", PbadBufMean(), _C->nC, _largestCluster, fullScore);
    }
    return success;
}

