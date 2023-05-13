#include "misc.h"
#include "sets.h"
#include "graph.h"
#include "hashmap.h"
#include "priority-queue.h"
#include "rand48.h"

unsigned _Gn; // number of nodes in the INPUT (user's) network
#define MAX_CLUSTERS 1000000 // maximum number of clusters in the cluster similarity graph

typedef struct _commType {
    unsigned n,m,nc2,k, index;
    double ED, score;
    SET *nodes;
	//hashmap_t *kin;
} CLUSTER;

static int ClusterScoreCompare(foint f1, foint f2)
{
    CLUSTER *c1 = (CLUSTER*)f1.v, *c2 = (CLUSTER*)f2.v;
    // abstract "smallest-at-top" wins, so reverse it for largest
    if(c1->score < c2->score) return +1;
    if(c1->score > c2->score) return -1;
    return 0;
}


// read the output of blant-clusters.sh, looks like this (after it's put through sed to remove non-digits)
// n  m   nc2 k ED%     list
// 19 154 171 5 90.0585 241 44 129 184 227 168 178 112 256 197 293 189 633 106 14 636 171 15 190
// returns number of nodes
CLUSTER *ReadCluster(FILE *fp)
{
    unsigned i,v;
    CLUSTER *c = (CLUSTER*) Calloc(1,sizeof(CLUSTER));
    c->nodes = SetAlloc(_Gn);
    assert(fscanf(fp, "%u%u%u%u%lf", &c->n, &c->m, &c->nc2, &c->k, &c->ED)==5);;
    c->ED /= 100; // blant outputs a percentage, we want a fraction
    for(i=0; i<c->n;i++) {assert(fscanf(fp, "%u", &v)==1); SetAdd(c->nodes,v);}
	//c->kin=hashmap_new();
    return c;
}

CLUSTER *_cluster[MAX_CLUSTERS]; // gets dynamically bigger as necessary, starting at _Gn
unsigned _numClus, *_finalMemberships; // array [_Gn] counting number of clusters each node in INPUT graph belongs to
double _overlapThresh = 0.5;
GRAPH *_clusterSimGraph;
SET *_finalComm, *_finalClusVisited;
SET **_clusterMemberships; // _clusterMemberships[u] = set of clusters node u \in G (input network) is in


enum Measure { undef, MEASURE_OMOD, MEASURE_EDN };
static enum Measure measure = undef;
double _stopT = 0; // stop threshold

PRIORITY_QUEUE *_PQ;
GRAPH *_inputNet;

// unsigned kin(CLUSTER *c, unsigned u){
// 	unsigned kinu;
// 	if (hashmap_get(c->kin, u, &kinu) == MAP_MISSING){
// 		kinu=0;
// 		SET *c_nodes=c->nodes;
// 		unsigned v, m;
// 		FOREACH(v,c_nodes){
// 			if(GraphAreConnected(_inputNet,u,v)){
// 				kinu++;
// 			}
// 		}
// 		hashmap_put(c->kin,u, kinu);
// 		printf("Putting kin for %d: %d", u, kinu);
// 	}
// 	else{
// 		printf("Getting kin for %d: %d", u, kinu);
// 	}
// 	return kinu;
// }
double scoreOfNodeInCommunity(CLUSTER *c, unsigned u, unsigned membershipCount)
{
    assert(measure != undef);
    switch(measure) {
    case MEASURE_EDN: return ( 1.0 / membershipCount ) * c->ED;
	break;
    case MEASURE_OMOD: Apology("Sorry"); 
		// unsigned kinu = kin(c,u);
		// unsigned ku= _inputNet->degree[u];
		// return ( (kinu - ( ku - kinu ) ) / ku ) * ( 1.0 / membershipCount ) * c->ED * ( 1.0 / c->n );
	break;
    default: Fatal("unknonw measure in scoreOfNodeInCommunity");
	return (-1); break;
    }
}


double communityScore(CLUSTER *c)
{
    double CS = 0;
    assert(c->n==SetCardinality(c->nodes));
    unsigned m; SET *s=c->nodes;
    FOREACH(m,s)
	CS += scoreOfNodeInCommunity(c, m, _finalMemberships[m]+1);
    return CS;
}

double _currentScore;

double potentialScore(CLUSTER *c)
{
    double P = _currentScore;
	P*=SetCardinality(_finalComm);
    SET *c_nodes=c->nodes;
    unsigned u, m;
    FOREACH_DECLARE(m,_finalComm);
    FOREACH(u,c_nodes) {
	if(_finalMemberships[u]!=0) {
	    FOREACH_LOOP(m,_finalComm) {
		CLUSTER *c2 = _cluster[m];
		if(SetIn(c2->nodes, u))
		    P=P-scoreOfNodeInCommunity(c2, u, _finalMemberships[u])+scoreOfNodeInCommunity(c2, u, _finalMemberships[u] + 1);
	    }
	}
    }
    P+=communityScore(c);
	P/=(SetCardinality(_finalComm)+1);
    return P;
}



void addToResult(CLUSTER *c)
{
    double P=potentialScore(c), diff=P-_currentScore;
    if (diff < _stopT) return;
    _currentScore=P;
    unsigned u; SET *cn=c->nodes;
    FOREACH(u,cn) _finalMemberships[u]++;
    SetAdd(_finalComm,c->index);
}

void expand(CLUSTER *c)
{
    if(GraphDegree(_clusterSimGraph, c->index)) {
	int j; 
	for(j=0;j<_clusterSimGraph->degree[c->index];j++) {
	    unsigned n=_clusterSimGraph->neighbor[c->index][j];
	    if(!SetIn(_finalClusVisited, n)) {
		double P=potentialScore(_cluster[n]), diff=P-_currentScore;
		if ( diff > _stopT ) {
		    _cluster[n]->score = P;
		    PriorityQueueInsert(_PQ,(foint)(void*)_cluster[n]);
		}
		SetAdd(_finalClusVisited,n);
	    }
	}
    }
}



// Compute the overlap amount (number of nodes of overlap) between this cluster and current list; result in globals below
SET *_overlapMatches; // set of other clusters that match c after calling ComputeClusterOverlap
static unsigned _overlap[MAX_CLUSTERS]; // the actual counts; only the members above are guaranteed valid

void ComputeClusterOverlap(const CLUSTER *c)
{
    static SET *dirty; // for every cluster pc in _overlap[], is the _overlap[pc] value dirty, or valid?
			// Use a bitvec since lots of deletions will occur and we don't want to keep qsorting it.
    if(dirty == NULL) {
	dirty = SetAlloc(MAX_CLUSTERS);
	assert(!_overlapMatches);
	_overlapMatches = SetAlloc(MAX_CLUSTERS);
	assert(_Gn > 0);
    } else
	SetReset(_overlapMatches);

    unsigned u; SET *cNodes=c->nodes;
    FOREACH(u,cNodes) {
	if(_clusterMemberships[u]) {
	    unsigned prevCluster; SET *myClusters = _clusterMemberships[u];
	    FOREACH(prevCluster, myClusters) {
		if(SetIn(dirty, prevCluster)) { _overlap[prevCluster] = 0; SetDelete(dirty,prevCluster);}
		++_overlap[prevCluster];
		SetAdd(_overlapMatches, prevCluster);
	    }
	}
    }
    // unsigned pc; FOREACH(pc, _overlapMatches) SetAdd(dirty, pc);
    SetCopy(dirty, _overlapMatches);
}

void init(int argc, char *argv[]){
	if(argc!=4) Fatal("USAGE: stopThresh measure inputNetwork.el");
    _stopT = atof(argv[1]); assert(_stopT>=0);
	printf("1");
	if(strcmp(argv[2], "OMOD") == 0) measure=MEASURE_OMOD;
	else if(strcmp(argv[2], "EDN") == 0) measure=MEASURE_EDN;
	assert(measure!=undef);
    printf("2");
	_finalMemberships = Calloc(_Gn, sizeof(_finalMemberships[0]));
    _clusterMemberships = Calloc(_Gn, sizeof(SET*));
	Boolean sparse = true, names=false;

	FILE* netFile = fopen(argv[3],"r");
	_inputNet = GraphReadEdgeList(netFile,sparse,false);
	_Gn=_inputNet->n;
	assert(_Gn>0);
	printf("Initialization done. Graph has %d nodes", _Gn);
	
	_clusterSimGraph = GraphAlloc(MAX_CLUSTERS, sparse, names);
    GraphMakeWeighted(_clusterSimGraph);	
}



int main(int argc, char *argv[])
{
	printf("Hola?");
    unsigned i, line=0;
    init(argc,argv);
    while(!feof(stdin) && _numClus < MAX_CLUSTERS) {
	++line; // numbering from 1
	CLUSTER *c = ReadCluster(stdin);

	fprintf(stderr, " %d(%d)", line, c->n);
	// Either remember, or forget, this cluster based on overlap with previous ones
	ComputeClusterOverlap(c);
	FOREACH(i,_overlapMatches) {
	    int maxCard = MAX(SetCardinality(_cluster[i]->nodes) , SetCardinality(c->nodes));
	    if(_overlap[i] / (1.0*maxCard) > _overlapThresh) {
		//printf("Skipping cluster %d\n", _numClus);
		SetFree(c->nodes);
		Free(c);
		c=NULL;
		break;
	    }
	}
	if(c) { // it was NOT disqualified based on too much overlap with previous cluster
	    fprintf(stderr, "A");
	    c->index = _numClus;
	    for(i=0;i<_numClus;i++) if(SetIn(_overlapMatches,i)) {
		//printf("sim %d %d = %d\n",i,_numClus,_overlap[i]);
		GraphSetWeight(_clusterSimGraph,i,_numClus,_overlap[i]);
	    }
	    unsigned u; SET *cNodes = c->nodes;
	    FOREACH(u,cNodes) { // Record new membership across the nodes of this cluster
		if(!_clusterMemberships[u]) _clusterMemberships[u] = SetAlloc(MAX_CLUSTERS);
		SetAdd(_clusterMemberships[u], _numClus);
	    }
	    _cluster[_numClus++] = c;
	    fscanf(stdin, " ");
	} else
	    fprintf(stderr, "R");
    }
    fprintf(stderr, "\n");
    if(!feof(stdin) && _numClus == MAX_CLUSTERS)
	Warning("cluster reading stopped at MAX_CLUSTERS %d; processing cluster graph anyway\n", MAX_CLUSTERS);

#if 0 // debugging just to check we read things correctly
    printf("Read %d clusters:\n", _numClus);
    for(i=0; i<_numClus; ++i) {
	CLUSTER *c=_cluster[i];
	printf("Cluster %d: n=%d m=%d nc2=%d k=%d ED=%g%%: ", i, c->n, c->m, c->nc2, c->k, c->ED);
	SetPrint(_cluster[i]->nodes);
	puts("");
    }
#endif

    _PQ = PriorityQueueAlloc(_numClus, ClusterScoreCompare, NULL);

    // For each connected component, mark it off and pick one node at random to put in the PQ
    int numVisited = 0, Varray[_numClus], rootNode=0;
    SET *visitedDFS = SetAlloc(_numClus);
    _finalClusVisited = SetAlloc(_numClus);
    for(rootNode=0;rootNode<_numClus;rootNode++) {
	if(!SetIn(visitedDFS, rootNode)) {
	    int prevNum = numVisited;
	    GraphVisitCC(_clusterSimGraph, rootNode, visitedDFS, Varray, &numVisited);
	    assert(numVisited == SetCardinality(visitedDFS));
	    // pick a node at random from the most recent CC
	    int CCsize = numVisited-prevNum, randOffset=CCsize*drand48(), ccRandNode=Varray[prevNum+randOffset];
	    _cluster[ccRandNode]->score = communityScore(_cluster[ccRandNode]);
	    PriorityQueueInsert(_PQ, (foint)(void*)_cluster[ccRandNode]);
	    SetAdd(_finalClusVisited, ccRandNode);
	}
    }
    printf("Cluster graph has %d nodes; PQ size is %d\n", _numClus, PriorityQueueSize(_PQ));
    _finalComm = SetAlloc(_numClus);

    while (PriorityQueueSize(_PQ)) {
	CLUSTER *c = PriorityQueueNext(_PQ).v;
	addToResult(c);
	expand(c);
    }

    if(measure==MEASURE_OMOD) printf("Qov=%g\n",_currentScore);
    else if(measure==MEASURE_EDN) printf("EDN=%g\n",_currentScore);


    unsigned m;
    FOREACH(m,_finalComm) {
	CLUSTER *c=_cluster[m];
	printf("%d nodes, %d of %d edges from k %d (%g%%):", c->n, c->m, c->nc2, c->k, 100*c->ED);
	unsigned u; SET *cn=c->nodes;
	FOREACH(u,cn) printf(" %d", u);
	puts("");
    }
    return 0;
}
