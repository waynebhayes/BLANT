#include "misc.h"
#include "sets.h"
#include "graph.h"
#include "queue.h"
#include "stats.h"
#include "rand48.h"
#include <stdio.h>

/* This should be an exact substitute for a similarly named bash+awk script in the blant-clusters.sh suite.
    The command-line arguments are:
	./build-clusters ED minClus edgeList.el sorted-blant-output
   where ED is the desired minimum edge density of the clusters to discover, minClus is the desired minimum cluster size;
	possible minClus values:
	    1) an integer>=3, that's the explicit smallest cluster to output (a good choice is k+1)
	    2) minClus in (0,1) defines a rough desired p-value, used to estimate the smallest "significant" cluster size
	    3) minClus = 0 means we attempt to figure out a decent p-value based on n+ED, and then apply #2.
   Finally, we expect two files: the edgeList (possibly weighted), and the sorted output from blant (usually stdin,
   specified by the filename "-"). The lines of the second file must be of the form:
     score node canonID nn v1 v2 v3... v_nn
   where "node" has "score" for canon/orbit ID, and has nn neigbors that appeared in said ID canons/orbits
*/

#define edgePredict 0 // try to predict this many missing edges from cliques/communities
#define minEdgeMean 0.5 // heuristic: avoid clusters connected by too-weak edges
#define DENSITY_LEEWAY 0.95 // we temporarily allow clusters + graphlets to be this factor below the requested density
#define OVERLAP 0.5


// These basically all have to be global in order to maintain a 1-to-1 correspondence with the AWK code
static double edgeDensity, *count;
static unsigned *line2node, *node2line; 
static SET *S, **graphletNeighbors, *visitedQ;
GRAPH *G;
QUEUE *Q;
// because I don't like the old names
#define QueueAdd QueuePut
#define QueueLength QueueSize
#define QueueNext QueueGet

unsigned EdgesIntoS(unsigned v) { // count edges from v into S, not including v (which can be is S)
    unsigned i, edgeHits=0, n1=SetCardinality(S),n2, Snodes[n1];
    n2=SetToArray(Snodes, S); assert(n1==n2);
    for(i=0; i<n1; i++) if(Snodes[i]!=v && GraphAreConnected(G,Snodes[i],v)) ++edgeHits;
    return edgeHits;
}
#if 0 // awk code hopefully not needed
function WgtEdgesIntoS(v,       edgeWgts,u) { // WEIGHT of edges from v into S, not including v (which can be is S)
    edgeWgts=0;
    for(u in S) if(u!=v && (v in edge[u])){
	ASSERT(edge[v][u],"edge error");
	edgeWgts+=edge[u][v];
    }
    return edgeWgts;
}
#endif

// Count the number of edges in the subgraph induced on nodes in T, and populate degreeInS with their degree
unsigned InducedEdges(SET *T, unsigned *array, unsigned *degreeInS) {
    unsigned i, j, u,v, m=0, n=SetToArray(array, T);
    for(i=0;i<n;i++) degreeInS[array[i]]=0; // set ONLY the relevant values to zero
    for(i=0;i<n;i++) {
	u=array[i];
	for(j=i+1;j<n;j++) {
	    v=array[j];
	    if(GraphAreConnected(G,u,v)) { ++degreeInS[u]; ++degreeInS[v]; ++m; }
	}
    }
    return m;
}

Boolean highRelClusCount(unsigned u, unsigned v) { // Heuristic... the order of (u,v) MATTERS
    if(count[u]==0) return 1;
    if(count[v]==0) return 0;
    else return count[v]*1.0/count[u]>=sqrt(edgeDensity)*0.7; // values in [0.5,0.9] work well
}

// Append the neighbors of u to the cluster being created around origin (which is the "seed" node of the cluster)
void AppendNeighbors(unsigned u, unsigned origin) {
    int n=SetCardinality(graphletNeighbors[u]);
    if(n) {
	// The following three lines will append the neighbors in order of most edges back into S
	// It seems reasonable, but is more expensive seems to have no significant difference. :-(
	//for (v in graphletNeighbors[u]) edgesIntoS[v] = EdgesIntoS(v);
	//PROCINFO["sorted_in"] = "@val_num_desc"; // for loop through edgesIntoS, largest first
	//for (v in edgesIntoS)
	// The default is to append them in random order... no significant difference on small networks
	// PROCINFO["sorted_in"]="randsort";
	unsigned array[n]; SetToArray(array, graphletNeighbors[u]);
	int i; unsigned which;
	for(i=n;i>0;i--) {
	    which = i*drand48(); // pick a random element from a shrinking array
	    unsigned v=array[which]; array[which] = array[i-1];
	    if(!SetIn(visitedQ,v) && (!node2line[v] || node2line[v] > node2line[origin]) && highRelClusCount(u, v)) {
		QueueAdd(Q, (foint)v);
		SetAdd(visitedQ,v);
	    }
	}
	assert(i==0 && which==0);
    }
}

void BuildClusters(double desiredEdgeDensity, double minClusArg, FILE *fp){
    unsigned seed = GetFancySeed(false), maxEdgesG = G->n*(G->n-1)/2, minClus;
    double eps=G->m*1.0/maxEdgesG;
    assert(minClusArg >= 0 && (minClusArg < 1 || (minClusArg>=3 && minClusArg == (int)minClusArg)));
    edgeDensity = desiredEdgeDensity;
    if(minClusArg>=1) {
	if(minClusArg<3 || minClusArg != (int)minClusArg) Fatal("when minClus>=1, it must be an integer >=3");
	minClus = (int)minClusArg;
    } else {
	double pVal;
	if(minClusArg>0) pVal = minClusArg;
	else pVal=1.0/(maxEdgesG*maxEdgesG); // heuristic to make pVal small enough to get significant clusters
	int e;
	for(e=1;e<G->m;e++) if(IntPow(eps,e)<pVal) break; // #of concurrently chosen edges to get < pVal
	for(minClus=2;minClus<G->n;minClus++)
	    if(edgeDensity*CombinChoose(minClus,2)>=e) break; // #nodes expected from the above to get those edges
	if(minClus < 3) minClus=3;
	fprintf(stderr, "minClus %d (%d edges = pVal %g < %g for ED %g)\n",minClus,e,IntPow(eps,e),pVal,eps);
    }

    char *env = getenv("ONLY_FIRST");
    Boolean onlyFirst = (env && strcmp(env,"1")==0);
    env = getenv("ONLY_BIGGEST");
    Boolean onlyBiggest = (env && strcmp(env,"1")==0);

    // Allocate the space for the globals
    count = Calloc(G->n, sizeof(double)); // the total count of graphlets above ED for each node
    unsigned line2nodeSize = G->n;
    line2node = Calloc(line2nodeSize, sizeof(unsigned)); // the node number indexed by the LINE NUMBER
    node2line = Calloc(G->n, sizeof(unsigned)); // opposite of the above: what node was on line FNR?
    graphletNeighbors = Calloc(G->n, sizeof(SET*));
    for(int i=0;i<G->n;i++) graphletNeighbors[i] = SetAlloc(G->n);
    S = SetAlloc(G->n);
    visitedQ = SetAlloc(G->n);
    Q = QueueAlloc(G->n);

    unsigned FNR=0; // line number like AWK except starting at 0 rather than 1
    while(!feof(fp)) {
	double score;
	char name[BUFSIZ]; strcpy(name, "<undefined>"); // because we might print it in the case of error
	unsigned canonID, numNeigh; // Note: we don't really USE the canonID for anything...
	// score node canonID nn v1 v2 v3... v_nn
	int numRead = fscanf(fp, "%lf%s%u%u", &score, name, &canonID, &numNeigh);
	if(numRead != 4) Fatal("couldn't find score,node,ID,nn; got only %d values: %g %s %u %u",
	    numRead, score, name, canonID, numNeigh);
	unsigned node = GraphNodeNameToInt(G, name);
	if(onlyFirst && count[node]) /* technically skip this node but we still have to read+discard all the neighbors */ ;
	else {
	    if(FNR>=line2nodeSize) {
		line2nodeSize *= 2;
		line2node = Realloc(line2node, line2nodeSize*sizeof(unsigned));
	    }
	    line2node[FNR]=node;
	    node2line[node]=FNR;
	    count[node] += score;
	}
	for(int i=0;i<numNeigh;i++) {
	    char name2[BUFSIZ]; strcpy(name2, "<undefined>"); // because we might print it in the case of error
	    if(fscanf(fp, "%s ", name2) != 1) Fatal("couldn't get neighbor %d of node %s",i,name);
	    unsigned neigh = GraphNodeNameToInt(G, name2);
	    if(!(onlyFirst && count[node])) {
		SetAdd(graphletNeighbors[node],neigh);
		SetAdd(graphletNeighbors[neigh],node); // do we really need both?
	    }
	}
	++FNR;
    }
    // END block
    double density_leeway = DENSITY_LEEWAY;
    if(edgeDensity==1) density_leeway=1; // no leeway from cliques
    // make started an empty set; started is a list of nodes we should NOT start a new BFS on
    SET *started = SetAlloc(G->n);
    unsigned numClus=0;
    unsigned Sorder[G->n]; // array that maintains ORDER in which nodes were added
    for(unsigned start=0; start<(int)(FNR*edgeDensity); start++) { // look for a cluster starting on line "start".
	unsigned origin=line2node[start];
	if(SetIn(started,origin)) continue;
	//printf "line %d origin %s;", start, origin > "/dev/stderr"
	SetAdd(started,origin);
	SetEmpty(S);
	SetEmpty(visitedQ);
	unsigned misses=0, Slen; // how many nodes have been skipped because they did not work?
	if(QueueLength(Q)>0) QueueEmpty(Q); // Empty the queue
	QueueAdd(Q, (foint)origin);
	SetAdd(visitedQ,origin);
	unsigned edgeCount = 0; // wgtEdgeCount = 0;
	while(QueueLength(Q)>0){
	    unsigned u = QueueNext(Q).ui;
	    assert(!SetIn(S,u));
	    unsigned newEdgeHits = EdgesIntoS(u); // wgtEdgeHits = WgtEdgesIntoS(u);
	    edgeCount += newEdgeHits; // wgtEdgeCount += wgtEdgeHits;
	    SetAdd(S,u);
	    Slen = SetCardinality(S);
	    //fprintf(stderr, "start %d u %u newEdgeHits %u edgeCount %u Slen %u\n", start,u,newEdgeHits,edgeCount,Slen);
	    assert(Slen > 0); // "Slen must be > 0 but is "Slen);
	    // if(Slen >= '$maxClus') break;
	    // CAREFUL: calling InducedEdges(S) every QueueNext() is VERY expensive; uncommenting the line below
	    // slows the program by more than 100x (NOT an exaggeration)
	    // WARN(edgeCount == InducedEdges(edge,S),"Slen "Slen" edgeCount "edgeCount" Induced(S) "InducedEdges(edge,S));
	    Sorder[Slen-1]=u; // no need to delete this element if u fails because Slen will go down by 1
	    if (Slen==1)
		AppendNeighbors(u, origin);
	    else {
		unsigned maxEdges = Slen*(Slen-1)/2;
		if(edgeCount < maxEdges*edgeDensity * density_leeway) {
		    SetDelete(S,u); // u drops the edge density too low, so nuke it
#if edgePredict
		    if(edgePredict) {
			if((edgeCount+edgePredict) >= maxEdges*edgeDensity) {
			    // just ONE edge would put us over the threshold; predict it should exist
			    unsigned numPredict=0;
			    for(uu in S)if(!(u in edge[uu]))
				//printf "predictEdge-%d\t%s\t%s\n",++numPredict, u,uu > "/dev/stderr"
			    WARN(numPredict<=edgePredict, "should only get "edgePredict" new edge predictions but got "numPredict);
			}
		    }
#endif
		    edgeCount -= newEdgeHits; // wgtEdgeCount -= wgtEdgeHits;
		    if(++misses > MAX(Slen, G->n/100.0)) break; // 1% of number of nodes is a heuristic...
		    // keep going until count decreases significantly; remove duplicate clusters in next awk
		} else {
		    misses=0;
		    AppendNeighbors(u, origin);
		}
	    }
	}
	Slen = SetCardinality(S);
	//printf " |S|=%d edgeMean %g", length(S), wgtEdgeCount/(edgeCount?edgeCount:1) > "/dev/stderr"
	STAT *stat = StatAlloc(0,0,0,0,0);
	if(QueueLength(Q)==0 && Slen > 1) {
	    // post-process to remove nodes with too low in-cluster degree---quite relevant for lower density
	    // communities where a node may be added because it does not
	    // reduce the *mean* degree of the cluster, but it does not really have sufficiently strong
	    // connections to the actual members of the community. We use two criteria:
	    // 1) the in-cluster degree is more than 3 sigma below the mean, or
	    // 2) the in-cluster degree is less than 1/3 the mode of the in-cluster degree distribution.
	    // The latter was added in response to our performance on the LFR graphs, but it does not appear
	    // to hurt performance anywhere else.
	    StatReset(stat); 
	    unsigned array[G->n], degreeInS[G->n], degFreq[G->n]; for(int i=0;i<G->n;i++) degFreq[i]=0;
	    unsigned tmpEdge = InducedEdges(S, array, degreeInS); // this call populates degreeInS
	    assert(tmpEdge == edgeCount);
	    for(int i=0; i<Slen; i++){ unsigned v=array[i];
		StatAddSample(stat, degreeInS[v]);
		++degFreq[degreeInS[v]];
	    }
	    unsigned maxFreq=0, degMode=0, dMin=MAX(1,StatMin(stat)), dMax=StatMax(stat);
	    for(unsigned d=dMax;d>=dMin;d--)
		if(degFreq[d] >= maxFreq) maxFreq=degFreq[degMode=d]; // use >= to extract SMALLEST mode
	    //printf " deg mean %g stdDev %g maxFreq %d degMode %d; pruning...", StatMean(""), StatStdDev(""), maxFreq, degMode > "/dev/stderr"
	    //PROCINFO["sorted_in"] = "@val_num_asc"; // for loop through in-degrees, smallest first
	    double sMean=StatMean(stat), sDev=StatStdDev(stat);
	    for(int i=0; i<Slen; i++){ unsigned v=array[i];
		if(degreeInS[v] < sMean - 3*sDev || degreeInS[v] < degMode/3.0) {
		    //printf " %s(%d)", v, degreeInS[v] > "/dev/stderr";
		    SetDelete(S,v);
		    edgeCount -= EdgesIntoS(v);
		    //wgtEdgeCount -= WgtEdgesIntoS(v);
		}
	    }
	    //printf " |S|=%d", _statN[""] > "/dev/stderr"
	}
	Slen = SetCardinality(S);
	unsigned array[G->n], degreeInS[G->n];
	unsigned tmpEdge = InducedEdges(S, array, degreeInS); // this call populates degreeInS
	assert(tmpEdge == edgeCount);
	//PROCINFO["sorted_in"] = "@val_num_asc"; // for loop through in-degrees, smallest first
	for(int i=0; i<Slen; i++) { unsigned v=array[i];
	    int tmpLen = SetCardinality(S);
	    if(tmpLen<2) break;
	    if(edgeCount*1.0/(tmpLen*(tmpLen-1)/2) >= edgeDensity) break; // break once ED is above threshsold
	    SetDelete(S,v);
	    edgeCount -= EdgesIntoS(v); // wgtEdgeCount -= WgtEdgesIntoS(v);
	    tmpEdge = InducedEdges(S, array, degreeInS); // this call populates degreeInS
	    assert(tmpEdge == edgeCount);
	}
	Slen=SetCardinality(S);
	if(Slen>1) {
	    double ed = edgeCount*1.0/(Slen*(Slen-1)/2);
	    assert(ed >= edgeDensity);
	}
	//printf " |S|=%d", _statN[""] > "/dev/stderr"
	if(Slen>=minClus /* && wgtEdgeCount*1.0/edgeCount>='$minEdgeMean'*/ ) {
	    unsigned maxEdges=Slen*(Slen-1)/2;
	    //print " ACCEPTED" > "/dev/stderr";
	    ++numClus; printf("%d %d %g", Slen, edgeCount, 1.0*edgeCount);
	    StatReset(stat);
	    tmpEdge = InducedEdges(S, array, degreeInS); // this call populates degreeInS
	    assert(tmpEdge == edgeCount);
	    for(int i=0;i<Slen;i++){unsigned u=array[i]; printf(" %s", G->name[u]); StatAddSample(stat, degreeInS[u]);}
	    //printf "final |S|=%d mean %g stdDev %g min %d max %d:\n", _statN[""], StatMean(""), StatStdDev(""), StatMin(""), StatMax("") > "/dev/stderr"
	    puts("");
	    if(onlyBiggest) break;
	    // now mark as "started" the first half of S that was filled... or more generally, some fraction.
	    // We choose the top (1-OVERLAP) fraction because OVERLAP is meant to be less stringent as it
	    // approaches 1, which in the case of THIS loop means that if OVERLAP is close to 1, we want
	    // to eliminate a SMALLER proportion (ie., allow more) of future BFS starts, so if OVERLAP=1,
	    // we eliminate NOTHING in this loop.
	    for(int i=1;i<Slen*(1-OVERLAP);i++) SetAdd(started,Sorder[i-1]);
	}
	//else print " REJECTED" > "/dev/stderr";
    }
    // end of AWK code
}

int main(int argc, char *argv[]) {
    if(argc!=5) Fatal("expecting 4 args: edgeDensity minClus edgeList.el blant-output");
    double ED = atof(argv[1]);
    double minClusArg = atof(argv[2]);
    FILE *graphFile = Fopen(argv[3],"r"); assert(graphFile);
    FILE *blantFile;
    if(strcmp(argv[4],"-")==0) blantFile=stdin;
    else blantFile = Fopen(argv[4],"r");
    assert(blantFile);
    Boolean self=false, directed=false, weighted=false;
    G = GraphReadEdgeList(graphFile, self, directed, weighted);
    fclose(graphFile);

    long seed;
    char *env = getenv("RANDOM_SEED");
    if(env) seed = atol(env);
    else seed = GetFancySeed(false);
    srand48(seed);

    BuildClusters(ED, minClusArg, blantFile);
}
