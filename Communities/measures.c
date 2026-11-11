#include "structs.h"
#include "measures.h"
#define TARGET_EDGE_DENSITY 0.5
// fakeN is the number of nodes inside the community I TELL the measure

// Because of the nature of how I implemented faster rejects, most of the time 
// its not the actual true number of nodes in the community

double IntraEdgeDensity(PARTITION * P, COMMUNITY *C, int fakeN){
    // Always does trivial solution

    if(fakeN < 2)
	return 0;
    return C->edgesIn / (fakeN *(fakeN - 1) / 2.0);
}

double InterEdgeDensity(PARTITION * P, COMMUNITY *C, int fakeN){
    // Always does trivial solution

    int tot = C->G->n - fakeN;
    //printf("tot = %d, edgesOut = %d\n", tot, edgesOut);
    
    return (fakeN == 0 || tot == 0) ? 0 : (double)C->edgesOut/(fakeN * tot);
}

double NewmanAndGirvan(PARTITION * P, COMMUNITY * C, int fakeN){  
    double firstTerm = C->edgesIn / (double)P->G->numEdges;
    double secondTerm = (2.0 * C->edgesIn + C->edgesOut) / (2.0 * P->G->numEdges);
    double ans = firstTerm - (secondTerm * secondTerm); 
    //printf("edgeCount = %d, P->G->numEdges = %d, firstTerm = %f, secondTerm = %f, ans = %f\n", C->edgesIn, P->G->numEdges, firstTerm, secondTerm, ans); 
    return ans;
}

double Expansion(PARTITION * P, COMMUNITY * C, int fakeN){
    // Always does trivial solution of making one big community
    // Horrible objective function
    int eOut = CommunityEdgeOutwards(P, C); 
    int eIn = CommunityEdgeCount(P, C);
    return eOut + eIn == 0 ? 9999999 : (double)(eOut / (eOut + eIn));  
}

double NormalizedCut(PARTITION * P, COMMUNITY * C, int fakeN){
    // Also always does trivial solution
    int eIn = C->edgesIn;
    int eOut = C->edgesOut;
    int denom1 = 2 * eIn + eOut; 
    int denom2 = 2 * (P->G->n - eIn) + eOut;
    if(denom1 == 0 || denom2 == 0){
	return 99999;
    }
    return (double)(eOut / (2 * eIn + eOut)) + (double)(eOut / (2 * (P->G->n - eIn) + eOut));
}

double CPM(PARTITION * P, COMMUNITY * C, int fakeN){
    // Constant Potts Model (Traag, Van Dooren & Nesterov, 2011)
    return C->edgesIn - 0.5 * (fakeN * (fakeN - 1)) / 2.0;
}

double Conductance(PARTITION * P, COMMUNITY * C, int fakeN){
    int denom = 2 * C->edgesIn + C->edgesOut; // Just the sum of degree of all nodes 
    if(denom == 0){
	// return infinity if denom is 0
	// Should never happen???
	return 99999;
    }
    double ans = (double)C->edgesOut/(double)denom;
    //printf("Conductance denom =  %d, edgesOut = %d, ans = %g\n", denom, C->edgesOut, ans);
    return ans;
}

double HayesScore(PARTITION * P, COMMUNITY *C, int fakeN){ 
#if DEBUG
    printf("Hayes Com %d, inEdges = %d, C->n = %d\n", C->id, C->edgesIn, fakeN);
#endif
    if(fakeN < 2){
	return 0;
    }
    // (double) first: edgesIn*edgesIn overflows int once a community has more than ~46k internal edges
    return (double)C->edgesIn * C->edgesIn / ((fakeN * (fakeN-1))/2.0);
    
}

double Modularity(PARTITION * P, COMMUNITY * C, int fakeN){
    // TODO: Actually implement modularity
    return 0; //IntraEdgeDensity(P, C, fakeN);
}


