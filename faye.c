// modelled after faye by Tuong Do
static SET *SampleGraphletFaye(SET *V, int *Varray, GRAPH *G, int k, int whichCC)
{
    /* Faye: Add a visited array to keep track of nodes. Initialize to 0 */
    int visited[G->n];
    static SET *outSet;
    static int numIsolatedNodes;
    if(!outSet)
       outSet = SetAlloc(G->n);  // we won't bother to free this since it's static.
    else if(G->n > outSet->n)
	SetResize(outSet, G->n);
    else
	SetEmpty(outSet);
    int v1, v2, i;
    int nOut = 0, outbound[G->n]; // vertices one step outside the boundary of V
    assert(V && V->n >= G->n);
    SetEmpty(V);
    int edge;
    do {
	edge = G->numEdges * RandomUniform();
	v1 = G->edgeList[2*edge];
    } while(!SetIn(_componentSet[whichCC], v1));
    v2 = G->edgeList[2*edge+1];
    SetAdd(V, v1); Varray[0] = v1;
    SetAdd(V, v2); Varray[1] = v2;

    /* Faye: Mark v1 and v2 as visited */
    visited[v1] = 1;
    visited[v2] = 1;    

    // The below loops over neighbors can take a long time for large graphs with high mean degree. May be faster
    // with bit operations if we stored the adjacency matrix... which may be too big to store for big graphs. :-(
    for(i=0; i < G->degree[v1]; i++)
    {
	int nv1 =  G->neighbor[v1][i];
	if(nv1 != v2)
	{
	    assert(!SetIn(V, nv1)); // assertion to ensure we're in line with faye
	    if (!visited[nv1]) { /* Faye: Check if it's visited */
            SetAdd(outSet, (outbound[nOut++] = nv1));
            visited[nv1] = 1;
        }
	}
    }
    for(i=0; i < G->degree[v2]; i++)
    {
	int nv2 =  G->neighbor[v2][i];
	if(nv2 != v1 && !SetIn(outSet, nv2))
	{
	    assert(!SetIn(V, nv2)); // assertion to ensure we're in line with faye
	    if (!visited[nv2]) { /* Faye: Check if it's visited */
            SetAdd(outSet, (outbound[nOut++] = nv2));
            visited[nv2] = 1;
        }
	}
    }
    for(i=2; i<k; i++)
    {
	int j;
	if(nOut == 0) // the graphlet has saturated it's connected component
	{
	    assert(SetCardinality(outSet) == 0);
	    assert(SetCardinality(V) < k);
#if ALLOW_DISCONNECTED_GRAPHLETS
	    /* Faye: check if the random node is visited instead 
        *while(SetIn(V, (j = G->n*RandomUniform()))    
        */
        while(visited[(j = G->n*RandomUniform())])
		; // must terminate since k <= G->n
	    outbound[nOut++] = j;
	    j = 0;
#else
	    static int depth;
	    depth++;
	    // must terminate eventually as long as there's at least one connected component with >=k nodes.
	    assert(depth < MAX_TRIES); // graph is too disconnected
	    V = SampleGraphletNodeBasedExpansion(V, Varray, G, k, whichCC);
	    depth--;
	    // Ensure the damn thing really *is* connected.
	    TINY_GRAPH *T = TinyGraphAlloc(k);
	    TinyGraphInducedFromGraph(T, G, Varray);
	    assert(NumReachableNodes(T,0) == k);
	    TinyGraphFree(T);
	    return V;
#endif
	}
	else
	    j = nOut * RandomUniform();
	v1 = outbound[j];
	SetDelete(outSet, v1);
	SetAdd(V, v1); Varray[i] = v1;
	outbound[j] = outbound[--nOut];	// nuke v1 from the list of outbound by moving the last one to its place
	for(j=0; j<G->degree[v1];j++) // another loop over neighbors that may take a long time...
	{
	    v2 = G->neighbor[v1][j];
        /* Faye: check if it's invisted instead
        * if(!SetIn(outSet, v2) && !SetIn(V, v2)) */
        if (!visited[v2]) {
		    SetAdd(outSet, (outbound[nOut++] = v2));
	        visited[v2] = 1;
        }
    }
    }
    assert(i==k);
#if PARANOID_ASSERTS
    assert(SetCardinality(V) == k);
    assert(nOut == SetCardinality(outSet));
#endif
    return V;
}

