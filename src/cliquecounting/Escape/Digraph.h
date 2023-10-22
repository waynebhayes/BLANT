#ifndef ESCAPE_DIGRAPH_H_
#define ESCAPE_DIGRAPH_H_

#include "Escape/ErrorCode.h"
#include "Escape/Graph.h"
#include <algorithm>

using namespace Escape;

// DAG structure has two pointers, one to the 
// adjacency list of outedge, one to the adjacency list of inedges
struct CDAG
{
    CGraph outlist;
    CGraph inlist;
};

void delCDAG(CDAG DAG)
{
  delCGraph(DAG.outlist);
  delCGraph(DAG.inlist);
}

// Structure for comparing nodes according to their degree.
// So u < v if degree of u less than that of v in graph g.

struct DegreeComp
{
    CGraph *g;
    DegreeComp(CGraph *g) { this->g = g;}

    bool operator () (VertexIdx u, VertexIdx v)
    {
        VertexIdx degu = g->offsets[u+1] - g->offsets[u];  // Degree of u
        VertexIdx degv = g->offsets[v+1] - g->offsets[v];  // Degree of v
    
        if (degu < degv || (degu == degv && u < v))    // Comparing degrees and breaking ties by id
            return true;
        else
            return false;
    }
};

//Construct DAG based on degree ordering
//
// Input: Pointer for CGraph g
// Output: CDAG for degree ordering in g
//         This is the CDAG for the DAG where each edge points from lower degree endpoint to higher degree endpoint.
//
//
//         The outlist in CDAG is guaranteed to be sorted by degrees. This means that the neighbors
//         of every vertex in the outlist are sorted by their degrees in g. This is quite useful in
//         further processing.
CDAG degreeOrdered(CGraph *g)
{
    CDAG ret;     // CDAG to be returned
    printf("In degreeOrdered. number of edges= %ld", g->nEdges);
    CGraph outdag = {g->nVertices, 0, new EdgeIdx[g->nVertices+1], new VertexIdx[g->nEdges+1]};  // Initialize DAG of out-edges
    CGraph indag = {g->nVertices, 0, new EdgeIdx[g->nVertices+1], new VertexIdx[g->nEdges+1]};   // Initialize DAG of in-edges
    EdgeIdx outcur = 0;
    EdgeIdx incur = 0;
    VertexIdx dest;
    VertexIdx degi;
    VertexIdx degdest;

    outdag.offsets[0] = 0;
    indag.offsets[0] = 0;
    for (VertexIdx i=0; i < g->nVertices; ++i)   // Looping over all vertices in g
    {
        for (EdgeIdx j = g->offsets[i]; j < g->offsets[i+1]; ++j)   // Looping over neighbors of i in g
        {
            dest = g->nbors[j];     // We are now looking at edge (i,dest)
            degi = g->offsets[i+1] - g->offsets[i];   // Degree of i
            degdest = g->offsets[dest+1]- g->offsets[dest];   // Degree of dest
            //printf("i=%ld dest=%ld degi=%ld degdest=%ld\n",i,dest,degi,degdest);

            //We now orient the edge depending of degi vs degdest.
            // We break ties according to vertex id.
            // In the output, the g-edge (i,dest) is either pointing to dest (in if condition) or pointing to i (in else condition).

            if (degi < degdest || (degi == degdest && i < dest))   
            {
                outdag.nbors[outcur] = dest;   // We want point edge from i to dest. So this directed edge is added to outdag.
                ++outcur;                      // Increment pointer in outdag.nbors and the number of edges in outdag.
                ++outdag.nEdges;
            }
            else
            {
                indag.nbors[incur] = dest;     // We point edge from dest to i. So this edge goes into indag.
                ++incur;                       // Pointer and number of edges incremented
                ++indag.nEdges;
            }
        }
        outdag.offsets[i+1] = outcur;         // We have finished all edges incident to i, so we can update offsets in DAGs.
        indag.offsets[i+1] = incur;
    }

    for (VertexIdx i=0; i < g->nVertices;++i)  // Loops over vertices
        std::sort(outdag.nbors+outdag.offsets[i], outdag.nbors+outdag.offsets[i+1], DegreeComp(g)); // In outdag, sort all neighbors of i according to their degree. Note that DegreeComp gives the desired comparator.

    ret.outlist = outdag;
    ret.inlist = indag;

    return ret;
}

#endif



