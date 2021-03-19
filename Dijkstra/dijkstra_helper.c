#include "../libwayne/include/misc.h"
#include "../libwayne/include/sets.h"
#include "../libwayne/include/graph.h"
#include "../libwayne/include/queue.h"
#include "../libwayne/include/rand48.h"
#include "../libwayne/include/Oalloc.h"
#include <assert.h>
#include <string.h>

int num_edges_back_to_subgraph(GRAPH *graph, int nodeid, int *aligned_nodes, int len)
{
    int edges = 0;
    int count = 0;
    while(count < len)
    {
        if (GraphAreConnected(graph, nodeid, aligned_nodes[count])) edges = edges + 1;
        count = count + 1;
    }
    return edges;
}

int num_edge_pairs_back_to_subgraph(GRAPH *g1, GRAPH *g2, int g1node, int g2node, int *aligned_pairs, int len)
{
    int edgepairs = 0;
    int count = 0;
    while(count < len)
    {
        if (GraphAreConnected(g1, g1node, aligned_pairs[count]) && GraphAreConnected(g2, g2node, aligned_pairs[count + 1])) edgepairs = edgepairs + 1;
        count = count + 2;
    }
    return edgepairs;
}

// this script seeks to replace the following functions
/*
def num_edges_back_to_subgraph(graph, node, aligned_nodes, trace=False):
    edges = 0
    for neighbor_node in aligned_nodes:
        if libCalc.GraphAreConnected(graph, node, neighbor_node):
            # if trace:
            #     print("Nnode in graph : " + str(neighbor_node))
            edges += 1
    return edges

def num_edge_pairs_back_to_subgraph(g1, g2, g1node, g2node, aligned_pairs):
    edgepairs = 0
    for n1, n2 in aligned_pairs:
        if libCalc.GraphAreConnected(g1, g1node, n1) and libCalc.GraphAreConnected(g2, g2node, n2):
            edgepairs += 1
    return edgepairs
*/