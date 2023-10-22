#ifndef ESCAPE_GRAPH_H_
#define ESCAPE_GRAPH_H_

#include <cstdlib>
#include <cstdio>

namespace Escape
{

using VertexIdx = int64_t;
using EdgeIdx   = int64_t;
using Count     = int64_t;

const EdgeIdx invalidEdge = -1;


//Graph in COO representation.  These are shallow POD objects and can be
//copied freely. Memory management is always done explicitly.
struct Graph
{
  VertexIdx  nVertices; //number of vertices in the graph
  EdgeIdx    nEdges;    //number of edges in this list
  VertexIdx *srcs;      //array of source vertices
  VertexIdx *dsts;      //array of destination vertices

  //Make a copy
  Graph copy() const;

  //For debugging.  Not for serialization.
  void print(FILE* f = stdout) const;
};

//Allocate memory to hold a graph.
Graph newGraph(VertexIdx nVertices, EdgeIdx nEdges);

//Release memory associated with the graph.  Should have been allocated
//with newGraph.
void delGraph(Graph g);

// structure holding pair of VertexIdx
struct Pair
{
    VertexIdx first;
    VertexIdx second;
};

//Basic binary search procedure
// Input: pointer array, index of last entry end, and val to search for
// Output: index if val is found, and -1 otherwise
VertexIdx binarySearch(EdgeIdx* array, VertexIdx end, EdgeIdx val);


// comparator that only compares the first in pair
bool pairCompareFirst(Pair firstPair, Pair nextPair);

// comparator that compares the second in pair. If they are equal, then compare first in pair
bool pairCompareSecond(Pair firstPair, Pair nextPair);


//Graph in CSR/CSC representation.
struct CGraph
{
  VertexIdx  nVertices;   //number of vertices in the graph
  EdgeIdx    nEdges;      //number of edges in the graph
  EdgeIdx   *offsets;     //vertex v's edges are [offset[v], offset[v + 1])
  VertexIdx *nbors;       //incoming or outgoing neighbors

  //Make a copy
  CGraph copy() const;

  //For debugging, not for serialization
  void print(FILE* f = stdout) const;

  // Checks if edge (v1, v2) is present
  int isEdge(VertexIdx v1, VertexIdx v2);

  void sortById() const;

  CGraph renameByDegreeOrder() const; 

  //Returns the index of the edge v1 -> v2 in the nbor list nbors.
  //Returns invalidEdge if v1 -> v2 does not exist
  EdgeIdx getEdgeBinary(VertexIdx v1, VertexIdx v2) const;

  bool isEdgeBinary(VertexIdx v1, VertexIdx v2) const;
  //Returns the index of the edge v1 -> v2 in the nbor list nbors.
  //Returns invalidEdge if v1 -> v2 does not exist
  EdgeIdx getEdge(VertexIdx v1, VertexIdx v2) const;

  EdgeIdx degree(VertexIdx v) const { return offsets[v + 1] - offsets[v]; }
	
	// Adding functions for moment-estimator
	VertexIdx randomNeighbor(VertexIdx v) const;
};

//Allocate memory for a CSR/CSC graph
CGraph newCGraph(VertexIdx nVertices, EdgeIdx nEdges);

//Release memory associated with a CSR/CSC graph.  Should have been allocated
//with newCGraph
void delCGraph(CGraph g);

//Make a CSR graph from a COO graph.  If inPlace is true, the input graph is
//destroyed, i.e. you should not call delGraph on it.
CGraph makeCSR(Graph g, bool inPlace = false);

//Make a CSC graph from a COO graph.  If inPlace is true, the input graph is
//destroyed, i.e. you should not call delGraph on it.
CGraph makeCSC(Graph g, bool inPlace = false);

//Calculates the degeneracy of the graph
double degeneracy(CGraph &g);

//Counts the number of wedges in the graph - can be directed or undirected
double countWedges(CGraph &g);

}
#endif
