#include "Escape/ErrorCode.h"

#ifndef GRAPHIO_H__
#define GRAPHIO_H__

//Some utilities for loading graph data

#include "Escape/ErrorCode.h"
#include "Escape/Graph.h"

#include <string>
#include <vector>


namespace Escape
{

enum class IOFormat
{
  none      //try to guess from file extension or probing
  , escape  //our own internal format
  , snap    //Stanford SNAP
  , matrix  //Matrix Market
  , bcsr    //Binary CSR
};


//Loads a graph from one of the file formats we support. The returned
//graph should be freed with delGraph when you are done with it.
ErrorCode loadGraph(const char *path, Graph& graph, int undirected, IOFormat fmt);

ErrorCode loadGraph_p(const char *path, Graph& graph, int undirected, IOFormat fmt, float p);

}
#endif
