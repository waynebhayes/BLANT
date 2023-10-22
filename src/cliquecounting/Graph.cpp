#include "Escape/Graph.h"
#include "Escape/JointSort.h"
#include <algorithm>
#include <map>
#include <unordered_set>
#include <vector>


using namespace Escape;
using namespace std;

//Basic binary search procedure
// Input: pointer array, index of last entry end, and val to search for
// Output: index if val is found, and -1 otherwise
VertexIdx Escape::binarySearch(EdgeIdx* array, VertexIdx end, EdgeIdx val)
{
    VertexIdx low = 0;
    VertexIdx high = end-1;
    VertexIdx mid;

    while (low <= high)
    {
        mid = (low+high)/2;
        if (array[mid] == val)
            return mid;
        if (array[mid] > val)
            high = mid-1;
        if (array[mid] < val)
            low = mid+1;
    }
    return -1;
}

// comparator that only compares the first in pair
bool Escape::pairCompareFirst(Pair firstPair, Pair nextPair)
{
    return firstPair.first < nextPair.first;
}

// comparator that compares the second in pair. If they are equal, then compare first in pair
bool Escape::pairCompareSecond(Pair firstPair, Pair nextPair)
{
    if (firstPair.second != nextPair.second)
        return firstPair.second < nextPair.second;
    return firstPair.first < nextPair.first;
}


Graph Escape::newGraph(VertexIdx nVertices, EdgeIdx nEdges)
{
  return {nVertices, nEdges, new VertexIdx[nEdges], new VertexIdx[nEdges]};
}


Graph Graph::copy() const
{
  Graph ret = newGraph(nVertices, nEdges);
  std::copy(srcs, srcs + nEdges, ret.srcs);
  std::copy(dsts, dsts + nEdges, ret.dsts);
  return ret;
}


void Graph::print(FILE* f) const
{
  fprintf(f, "nVertices=%ld nEdges=%ld\n", (int64_t) nVertices, (int64_t) nEdges);
  for (EdgeIdx i = 0; i < nEdges; ++i)
    fprintf(f, "  %ld   %ld -> %ld\n", (int64_t) i, (int64_t) srcs[i], (int64_t) dsts[i]);
}


void Escape::delGraph(Graph g)
{
  delete[] g.srcs;
  delete[] g.dsts;
}


CGraph Escape::newCGraph(VertexIdx nVertices, EdgeIdx nEdges)
{
  return {nVertices, nEdges, new EdgeIdx[nVertices + 1], new VertexIdx[nEdges]};
}


CGraph CGraph::copy() const
{
  CGraph ret = newCGraph(nVertices, nEdges);
  std::copy(offsets, offsets + nVertices + 1, ret.offsets);
  std::copy(nbors, nbors + nEdges, ret.nbors);
  return ret;
}


void CGraph::print(FILE *f) const
{
  fprintf(f, "nVertices=%ld nEdges=%ld\n", (int64_t) nVertices, (int64_t) nEdges);
  for (VertexIdx i = 0; i < nVertices; ++i)
  {
    fprintf(f, "%ld: ", (int64_t) i);
    for (EdgeIdx j = offsets[i]; j < offsets[i + 1]; ++j)
      fprintf(f, "%ld ", (int64_t) nbors[j]);
    fprintf(f, "\n");
  }
}


void Escape::delCGraph(CGraph g)
{
  delete[] g.offsets;
  delete[] g.nbors;
}


// Checks if edge (v1, v2) is present in CGraph
// If edge is present: Return index (in nbors) of v2 in neighbors of v1
// If edge is not present: Return -1
int CGraph::isEdge(VertexIdx v1, VertexIdx v2)
{
    if (v1 >= nVertices)
        return -1;
    for (EdgeIdx i=offsets[v1]; i < offsets[v1+1]; ++i)
        if (nbors[i] == v2)
            return i;
    return -1;
}

//Checks if edge (v1, v2) is present using binary search
//Assumes CGraph is sorted by ID
//If edge is present: Return index (in nbors) of v2 in neighbors of v2
//If edge is not present: return -1
EdgeIdx CGraph::getEdgeBinary(VertexIdx v1, VertexIdx v2) const
{
    if (v1 >= nVertices)
        return -1;
    EdgeIdx low = offsets[v1];
    EdgeIdx high = offsets[v1+1]-1;
    EdgeIdx mid;

    while(low <= high)
    {
        mid = (low+high)/2;
        if (nbors[mid] == v2)
            return mid;
        if (nbors[mid] > v2)
            high = mid-1;
        else
            low = mid+1;
    }
    return -1;
}

//Checks if edge (v1, v2) is present using binary search
//Assumes CGraph is sorted by ID
//If edge is present: return true
//If edge is not present: return false
bool CGraph::isEdgeBinary(VertexIdx v1, VertexIdx v2) const
{
    VertexIdx deg1 = offsets[v1+1] - offsets[v1];
    VertexIdx deg2 = offsets[v2+1] - offsets[v2];

//     if(deg2 < deg1)
//     {
//         VertexIdx swp = v1;
//         v1 = v2;
//         v2 = swp;
//     }
// 
    EdgeIdx low = offsets[v1];
    EdgeIdx high = offsets[v1+1]-1;
    EdgeIdx mid;

    while(low <= high)
    {
        mid = (low+high)/2;

        if (nbors[mid] == v2)
            return true;
        if (nbors[mid] > v2)
            high = mid-1;
        else
            low = mid+1;
    }
    return false;
}

//Same as above, changing name since we retrieve the edge index as well.
//There is a good reason to have both isEdge and getEdge, since it's easier
//to answer the first question.  Recommend changing isEdge to return bool - vv
EdgeIdx CGraph::getEdge(VertexIdx v1, VertexIdx v2) const
{
    if (v1 >= nVertices)
        return -1;
    for (EdgeIdx i=offsets[v1]; i < offsets[v1+1]; ++i)
        if (nbors[i] == v2)
            return i;
    return -1;
}

VertexIdx CGraph::randomNeighbor(VertexIdx v) const
{
	VertexIdx nbr = nbors[offsets[v] + (rand() % (offsets[v+1] - offsets[v]))];  	
    return nbr;
}

// This sorts each individual adjacency list by vertex ID. This is useful for
// doing a binary search, or for merging neighbor lists to find common neighbors.
void CGraph::sortById() const
{
    for (VertexIdx i=0; i < nVertices; i++)
        std::sort(nbors+offsets[i],nbors+offsets[i+1]);
}

// This outputs a new, isomorphic CGraph where vertex labels are in increasing order corresponding to degree.
// Thus, (after the relabeling), for all i < j, the degree of i is less than that of j.

CGraph CGraph::renameByDegreeOrder() const
{
    CGraph ret = newCGraph(nVertices, nEdges);
    Pair *deg_info = new Pair[nVertices];

    VertexIdx *mapping = new VertexIdx[nVertices];
    VertexIdx *inverse = new VertexIdx[nVertices];

    
    // Construct array of pairs, storing old vertex label and degree
    for (VertexIdx i=0; i < nVertices; i++)
    {
        deg_info[i].first = i;
        deg_info[i].second = offsets[i+1]-offsets[i];
    }

    // sort the pairs by degree (if degree is same, sort by old vertex label)
    std::sort(deg_info,deg_info+nVertices,Escape::pairCompareSecond);

    // Construct the mapping of old vertex label to new vertex label
    // So mapping[i] is what i is mapped to
    // And inverse[i] is what maps to i
    for (VertexIdx i=0; i < nVertices; i++)
    {
        mapping[deg_info[i].first] = i;
        inverse[i] = deg_info[i].first;
    }

    // Initialize offsets of output CGraph
    ret.offsets[0] = 0;
    EdgeIdx current = 0;


    // Loop over new vertices
    for (VertexIdx new_label=0; new_label < nVertices; new_label++)
    {
        VertexIdx old_label = inverse[new_label]; // corresponding old label for new vertices
        for (EdgeIdx pos = offsets[old_label]; pos < offsets[old_label+1]; pos++) // loop over neighbors of old label
        {
            VertexIdx old_nbr = nbors[pos];
            VertexIdx new_nbr = mapping[old_nbr]; //corresponding new neighbor
            ret.nbors[current] = new_nbr; // insert new neighbor in nbors of output
            current++;
        }
        ret.offsets[new_label+1] = current; // all neighbors of new_label have been added, so we set offset for new_label+1
    }

    return ret;
}


CGraph Escape::makeCSR(Graph g, bool inPlace)
{
  //create a temporary graph for sorting unless the user requests in-place
  //operation.
  printf("In makeCSR. number of edges = %ld \n", g.nEdges);
  auto tmpG = inPlace ? g : g.copy();

  auto begin = JSIterator<VertexIdx, VertexIdx>{tmpG.srcs, tmpG.dsts};
  auto end   = begin + tmpG.nEdges;
  std::sort(begin, end);

  //We steal the dsts array from tmpG.
  CGraph ret = {g.nVertices, g.nEdges, new EdgeIdx[g.nVertices + 1], tmpG.dsts};

  //Now we have everything sorted by src, compress:
  VertexIdx cv = 0;
  for (EdgeIdx i = 0; i < tmpG.nEdges; ++i)
  {
    auto src = tmpG.srcs[i];
    while (cv <= src)
      ret.offsets[cv++] = i;
  }
  while (cv <= g.nVertices)
    ret.offsets[cv++] = g.nEdges;

  delete[] tmpG.srcs; //we retain tmpG.dsts in the output

  ret.sortById();

  return ret; 
}


CGraph Escape::makeCSC(Graph g, bool inPlace)
{
  return makeCSR({g.nVertices, g.nEdges, g.dsts, g.srcs}, inPlace);
}

double Escape::degeneracy(CGraph &g)
{
  // CGraph cg = g.copy();
  VertexIdx n = g.nVertices;
  VertexIdx min_deg = n;
  double degen = 0;

  map <VertexIdx, unordered_set<VertexIdx> > deg_list;    
  map <VertexIdx, bool> touched;
  vector<VertexIdx> cur_degs;
  cur_degs.resize(n);
	//deg_list.resize(40000);
	//printf("before for\n");
  for (int i=0; i<n; i++)
  {
    VertexIdx deg = g.degree(i);
    deg_list[deg].insert(i);
    cur_degs[i] = deg;
    if (min_deg > deg) min_deg = deg;
  }
	//printf("after for\n");
  for(int i=0; i<n; i++)
  {
	//printf("Inside first for\n");
    while(deg_list[min_deg].size() == 0)
      min_deg++;
    unordered_set<VertexIdx>::iterator it = deg_list[min_deg].begin();
    VertexIdx source = *it;
    touched[source] = true;
    if (degen < min_deg) 
    {
      degen = min_deg;
      // printf("degen = %f\n", degen);
    }
    deg_list[min_deg].erase(it);
   //printf("Calling second for\n"); 
    for (int j=g.offsets[source]; j<g.offsets[source+1]; j++)
    {
      VertexIdx nbr = g.nbors[j];
      if (touched[nbr] == true)
          continue;
      VertexIdx deg = cur_degs[nbr];
      deg_list[deg].erase(nbr);
      deg_list[deg-1].insert(nbr);
      if (deg-1 < min_deg)
          min_deg = deg-1;
      cur_degs[nbr] = deg - 1;
            
    }

  }
  printf("returning degen = %f\n", degen);
  return degen;
}

double Escape::countWedges(CGraph &g)
{
  double wedges = 0;
  for (VertexIdx i = 0; i < g.nVertices; i++)
  {
    double deg = g.degree(i);
    wedges += deg * (deg - 1) / 2;
  }
  return wedges;
}
