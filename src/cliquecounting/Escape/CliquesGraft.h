#ifndef ESCAPE_CLIQUES_EDGE_GRAFT_H_
#define ESCAPE_CLIQUES_EDGE_GRAFT_H_

#include <algorithm>
#include <chrono>
#include <random>
#include "Escape/ErrorCode.h"
#include "Escape/Graph.h"
#include "Escape/Digraph.h"
#include "Escape/Utils.h"
#include "Escape/nCr.h"
#include "Escape/CliqueHelper.h"
#include "JointSort.h"

using namespace Escape;
using namespace std;
using namespace std::chrono;

double get_edge_cliques(Graph &g, CGraph &CG, CDAG &DG, EdgeIdx edge, unsigned long long clique_size)
{
    double num_cliques = 0;
    VertexIdx v1 = g.srcs[edge];
    VertexIdx v2 = g.dsts[edge];
    
    Stack* stack = newStack();
    VertexSet *path = NULL, *nbrs = NULL;

    path = newVertexSet(0);
    nbrs = newVertexSet(CG.offsets[v1+1] - CG.offsets[v1]);

    int c = 0;
    for (int i=CG.offsets[v1]; i<CG.offsets[v1+1]; i++)
    {
        for (int j=CG.offsets[v2]; j<CG.offsets[v2+1]; j++)
        {
            if (CG.nbors[i] == CG.nbors[j])
            {
                nbrs->vertices[c] = CG.nbors[i];
                c++;
            }
        }
    }
    nbrs->nVertices = c;
    PartialClique *pc1 = newPartialClique(path, nbrs);
    stack->push(pc1);
    StackItem* si = stack->pop();

    while (si != NULL)
    {
        path = si->pc->path;
        nbrs = si->pc->nbrs;  // get the current path and its nbrs
        VertexIdx psize = path->nVertices, nsize = nbrs->nVertices;

        int flag = 0;
        if (psize + nsize < clique_size)
        {
        }
        else if (psize == clique_size)
        {   
            num_cliques++;
        }
        else if (psize == clique_size - 1)
        {
            num_cliques += nsize;
        }
        else if (psize == clique_size - 2)
        {
            for (int i=0; i<nsize-1; i++)
                for (int j=i+1; j<nsize; j++)
                    if (CG.isEdge(nbrs->vertices[i], nbrs->vertices[j]) != -1)
                        num_cliques++;
        }
        else
        {
            for (VertexIdx i=0; i<nsize; i++)
            {
                VertexIdx v = nbrs->vertices[i];
                VertexSet *new_path = newVertexSet(psize+1);

                int s = DG.outlist.offsets[v+1]-DG.outlist.offsets[v], l = nsize;
                if (l < s)
                {
                    l = s;
                    s = nsize;
                }

                VertexSet *common = newVertexSet(s);
                int k = 0, offset = DG.outlist.offsets[v];
                for (int l =0; l < nbrs->nVertices; l++) {
                    for (int j=0; j<DG.outlist.offsets[v+1]-DG.outlist.offsets[v]; j++) {
                        if (nbrs->vertices[l] == DG.outlist.nbors[offset + j])
                        {
                            common->vertices[k] = nbrs->vertices[l];
                            k++;
                            break;
                        }
                    }
                }
                if (k == 0) 
                {
                    delete[] common->vertices;
                    delete common;
                    delete[] new_path->vertices;
                    continue;
                }
                else
                {
                    common->nVertices = k; 
                    PartialClique *pc = newPartialClique(new_path, common);
                    stack->push(pc);
                }
            }
        }
        
        if (nbrs->nVertices > 0) delete[] nbrs->vertices;
        delete nbrs;
        if (path->nVertices > 0) delete[] path->vertices;
        delete path;
        delete si->pc;
        delete si;
        si = stack->pop();
    }

    delStack(stack);
    double ret = num_cliques; 
    return ret;
}

void cliques_graft(Graph &g, CGraph &CG, CDAG &DG, ofstream &of, string gname, string fname, unsigned long long clique_size, int num_samples)
{    
    high_resolution_clock::time_point t_elim_s = high_resolution_clock::now();
    Stack* stack = newStack();
    VertexSet *path = NULL, *nbrs = NULL;
 
    double num_cliques = 0, est = 0;

    std::random_device rd{}; // use to seed the rng 
    std::mt19937 rng{rd()}; // rng

    std::uniform_int_distribution<int> distribution(0,floor((CG.nEdges-1)/2));

	cout << "Starting sampling" << endl;
    for (int i=0; i<num_samples; i++)
    {
        EdgeIdx edge = distribution(rng);
        num_cliques += get_edge_cliques(g, CG, DG, 2*edge, clique_size - 2);
    }
    est = num_cliques * CG.nEdges/ (2 * num_samples * nCr[clique_size][2]);

	high_resolution_clock::time_point t_elim_e = high_resolution_clock::now();
    auto duration_program = std::chrono::duration_cast<std::chrono::microseconds>( t_elim_e - t_elim_s ).count();
    
    cout << "time reqd: " << duration_program << endl;
	cout << "number of cliques = " << num_cliques << endl;
    cout << "est = " << est << endl;
    
    of << clique_size << "," << est << "," << duration_program << "," << endl;
}

#endif
