#ifndef ESCAPE_CLIQUES_BRUTE_FORCE_H_
#define ESCAPE_CLIQUES_BRUTE_FORCE_H_

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

void cliques_brute_force(CGraph &CG, CDAG &DG, ofstream &of, string gname, string fname, unsigned long long clique_size)
{    
    high_resolution_clock::time_point t_elim_s = high_resolution_clock::now();
    
    Stack* stack = newStack();
    VertexSet *path = NULL, *nbrs = NULL;
    double num_cliques = 0;

    for (VertexIdx c=0; c<CG.nVertices; c++)
    {
        path = newVertexSet(1);
        path->vertices[0] = c;
        nbrs = newVertexSet(DG.outlist.offsets[c+1] - DG.outlist.offsets[c]);
        std::copy(DG.outlist.nbors+DG.outlist.offsets[c], DG.outlist.nbors+DG.outlist.offsets[c+1], nbrs->vertices);
        PartialClique* pc = newPartialClique(path, nbrs);
        stack->push(pc);

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
    }

    delStack(stack);

    high_resolution_clock::time_point t_elim_e = high_resolution_clock::now();
    auto duration_program = std::chrono::duration_cast<std::chrono::microseconds>( t_elim_e - t_elim_s ).count();
    
    cout << "time reqd: " << duration_program << endl;
	cout << "number of cliques = " << num_cliques << endl;
	cout << "clique size = " << clique_size << endl;

	of << clique_size << "," << num_cliques << "," << duration_program << "," << endl;
}

#endif
