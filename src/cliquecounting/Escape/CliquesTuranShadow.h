#ifndef ESCAPE_CLIQUES_H_
#define ESCAPE_CLIQUES_H_

#include "Escape/CliqueHelper.h"
#include <cassert>
#include <vector>
#include <random>
#include <ctgmath>
#include <iomanip>

using namespace Escape;
using namespace std;

bool DEBUG = false;


/* calculates the density of the nbrs. also gets degree of every vertex, 
does degree ordering of the vertices, and gets child, grandchild pair 
*/
double get_density(CGraph CG, VertexSet *nbrs, EdgeIdx& num_edges)
{
    double n = nbrs->nVertices;
    double total = 0;
    double actual = 0;
    if (n==0) return 0;
    if (n==1) return 1;

    int *degrees = (int *) calloc(n, sizeof(int));

    for (int i=0; i<n-1; i++)
    {
        for (int j=i+1; j<n; j++)
        {
            total++;
            if (CG.isEdge(nbrs->vertices[i], nbrs->vertices[j]) != -1)
            {
                degrees[i]++;
                degrees[j]++;
                actual++;
            }
        }
    }

    double density = actual/total;

    VertexIdx tmp;

    for (int i=n-2; i>=0; i--)
    {
        for (int j=0; j<=i; j++)
        {
            if (degrees[j] > degrees[j+1])
            {
                tmp = degrees[j];
                degrees[j] = degrees[j+1];
                degrees[j+1] = tmp;
                tmp = nbrs->vertices[j];
                nbrs->vertices[j] = nbrs->vertices[j+1];
                nbrs->vertices[j+1] = tmp;
            }
        }
    }
    num_edges = actual;
	free(degrees);
    return density;
}

void cliques_turan_shadow(CGraph &CG, CDAG &DG, ofstream &of, string gname, string fname, int runs, int num_samples, int clique_size)
{    
    high_resolution_clock::time_point treegen_s = high_resolution_clock::now();

    double num_cliques = 0, num_leaves = 0;
    double non_leaf_cliques = 0;
    double discarded_leaves = 0;
    int k = 0; 
    double succ_ratio = 0, leaf_samples = 0, leaf_size = 0, density = 0, leaf_cliques = 0;
    VertexIdx psize = 0, nsize = 0;
    PartialClique *pc = NULL;
    VertexIdx arr[10000];  // array for storing numbers 0 through leaf_size, which will be used by the random number generator
    VertexIdx rng[10000]; // will be used for the random number generator
    double N = 0;    

    vector<int> leafk;
    vector<int> leafn;
    int leaves_inc = 10000000;
    int max_leaves = leaves_inc;
    leafk.resize(leaves_inc);
    leafn.resize(leaves_inc);
    double internal = 0;


    for (int i=0; i<10000; i++)
        arr[i] = i;

    double internal_nodes = 0;
	double shadow_size = 0;
    srand(time(NULL));
    
    Stack *stack = newStack();
    vector <VertexSet*> leaves;
    vector <double> leafExperiments;
    leafExperiments.resize(leaves_inc);
    leaves.resize(leaves_inc);

    for (VertexIdx i=0; i<CG.nVertices; i++)
    {
    	VertexSet *path = newVertexSet(1);
    	path->vertices[0] = i;
    	VertexSet *nbrs = newVertexSet(DG.outlist.offsets[i+1] - DG.outlist.offsets[i]);
    	std::copy(DG.outlist.nbors+DG.outlist.offsets[i], DG.outlist.nbors+DG.outlist.offsets[i+1], nbrs->vertices);
    	pc = newPartialClique(path, nbrs);
    	stack->push(pc);

        StackItem* si = stack->pop();
        int flag = 0;
    	EdgeIdx num_edges = 0;
        while (si != NULL)  // unraveling the stack
        {
            VertexSet *path = si->pc->path;
            VertexSet *nbrs = si->pc->nbrs;  // get the current path and its nbrs
            psize = path->nVertices;
            nsize = nbrs->nVertices;
            flag = 0;

    	    k = clique_size - psize;
            density = get_density(CG, nbrs,num_edges);
    				
            if (psize + nsize < clique_size) 
            {
                discarded_leaves++;
                flag = 1;
            }
            else if (psize == clique_size)
            {
                num_cliques++;
                non_leaf_cliques++;
                flag = 1;
            }
            else if (psize == clique_size - 1)
            {
                num_cliques += nbrs->nVertices;
                non_leaf_cliques += nbrs->nVertices;
                flag = 1;
            }
            else if (k == 2)
    		{
    				non_leaf_cliques += num_edges;
    				flag = 1;
    		}
    		else
    		{
                double thresh_dens = 1.0 - pow((k-1), -1);
                if (density >= thresh_dens)
                {
                    leaves[num_leaves] = nbrs;
                    leafExperiments[num_leaves] = nCr[nsize][k];
                    leafk[num_leaves] = k;
                    leafn[num_leaves] = nsize;
                    N += nCr[nsize][k];
                  
                    shadow_size += nsize;

                    num_leaves++;
                    if (num_leaves == max_leaves) 
                    {
                        max_leaves += leaves_inc;
                        leaves.resize(max_leaves);
                        leafExperiments.resize(max_leaves);
                        leafk.resize(max_leaves);
                        leafn.resize(max_leaves);
                    }
                }
                else
                {
                    internal++;
                    for (int i=0; i<nbrs->nVertices-1; i++)
                    {
                        VertexSet *new_path = newVertexSet(1);
                        new_path->nVertices = path->nVertices + 1;
                        VertexSet *new_nbrs = newVertexSet(nsize);

                        int c = 0;
                        for (int j=i+1; j<nbrs->nVertices; j++)
                        {
                            if (CG.isEdge(nbrs->vertices[i], nbrs->vertices[j]) != -1)
                            {
                                new_nbrs->vertices[c] = nbrs->vertices[j];
                                c++;
                            }
                        }
                        new_nbrs->nVertices = c;
    			
                        PartialClique *pc1 = newPartialClique(new_path, new_nbrs);
                        stack->push(pc1);
                    }

                    flag = 1;
                }
            }
            if (flag == 1)
            {
                if (nbrs->nVertices >= 0) delete[] nbrs->vertices;
                delete nbrs;
                
            }
            if (path->nVertices > 0) delete[] path->vertices;
            delete path;
            delete si->pc;
            delete si;
            si = stack->pop();
            continue;
        }

    }
	high_resolution_clock::time_point treegen_e = high_resolution_clock::now();
    auto duration_treegen = std::chrono::duration_cast<std::chrono::microseconds>( treegen_e - treegen_s ).count();

    cout << "internal = " << internal << endl;
    cout << "leaves = " << num_leaves << endl;
    cout << "discarded_leaves = " << discarded_leaves << endl;
    cout << "internal + leaves + discarded_leaves = " << internal + num_leaves + discarded_leaves << endl;

    string datfile(gname + "_" + to_string(clique_size) + "_" + to_string(num_samples) + "_" + to_string(runs) + "_data");
    std::string dfname = "../results/cliques/" + datfile;
    std::cout << dfname << std::endl;
    std::ofstream df;
    df.open(dfname);
    if (!df.is_open())
    {
        std::cout << "Could not open output file." << std::endl;
        exit(1);
    }
    printf("Output File: %s\n", dfname.c_str());
    df << "estimate,succ_ratio,time," << endl;
    df.precision(2);
    double est;

	if (num_leaves != 0)
	{
		std::default_random_engine generator;
		std::discrete_distribution<int> distribution (leafExperiments.begin(), leafExperiments.end());
	
    	for (int r=0; r<runs; r++)
    	{
        
        	high_resolution_clock::time_point sampling_s = high_resolution_clock::now();
        	auto duration_sampling = std::chrono::duration_cast<std::chrono::microseconds>(sampling_s - sampling_s).count();
        	leaf_cliques = 0;
        	for (int i=0; i<num_samples; ++i) 
        	{
            	int leaf = distribution(generator);
            	k = leafk[leaf];
            	VertexSet *l = leaves[leaf];
            	std::copy(arr, arr+l->nVertices, rng);
            	std::random_shuffle ( rng, rng+l->nVertices);
            	VertexIdx *sample = (VertexIdx *) calloc(k,sizeof(VertexIdx)); 
            	for (int j = 0; j < k; j++)
                	sample[j] = l->vertices[rng[j]];

            	if (isClique(sample, k, CG)) 
            	{
                	leaf_cliques++; //experiment was a success. Sampled set of vertices do form a clique
            	}

            	free(sample);
        	}

        	high_resolution_clock::time_point sampling_e = high_resolution_clock::now();
        	duration_sampling = std::chrono::duration_cast<std::chrono::microseconds>( sampling_e - sampling_s ).count();
        	modf((leaf_cliques * N / num_samples)+non_leaf_cliques, &est);
        	succ_ratio = leaf_cliques / num_samples;
        	df << est << "," <<  succ_ratio << "," << duration_sampling << "," << endl;
    	}
    }
	else 
     	df << "0,0,0," << endl;
    		
	df.close();

    string outfile(gname + "_" + to_string(clique_size) + "_" + to_string(num_samples) + "_params");
    std::string ofname = "../results/cliques/" + outfile;
    std::cout << ofname << std::endl;
    std::ofstream ofparams;
    ofparams.open(ofname);
    if (!ofparams.is_open())
    {
        std::cout << "Could not open output file." << std::endl;
        exit(1);
    }
    printf("Output File: %s\n", ofname.c_str());

    ofparams << "internal," << "leaves," << "discarded leaves,"  << "shadow size," << "N,";
    ofparams << "tree_t," << "num_samples," << endl;
    ofparams << internal << "," << num_leaves << "," << discarded_leaves << "," << shadow_size << "," << N << ",";
    ofparams << duration_treegen << "," << num_samples << "," << endl;
    ofparams << endl;
    ofparams.close();
   
	for (int i=0; i<num_leaves; i++)
	{
		if (leaves[i]->nVertices >= 0) delete[] leaves[i]->vertices;
            delete leaves[i];
            
	} 
    cout << "clique_size = " << clique_size << endl;;
    cout << "Num of leaves = " << num_leaves << endl;
    cout << "Num of discarded leaves = " << discarded_leaves << endl;
    cout << "Number of non_leaf_cliques = " << non_leaf_cliques << endl;
    cout << "Estimated number of cliques = " << est << endl;
    cout << "N = " << N << endl;
    cout << "time reqd to build tree: " << duration_treegen << endl;
}

#endif
