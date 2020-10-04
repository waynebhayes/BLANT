#ifndef BLANT_KOVACS_H
#define BLANT_KOVACS_H

#include "../blant.h"

extern int _kovacsOrdinal;
extern int _kovacsOrbit1, _kovacsOrbit2; // the columns (ie orbits) you want
extern float **_KovacsScore, **_KovacsNorm; // will be allocated lower triangle of n by n matrix of node-pair counts
extern int _KovacsMotifCount[MAX_CANONICALS][maxK][maxK]; //pre-computed matrices of 64-bit ints take only about 6MB
extern float _KovacsMotifNorm[MAX_CANONICALS][maxK][maxK]; 
extern unsigned _kovacsOrbitPairSeen[MAX_CANONICALS][maxK][maxK];
extern unsigned _kovacsOrbitPairEdge[MAX_CANONICALS][maxK][maxK];

void PreComputeKovacs(TINY_GRAPH *g, int topOrdinal, TINY_GRAPH *gTop, char perm[maxK]);
void ProcessKovacsNorm(TINY_GRAPH *g, GRAPH *G, unsigned Varray[]);
void ProcessKovacsPreComputed(TINY_GRAPH *g, unsigned Varray[]);
void ProcessKovacsAllOrbits(TINY_GRAPH *g, GRAPH *G, unsigned Varray[]);

#endif