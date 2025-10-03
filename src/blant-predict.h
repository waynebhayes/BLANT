// This software is part of github.com/waynebhayes/BLANT, and is Copyright(C) Wayne B. Hayes 2025, under the GNU LGPL 3.0
// (GNU Lesser General Public License, version 3, 2007), a copy of which is contained at the top of the repo.
#ifndef BLANT_PREDICT_H
#define BLANT_PREDICT_H

void Predict_Init(GRAPH *G); // must be called once before any prediction starts
void Predict_ProcessGraphlet(GRAPH *G, unsigned Varray[], TINY_GRAPH *g, Gint_type Gint, Gordinal_type GintOrdinal);
void Predict_Flush(GRAPH *G); // call this near the end to flush outputs
void Predict_Shutdown(GRAPH *G); // should be called at end to free memory (to avoid valgrind warnings)
void Predict_ProcessLine(GRAPH *G, unsigned long lineNum, char line[]); // parent calls this to process outputline from child
int  Predict_Merge(GRAPH *G);
extern Gint_type Predict_canon_map(Gint_type num, int k, unsigned char* return_permutation);
extern Gordinal_type Predict_canon_to_ordinal(Gint_type canon, int k);
Gordinal_type L_K_Func(Gint_type Gint);

#endif
