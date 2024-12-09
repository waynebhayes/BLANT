#ifndef BLANT_PREDICT_H
#define BLANT_PREDICT_H

void Predict_Init(GRAPH *G);
void Predict_Shutdown(GRAPH *G);
void Predict_Flush(GRAPH *G);
void Predict_ProcessLine(GRAPH *G, unsigned long lineNum, char line[]);
void Predict_ProcessGraphlet(GRAPH *G, unsigned Varray[], TINY_GRAPH *g, Gint_type Gint, Gordinal_type GintOrdinal);
int  Predict_Merge(GRAPH *G);
extern Gint_type smaller_canon_map(Gint_type num, int k, unsigned char* return_permutation);
extern Gordinal_type canon_to_ordinal(Gint_type canon, int k);
Gordinal_type L_K_Func(Gint_type Gint);

#endif
