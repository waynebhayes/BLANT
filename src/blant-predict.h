#ifndef BLANT_PREDICT_H
#define BLANT_PREDICT_H

void Predict_Init(GRAPH *G);
void Predict_FlushMotifs(GRAPH *G);
void Predict_ProcessLine(GRAPH *G, char line[]);

void Predict_AccumulateMotifs(GRAPH *G, unsigned Varray[], TINY_GRAPH *g, Gint_type Gint, int GintOrdinal);
int  Predict_Merge(GRAPH *G);

#endif
