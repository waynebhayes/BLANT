#ifndef BLANT_INPUT_H
#define BLANT_INPUT_H

#include "graph.h"

GRAPH* ReadGraph(int argc, char *argv[], int *optind, char supportNodeNames, Boolean GRAPH_GEN, Boolean weighted, Boolean useComplement);

#endif 