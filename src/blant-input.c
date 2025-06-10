#include <stdbool.h>
#include <stdio.h>
#include "blant.h"
#include "blant-input.h"
#include "graph.h"
#include "misc.h"
#include "sets.h"
#include "bintree.h"

GRAPH* ReadGraph(int argc, char *argv[], int *optind, char supportNodeNames, Boolean GRAPH_GEN, Boolean weighted, Boolean useComplement) {
    FILE *fpGraph;
    int piped = 0;
    if(!argv[*optind]) {
        fpGraph = stdin;
        if(isatty(0) && !_quiet) Warning("reading graph input file from terminal, press ^D to finish");
    } else {
        char *graphFileName = argv[*optind];
        // readFile doesn't actually READ anything, it only prepares the file for reading by, for example, uncompressing it if it's name ends in ".gz"
        fpGraph = readFile(graphFileName, &piped);
        if(!fpGraph) Fatal("cannot open graph input file '%s'\n", argv[*optind]);
        (*optind)++;
    }
    assert(*optind == argc || GRAPH_GEN || _windowSampleMethod == WINDOW_SAMPLE_DEG_MAX);

    // Read network using native Graph routine.
    GRAPH *G = GraphReadEdgeList(fpGraph, SPARSE, supportNodeNames, weighted);
    if(useComplement) G->useComplement = true;

    if(supportNodeNames) {
        assert(G->name);
        _nodeNames = G->name;
    }
    if(fpGraph != stdin) closeFile(fpGraph, &piped);

    return G;
}
