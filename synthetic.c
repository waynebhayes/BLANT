#include <sys/file.h>
#include <unistd.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "misc.h"
#include "tinygraph.h"
#include "graph.h"
#include "blant.h"
#include "sets.h"

int main(int argc, char *argv[])
{
    int i, opt, numSamples=0;

    if(argc == 1)
    {
	fprintf(stderr, "%s\n", USAGE);
	exit(1);
    }

    _k = 0;

    while((opt = getopt(argc, argv, "s:k:r:")) != -1)
    {
	switch(opt)
	{
	case 'r': _seed = atoi(optarg);
	    break;
	case 's': numSamples = atoi(optarg);
	    if(numSamples < 0) Fatal("numSamples must be non-negative\n%s", USAGE);
	    break;
	case 'k': _k = atoi(optarg);
	    if(!(3 <= _k && _k <= 8)) Fatal("k must be between 3 and 8\n%s", USAGE);
	    break;
	default: Fatal("unknown option %c\n%s", opt, USAGE);
	}
    }

    if(!argv[optind]) Fatal("no input graph file specified\n%s", USAGE);
    char *graphFileName = argv[optind];
    FILE *fpGraph = fopen(argv[optind], "r");
    if(!fpGraph) Fatal("cannot open graph input file '%s'\n", argv[optind]);
    optind++;
    assert(optind == argc);

    SetGlobalCanonMaps(); // needs _k to be set

    // Read it in using native Graph routine.
    GRAPH *G = GraphReadEdgeList(fpGraph, true, _supportNodeNames); // sparse=true
    if(_supportNodeNames)
    {
	assert(G->name);
	_nodeNames = G->name;
    }
    fclose(fpGraph);
    return RunBlantInThreads(_k, numSamples, G);
}
