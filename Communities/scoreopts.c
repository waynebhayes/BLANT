#include "structs.h"
#include "measures.h"


typedef struct 	{
    const char * name;
    double (*func)(PARTITION * P, COMMUNITY * C, int fakeN);
    int optimizingForMax;
} SCORINGOPTION;

static const SCORINGOPTION scoreOpts[] = { 
    {"IntraEdgeDensity", IntraEdgeDensity, 1},
    {"InterEdgeDensity", InterEdgeDensity, -1},
    {"Conductance", Conductance, -1}, 
    {"HayesScore", HayesScore, 1}, 
    {"NewmanAndGirvan", NewmanAndGirvan, 1},
    {"Expansion", Expansion, -1},
    {"NormalizedCut", NormalizedCut, -1},
    {"CPM", CPM, 1},

};

static const int numScoreOpts = sizeof(scoreOpts) / sizeof(SCORINGOPTION);

