#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "uthash.h"
#include "graph.h"
#include "sets.h"

typedef struct keyvalue{
	int key, value;
	UT_hash_handle hh;
} KeyValue;

typedef struct dictionary{
	KeyValue* hashTable;
} Dictionary;

typedef struct change{
    int k, linenum, original, new;
} Change;

typedef struct revertstack{
    int tos, size;
    Change* space;
} RevertStack;

typedef struct int_tuple{
	int nums[2];
} IntTuple;

typedef struct smallworld{
	int make;  // +1 for make more small world, 0 for do nothing, -1 for make less small world
	int khop_interval;
	int lowerHops, upperHops;
} SmallWorld;

typedef struct clustcoffstate{
	int* histograms[2];  // stores num of elts in a bin
	int histograms_size[2];  // stores num of bins
	double histograms_bin_size[2];  // stores bin size/width
} ClustCoffState;

typedef struct gkstate{
	long udotv, sq_length_u, sq_length_v;
} GKState;

// key:value hash helpers
int dictionary_create(Dictionary* dictionary);
int dictionary_get(Dictionary* dictionary, int key, int default_value);
void dictionary_set(Dictionary* dictionary, int key, int value);
KeyValue* getIterator(Dictionary* dictionary);
int getNext(KeyValue** iterator, int* key, int* value);

// revert stack helpers
int create_stack(RevertStack* stack, int size);
int init_stack(RevertStack* stack);
int push(RevertStack* stack, Change elt);
int pop(RevertStack* stack, Change* elt);

// math helpers
int getIntMedian(int* nums, int start, int end);
double getDoubleMedian(double* nums, int start, int end);
double PoissonDistribution(double l, int k);

// optimum histogram bin size
double getDoubleBinSize(int n, double localClustCoff[n], double* scratchspace);
int getIntegerBinSize(int n, int GDVcolumn[n], int* scratchspace);

// node selection helpers
int getRandomNodeAtHops(GRAPH* G, int src, int hops);
int getRandomConnectedNode(GRAPH* G, int src);

// k-hop helpers
void sampleKHop(GRAPH* G, Dictionary* khop, double quality, int nodesBySp[G->n]);
int compareKHopByMedian(Dictionary* khop[2], int medians[2], int MAX_Keys[2]);
void print_khop_sample(Dictionary* khop);
