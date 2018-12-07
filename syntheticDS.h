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

typedef struct hops{
	int valid, lower, upper;
} Hops;

int dictionary_create(Dictionary* dictionary);
int dictionary_get(Dictionary* dictionary, int key, int default_value);
void dictionary_set(Dictionary* dictionary, int key, int value);
KeyValue* getIterator(Dictionary* dictionary);
int getNext(KeyValue** iterator, int* key, int* value);
int create_stack(RevertStack* stack, int size);
int init_stack(RevertStack* stack);
int push(RevertStack* stack, Change elt);
int pop(RevertStack* stack, Change* elt);
int compare_ints(const void* a, const void* b);
int compare_doubles(const void* a, const void* b);
int getIntMedian(int* nums, int start, int end);
double getDoubleMedian(double* nums, int start, int end);
double PoissonDistribution(double l, int k);
double getDoubleBinSize(int n, double localClustCoff[n], double* scratchspace);
int getIntegerBinSize(int n, int GDVcolumn[n], int* scratchspace);
int getRandomNode(GRAPH* G, int src, int hops);
void sampleKHop(GRAPH* G, Dictionary* khop, double quality);
void print_khop(Dictionary* khop);