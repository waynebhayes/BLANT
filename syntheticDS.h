#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "uthash.h"

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

int dictionary_create(Dictionary* dictionary);
int dictionary_get(Dictionary* dictionary, int key, int default_value);
void dictionary_set(Dictionary* dictionary, int key, int value);
KeyValue* getIterator(Dictionary* dictionary);
int getNext(KeyValue** iterator, int* key, int* value);
int create_stack(RevertStack* stack, int size);
int init_stack(RevertStack* stack);
int push(RevertStack* stack, Change elt);
<<<<<<< HEAD
int pop(RevertStack* stack, Change* elt);
int compare_ints(const void* a, const void* b);
int getMedian(int* nums, int start, int end);
=======
int pop(RevertStack* stack, Change* elt);
>>>>>>> b5043c2c403bdcff6e100cc40eee02f2ad41e0d5
