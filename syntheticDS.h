#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "uthash.h"

typedef struct keyvalue{
	int key, value;
	UT_hash_handle hh;
} KeyValue;

typedef struct storage{
	struct storage* next;
	KeyValue node;
} Storage;

typedef struct dictionary{
	// int:int as key:value pairs
	KeyValue* hashTable;
	Storage* storage_head;
	Storage* storage_tail;
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
Storage* getKeyValueIterator(Dictionary* dictionary);
void getNext(Storage* iterator, int* key, int* value);
int create_stack(RevertStack* stack, int size);
int init_stack(RevertStack* stack);
int push(RevertStack* stack, Change elt);
int pop(RevertStack* stack, Change* elt);




