#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "uthash.h"
#include "syntheticDS.h"

// DICTIONARY
// wrapper around uthash
// https://raw.githubusercontent.com/troydhanson/uthash/master/src/uthash.h
Storage* addToStore(Storage* head, Storage* tail, KeyValue* elt){
	if (tail == NULL){
		tail = (Storage*) malloc(sizeof(Storage));
		head = tail;
	}else{
		tail->next = (Storage*) malloc(sizeof(Storage));
		tail = tail->next;
	}
	memcpy(&(tail->node), elt, sizeof(KeyValue));
	tail->next = NULL;
	return tail;  // this pointer should not be pointed to another mem location
}

int dictionary_create(Dictionary* dictionary){
	dictionary->hashTable = NULL;
	dictionary->storage_head = NULL;
	dictionary->storage_tail = NULL;
}

int dictionary_get(Dictionary* dictionary, int key, int default_value){
	KeyValue* s;
	HASH_FIND_INT(dictionary->hashTable, &key, s);
	if (s == NULL){
		return default_value;
	}else{
		return s->value;
	}
}

void dictionary_set(Dictionary* dictionary, int key, int value){
	KeyValue* s;
	HASH_FIND_INT(dictionary->hashTable, &key, s);
	if (s == NULL){
		// no such key exists
		s = (KeyValue*) malloc(sizeof(KeyValue));
		s->key = key;
		s->value = value;
		KeyValue* t = addToStore(dictionary->storage_head, dictionary->storage_tail, s);
		HASH_ADD_INT(dictionary->hashTable, key, t);
	}else{
		// overwrite
		s->value = value;
	}
}

// iterator for key:value pairs
Storage* getKeyValueIterator(Dictionary* dictionary){
	Storage* iterator;
	if (dictionary->storage_head != NULL)
		memcpy(&iterator, &(dictionary->storage_head), sizeof(Storage*));
	return iterator;
}

void getNext(Storage* iterator, int* key, int* value){
	if(iterator == NULL)
		return;
	memcpy(key, iterator->node.key, sizeof(key));
	memcpy(value, iterator->node.value, sizeof(value));
	iterator = iterator->next;
}


// STACK
int create_stack(RevertStack* stack, int size){
    stack->tos = -1;
    stack->size = size;  // usually size = k*_numSamples (atmost every line in BLANT can change)
    stack->space = (Change*) Malloc(size * sizeof(Change));
    if (stack->space != NULL)
        return 0;
    else
        return -1;
}

int init_stack(RevertStack* stack){
    stack->tos = -1;
    if (stack->space != NULL)
        return 0;
    else
        return -1;
} 

int push(RevertStack* stack, Change elt){
    if (stack->tos == (stack->size-1))
        return -1;
    (stack->tos) += 1;
    memcpy(&(stack->space)[stack->tos], &elt, sizeof(Change));
    return 0;
}

int pop(RevertStack* stack, Change* elt){
    if (stack->tos == -1)
        return -1;
    memcpy(elt, &(stack->space)[stack->tos], sizeof(Change));
    stack->tos -= 1;
    return 0;
}
