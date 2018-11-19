#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "uthash.h"
#include "syntheticDS.h"

// DICTIONARY
// wrapper around 'uthash'
// https://raw.githubusercontent.com/troydhanson/uthash/master/src/uthash.h

int dictionary_create(Dictionary* dictionary){
	dictionary->hashTable = NULL;
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
		HASH_ADD_INT(dictionary->hashTable, key, s);
	}else{
		// overwrite
		s->value = value;
	}
}

// iterator

KeyValue* getIterator(Dictionary* dictionary){
	KeyValue* iterator = dictionary->hashTable;
	return iterator;
}

int getNext(KeyValue** iterator, int* key, int* value){
	if((iterator==NULL) || (*iterator == NULL))
		return -1;
	memcpy(key, &((*iterator)->key), sizeof(int));
	memcpy(value, &((*iterator)->value), sizeof(int));
	*iterator = (*iterator)->hh.next;
	return 0;
}


// STACK
int create_stack(RevertStack* stack, int size){
    stack->tos = -1;
    stack->size = size;
    stack->space = (Change*) malloc(size * sizeof(Change));
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
