#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <ctype.h>
#include "hashmap.h"
/*
#define MAP_MISSING -3
#define MAP_FULL -2
#define MAP_OMEM -1
#define MAP_OK 0
*/

int PrintHashEntry(int key, any_t data)
{
    char *s=(char*)data;
    printf("entry %d is \"%s\"\n", key, s);
    return MAP_OK;
}

int main(int argc, char *argv[])
{
    FILE *fp = stdin;
    if(argc==2) fp=fopen(argv[1],"r");
    assert(fp);
    int key;
    char *value, line[BUFSIZ];
    hashmap_t H = hashmap_new();
    while(fgets(line,sizeof(line),fp) && strlen(line)>1){
	assert(1==sscanf(line, "%d", &key));
	line[strlen(line)-1]='\0';
	value=line;
	while(isspace(*value))++value;
	while(isdigit(*value))++value;
	while(isspace(*value))++value;
	//printf("Inserting key %d value <%s>\n",key,value);
	hashmap_put(H, key, strdup(value));
    }

    printf("Finished entering. Now iterate through all elements:\n");
    hashmap_iterate(H, PrintHashEntry);
    printf("Finished enumerating elements\n");

    int hash_status;
    while(fgets(line,sizeof(line),stdin) && 1==sscanf(line, "%d", &key)) {
	hash_status = hashmap_get(H, key, (void**)&value);
	if(hash_status == MAP_MISSING)
	    printf("No such element %d\n", key);
	else if(hash_status == MAP_OK)
	    printf("key %d gave <%s>\n", key, value);
	else assert(0);
    }
    return 0;
}

#if 0

/*
 * Get an element from the hashmap. Return MAP_OK or MAP_MISSING.
 */
extern int hashmap_get(hashmap_t in, int key, any_t *arg);

/*
 * Remove an element from the hashmap. Return MAP_OK or MAP_MISSING.
 */
extern int hashmap_remove(hashmap_t in, int key);

/*
 * Get any element. Return MAP_OK or MAP_MISSING.
 * remove - should the element be removed from the hashmap
 */
extern int hashmap_get_one(hashmap_t in, any_t *arg, int remove);

/*
 * Free the hashmap
 */
extern void hashmap_free(hashmap_t in);

/*
 * Get the current size of a hashmap
 */
extern int hashmap_length(hashmap_t in);

#endif
