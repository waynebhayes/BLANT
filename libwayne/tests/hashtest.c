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

int main(int argc, char *argv[])
{
    FILE *fp = stdin;
    if(argc==2) fp=fopen(argv[1],"r");
    assert(fp);
    int key;
    char *value, line[BUFSIZ];
    map_t H = hashmap_new();
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

    /*
     * Iteratively call f with argument (item, data) for
     * each element data in the hashmap. The function must
     * return a map status code. If it returns anything other
     * than MAP_OK the traversal is terminated. f must
     * not reenter any hashmap functions, or deadlock may arise.
     */
    //extern int hashmap_iterate(map_t in, PFany f, any_t item);

    int hash_status;
    while(fgets(line,sizeof(line),stdin) && 1==sscanf(line, "%d", &key)) {
	hash_status = hashmap_get(H, key, &value);
	if(hash_status == MAP_MISSING)
	    printf("No such element %d\n", key);
	else if(hash_status == MAP_OK)
	    printf("key %d gave <%s>\n", key, value);
	else assert(0);
    }
}

#if 0

/*
 * Get an element from the hashmap. Return MAP_OK or MAP_MISSING.
 */
extern int hashmap_get(map_t in, int key, any_t *arg);

/*
 * Remove an element from the hashmap. Return MAP_OK or MAP_MISSING.
 */
extern int hashmap_remove(map_t in, int key);

/*
 * Get any element. Return MAP_OK or MAP_MISSING.
 * remove - should the element be removed from the hashmap
 */
extern int hashmap_get_one(map_t in, any_t *arg, int remove);

/*
 * Free the hashmap
 */
extern void hashmap_free(map_t in);

/*
 * Get the current size of a hashmap
 */
extern int hashmap_length(map_t in);

#endif 
