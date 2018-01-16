#include "misc.h"

/* Dynamic arrays: arrays that expand to fill the size needed
*/
typedef struct _dynamicArray {
    foint *array;
    int maxSize, objSize;
} ARRAY;

ARRAY *ArrayAlloc(int objSize, int initNumObjs);
void *ArrayWrite(ARRAY *, int element, foint value);
void *ArrayReadElement(ARRAY *, int element);
