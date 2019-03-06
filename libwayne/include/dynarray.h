#include "misc.h"

/* Dynamic arrays: arrays that expand to fill the size needed
*/
typedef struct _dynamicArray {
    foint *array;
    int size, maxSize;
} ARRAY;

ARRAY *ArrayAlloc(int initNumObjs);
void ArrayFree(ARRAY *A);
// Add an element to the end of the dynArray.
foint ArrayAdd(ARRAY *A, foint new);
// Change an item at a specified position.
foint ArraySet(ARRAY *A, int pos, foint new);
// Get an item at specified position from the dynarray.
foint ArrayAt(ARRAY *A, int pos);
// Remove item at a specified postion. 
foint ArrayRemoveAt(ARRAY *A, int pos);
// Remove the first occurence of a specified element.
foint ArrayRemove(ARRAY *A, foint elem, pCmpFcn cmp);
// Remove all occurences of a specified element.
foint ArrayRemoveAll(ARRAY *A, foint elem, pCmpFcn cmp);
// Append B to the end of A.
ARRAY *ArrayAppend(ARRAY *C, ARRAY *A, ARRAY *B);
int ArraySize(ARRAY *A);

