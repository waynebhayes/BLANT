#define MAX_DIMENSIONS  8

typedef struct _checked_array
{
    int dim;    /* number of dimensions */
    int upper[MAX_DIMENSIONS]; /* upper[dim] of upper array bounds */
    void *array;    /* the actual array */
} CHECKED_ARRAY;

/* Use the C sizeof() operator for the first argument, then the
** number of dimentions, followed by the length in each dimension
*/
CHECKED_ARRAY *ArrayAlloc(int elementSize, int dim, ...);
void *ArrayAccess(CHECKED_ARRAY *, ...);
