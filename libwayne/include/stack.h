#ifndef _STACK_H
#define _STACK_H

/* array implementation of stacks of ints, from Lewis & Denenberg, p 76. */

#include "misc.h"

typedef struct _stackStruct {
    int maxSize;
    foint *stack;  /* will be an array[size+1], with sizeof stack in stack[0] */
} STACK;

STACK *StackAlloc(int maxSize);
void StackFree(STACK *s);
foint StackTop(STACK *s);

/* also returns old top */
foint StackPop(STACK *s);

/* returns element pushed */
foint StackPush(STACK *s, foint new);

/* returns the element n items below the top.  StackBelowTop(s, 0)
** gives the same result as StackTop(s)
*/
foint StackBelowTop(STACK *s, int n);
int StackSize(STACK *s);

#endif  /* _STACK_H */
