#include <stdio.h>
#include <string.h>
#include "misc.h"
#include "stack.h"

int main(void)
{
    STACK *lineStack = StackAlloc(1000);	/* grows automagically */
    char buf[102400];

    while(fgets(buf, sizeof(buf), stdin))
	StackPush(lineStack, (foint)(void*)strdup(buf));

    while(StackSize(lineStack))
	fputs((char*)StackPop(lineStack).v, stdout);
    return 0;
}
