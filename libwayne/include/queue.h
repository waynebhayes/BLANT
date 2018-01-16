/* array implementation of queues of foints, from Lewis & Denenberg, p 77. */

#ifndef _QUEUE_H
#define _QUEUE_H
#include "misc.h"   /* for foints */

typedef struct _queueStruct {
    int maxSize, front, length;
    foint *queue;  /* will be an array[maxSize] */
} QUEUE;

/* maxSize is only an estimate; it will be dynamically increased as necessary */
QUEUE *QueueAlloc(int maxSize);

void QueueFree(QUEUE *q);

foint QueueFront(QUEUE *q); /* peek at front */
foint QueueGet(QUEUE *q);       /* get and delete front of queue */

/* returns element pushed */
foint QueuePut(QUEUE *q, foint new);

/* peek at the element n below the front */
foint QueueBelowTop(QUEUE *q, int n);

int QueueSize(QUEUE *q);

#endif  /* _QUEUE_H */
