/*
** mem-debug.c
*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "misc.h"
#include "mem-debug.h"


void Fatal_fl( const char *message, const char *file, const int line )
{
    fprintf( stderr, "[%s:%d] Error: %s\n", file, line, message );

#ifdef DEBUG
    assert( false );
#endif

    exit( 1 );
}


/*
** debugprintf:  empty function to allow removal of DEBUGPRINTF statements
** when debugging is turned off by undefining the constant DEBUG.
**
** Do not call this function directly.  Use the macro DEBUGPRINTF
** instead.
*/
void debugprintf( const char *format, ... ) {}


/*
** MEMORY_BLOCK_HEADER: structure used by Malloc (MALLOC) and Free (FREE)
** to facilitate memory tracking and reporting.
*/
#define ALIGN_PADDING 2
typedef struct _memoryBlockHeader {
    struct _memoryBlockHeader *next; /* pointer to next memory block */
    struct _memoryBlockHeader *prev; /* pointer to previous memory block */
    const char *file; /* pointer to constant string containing the name of
			 the source file in which the block was allocated */
    int line;	/* line number file where the block was allcoated */
    int size;          /* size of the user portion of the block, in bytes */
    const char *fileFree; /* pointer to constant string containing the name of
			 the source file in which the block was freed */
    int lineFree;
    double align[ALIGN_PADDING];/*force alignement of the block and to provide
    				** padding between this header and the user
				** portion of the memory block */
} MEMORY_BLOCK_HEADER;


/*
** memoryBlockList: list of MEMORY_BLOCK_HEADER structures for all blocks
** of memory currently allocated by Malloc.
*/
static MEMORY_BLOCK_HEADER memoryBlockList = {
    NULL,
    NULL,
    __FILE__,
    __LINE__,
    0
};


/*
** Total number of memory blocks currently allocated.
*/
static int memoryBlockCount = 0;


/*
** Total number of bytes of memory currently allocated.
*/
static int memoryTotalUsage = 0;


/*
** true if memory tracking and reporting is enabled, false otherwise.
*/
static int memoryTrackingEnabled = false;


void *Malloc_fl( const int size, const char *file, const int line )
{
    MEMORY_BLOCK_HEADER *memoryBlockHeader;
    void *memoryBlock;
    int totalSize;

    assert( size >= 0 );
    totalSize=size + (memoryTrackingEnabled ? sizeof(MEMORY_BLOCK_HEADER) : 0);

    memoryBlock = malloc( totalSize );
    if( memoryBlock == NULL )
	Fatal_fl( "out of memory", file, line );

    if( memoryTrackingEnabled )
    {
	int *align;
	memoryBlockHeader = (MEMORY_BLOCK_HEADER *)memoryBlock;
	memoryBlock = (void *)( (unsigned char *)memoryBlockHeader + sizeof(MEMORY_BLOCK_HEADER) );

	memoryBlockHeader->file = file;
	memoryBlockHeader->line = line;
	memoryBlockHeader->size = size;
	memoryBlockHeader->fileFree = "(not freed yet)";
	memoryBlockHeader->lineFree = 0;

	/* use DEADBEEF as a flag that this is actually malloc'd memory */
	assert(sizeof(double) == 2*sizeof(int));
	align = (int*)memoryBlockHeader->align;
	align[0] = 0xDEAD;
	align[1] = 0xBEEF;

	memoryBlockHeader->prev = &memoryBlockList;
	memoryBlockHeader->next = memoryBlockList.next;
	if( memoryBlockHeader->next != NULL )
		memoryBlockHeader->next->prev = memoryBlockHeader;
	memoryBlockList.next = memoryBlockHeader;

	assert( memoryTotalUsage >= 0 );
	memoryTotalUsage += size;

	assert( memoryBlockCount >= 0 );
	memoryBlockCount++;
    }

    return memoryBlock;
}


void Free_fl( void *memoryBlock, const char *file, const int line )
{
    MEMORY_BLOCK_HEADER *memoryBlockHeader;

    if( memoryBlock == NULL )
        Fatal_fl( "attempt to free NULL pointer", file, line );

    if( memoryTrackingEnabled )
    {
	int *align;
	memoryBlockHeader = (MEMORY_BLOCK_HEADER*)
	    ((void*)
		((unsigned char*)memoryBlock - sizeof(MEMORY_BLOCK_HEADER)));

	/* verify that this was actually malloc'd */
	align = (int*)memoryBlockHeader->align;
	if(align[0] == 0xBEEF && align[1] == 0xDEAD)
	{
	    fprintf(stderr, "memory allocated at %s line %d was already freed at %s line %d; ",
		memoryBlockHeader->file, memoryBlockHeader->line,
		memoryBlockHeader->fileFree, memoryBlockHeader->lineFree);
	    Fatal_fl("Tried to free a second time at %s line %d", file, line);
	}
	if(align[0] != 0xDEAD || align[1] != 0xBEEF)
	    Fatal_fl("call to FREE on memory which is either non-MALLOC'd, already freed and re-used, or corrupt", file, line);

	/* mark this as FREE'd */
	align[0] = 0xBEEF;
	align[1] = 0xDEAD;
	memoryBlockHeader->fileFree = file;
	memoryBlockHeader->lineFree = line;

	memoryBlockHeader->prev->next = memoryBlockHeader->next;
	if( memoryBlockHeader->next != NULL )
	    memoryBlockHeader->next->prev = memoryBlockHeader->prev;

	memoryBlockCount--;
	assert( memoryBlockCount >= 0 );

	memoryTotalUsage -= memoryBlockHeader->size;
	assert( memoryTotalUsage >= 0 );

	free( memoryBlockHeader );
    }
    else
	free( memoryBlock );
}


void MemoryAllocationReport( const char *file, const int line )
{
    MEMORY_BLOCK_HEADER *memoryBlockHeader = memoryBlockList.next;

    /* do nothing if memory tracking and leak detection is not enabled */
    if( !memoryTrackingEnabled )
	return;

    fprintf( stderr, "MemoryAllocationReport from file %s, line %d:\n",
	file, line);

    if( memoryBlockCount > 0 )
    {
	fprintf( stderr, "Total leaked memory: %d byte%s in %d block%s.\n\n",
	    memoryTotalUsage, ((memoryTotalUsage > 1)?"s":""),
	    memoryBlockCount, ((memoryBlockCount > 1)?"s":""));
	fprintf( stderr, "Address:  Bytes:      Line:  File:\n" );
	fprintf( stderr, "--------  ----------  -----  ----------\n" );

	while( memoryBlockHeader != NULL )
	{
	    fprintf( stderr, "%8p  %10.1d  %5.1d  %s\n",
		(void *)((unsigned char*)memoryBlockHeader +
		    sizeof(MEMORY_BLOCK_HEADER) ),
		memoryBlockHeader->size,
		memoryBlockHeader->line,
		memoryBlockHeader->file );

	    memoryBlockHeader = memoryBlockHeader->next;
	}
    }
    else
	fprintf( stderr, "No memory is currently allocated.\n" );
}


/*
** MemoryLeakDetect: detect memory leaks (if any) when your program terminates.
** If memory leaks are detected, MemoryAllocationReport is called to display a
** list of allocated memory blocks.  After displaying the list, the function
** released all memory blocks before the program terminates.
**
** IMPORTANT NOTE:
** Under no circumstances should you call this function directly!  It is
** "registered" by EnableMemoryTracking and is called automatically when
** the program terminates.
*/
static void MemoryLeakDetect( void )
{
    MEMORY_BLOCK_HEADER *memoryBlockHeader;

    if( !memoryTrackingEnabled )
	FATAL( "memory reporting and leak detection is not enabled");

    if( memoryBlockCount > 0 || memoryBlockList.next != NULL )
    {
	fprintf( stderr, "---\nYour program has exited, and has leaked the following memory:\n" );
	MemoryAllocationReport(__FILE__, __LINE__);
	while( memoryBlockList.next != NULL )
	{
	    memoryBlockHeader = memoryBlockList.next;
	    memoryBlockList.next = memoryBlockHeader->next;
	    free( (void *)memoryBlockHeader );
	}
    }
}


void EnableMemoryTracking( const char *file, const int line )
{
    memoryTrackingEnabled = true;
    if( atexit( MemoryLeakDetect  ) != 0 )
	Fatal_fl( "could not initialize memory tracking and leak detection",
	    file, line);
}
