/*
** mem-debug.h: note that memory tracking is available regardless of whether
** DEBUG is defined or not; it's dependent only upon whether or not you
** call EnableMemoryTracking().  DEBUG just gives you the chance for
** more verbose output.
**
** If you explicitly want memory tracking NEVER to be on, you must define
** NDEBUG.
*/

#ifndef _DEBUG_H_
#define _DEBUG_H_ /* to ensure we don't include it more than once */


/*
** Fatal (FATAL): Display an error message and terminate the program.  If
** you are debugging, make sure you compile with WITH the constant
** DEBUG defined.  Doing so will cause and "assert(false)" statement to
** be inserted before the program is terminated by calling "exit(1)".
** This allows the debugger to catch your program before it exits,
** allowing you to debug the problem.  If DEBUG is not defined, the program
** will simply terminate when this function is called (i.e., the program
** will terminate before the debugger has a chance to catch the error).
**
** IMPORTANT NOTE: You should not call this function directly.  Use the
** macro FATAL instead, because FATAL fills in the 'file' and 'line'
** parameters for you.
*/
void Fatal_fl( const char *message, const char *file, const int line );
#define FATAL(message) Fatal_fl( (message), __FILE__ , __LINE__ )


/*
** EnableMemoryTracking:  Turn on memory tracking and reporting for memory
** blocks allocated by Malloc (MALLOC) and released by Free (FREE).  To
** use memory tracking, this function must be called once (and only once)
** at the beginning of your program BEFORE any memory is allocated.
**
** IMPORTANT NOTES:
** 1. Only memory blocks allocated with Malloc (MALLOC) and released with
**    Free (FREE) are tracked and reported.
** 2. If you wish to use the memory tracking and reporting features of
**    Malloc (MALLOC), Free (FREE), and MemoryAllocationReport, you must call
**    this function ONCE and ONLY ONCE at the beginning of your program
**    BEFORE any memory is allocated.
*/
#ifndef NDEBUG
void EnableMemoryTracking( const char *file, const int line );
#define ENABLE_MEMORY_TRACKING() EnableMemoryTracking( __FILE__, __LINE__)
#else
#define ENABLE_MEMORY_TRACKING() /* nothing */
#endif

/*
** Malloc (MALLOC):  Allocate and return a pointer to a block of 'size'
** bytes of memory.  If the allocation fails (due to an out of memory
** condition), the program is terminated.  If memory tracking and
** reporting is enabled (by calling EnableMemoryTracking before any calls
** to Malloc), a list of allocated memory blocks will be maintained.  If
** you fail to release allocated memory blocks before your program
** terminates, you will receive an error message saying 'Memory leaks
** detected' when your program exits, followed by a list of the allocated
** blocks, their memory address, their size in bytes, and the location in
** the source code where they were allocated (i.e., the file name and line
** number).  You can also manually view this list at any time by calling
** calling the function MemoryAllocationReport.
**
** IMPORTANT NOTES:
** 1. You should not call this function directly.  Instead use the macro MALLOC
**    which fills in the 'file' and 'line' parameters for you.
** 2. All memory allocated by Malloc (i.e., MALLOC) must be freed by
**    calling Free (i.e., FREE).
** 3. DO NOT mix calls to Malloc and Free with calls to the standard C library
**    functions malloc and free.
*/
#ifndef NDEBUG
void *Malloc_fl( const int size, const char *file, const int line );  /* Do not call this function directly.  Use MALLOC instead. */
#define MALLOC(size) Malloc_fl( (size), __FILE__, __LINE__ )
#else
#define MALLOC(size) malloc(size)
#endif


/*
** Free (FREE):  Release a block of memory allocated by Malloc (i.e., MALLOC).
**
** IMPORTANT NOTES:
** 1. You should not call this function directly.  Use the macro FREE instead.
** 2. All memory allocated by Malloc (i.e., MALLOC) must be freed by
**    calling Free (i.e., FREE).
** 3. DO NOT mix calls to Malloc and Free with calls to the standard C
**    library functions malloc and free.
*/
#ifndef NDEBUG
void Free_fl( void *memoryBlock, const char *file, const int line );
#define FREE(pointer) Free_fl( (pointer), __FILE__, __LINE__ )
#else
#define FREE(p) free(p)
#endif


/*
** MemoryAllocationReport:  Displays a list of all memory blocks allocated
** by Malloc (MALLOC) that have not been deallocated by calling Free
** (FREE).  The list displays the following information about each memory
** block: the memory address of the block, the size of the memory block in
** bytes, and the location (i.e., the name of source file and the line
** number in that file) where the blocks was allocated.  The function
** EnableMemoryTracking must be called BEFORE any calls to Malloc (MALLOC)
** for MemoryAllocationReport to function correctly.
**
** IMPORTANT NOTES:
** 1. Only memory blocks allocated with Malloc (MALLOC) and released with
**    Free (FREE) are tracked and reported.
** 2. If you call this function without first enabling memory tracking and
**    reporting (i.e., by calling EnableMemoryTracking), then the function
**    will not display any output.
*/
#ifndef NDEBUG
void MemoryAllocationReport( const char *file, const int line );
#define MEMORY_ALLOCATION_REPORT() MemoryAllocationReport(__FILE__,__LINE__)
#else
#define MEMORY_ALLOCATION_REPORT() /* nothing */
#endif


/*
** DEBUGPRINTF:  Execute a printf statement if in debugging mode (i.e., if
** the constant DEBUG is defined).  If DEBUG is not defined, the
** DEBUGPRINTF statement is automatically skipped by the compiler (i.e.,
** you do not need to delete it from the source code).
**
** NOTE:  Do not call 'debugprintf' directly, call DEBUGPRINTF instead.
*/
#ifdef DEBUG
#define DEBUGPRINTF printf
#else
#define DEBUGPRINTF debugprintf
#endif
void debugprintf( const char *format, ... ); /* do not call this function directly, call DEBUGPRINTF instead */
#define TRACE DEBUGPRINTF

/*
** HERE:  a simple macro that prints the current soure file name and line
** number when debugging (i.e., when the constant DEBUG is defined).
** If DEBUG is not defined, the HERE statement is automatically skipped by
** the compiler (i.e., you do not have to removed it from the code).  This
** macro can be very useful for determining where a program fails.
*/
#ifdef DEBUG
#define HERE fprintf( stderr, "[%s:%d]", __FILE__, __LINE__ )
#else
#define HERE ((void *)0)
#endif

#endif /* _DEBUG_H_ */
