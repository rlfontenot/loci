/*
  memory.h
  defines memory allocation routines that include out of memory checking.
*/
#ifndef MEMORY_HEADER

#define MEMORY_HEADER

#include <malloc.h>
#include "const.h"

#ifdef ansi
void *ecalloc(int) ;
#else
void *ecalloc() ;
#endif

#define ALLOC(size) ecalloc(size)
#define NEW(type) ((type *) ALLOC(sizeof(type)))
#define FREE(ptr)  free(ptr)

#define NEWSTR(str) strcpy(ALLOC(strlen(str)+1),str)
#endif
