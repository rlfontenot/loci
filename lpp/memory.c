/*
 * memory.c
 * contains error checking allocation routine.
 */

#include <stdio.h>
#include "memory.h"

void *ecalloc(size)
register int size ;
{
    register void *vp ;

    if(!(vp = (void *) malloc(size))) {
        fprintf(stderr,"Out of Memory\n") ;
        exit(-1) ;
    }
    return vp ;
}

