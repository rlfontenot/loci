/*
 * hash.c
 *  This file includes a set of simple hash table creation and lookup
 *  routines:
 *
 *  makehash(size)
 *   makehash returns a hashtable of size ~2^size, predefined arguments
 *   of DEFAULT_HASH_SIZE, BIG_HASH_SIZE and SMALL_HASH_SIZE may
 *   be used.
 *  findhashent(hashtable,name)
 *   findhash returns the void pointer that was previously stored in
 *   the hash under the character string pointed to by name.  If there
 *   is no entry under name then a NULL pointer is returned.
 *  addptrtohash(hashtable,name,pointer)
 *   addptrtohash installs a new entry into the hashtable.  This associates
 *   the string name with pointer.
 *  printhash()
 *   printhash outputs the current hashtable and names to stdout.
 *   Used for debugging purposes.
 */

#include <stdio.h>

#include "memory.h"
#include "debug.h"
#include "hash.h"

static int primes[] = { 1,2,3,7,13,31,61,127,251,509,1021,2039,4093 } ;

static int numprimes = sizeof(primes)/sizeof(int) ;
     
static struct hashents *
  makehashent(hashentname)
char *hashentname ;
{
    register struct hashents *hep ;
    
    fatal(!hashentname) ;
    hep = NEW(struct hashents) ;
    hep->name = (char *) ALLOC(strlen(hashentname)+1) ;
    strcpy(hep->name,hashentname) ;
    hep->next = NULL ;
    hep->entry = NULL ;
    return hep ;
}

static int
  hashfunc(name,currenthash)
char *name ;
struct hashtable *currenthash ;
{
    register int tot=0 ;
    register char *s = name ;
    register int hashsize ;
    
    fatal(!name) ;
    fatal(!currenthash) ;
    hashsize =currenthash->hashsize ;
    while(*s != '\0')
      tot = ((tot<<6)+((int)*s++))%hashsize ;
    tot = (tot<0)?-tot:tot ;
    return tot ;
}



/* make a new hashtable.  Size is the approximate size of the hash table log2.
   i.e. a for a hashtable with about 4000 entries size would be 12.  If size is
   less than zero then the DEFAULT_HASH_SIZE will be used. */

struct hashtable *
makehash(size)
int size ;
{
    register struct hashtable *hp ;
    register int i ;
    int hashsize ;
    
    size = (size<0)?DEFAULT_HASH_SIZE:size ;
    hashsize = primes[(size<numprimes-1)?size:numprimes-1] ;
    hp = (struct hashtable *) ALLOC(sizeof(*hp)+
                                    sizeof(struct hashents *)*(hashsize-1)) ;
    hp->hashsize = hashsize ;
    for(i=0;i<hashsize;i++)
      hp->h_table[i] = NULL ;
    return hp ;
}

/* Get the data for hash table entry (name) in hashtable pointed to by hp. */

void *
findhashent(hp,name)
struct hashtable *hp ;
char *name ;
{
    register struct hashents *hep ;
    register int hindex ;
    
    fatal(!name) ;
    fatal(!hp) ;
    hindex = hashfunc(name,hp) ;
    for(hep = hp->h_table[hindex];hep;hep=hep->next) {
        fatal(!hep->name) ;
        if(!strcmp(name,hep->name))
          break ;
    }
    return hep?hep->entry:NULL ;
}


/* Add an entry to the hash under (name) */
void
addptrtohash(hp,name,ptr)
struct hashtable *hp ;
char *name ;
void *ptr ;
{
    register struct hashents *hep = makehashent(name) ;
    register int hindex ;
    
    
    fatal(!name) ;
    fatal(!ptr) ;
    fatal(!hp) ;
    warn(findhashent(hp,name)) ;
    
    hep->entry = ptr ;
    hindex = hashfunc(hep->name,hp) ;
    hep->next = hp->h_table[hindex] ;
    hp->h_table[hindex] = hep ;
}

/* Display the hashtable */
void
printhash(hp)
struct hashtable *hp ;
{
    int i ;
    struct hashents *hep ;

    fatal(!hp) ;
    
    for(i=0;i<hp->hashsize;i++) {
        printf("[%d] ",i) ;
        if(hp->h_table[i])
          for(hep=hp->h_table[i];hep;hep=hep->next)
            printf(" %s",hep->name) ;
        else
          printf("(nil)") ;
        printf("\n") ;
    }
}


void
operateOnHash(hp,fp)
struct hashtable *hp ;
#ifdef ansi
void (*fp)(char *,void *) ;
#else
int (*fp)() ;
#endif
{
    int i ;
    struct hashents *hep ;
    
    fatal(!hp) ;
    fatal(!fp) ;
    
    for(i=0;i<hp->hashsize;i++)
      if(hp->h_table[i])
        for(hep=hp->h_table[i];hep;hep=hep->next)
          (*fp)(hep->name,hep->entry) ;
}
