#ifndef HASH_HEADER
#define HASH_HEADER
#include "const.h"

#define DEFAULT_HASH_SIZE (7)
#define BIG_HASH_SIZE (12)
#define SMALL_HASH_SIZE (5)

struct hashents {
    char *name ;
    struct hashents *next ;
    void *entry ;
} ;

struct hashtable {
    int hashsize ;
    struct hashents *h_table[1] ;
} ;

#ifdef ansi
struct hashtable *makehash(int) ;
void             *findhashent(struct hashtable *, char *) ;
void              addptrtohash(struct hashtable *, char *, void *) ;
void              printhash(struct hashtable *) ;
void              operateOnHash(struct hashtable *,void (*)(char *,void *)) ;
#else
struct hashtable *makehash() ;
void             *findhashent() ;
void              addptrtohash() ;
void              printhash() ;
void              operateOnHash() ;

#endif

#endif
