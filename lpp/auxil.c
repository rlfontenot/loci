#include <stdio.h>
#include "hash.h"

extern FILE* yyin ;
struct hashtable* global_def ;

main(int argc, char* argv[])
{
   ++argv, --argc ;
   if(argc > 0)
     yyin = fopen(argv[0], "r") ;
   else
     yyin = stdin ;
   
   global_def = makehash(DEFAULT_HASH_SIZE) ;
   yyparse();
   return 0;
}

yyerror(msg)
   char *msg;
{
   extern long yypos;

   printf("line %d: %s\n", yypos, msg);
   exit(1);
}

yywrap()
{
   return 1;
}
