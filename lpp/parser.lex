
D			[0-9]
L			[a-zA-Z_]
H			[a-fA-F0-9]
E			[Ee][+-]?{D}+
FS			(f|F|l|L)
IS			(u|U|l|L)*


/* the "incl" state is used for picking up the name
 * of an include file
 */
%x incl

%{
#include <stdio.h>
#include <strings.h>
#include "TYPE.h"
#include "yygrammar.h"

#define MAX_INCLUDE_DEPTH 10
YY_BUFFER_STATE include_stack[MAX_INCLUDE_DEPTH];
int include_stack_ptr = 0;

void count() ;
char* tmp ;
%}

%%
"/*"			{ comment() ; }
"//"                    { comment2() ; }

"auto"			{ count() ; return(AUTO) ; }
"bool"                  { count() ; return(BOOL) ; }
"break"			{ count() ; return(BREAK) ; }
"case"			{ count() ; return(CASE) ; }
"char"			{ count() ; return(CHAR) ; }
"const"			{ count() ; return(CONST) ; }
"continue"		{ count() ; return(CONTINUE) ; }
"default"		{ count() ; return(DEFAULT) ; }
"do"			{ count() ; return(DO) ; }
"double"		{ count() ; return(DOUBLE) ; }
"else"			{ count() ; return(ELSE) ; }
"enum"			{ count() ; return(ENUM) ; }
"extern"		{ count() ; return(EXTERN) ; }
"float"			{ count() ; return(FLOAT) ; }
"for"			{ count() ; return(FOR) ; }
"goto"			{ count() ; return(GOTO) ; }
"if"			{ count() ; return(IF) ; }
"int"			{ count() ; return(INT) ; }
"long"			{ count() ; return(LONG) ; }
"register"		{ count() ; return(REGISTER) ; }
"return"		{ count() ; return(RETURN) ; }
"short"			{ count() ; return(SHORT) ; }
"signed"		{ count() ; return(SIGNED) ; }
"sizeof"		{ count() ; return(SIZEOF) ; }
"static"		{ count() ; return(STATIC) ; }
"struct"		{ count() ; return(STRUCT) ; }
"switch"		{ count() ; return(SWITCH) ; }
"typedef"		{ count() ; return(TYPEDEF) ; }
"union"			{ count() ; return(UNION) ; }
"unsigned"		{ count() ; return(UNSIGNED) ; }
"void"			{ count() ; return(VOID) ; }
"volatile"		{ count() ; return(VOLATILE) ; }
"while"			{ count() ; return(WHILE) ; }

"datatypes"             { count() ; return(DATATYPES) ; }
"pointwise rule"        { count() ; return(POINTWISE_RULE) ; }
"unit rule"             { count() ; return(UNIT_RULE) ; }
"apply rule"            { count() ; return(APPLY_RULE) ; }
"singleton rule"        { count() ; return(SINGLETON_RULE) ; }
"$"                     { count() ; return(ITER_VAR_OP) ; }
"$="			{ count() ; return(JOIN) ; }
"conditional"           { count() ; return(CONDITIONAL) ; }
"$include"              { BEGIN(incl); }
<incl>[ \t]*            /* eat the whitespace */
<incl>[^ \t\n]+         { /* got the include file name */
        char* file_name = (char*)malloc(strlen(yytext)) ;

        if(include_stack_ptr >= MAX_INCLUDE_DEPTH)
        {
           fprintf(stderr, "Includes nested too deeply\n") ;
           exit(1) ;
        }

        include_stack[include_stack_ptr++] = YY_CURRENT_BUFFER ;

        strncpy(file_name, yytext+1, strlen(yytext)-2) ;
        yyin = fopen(file_name, "r") ;

        if(!yyin)
        {
           fprintf(stderr, "error open file: %s\n", file_name) ;
           exit(1) ;
        }

        yy_switch_to_buffer(yy_create_buffer(yyin,YY_BUF_SIZE)) ;

        BEGIN(INITIAL) ;
}

<<EOF>>                 {
        if(--include_stack_ptr < 0)
        {
           yyterminate() ;
        }

        else
        {
           yy_delete_buffer(YY_CURRENT_BUFFER) ;
           yy_switch_to_buffer(include_stack[include_stack_ptr]) ;
        }
}


{L}({L}|{D})*		{
  count() ;
  tmp = (char *)malloc(strlen(yytext)+1) ;
  strcpy(tmp, yytext) ;
  yylval.stringval = tmp ;
  return(check_type()) ;
}

0[xX]{H}+{IS}?		{ 
  count() ;
  tmp = (char *)malloc(strlen(yytext)+1) ;
  strcpy(tmp, yytext) ;
  yylval.stringval = tmp ;
  return(CONSTANT) ;
}

0{D}+{IS}?		{ 
  count() ; 
  tmp = (char *)malloc(strlen(yytext)+1) ;
  strcpy(tmp, yytext) ;
  yylval.stringval = tmp ;
  return(CONSTANT) ;
}

{D}+{IS}?		{ 
  count() ; 
  tmp = (char *)malloc(strlen(yytext)+1) ; 
  strcpy(tmp, yytext) ; 
  yylval.stringval = tmp ; 
  return(CONSTANT) ;
}

'(\\.|[^\\'])+'		{ 
  count() ; 
  tmp = (char *)malloc(strlen(yytext)+1) ; 
  strcpy(tmp, yytext) ; 
  yylval.stringval = tmp ; 
  return(CONSTANT) ;
}

{D}+{E}{FS}?		{ 
  count() ; 
  tmp = (char *)malloc(strlen(yytext)+1) ; 
  strcpy(tmp, yytext) ; 
  yylval.stringval = tmp ; 
  return(CONSTANT) ;
}

{D}*"."{D}+({E})?{FS}?	{ 
  count() ; 
  tmp = (char *)malloc(strlen(yytext)+1) ; 
  strcpy(tmp, yytext) ; 
  yylval.stringval = tmp ; 
  return(CONSTANT) ;
}

{D}+"."{D}*({E})?{FS}?	{ 
  count() ; 
  tmp = (char *)malloc(strlen(yytext)+1) ; 
  strcpy(tmp, yytext) ; 
  yylval.stringval = tmp ; 
  return(CONSTANT) ;
}

\"(\\.|[^\\"])*\"	{ 
  count() ; 
  tmp = (char *)malloc(strlen(yytext)+1) ; 
  strcpy(tmp, yytext) ; 
  yylval.stringval = tmp ; 
  return(STRING_LITERAL) ;
}

">>="			{ count() ; return(RIGHT_ASSIGN) ; }
"<<="			{ count() ; return(LEFT_ASSIGN) ; }
"+="			{ count() ; return(ADD_ASSIGN) ; }
"-="			{ count() ; return(SUB_ASSIGN) ; }
"*="			{ count() ; return(MUL_ASSIGN) ; }
"/="			{ count() ; return(DIV_ASSIGN) ; }
"%="			{ count() ; return(MOD_ASSIGN) ; }
"&="			{ count() ; return(AND_ASSIGN) ; }
"^="			{ count() ; return(XOR_ASSIGN) ; }
"|="			{ count() ; return(OR_ASSIGN) ; }
">>"			{ count() ; return(RIGHT_OP) ; }
"<<"			{ count() ; return(LEFT_OP) ; }
"++"			{ count() ; return(INC_OP) ; }
"--"			{ count() ; return(DEC_OP) ; }
"->"			{ count() ; return(PTR_OP) ; }
"&&"			{ count() ; return(AND_OP) ; }
"||"			{ count() ; return(OR_OP) ; }
"<="			{ count() ; return(LE_OP) ; }
">="			{ count() ; return(GE_OP) ; }
"=="			{ count() ; return(EQ_OP) ; }
"!="			{ count() ; return(NE_OP) ; }
";"			{ count() ; return(';') ; }
"{"			{ count() ; return('{') ; }
"}"			{ count() ; return('}') ; }
","			{ count() ; return(',') ; }
":"			{ count() ; return(':') ; }
"="			{ count() ; return('=') ; }
"("			{ count() ; return('(') ; }
")"			{ count() ; return(')') ; }
"["			{ count() ; return('[') ; }
"]"			{ count() ; return(']') ; }
"."			{ count() ; return('.') ; }
"&"			{ count() ; return('&') ; }
"!"			{ count() ; return('!') ; }
"~"			{ count() ; return('~') ; }
"-"			{ count() ; return('-') ; }
"+"			{ count() ; return('+') ; }
"*"			{ count() ; return('*') ; }
"/"			{ count() ; return('/') ; }
"%"			{ count() ; return('%') ; }
"<"			{ count() ; return('<') ; }
">"			{ count() ; return('>') ; }
"^"			{ count() ; return('^') ; }
"|"			{ count() ; return('|') ; }
"?"			{ count() ; return('?') ; }

[ \t\v\n\f]		{ count() ; }
.			{ /* ignore bad characters */ }

%%
/*
yywrap()
{
  return(1) ;
}
*/

comment()
{
  char c, c1 ;

loop:
  while ((c=input())!='*' && c!=0)
    putchar(c) ;

  if ((c1=input())!='/' && c!=0)
  {
    unput(c1) ;
    goto loop ;
  }

  if (c != 0)
    putchar(c1) ;
}

comment2()
{
  char c ;
  while ((c=input())!='\n' && c!=0)
    ;
}


int column = 0 ;

void count()
{
  int i ;

  for (i = 0 ; yytext[i] != '\0' ; i++)
    if (yytext[i] == '\n')
    {
      column = 0 ;
      ++yypos  ;
    }
    else if (yytext[i] == '\t')
      column += 8-(column%8) ;
    else
      column++ ;

//  ECHO;
}


int check_type()
{
/*
*  pseudo code --- this is what it should check
*
*  if (yytext == type_name)
*    return(TYPE_NAME) ;
*
*  return(IDENTIFIER) ;
*/

/*
*  it actually will only return IDENTIFIER
*/

  return(IDENTIFIER) ;
}
