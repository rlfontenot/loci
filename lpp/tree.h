#ifndef TREE_HEAD
#define TREE_HEAD
#include <stdio.h>
#include <stdlib.h>
#include "hash.h"

typedef enum OP_ID { 
  OP_UADD, OP_USUB,OP_PLUS_PLUS, OP_MINUS_MINUS, 
  OP_ADD, OP_SUB, OP_MUL, OP_DIV, OP_MOD, OP_EQU, OP_PLUS_EQU,
  OP_MINUS_EQU, OP_MUL_EQU, OP_DIV_EQU, OP_MOD_EQU,
  OP_NOT, 
  OP_AND, OP_OR, OP_XOR, OP_LOGICAL_AND, OP_LOGICAL_OR,
  OP_AND_EQU, OP_OR_EQU, OP_XOR_EQU, OP_EQU_EQU,
  OP_GREATER_EQU, OP_GREATER, OP_LESS, OP_LESS_EQU, OP_NOT_EQU,
  OP_LEFT_SHIFT, OP_RIGHT_SHIFT, OP_LEFT_SHIFT_EQU,
  OP_RIGHT_SHIFT_EQU,
  OP_TILDE, OP_ARROW, OP_ARROW_STAR, OP_COLON, OP_COLON_COLON, 
  OP_TERMINAL, OP_EARMARK, OP_DOT, OP_DOT_STAR, OP_SEMI_COLON,
  OP_COMMA, OP_PARENTHESIS_OPEN, OP_PARENTHESIS_CLOSE,
  OP_BRACKET_OPEN, OP_BRACKET_CLOSE, OP_BRACE_OPEN,
  OP_BRACE_CLOSE, OP_PUNCH, OP_PARENTHESIS, OP_BRACKET,
  OP_BRACE, OP_LESS_GREATER, OP_ARRAY, OP_ARRAY_INDEX,
  OP_VAR_DEF, OP_CONTROL_STRUCTURE, OP_SWITCH_STATEMENT_LIST,
  OP_LABLE, OP_STATEMENT_LIST, OP_TYPEDEF, OP_TD_C, OP_TD_CPP,
  OP_ARGUMENT_DEF, OP_ELIPSIS, OP_RANGE, OP_CASE, OP_DEFAULT,
  OP_IF, OP_ELSE, OP_SWITCH, OP_WHILE, OP_DO, OP_FOR, OP_GOTO,
  OP_CONTINUE, OP_BREAK, OP_RETURN, OP_ID, OP_NUM, OP_STRING,
  OP_FUNC,

  OP_MAP, OP_PARAMETRIC_ID, OP_PARAMETRIC_ARG, /*OP_ITERWARG_ID,*/ OP_ITER_ARG, OP_DOLLAR, OP_JOIN
} op_id ;

typedef enum DP_ID {
  DP_TYPE_NAME, DP_TYPEDEF, DP_EXTERN, DP_STATIC, DP_AUTO, 
  DP_REGISTER, DP_CHAR, DP_SHORT, DP_INT, DP_LONG, DP_SIGNED, 
  DP_UNSIGNED, DP_FLOAT, DP_DOUBLE, DP_CONST, DP_VOLATILE, 
  DP_VOID, DP_STRUCT, DP_UNION, DP_ENUM, DP_BOOL,

  DP_PARAM, DP_STORE, DP_STOREVEC, DP_STOREMAT, DP_MULTISTORE,
  DP_DSTORE, DP_CONSTRAINT, DP_MAP, DP_MAPVEC, DP_MULTIMAP,
  DP_DMAP, DP_DMAPVEC, DP_DMULTIMAP, DP_CONDITIONAL,

  DP_USER
} dp_id ;

typedef enum RULE_ID {
  POINTWISE_T, SINGLETON_T, UNIT_T, APPLY_T
} rule_id ;

typedef enum BOOL {
  false, true
} bool ;

typedef union ATTR {
  struct Var {
    bool is_loci_var ;
    bool is_loci_par_var ;
    /*04/07*/
    bool is_loci_iter_var ;
    /*the var used around
      the map op ->*/
    bool is_map_governed ;
    bool is_output ;
  } var ;
  struct Map {
     /*the -> used in loci_map or
      the common -> used for pointer*/
     bool is_loci_map ;
  } map ;
  /* others to be add later */
} Attr ;

typedef struct EXPTREENODE {
  op_id op ;
  char* name ;
  struct EXPTREENODE* parent ;
  struct EXPTREELIST* list1 ;
  Attr* analysis ;
} ExpTreeNode ;

typedef struct EXPTREELIST {
  ExpTreeNode* node ;
  struct EXPTREELIST* next ;
} ExpTreeList ;

typedef struct FUNCLISTNODE {
  ExpTreeNode* exp ;
  struct FUNCLISTNODE* prev ;
  struct FUNCLISTNODE* next ;
} FuncListNode ;

typedef struct FUNCLISTHEAD {
  /* some info area*/
  /* to be defined later  */
  bool is_rule ;
  bool is_cond ;/*04/24*/
  rule_id rule_type ;
  ExpTreeList* unit_arg ;
  struct TYPE_LIST* apply_arg ;
  ExpTreeList* cond_arg ;/*04/24*/
/*  char* singleton_arg ;*/
  FuncListNode* list ;
} FuncListHead ;

typedef struct FILELISTNODE {
  FuncListHead* func ;
  struct FILELISTNODE* prev ;
  struct FILELISTNODE* next ;
} FileListNode ;

typedef struct FILELISTHEAD {
  /* some info area*/
  /* to be defined later  */
  struct hashtable* def ;
  FileListNode* list ;
} FileListHead ;

/****************************/
typedef struct ID_LIST {
  char* name ;
  char* init ;
  struct ID_LIST* next ;
} id_list ;

typedef struct TYPE_INFO {
  dp_id type ;
  char* user_tname ;
  bool is_elabt ;
  bool begin_temp ;
  bool close_temp ;
} type_info ;

typedef struct TYPE_LIST {
  struct TYPE_INFO* info ;
  struct TYPE_LIST* next ;
} type_list ;

typedef struct TABLE_ENTRY {
  type_list* tlist ;
  char* init ;
} table_entry ;



#endif
