#ifndef PROC_HEAD
#define PROC_HEAD

#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include "tree.h"
#include "hash.h"

void createExpNode(ExpTreeNode **, op_id, char *) ;
void createExpNode_UnaryOp(ExpTreeNode **, op_id, ExpTreeNode *) ;
void createExpNode_BinaryOp(ExpTreeNode **, op_id, ExpTreeNode *, ExpTreeNode *) ;
void create_tlist(type_list **, dp_id, char*, bool, bool, bool, type_list *) ;
void filltable(struct hashtable *, type_list *, id_list *) ;
/*--------------------------------------------------------------------+
| itoa() - manage the sign, compute the string equivalent, and call   |
|          memcpy().                                                  |
+--------------------------------------------------------------------*/
#define INTSIZE 10
char *itoa(int) ;

static int name_count = 1 ;

void Input_Map(ExpTreeList *) ;

void collapse_var_list(ExpTreeList *) ;

void Set_Par_Arg(ExpTreeList *) ;

void Set_Iter_Arg(ExpTreeList *) ;

dp_id Set_dp_id(char *) ;

void TransExp2(ExpTreeNode *, rule_id) ;

void Print_OP_FUNC(ExpTreeNode *, rule_id) ;

void Print_Par_Var(ExpTreeNode *, rule_id) ;

void Print_map_var(ExpTreeNode *) ;

void Print_Iter_Arg(ExpTreeList *) ;

void Print_Type(type_info *) ;

void Print_tlist(type_list *) ;

void PrintOP(op_id) ;

void TransExp(ExpTreeNode *, rule_id) ;

bool Is_Loci_var(char *, struct hashtable *) ;

bool Is_Loci_map_var(char *, struct hashtable *) ;

void ExpNodeAnalysis(ExpTreeNode *, ExpTreeList **, ExpTreeList **) ;

void process_pw_exp(ExpTreeNode *, ExpTreeList **, ExpTreeList **) ;

void create_classname(ExpTreeList *, FuncListHead *) ;

void Print_classname(ExpTreeList *, FuncListHead *) ;

void Print_namestore(ExpTreeList *, rule_id) ;

void Print_input(ExpTreeList *, ExpTreeList *, rule_id) ;

void Print_output(ExpTreeList *, rule_id) ;

void process_pointwise_rule(FuncListHead *) ;

void process_unit_rule(FuncListHead *) ;

void process_apply_rule(FuncListHead *) ;

void process_singleton_rule(FuncListHead *) ;

void process_function(FuncListHead *) ;

void process_loci_rule(FuncListHead *) ;

void process_file_node(FuncListHead *) ;

void process_file(FileListHead *) ;



#endif
