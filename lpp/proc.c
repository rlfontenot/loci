#include "proc.h"
extern struct hashtable* global_def ;


void warning(char* s)
{
   printf("%s\n", s) ;
}

void createExpNode
(ExpTreeNode** enode,
 op_id op, char* s)
{
   ExpTreeNode* tmp ;
 
   tmp = (ExpTreeNode *)malloc
     (sizeof(ExpTreeNode)) ;
   tmp->op = op ;
   tmp->name = s ;
   tmp->parent = NULL ;
   tmp->list1 = NULL ;
   tmp->analysis = NULL ;

   *enode = tmp ;
}


void createExpNode_UnaryOp
(ExpTreeNode** enode,
 op_id op, 
 ExpTreeNode* subnode)
{
   ExpTreeNode* ret ;
 
   ExpTreeList* tmp ;
   tmp = (ExpTreeList *)malloc
     (sizeof(ExpTreeList)) ;
   tmp->node = subnode ;
   tmp->next = NULL ;
   ret = (ExpTreeNode *)malloc
     (sizeof(ExpTreeNode)) ;
   ret->op = op ;
   ret->name = NULL ;
   ret->parent = NULL ;
   ret->list1 = tmp ;
   ret->analysis = NULL ;
   subnode->parent = ret ;

   *enode = ret ;
}


void createExpNode_BinaryOp
(ExpTreeNode** enode,
 op_id op, 
 ExpTreeNode* subnode1,
 ExpTreeNode* subnode2)
{
   ExpTreeNode* ret ;

   ExpTreeList* tmp1 ;
   ExpTreeList* tmp2 ;
   tmp1 = (ExpTreeList *)malloc
     (sizeof(ExpTreeList)) ;
   tmp2 = (ExpTreeList *)malloc
     (sizeof(ExpTreeList)) ;
   tmp1->node = subnode1 ;
   tmp1->next = tmp2 ;
   tmp2->node = subnode2 ;
   tmp2->next = NULL ;
   ret = (ExpTreeNode *)malloc
     (sizeof(ExpTreeNode)) ;
   ret->op = op ;
   ret->name = NULL ;
   ret->parent = NULL ;
   ret->list1 = tmp1 ;
   ret->analysis = NULL ;
   subnode1->parent = ret ;
   subnode2->parent = ret ;

   *enode = ret ;
}


void create_tlist
(type_list** tlist,
 dp_id type,
 char* utname,
 bool e,
 bool b,
 bool c,
 type_list* next)
{
   type_info* tmp ;
   tmp = (type_info *)malloc
     (sizeof(type_info)) ;
   tmp->type = type ;
   tmp->user_tname = utname ;
   tmp->is_elabt = e ;
   tmp->begin_temp = b ;
   tmp->close_temp = c ;
   (*tlist) = (type_list *)
     malloc(sizeof(type_list)) ;
   (*tlist)->info = tmp ;
   (*tlist)->next = next ;
}

void filltable
(struct hashtable* table,
 type_list* tlist,
 id_list* idlist)
{
   table_entry* entry ;
   id_list* curr ;
  
   entry = (table_entry *)
     malloc(sizeof(table_entry)) ;
   entry->tlist = tlist ;
   entry->init = idlist->init ;
  
   curr = idlist ;
   while(curr != NULL)
   {
      if((findhashent(table,curr->name))
         == NULL)
        addptrtohash(table,
                     curr->name,
                     entry) ;
      else
      {
         printf("WARNING: %s redefined globally!\n",
                 curr->name) ;
         printf("Exit PreProcessing!\n") ;
         exit(-1) ;
      }
      curr=curr->next ;
   }

}


/***************Processing*********************/

/*-----------------------------------------------+                     
| itoa() - manage the sign,                      |
|          compute the string equivalent,        |
|          and call memcpy().                    |
+-----------------------------------------------*/
char *itoa(int value)
{
   int count,  /* number of characters in string         */
       i,      /* loop control variable                  */
       sign ;  /* determine if the value is negative     */
   char *ptr,    /* temporary pointer, index into string */
        *string, /* return value                         */
        *temp ;  /* temporary string array               */

   count = 0 ;
   if ((sign=value) < 0)
   /* assign value to sign, if negative,                 */
   /* keep track and invert value                        */
   {                          
      value = -value ;
      count++ ;        /* increment count                */
   }

   /* allocate INTSIZE plus 2 bytes (sign and NULL)      */
   temp = (char *) malloc(INTSIZE+2) ;
   if (temp == NULL)
   {
      return(NULL) ;
   }
   memset(temp, '\0', INTSIZE+2) ;

   string = (char *)malloc(INTSIZE+2) ;
   if (string == NULL)
   {
      return(NULL) ;
   }
   memset(string, '\0', INTSIZE+2) ;
   ptr = string ; /* set temporary ptr to string         */

/*-----------------------------------------------+
| NOTE: This process reverses the order of an    |
|       integer, ie: value = -1234 equates to:   |
|       char [4321-]                             |
|       Reorder the values using for loop below  |
+-----------------------------------------------*/
   do {
     *temp++ = value%10+'0' ;    /* obtain modulus and or with '0'   */
     count++ ;                   /* increment count, track iterations*/
   }  while ((value/=10) >0) ;

   if (sign < 0)              /* add '-' when sign is negative        */
     *temp++ = '-' ;

   *temp-- = '\0' ;           /* ensure null terminated and point     */
                             /* to last char in array                */

/*-----------------------------------------------+
| reorder the resulting char *string:            |
| temp - points to the last char                 |
|        in the temporary array                  |
| ptr  - points to the first element             |
|        in the string array                     |
+-----------------------------------------------*/
   for (i = 0; i < count; i++, temp--, ptr++)
   {
      memcpy(ptr, temp, sizeof(char)) ;
   }

   return(string) ;
}


void process_file
(FileListHead* file)
{
   FileListNode* func ;
 
   func = file->list ;
   while(func != NULL)
   {
      process_file_node(func->func) ;
      func = func->next ;
   }
}

void process_file_node
(FuncListHead* func)
{
   if((func->is_rule) == true)
     process_loci_rule(func) ;
   else
     process_function(func) ;
}

void process_loci_rule
(FuncListHead* func)
{
   switch(func->rule_type)
   {
      case POINTWISE_T:
        process_pointwise_rule(func) ;
        break ;
      case UNIT_T:
        process_unit_rule(func) ;
        break ;
      case APPLY_T:
        process_apply_rule(func) ;
        break ;
      case SINGLETON_T:
        process_singleton_rule(func) ;
        break ;
      default:
        break ;
   }
}

void process_function
(FuncListHead* func)
{
}


void process_singleton_rule
(FuncListHead* func)
{
   rule_id rtype ;
   FuncListNode* node ;
   ExpTreeList* loci_var ;
   ExpTreeList* loci_curr_var ;
   ExpTreeList* loci_map ;
   ExpTreeList* loci_curr_map ;
   table_entry* entry ;
  
   rtype = func->rule_type ;

   loci_var = (ExpTreeList *)
     malloc(sizeof(ExpTreeList)) ;
   loci_var->node = NULL ;
   loci_var->next = NULL ;
   loci_curr_var = loci_var ;

   loci_map = (ExpTreeList *)
     malloc(sizeof(ExpTreeList)) ;
   loci_map->node = NULL ;
   loci_map->next = NULL ;
   loci_curr_map = loci_map ;

   node = func->list ;
   while(node != NULL)
   {
      process_pw_exp(node->exp,
                     &loci_curr_var,
                     &loci_curr_map) ;
      node = node->next ;
   }
   /* find out the same repeated variables in the list */
   collapse_var_list(loci_var) ;

   /*
    *****************************************
    ******* write out the rule class ********
   *
   printf("\n") ;
   printf("class %s : public singleton_rule {\n", classname) ;
   */
   Print_classname(loci_var, func) ;
   Print_namestore(loci_var, func->rule_type) ;
   Print_input(loci_var, loci_map, func->rule_type) ;
   Print_output(loci_var, func->rule_type) ;

   if(func->is_cond == true)
   {
      entry = (table_entry *)findhashent
        (global_def, func->cond_arg->node->name) ;
      if((entry->tlist->info->type==DP_PARAM) &&
         (entry->tlist->next->info->type==DP_BOOL))
      {
         printf("\t") ;
         printf("conditional(\"%s{n}\") ;\n", 
                func->cond_arg->node->name) ;

         printf("}\n") ;
      }
      else printf("This is not a conditional rule!\n") ;
   }

   printf("\t") ; 
   printf("virtual void compute(const sequence &seq) {\n") ;
   node = func->list ;
   while(node != NULL)
   {
      printf("\t\t") ;
      TransExp(node->exp, rtype) ;
      printf(" ;\n") ;
      node = node->next ;
   }
   printf("\t") ; printf("}\n") ;


   printf("} ;\n") ;


   printf("\n") ;
   printf("register_rule<") ;
   create_classname(loci_var, func) ;
   printf("> register_") ;
   create_classname(loci_var, func) ;
   printf(" ;\n") ;
   printf("\n") ;
}

void process_unit_rule
(FuncListHead* func)
{
   rule_id rtype ;
   FuncListNode* node ;
   ExpTreeList* loci_var ;
   ExpTreeList* loci_curr_var ;
   ExpTreeList* loci_map ;
   ExpTreeList* loci_curr_map ;
   table_entry* entry ;
  
   rtype = func->rule_type ;

   loci_var = (ExpTreeList *)
     malloc(sizeof(ExpTreeList)) ;
   loci_var->node = NULL ;
   loci_var->next = NULL ;
   loci_curr_var = loci_var ;

   loci_map = (ExpTreeList *)
     malloc(sizeof(ExpTreeList)) ;
   loci_map->node = NULL ;
   loci_map->next = NULL ;
   loci_curr_map = loci_map ;


   node = func->list ;
   while(node != NULL)
   {
      process_pw_exp(node->exp,
                     &loci_curr_var,
                     &loci_curr_map) ;
      node = node->next ;
   }
   /* find out the same repeated variables in the list */
   collapse_var_list(loci_var) ;

   /*
    *****************************************
    ******* write out the rule class ********
   *
   printf("\n") ;
   printf("class %s : public unit_rule {\n", classname) ;
*/
   Print_classname(loci_var, func) ;
   Print_namestore(loci_var, func->rule_type) ;

   /*   Print_input(loci_var, loci_map, func->rule_type) ;*/
   printf("\t\t") ; printf("constraint(\"%s\") ;\n",
                           func->unit_arg->node->name) ;

   Print_output(loci_var, func->rule_type) ;

   printf("\t") ; printf("void calculate(Entity e) {\n") ;
   node = func->list ;
   while(node != NULL)
   {
      printf("\t\t") ;
      TransExp(node->exp, rtype) ;
      printf(" ;\n") ;
      node = node->next ;
   }
   printf("\t") ; printf("}\n") ;


   printf("\t") ; 
   printf("virtual void compute(const sequence &seq) {\n") ;
   printf("\t\t") ;
   printf("do_loop(seq,this,&") ;
   create_classname(loci_var, func) ;
   printf("::calculate) ;\n") ;
   printf("\t") ; printf("}\n") ;


   printf("} ;\n") ;


   printf("\n") ;
   printf("register_rule<") ;
   create_classname(loci_var, func) ;
   printf("> register_") ;
   create_classname(loci_var, func) ;
   printf(" ;\n") ;
   printf("\n") ;
}

void process_apply_rule
(FuncListHead* func)
{
   rule_id rtype ;
   FuncListNode* node ;
   ExpTreeList* loci_var ;
   ExpTreeList* loci_curr_var ;
   ExpTreeList* loci_map ;
   ExpTreeList* loci_curr_map ;
   table_entry* entry ;
  
   rtype = func->rule_type ;

   loci_var = (ExpTreeList *)
     malloc(sizeof(ExpTreeList)) ;
   loci_var->node = NULL ;
   loci_var->next = NULL ;
   loci_curr_var = loci_var ;

   loci_map = (ExpTreeList *)
     malloc(sizeof(ExpTreeList)) ;
   loci_map->node = NULL ;
   loci_map->next = NULL ;
   loci_curr_map = loci_map ;

   node = func->list ;
   while(node != NULL)
   {
      process_pw_exp(node->exp,
                     &loci_curr_var,
                     &loci_curr_map) ;
      node = node->next ;
   }
   /* find out the same repeated variables in the list */
   collapse_var_list(loci_var) ;

   /*
    *****************************************
    ******* write out the rule class ********
   *
   printf("\n") ;
   printf("class %s : public apply_rule<", classname) ;
   */
   Print_classname(loci_var, func) ;
   Print_namestore(loci_var, func->rule_type) ;
   Print_input(loci_var, loci_map, func->rule_type) ;

   loci_curr_var = loci_var ;
   printf("\t\t") ; printf("output(\"") ;
   while(loci_curr_var->next != NULL)
   {
      if(loci_curr_var->node->analysis->var.is_output
         == true)
      {
         if(loci_curr_var->node->analysis->var.is_map_governed
            == true)
           Print_map_var(loci_curr_var->node->parent) ;

         if(loci_curr_var->node->analysis->var.is_loci_par_var
            == true)
         {
            Print_OP_FUNC(loci_curr_var->node, rtype) ;
            printf("\") ;") ; 
         }
         /*04/07*/
         else
           if(loci_curr_var->node->analysis->var.is_loci_iter_var
              == true)
             printf("$%s\") ;", loci_curr_var->node->name) ;
           else
             printf("%s\") ;", loci_curr_var->node->name) ;
         break ;
      }
      loci_curr_var = loci_curr_var->next ;
   }
   printf("\n") ;

   loci_curr_var = loci_var ;
   printf("\t\t") ; printf("constraint(\"") ;
   while(loci_curr_var->next != NULL)
   {
      if(loci_curr_var->node->analysis->var.is_output
         == true)
      {
         if(loci_curr_var->node->analysis->var.is_map_governed
            == true)
            Print_map_var(loci_curr_var->node->parent) ;

         if(loci_curr_var->node->analysis->var.is_loci_par_var
            == true)
         {
            Print_OP_FUNC(loci_curr_var->node, rtype) ;
            printf("\") ;") ; 
         }
         else
            printf("%s\") ;", loci_curr_var->node->name) ;
         break ;
      }
      loci_curr_var = loci_curr_var->next ;
   }
   printf("\n") ;
   printf("\t") ; printf("}\n") ;


   printf("\t") ; printf("void calculate(Entity e) {\n") ;
   node = func->list ;
   while(node != NULL)
   {
      printf("\t\t") ;
      TransExp(node->exp, rtype) ;
      printf(" ;\n") ;
      node = node->next ;
   }
   printf("\t") ; printf("}\n") ;


   printf("\t") ; 
   printf("virtual void compute(const sequence &seq) {\n") ;
   printf("\t\t") ;
   printf("do_loop(seq,this,&") ;
   create_classname(loci_var, func) ;
   printf("::calculate) ;\n") ;
   printf("\t") ; printf("}\n") ;


   printf("} ;\n") ;


   printf("\n") ;
   printf("register_rule<") ;
   create_classname(loci_var, func) ;
   printf("> register_") ;
   create_classname(loci_var, func) ;
   printf(" ;\n") ;
   printf("\n") ;
}

void process_pointwise_rule
(FuncListHead* func)
{
   rule_id rtype ;
   FuncListNode* node ;
   ExpTreeList* loci_var ;
   ExpTreeList* loci_curr_var ;
   ExpTreeList* loci_map ;
   ExpTreeList* loci_curr_map ;
   table_entry* entry ;
  
   rtype = func->rule_type ;

   loci_var = (ExpTreeList *)
     malloc(sizeof(ExpTreeList)) ;
   loci_var->node = NULL ;
   loci_var->next = NULL ;
   loci_curr_var = loci_var ;

   loci_map = (ExpTreeList *)
     malloc(sizeof(ExpTreeList)) ;
   loci_map->node = NULL ;
   loci_map->next = NULL ;
   loci_curr_map = loci_map ;

   node = func->list ;
   while(node != NULL)
   {
      process_pw_exp(node->exp,
                     &loci_curr_var,
                     &loci_curr_map) ;
      node = node->next ;
   }
   /* find out the same repeated variables in the list */
   collapse_var_list(loci_var) ;

   /*
    *****************************************
    ******* write out the rule class ********
   *
   printf("\n") ;
   printf("class %s : public pointwise_rule {\n", classname) ;
*/
   Print_classname(loci_var, func) ;
   Print_namestore(loci_var, func->rule_type) ;
   Print_input(loci_var, loci_map, func->rule_type) ;
   Print_output(loci_var, func->rule_type) ;

   printf("\t") ; printf("void calculate(Entity e) {\n") ;
   node = func->list ;
   while(node != NULL)
   {
      printf("\t\t") ;
      TransExp(node->exp, rtype) ;
      printf(" ;\n") ;
      node = node->next ;
   }
   printf("\t") ; printf("}\n") ;


   printf("\t") ; 
   printf("virtual void compute(const sequence &seq) {\n") ;
   printf("\t\t") ;
   printf("do_loop(seq,this,&") ;
   create_classname(loci_var, func) ;
   printf("::calculate) ;\n") ;
   printf("\t") ; printf("}\n") ;


   printf("} ;\n") ;


   printf("\n") ;
   printf("register_rule<") ;
   create_classname(loci_var, func) ;
   printf("> register_") ;
   create_classname(loci_var, func) ;
   printf(" ;\n") ;
   printf("\n") ;
}

void process_pw_exp
(ExpTreeNode* exp,
 ExpTreeList** var_list,
 ExpTreeList** map_list)
{
   ExpNodeAnalysis(exp, var_list, map_list) ;
}

void create_classname
(ExpTreeList* loci_var,
 FuncListHead* func)
{
   rule_id rtype ;
   ExpTreeList* loci_curr_var ;

   rtype = func->rule_type ;
   loci_curr_var = loci_var ;

   while(loci_curr_var->next != NULL)
   {
      if(loci_curr_var->node->analysis->var.is_loci_par_var
         == true)
         Print_Par_Var(loci_curr_var->node, rtype) ;
      else
      {
         if(loci_curr_var->node->analysis->var.is_loci_iter_var
            == true)
           printf("$%s", loci_curr_var->node->name) ;
         else
           printf("%s", loci_curr_var->node->name) ;
      }
      if(loci_curr_var->next->next != NULL)
         printf("_") ; 
      loci_curr_var = loci_curr_var->next ;
   }
}

void Print_classname
(ExpTreeList* loci_var,
 FuncListHead* func)
{
   rule_id rtype ;
   ExpTreeList* loci_curr_var ;
   table_entry* entry ;
   rtype = func->rule_type ;

   printf("\n") ;
   printf("class ") ;
   create_classname(loci_var, func) ;

   switch(rtype)
   {
      case POINTWISE_T:
        printf(" : public pointwise_rule {\n") ;
        break ;
      case UNIT_T:
        printf(" : public unit_rule {\n") ;
        break ;
      case APPLY_T:
        printf(" : public apply_rule <") ;/*add apply_arg later*/
        loci_curr_var = loci_var ;
        while(loci_curr_var->next != NULL)
        {
           if(loci_curr_var->node->analysis->var.is_output
              == true)
           {
              entry = (table_entry *)findhashent
                (global_def, loci_curr_var->node->name) ;
              Print_tlist(entry->tlist) ;
              break ;
           }
           loci_curr_var = loci_curr_var->next ;
        }
        if(func->apply_arg != NULL)
	{
           printf(", ") ;
           Print_tlist(func->apply_arg) ;
        }
        printf(" > {\n") ;
        break ;
      case SINGLETON_T:
        printf(" : public singleton_rule {\n") ;
        break ;
      default:
        printf(" : public {\n") ;
        break ;
   }

   loci_curr_var = loci_var ;
   while(loci_curr_var->next != NULL)
   {
      printf("\t") ;
      if(loci_curr_var->node->analysis->var.is_output
         == false)
        printf("const_") ;
      entry = (table_entry *)findhashent
        (global_def, loci_curr_var->node->name) ;
      Print_tlist(entry->tlist) ; printf(" ") ;
      if(loci_curr_var->node->analysis->var.is_loci_par_var
         == true)
      {
         Print_Par_Var(loci_curr_var->node, rtype) ;
         printf(" ;\n") ;
      }
      else
         printf("%s ;\n", loci_curr_var->node->name) ;
      loci_curr_var = loci_curr_var->next ;
   }
   printf("public:\n") ;
   printf("\t") ;
   create_classname(loci_var, func) ;
   printf("() {\n") ;
}

void Print_namestore
(ExpTreeList* loci_var,
 rule_id rtype)
{
   ExpTreeList* loci_curr_var ;

   loci_curr_var = loci_var ;
   while(loci_curr_var->next != NULL)
   {
      printf("\t\t") ; printf("name_store") ;
      if(loci_curr_var->node->analysis->var.is_loci_par_var
         == true)
      {
         printf("(\"") ; Print_OP_FUNC(loci_curr_var->node, rtype) ;
         printf("\",") ; Print_Par_Var(loci_curr_var->node, rtype) ;
         printf(") ;\n") ;
      }
      /*04/07*/
      else
      {
         switch(rtype)
         {
            case POINTWISE_T:
            case UNIT_T:
            case APPLY_T:
              if(loci_curr_var->node->analysis->var.is_loci_iter_var
                 == true)
                printf("(\"$%s\",%s) ;\n", loci_curr_var->node->name,
                                           loci_curr_var->node->name) ;
              else
                printf("(\"%s\",%s) ;\n", loci_curr_var->node->name,
                                          loci_curr_var->node->name) ;
              break ;
            /*04/08*/
            case SINGLETON_T:
              if(loci_curr_var->node->analysis->var.is_loci_iter_var
                 == true)
                printf("(\"$") ;
              else
                printf("(\"") ;

              printf("%s", loci_curr_var->node->name,
                           loci_curr_var->node->name) ;
              if((loci_curr_var->node->list1!=NULL) &&
                 (loci_curr_var->node->list1->node->op==OP_ITER_ARG))
	      {
                 printf("{") ;
                 Print_Iter_Arg(loci_curr_var->node->list1) ;
                 printf("}") ;
              }
              printf("\",%s) ;\n", loci_curr_var->node->name,
                                   loci_curr_var->node->name) ;
              break ;
            default:
              break ;
         }
      }
      loci_curr_var = loci_curr_var->next ;
   }
}

void Print_input
(ExpTreeList* loci_var,
 ExpTreeList* loci_map,
 rule_id rtype)
{
   ExpTreeList* loci_curr_var ;
   ExpTreeList* loci_curr_map ;

   loci_curr_var = loci_var ;
   loci_curr_map = loci_map ;
   printf("\t\t") ; printf("input(\"") ;
   while(loci_curr_var->next != NULL)
   {
      if(loci_curr_var->node->analysis->var.is_map_governed
         == false)
      {
         if(loci_curr_var->node->analysis->var.is_output
            == false)
         {
	    if(loci_curr_var->node->analysis->var.is_loci_par_var
               == true)
            {
	       Print_OP_FUNC(loci_curr_var->node, rtype) ;
               if(loci_curr_var->next->next != NULL)
                 printf(",") ;
	       /*	       printf(",") ;*/
            }
            /*04/07*/
            else
            {
               switch(rtype)
               {
                  case POINTWISE_T:
		    /*case UNIT_T:*/
                  case APPLY_T:
                    if(loci_curr_var->node->analysis->var.is_loci_iter_var
                       == true)
                      printf("$%s", loci_curr_var->node->name) ;
		    /*                      printf("$%s,", loci_curr_var->node->name) ;*/
                    else
                      printf("%s", loci_curr_var->node->name) ;
		    /*                      printf("%s,", loci_curr_var->node->name) ;*/
                    break ;
                  /*04/08*/
                  case SINGLETON_T:
                    if(loci_curr_var->node->analysis->var.is_loci_iter_var == true)
                      printf("$") ;

                    printf("%s", loci_curr_var->node->name,
                                 loci_curr_var->node->name) ;
                    if((loci_curr_var->node->list1!=NULL) &&
                       (loci_curr_var->node->list1->node->op==OP_ITER_ARG))
                    {
                       printf("{") ;
                       Print_Iter_Arg(loci_curr_var->node->list1) ;
                       printf("}") ;
                    }
		    /*                    printf(",") ;*/
                    break ;
                  default:
                    break ;
               }
               if(loci_curr_var->next->next != NULL)
                 printf(",") ;
            }
         }
      }
      loci_curr_var = loci_curr_var->next ;
   }
   Input_Map(loci_map) ;
   printf("\") ;") ; printf("\n") ;
}

void Print_output
(ExpTreeList* loci_var,
 rule_id rtype)
{
   ExpTreeList* loci_curr_var ;

   loci_curr_var = loci_var ;
   printf("\t\t") ; printf("output(\"") ;
   while(loci_curr_var->next != NULL)
   {
      if(loci_curr_var->node->analysis->var.is_output
         == true)
      {
         if(loci_curr_var->node->analysis->var.is_loci_par_var
            == true)
         {
            Print_OP_FUNC(loci_curr_var->node, rtype) ;
	    printf("\") ;") ; 
         }
         /*04/07*/
         else
         {
            switch(rtype)
            {
               case POINTWISE_T:
               case UNIT_T:
                 if(loci_curr_var->node->analysis->var.is_loci_iter_var
                    == true)
                   printf("$%s\") ;", loci_curr_var->node->name) ;
                 else
                   printf("%s\") ;", loci_curr_var->node->name) ;
                 break ;
               /*04/08*/
               case SINGLETON_T:
                 if(loci_curr_var->node->analysis->var.is_loci_iter_var == true)
                   printf("$") ;

                 printf("%s", loci_curr_var->node->name,
                              loci_curr_var->node->name) ;
                 if((loci_curr_var->node->list1!=NULL) &&
                    (loci_curr_var->node->list1->node->op==OP_ITER_ARG))
                 {
                    printf("{") ;
                    Print_Iter_Arg(loci_curr_var->node->list1) ;
                    printf("}") ;
                 }
                 printf("\") ;") ;
                 break ;
               default:
                 break ;
            }
         }
         break ;
      }
      loci_curr_var = loci_curr_var->next ;
   }
   printf("\n") ;
   printf("\t") ; printf("}\n") ;
}

/* 
 * find out loci specific variables 
 * and operators in the expression tree,
 * and check the "in/out put" of a variable
 */
void ExpNodeAnalysis
(ExpTreeNode* exp,
 ExpTreeList** var_list,
 ExpTreeList** map_list)
{
   ExpTreeList* curr ;
   Attr* att ;
   ExpTreeList* advance ;
   static bool output = false ;

   exp->analysis = (Attr *)
     malloc(sizeof(Attr)) ;
   att = exp->analysis ;

   if(exp->list1 != NULL)
     curr = exp->list1 ;
   else
     curr = NULL ;

   switch(exp->op)
   {
      case OP_ID:
        if(Is_Loci_var(exp->name,global_def)
           == true)
        {
           att->var.is_loci_var = true ;

           /*04/07*/
           if((exp->parent!=NULL) &&
              (exp->parent->op==OP_DOLLAR))
             att->var.is_loci_iter_var = true ;

           if((output==true) &&
              (Is_Loci_map_var(exp->name,global_def)
               ==false))
           {
              att->var.is_output = true ;
              output = false ;
           }
           else
              att->var.is_output = false ;

           if((exp->parent!=NULL) &&
              (exp->parent->op==OP_MAP))
	     att->var.is_map_governed = true ;
           else 
             att->var.is_map_governed = false ;

           advance = (*var_list) ;
           advance->node = exp ;
           advance->next = (ExpTreeList *)
             malloc(sizeof(ExpTreeList)) ;
           *var_list = advance->next ;
           (*var_list)->node = NULL ;
           (*var_list)->next = NULL ;
        }
        break ;
      case OP_PARAMETRIC_ID:
        att->var.is_loci_var = true ;
        att->var.is_loci_par_var = true ;
        if(output == true)
        {
           att->var.is_output = true ;
           output = false ;
        }
        else
          att->var.is_output = false ;
        if((exp->parent!=NULL) &&
           (exp->parent->op==OP_MAP))
          att->var.is_map_governed = true ;
        else
          att->var.is_map_governed = false ;

        advance = (*var_list) ;
        advance->node = exp ;
        advance->next = (ExpTreeList *)
          malloc(sizeof(ExpTreeList)) ;
        *var_list = advance->next ;
        (*var_list)->node = NULL ;
        (*var_list)->next = NULL ;
        break ;
      case OP_MAP:
        att->map.is_loci_map = true ;

        advance = (*map_list) ;
        advance->node = exp ;
        advance->next = (ExpTreeList *)
          malloc(sizeof(ExpTreeList)) ;
        *map_list = advance->next ;
        (*map_list)->node = NULL ;
        (*map_list)->next = NULL ;
        break ;
      case OP_EQU:
      case OP_PLUS_EQU:
      case OP_MINUS_EQU:
      case OP_MUL_EQU:
      case OP_DIV_EQU:
      case OP_MOD_EQU:
      case OP_AND_EQU:
      case OP_OR_EQU:
      case OP_XOR_EQU:
      case OP_LEFT_SHIFT_EQU:
      case OP_RIGHT_SHIFT_EQU:
      case OP_JOIN:
        output = true ;
        break ;
      default:
        break ;
   }
   /* recursive analyze the exp tree */
   while(curr != NULL)
   {
      ExpNodeAnalysis(curr->node,
                      var_list,
                      map_list) ;
      curr = curr->next ;
   }
}

/* need to be add types later */
bool Is_Loci_var
(char* s, struct hashtable* table)
{
  table_entry* entry ;
  dp_id type ;

  entry = (table_entry *)
    findhashent(table, s) ;
  if(entry == NULL)
    return false ;

  type = entry->tlist->info->type ;
  if((type==DP_PARAM) || (type==DP_STORE) ||
     (type==DP_STOREVEC) || (type==DP_STOREMAT) ||
     (type==DP_MULTISTORE) || (type==DP_DSTORE) ||
     (type==DP_CONSTRAINT) || (type==DP_MAP) ||
     (type==DP_MAPVEC) || (type==DP_MULTIMAP) ||
     (type==DP_DMAP) || (type==DP_DMAPVEC) ||
     (type==DP_DMULTIMAP))
    return true ;
  else
    return false ;
}

bool Is_Loci_map_var
(char* s, struct hashtable* table)
{
  table_entry* entry ;
  dp_id type ;

  entry = (table_entry *)
    findhashent(table, s) ;
  if(entry == NULL)
    return false ;

  type = entry->tlist->info->type ;
  if((type==DP_MAP) || (type==DP_MAPVEC) ||
     (type==DP_MULTIMAP) || (type==DP_DMAP) ||
     (type==DP_DMAPVEC) || (type==DP_DMULTIMAP))
    return true ;
  else
    return false ;
}

void TransExp
(ExpTreeNode* exp,
 rule_id rtype)
{
  table_entry* entry ;
  dp_id type ;

  ExpTreeList* list = exp->list1 ;

  switch(exp->op)
  {
     case OP_UADD:
     case OP_USUB:
     case OP_NOT:
       printf("( ") ;
       PrintOP(exp->op) ;
       TransExp(list->node, rtype) ;
       printf(" )") ;
       break ;
     case OP_ADD:
     case OP_SUB:
     case OP_MUL:
     case OP_DIV:
     case OP_MOD:
     case OP_AND:
     case OP_OR:
     case OP_XOR:
     case OP_LOGICAL_AND: 
     case OP_LOGICAL_OR:
     case OP_LESS:
     case OP_GREATER:
     case OP_EQU_EQU:
     case OP_NOT_EQU:
     case OP_LESS_EQU:
     case OP_GREATER_EQU:
     case OP_LEFT_SHIFT:
     case OP_RIGHT_SHIFT:
       printf("( ") ;
       TransExp(list->node, rtype) ;
       PrintOP(exp->op) ;
       TransExp(list->next->node, rtype) ;
       printf(" )") ;
       break ;
     case OP_EQU:
     case OP_PLUS_EQU:
     case OP_MINUS_EQU:
     case OP_MUL_EQU:
     case OP_DIV_EQU:
     case OP_MOD_EQU:
     case OP_AND_EQU:
     case OP_OR_EQU:
     case OP_XOR_EQU:
     case OP_LEFT_SHIFT_EQU:
     case OP_RIGHT_SHIFT_EQU:
       TransExp(list->node, rtype) ;
       PrintOP(exp->op) ;
       TransExp(list->next->node, rtype) ;
       break ;
     case OP_ID:
       if(Is_Loci_var(exp->name,global_def)
          == true)
       {
         /*04/07*/
         if(rtype == SINGLETON_T)
           printf("*%s", exp->name) ;
         else
           printf("%s[e]", exp->name) ;
       }
       else
         printf("%s", exp->name) ;
       break ;
     case OP_NUM:
       printf("%s", exp->name) ;
       break ;
     case OP_PARAMETRIC_ID:
       Print_Par_Var(exp, rtype) ;
       printf("[e]") ;
       break ;
     case OP_FUNC:
       Print_OP_FUNC(exp, rtype) ;
       break ;    
     case OP_MAP:
       TransExp2(list->next->node, rtype) ;
       printf("[") ;
       TransExp2(list->node, rtype) ;
       printf("[e]]") ;
       break ;
     case OP_JOIN:
       printf("join(") ;
       TransExp(list->node, rtype) ;
       printf(", ") ;
       TransExp(list->next->node, rtype) ;
       printf(")") ;
       break ;
     case OP_DOLLAR:
       entry = (table_entry *)findhashent
         (global_def, exp->list1->node->name) ;
       type = entry->tlist->next->info->type ;

       entry = (table_entry *)findhashent
         (global_def, exp->parent->list1->next->node->name) ;
       if(type == entry->tlist->next->info->type)
         TransExp(list->node, rtype) ;
       else
       {
          Print_Type(entry->tlist->next->info) ;
          printf("( ") ;
          TransExp(list->node, rtype) ;
          printf(" )") ;
       }
       break ;
     default:
       break ;
  }
}

void TransExp2
(ExpTreeNode* exp,
 rule_id rtype)
{
  table_entry* entry ;
  dp_id type ;

  ExpTreeList* list = exp->list1 ;

  switch(exp->op)
  {
    case OP_UADD:
    case OP_USUB:
    case OP_NOT:
      printf("( ") ;
      PrintOP(exp->op) ;
      TransExp(list->node, rtype) ;
      printf(" )") ;
      break ;
    case OP_ADD:
    case OP_SUB:
    case OP_MUL:
    case OP_DIV:
    case OP_MOD:
    case OP_AND:
    case OP_OR:
    case OP_XOR:
    case OP_LOGICAL_AND: 
    case OP_LOGICAL_OR:
    case OP_LESS:
    case OP_GREATER:
    case OP_EQU_EQU:
    case OP_NOT_EQU:
    case OP_LESS_EQU:
    case OP_GREATER_EQU:
    case OP_LEFT_SHIFT:
    case OP_RIGHT_SHIFT:
      printf("( ") ;
      TransExp(list->node, rtype) ;
      PrintOP(exp->op) ;
      TransExp(list->next->node, rtype) ;
      printf(" )") ;
      break ;
    case OP_EQU:
    case OP_PLUS_EQU:
    case OP_MINUS_EQU:
    case OP_MUL_EQU:
    case OP_DIV_EQU:
    case OP_MOD_EQU:
    case OP_AND_EQU:
    case OP_OR_EQU:
    case OP_XOR_EQU:
    case OP_LEFT_SHIFT_EQU:
    case OP_RIGHT_SHIFT_EQU:
      TransExp(list->node, rtype) ;
      PrintOP(exp->op) ;
      TransExp(list->next->node, rtype) ;
      break ;
    case OP_ID:
    case OP_PARAMETRIC_ARG:
      printf("%s", exp->name) ;
      break ;
    case OP_PARAMETRIC_ID:
      Print_OP_FUNC(exp, rtype) ;
      break ;
    case OP_JOIN:
      printf("join(") ;
      TransExp(list->node, rtype) ;
      printf(", ") ;
      TransExp(list->next->node, rtype) ;
      printf(")") ;
      break ;
    case OP_DOLLAR:
      entry = (table_entry *)findhashent
        (global_def, exp->list1->node->name) ;
      type = entry->tlist->next->info->type ;

      entry = (table_entry *)findhashent
        (global_def, exp->parent->list1->next->node->name) ;
      if(type == entry->tlist->next->info->type)
        TransExp(list->node, rtype) ;
      else
      {
         Print_Type(entry->tlist->next->info) ;
         printf("( ") ;
         TransExp(list->node, rtype) ;
         printf(" )") ;
      }
      break ;
    default:
      break ;
  }
}

void PrintOP(op_id op)
{
  switch(op)
  {
     case OP_NOT: printf( " ! " ) ; break ;
     case OP_UADD:
     case OP_ADD: printf( " + " ) ; break ;
     case OP_USUB:
     case OP_SUB: printf( " - " ) ; break ;
     case OP_MUL: printf( " * " ) ; break ;
     case OP_DIV: printf( " / " ) ; break ;
     case OP_MOD: printf( " % " ) ; break ;
     case OP_EQU: printf( " = " ) ; break ;
     case OP_PLUS_EQU:  printf( " += " ) ; break ;
     case OP_MINUS_EQU: printf( " -= " ) ; break ;
     case OP_MUL_EQU:   printf( " *= " ) ; break ;
     case OP_DIV_EQU:   printf( " /= " ) ; break ;
     case OP_MOD_EQU:   printf( " %= " ) ; break ;
     case OP_AND: printf( " & " ) ; break ;
     case OP_OR:  printf( " | " ) ; break ;
     case OP_XOR: printf( " ^ " ) ; break ;
     case OP_LOGICAL_AND: printf( " && " ) ; break ;
     case OP_LOGICAL_OR:  printf( " || " ) ; break ;
     case OP_AND_EQU:  printf( " &= " ) ; break ;
     case OP_OR_EQU:   printf( " |= " ) ; break ;
     case OP_XOR_EQU:  printf( " ^= " ) ; break ;
     case OP_LESS:     printf( " < " ) ; break ;
     case OP_GREATER:  printf( " > " ) ; break ;
     case OP_EQU_EQU:  printf( " == " ) ; break ;
     case OP_NOT_EQU:  printf( " != " ) ; break ;
     case OP_LESS_EQU: printf( " <= " ) ; break ;
     case OP_GREATER_EQU: printf( " >= " ) ; break ;
     case OP_LEFT_SHIFT:  printf( " << " ) ; break ;
     case OP_RIGHT_SHIFT: printf( " >> " ) ; break ;
     case OP_LEFT_SHIFT_EQU:  printf( " <<= " ) ; break ;
     case OP_RIGHT_SHIFT_EQU: printf( " >>= " ) ; break ;
     default: break ;
  }
}

void Print_tlist
(type_list* tlist)
{
  int temp_close = 0 ;
  int i ;

  while(tlist != NULL)
  {
     Print_Type(tlist->info) ;
     if(tlist->info->is_elabt == true)
       printf("::") ;
     else
     {
        if(tlist->info->begin_temp
           == true)
        {
           printf("<") ;
           ++temp_close ;
        }
        else
        {
           if(tlist->info->close_temp
              == true)
           {
              printf(">") ;
              --temp_close ;
           }
           if(tlist->next != NULL)
             printf(", ") ;
        }
     }
     tlist = tlist->next ;
  }

  if(temp_close > 0)
  {
     for(i=0; i<temp_close; ++i)
       printf(" >") ;
  }
}

void Print_Type
(type_info* info)
{
  if(info->type != DP_USER)
  {
     switch(info->type)
     {
        case DP_PARAM:
          printf("param") ; break ;
        case DP_STORE:
          printf("store") ; break ;
        case DP_STOREVEC:
          printf("storeVec") ; break ;
        case DP_STOREMAT:
          printf("storeMat") ; break ;
        case DP_MULTISTORE:
          printf("multiStore") ; break ;
        case DP_DSTORE:
          printf("dstore") ; break ;
        case DP_CONSTRAINT:
          printf("constraint") ; break ;
        case DP_MAP:
          printf("Map") ; break ;
        case DP_MAPVEC:
          printf("MapVec") ; break ;
        case DP_MULTIMAP:
          printf("multiMap") ; break ;
        case DP_DMAP:
          printf("dMap") ; break ;
        case DP_DMAPVEC:
          printf("dMapVec") ; break ;
        case DP_DMULTIMAP:
          printf("dmultiMap") ; break ;
        case DP_VOID:
          printf("void") ; break ;
        case DP_CHAR:
          printf("char") ; break ;
        case DP_SHORT:
          printf("short") ; break ;
        case DP_INT:
          printf("int") ; break ;
        case DP_FLOAT:
          printf("float") ; break ;
        case DP_DOUBLE:
          printf("double") ; break ;
        case DP_LONG:
          printf("long") ; break ;
        case DP_SIGNED:
          printf("signed") ; break ;
        case DP_UNSIGNED:
          printf("unsigned") ; break ;
        case DP_BOOL:
          printf("bool") ; break ;
        default: break ;
     }
  }
  else
  {
     printf("%s", info->user_tname) ;
  }
}

void Print_OP_FUNC
(ExpTreeNode * node,
 rule_id rtype)
{
  ExpTreeList *args = node->list1 ;
  
  printf("%s(", node->name) ;
  while(args->next != NULL)
  {
    TransExp2(args->node, rtype) ;
    printf(",") ;
    args = args->next ;
  }
  TransExp2(args->node, rtype) ;
  printf(")") ;
}

void Print_Par_Var
(ExpTreeNode * node,
 rule_id rtype)
{
  ExpTreeList *args = node->list1 ;

  printf("%s", node->name) ;
  while(args != NULL)
  {
    printf("_") ;
    TransExp2(args->node, rtype) ;
    args = args->next ;
  }
}

void Print_map_var
(ExpTreeNode *node)
{
  ExpTreeNode* subnode1 ;
  ExpTreeNode* subnode2 ;

  if(node->list1 == NULL)
  {
    subnode1 = NULL ;
    subnode2 = NULL ;
  }
  else
  {
    subnode1 = node->list1->node ;
    subnode2 = node->list1->next->node ;
  }
  
  if(subnode1 != NULL)
    Print_map_var(subnode1) ;

    if((node->op == OP_ID) &&
       (node->analysis->var.is_output
        == false))
      printf("%s", node->name) ;
    if(node->op == OP_MAP)
      printf("->") ;

  if(subnode2 != NULL)
    Print_map_var(subnode2) ;
}


void Print_Iter_Arg
(ExpTreeList *arg_list)
{
  ExpTreeList *curr = arg_list->next ;

  printf("%s", arg_list->node->name) ;
  if(arg_list->next != NULL)
    printf(",") ;
  while(curr != NULL)
  {
     Print_Iter_Arg(curr) ;
     curr = curr->next;
  }
}


void Set_Par_Arg
(ExpTreeList *arg_list)
{
  ExpTreeList *curr = arg_list->node->list1 ;

  if(arg_list->node->op == OP_ID)
    arg_list->node->op = OP_PARAMETRIC_ARG ;

  while(curr != NULL)
  {
     Set_Par_Arg(curr) ;
     curr = curr->next ;
  }
}


void Set_Iter_Arg
(ExpTreeList *arg_list)
{
  ExpTreeList *curr = arg_list->next ;

  if(arg_list->node->op == OP_ID)
    arg_list->node->op = OP_ITER_ARG ;

  while(curr != NULL)
  {
     Set_Iter_Arg(curr) ;
     curr = curr->next;
  }
}


/* remove the same repeated variables from the list
 * this is an N square algorithm, probably not good
 * enough.
 */
void collapse_var_list
(ExpTreeList* var_list)
{
  ExpTreeList* prev, *curr, *start ;

  start = var_list ;
  while(start->next != NULL)
  {
    prev = start ;
    curr = start->next ;
    while(curr->next != NULL)
    {
      if(strcmp(start->node->name,curr->node->name)
         == 0)
      {
        prev->next = curr->next ;
        curr = curr->next ;
      }
      else
      {
        curr = curr->next ;
        prev = prev->next ;
      }
    }
    start = start->next ;
  }
}

void Input_Map
(ExpTreeList* map_list)
{
  ExpTreeList* prev, *curr, *start ;
  char* name1, *name2 ;

  start = map_list ;
  while(start->next != NULL)
  {
    name1 = start->node->list1->next->node->name ;
    prev = start ;
    curr = start->next ;
    printf("(%s", start->node->list1->node->name) ;
    while(curr->next != NULL)
    {
      name2 = curr->node->list1->next->node->name ;
      if(strcmp(name1,name2) == 0)
      {
         printf(",") ;
         printf("%s", curr->node->list1->node->name) ;
         prev->next = curr->next ;
         curr = curr->next ;
      }
      else
      {
         curr = curr->next ;
         prev = prev->next ;
      }
    }
    printf(")->%s", name1) ;
    if(start->next->next != NULL)
      printf(",") ;
    /*    printf(")->%s,", name1) ;*/
    start = start->next ;
  }
  /*  printf("\b") ;*/
}

int stringcmp(char s1[], char s2[])
{
  int i ;
  int l1 ; int l2 ;
  l1 = strlen(s1) ;
  l2 = strlen(s2) ;
  if(l1 != l2) return false ;
  
  for (i=0; i<l1; i++)
    if (s1[i] != s2[i])
      return false ;
  return true ;
}

dp_id Set_dp_id(char* utname)
{
  if(stringcmp(utname, "param") == true)
    return DP_PARAM ;
  else if(stringcmp(utname, "store") == true)
    return DP_STORE ;
  else if(stringcmp(utname, "storeVec") == true)
    return DP_STOREVEC ;
  else if(stringcmp(utname, "storeMat") == true)
    return DP_STOREMAT ;
  else if(stringcmp(utname, "multiStore") == true)
    return DP_MULTISTORE ;
  else if(stringcmp(utname, "dstore") == true)
    return DP_DSTORE ;
  else if(stringcmp(utname, "constraint") == true)
    return DP_CONSTRAINT ;
  else if(stringcmp(utname, "Map") == true)
    return DP_MAP ;
  else if(stringcmp(utname, "MapVec") == true)
    return DP_MAPVEC ;
  else if(stringcmp(utname, "multiMap") == true)
    return DP_MULTIMAP ;
  else if(stringcmp(utname, "dMap") == true)
    return DP_DMAP ;
  else if(stringcmp(utname, "dMapVec") == true)
    return DP_DMAPVEC ;
  else if(stringcmp(utname, "dmultiMap") == true)
    return DP_DMULTIMAP ;
  else
    return DP_USER ;
}

/**** end-of-file ****/


