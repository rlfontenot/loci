#ifndef VISITORABS_H
#define VISITORABS_H

// the abstract base class of visitor
// used as interface for other files (comp_dag.cc etc)
namespace Loci {
  
  class loop_compiler ;
  class dag_compiler ;
  class conditional_compiler ;
  
  class visitor {
  public:
    virtual ~visitor() {}
    virtual void visit(loop_compiler&) = 0 ;
    virtual void visit(dag_compiler&) = 0 ;
    virtual void visit(conditional_compiler&) = 0 ;
  } ;
}

#endif
