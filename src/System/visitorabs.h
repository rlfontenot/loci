#ifndef VISITORABS_H
#define VISITORABS_H

// the abstract base class of visitor
// used as interface for other files (comp_dag.cc etc)
namespace Loci {
  
  class loop_compiler ;
  class dag_compiler ;
  class conditional_compiler ;
  class chomp_compiler ;
  class impl_recurse_compiler ;
  class recurse_compiler ;
  
  class visitor {
  public:
    virtual ~visitor() {}
    // the following three are mandatory methods
    virtual void visit(loop_compiler&) = 0 ;
    virtual void visit(dag_compiler&) = 0 ;
    virtual void visit(conditional_compiler&) = 0 ;
    // the recursive compiler may not needed in the
    // visit, therefore we supply two empty defaults
    virtual void visit(impl_recurse_compiler&) {}
    virtual void visit(recurse_compiler&) {}
  } ;
}

#endif
