// This is an example file for using the constraint_rule type.
// A constraint rule is mainly used to compute a constraint
// at runtime. The computed constraint can be used in further
// Loci rules. A constraint_rule is similar to the conventional
// pointwise_rule etc. except for that the targets are typically
// all constraint. Note: in specifying the target type, we will 
// need to use "Constraint" (notice the capital "C") instead of
// the lower case "constraint." This is due to historic reason. 
// The lower case "constraint" has been used in other context
// in Loci. 

// Below is an example constraint_rule. In this rule, we compute
// different limiters based on the input param, which is a string
// that should be specified in the .vars configuration file. 

class constraint_limiter: public constraint_rule {
  const_param<std::string> limiter ;
  Constraint V_limiter, B_limiter, N_limiter, Z_limiter ;
public:
  constraint_limiter() {
    name_store("limiter",limiter) ;
    name_store("V_limiter",V_limiter) ;
    name_store("B_limiter",B_limiter) ;
    name_store("N_limiter",N_limiter) ;
    name_store("Z_limiter",Z_limiter) ;
    input("limiter") ;
    output("V_limiter,B_limiter,N_limiter,Z_limiter") ;
  }
  virtual void compute(const sequence& seq) {
    if(*limiter == "venkatakrishnan" || *limiter == "V") {
      V_limiter = ~EMPTY ;
      /////
      B_limiter = EMPTY ;
      N_limiter = EMPTY ;
      Z_limiter = EMPTY ;
    } else if(*limiter == "barth" || *limiter == "B") {
      B_limiter = ~EMPTY ;
      /////
      V_limiter = EMPTY ;
      N_limiter = EMPTY ;
      Z_limiter = EMPTY ;
    } else if(*limiter == "none") {
      N_limiter = ~EMPTY ;
      /////
      V_limiter = EMPTY ;
      B_limiter = EMPTY ;
      Z_limiter = EMPTY ;
    } else if(*limiter == "zero") {
      Z_limiter = ~EMPTY ;
      //////
      V_limiter = EMPTY ;
      B_limiter = EMPTY ;
      N_limiter = EMPTY ;
    } else {
      cerr << "limiter " << *limiter
           << " not supported for generalized grids" << endl ;
      cerr << "defaulting to venkatakrishnan limiter" << endl ;
      V_limiter = ~EMPTY ;
      //////
      B_limiter = EMPTY ;
      N_limiter = EMPTY ;
      Z_limiter = EMPTY ;
    }
  }
} ;
register_rule<constraint_limiter> register_constraint_limiter ;

// the end

