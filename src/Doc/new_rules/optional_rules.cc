// This is an example file for using the optional_rule type.
// An optional rule is mainly used to register the type info.
// of a parameter in the fact_db. If the parameter is specified
// in the .vars configuration file, then it will be created
// in the fact_db and its value will be read in. Otherwise if
// the user does not specify it in the .vars file, then it
// will NOT appear in the fact_db.

// Here are some examples to use the optinal_rule.
// NOTE: an optional_rule should only have target facts and NO source
// facts. In other words, no "input" methods in the constructor
// should be specified. The targets should be of type "param."
// Since the purpose of optional_rule is to register the type, we don't
// need the compute method. But since "compute" is pure virtual in the
// base class, it'll have to be defined. We just need to set it to
// be empty.

class optional_restart_freq: public optional_rule {
  param<int> restart_freq ;
public:
  optional_restart_freq() {
    name_store("restart_freq",restart_freq) ;
    output("restart_freq") ;
  }
  virtual void compute(const sequence& seq) {}
} ;
register_rule<optional_restart_freq> register_optional_restart_freq ;

class optional_bc_info: public optional_rule {
  param<bc_options> bc_info ;
public:
  optional_bc_info() {
    name_store("boundary_conditions",bc_info) ;
    output("boundary_conditions") ;
  }
  virtual void compute(const sequence& seq) {}
} ;
register_rule<optional_bc_info> register_optional_bc_info ;


// the end
