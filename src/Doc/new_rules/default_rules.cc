// This is an example file for using the default_rule type.
// A default rule is mainly used to set up a default parameter
// in the fact_db. If the parameter is specified in the .vars
// configuration file, then its value will be updated during
// the read in process. Otherwise its value will remain unchanged.

// Here are some examples to use the default_rule.
// NOTE: a default_rule should only have target facts and NO source
// facts. In other words, no "input" methods in the constructor
// should be specified. The targets should be of type "param."
// A compute method is usually required to set up the default value.

class default_plot_modulo : public default_rule {
  param<int> plot_modulo ;
public:
  default_plot_modulo() {
    name_store("plot_modulo",plot_modulo) ;
    output("plot_modulo") ;
  }
  virtual void compute(const sequence &seq) {
    *plot_modulo = 0 ;
  }
} ;
register_rule<default_plot_modulo> register_default_plot_modulo ;

class default_time_integration: public default_rule {
  param<std::string> time_integration ;
public:
  default_time_integration() {
    name_store("time_integration",time_integration) ;
    output("time_integration") ;
  }
  virtual void compute(const sequence& seq) {
    *time_integration = "second_order" ;
  }
} ;
register_rule<default_time_integration> register_default_time_integration ;


// the end
