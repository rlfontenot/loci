#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <rule.h>

namespace Loci {
  rule prepend_time_level(rule r, std::string prepend) ;
  rule rename_rule(rule r, std::map<variable, variable> &vm) ; 
  rule_db parametric_rdb(rule_db &rdb,variableSet query_vars) ;
}
