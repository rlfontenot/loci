#include <rule.h>

namespace Loci {
  rule_db local_modify_time_vars(rule_db& rdb, const std::string &sn) ;
  rule rename_rule(rule r, std::map<variable, variable> &vm) ; 
  rule_db parametric_rdb(rule_db &rdb) ;
}
