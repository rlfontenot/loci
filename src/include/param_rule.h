#ifndef PARAM_RULE_H
#define PARAM_RULE_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>


#include <rule.h>

namespace Loci {
  rule prepend_time_level(rule r, std::string prepend) ;
  rule rename_rule(rule r, std::map<variable, variable> &vm) ; 
  rule_db parametric_rdb(const rule_db &rdb,variableSet query_vars) ;
}

#endif
