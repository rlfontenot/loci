
#ifndef SCHEDULER_H
#define SCHEDULER_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>


#include <rule.h>
#include <fact_db.h>
#include <execute.h>


namespace Loci {

  // function prototype for internal queries
  executeP create_internal_execution_schedule(rule_db& par_rdb,
                                              fact_db &facts,
                                              const variableSet& target,
                                              int nth=1) ;
  bool internalQuery(rule_db& par_rdb, fact_db& facts,
                     const variableSet& query) ;
    

  extern executeP create_execution_schedule(const rule_db &rdb,
                                            fact_db &facts,
                                            const variableSet& target,
                                            int nth=1) ;

  extern bool makeQuery(const rule_db &rdb, fact_db &facts,
                        const std::string& target) ;
}

#endif
