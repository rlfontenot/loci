
#ifndef SCHEDULER_H
#define SCHEDULER_H

#include <rule.h>
#include <fact_db.h>
#include <execute.h>


namespace Loci {

  extern executeP create_execution_schedule(rule_db &rdb, fact_db &facts,
                                            std::string target) ;
}

#endif
