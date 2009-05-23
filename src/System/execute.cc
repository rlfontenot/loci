#include <ostream>
#include <execute.h>
#include "loci_globs.h"
#include <distribute.h>

using std::cerr ;
using std::endl ;
using std::ostream ;


namespace Loci {
  struct exec_info  {
    Loci::executeP exec_routine;
    Loci::fact_db *current_fact_db ;
    exec_info() {} ;
    exec_info(Loci::executeP &ep, Loci::fact_db &facts) {
      exec_routine = ep ; current_fact_db = &facts ;
    }
  } ;

  Loci::fact_db *current_fact_db ;

  void execute_list::execute(fact_db &facts, sched_db& scheds) {
    std::vector<executeP>::iterator eli ;
    for(eli=elist.begin();eli!=elist.end();++eli)
      (*eli)->execute(facts, scheds) ;
  }

  void execute_list::Print(std::ostream &s) const {
    std::vector<executeP>::const_iterator eli ;
    for(eli=elist.begin();eli!=elist.end();++eli)
      (*eli)->Print(s) ;
  }

  void execute_list::dataCollate(collectData &data_collector) const {
    std::vector<executeP>::const_iterator eli ;

    for(eli=elist.begin();eli!=elist.end();++eli)
      (*eli)->dataCollate(data_collector) ;
  }


  void execute_sequence::execute(fact_db &facts, sched_db& scheds) {
    std::vector<executeP>::iterator eli ;

    for(eli=elist.begin();eli!=elist.end();++eli)
      (*eli)->execute(facts, scheds) ;
  }

  void execute_sequence::Print(std::ostream &s) const {
    std::vector<executeP>::const_iterator eli ;
    for(eli=elist.begin();eli!=elist.end();++eli)
      (*eli)->Print(s) ;
  }

  void execute_sequence::dataCollate(collectData &data_collector) const {
    std::vector<executeP>::const_iterator eli ;

    for(eli=elist.begin();eli!=elist.end();++eli)
      (*eli)->dataCollate(data_collector) ;
  }

}
