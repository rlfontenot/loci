#include "execute.h"

namespace Loci {
  void execute_list::execute(fact_db &facts) {
    std::list<executeP>::iterator eli ;
    for(eli=elist.begin();eli!=elist.end();++eli)
      (*eli)->execute(facts) ;
  }

  void execute_list::Print(std::ostream &s) const {
    std::list<executeP>::const_iterator eli ;
    for(eli=elist.begin();eli!=elist.end();++eli)
      (*eli)->Print(s) ;
  }
  
}
