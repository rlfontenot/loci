#include <Tools/eventNotify.h>

#include <vector>

namespace Loci {
  eventDispatcher::eventDispatcher() {}

  eventDispatcher::~eventDispatcher() {
    warn(notify_group.size()!=0) ;
  }

  void eventDispatcher::engage(eventNotify *p) {
    bmutex l(mutex) ;
    warn(!p) ;
    notify_group.push_back(p) ;
  }

  void eventDispatcher::disengage(eventNotify *p) {
    bmutex l(mutex) ;
    
    warn(!p) ;

    if(notify_group.begin() == notify_group.end()) {
      std::cerr << "disengage with empty list" << std::endl ;
      std::cerr << "p = " << p << std::endl ;
      return ;
    }
    notify_list::iterator nlp = notify_group.end() ;
    for(--nlp;nlp != notify_group.begin() && *nlp != p;--nlp)
      /* NULL STATEMENT */ ;

    warn(*nlp != p) ;
        
    if(*nlp == p)
      notify_group.erase(nlp) ;
  }    

  void eventDispatcher::dispatch_notify() {
    mutex.lock() ;
    notify_list copy ;
    notify_list::iterator nlp ;
    for(nlp=notify_group.begin();nlp!=notify_group.end();++nlp)
      copy.push_back(*nlp) ;
    mutex.unlock() ;
    for(nlp=copy.begin();nlp!=copy.end();++nlp) 
      (*nlp)->notification() ;
  }

  eventNotify::~eventNotify() {}

}
