#include <Tools/eventNotify.h>

#include <vector>

namespace Loci {
  eventDispatcher::eventDispatcher() {}

  eventDispatcher::~eventDispatcher() {
    warn(notify_group.size()!=0) ;
  }

  void eventDispatcher::engage(eventNotify *p) {
    //    bmutex l(mutex) ;
    mutex.lock() ;
    warn(!p) ;
    notify_group.push_back(p) ;
    mutex.unlock() ;
  }

  void eventDispatcher::disengage(eventNotify *p) {
    //    bmutex l(mutex) ;
    using std::cerr ;
    using std::endl ;
    
    mutex.lock() ;
    warn(!p) ;

    if(notify_group.begin() == notify_group.end()) {
      mutex.unlock() ;
      cerr << "disengage with empty list" << endl ;
      cerr << "p = " << p << endl ;
      return ;
    }
    notify_list::iterator nlp = notify_group.end() ;
    for(--nlp;nlp != notify_group.begin() && *nlp != p;--nlp)
      /* NULL STATEMENT */ ;

    warn(*nlp != p) ;
    if(*nlp != p) {
      cerr << "list = " << endl ;
      for(nlp=notify_group.begin();nlp!=notify_group.end();++nlp)
        cerr << *nlp << " " ;
      cerr << endl << "p = " << p << endl ;
    }
        
        
    if(*nlp == p)
      notify_group.erase(nlp) ;
    mutex.unlock() ;
  }    

  void eventDispatcher::dispatch_notify() {
    //    bmutex l(mutex) ;
    mutex.lock() ;

    notify_list copy ;
    notify_list::iterator nlp ;
    for(nlp=notify_group.begin();nlp!=notify_group.end();++nlp)
      copy.push_back(*nlp) ;
    for(nlp=copy.begin();nlp!=copy.end();++nlp) 
      (*nlp)->notification() ;
    mutex.unlock() ;
  }

  eventNotify::~eventNotify() {}

}
