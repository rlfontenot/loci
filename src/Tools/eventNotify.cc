#include <Tools/eventNotify.h>

#include <vector>

namespace Loci {
eventDispatcher::eventDispatcher()
{
}

eventDispatcher::~eventDispatcher()
{
  warn(notify_group.begin() != notify_group.end()) ;
}

void eventDispatcher::engage(eventNotify *p)
{
  bmutex l(mutex) ;
  warn(!p) ;
  unsigned long notify_tag =  p->ident() ;
  warn(notify_group.find(notify_tag) != notify_group.end()) ;
  notify_group[notify_tag] = p ;
}

void eventDispatcher::disengage(eventNotify *p)
{
  bmutex l(mutex) ;
  warn(!p) ;
  notify_list::iterator di = notify_group.find(p->ident()) ;
  fatal(di == notify_group.end()) ;
  notify_group.erase(di) ;
}    

void eventDispatcher::dispatch_notify()
{
  mutex.lock() ;
  std::vector<eventNotify *> ln ;
  notify_list::iterator nlp ;
  for(nlp=notify_group.begin();nlp!=notify_group.end();++nlp) 
    ln.push_back((*nlp).second) ;
  mutex.unlock() ;
  for(int i=0;i<ln.size();++i)
    ln[i]->notification() ;
}

eventNotify::~eventNotify()
{
}

}
