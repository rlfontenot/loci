#include <Tools/eventNotify.h>

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
    warn(!p) ;
    unsigned long notify_tag =  p->ident() ;
    warn(notify_group.find(notify_tag) != notify_group.end()) ;
    notify_group[notify_tag] = p ;
}

void eventDispatcher::disengage(eventNotify *p)
{
    warn(!p) ;

    notify_list::iterator di = notify_group.find(p->ident()) ;
    fatal(di == notify_group.end()) ;
    notify_group.erase(di) ;
//    int num_erase = notify_group.erase(p->ident()) ;
//    warn(num_erase != 1) ;
}    

void eventDispatcher::dispatch_notify()
{
    notify_list::iterator nlp ;

    for(nlp=notify_group.begin();nlp!=notify_group.end();++nlp) 
      (*nlp).second->notification() ;
}

eventNotify::~eventNotify()
{
}

}
