#ifndef EVENTNOTIFY_H
#define EVENTNOFITY_H

#include <Tools/debug.h>
#include <Tools/identity.h>


#include <map>

namespace Loci {
    
class eventNotify : public Identity {
  public:
    eventNotify() { }
    virtual ~eventNotify() ;
    virtual void notification()  = 0 ;
} ;


class eventDispatcher {
    typedef std::map<unsigned long,eventNotify *> notify_list ;
    notify_list notify_group ;
  public:
    eventDispatcher() ;
    ~eventDispatcher() ;
    void engage(eventNotify *) ;
    void disengage(eventNotify *) ;
    void dispatch_notify() ;
} ;
}

#endif
