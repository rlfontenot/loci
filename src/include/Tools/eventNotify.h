#ifndef EVENTNOTIFY_H
#define EVENTNOFITY_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Tools/debug.h>
#include <Tools/identity.h>
#include <Tools/lmutex.h>

#include <list>

namespace Loci {
    
  class eventNotify : public Identity {
  public:
    eventNotify() { }
    virtual ~eventNotify() ;
    virtual void notification()  = 0 ;
  } ;


  class eventDispatcher {
    typedef std::list<eventNotify *> notify_list ;
    notify_list notify_group ;
    lmutex mutex ;
  public:
    eventDispatcher() ;
    ~eventDispatcher() ;
    void engage(eventNotify *) ;
    void disengage(eventNotify *) ;
    void dispatch_notify() ;
  } ;
}

#endif
