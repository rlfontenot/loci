#ifndef LMUTEX_H
#define LMUTEX_H

#ifdef PTHREADS
#include <pthread.h>
#include <errno.h>

namespace Loci {

  class lmutex {
    pthread_mutex_t mutex ;
  public:
    lmutex() {
      pthread_mutexattr_t mattr ;
      pthread_mutexattr_init(&mattr) ;
      pthread_mutexattr_setpshared(&mattr, PTHREAD_PROCESS_SHARED) ;
      pthread_mutex_init(&mutex,&mattr) ;
      pthread_mutexattr_destroy(&mattr) ;
    }
    ~lmutex() {
      pthread_mutex_destroy(&mutex) ;
    }
    void lock() {
      const int err = pthread_mutex_lock(&mutex) ;
      fatal(err==EDEADLK || err==EINVAL) ;
    }
    void unlock() {
      const int err = pthread_mutex_unlock(&mutex) ;
      fatal(err==EPERM) ;
    }
    
  } ;

  class bmutex {
    lmutex &mutex ;
  public:
    bmutex(lmutex &m) : mutex(m) { mutex.lock() ; }
    ~bmutex() { mutex.unlock() ; }
  } ;
}
#else
namespace Loci {

  class lmutex {
  public:
    lmutex() {}
    void lock() {}
    void unlock() {}
  } ;

  class bmutex {
    lmutex &mutex ;
  public:
    bmutex(lmutex &m) : mutex(m) { mutex.lock() ; }
    ~bmutex() { mutex.unlock() ; }
  } ;
}
#endif    

#endif
