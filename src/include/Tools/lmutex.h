#ifndef LMUTEX_H
#define LMUTEX_H
#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif

#include <Config/conf.h>

#ifdef PTHREADS
#include <pthread.h>
#include <errno.h>


namespace Loci {

  class lmutex {
    pthread_mutex_t mutex ;
    pthread_mutexattr_t mattr ;
  public:
    lmutex() {
      pthread_mutexattr_init(&mattr) ;
#ifdef LINUX
      // Linux doesn't have shared attribute for mutexes
#else
      pthread_mutexattr_setpshared(&mattr, PTHREAD_PROCESS_SHARED) ;
#endif
      //      pthread_mutexattr_settype(&mattr, PTHREAD_MUTEX_RECURSIVE) ;
      pthread_mutex_init(&mutex,&mattr) ;
    }
    ~lmutex() {
#ifdef DEBUG
      const int err =
#endif
        pthread_mutex_destroy(&mutex) ;
      fatal(err==EBUSY) ;
      pthread_mutexattr_destroy(&mattr) ;
    }
    void lock() {
#ifdef DEBUG
      const int err =
#endif
        pthread_mutex_lock(&mutex) ;
      fatal(err==EDEADLK || err==EINVAL) ;
    }
    void unlock() {
#ifdef DEBUG
      const int err =
#endif
        pthread_mutex_unlock(&mutex) ;
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
