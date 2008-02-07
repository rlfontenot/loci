//#############################################################################
//#
//# Copyright 2008, Mississippi State University
//#
//# This file is part of the Loci Framework.
//#
//# The Loci Framework is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The Loci Framework is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# Lesser GNU General Public License for more details.
//#
//# You should have received a copy of the Lesser GNU General Public License
//# along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################
#ifndef LMUTEX_H
#define LMUTEX_H

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
