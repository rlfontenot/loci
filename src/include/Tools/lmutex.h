//#############################################################################
//#
//# Copyright 2008, 2015, Mississippi State University
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
namespace Loci {
  // class lmutex {
  // public:
  //   lmutex() {}
  //   void lock() {}
  //   void unlock() {}
  // } ;

  // class bmutex {
  // public:
  //   bmutex(lmutex &m) {}
  // } ;

  // class lmutex {
  //   pthread_mutex_t m;
  //   pthread_mutexattr_t at;
  // public:
  //   lmutex() { pthread_mutexattr_init(&at); pthread_mutex_init(&m,NULL); }
  //   void lock() { pthread_mutex_lock(&m); }
  //   void unlock() { pthread_mutex_unlock(&m); }
  //   ~lmutex() { pthread_mutex_destroy(&m); pthread_mutexattr_destroy(&at); }
  // } ;
  
  // class bmutex {
  //   lmutex &mutex ;
  // public:
  //   bmutex(lmutex &m) : mutex(m) { mutex.lock() ; }
  //   ~bmutex() { mutex.unlock() ; }
  // } ;

  class lmutex {
    pthread_spinlock_t m;
  public:
    lmutex() { pthread_spin_init(&m,0); }
    void lock() { pthread_spin_lock(&m); }
    void unlock() { pthread_spin_unlock(&m); }
    ~lmutex() { pthread_spin_destroy(&m); }
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
  public:
    bmutex(lmutex &m) {}
  } ;
}
#endif

#endif
