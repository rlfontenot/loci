#include <execute.h>

#ifdef PTHREADS
#include <pthread.h>
#include <semaphore.h>
#endif

#include <Tools/stream.h>

#ifdef PTHREADS
sem_t thread_barrier, barrier_ack ;
#endif

namespace Loci {
  struct exec_info  {
    Loci::executeP exec_routine;
    Loci::fact_db *current_fact_db ;
    exec_info() {} ;
    exec_info(Loci::executeP &ep, Loci::fact_db &facts) {
      exec_routine = ep ; current_fact_db = &facts ;
    }
  } ;
  std::vector<exec_info> thread_schedule[Loci::max_threads] ;
  int num_created_threads = 1 ;
  bool work_in_queue = false ;
  
  int thread_num[Loci::max_threads] ;

  std::vector<Loci::executeP> *current_execute_list ;
  Loci::fact_db *current_fact_db ;
}

void process_thread(int i) {
  using namespace Loci ;
  for(int w=0;w<thread_schedule[i].size();++w) {
    exec_info &ei = thread_schedule[i][w] ;
    ei.exec_routine->execute(*ei.current_fact_db) ;
  }
}

#ifdef PTHREADS

extern "C" {

  void *worker_thread( void *ptr) {
    using namespace Loci ;
    int tnum = *(int *)(ptr) ;
    (*current_execute_list)[tnum]->execute(*current_fact_db) ;
    pthread_exit(0) ;
    return NULL ;
  }

  void *worker_thread2( void *ptr) {
    int tnum = *(int *)(ptr) ;
    process_thread(tnum) ;
    pthread_exit(0) ;
    return NULL ;
  }

#ifdef TRYOUT
  void *worker_thread3( void *ptr) {
    int tnum = *(int *)(ptr) ;
    do {
      while(sem_wait(&thread_barrier) != 0)
        if(errno != EINTR) {
          perror("sem_wait") ;
        }
      process_thread(tnum) ;
      if(sem_post(&barrier_ack) != 0)
        perror("sem_post") ;
      
    } while(true) ;
  }
#endif
}


#endif

namespace Loci {
  void execute_list::execute(fact_db &facts) {
    static int round_robin_allocate = 0 ;
    std::vector<executeP>::iterator eli ;
#ifdef PTHREADS
    for(eli=elist.begin();eli!=elist.end();++eli) {
      if(!(*eli)->is_control_thread() && num_created_threads>1) {
        work_in_queue = true ;
        thread_schedule[round_robin_allocate].push_back(exec_info(*eli,facts)) ;
        if(num_created_threads > 1)
          round_robin_allocate = (round_robin_allocate+1)%num_created_threads ;
      } else
        (*eli)->execute(facts) ;
    }
#else
    for(eli=elist.begin();eli!=elist.end();++eli)
      (*eli)->execute(facts) ;
#endif    
  }

  void execute_list::Print(std::ostream &s) const {
    std::vector<executeP>::const_iterator eli ;
    for(eli=elist.begin();eli!=elist.end();++eli)
      (*eli)->Print(s) ;
  }

  void execute_sequence::execute(fact_db &facts) {
    std::vector<executeP>::iterator eli ;

    for(eli=elist.begin();eli!=elist.end();++eli)
      (*eli)->execute(facts) ;
  }

  void execute_sequence::Print(std::ostream &s) const {
    std::vector<executeP>::const_iterator eli ;
    for(eli=elist.begin();eli!=elist.end();++eli)
      (*eli)->Print(s) ;
  }

  void execute_par::execute(fact_db &facts) {
#ifdef PTHREADS
    for(int i=0;i!=elist.size();++i) {
      work_in_queue = true ;
      thread_schedule[i%num_created_threads].
        push_back(exec_info(elist[i],facts)) ;
    }
#else
    std::vector<executeP>::iterator eli ;
    for(eli=elist.begin();eli!=elist.end();++eli)
      (*eli)->execute(facts) ;
#endif
  }

  void execute_par::Print(std::ostream &s) const {
    s << "ParBegin {" << std::endl ;
    std::vector<executeP>::const_iterator eli ;
    for(eli=elist.begin();eli!=elist.end();++eli)
      (*eli)->Print(s) ;
    s << "} ParEnd" << std::endl ;
  }


  void execute_create_threads::execute(fact_db &facts) {
    num_created_threads = num_threads ;

#ifdef TRYOUT
    sem_init(&thread_barrier,0,0) ;
    sem_init(&barrier_ack,0,0) ;
    
    static pthread_t tids[max_threads] ;
    pthread_attr_t attr ;

    pthread_attr_init(&attr) ;
    pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM) ;

    for(int i=0;i<num_created_threads;++i) {
      thread_num[i] = i ;
      pthread_create(&tids[i], &attr, &worker_thread3, &thread_num[i]) ;
    }
#endif
  }

  void execute_create_threads::Print(std::ostream &s) const {
    s << "create " << num_threads << " threads" << std::endl ;
  }

  void execute_destroy_threads::execute(fact_db &facts) {
  }

  void execute_destroy_threads::Print(std::ostream &s) const {
    s << "destroy threads" << std::endl ;
  }
  
  void execute_thread_sync::execute(fact_db &facts) {
#ifdef TRYOUT
    static bool in_barrier = false ;
    if(in_barrier) {
      cerr << "nested barrier calls, does not make sense."<<endl ;
    }
    in_barrier = true ;
    if(work_in_queue) {
      for(int i=0;i<num_created_threads;++i)
        sem_post(&thread_barrier) ;
      for(int i=0;i<num_created_threads;++i)
        sem_wait(&barrier_ack) ;
      for(int i=0;i<num_created_threads;++i)
        thread_schedule[i].clear() ;
      work_in_queue = false ;
    }
    in_barrier = false ;
    
#else
#ifdef PTHREADS
    static bool in_barrier = false ;

    if(in_barrier) {
      cerr << "nested barrier calls, does not make sense."<<endl ;
    }
    in_barrier = true ;
    if(work_in_queue) {
      pthread_t tids[max_threads] ;
      pthread_attr_t attr ;

      pthread_attr_init(&attr) ;
      pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM) ;
    
      for(int i=0;i<num_created_threads;++i) {
        thread_num[i] = i ;
        pthread_create(&tids[i], &attr, &worker_thread2, &thread_num[i]) ;
      }

      for(int i=0;i<num_created_threads;++i) {
        pthread_join(tids[i],NULL) ;
      }
      for(int i=0;i<num_created_threads;++i) {
        thread_schedule[i].clear() ;
      }
      work_in_queue = false ;
    }
    in_barrier = false ;
#endif
#endif
  }

  void execute_thread_sync::Print(std::ostream &s) const {
    s << "thread barrier " << note << std::endl ;
  }

}
