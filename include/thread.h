//#############################################################################
//#
//# Copyright 2008-2019, Mississippi State University
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
#ifndef THREAD_H
#define THREAD_H

// disable this file by default for now...
#ifdef PTHREADS

#include <variable.h>
#include <rule.h>
#include <execute.h>
#include <distribute.h>
#include <Tools/cptr.h>

#include <pthread.h>
#include <semaphore.h>
#include <sched.h>

#include <vector>
#include <map>
#include <set>
#include <deque>
#include <string>
#include <exception>
#include <utility>
#include <cctype>

#include "comp_tools.h"
#include "sched_tools.h"
#include "loci_globs.h"

#ifndef PAPI_DEBUG
#define long_long long
#endif

namespace Loci {

  extern int current_rule_id;

  // a simple thread exception
  struct ThreadException : public std::exception {
    ThreadException(const std::string& s="thread exception"):msg(s) {}
    ~ThreadException() throw() {} // make g++ happy
    virtual const char* what() const throw()
    { return msg.c_str(); }
  private:
    std::string msg;
  };

  // this is the abstract interface for all the thread execution modules
  // all original sequential execution modules that could be made
  // multi-threaded would generate corresponding threaded version as
  // described later in this file
  class ThreadedExecutionModule : public CPTR_type {
  public:
    virtual void execute(fact_db& facts, sched_db& scheds) = 0;
    // print execution sequence for the current thread
    virtual void print_seq(std::ostream& s) const = 0;
    // return the timing information of this execution module
    virtual double get_time() const = 0;
    // get cache misses
    virtual long_long get_l1_dcm() const = 0;
    virtual long_long get_l2_dcm() const = 0;
    virtual ~ThreadedExecutionModule() {}
  };
  typedef CPTR<ThreadedExecutionModule> ThreadedEMP;

  // an abstract interface for thread control.
  // this is used to globally control threads behavior, such as
  // creating/starting/stopping threads, do thread sync, etc.
  //
  // this interface is only used by the main control thread, which
  // does NOT participate in the work performed by all the work
  // threads created through this interface.
  class ThreadControl {
  public:
    // this function creates all the work threads
    // (the number of the threads is determined at the construction time)
    // Each work thread will be in a suspended state
    // upon creation (essentially sits there idly waiting to be fed
    // a thread execution module). 
    virtual void create_threads() = 0;
    // this function terminates all created work threads.
    virtual void shutdown() = 0;
    // this function restarts all suspended threads, feeding each
    // a new thread execution module.
    virtual void restart(std::vector<ThreadedEMP>& tem) = 0;
    // this function causes all the work threads to be sequentialized
    // i.e., thread 0 executes first, and then thread 1 executes,
    // and so on... this is mainly supplied as a debugging tool.
    virtual void sequential_restart(std::vector<ThreadedEMP>& tem) = 0;
    // this function causes the main control thread to suspend
    // waiting for all the work-thread to finish their current
    // execution module, and then wakes up.
    // It also means that after this function returns, all working
    // threads will be in a suspended state.
    virtual void wait_threads() = 0;
    // returns the total number of work threads created so far
    virtual int num_threads() const = 0;
    // sets up the the global fact and sched database for use in execution.
    virtual void setup_facts_scheds(fact_db&, sched_db&) = 0;
    // this just returns the minimum work per thread (as entity size)
    // i.e., if the total work is less than this value times the
    // number of threads, then it is not worth to lauch the work on
    // multiple threads.
    virtual int min_work_per_thread() const = 0;
    // given a thread id, this function will return the reference
    // to a vector of thread id that represents the (global) reduction
    // partner of the given thread
    virtual const std::vector<int>& get_reduction_partners(int id) const = 0;
    // returns the thread id of the reduction tree root
    virtual int get_reduction_root() const = 0;

    // this one determines if a calling thread is the leading thread
    virtual bool is_leading_thread() const = 0;

    // these provide a global lock on all threads managed by 
    // the thread control
    virtual void atomic_begin() = 0;
    virtual void atomic_end() = 0;

    // this method provides the local ID for each calling thread
    virtual std::string get_local_id() const = 0;

    virtual ~ThreadControl() {}
  };

  extern ThreadControl* thread_control;

  // this function partitions a sequence into n parts
  std::vector<sequence> partition_seq(const sequence& s, int n);
  
  // this is an implementation of the ThreadControl interface
  // based on the POSIX thread and semaphore implementation.
  class ThreadControl_pthread : public ThreadControl {
  public:
    ThreadControl_pthread(int n);
    ~ThreadControl_pthread();
    void create_threads();
    void shutdown();
    void restart(std::vector<ThreadedEMP>&);
    void sequential_restart(std::vector<ThreadedEMP>&);
    void wait_threads();
    int num_threads() const { return tnum; }
    void setup_facts_scheds(fact_db& facts, sched_db& scheds)
    { factsP = &facts; schedsP = &scheds; }
    // 10 is just an arbitrary value for now...
    // we may want to define a better measure in the future...
    int min_work_per_thread() const { return 10; }
    const std::vector<int>& get_reduction_partners(int id) const;
    int get_reduction_root() const { return 0; }
    bool is_leading_thread() const;
    void atomic_begin() { global_spin.lock(); }
    void atomic_end()  { global_spin.unlock(); }
    std::string get_local_id() const;
  private:
    // disable copy and assignment
    ThreadControl_pthread(const ThreadControl_pthread&);
    ThreadControl_pthread& operator=(const ThreadControl_pthread&);
    // the thread function assigned to eacher worker
    static void* thread_fun(void* arg);
    // whether the work threads are active or not.
    bool active;
    // flag to indicate shutdown
    bool stop;
    // this would point to the actual work for each work thread
    std::vector<ThreadedEMP>* work;
    // pointers to the current fact and sched database
    fact_db* factsP;
    sched_db* schedsP;
    // a small structure to build the argument to the thread function
    struct TArg {
      int tid;                // internal thread id, from 0 to n-1
      pthread_t ptid;           // thread id assigned by the pthread lib
      ThreadControl_pthread* self; // pointer to the class object created
                                   // used to access the object internals
                                   // from each thread
    };
    // records the number of threads created
    int tnum;
    std::vector<TArg> thread_args;
    // the pthreads attributes when creating threads
    pthread_attr_t pattr;
    // semaphores used to block the work threads
    std::vector<sem_t> worker_block;
    // semaphore used to block the main thread
    sem_t main_block;
    // this is used to parallel start all threads
    std::vector<std::vector<int> > start_partners;
    // these types and data structures are used to signal the
    // finishing of all the work threads
    struct Terminate {
      bool done;           // finish flag
      pthread_spinlock_t done_lock; // spin locks that guard "done"
      int p; // the parent node
      std::vector<int> c; // all the child nodes, could be 0, 1, or 2
    };
    std::vector<Terminate> finish_notification;
    std::vector<vector<int> > reduction_partners;
    // the pthread id of the main and the leading work thread
    pthread_t main_pid;
    pthread_t lead_work_pid;
    lmutex global_spin;
    std::map<pthread_t,std::string> id_map;
    void build_term_tree(int pid, int depth, size_t& myid);
  };

  // this is a undirected graph type used in the graph coloring algorithm
  // the graph is assumed to have only one edge between each vertices pair
  // vertices are represented using integer sets.
  class UDG {
  public:
    // add an edge between vertices x and y. if x and/or y do not exist,
    // first create those vertices in the graph.
    void add_edge(int x, int y);
    // add an edge between x and every vertex in ys.
    void add_edges(int x, const entitySet& ys);
    // add edges to all pairs in the input vertices
    void add_edges(const entitySet& xs);
    // delete the link between x and y.
    void delete_edge(int x, int y);
    // delete all the edges between x and every vertex in ys.
    void delete_edges(int x, const entitySet& ys);
    // add a single vertex with no edges attached
    // if the vertex is already in the graph, then this call has
    // no effects.
    void add_vertex(int x);
    // add a set of vertices with no edges attached to them.
    // if any of them is already in the graph, then nothing is done
    // to that vertex.
    void add_vertices(const entitySet& xs);
    // remove the vertex x from the graph, also deletes
    // all edges connecting to x.
    void delete_vertex(int x);
    // remove all vertices in the set xs.
    void delete_vertices(const entitySet& xs);
    // contract vertex y into vertex x, i.e., remove y first, and then
    // make all neighbors of y the neighbor of x.
    void contract_vertex(int x, int y);
    // returns the vertices connected to vertex x.
    const entitySet& neighbors(int x) const;
    // returns the vertices that are NOT directly connected to x.
    entitySet non_neighbors(int x) const;
    // returns all of the vertices in the graph.
    const entitySet& all_vertices() const { return vertices; }
    // return all of the edges in the graph.
    std::vector<pair<int,int> > all_edges() const;
    // return all vertices with no connection to others (dangling vertices)
    entitySet dangling_vertices() const;
    // this function will compute the vertex with the maximum degree
    // i.e., the vertex that has the maximum number of neighbors.
    // it returns a reference to the neighbor set and also sets the
    // vertex id in the passed in parameter.
    const entitySet& max_degree(int& x) const;
    // this is an overloaded version of max_degree, in which the vertex
    // with the maximum degree is chosen from the passed in vertices set.
    // it returns a reference to the vertex with maximum degree and also
    // sets the vertex id in the passed in parameter.
    const entitySet& max_degree(const entitySet& xs, int& x) const;
    // this function returns the set of common neighbors between x and y.
    entitySet common_neighbors(int x, int y) const;
    // produce some statistical information about the graph
    std::ostream& show_stats(std::ostream& s) const;
  private:
    entitySet vertices;
    std::map<int,entitySet> g;
    friend std::ostream& operator<<(std::ostream& s, const UDG& g);
    friend std::map<int,int> lf_color(const UDG& g);
  };
  // io functions
  std::istream& operator>>(std::istream& s, UDG& g);
  std::ostream& operator<<(std::ostream& s, const UDG& g);
  
  class UDGRangeError: public std::exception {
  public:
    virtual const char* what() const throw()
    { return "UDG: vertices range error!"; }
  };
  
  class UDGFormatError: public std::exception {
  public:
    virtual const char* what() const throw()
    { return "UDG: input format error!"; }
  };

  // the "Recursive Largest First" graph color algorithm
  // Notice: the input graph is modified!
  std::map<int,int> rlf_color(UDG& g);

  // "largest first" graph color algorithm
  std::map<int,int> lf_color(const UDG& g);

  // function to verify the correctness of graph coloring
  bool verify_color(const std::map<int,int>& color, const UDG& g);

  // this function generates an interference graph for
  // a local reduction rule.
  void build_interference_graph(const rule& r,
                                fact_db& facts, sched_db& scheds,
                                const entitySet& context, UDG& g);
  
  // partition entity set based on colors
  std::vector<entitySet> partition_color(const std::map<int,int>& color);

  // this execution module creates all work threads
  class StartThreads : public execute_modules {
  public:
    StartThreads() {}
    void execute(fact_db& facts, sched_db& scheds)
    {
      stopWatch s;
      s.start();
      thread_control->create_threads();
      timer.addTime(s.stop(),1);
    }
    void Print(std::ostream& s) const
    {
      printIndent(s);
      s << "starting threads" << std::endl;
    }
    std::string getName() { return "StartThreads"; }
    void dataCollate(collectData& data_collector) const
    {data_collector.accumulateTime(timer,EXEC_CONTROL,"start_threads");}
  protected:
    timeAccumulator timer;
  };
  
  // this execution module causes all work threads to be shutdown,
  // and the thread control object to be deleted.
  class ShutDownThreads : public execute_modules {
  public:
    void execute(fact_db& facts, sched_db& scheds)
    {
      stopWatch s;
      s.start();
      thread_control->shutdown();
      delete thread_control;
      thread_control = 0;
      timer.addTime(s.stop(),1);
    }
    void Print(std::ostream& s) const
    {
      printIndent(s);
      s << "shutting down threads" <<std::endl;
    }
    std::string getName() { return "ShutDownThreads"; }
    void dataCollate(collectData& data_collector) const
    {data_collector.accumulateTime(timer,EXEC_CONTROL,"shutdown_threads");}
  protected:
    timeAccumulator timer;
  };

  // we will then create threaded execution modules for the various
  // Loci computations that can be made multithreaded.

  // this is the module for threading pointwise rules
  class ExecuteThreaded_pointwise: public ThreadedExecutionModule {
  public:
    ExecuteThreaded_pointwise(const sequence& s, executeP er)
      :seq(s), exec_size(s.size()), exec_rule(er) {}
    void execute(fact_db& facts, sched_db& scheds)
    { 
      stopWatch s;
      s.start();
      exec_rule->execute_kernel(seq);
      timer.addTime(s.stop(), exec_size);
    }
    void print_seq(std::ostream& s) const
    {
      if(verbose || seq.num_intervals() < 4) {
        s << seq;
      } else {
        s << "{" << seq[0] << "...}";
      }
    }
    double get_time() const { return timer.getTime(); }
    long_long get_l1_dcm() const { return 0; } // to be implemented!!
    long_long get_l2_dcm() const { return 0; } // to be implemented!!
  private:
    // the execution sequence
    sequence seq;
    // the size of the execution sequence
    size_t exec_size;
    // we will just delegate the execution to the underlying concrete rule
    executeP exec_rule;
    // a timer to record the execution time
    timeAccumulator timer;
  };

  // this is a module that hooks up the compiler and the threaded
  // execution module for threaded pointwise rules defined above
  class Threaded_execute_rule: public execute_modules {
  public:
    Threaded_execute_rule(rule r, const sequence& s,
                          fact_db& facts, sched_db& scheds)
      :rule_name(r), seq(s)
    {
      exec_rule = new execute_rule(r, s, facts, scheds);
      int tnum = thread_control->num_threads();
      std::vector<sequence> partition = partition_seq(seq, tnum);
      tot_mem = 0;
      for(int i=0;i<tnum;++i) {
        texec.push_back(new ExecuteThreaded_pointwise(partition[i], 
                                                      exec_rule));
        tot_mem += 16*partition[i].num_intervals();
      }
      tot_mem += 8*seq.num_intervals();
      // DEBUG
      // Loci::debugout << "--------" << r << std::endl;
      // Loci::debugout << "total seq intervals: "
      //                << seq.num_intervals()
      //                << " size: " << seq.size() << std::endl;
      // for(int i=0;i<tnum;++i)
      //   Loci::debugout << "thread " << i << " intervals: "
      //                  << partition[i].num_intervals()
      //                  << " size: " << partition[i].size() << std::endl;
      ////////
    }
    void execute(fact_db& facts, sched_db& scheds)
    {
      current_rule_id = rule_name.ident();
      stopWatch s;
      s.start();

      exec_rule->execute_prelude(seq);
      timer_comp.addTime(s.stop(), seq.size());

      s.start();
      thread_control->setup_facts_scheds(facts,scheds);
      thread_control->restart(texec);
      timer_ctrl.addTime(s.stop(),seq.size());

      thread_control->wait_threads();

      s.start();
      exec_rule->execute_postlude(seq);
      timer_comp.addTime(s.stop(), seq.size());
      current_rule_id = 0;
    }
    void Print(std::ostream& s) const
    {
      printIndent(s);
      s << rule_name << " over sequence ";
      if(verbose || seq.num_intervals() < 4)
        s << seq;
      else
        s << "{" << seq[0] << "...}";
      s << " using " << texec.size() << " threads";
      if(verbose) {
        s << " (";
        s << "thread " << 0 << ": ";
        texec[0]->print_seq(s);
        for(size_t i=1;i<texec.size();++i) {
          s << ", thread " << i << ": " ;
          texec[i]->print_seq(s);
        }
        s << ")";
      }
      s << std::endl;
    }

    std::string getName() { return "Threaded_execute_rule"; }

    void dataCollate(collectData& data_collector) const
    {
      std::ostringstream oss;
      oss << "rule: " << rule_name;
      // control time
      data_collector.accumulateTime(timer_ctrl, EXEC_CONTROL, oss.str());
      // we also need to add the computation time, which is the max
      // time from all threads plus the compute time at this level
      double mt = 0;
      for(size_t i=0;i<texec.size();++i) {
        double lt = texec[i]->get_time();
        if (lt > mt)
          mt = lt;
      }
      timeAccumulator new_comp = timer_comp;
      new_comp.addTime(mt, seq.size());
      data_collector.accumulateTime(new_comp, 
                                    EXEC_COMPUTATION, oss.str());
      // collect scheduling memory footprint size
      data_collector.accumulateSchedMemory(oss.str(),tot_mem);
      // DEBUG (report individual thread timing in debug file
      Loci::debugout << "--------" << rule_name << std::endl;
      for(size_t i=0;i<texec.size();++i) {
        Loci::debugout << "thread " << i << " time: "
                       << texec[i]->get_time() << std::endl;
      }
      ////////
      // collect cache misses
      long_long total_l1_dcm = 0;
      long_long total_l2_dcm = 0;
      for(size_t i=0;i<texec.size();++i) {
        total_l1_dcm += texec[i]->get_l1_dcm();
        total_l2_dcm += texec[i]->get_l2_dcm();
      }
      data_collector.accumulateDCM(oss.str(),total_l1_dcm,total_l2_dcm);
      // report individual thread cache misses in debug file
      Loci::debugout << "********" << rule_name << std::endl;
      for(size_t i=0;i<texec.size();++i) {
        long_long l1 = texec[i]->get_l1_dcm();
        long_long l2 = texec[i]->get_l2_dcm();
        Loci::debugout << "thread " << i << " L1 DCM: " << l1
                       << ", L2 DCM: " << l2 << std::endl;
      }
    }
  protected:
    timeAccumulator timer_comp; // computatation timer
    timeAccumulator timer_ctrl; // control timer
    rule rule_name;
    sequence seq;
    std::vector<ThreadedEMP> texec;
    double tot_mem;             // tabulating the total memory
    executeP exec_rule; // the serial execution module
  };

  // a type used in threaded global reduction
  struct ParamReductionSignal {
    bool done;
    pthread_spinlock_t done_lock;
  };

  // threaded execution module for global reduction computation
  class ExecuteThreaded_param_reduction: public ThreadedExecutionModule {
  public:
    ExecuteThreaded_param_reduction
    (int i, const sequence& s, executeP er,
     CPTR<joiner> jop, const std::vector<int>& rp,
     std::vector<ParamReductionSignal>& ds,
     std::vector<storeRepP>& pr, storeRepP f)
      :tid(i),seq(s),exec_size(s.size()),exec_rule(er),
       join_op(jop),partners(rp),done_signals(ds),partial(pr),
       final_s(f) {}
    
    void execute(fact_db& facts, sched_db& scheds)
    {
      stopWatch s;
      s.start();
      // it is important that we should copy the unit (initial) value
      // to the local partial result before executing the rule
      partial[tid]->fast_copy(final_s, EMPTY);//final_s->domain());//entitySet(seq));
      //partial[tid]->copy(final_s, final_s->domain());//entitySet(seq));
      // we will then execute the rule
      exec_rule->execute_kernel(seq);
      // then we will need to perform the reduction in parallel.
      // we will just need to reverse the termination tree order
      // to perform the reduction
      for(size_t i=0;i<partners.size();++i) {
        int p = partners[i];
        bool done = false;
        while(!done) {
          pthread_spin_lock(&(done_signals[p].done_lock));
          done = done_signals[p].done;
          pthread_spin_unlock(&(done_signals[p].done_lock));
        }
        // we will then reset the done flag in the child
        // we don't need to acquire the lock at this time since
        // this child should have finished by now.
        // the last finished thread's done signal is reset
        // by the main thread when the final result is read
        done_signals[p].done = false;
        // perform the join operation
        join_op->SetArgs(partial[tid], partial[p]);
        join_op->Join(seq);
      }
      // after all join is done, we will then write the local done signal
      pthread_spin_lock(&(done_signals[tid].done_lock));
      done_signals[tid].done = true;
      pthread_spin_unlock(&(done_signals[tid].done_lock));
      timer.addTime(s.stop(),exec_size);
    }
    void print_seq(std::ostream& s) const
    {
      if(verbose || seq.num_intervals() < 4) {
        s << seq;
      } else {
        s << "{" << seq[0] << "...}";
      }
    }
    double get_time() const { return timer.getTime(); }
    long_long get_l1_dcm() const { return 0; } // not implemented yet
    long_long get_l2_dcm() const { return 0; }
  private:
    // thread id
    int tid;
    // the execution sequence (for printing purpose only)
    sequence seq;
    size_t exec_size;
    // we will just delegate the main computation to the underlying
    // concrete rule.
    executeP exec_rule;
    // the join operator
    CPTR<joiner> join_op;
    // the reduction partners of this thread
    const std::vector<int>& partners;
    // this is the done signal pointer
    std::vector<ParamReductionSignal>& done_signals;
    // pointer to the partial results of every thread
    std::vector<storeRepP>& partial;
    // this is the pointer to the store that the final result should go
    // here it is used to copy the initial value
    storeRepP final_s;
    timeAccumulator timer;
  };  

  // this is a module that hooks up the compiler and the
  // threaded execution module for param reduction rule defined above
  class Threaded_execute_param_reduction: public execute_modules {
  public:
    Threaded_execute_param_reduction(rule r, const sequence& s,
                                     // target is the name of the final result
                                     const variable& target,
                                     fact_db& facts, sched_db& scheds)
      :rule_name(r), seq(s)
    {
      exec_rule = new execute_rule(r, s, facts, scheds);
      int tnum = thread_control->num_threads();
      done_signals.resize(tnum);
      for(int i=0;i<tnum;++i) {
        done_signals[i].done = false;
        if(pthread_spin_init(&(done_signals[i].done_lock), 0)) {
          // need better error handling later...
          Loci::Abort();
        }
      }

      final_s = facts.get_variable(target);
      std::vector<sequence> partition = partition_seq(seq, tnum);
      rule_implP ar = rule_name.get_info().rule_impl;
      CPTR<joiner> jop = ar->get_joiner();
      
      for(int i=0;i<tnum;++i) {
        storeRepP sp = final_s->new_store(EMPTY);
        //s->allocate(entitySet(partition[i]));
        sp->allocate(final_s->domain());
        partial.push_back(sp);
      }

      tot_mem = 0;
      for(int i=0;i<tnum;++i) {
        texec.push_back
          (new ExecuteThreaded_param_reduction
           (i, partition[i], //exec_rule,
            new execute_rule(rule_name,partition[i],
                             facts,target,partial[i],scheds),
            jop,thread_control->get_reduction_partners(i),
            done_signals,partial,final_s));
        tot_mem += 16*partition[i].num_intervals();
      }
      tot_mem += 8*seq.num_intervals();

      reduction_root = thread_control->get_reduction_root();
    }

    ~Threaded_execute_param_reduction()
    {
      for(size_t i=0;i<done_signals.size();++i)
        pthread_spin_destroy(&(done_signals[i].done_lock));
      // destroy the local storage
      for(size_t i=0;i<partial.size();++i)
        partial[i]->allocate(EMPTY);
    }
    
    void execute(fact_db& facts, sched_db& scheds)
    {
      current_rule_id = rule_name.ident();
      stopWatch s;
      s.start();

      exec_rule->execute_prelude(seq);
      timer_comp.addTime(s.stop(), seq.size());

      s.start();
      thread_control->setup_facts_scheds(facts,scheds);
      thread_control->restart(texec);
      timer_ctrl.addTime(s.stop(), seq.size());

      thread_control->wait_threads();
      // there are two things to do upon resuming the main thread:
      // 1) reset the done signal of the first work thread so that
      //    it can be reused the next time. (this is not really needed
      //    anyway since no one is going to examine the reduction root's
      //    done flag)
      // 2) store the reduction result into the final place
      s.start();
      done_signals[reduction_root].done = false;
      final_s->copy(partial[reduction_root], final_s->domain());
      // finally execute the postlude
      exec_rule->execute_postlude(seq);
      timer_comp.addTime(s.stop(), seq.size());
      current_rule_id = 0;
    }

    void Print(std::ostream& s) const
    {
      printIndent(s);
      s << rule_name << " over sequence ";
      if(verbose || seq.num_intervals() < 4)
        s << seq;
      else
        s << "{" << seq[0] << "...}";
      s << " using " << texec.size() << " threads";
      if(verbose) {
        s << " (";
        s << "thread " << 0 << ": ";
        texec[0]->print_seq(s);
        for(size_t i=1;i<texec.size();++i) {
          s << ", thread " << i << ": " ;
          texec[i]->print_seq(s);
        }
        s << ")";
      }
      s << std::endl;
    }

    std::string getName() { return "Threaded_execute_param_reduction"; }

    void dataCollate(collectData& data_collector) const
    {
      std::ostringstream oss;
      oss << "rule: " << rule_name;
      // we will need to collect the control time, which is
      // the control time at this level and the maximum control
      // time spent in any thread
      data_collector.accumulateTime(timer_ctrl, EXEC_CONTROL, oss.str());

      double mt = 0;
      for(size_t i=0;i<texec.size();++i) {
        double lt = texec[i]->get_time();
        if (lt > mt)
          mt = lt;
      }
      timeAccumulator new_comp = timer_comp;
      new_comp.addTime(mt, seq.size());
      data_collector.accumulateTime(new_comp,
                                    EXEC_COMPUTATION, oss.str());
      //
      data_collector.accumulateSchedMemory(oss.str(), tot_mem);
    }
  protected:
    timeAccumulator timer_comp;
    timeAccumulator timer_ctrl;
    rule rule_name;
    sequence seq;
    std::vector<ThreadedEMP> texec;
    // here is a vector of partial results for the threads to fill
    // after the thread excution, partial[root] will contain the final result
    std::vector<storeRepP> partial;
    // a vector of reduction signals for the work threads to use
    std::vector<ParamReductionSignal> done_signals;
    // the store in the fact database that the final result will go
    storeRepP final_s;
    // the root thread id of the reduction tree
    int reduction_root;
    double tot_mem;             // tabulate the scheduler memory
    executeP exec_rule;
  };

#ifdef THREAD_CHOMP
  // this is the module for threading chomped rules
  class ExecuteThreaded_chomp: public ThreadedExecutionModule {
  public:
    ExecuteThreaded_chomp(const sequence& s, rule_implP rp)
      :seq(s), exec_size(s.size()), exec_chomp(rp) {}
    void execute(fact_db& facts, sched_db& scheds)
    { 
      stopWatch s;
      s.start();
      exec_chomp->compute(seq); 
      timer.addTime(s.stop(), exec_size);
    }
    void print_seq(std::ostream& s) const
    {
      if(verbose || seq.num_intervals() < 4) {
        s << seq;
      } else {
        s << "{" << seq[0] << "...}";
      }
    }
    double get_time() const { return timer.getTime(); }
    long_long get_l1_dcm() const { 0; }
    long_long get_l2_dcm() const { 0; }
  private:
    // the execution sequence (for printing purpose only)
    sequence seq;
    size_t exec_size;
    // we will just delegate the execution to the underlying
    // concrete rule.
    rule_implP exec_chomp;
    timeAccumulator timer;
  };

  // this is a module that hooks up the compiler and the threaded
  // execution module for threaded pointwise rules defined above
  class Threaded_execute_chomp: public execute_modules {
  public:
    Threaded_execute_chomp
    (const sequence& s,
     const std::vector<rule>& rs,
     const std::vector<std::pair<rule,rule_compilerP> >& comp,
     const std::deque<entitySet>& rseq,
     const variableSet& cv,fact_db& facts,sched_db& scheds)
      :chomp_rules(rs),rule_compilers(comp),rule_seqs(rseq),
       chomp_vars(cv),seq(s)
    {
      seq_size = s.size();
      // we will then need to split the rule sequences.
      // Note: we need to partition the entire collection of
      // sequence (the passed in parameter "s") instead of splitting
      // the individual rule sequences since we want to assign
      // each entity to only one thread
      int threads = thread_control->num_threads();
      std::vector<std::deque<entitySet> > rule_seqs_split(threads);
      std::vector<sequence> partitions = partition_seq(s,threads);
      
      for(int t=0;t<threads;++t) {
        // p is the partition given to thread t
        entitySet p = entitySet(partitions[t]);
        // rst is the rule sequence deque for thread t
        std::deque<entitySet>& rst = rule_seqs_split[t];
        // split the rule sequences according to the partition
        for(size_t i=0;i<rule_seqs.size();++i) {
          entitySet s = rule_seqs[i] & p;
          rst.push_back(s);
        }
      }

      // create threaded execute chomp modules
      for(int t=0;t<threads;++t) {
        const std::deque<entitySet>& trs = rule_seqs_split[t];
        texec.push_back
          (new ExecuteThreaded_chomp
           (partitions[t],
            new execute_chomp
            (entitySet(partitions[t]),rule_compilers,trs,chomp_vars,facts)));
      }
    }
    void execute(fact_db& facts, sched_db& scheds)
    {
      stopWatch s;
      s.start();
      thread_control->setup_facts_scheds(facts,scheds);
      // thread_control->restart(texec);
      // thread_control->wait_threads();
      thread_control->sequential_restart(texec);
      timer.addTime(s.stop(),seq_size);
    }
    void Print(std::ostream& s) const
    {
      printIndent(s);
      s << "--Start chomping over sequence ";
      if(verbose || seq.num_intervals() < 4)
        s<< seq;
      else
        s << "{" << seq[0] << "...}";
      s << " using " << texec.size() << " threads" << std::endl;
      if(verbose) {
        s << " (";
        s << "thread " << 0 << ": ";
        texec[0]->print_seq(s);
        for(size_t i=1;i<texec.size();++i) {
          s << ", thread " << i << ": " ;
          texec[i]->print_seq(s);
        }
        s << ")" << std::endl;
      }
      printIndent(s);
      s << "--Perform chomping on these rules: " << std::endl;
      for(size_t i=0;i<chomp_rules.size();++i) {
        printIndent(s);
        s << "-- " << chomp_rules[i] << std::endl;
      }
      printIndent(s);
      s << "--End chomping" << std::endl;
    }
    std::string getName() { return "Threaded_execute_chomp"; }

    void dataCollate(collectData& data_collector) const
    { // TO BE REVISED...
      std::ostringstream oss ;
      oss <<"chomp: " ;
      for(size_t i=0;i<chomp_rules.size();++i) 
        oss << "[" << chomp_rules[i] << "]";
      data_collector.accumulateTime(timer,EXEC_COMPUTATION,oss.str()) ;
    }
  private:
    timeAccumulator timer_comp;
    timeAccumulator timer_ctrl;
    std::vector<rule> chomp_rules;
    std::vector<std::pair<rule,rule_compilerP> > rule_compilers;
    std::deque<entitySet> rule_seqs;
    variableSet chomp_vars;
    sequence seq;               // union of all rules' sequences
    int seq_size;
    std::vector<ThreadedEMP> texec;
  };
#endif
  
  // threaded execution module for local reduction computation
  // this is the computation part
  class ExecuteThreaded_local_reduction_comp: 
    public ThreadedExecutionModule {
  public:
    ExecuteThreaded_local_reduction_comp
    (int i, const sequence& s, executeP ue, executeP er)
    :tid(i),seq(s),exec_size(s.size()),unit_exec(ue),exec_rule(er) {}
    
    void execute(fact_db& facts, sched_db& scheds)
    {
      stopWatch s;
      s.start();
      // first execute the unit rule to fill the unit values
      unit_exec->execute(facts,scheds);
      // we will then execute the rule
      exec_rule->execute_kernel(seq);
      timer.addTime(s.stop(),exec_size);
    }
    void print_seq(std::ostream& s) const
    {
      if(verbose || seq.num_intervals() < 4) {
        s << seq;
      } else {
        s << "{" << seq[0] << "...}";
      }
    }
    double get_time() const { return timer.getTime(); }
    long_long get_l1_dcm() const { return 0; } // not implemented yet
    long_long get_l2_dcm() const { return 0; }
  private:
    // thread id
    int tid;
    // the execution sequence (note it is a reference)
    const sequence& seq;
    size_t exec_size;
    // this is the unit rule's execution module, we need it 
    // to set the unit value in the private storage
    executeP unit_exec;
    // we will just delegate the main computation to the underlying
    // concrete rule.
    executeP exec_rule;
    timeAccumulator timer;
  };  

  // threaded execution module for local reduction computation
  // this is the reduction part
  class ExecuteThreaded_local_reduction_combine: 
    public ThreadedExecutionModule {
  public:
    ExecuteThreaded_local_reduction_combine
    (int i, const std::vector<sequence>& rdom,
     CPTR<joiner> jop, std::vector<storeRepP>& pr, storeRepP t)
      :tid(i),reduction_dom(rdom),join_op(jop),partial(pr),target(t) {}
    
    void execute(fact_db& facts, sched_db& scheds)
    {
      stopWatch s;
      s.start();
      // we will need to perform the reduction in parallel.
      for(size_t i=0;i<partial.size();++i) {
        // perform the join operation
        join_op->SetArgs(target, partial[i]);
        join_op->Join(reduction_dom[tid]);
      }
      timer.addTime(s.stop(),reduction_dom[tid].size());
    }
    void print_seq(std::ostream& s) const
    { /* no need to print in the combine module... */ }
    double get_time() const { return timer.getTime(); }
    long_long get_l1_dcm() const { return 0; } // not implemented yet
    long_long get_l2_dcm() const { return 0; }
  private:
    // thread id
    int tid;
    // this is sequence used in reduction
    const std::vector<sequence>& reduction_dom;
    // the join operator
    CPTR<joiner> join_op;
    // pointer to the partial results of every thread
    std::vector<storeRepP>& partial;
    // this is the original target of the rule
    storeRepP target;
    timeAccumulator timer;
  };  

  // this is a module that hooks up the compiler and the
  // threaded execution module for param reduction rule defined above
  class Threaded_execute_local_reduction: public execute_modules {
  public:
    Threaded_execute_local_reduction(rule r, rule u, const sequence& s,
                                     // target is the name of 
                                     // the final result
                                     const variable& target,
                                     fact_db& facts, sched_db& scheds)
      :rule_name(r), unit_rule(u), seq(s), tvar(target)
    {
      // creates a real rule execution module
      exec_rule = new execute_rule(r, s, facts, scheds);
      int tnum = thread_control->num_threads();
      // partition the execution sequence among the threads
      thread_seqs = partition_seq(seq, tnum);
      // obtain the rep of the target of the rule
      final_s = facts.get_variable(target);
      // obtain the join operator of the rule
      rule_implP ar = rule_name.get_info().rule_impl;
      CPTR<joiner> jop = ar->get_joiner();
      // figuring out the computed domain of the target
      entitySet context;
      for (sequence::const_iterator si=seq.begin();si!=seq.end();++si)
        context += *si;
      const rule_impl::info& rinfo = rule_name.get_info().desc;
      tvar_exist = vmap_target_exist(*rinfo.targets.begin(),
          facts,context,scheds);
      
      // create all the private storage for work threads
      for(int i=0;i<tnum;++i) {
        storeRepP sp = final_s->new_store(EMPTY);
        partial.push_back(sp);
      }

      tot_mem = seq.num_intervals(); 
      // create the computation and reduction thread module
      for(int i=0;i<tnum;++i) {
        texec_comp.push_back
          (new ExecuteThreaded_local_reduction_comp
           (i, thread_seqs[i],
            new execute_rule(unit_rule,tvar_exist,
                             facts,target,partial[i],scheds),
            new execute_rule(rule_name,thread_seqs[i],
                             facts,target,partial[i],scheds)));
        texec_combine.push_back
          (new ExecuteThreaded_local_reduction_combine
           (i, reduction_dom, jop->clone(), partial, final_s));
      }
    }

    ~Threaded_execute_local_reduction()
    {
      // destroy the local storage
      for(size_t i=0;i<partial.size();++i)
        partial[i]->allocate(EMPTY);
    }
    
    void execute(fact_db& facts, sched_db& scheds)
    {
      current_rule_id = rule_name.ident();
      stopWatch s;
      s.start();
      exec_rule->execute_prelude(seq);
      timer_comp.addTime(s.stop(), seq.size());

      // first start the computation part
      for (size_t i=0;i<partial.size();++i) {
        partial[i]->allocate(final_s->domain());
      }

      s.start();
      thread_control->setup_facts_scheds(facts,scheds);
      thread_control->restart(texec_comp);
      timer_ctrl.addTime(s.stop(), seq.size());

      thread_control->wait_threads();

      // then start the reduction part
      s.start();
      // this step is to generate the parallel reduction sequence
      // for each of the work thread (before starting the threads)
      // we need to do it here because we only get to know the
      // domain of the target at execution time
      if (reduction_dom.empty()) {
        int tnum = thread_control->num_threads();
        reduction_dom = partition_seq(tvar_exist, tnum);
      }

      thread_control->setup_facts_scheds(facts,scheds);
      thread_control->restart(texec_combine);
      timer_ctrl.addTime(s.stop(), seq.size());

      thread_control->wait_threads();

      // destroy the local storage
      for(size_t i=0;i<partial.size();++i)
        partial[i]->allocate(EMPTY);

      s.start();
      // finally execute the postlude
      exec_rule->execute_postlude(seq);
      timer_comp.addTime(s.stop(), seq.size());
      current_rule_id = 0;
    }

    void Print(std::ostream& s) const
    {
      printIndent(s);
      s << rule_name << " over sequence ";
      if(verbose || seq.num_intervals() < 4)
        s << seq;
      else
        s << "{" << seq[0] << "...}";
      s << " using " << texec_comp.size() << " threads";
      if(verbose) {
        s << " (";
        s << "thread " << 0 << ": ";
        texec_comp[0]->print_seq(s);
        for(size_t i=1;i<texec_comp.size();++i) {
          s << ", thread " << i << ": " ;
          texec_comp[i]->print_seq(s);
        }
        s << ")";
      }
      s << std::endl;
    }

    std::string getName() { return "Threaded_execute_local_reduction"; }

    void dataCollate(collectData& data_collector) const
    {
      std::ostringstream oss;
      oss << "rule: " << rule_name;
      // we will need to collect the control time, which is
      // the control time at this level and the maximum control
      // time spent in any thread
      data_collector.accumulateTime(timer_ctrl, EXEC_CONTROL, oss.str());

      double mt = 0;
      for(size_t i=0;i<texec_comp.size();++i) {
        double lt = texec_comp[i]->get_time();
        if (lt > mt)
          mt = lt;
      }
      timeAccumulator new_comp = timer_comp;
      new_comp.addTime(mt, seq.size());

      mt = 0;
      for(size_t i=0;i<texec_combine.size();++i) {
        double lt = texec_combine[i]->get_time();
        if (lt > mt)
          mt = lt;
      }
      new_comp.addTime(mt, seq.size());

      data_collector.accumulateTime(new_comp,
                                    EXEC_COMPUTATION, oss.str());
      //
      data_collector.accumulateSchedMemory(oss.str(), tot_mem);
    }
  protected:
    timeAccumulator timer_comp;
    timeAccumulator timer_ctrl;
    // this is the unit rule associated with this apply rule
    rule rule_name;
    rule unit_rule;
    sequence seq;
    // the computation thread module
    std::vector<ThreadedEMP> texec_comp;
    // the reduction thread module
    std::vector<ThreadedEMP> texec_combine;
    // the partition of the execution sequence for threads
    std::vector<sequence> thread_seqs;
    // here is a vector of partial results for the threads to fill
    std::vector<storeRepP> partial;
    // the parallel reduction domain for each thread
    // this is to be filled at execution time
    std::vector<sequence> reduction_dom;
    // the store in the fact database that the final result will go
    storeRepP final_s; 
    double tot_mem;             // tabulate the scheduler memory
    executeP exec_rule;
    variable tvar;
    // the computed domain of the target variable in this rule
    entitySet tvar_exist;
  };
  
} // end of namespace Loci

#endif  // end #ifdef PTHREADS

#endif
