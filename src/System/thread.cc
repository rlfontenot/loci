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

// disable this file by default for now
#ifdef PTHREADS

#include "thread.h"
#include <sstream>
using std::stringstream;

using std::cout;
using std::cerr;
using std::endl;
using std::vector;
using std::string;
using std::map;
using std::set;
using std::pair;
using std::istream;
using std::ostream;
using std::deque;

namespace Loci {

  ThreadControl* thread_control = 0;

  // this implements the leading execution predicate
  bool is_leading_execution()
  {
    bool lp = MPI_processes>1 ? MPI_rank==0 : true;
    bool lt = thread_control ? thread_control->is_leading_thread() : true;
    return lp && lt;
  }

  // these two functions are intended to provide a global atomic
  // region of execution (the ideal use of them is to use the 
  // Loci preprocessor to hide these calls from the users)
  void global_atomic_region_begin()
  {
    if (thread_control)
      thread_control->atomic_begin();
  }

  void global_atomic_region_end()
  {
    if (thread_control)
      thread_control->atomic_end();
  }
  
  // this function generates a unique name on each process/thread
  std::string gen_local_name(const string& prefix,
                             const string& suffix)
  {
    string mpi_id = "";
    string thread_id = "";
    if (MPI_processes > 1) {
      stringstream ss; ss << MPI_rank;
      mpi_id += ss.str();
    }
    if (thread_control) {
      thread_id += thread_control->get_local_id();
    }
    string name = prefix;
    if (!mpi_id.empty()) {
      name += ".";
      name += mpi_id;
    }
    if (!thread_id.empty()) {
      name += ".";
      name += thread_id;
    }
    if (!suffix.empty()) {
      name += ".";
      name += suffix;
    }
    return name;
  }

  // this function partitions a sequence into n parts
  std::vector<sequence> partition_seq(const sequence& s, int n)
  {
    int psize = s.size() ;
    int div = psize/n ;
    int rem = psize%n ;
    vector<sequence> partition(n) ;

    sequence::const_iterator si = s.begin();
    for(int p=0;p<n;++p) {
      sequence& sp = partition[p];
      for(int k=0;k<div;++k,++si)
        sp += *si;
      if(rem > 0) {
        sp += *si;
        ++si;
        --rem;
      }
    }

    return partition ;    
  }

  ThreadControl_pthread::~ThreadControl_pthread()
  {
    for(int i=0;i<tnum;++i) {
      sem_destroy(&worker_block[i]);
    }
    pthread_attr_destroy(&pattr);
    sem_destroy(&main_block);
    for(size_t i=0;i<finish_notification.size();++i) {
      pthread_spin_destroy(&(finish_notification[i].done_lock));
    }
  }

  // these two functions build the termination tree
  
  void ThreadControl_pthread::
  build_term_tree(int pid, int depth, size_t& myid)
  {
    if(myid >= finish_notification.size())
      return;

    size_t me = myid;
    finish_notification[me].done = false;
    if(pthread_spin_init(&(finish_notification[me].done_lock), 0)) {
      // need better error handling later...
      Loci::Abort();
    }
    finish_notification[me].p = pid;
    if(pid >= 0) {
      finish_notification[pid].c.push_back(me);
    }
    ++myid;

    if(depth>0) {
      build_term_tree(me, depth-1, myid);
      build_term_tree(me, depth-1, myid);
    }
  }

  // a function that computes the nearest power of 2 for a given integer
  int pow2_roundup(int n)
  {
    int start = 1;
    while(start < n)
      start <<= 1;
    return start;
  }

  // build a list of partners for the parallel starting process
  void build_partner(int me, int limit, int stride,
                     vector<vector<int> >& partner)
  {
    stride >>= 1;
    if(stride < 1)
      return;
    int p = me + stride;
    if(p < limit)
      partner[me].push_back(p);
    build_partner(me, limit, stride, partner);
    if(p < limit)
      build_partner(p, limit, stride, partner);
  }
  
  ThreadControl_pthread::ThreadControl_pthread(int n)
    :active(false),stop(false),work(0),factsP(0),schedsP(0),tnum(n)
  {
    thread_args.resize(n);
    worker_block.resize(n);
    if(sem_init(&main_block, 0, 0)) {
      // need to better handle errors...
      Loci::Abort();
    }

    // build the parallel start partner
    int stride = pow2_roundup(tnum);
    start_partners.resize(tnum);
    build_partner(0, tnum, stride, start_partners);

    // build the termination notification tree
    int size = tnum;
    finish_notification.resize(size);
    int depth = 0;
    while(size > 1) {
      ++depth;
      size /= 2;
    }
    size_t init_id = 0;
    build_term_tree(-1, depth, init_id);

    // build the reduction partner (just the reverse of the start_partners)
    reduction_partners.resize(tnum);
    for(int i=0;i<tnum;++i) {
      vector<int>::const_reverse_iterator b,e;
      b = start_partners[i].rbegin();
      e = start_partners[i].rend();
      reduction_partners[i] = vector<int>(b,e);
    }

    // system contention scope
    pthread_attr_init(&pattr);
    pthread_attr_setscope(&pattr, PTHREAD_SCOPE_SYSTEM);
    for(int i=0;i<n;++i) {
      if(sem_init(&worker_block[i], 0, 0)) {
        // need to better handle errors...
        Loci::Abort();
      }
    }

    stop = false;
  }
  
  void ThreadControl_pthread::create_threads()
  {
    // first record the main control thread id
    main_pid = pthread_self();
    id_map[main_pid] = "tmain";
    // creating threads
    for(int i=0;i<tnum;++i) {
      thread_args[i].tid = i;
      thread_args[i].self = this;
      int rc = pthread_create(&(thread_args[i].ptid), &pattr, thread_fun,
                              (void*)(&thread_args[i]));
      if(rc) {
        // need to handle errors better later...
        Loci::Abort();
      }
      string id = "t";
      stringstream ss; ss << i;
      id += ss.str();
      id_map[thread_args[i].ptid] = id;
    }
    // record the leading work thread id
    lead_work_pid = thread_args[0].ptid;
  }

  string ThreadControl_pthread::get_local_id() const
  {
    pthread_t self = pthread_self();
    map<pthread_t,string>::const_iterator mi = id_map.find(self);
    if (mi != id_map.end())
      return mi->second;
    else
      return "tUNKNOWN";
  }

  // determine if the calling thread is leading thread or not.
  // the leading thread in the pthread implementation is defined
  // to be either the main control thread or the leading work 
  // thread.  This is valid since in our design, the main control
  // thread will never work together with any of the work threads.
  bool ThreadControl_pthread::is_leading_thread() const
  {
    pthread_t sid = pthread_self();
    return (pthread_equal(sid,main_pid) || 
            pthread_equal(sid,lead_work_pid));
  }

  const vector<int>& ThreadControl_pthread::
  get_reduction_partners(int id) const
  {
    if(id < 0 || id >= tnum)
      throw ThreadException("get_reduction_partners: thread id not valid");
    return reduction_partners[id];
  }

  // this is the thread function to execute on each thread upon creation.
  // basically it just waits until there is work to do.
  void* ThreadControl_pthread::thread_fun(void* arg)
  {
    // retrieve the argument
    TArg* targ = (TArg*)arg;
    int tid = targ->tid;
    ThreadControl_pthread* self = targ->self;

    // // we want to set the CPU affinity of threads so that
    // // each thread is attached to a specific CPU and wont migrate (hopefully)
    // cpu_set_t cpuset;
    
    // CPU_ZERO(&cpuset);
    // CPU_SET(Loci::MPI_rank*self->tnum + tid, &cpuset);
    
    // if( (pthread_setaffinity_np(pthread_self(),
    //                             sizeof(cpu_set_t), &cpuset)) != 0)
    //   cerr << "WARNING: pthread CPU affinity setting failed!" << endl;

    // this holds the pointer to the next work module
    ThreadedEMP w = 0;

    while(true) {
      // first block on the semaphore waiting for the
      // go signal from the main thread
      sem_wait(&(self->worker_block[tid]));

      // parallel start...
      for(size_t i=0;i<(self->start_partners)[tid].size();++i) {
        int p = (self->start_partners)[tid][i];
        sem_post(&((self->worker_block)[p]));
      }

      if(self->stop)
        break;

      // do the module work
      w = (*(self->work))[tid];
      fact_db* f = self->factsP;
      sched_db* s = self->schedsP;
      if(w != 0 && f && s) {
        w->execute(*f, *s);
      }

      // propagate the done flag so that eventually the main thread
      // is notified that all are finished.

      // first of all, we will check are all children in the tree done yet?
      vector<Terminate>& fn = self->finish_notification;
      for(size_t i=0;i<fn[tid].c.size();++i) {
        int cid = fn[tid].c[i];
        bool done = false;
        while(!done) {
          pthread_spin_lock(&(fn[cid].done_lock));
          done = fn[cid].done;
          pthread_spin_unlock(&(fn[cid].done_lock));
        }
        // we will then reset the done flag in the child
        // we don't need to acquire the lock at this time since
        // this child should have finished by now.
        fn[cid].done = false;
      }
      // then we need to write the done flag
      // the root thread needn't do this since it just needs to
      // notify the main thread
      if(fn[tid].p < 0) {
        sem_post(&self->main_block);
      } else {
        pthread_spin_lock(&(fn[tid].done_lock));
        fn[tid].done = true;
        pthread_spin_unlock(&(fn[tid].done_lock));
      }
    }
    return 0;
  }

  void ThreadControl_pthread::shutdown()
  {
    // only the main thread is able to call this function
    // also, if the main thread is active, then all the
    // work threads must be blocked, so it is safe to update
    // any data at this moment.
    stop = true;
    vector<ThreadedEMP> w(tnum, ThreadedEMP(0));
    restart(w);
    // joining all the threads
    for(int i=0;i<tnum;++i) {
      if(pthread_join(thread_args[i].ptid,NULL)) {
        // better error handling later...
        Loci::Abort();
      }
    }
    active = false;
  }
  
  void ThreadControl_pthread::
  restart(std::vector<ThreadedEMP>& w)
  {
    if( (int)w.size() != tnum) {
      throw ThreadException("restart work size wrong!");
    }
    work = &w;
    active = true;
    sem_post(&worker_block[0]);
  }
  
  void ThreadControl_pthread::wait_threads()
  {
    if(!active)
      return;

    sem_wait(&main_block);
    active = false;
  }
  
  void ThreadControl_pthread::
  sequential_restart(std::vector<ThreadedEMP>& w)
  {
    if( (int)w.size() != tnum) {
      throw ThreadException("sequential_restart work size wrong!");
    }
    std::vector<ThreadedEMP> nw(tnum);
    for(int t=0;t<tnum;++t) {
      for(int i=0;i<tnum;++i)
        nw[i] = ThreadedEMP(0);
      nw[t] = w[t];
      restart(nw);
      wait_threads();
    }
  }

  // the implementation of the UDG graph type
  // and the graph coloring algorithm
  void UDG::add_edge(int x, int y)
  {
    // since it is undirected graph, we need to add y to the record of x,
    // and also add x to the associated set of y

    // if x and/or y do not exist, these will automatically create them.
    entitySet& xn = g[x];
    entitySet& yn = g[y];
    xn += y;
    yn += x;
    vertices += x;
    vertices += y;
  }

  void UDG::add_edges(int x, const entitySet& ys)
  {
    entitySet& xn = g[x];
    xn += ys;
    vertices += x;
    for(entitySet::const_iterator ei=ys.begin();ei!=ys.end();++ei) {
      entitySet& yn = g[*ei];
      yn += x;
      vertices += *ei;
    }
  }

  void UDG::add_edges(const entitySet& xs)
  {
    for(entitySet::const_iterator ei=xs.begin();ei!=xs.end();++ei) {
      entitySet others = xs; others -= *ei;
      add_edges(*ei, others);
    }
  }
  
  void UDG::delete_edge(int x, int y)
  {
    map<int,entitySet>::iterator mi = g.find(x);
    if(mi != g.end()) {
      entitySet& xn = mi->second;
      xn -= y;
    }
    mi = g.find(y);
    if(mi != g.end()) {
      entitySet& yn = mi->second;
      yn -= x;
    }
  }

  void UDG::delete_edges(int x, const entitySet& ys)
  {
    map<int,entitySet>::iterator mi = g.find(x);
    if(mi != g.end()) {
      entitySet& xn = mi->second;
      xn -= ys;
    }
    for(entitySet::const_iterator ei=ys.begin();ei!=ys.end();++ei) {
      mi = g.find(*ei);
      if(mi != g.end()) {
        entitySet& yn = mi->second;
        yn -= x;
      }
    }
  }

  void UDG::add_vertex(int x)
  {
    vertices += x;
    g[x];
  }

  void UDG::add_vertices(const entitySet& xs)
  {
    vertices += xs;
    for(entitySet::const_iterator ei=xs.begin();ei!=xs.end();++ei)
      g[*ei];
  }

  void UDG::delete_vertex(int x)
  {
    map<int,entitySet>::iterator mi = g.find(x);
    if(mi != g.end()) {
      entitySet& xn = mi->second;
      for(entitySet::const_iterator ei=xn.begin();ei!=xn.end();++ei) {
        map<int,entitySet>::iterator mi2 = g.find(*ei);
        if(mi2 != g.end()) {
          entitySet& vn = mi2->second;
          vn -= x;
        }
      }
      g.erase(mi);
      vertices -= x;
    }
  }

  void UDG::delete_vertices(const entitySet& xs)
  {
    for(entitySet::const_iterator ei=xs.begin();ei!=xs.end();++ei)
      delete_vertex(*ei);
  }

  void UDG::contract_vertex(int x, int y)
  {
    // remove y and make all neighbors of y neighbors of x
    map<int,entitySet>::iterator mi = g.find(y);
    if(mi != g.end()) {
      entitySet yn = mi->second;
      for(entitySet::const_iterator ei=yn.begin();ei!=yn.end();++ei) {
        map<int,entitySet>::iterator mi2 = g.find(*ei);
        if(mi2 != g.end()) {
          entitySet& vn = mi2->second;
          vn -= y;
        }
      }
      g.erase(mi);
      vertices -= y;
      // we want to remove (possible) x from cyn so that when we
      // connect these vertices to x, we don't have a self-loop.
      yn -= x;
      add_edges(x, yn);
    }
  }

  const entitySet& UDG::neighbors(int x) const
  {
    map<int,entitySet>::const_iterator mi = g.find(x);
    if(mi == g.end())
      throw UDGRangeError();

    return mi->second;
  }

  entitySet UDG::non_neighbors(int x) const
  {
    map<int,entitySet>::const_iterator mi = g.find(x);
    if(mi == g.end())
      throw UDGRangeError();

    const entitySet& xn = mi->second;

    entitySet nn = vertices - xn;
    nn -= x;

    return nn;
  }

  vector<pair<int,int> > UDG::all_edges() const
  {
    vector<pair<int,int> > edges;
    entitySet seen;
    for(map<int,entitySet>::const_iterator mi=g.begin();mi!=g.end();++mi) {
      int x = mi->first;
      const entitySet& n = mi->second;
      seen += x;
      entitySet o = n - seen;
      for(entitySet::const_iterator ei=o.begin();ei!=o.end();++ei) {
        edges.push_back(pair<int,int>(x, *ei));
      }
    }
    return edges;
  }

  entitySet UDG::dangling_vertices() const
  {
    entitySet dv;
    for(map<int,entitySet>::const_iterator mi=g.begin();mi!=g.end();++mi) {
      if(mi->second == EMPTY)
        dv += mi->first;
    }
    return dv;
  }
  
  const entitySet& UDG::max_degree(int& x) const
  {
    if(vertices == EMPTY)
      throw UDGRangeError();

    map<int,entitySet>::const_iterator mx, mi = g.begin();
    mx = mi;
    int md = mx->second.size();
  
    for(++mi;mi!=g.end();++mi) {
      int d = mi->second.size();
      if(d > md) {
        md = d;
        mx = mi;
      }
    }
    x = mx->first;
    return mx->second;
  }

  const entitySet& UDG::max_degree(const entitySet& xs, int& x) const
  {
    if(vertices == EMPTY)
      throw UDGRangeError();

    map<int,entitySet>::const_iterator mx, mi = g.begin();
    entitySet::const_iterator ei = xs.begin();
    int md = -1;
    while(mi!=g.end() && ei!=xs.end()) {
      if(mi->first == *ei) {
        int d = mi->second.size();
        if(d > md) {
          md = d;
          mx = mi;
        }
        ++mi; ++ei;
      } else if(mi->first < *ei) {
        ++mi;
      } else {
        ++ei;
      }
    }

    if(md == -1)
      throw UDGRangeError();

    x = mx->first;
    return mx->second;
  }

  entitySet UDG::common_neighbors(int x, int y) const
  {
    map<int,entitySet>::const_iterator mi = g.find(x);
    if(mi == g.end())
      throw UDGRangeError();
    const entitySet& xn = mi->second;

    mi = g.find(y);
    if(mi == g.end())
      throw UDGRangeError();
    const entitySet& yn = mi->second;

    return (xn & yn);
  }

  ostream& UDG::show_stats(std::ostream& s) const
  {
    s << "total number of vertices: " << vertices.size() << endl;
    int en = 0;
    int max_d = 0;
    int min_d = vertices.size() + 1;
    for(map<int,entitySet>::const_iterator mi=g.begin();mi!=g.end();++mi) {
      int es = mi->second.size();
      en += es;
      if(es > max_d)
        max_d = es;
      if(es < min_d)
        min_d = es;
    }
    en /= 2;
    s << "total number of edges: " << en << endl;
    s << "maximum degree: " << max_d << endl;
    s << "minimum degree: " << min_d << endl;
    s << "total number of dangling vertices: "
      << dangling_vertices().size() << endl;

    return s;
  }

  istream& remove_space(istream& s)
  {
    int c = s.peek();
    while(isspace(c)) {
      s.get();
      c = s.peek();
    }
    return s;
  }

  istream& operator>>(istream& s, UDG& g)
  {
    while(s) {
      remove_space(s);
      char p = s.peek();
      if(p != '(')
        return s;

      s.get();

      int x,y;
      s >> x;
      if(s.eof() || x<0)
        throw UDGFormatError();

      // peek to see if there is any edges attached
      p = s.peek();
      if(p == ')') {
        s.get();
        g.add_vertex(x);
        continue;
      }
      
      s >> p;
      if(s.eof() || p!=',')
        throw UDGFormatError();
    
      s >> y;
      if(s.eof() || y<0)
        throw UDGFormatError();
    
      s >> p;
      if(s.eof())
        throw UDGFormatError();

      if(p == ')')
        g.add_edge(x,y);
      else if(p == '|') {
        entitySet vs;
        vs += y;
        while(s) {
          s >> y;
          if(s.eof() || y<0)
            throw UDGFormatError();
          vs += y;
          s >> p;
          if(p == ')')
            break;
          if(p != '|')
            throw UDGFormatError();
        }
        g.add_edges(x,vs);
        // consecutive vertices list
      } else
        throw UDGFormatError();
    }
    return s;
  }

  ostream& operator<<(ostream& s, const UDG& g)
  {
    for(map<int,entitySet>::const_iterator mi=g.g.begin();mi!=g.g.end();++mi) {
      if(mi->second == EMPTY) {
        s << "(" << mi->first << ")" << endl;
      } else if(mi->second.size() == 1) {
        s << "(" << mi->first << "," << *(mi->second.begin()) << ")" << endl;
      } else {
        entitySet::const_iterator ei = mi->second.begin();
        s << "(" << mi->first << "," << *ei;
        for(++ei;ei!=mi->second.end();++ei)
          s << "|" << *ei;
        s << ")" << endl;
      }
    }
    return s;
  }

  map<int,int> rlf_color(UDG& g)
  {
    map<int,int> color;
    int c = 0;
    while(g.all_vertices() != EMPTY) {
      int x;
      g.max_degree(x);
      color[x] = c;
      entitySet nn = g.non_neighbors(x);
      while(nn != EMPTY) {
        int maxcn = -1;
        int yd = -1;
        int y;
        for(entitySet::const_iterator ni=nn.begin();ni!=nn.end();++ni) {
          int z = *ni;
          entitySet cn = g.common_neighbors(x, z);
          int cns = cn.size();
          int zd = g.neighbors(z).size();
          if(cns>maxcn || (cns==maxcn && zd<yd)) {
            y = z;
            yd = zd;
            maxcn = cns;
          }
        }
        if(maxcn == 0)
          g.max_degree(nn, y);
        color[y] = c;
        g.contract_vertex(x, y);
        nn = g.non_neighbors(x);
      }
      g.delete_vertex(x);
      ++c;
    }
    return color;
  }

  // small type used to sort the graph vertices
  struct VDT {
    int v;
    int d;
    VDT(int vv, int dd):v(vv),d(dd) {}
  };
  bool operator<(const VDT& vd1, const VDT& vd2)
  {
    return vd1.d < vd2.d;
  }

  map<int,int> lf_color(const UDG& g)
  {
    // first we will need to sort all the vertices according to
    // their degrees in a non-increasing order
    vector<VDT> vertices;
    for(map<int,entitySet>::const_iterator mi=g.g.begin();mi!=g.g.end();++mi) {
      vertices.push_back(VDT(mi->first,mi->second.size()));
    }
    sort(vertices.begin(), vertices.end());
    // then we will start color
    map<int,int> color;
    for(vector<VDT>::const_iterator
          vi=vertices.begin();vi!=vertices.end();++vi) {
      // get the neighbors;
      int v = vi->v;
      const entitySet& nn = g.neighbors(v);
      // find out the smallest color not yet used in the neighborhood
      vector<int> cc;
      map<int,int>::const_iterator mi;
      for(entitySet::const_iterator ei=nn.begin();ei!=nn.end();++ei) {
        mi = color.find(*ei);
        if(mi != color.end())
          cc.push_back(mi->second);
      }
      sort(cc.begin(), cc.end());
      // search cc to find out the smallest color not yet used
      int c = 0;
      for(vector<int>::const_iterator cci=cc.begin();cci!=cc.end();++cci) {
        if(c < *cci)
          break;
        else
          c = *cci + 1;
      }
      color[v] = c;
    }
    return color;
  }

  bool verify_color(const map<int,int>& color, const UDG& g)
  {
    for(map<int,int>::const_iterator mi=color.begin();mi!=color.end();++mi) {
      int v = mi->first;
      int c = mi->second;
      const entitySet& nn = g.neighbors(v);
      for(entitySet::const_iterator ei=nn.begin();ei!=nn.end();++ei) {
        map<int,int>::const_iterator ci=color.find(*ei);
        if(ci != color.end()) {
          if(ci->second == c)
            return false;
        } else
          return false;
      }
    }
    return true;
  }
  

  // this function generates an interference graph for
  // a local reduction rule.
  void build_interference_graph(const rule& r,
                                fact_db& facts, sched_db& scheds,
                                const entitySet& context, UDG& g)
  {
    // DEBUG
    // double s = 0, t1 = 0, t2 = 0, t3 = 0;
    ////////

    // DEBUG
    // s = MPI_Wtime();
    ////////
    
    // first we need to get the mapping info in the targets of the rule
    const rule_impl::info& rinfo = r.get_info().desc;
    // build a mapping from each entity in the context to all the
    // target entity sets, we assume that there is just one target
    set<vmap_info>::const_iterator vmi = rinfo.targets.begin();
    map<Entity,entitySet> s2t;
    for(entitySet::const_iterator ei=context.begin();ei!=context.end();++ei) {
      entitySet input; input += *ei;
      s2t[*ei] = vmap_target_exist(*vmi, facts, input, scheds);
    }

    // DEBUG
    // t1 = MPI_Wtime() - s;
    // s = MPI_Wtime();
    ////////
    
    // we then build an inverse map of "s2t", that is, for each entity (e)
    // in the target's domain, we will get what the input entities are
    // responsible to write to "e"
    // the map "s2t" tells us for each of the input entity, which portion of
    // the target is being written, so the inverse of "s2t" tells us
    // the set of input entities that are writing to the same location,
    // i.e., the interference map.
    map<Entity,entitySet> t2s;
    for(map<Entity,entitySet>::const_iterator
          mi=s2t.begin();mi!=s2t.end();++mi) {
      Entity ss = mi->first;
      const entitySet& tt = mi->second;
      for(entitySet::const_iterator ei=tt.begin();ei!=tt.end();++ei) {
        t2s[*ei] += ss;
      }
    }

    // DEBUG
    // t2 = MPI_Wtime() - s;
    // s = MPI_Wtime();
    ////////
    
    // we are now ready to generate the interference graph.
    // all the clusters in each entry in the mapping "t2s" are connected
    // in the final graph.

    // DEBUG
    // int max_s = 0;
    //
    
    for(map<Entity,entitySet>::const_iterator
          mi=t2s.begin();mi!=t2s.end();++mi) {
      g.add_edges(mi->second);
      // DEBUG
      //if(mi->second.size() > max_s)
      //  max_s = mi->second.size();
      //
    }

    // DEBUG
    // t3 = MPI_Wtime() - s;
    // cout << "max partial reduction domain size: " << max_s << endl;
    // cout << "t1 = " << t1 << endl;
    // cout << "t2 = " << t2 << endl;
    // cout << "t3 = " << t3 << endl;
    // 
  }

  // partition entity set based on colors
  vector<entitySet> partition_color(const map<int,int>& color)
  {
    map<int,entitySet> part_map;
    for(map<int,int>::const_iterator mi=color.begin();mi!=color.end();++mi)
      part_map[mi->second] += mi->first;
    vector<entitySet> part;
    for(map<int,entitySet>::const_iterator
          mi=part_map.begin();mi!=part_map.end();++mi)
      part.push_back(mi->second);

    return part;
  }

  //pthread_mutex_t copy_mutex = PTHREAD_MUTEX_INITIALIZER; 
  
  
} // end of namespace Loci

#else
#include "distribute.h"
namespace Loci {
  // this implements the leading execution predicate
  bool is_leading_execution()
  {
    return MPI_rank==0 ;
  }

  // these two functions are intended to provide a global atomic
  // region of execution (the ideal use of them is to use the 
  // Loci preprocessor to hide these calls from the users)
  void global_atomic_region_begin()
  {
  }

  void global_atomic_region_end()
  {
  }
}
#endif
