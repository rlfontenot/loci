#ifndef RULE_H
#define RULE_H

#include <Tools/debug.h>
#include <entitySet.h>
#include <store_rep.h>

#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <set>

#include <Tools/cptr.h>
#include <Tools/expr.h>
#include <variable.h>
#include <Map.h>
#include <parameter.h>
#include <storeVec.h>

namespace Loci {

  class fact_db ;
  
  class joiner : public CPTR_type {
  public:
    virtual CPTR<joiner> clone() = 0 ;
    virtual storeRepP getTargetRep() = 0 ;
    virtual void SetArgs(storeRepP &target, storeRepP &source) = 0 ;
    virtual void Join(const sequence &seq) = 0 ;
    virtual void Join(Map &t2s, const sequence &seq) = 0 ;
  } ;  

  class rule_impl : public CPTR_type {
  public:
    struct info {
      std::set<vmap_info> sources, targets, constraints ;
      variableSet conditionals ;
      std::string rule_identifier() const ;

      variableSet input_vars() const ;
      variableSet output_vars() const ;
    } ;
    typedef CPTR<rule_impl> rule_implP ;
    enum rule_impl_type {POINTWISE,SINGLETON,UNIT,APPLY,JOIN, UNKNOWN} ;
  private:
    rule_impl_type rule_impl_class ;
    bool rule_threading ;
    mutable std::string name ;

    info rule_info ;
        
    typedef std::map<variable, store_instance *> storeIMap ;
    storeIMap var_table ;
    void source(const std::string &invar) ;
    void target(const std::string &outvar) ;

  protected:
    rule_impl(rule_impl &f) { fatal(true) ; }
    void rule_class(rule_impl_type ft) { rule_impl_class = ft ; }
    void disable_threading() { rule_threading = false ; }
    void rule_name(const std::string &name) ;
    void name_store(const std::string &name,store_instance &si) ;
    void input(const std::string &invar) { source(invar) ; }
    void output(const std::string &outvar) { target(outvar) ; }
    void constraint(const std::string &constrain) ;
    void conditional(const std::string &cond) ;

    void base_copy(const rule_impl &f) {
      rule_impl_class = f.rule_impl_class ;
      rule_threading = f.rule_threading ;
      name = f.get_name() ;
      rule_info = f.rule_info ;
      var_table = f.var_table ;
    }

  public:
    rule_impl() ;
    bool check_perm_bits() const ;
    bool thread_rule() const { return rule_threading; }
    void initialize(fact_db &facts) ;
    std::string get_name() const ;
    rule_impl_type get_rule_class() const { return rule_impl_class ; }
    const info &get_info() const { return rule_info ; }
    void set_store(variable v, const storeRepP &p) ;
    void set_store(const std::string &name, const storeRepP &p) 
      { set_store(variable(expression::create(name)),p) ; }
    storeRepP get_store(variable v) const ;
    storeRepP get_store(const std::string &name) const
      { return get_store(variable(expression::create(name))) ; }

    void set_variable_times(time_ident tl) ;
    void copy_store_from(rule_impl &f) ;
    void Print(std::ostream &s) const ;

    virtual rule_implP new_rule_impl() const ;
    virtual void compute(const sequence &) = 0 ;
    virtual CPTR<joiner> get_joiner() = 0 ;
  } ;

  typedef rule_impl::rule_implP rule_implP ;

  template <class T> class copy_rule_impl : public rule_impl {
    rule_implP realrule_impl ;
  public:
    copy_rule_impl() {realrule_impl = new T ; base_copy(*realrule_impl) ;}
    virtual rule_implP new_rule_impl() const { return new copy_rule_impl<T> ; }
    virtual void compute(const sequence &s)  { realrule_impl->compute(s) ; }
    virtual CPTR<joiner> get_joiner() { return realrule_impl->get_joiner() ; }
  } ;
        
  class pointwise_rule : public rule_impl {
   protected:
    pointwise_rule() { rule_class(POINTWISE) ; }
    void name_store(const std::string &name, store_instance &si)
    { rule_impl::name_store(name,si) ; }
    void input(const std::string &invar)
    { rule_impl::input(invar) ; }
    void output(const std::string &outvar)
    { rule_impl::output(outvar) ; }
    void constraint(const std::string &constrain)
    { rule_impl::constraint(constrain) ; }
    void conditional(const std::string &cond)
    { rule_impl::conditional(cond) ; }
    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0) ; }
  } ;

  class singleton_rule : public rule_impl {
  protected:
    singleton_rule() { rule_class(SINGLETON) ; }
    void name_store(const std::string &name, store_instance &si)
    { rule_impl::name_store(name,si) ; }
    void input(const std::string &invar)
    { rule_impl::input(invar) ; }
    void output(const std::string &outvar)
    { rule_impl::output(outvar) ; }
    void constraint(const std::string &constrain)
    { rule_impl::constraint(constrain) ; }
    void conditional(const std::string &cond)
    { rule_impl::conditional(cond) ; }
    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0) ; }
  } ;
  
  class unit_rule : public rule_impl {
   protected:
    unit_rule() { rule_class(UNIT) ; }
    void name_store(const std::string &name, store_instance &si)
    { rule_impl::name_store(name,si) ; }
    void input(const std::string &invar)
    { rule_impl::input(invar) ; }
    void output(const std::string &outvar)
    { rule_impl::output(outvar) ; }
    void constraint(const std::string &constrain)
    { rule_impl::constraint(constrain) ; }
    void conditional(const std::string &cond)
    { rule_impl::conditional(cond) ; }
    virtual CPTR<joiner> get_joiner() { return CPTR<joiner>(0) ; }
  } ;

  template <class T, class Op> class joinOp : public joiner {
    Op join ;
    T t,s ;
  public:
    virtual CPTR<joiner> clone()
    { return CPTR<joiner>(new joinOp<T,Op> ); }

    virtual storeRepP getTargetRep()
    { T st ; storeRepP rep = st.Rep(); return rep; }

    virtual void SetArgs(storeRepP &target, storeRepP &source)
    { s.setRep(source) ; t.setRep(target) ; }

    virtual void Join(const sequence &seq) {
      for(sequence::const_iterator i=seq.begin();i!=seq.end();++i) 
        join(t[*i],s[*i]) ;
    }

    virtual void Join(Map &t2s, const sequence &seq) { 
      for(sequence::const_iterator i=seq.begin();i!=seq.end();++i) 
        join(t[*i],s[t2s[*i]]) ;
    }
  } ;

  template<class Type,class Op> class joinOp<param<Type>,Op> : public joiner {
    Op join ;
    param<Type> s,t ;
  public:
    virtual CPTR<joiner> clone()
    { return CPTR<joiner>(new joinOp<param<Type>,Op> ); }

    virtual void SetArgs(storeRepP &target, storeRepP &source)
    { s.setRep(source) ; t.setRep(target) ; }

    virtual storeRepP getTargetRep()
    { param<Type> st ; storeRepP rep = st.Rep(); return rep; }

    virtual void Join(const sequence &seq)
    {join(*t,*s) ;}

    virtual void Join(Map &t2s, const sequence &seq)
    {join(*t,*s) ;}
  } ;

  template<class Type,class Op> class joinOp<storeVec<Type>,Op> : public joiner {
    Op join ;
    storeVec<Type> s,t ;
  public:
    virtual CPTR<joiner> clone()
    { return CPTR<joiner>(new joinOp<storeVec<Type>,Op> ); }

    virtual void SetArgs(storeRepP &target, storeRepP &source)
    { s.setRep(source) ; t.setRep(target) ; }

    virtual storeRepP getTargetRep()
    { param<Type> st ; storeRepP rep = st.Rep(); return rep; }

    virtual void Join(const sequence &seq) {
      for(sequence::const_iterator i=seq.begin();i!=seq.end();++i) {
        Vect<Type> m = t[*i] ;
        join(m,s[*i]) ;
      }
    }

    virtual void Join(Map &t2s, const sequence &seq) {
      for(sequence::const_iterator i=seq.begin();i!=seq.end();++i) {
        Vect<Type> m = t[*i] ;
        join(m,s[t2s[*i]]) ;
      }
    }
  } ;

  template<class Type,class Op> class joinOp<storeMat<Type>,Op> : public joiner {
    Op join ;
    storeMat<Type> s,t ;
  public:
    virtual CPTR<joiner> clone()
    { return CPTR<joiner>(new joinOp<storeMat<Type>,Op> ); }

    virtual void SetArgs(storeRepP &target, storeRepP &source)
    { s.setRep(source) ; t.setRep(target) ; }

    virtual storeRepP getTargetRep()
    { param<Type> st ; storeRepP rep = st.Rep(); return rep; }

    virtual void Join(const sequence &seq) {
      for(sequence::const_iterator i=seq.begin();i!=seq.end();++i) {
        Mat<Type> m = t[*i] ;
        join(m,s[*i]) ;
      }
    }

    virtual void Join(Map &t2s, const sequence &seq) {
      for(sequence::const_iterator i=seq.begin();i!=seq.end();++i) {
        Mat<Type> m = t[*i] ;
        join(m,s[t2s[*i]]) ;
      }
    }
  } ;

  template<class Type, class Op> class joinOp<multiStore<Type>,Op> :
  public joiner {
    Op join ;
    multiStore<Type> s,t ;
  public:
    virtual CPTR<joiner> clone()
    { return CPTR<joiner>(new joinOp<multiStore<Type>,Op> ); }

    virtual void SetArgs(storeRepP &target, storeRepP &source)
    { s.setRep(source) ; t.setRep(target) ; }

    virtual storeRepP getTargetRep()
    { multiStore<Type> st ; storeRepP rep = st.Rep(); return rep; }

    virtual void Join(const sequence &seq) {
      for(sequence::const_iterator i=seq.begin();i!=seq.end();++i) {
        int szt = t.end(*i)-t.begin(*i) ;
        fatal(szt != (s.end(*i)-s.begin(*i))) ;
        for(int j=0;j<szt;++j)
          join(t[*i][j],s[*i][j]) ;
      }
    }

    virtual void Join(Map &t2s, const sequence &seq) {
      for(sequence::const_iterator i=seq.begin();i!=seq.end();++i) {
        int szt = t.end(*i)-t.begin(*i) ;
        fatal(szt != (s.end(t2s[*i])-s.begin(t2s[*i]))) ;
        for(int j=0;j<szt;++j)
          join(t[*i][j],s[t2s[*i]][j]) ;
      }
    }
  } ;
    
  template <class T, class Op > class apply_rule : public rule_impl {
  protected:
    apply_rule() { rule_class(APPLY) ; }
    void name_store(const std::string &name, store_instance &si)
    { rule_impl::name_store(name,si) ; }
    void input(const std::string &invar)
    { rule_impl::input(invar) ; }
    void output(const std::string &outvar)
    { rule_impl::output(outvar) ; }
    void constraint(const std::string &constrain)
    { rule_impl::constraint(constrain) ; }
    void conditional(const std::string &cond)
    { rule_impl::conditional(cond) ; }
    void join(typename T::containerType &t1,
              const typename T::containerType &t2) {
      Op f ;
      f(t1,t2) ;
    }
  public:
    virtual CPTR<joiner> get_joiner() {
      CPTR<joinOp<T,Op> > jp = new joinOp<T,Op> ;
      return CPTR<joiner>(jp) ;
    }

  } ;
  template <class T> struct Summation {
    void operator()(T &res, const T &arg)
    { res += arg ; }
  } ;

  template <class T> struct Product {
    void operator()(T &res, const T &arg)
    { res *= arg ; }
  } ;

  template <class T> struct Maximum {
    void operator()(T &res ,const T &arg)
    { res = max(res,arg) ; }
  } ;

  template <class T> struct Minimum {
    void operator()(T &res, const T &arg)
    { res = min(res,arg) ; }
  } ;
  
  class rule {
  public:
    enum rule_type {BUILD=0,COLLAPSE=1,GENERIC=2,TIME_SPECIFIC=3,INTERNAL=4} ;
    struct info {
      rule_implP rule_impl ;
        
      rule_impl::info desc ;
        
      std::string rule_ident ;
        
      time_ident source_level, target_level ;
      variableSet source_vars, target_vars, map_vars, constraint_vars ;
        
      rule_type rule_class ;
      bool output_is_parameter ;
      bool time_advance ;
      string internal_qualifier ;
        
      const std::string &name() const { return rule_ident ; }
      info() { rule_ident = "NO_RULE" ;}
      info(const rule_implP &fp) ;
      info(const info &fi, time_ident tl) ;
      info(const std::string &s) ;

      const info &get_info() const { return *this ; }
      const variableSet &sources() const { return source_vars ; }
      const variableSet &targets() const { return target_vars ; }
      const variableSet &maps() const { return map_vars ; }
      const variableSet &constraints() const {return constraint_vars; }
      
      const string &qualifier() const { return internal_qualifier ; }
      rule_type type() const { return rule_class ; }
      time_ident target_time() const { return target_level ; }
      time_ident source_time() const { return source_level ; }
      time_ident time() const { return source_level ; }
      int ident() const { return rule::rdb->get_id(*this) ; }
      rule_implP get_rule_implP() const ;
    } ;
  private:
    friend class rule::info ;
    struct rule_db {
      std::vector<info> fiv ;
      std::map<std::string,int> fmap ;
      int get_id(const info &fi) {
        std::map<std::string,int>::iterator fmi ;
        if((fmi = fmap.find(fi.name())) == fmap.end()) {
          fiv.push_back(fi) ;
          int id = - fiv.size() ;
          fmap[fi.name()] = id ;
          return id ;
        }
        return fmi->second ;
      }
      const info &get_info(int id) const { return fiv[-(id+1)] ; }
    } ;
    static rule_db *rdb ;
    int id ;
    void create_rdb() {if(0==rdb) rdb = new rule::rule_db ; }
  public:
    rule() { create_rdb() ; id = rdb->get_id(info()) ;}
    explicit rule(int i)
      { create_rdb(); id = i ; }
    rule(const rule_implP &fp)
      { create_rdb(); id = rdb->get_id(info(fp)) ; }
    rule(rule f, time_ident tl)
      { create_rdb(); id = rdb->get_id(info(rdb->get_info(f.id),tl)) ; }
    rule(const std::string &s)
      { create_rdb(); id = rdb->get_id(info(s)) ; }
      
    rule parent() const { return rule(*this,time().parent()) ; }
      
    std::ostream &Print(std::ostream &s) const
      { s << rdb->get_info(id).name() ; return s ; }
      
    bool operator<(const rule &f) const { return id < f.id ; }
    bool operator==(const rule &f) const { return id == f.id ; }
    bool operator!=(const rule &f) const { return id != f.id ; }

    int ident() const { return id ; }
    const rule::info &get_info() const { return rdb->get_info(id) ; }

    const variableSet &sources() const { return rdb->get_info(id).sources(); }
    const variableSet &targets() const { return rdb->get_info(id).targets(); }
      
    rule_type type() const { return rdb->get_info(id).type() ; }
    time_ident target_time() const { return rdb->get_info(id).target_time() ;}
    time_ident source_time() const { return rdb->get_info(id).source_time() ;}
    time_ident time() const { return rdb->get_info(id).time() ; }
    rule_implP get_rule_implP() const { return rdb->get_info(id).get_rule_implP() ;}
  } ;

  inline std::ostream &operator<<(std::ostream &s, const rule &f)
    { return f.Print(s) ; }

  class ruleSet : public intervalSet {
  public:
    ruleSet() {}
    explicit ruleSet(const exprP &e) ;
    explicit ruleSet(const intervalSet &v)
      {*(static_cast<intervalSet *>(this)) = v ;}
    ruleSet &operator=(const intervalSet &v)
      {*(static_cast<intervalSet *>(this)) = v ; return *this ;}
    ruleSet &operator+=(const rule &v)
      { *this += v.ident() ; return *this ; }
    ruleSet &operator-=(const rule &v)
      { *this -= v.ident() ; return *this ; }
    bool inSet(const rule &v) const
      { return intervalSet::inSet(v.ident()) ; }
    class ruleSetIterator {
      intervalSet::const_iterator ii ;
    public:
      ruleSetIterator() {}
      ruleSetIterator(const intervalSet::const_iterator &i)
        { ii = i ; }
      rule operator*() const { return rule(*ii) ; }
      const rule::info *operator->() const
        { return &(rule(*ii).get_info()) ; }
      ruleSetIterator &operator++() { ++ii ; return *this ;}
      ruleSetIterator operator++(int )
        { return ruleSetIterator(ii++); }
      bool operator==(const ruleSetIterator &i) const 
        { return ii == i.ii ; }
      bool operator!=(const ruleSetIterator &i) const
        { return ii != i.ii ; } ;
    } ;
    typedef ruleSetIterator const_iterator ;
    const_iterator begin() const {
      return const_iterator(intervalSet::begin()) ; }
    const_iterator end() const {
      return const_iterator(intervalSet::end()) ; }
    std::ostream &Print(std::ostream &s) const ;

  } ;

  inline std::ostream &operator<<(std::ostream &s, const ruleSet v)
    { return v.Print(s) ; }

  class register_rule_type {
  public:
    virtual rule_implP get_func() const = 0 ;
  } ;

  class global_rule_impl_list {
  public:
    class rule_list_iterator ;
    friend class rule_list_iterator ;
  private:
    class rule_list_ent {
    public:
      rule_list_ent(register_rule_type *p, rule_list_ent *nxt) :
        rr(p), next(nxt) {}
      register_rule_type *rr ;
      rule_list_ent *next ;
    } ;
    static rule_list_ent *list ;
  public:
    class rule_list_iterator {
      rule_list_ent *p ;
    public:
      rule_list_iterator(rule_list_ent *ptr) : p(ptr) {}
      rule_implP operator*() { return p->rr->get_func() ; }
      rule_list_iterator &operator++() {
        p = p->next ;
        return *this ;
      }
      rule_list_iterator operator++(int ) {
        rule_list_iterator tmp(p) ;
        p = p->next ;
        return tmp ;
      }
      bool operator==(const rule_list_iterator &i) { return i.p == p ; }
      bool operator!=(const rule_list_iterator &i) { return i.p != p ; }
    } ;
    typedef rule_list_iterator iterator ;
        
    global_rule_impl_list() {}
    ~global_rule_impl_list() ;
    void push_rule(register_rule_type *rr) ;
    iterator begin() { return iterator(list) ; }
    iterator end() { return iterator(0) ; }
  } ;
    
  extern global_rule_impl_list global_rule_list ;    
    
  template<class T> class register_rule : public register_rule_type {
    //    copy_rule_impl<T> f ;
  public:
    register_rule() { global_rule_list.push_rule(this) ; }
    virtual rule_implP get_func() const { return new copy_rule_impl<T> ; }
  } ;
  
    
  class rule_db {
    typedef std::map<variable,ruleSet> varmap ;
    typedef varmap::const_iterator vc_iterator ;
    typedef std::map<std::string,rule> rule_map_type ;
      
    static const ruleSet EMPTY_RULE ;

    ruleSet known_rules ;

    rule_map_type name2rule ;
    varmap srcs2rule,trgt2rule ;

  public:
      
    void add_rule(const rule_implP &fp) ;
    void add_rule(rule f) ;
    void add_rules(global_rule_impl_list &gfl) ;
      
    rule_implP get_rule(const std::string &name) {
      return name2rule[name].get_info().rule_impl ;
    }
    const ruleSet &all_rules() const { return known_rules ; }
    const ruleSet &rules_by_source(variable v) const {
      vc_iterator vmi = srcs2rule.find(v) ;
      if(vmi == srcs2rule.end())
        return EMPTY_RULE ;
      else
        return vmi->second;
    }
    const ruleSet &rules_by_target(variable v) const {
      vc_iterator vmi = trgt2rule.find(v) ;
      if(vmi == trgt2rule.end())
        return EMPTY_RULE ;
      else
        return vmi->second;
    }
      
  } ;

}

#endif
