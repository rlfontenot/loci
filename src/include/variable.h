#ifndef VARIABLE_H
#define VARIABLE_H

#include <ostream>
#include <iostream>
#include <sstream>

#include <map>
#include <set>
#include <vector>
#include <list>

#include <Tools/expr.h>
#include <Tools/intervalSet.h>

namespace Loci {
  class time_ident {
    struct time_hierarchy {
      struct time_info {
        std::string level_name ;
        int parent_id ;
        std::vector<int> children ;
        time_info() { level_name = "*"; parent_id = 0 ; }
        time_info(const std::string &name, int parent)
          { level_name = name; parent_id = parent ; }
      } ;
      std::vector<time_info> time_db ;
      int add_level(const std::string &lname,int level) ;
      time_hierarchy()
        { time_db.push_back(time_info()) ; }
      int parent(int i) const{ return time_db[i].parent_id ; }
      const std::string &name(int i) const { return time_db[i].level_name ; }
      const std::vector<int> get_children(int i) const
        { return time_db[i].children; } ;
    } ;
  private:
    static time_hierarchy *thp ;
    int id ;
    void create_thp() {if(0 == thp) thp = new time_hierarchy ; }
  public:
    time_ident() { create_thp() ; id = 0 ; }
    explicit time_ident(int i) { create_thp() ; id = i ; }
    time_ident(const time_ident &t) { create_thp() ; id = t.id ; }
    time_ident(const exprP &exp) ;
    time_ident parent() const { return time_ident(thp->parent(id)) ; }
    int ident() const { return id ; }
    std::vector<time_ident> children() ;
    const std::string &level_name() const { return thp->name(id) ; }
    bool before(const time_ident &t) const ;
    bool operator<(const time_ident &t) const { return id < t.id ; }
    bool operator==(const time_ident &t) const { return id == t.id ; }
    bool operator!=(const time_ident &t) const { return id != t.id ; }
    std::ostream &Print(std::ostream &s) const ;
  } ;

  inline std::ostream & operator<<(std::ostream &s, const time_ident &ti) {
    ti.Print(s) ;
    return s;
  }




  template <class Key> class key_ident {
    typedef typename std::map<Key,int>::iterator map_iterator ;
    std::map<Key, int> key_map ;
    std::vector<Key> key_vec ;
  public:
    int get_id(const Key &k) {
      map_iterator i ;
      if((i=key_map.find(k)) != key_map.end()) 
        return i->second ;
      else {
        int val = (key_map[k] = key_vec.size()) ;
        key_vec.push_back(k) ;
        return val ;
      }
    }
    int size() const { return key_vec.size() ; }
    const Key &operator[](int id) const { return key_vec[id] ; }
    int operator[](const Key &k) const { return key_map[k] ; }
  } ;

  class variable {
  public:
    struct info {
      bool         tvar ;
      bool         assign ;
      std::string  name ;
      time_ident   time_id ;
      int          offset ;
      std::vector<std::string> priority ;
            
      bool operator<(const info &v) const ;
      bool operator==(const info &v) const ;
            
      info() {
        tvar    = false;
        assign  = false;
        name    = "*NONAME*" ;
        offset  = 0;
      }
      ostream &Print(ostream &s) const ;
      const time_ident & time() const { return time_id ; }
      const info & get_info() const { return *this ; }
      variable parent() const ;
      variable drop_assign() const ;
      variable drop_priority() const ;
      variable new_offset(int o) const ;
      int ident() const { return variable::vdb->vars.get_id(*this) ; }
    };
  private:
    friend class variable::info ;
    struct variable_db {
      key_ident<info> vars ;
    } ;
    static variable_db *vdb ;
    int id ;
    void create_vdb() {if(0 == vdb) vdb = new variable_db ; }
    variable(const info &v) { create_vdb() ; id = vdb->vars.get_id(v) ; }
  public:
    variable() { create_vdb() ; id = vdb->vars.get_id(info()) ; }
    explicit variable(int i) { create_vdb() ; id = i ; }
    explicit variable(const exprP &p) ;
    explicit variable(string s)
    { id = variable(expression::create(s)).ident() ; }
    explicit variable(const time_ident &t) ;
    explicit variable(const variable &v, const time_ident &t) ;
        
    ostream &Print(std::ostream &s) const { return vdb->vars[id].Print(s) ; }
    std::string str() const
      { std::ostringstream ss ; Print(ss) ; return ss.str() ; }

    int ident() const { return id ; }
    bool operator<(const variable &v) const { return id < v.id; }
    bool operator==(const variable &v) const { return id == v.id; }
    bool operator!=(const variable &v) const { return id != v.id; }
    const time_ident & time() const { return vdb->vars[id].time_id ; }
    const info & get_info() const { return vdb->vars[id]; }

    variable parent() const { return vdb->vars[id].parent() ;}
    variable drop_assign() const { return vdb->vars[id].drop_assign() ; }
    variable drop_priority() const { return vdb->vars[id].drop_priority() ; }
    variable new_offset(int o) const { return vdb->vars[id].new_offset(o) ; }
  } ;

  inline ostream &operator<<(ostream &s, const variable &v)
    { return v.Print(s) ; }

  class variableSet : public intervalSet {
  public:
    variableSet() {}
    explicit variableSet(const exprP &e) ;
    explicit variableSet(const intervalSet &v)
      {*(static_cast<intervalSet *>(this)) = v ;}
    variableSet &operator=(const intervalSet &v)
      {*(static_cast<intervalSet *>(this)) = v ; return *this ;}
    variableSet &operator+=(const variable &v)
      { *this += v.ident() ; return *this ; }
    variableSet &operator-=(const variable &v)
      { *this -= v.ident() ; return *this ; }
    bool inSet(const variable &v) const
      { return intervalSet::inSet(v.ident()) ; }
    class variableSetIterator {
      intervalSet::const_iterator ii ;
    public:
      variableSetIterator() {}
      variableSetIterator(const intervalSet::const_iterator &i)
        { ii = i ; }
      variable operator*() const { return variable(*ii) ; }
      const variable::info *operator->() const
        { return &(variable(*ii).get_info()) ; }
      variableSetIterator &operator++() { ++ii ; return *this ;}
      variableSetIterator operator++(int )
        { return variableSetIterator(ii++); }
      bool operator==(const variableSetIterator &i) const 
        { return ii == i.ii ; }
      bool operator!=(const variableSetIterator &i) const
        { return ii != i.ii ; } ;
    } ;
    typedef variableSetIterator const_iterator ;
    const_iterator begin() const {
      return const_iterator(intervalSet::begin()) ; }
    const_iterator end() const {
      return const_iterator(intervalSet::end()) ; }
    std::ostream &Print(std::ostream &s) const ;

  } ;

  inline std::ostream &operator<<(std::ostream &s, const variableSet v)
    { return v.Print(s) ; }


  struct vmap_info {
    std::vector<variableSet> mapping ;
    variableSet var ;
    std::vector<std::pair<variable,variable> > assign ;
    bool operator<(const vmap_info &di) const {
      if(mapping < di.mapping)
        return true ;
      if(!(di.mapping < mapping) && var < di.var)
        return true ;
      return false ;
    }
    vmap_info() {}
    vmap_info(const exprP &e) ;
    std::ostream &Print(std::ostream &s) const ;
  } ;

  inline std::ostream &operator<<(std::ostream &s, const vmap_info &di) {
    return di.Print(s) ;
  }


  inline std::ostream &operator<<(std::ostream &s, const std::set<vmap_info> &v) {
    std::set<Loci::vmap_info>::const_iterator i ;
    for(i = v.begin();i!=v.end();) {
      s << (*i) ;
      ++i ;
      if(i!=v.end())
        s << "," ;
    }
    return s ;
  }
}    

#endif
