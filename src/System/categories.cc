#include "dist_tools.h"
#include <Tools/debug.h>
#include <entitySet.h>


#include "Tools/debugger.h"

#include <vector>
using std::vector ;
#include <set>
using std::set ;
#include <map>
using std::map ;
#include <list>
using std::list ;
#include <string>
using std::string ;

#include <iostream>
using std::cout ;
using std::cerr ; 
using std::endl ;
using std::ios ;
using std::ostream ;

#include <algorithm>
using std::sort ;

#define VERBOSE


namespace Loci {
 
  inline bool compare_var_sizes(const pair<variable,int> &p1, const pair<variable,int> &p2) {
    return p1.second > p2.second ;
  }

  inline bool compare_groups(const pair<variableSet,string> &p1,
                             const pair<variableSet,string> &p2) {
    return p1.second < p2.second ;
  }


  vector<entitySet> new_categories(map<variable, entitySet> &vm, vector<interval> &pvec) {

    // Create groups of variables that form categories
    // This is stored in groupings
    map<variable,entitySet>::const_iterator vmi ;
    map<variableSet, entitySet> groupings ;
    for(size_t i = 0;i<pvec.size();++i) {
      variableSet vss ;
      for(vmi=vm.begin();vmi!=vm.end();++vmi) {
        entitySet isect = vmi->second & entitySet(pvec[i]);
        if(isect != EMPTY && isect != entitySet(pvec[i])) {
          cerr << "something wrong with pvec" << endl ;
        }
        if(isect != EMPTY) 
          vss += vmi->first ;
      }
      if(vss != EMPTY) {
        groupings[vss] += pvec[i] ;
      }
    }


#ifdef VERBOSE
    Loci::debugout << "groupings.size() = " << groupings.size() << endl ;
#endif
    map<variableSet, entitySet>::const_iterator miter ;

    // Now we need to define an ordering of groupings that enables best memory utilization
    // We do this first by finding out how many entities are associated with each variable
    variableSet all_vars ;
    for(miter = groupings.begin(); miter != groupings.end(); ++miter) 
      all_vars += miter->first ;
    // all_vars gets all of the variables that are used in the groupings

    // We now get the number of entities associated with each variable in vsizes
    map<variable,int> vsizes ;
    map<variable,int>::const_iterator vsize_iter ;
    for(variableSet::const_iterator vi=all_vars.begin();vi!=all_vars.end();++vi)
      vsizes[*vi] = 0 ;
    for(miter = groupings.begin(); miter != groupings.end(); ++miter) 
      for(variableSet::const_iterator vi=miter->first.begin();vi!=miter->first.end();++vi)
        vsizes[*vi] += miter->second.size() ;

    // Now we create a sorted list of variables based on the number of entities
    vector<pair<variable,int> > vsort ;
    for(vsize_iter = vsizes.begin();vsize_iter!=vsizes.end();++vsize_iter) {
      vsort.push_back(*vsize_iter) ;
    }

    sort(vsort.begin(),vsort.end(),compare_var_sizes) ;


    // Now we sort groups such that the largest number of entities get grouped together first
    // We do this by associating a bit string with each variable, where the position in the
    // bit string relates to the size of the variable
    vector<pair<variableSet,string> > gsort ;

    for(miter = groupings.begin(); miter != groupings.end(); ++miter) {
      string s ;
      const variableSet &vs = miter->first ;
      for(size_t i=0;i<vsort.size();++i) {
        s += vs.inSet(vsort[i].first)?'1':'0' ;
      }
      gsort.push_back(make_pair(vs,s)) ;
    }
    // After sorting based on this generated key, the groups will have an ordering that
    // is good for memory layout.
    sort(gsort.begin(),gsort.end(),compare_groups) ;

    // Now create a list of entities in the order that will produce the best memory layout
    vector<entitySet> tmp_pvec ;
    for(size_t i=0;i<gsort.size();++i) {
      variableSet v = gsort[i].first ;
      entitySet e = groupings[v] ;
      tmp_pvec.push_back(e) ;
#ifdef VERBOSE
      Loci::debugout <<  "******************************************************" << endl ;
      Loci::debugout << " grouping variables " << v << endl ;
      Loci::debugout << " Entities shared = " << e << endl ;
      Loci::debugout << " Total Entities causing the grouping = " << e.size() << endl ;  
      Loci::debugout << "******************************************************" << endl ;
#endif
    }

    return tmp_pvec ;

  }
 
  class cat_tree : public CPTR_type {
  public:
    typedef CPTR<cat_tree>     treeP ;
    typedef list<treeP> treeList ;
    variableSet                generalization ;
    treeList                    before, after ;
    variableSet get_below_vars() ;
    int is_swapped ;
    cat_tree() { is_swapped = 0; }
    void Print(ostream &s) ;
    void fill_list(list<treeP> &list_tp) ;
    list<variableSet> get_categories(variableSet &vs) ;
    void recursive_reorder(list<variableSet> &lvs) ;
    void reorder(list<variableSet> &lvs) ;
    void down(int &i) ;
    bool operator==(const treeP &tp) const {
      if((generalization == tp->generalization) && (before.size() == tp->before.size()) && (after.size() == tp->after.size()))
        return 1 ;
      else
        return 0 ;
    }
    void before2after() ;
  } ;

  typedef cat_tree::treeP     treeP ;
  typedef cat_tree::treeList  treeList ;

  inline bool operator<(const treeP &tp1, const treeP &tp2) {
  
    if(tp1->generalization != tp2->generalization)
      return 0 ;
    else
      return  tp1->generalization < tp2->generalization ||
        (tp1->generalization==tp2->generalization &&
         (tp1->before.size()< tp2->before.size())) ; 
  }

  variableSet cat_tree::get_below_vars() {
    variableSet vs = generalization ;
    for(list<treeP>::iterator li = before.begin(); li != before.end(); ++li)
      vs += (*li)->get_below_vars() ;
    for(list<treeP>::iterator li = after.begin(); li != after.end(); ++li)
      vs += (*li)->get_below_vars() ;
    return vs ;
  }

  void cat_tree::Print(ostream &s) {
    for(list<treeP>::iterator li = before.begin(); li != before.end(); ++li)
      (*li)->Print(s) ;
    s << "  generalization = " << generalization << "  before.size = " << before.size() << "  after.size()  =  " << after.size() << endl ;
    for(list<treeP>::iterator li = after.begin(); li != after.end(); ++li)
      (*li)->Print(s) ;
  }

  void cat_tree::down(int &i ) {
    i++ ;
    for(list<treeP>::iterator li = before.begin(); li != before.end(); ++li)
      (*li)->down(i) ;
  }

  list<variableSet> cat_tree::get_categories(variableSet &vs) {
    vs += generalization ;
    variableSet tvs  = vs ;
    list<variableSet> lvs ;
  
    if(before.size() > 1) {
      for(list<treeP>::iterator li = before.begin(); li != before.end(); ++li) {
        variableSet tmp_vs = vs;
        tmp_vs += (*li)->generalization ;
        list<variableSet> tmp_lvs = (*li)->get_categories(tmp_vs) ;
        for(list<variableSet>::iterator lvi = tmp_lvs.begin(); lvi != tmp_lvs.end(); ++lvi)
          lvs.push_back(*lvi) ;
      }
    }
    else {
      for(list<treeP>::iterator li = before.begin(); li != before.end(); ++li) {
        variableSet tmp_vs = vs;
        list<variableSet> tmp_lvs ;
        if(after.size()) {
          tmp_vs += (*li)->generalization ;
          tmp_lvs = (*li)->get_categories(tmp_vs) ;
          for(list<variableSet>::iterator lvi = tmp_lvs.begin(); lvi != tmp_lvs.end(); ++lvi)
            lvs.push_back(*lvi) ;
        }
        else {
          vs += (*li)->generalization ;
          tmp_lvs = (*li)->get_categories(vs) ;
          for(list<variableSet>::iterator lvi = tmp_lvs.begin(); lvi != tmp_lvs.end(); ++lvi)
            lvs.push_back(*lvi) ;
        }  
      }
    }
    if(tvs != EMPTY)
      lvs.push_back(tvs) ;
    if(after.size() > 1) {
      for(list<treeP>::iterator li = after.begin(); li != after.end(); ++li) {
        variableSet tmp_vs = tvs;
        tmp_vs += (*li)->generalization ;
        list<variableSet> tmp_lvs = (*li)->get_categories(tmp_vs) ;
        for(list<variableSet>::iterator lvi = tmp_lvs.begin(); lvi != tmp_lvs.end(); ++lvi)
          lvs.push_back(*lvi) ;
      }
    }
    else {
      for(list<treeP>::iterator li = after.begin(); li != after.end(); ++li) {
        variableSet tmp_vs = vs ;
        list<variableSet> tmp_lvs ;
        if(before.size()) {
          tmp_vs += (*li)->generalization ;
          tmp_lvs = (*li)->get_categories(tmp_vs) ;
        }
        else {
          vs += (*li)->generalization ;
          tmp_lvs = (*li)->get_categories(vs) ;
        }
        for(list<variableSet>::iterator lvi = tmp_lvs.begin(); lvi != tmp_lvs.end(); ++lvi)
          lvs.push_back(*lvi) ;
      }  
    }
  
    return lvs ;
  }

  void cat_tree::fill_list(list<treeP> &list_tp) {
    for(list<treeP>::iterator li = before.begin(); li != before.end(); ++li)
      (*li)->fill_list(list_tp) ;
    list_tp.push_front(this) ;
    for(list<treeP>::iterator li = after.begin(); li != after.end(); ++li)
      (*li)->fill_list(list_tp) ;
  }

  void cat_tree::before2after() {
    list<treeP> vt ;
    for(list<treeP>::iterator li = before.begin(); li != before.end(); ++li) {
      if((*li)->before.size())
        vt.push_back(*li) ;
      else
        after.push_back(*li) ;
    }
    before = vt ;
  }

  void cat_tree::reorder(list<variableSet> &lvs) {
    std::set<variableSet> tmp_vars ;
    list<treeP> lone_ltp, tmp_buf, tmp_before ;
    for(list<treeP>::iterator li = before.begin(); li != before.end(); ++li)
      if((*li)->before.size() || (*li)->after.size())
        tmp_buf.push_back(*li) ;
      else
        lone_ltp.push_back(*li) ;
    for(list<treeP>::iterator li = tmp_buf.begin(); li != tmp_buf.end(); ++li) {
      tmp_before.push_back(*li) ;
      tmp_vars.insert((*li)->generalization) ;
      variableSet vset = (*li)->get_below_vars() ;
      for(list<treeP>::iterator tmp_li = lone_ltp.begin(); tmp_li != lone_ltp.end(); ++tmp_li)
        if((((*tmp_li)->generalization) & vset) != EMPTY)
          if(tmp_vars.find((*tmp_li)->generalization) == tmp_vars.end()) {
            tmp_before.push_back(*tmp_li) ;
            tmp_vars.insert((*tmp_li)->generalization) ;
          }
    }
    for(list<treeP>::iterator tmp_li = lone_ltp.begin(); tmp_li != lone_ltp.end(); ++tmp_li)
      if(tmp_vars.find((*tmp_li)->generalization) == tmp_vars.end()) {
        tmp_vars.insert((*tmp_li)->generalization) ;
        tmp_before.push_back(*tmp_li) ;
      }
    before = tmp_before ;
    variableSet int_vset ;
    for(list<variableSet>::iterator li = lvs.begin(); li != lvs.end(); ++li) {
      for(variableSet::const_iterator vi = (*li).begin(); vi != (*li).end(); ++vi)
        int_vset += *vi ;
    }
    for(list<treeP>::iterator lii = before.begin(); lii != before.end(); ++lii) {
      entitySet vs = entitySet((*lii)->get_below_vars()) ;
      for(list<treeP>::iterator lij = lii; lij != before.end(); ++lij) {
        if(lii != lij) {
          entitySet next_vs =  entitySet((*lij)->get_below_vars()) ;
          entitySet tmp_vs = vs & next_vs ;
          tmp_vs -= (tmp_vs & entitySet(int_vset)) ;
          if(tmp_vs != EMPTY) {
            lvs.push_front(variableSet(tmp_vs)) ;
            int_vset += variableSet(tmp_vs) ;
            list<treeP> tlp, tmp_tlp ;
            for(list<treeP>::iterator tlpi = (*lij)->before.begin(); tlpi != (*lij)->before.end(); ++tlpi)
              if(((*tlpi)->get_below_vars() & variableSet(tmp_vs)) == EMPTY)
                (*lij)->after.push_back(*tlpi) ;
              else
                tlp.push_back(*tlpi) ;
            (*lij)->before = tlp ; 
            (*lij)->is_swapped = 1 ;
          }
        }
      }
      for(list<treeP>::iterator lij = before.begin(); lij != lii; ++lij) {
        if(!(*lij)->is_swapped) {
          entitySet vs = entitySet((*lij)->get_below_vars()) ;
          for(list<variableSet>::iterator lvi = lvs.begin(); lvi != lvs.end(); ++lvi) {
            entitySet next_vs =  entitySet(*lvi) ;
            entitySet tmp_vs = vs & next_vs ;
            if(tmp_vs != EMPTY) {
              list<treeP> tlp, tmp_tlp ;
	    
              for(list<treeP>::iterator tlpi = (*lij)->before.begin(); tlpi != (*lij)->before.end(); ++tlpi)
                if(((*tlpi)->get_below_vars() & variableSet(tmp_vs)) == EMPTY)
                  tmp_tlp.push_back(*tlpi) ;
                else
                  (*lij)->after.push_front(*tlpi) ;
              (*lij)->before = tmp_tlp ; 
	  
            }
          }
        }
      }
    }
  }


  void cat_tree::recursive_reorder(list<variableSet> &tmp_lvs) {
    for(list<treeP>::iterator li = before.begin(); li != before.end(); ++li)
      (*li)->recursive_reorder(tmp_lvs) ;
    reorder(tmp_lvs) ;
    for(list<treeP>::iterator li = after.begin(); li != after.end(); ++li)
      (*li)->recursive_reorder(tmp_lvs) ;
  } 
 
  vector<vector<variableSet> > create_orig_matrix(vector<variableSet> &ves) {
    vector<variableSet>::const_iterator vesi, viter ;
    vector<vector<variableSet> > vvvs(ves.size()) ;
    unsigned int tmp = 0 ;
    viter = ves.begin() ;
    while(tmp < ves.size()) { 
      entitySet vs ;
      for(vesi = ves.begin(); vesi != ves.end(); ++vesi) { 
        vs = (entitySet(*viter) & entitySet(*vesi)) ;
        vvvs[tmp].push_back(variableSet(vs));
      }
      tmp++ ;
      viter++ ;
    }
    /*
      if(Loci::MPI_rank == 0) {
      cout << " original_matrix =  " << endl ;
      for(int i = 0; i < vvvs.size(); ++i) {
      for(int j = 0; j < vvvs[i].size(); ++j) 
      cout << vvvs[i][j] << "   "  ; 
      cout << endl ;
      }
      }
    */
    return vvvs ;
  }

  vector<vector<variableSet> > create_sub_matrix(vector<variableSet> &ves, variableSet v) {
    vector<variableSet>::const_iterator vesi, viter ;
    vector<vector<variableSet> > vvvs ;
    vector<variableSet>  tmp_ves(ves.size()) ;
    unsigned int tmp = 0 ;
    for(unsigned int i = 0; i < ves.size(); ++i) 
      tmp_ves[i] = variableSet(entitySet(ves[i])- entitySet(v)) ;
    viter = tmp_ves.begin() ;
    while(tmp < tmp_ves.size()) { 
      entitySet vs ;
      vector<variableSet> tmp_vvs ;
      for(vesi = tmp_ves.begin(); vesi != tmp_ves.end(); ++vesi) {
        vs = (entitySet(*viter) & entitySet(*vesi)) ;
        tmp_vvs.push_back(variableSet(vs));
      }
      vvvs.push_back(tmp_vvs) ;
      tmp++ ;
      viter++ ;
    }
    /*
      if(Loci::MPI_rank == 0) {
      cout << " Sub matrix =  " << endl ;
      for(int i = 0; i < vvvs.size(); ++i) {
      for(int j = 0; j < vvvs[i].size(); ++j) 
      cout << vvvs[i][j] << "   "  ; 
      cout << endl ;
      }
      }
    */
    return vvvs ;
  }

  treeP recursive_matrix_tree(vector<vector<variableSet> > &tmp_vvvs, variableSet vset, std::set<variableSet> total_varset) {
    vector<vector<variableSet> > vvvs ;
    treeP tp = new cat_tree ;
    variableSet vs ;
    tp->generalization = vset ;
    tp->is_swapped = 0 ;
    if(tmp_vvvs.size() == 0) {
      return tp ;
    }
    vector<variableSet> tmp_vvs ;
    for(unsigned int i = 0; i < tmp_vvvs.size(); ++i) {
      if(tmp_vvvs[i][i] != EMPTY) {
        tmp_vvs.push_back(tmp_vvvs[i][i]) ;
      }
    }
    if(tmp_vvs.size() < tmp_vvvs.size())
      vvvs = create_orig_matrix(tmp_vvs) ;
    else
      vvvs = tmp_vvvs ;

    for(unsigned int i = 0; i < vvvs.size(); ++i) {
      std::set<variableSet> set_vars ;
      vector<variableSet> vvs ;
      entitySet es ;
      es = entitySet(vvvs[i][i]) ;
      if(set_vars.find(vvvs[i][i]) == set_vars.end()) {
        vvs.push_back(vvvs[i][i]) ;
        set_vars.insert(vvvs[i][i]) ;
      }
      vector<vector<variableSet> > sub_vvvs ;
      for(unsigned int j = 0; j < vvvs[i].size(); ++j)
        if(j != i) {
          if(vvvs[i][j] != EMPTY) {
            es &= entitySet(vvvs[i][j]) ;
            vs = variableSet(es) ;
            if(vs != EMPTY) {
              if(set_vars.find(vvvs[j][j]) == set_vars.end()) {
                vvs.push_back(vvvs[j][j]) ;
                set_vars.insert(vvvs[j][j]) ;
              }
            }
          }
        }
      if(vvs.size()) {
        es = vvs[0] ;
        for(unsigned int j = 1; j < vvs.size(); ++j)
          es &= entitySet(vvs[j]) ;
        if(total_varset.find(variableSet(es)) == total_varset.end()) {
          total_varset.insert(variableSet(es)) ;
          sub_vvvs = create_sub_matrix(vvs, variableSet(es)) ;
          tp->before.push_back(recursive_matrix_tree(sub_vvvs,
                                                     variableSet(es), total_varset));
        }
      }
      if(set_vars.size() == vvvs.size() )
        break ;
    }
    return tp ;
  }

  vector<entitySet> modified_categories(fact_db &facts, map<variable, entitySet> &vm, vector<interval> &pvec) {
    vector<entitySet> tmp_pvec ;
    vector<variableSet> vvs(pvec.size()) ;
    map<variableSet, entitySet> mve ;
    //map<variableSet, vector<int> > mv_vec ;
    //map<variableSet, vector<int> >::iterator mv_iter ;
    map<variableSet, entitySet>::const_iterator miter ;
    map<variable, entitySet>::const_iterator svi ;
    variableSet initial_varset ;
    //  double start_time = MPI_Wtime() ;
#ifdef VERBOSE
    Loci::debugout << " Size of the vector passed to modified categories = " << pvec.size() << endl ;
#endif
    for(unsigned int i = 0; i < pvec.size(); ++i) {
      for(svi = vm.begin(); svi != vm.end(); ++svi) {
        if((svi->second & entitySet(pvec[i])) != EMPTY) {
          vvs[i] += svi->first ;
          initial_varset += svi->first ;
        }
      }
      if(vvs[i] != EMPTY) {
        mve[vvs[i]] += pvec[i]; 
      }
    }

    vector<variableSet> tmp_vvs ;
  
    for(miter = mve.begin(); miter != mve.end(); ++miter) {
      tmp_vvs.push_back(miter->first) ;
#ifdef VERBOSE
      Loci::debugout << "******************************************************" << endl ;
      Loci::debugout << " grouping variables " << miter->first << endl ;
      Loci::debugout << " Entities shared = " << miter->second << endl ;
      Loci::debugout << " Total Entities causing the grouping = " << miter->second.size() << endl ;  
      Loci::debugout << "******************************************************" << endl ;
#endif
    }

#ifdef VERBOSE
    Loci::debugout << " The number of variable sets grouped due to common categories = " << tmp_vvs.size() << endl ;
#endif
    vector<vector<variableSet> > vvvs = create_orig_matrix(tmp_vvs) ;
    vvs.clear() ;
    vector<variableSet> indep ;
    std::set<variableSet> total_vars, itotal_vars ;
  
    for(unsigned int i = 0; i < vvvs.size(); ++i) {
      entitySet vs, tmp_vs ;
      for(unsigned int j = 0; j < vvvs[i].size(); ++j) {
        if(j != i) {
          vs = entitySet(vvvs[i][i]) & entitySet(vvvs[i][j]) ;
          tmp_vs += vs ;
          if(vs != EMPTY) {
            if(total_vars.find(vvvs[j][j]) == total_vars.end()) {
              vvs.push_back(vvvs[j][j]) ;
              total_vars.insert(vvvs[j][j]) ;
            }
          }
        }
      }
      if(tmp_vs == EMPTY) 
        indep.push_back(vvvs[i][i]) ;
    }
    vector<vector<variableSet> > tmp_vvvs ;
    vvvs = create_orig_matrix(vvs) ;
    for(unsigned int i = 0; i < vvvs.size(); ++i) {
      vvs.clear() ;
      entitySet es = entitySet(vvvs[i][i]);
      for(unsigned int j = 0; j < vvvs[i].size(); ++j) {
        entitySet tmp_es ;
        if((vvvs[i][j] != EMPTY))  {
          if(itotal_vars.find(vvvs[j][j]) == itotal_vars.end()) {
            tmp_es = es & entitySet(vvvs[i][j]) ;
            if(tmp_es != EMPTY) {
              es = tmp_es ;
              itotal_vars.insert(vvvs[j][j]) ;
              vvs.push_back(vvvs[j][j]) ;
            }
          }
          else
            j++ ;
        }
      }
      if(vvs.size())
        tmp_vvvs.push_back(vvs) ;
      if(itotal_vars.size() == vvvs.size())
        break ; 
    }
    treeP super_tp = new cat_tree ;
    variableSet vset ;
    super_tp->generalization = vset ;
    for(unsigned int i = 0; i < tmp_vvvs.size(); ++i) {
      treeP tp = new cat_tree ;
      variableSet vs ;
      vvvs = create_orig_matrix(tmp_vvvs[i]) ;
      std::set<variableSet> total_varset ;
      tp = recursive_matrix_tree(vvvs, vs, total_varset) ;
      super_tp->before.push_back(tp) ;
    }
    treeP tp = new cat_tree ;
    tp->generalization = vset ;
    for(unsigned int i = 0 ; i < indep.size(); ++i) {
      treeP tmp_tp = new cat_tree ;
      tmp_tp->generalization = indep[i] ;
      tp->before.push_back(tmp_tp) ;
    } 
    for(list<treeP>::iterator li = super_tp->before.begin(); li != super_tp->before.end(); ++li) 
      for(list<treeP>::iterator lj = (*li)->before.begin(); lj != (*li)->before.end(); ++lj)
        tp->before.push_back(*lj) ;
  
    list<variableSet> lvs = tp->get_categories(vset) ; 
    list<variableSet> tlist ;
    tp->recursive_reorder(tlist) ;
  
    lvs = tp->get_categories(vset) ;
    entitySet total_entities ;
    HASH_MAP(int, entitySet) map_vec ;
    vector<int> sort_vec ;
    for(list<variableSet>::iterator lvi = lvs.begin(); lvi != lvs.end(); ++lvi) {
      if(mve.find(*lvi) != mve.end()) {
        entitySet tmp = mve[*lvi] ; 
        if(tmp!= EMPTY) {
          //map_vec[*ei] = tmp ;
          //sort_vec.push_back(*ei) ;
          //Loci::debugout << " Correspoding to " << *lvi << " pushing back " << tmp << endl ;
          tmp -= total_entities ;
          if(tmp!= EMPTY) {
            tmp_pvec.push_back(tmp) ;
            //ei = tmp.begin() ;
            //map_vec[*ei] = tmp ;
            //sort_vec.push_back(*ei) ;
          }
          total_entities += tmp ;
        }
      }
    }

#ifdef VERBOSE
    for(size_t i = 0; i < tmp_pvec.size(); ++i)
      Loci::debugout << " tmp_pvec[" << i << " ] = " << tmp_pvec[i] << endl ;
    Loci::debugout << " total_entities = " << total_entities << endl ;
#endif
    //vector<int>::const_iterator svci = sort_vec.begin();
    /*
      while(svci != sort_vec.end()) {
      entitySet tmp = map_vec[*svci] ; 
      ++svci ;
      while(svci != sort_vec.end()) 
      if(tmp.Max() > map_vec[*svci].Min()) {
      tmp += map_vec[*svci] ;
      ++svci ;
      }
      else 
      break ;
      tmp_pvec.push_back(tmp) ;
      }
    */
    //for(vector<int>::const_iterator svci = sort_vec.begin(); svci != sort_vec.end(); ++svci)
    //tmp_pvec.push_back(map_vec[*svci]) ;
  
    variableSet vars = facts.get_typed_variables() ;
    for(variableSet::const_iterator vsi = vars.begin(); vsi != vars.end(); ++vsi) {
      storeRepP p = facts.get_variable(*vsi) ;
      if((p->RepType() == MAP)) {
        MapRepP mp = MapRepP(p->getRep()) ;
        entitySet image_dom = mp->image(p->domain()) ;
        entitySet out_of_set = image_dom - total_entities ;
        entitySet left_out_categories = all_collect_entitySet(out_of_set) ;
        if(left_out_categories != EMPTY) {
          Loci::debugout << " left out stuff  = " << left_out_categories  << endl ;
          //ei = left_out_categories.begin() ;
          //map_vec[*ei] = left_out_categories ;
          //sort_vec.push_back(*ei) ;
          tmp_pvec.push_back(left_out_categories) ;
          total_entities += left_out_categories ;
        }
      }
    }
    //std::sort(sort_vec.begin(), sort_vec.end()) ; 
    //for(vector<int>::const_iterator svci = sort_vec.begin(); svci != sort_vec.end(); ++svci)
    //tmp_pvec.push_back(map_vec[*svci]) ;
  
    return tmp_pvec ;
  }

  void categories(fact_db &facts,vector<entitySet> &pvec) {
    entitySet active_set ;
    set<entitySet> set_of_sets ;
    set_of_sets.insert(~EMPTY) ;
    variableSet vars = facts.get_typed_variables() ;
    for(variableSet::const_iterator vi=vars.begin();vi!=vars.end();++vi) {
      storeRepP p = facts.get_variable(*vi) ;
      if(p->RepType() == MAP) {
        MapRepP mp = MapRepP(p->getRep()) ;
        active_set += p->domain() ;
        active_set += mp->image(p->domain()) ;
        set_of_sets.insert(p->domain()) ;
      } else if(p->RepType() == STORE) {
        active_set += p->domain() ;
        set_of_sets.insert(p->domain()) ;
      } else {
        if(p->domain() != ~EMPTY) {
	  entitySet tmp = p->domain() ;
	  if(facts.is_distributed_start())
	    tmp = all_collect_entitySet(tmp) ;
	  active_set += tmp ;
	  set_of_sets.insert(p->domain()) ;
	}
      }
    }
    vector<int> vals,vals2 ;
    set<entitySet>::iterator si ;
    for(si=set_of_sets.begin();si!=set_of_sets.end();++si) {
      entitySet s = *si & active_set ;
      for(int i = 0;i < s.num_intervals(); ++i) {
	vals.push_back(s[i].first-1) ;
	vals.push_back(s[i].first) ;
	vals.push_back(s[i].second) ;
	vals.push_back(s[i].second+1) ;
      }
    }
    
    std::sort(vals.begin(),vals.end()) ;
    
    vals2.push_back(vals[0]) ;
    for(unsigned int i=1;i<vals.size();++i) {
      if(vals[i-1] == vals[i])
        continue ;
      vals2.push_back(vals[i]) ;
    }
    vector<interval> tmp_pvec ;
    FATAL(vals2.size() < 4) ;
    for(int i=1;i<int(vals2.size())-1;) {
      int i1 = vals2[i] ;
      int i2 = vals2[i+1] ;
      if(i1+1 == i2)
        i2 = i1 ;
      else 
        ++i ;
      ++i ;
      interval iv(i1,i2) ;
      if(facts.is_distributed_start())
	tmp_pvec.push_back(iv) ;
      else
	pvec.push_back(entitySet(iv)) ;
    }
    
    if(facts.is_distributed_start()) {
      map<variable, entitySet> vm ;
      entitySet total_entities ;
      for(variableSet::const_iterator vi=vars.begin();vi!=vars.end();++vi) {
	entitySet active_set ;
	storeRepP p = facts.get_variable(*vi) ;
	if((p->RepType() == MAP)) {
	  entitySet tmp = p->domain() ;
	  entitySet tmp_all_collect = Loci::all_collect_entitySet(tmp) ;
	  vm[*vi] = tmp_all_collect ; 
	  total_entities += tmp_all_collect ;
	}
	else if((p->RepType() == STORE)) {
	  entitySet tmp = p->domain() ;
	  entitySet all_collect = Loci::all_collect_entitySet(tmp) ;
	  vm[*vi] = all_collect ;
	  total_entities += all_collect ;
	}
	else {
	  if(p->domain() != ~EMPTY) {
	    entitySet tmp = p->domain() ;
	    entitySet all_collect = Loci::all_collect_entitySet(tmp) ;
	    vm[*vi] = all_collect ; 
	    total_entities += all_collect ;
	  }
	}
      }
      for(variableSet::const_iterator vsi = vars.begin(); vsi != vars.end(); ++vsi) {
	storeRepP p = facts.get_variable(*vsi) ;
	if((p->RepType() == MAP)) {
	  MapRepP mp = MapRepP(p->getRep()) ;
	  entitySet image_dom = mp->image(p->domain()) ;
	  entitySet out_of_set = image_dom - total_entities ;
	  entitySet left_out_categories = all_collect_entitySet(out_of_set) ;
	  if(left_out_categories != EMPTY) {
	    std::string name = "image_" ;
	    name.append(vsi->get_info().name) ;
	    variable v = variable(name) ;
	    vm[v] = left_out_categories ; 
	  }
	}
      } 
      //      pvec = modified_categories(facts, vm, tmp_pvec) ;
      pvec = new_categories(vm, tmp_pvec) ;      
    }
  }
}
