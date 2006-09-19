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
using std::pair ;
using std::make_pair ;

using std::ostringstream ;
//#define VERBOSE


namespace Loci {

  // Collect entitities to a unified entitySet that is distributed across processors according to the partition ptn.
  entitySet dist_collect_entitySet(entitySet inSet, const vector<entitySet> &ptn) {
    const int p = MPI_processes ;
    vector<entitySet> distSet(p) ;
    for(int i=0;i<p;++i) 
      distSet[i] = inSet & ptn[i] ;

    vector<int> send_sz(p) ;
    for(int i=0;i<p;++i)
      if(i!=MPI_rank)
        send_sz[i] = distSet[i].num_intervals()*2 ;
      else
        send_sz[i] = 0 ;
    
    vector<int> recv_sz(p) ;
    MPI_Alltoall(&send_sz[0],1,MPI_INT,
                 &recv_sz[0],1,MPI_INT,
                 MPI_COMM_WORLD) ;

    int size_send = 0 ;
    int size_recv = 0 ;
    for(int i=0;i<p;++i) {
      size_send += send_sz[i] ;
      size_recv += recv_sz[i] ;
    }

    int *send_store = new int[size_send] ;
    int *recv_store = new int[size_recv] ;
    int *send_displacement = new int[p] ;
    int *recv_displacement = new int[p] ;

    send_displacement[0] = 0 ;
    recv_displacement[0] = 0 ;
    for(int i = 1; i <  p; ++i) {
      send_displacement[i] = send_displacement[i-1] + send_sz[i-1] ;
      recv_displacement[i] = recv_displacement[i-1] + recv_sz[i-1] ;
    }
    for(int i = 0; i <  p; ++i)
      if(i != MPI_rank)
        for(int j=0;j<distSet[i].num_intervals();++j) {
          send_store[send_displacement[i]+j*2] = distSet[i][j].first ;
          send_store[send_displacement[i]+j*2+1] = distSet[i][j].second ;
        }
    
    
    MPI_Alltoallv(send_store,&send_sz[0], send_displacement , MPI_INT,
		  recv_store, &recv_sz[0], recv_displacement, MPI_INT,
		  MPI_COMM_WORLD) ;  

    entitySet retval = distSet[MPI_rank] ;
    for(int i = 0; i <  p; ++i) 
      for(int j=0;j<recv_sz[i]/2;++j) {
        int i1 = recv_store[recv_displacement[i]+j*2]  ;
        int i2 = recv_store[recv_displacement[i]+j*2+1] ;
        retval += interval(i1,i2) ;
      }
    delete[] recv_displacement ;
    delete[] send_displacement ;
    delete[] recv_store ;
    delete[] send_store ;

    return retval ;
  }
  
  
  inline bool compare_var_sizes(const pair<variable,int> &p1, const pair<variable,int> &p2) {
    return p1.second > p2.second ;
  }

  inline bool compare_groups(const pair<variableSet,string> &p1,
                             const pair<variableSet,string> &p2) {
    return p1.second < p2.second ;
  }


  // Compute associations between variables and their domains
  // Note, we also create image variables to represent the images of maps
  //   over their domains.
  void getVariableAssociations(map<variable,entitySet> &vm, fact_db &facts) {
    vector<entitySet> ptn = facts.get_init_ptn() ;
    entitySet total_entities ;
    vm.clear() ;
    variableSet vars = facts.get_typed_variables() ;
    for(variableSet::const_iterator vi=vars.begin();vi!=vars.end();++vi) {
      entitySet active_set ;
      storeRepP p = facts.get_variable(*vi) ;
      if((p->RepType() == MAP)) {
        entitySet tmp = dist_collect_entitySet(p->domain(), ptn) ;
        vm[*vi] = tmp ;
        total_entities += tmp ;
        MapRepP mp = MapRepP(p->getRep()) ;
        entitySet image_dom = mp->image(p->domain()) ;
        std::string name = "image_" ;
        name.append(vi->get_info().name) ;
        variable v = variable(name) ;
        image_dom = dist_collect_entitySet(image_dom, ptn) ;
        vm[v] = image_dom ;
        total_entities += image_dom ;
      } else if((p->RepType() == STORE)) {
        entitySet tmp = dist_collect_entitySet(p->domain(), ptn) ;
        vm[*vi] = tmp ;
        total_entities += tmp ;
      } else {
        if(p->domain() != ~EMPTY) {
          entitySet all_collect = dist_collect_entitySet(p->domain(),ptn) ;
          vm[*vi] = all_collect ; 
          total_entities += all_collect ;
        }
      }
    }
    variable v_leftout("SENTINEL_LEFTOUT") ;
    vm[v_leftout] = ptn[MPI_rank]-total_entities ;
  }

  // Compute the distict intervals that are required to describe all of
  // the variable domain sets.
  void getLocalIntervals(vector<interval> &pvec,
                         const map<variable,entitySet> &vm) {
    pvec.clear() ;
    vector<int> vals ;
    set<entitySet>::iterator si ;
    map<variable,entitySet>::const_iterator mi ;
    entitySet totSet ;
    for(mi=vm.begin();mi!=vm.end();++mi) {
      entitySet s = mi->second ; 
      totSet += s ;
      for(int i = 0;i < s.num_intervals(); ++i) {
	vals.push_back(s[i].first-1) ;
	vals.push_back(s[i].first) ;
	vals.push_back(s[i].second) ;
	vals.push_back(s[i].second+1) ;
      }
    }
    
    std::sort(vals.begin(),vals.end()) ;

    vector<int>::iterator uend = std::unique(vals.begin(),vals.end()) ;

    vector<int>::iterator ii = vals.begin() ;

    ++ii ;
    while(ii+1 != uend) {
      int i1 = *ii ;
      int i2 = *(ii+1)  ;
      if(i1+1 == i2)
        i2 = i1 ;
      else
        ++ii ;
      ++ii ;
      interval iv(i1,i2) ;
      if((entitySet(iv)&totSet) != EMPTY)
        pvec.push_back(iv) ;
    }
  }

  // Identify categories by the association of attributes (variables)
  // with entities.  These categories are only for entities owned by
  // the current processors (i.e. local)
  void getLocalCategories(map<variableSet, entitySet> &cmap,
                          const map<variable, entitySet> &vm) {
    cmap.clear() ;
    vector<interval> pvec ;
    getLocalIntervals(pvec,vm) ;
    vector<variableSet> cat_names(pvec.size()) ;
    map<variable,entitySet>::const_iterator mi ;
    for(mi=vm.begin();mi!=vm.end();++mi) {
      variable v = mi->first ;
      entitySet e = mi->second ;
      for(size_t i=0;i<pvec.size();++i) {
        entitySet isect = entitySet(pvec[i]) & e ;
        if(isect != EMPTY) {
          WARN(isect != pvec[i]) ;
          cat_names[i] += v ;
        }
      }
    }
    
    for(size_t i=0;i<pvec.size();++i) 
      cmap[cat_names[i]] += pvec[i] ;

  }

  

  // Unify local categories across all processors, return in a universal order
  // variables and categories.
  void getGlobalCategories(vector<variableSet> &clist,
                           vector<variable> &vlist,
                           const map<variableSet,entitySet> &cmap,
                           const map<variable,entitySet> &vm) {
    vector<string> vnlist ;
    map<variable,entitySet>::const_iterator ii ;
    for(ii=vm.begin();ii!=vm.end();++ii) {
      ostringstream oss ;
      oss << ii->first ;
      vnlist.push_back(oss.str()) ;
    }
    sort(vnlist.begin(),vnlist.end()) ;
    vlist.clear() ;
    map<variable,int> vident ;
    for(size_t i=0;i<vnlist.size();++i) {
      variable v(vnlist[i]) ;
      vlist.push_back(v) ;
      vident[v] = i ;
    }
    
    set<entitySet> cat_set ;
    map<variableSet,entitySet>::const_iterator ci ;
    for(ci=cmap.begin();ci!=cmap.end();++ci) {
      entitySet tmp ;
      variableSet s = ci->first ;
      for(variableSet::const_iterator vi=s.begin();vi!=s.end();++vi) 
        tmp += vident[*vi] ;
      cat_set.insert(tmp) ;
    }

    set<entitySet> cat_set_unified ;
    int comm_root = 0 ;
    while(comm_root != -1) {
      vector<int> buffer ;
      if(comm_root == MPI_rank) {
        set<entitySet>::const_iterator si ;
        for(si=cat_set.begin();si!=cat_set.end();++si) {
          entitySet s = *si ;
          buffer.push_back(s.size()) ;
          for(entitySet::const_iterator ei=s.begin();ei!=s.end();++ei)
            buffer.push_back(*ei) ;
        }
      }
      int buf_size = buffer.size() ;
      MPI_Bcast(&buf_size,1,MPI_INT,comm_root,MPI_COMM_WORLD) ;
      if(comm_root != MPI_rank)
        buffer = vector<int>(buf_size) ;
      MPI_Bcast(&buffer[0],buf_size,MPI_INT,comm_root,MPI_COMM_WORLD) ;

      for(size_t ex=0;ex<buffer.size();) {
        int num_ents = buffer[ex] ;
        ex++ ;
        entitySet s ;
        for(int i=0;i<num_ents;++i) {
          s += buffer[ex] ;
          ex++ ;
        }
        cat_set_unified.insert(s) ;
        cat_set.erase(s) ;
      }
      int r = -1 ;
      if(cat_set.size() != 0)
        r = MPI_rank ;
      MPI_Allreduce(&r,&comm_root,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD) ;
    }
    clist.clear() ;
    set<entitySet>::const_iterator si ;
    for(si=cat_set_unified.begin();si!=cat_set_unified.end();++si) {
      entitySet s = *si ;
      variableSet vs ;
      for(entitySet::const_iterator ei=s.begin();ei!=s.end();++ei)
        vs += vlist[*ei] ;
      clist.push_back(vs) ;
    }
  }

  // Collect the sizes of the distributed variable sets
  void getVarSizes(vector<pair<variable,int> > &var_sizes,
                   const vector<variable> &var_list,
                   const map<variable,entitySet> &vm) {
    var_sizes.clear() ;
    for(size_t i=0;i<var_list.size();++i) {
      variable v = var_list[i] ;
      map<variable,entitySet>::const_iterator vinfo = vm.find(v) ;
      int local_size = vinfo->second.size() ;
      int vsize = 0 ;
      MPI_Allreduce(&local_size,&vsize,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD) ;
      var_sizes.push_back(pair<variable,int>(v,vsize)) ;
    }
    std::stable_sort(var_sizes.begin(),var_sizes.end(),compare_var_sizes) ;
  }

  void getCategoryKeys(vector<pair<variableSet,string> > &cat_keys,
                  const vector<pair<variable,int> > &var_sizes,
                  const vector<variableSet> &cat_list) {
    cat_keys.clear() ;
    for(size_t i=0;i<cat_list.size();++i) {
      variableSet vset = cat_list[i] ;
      string key ;
      for(size_t j=0;j<var_sizes.size();++j) {
        key += vset.inSet(var_sizes[j].first)?'1':'0' ;
      }
      cat_keys.push_back(pair<variableSet,string>(vset,key)) ;
    }
    std::stable_sort(cat_keys.begin(),cat_keys.end(),compare_groups) ;
  }
                  
  void categories(fact_db &facts,vector<entitySet> &pvec) {
    map<variable, entitySet> vm ;
    getVariableAssociations(vm,facts) ;
    map<variableSet, entitySet> cmap ;
    getLocalCategories(cmap,vm) ;
    vector<variableSet> cat_list ;
    vector<variable> var_list ;
    getGlobalCategories(cat_list,var_list,cmap,vm) ;
    vector<pair<variable,int> > var_sizes ;
    getVarSizes(var_sizes,var_list,vm) ;
    vector<pair<variableSet,string> > cat_keys ;
    getCategoryKeys(cat_keys,var_sizes,cat_list) ;

    pvec.clear() ;
    debugout << "num_categories = " << cat_list.size() <<endl ;
    for(size_t i=0;i<cat_keys.size();++i) {
      variableSet vset = cat_keys[i].first ;
      entitySet eset = cmap[vset] ;
      pvec.push_back(all_collect_entitySet(eset)) ;
      debugout << "category " << cat_keys[i].first << " = " << cat_keys[i].second << ", pvec_n = " << pvec[i].num_intervals() <<  endl ;

    }
  }
}
