#include "comp_tools.h"

namespace Loci {
  void impl_recurse_compiler::set_var_existence(fact_db &facts) {
    variableSet::const_iterator vi ;
    variableSet tvars ;
    tvars = impl.targets() ;
    fcontrol &fctrl = control_set ;
    const rule_impl::info &finfo = impl.get_info().desc ;
    warn(impl.type() == rule::INTERNAL) ;
    entitySet sources = ~EMPTY ;

    entitySet constraints = ~EMPTY ;
    set<vmap_info>::const_iterator si ;
    for(si=finfo.sources.begin();si!=finfo.sources.end();++si) {
      if((si->var & tvars) == EMPTY)
        sources &= vmap_source_exist(*si,facts) ;
      else {
        const int num_maps = si->mapping.size() ;
        if(num_maps != 0) {
          for(vi=si->mapping[0].begin();vi!=si->mapping[0].end();++vi) {
            sources &= facts.variable_existence(*vi) ;
          }
        } 
        
        tvars += si->var ;
        vector<variableSet::const_iterator> miv(num_maps) ;
        for(int i=0;i<num_maps;++i)
          miv[i] = si->mapping[i].begin() ;
        
        for(vi=si->var.begin();vi!=si->var.end();++vi) {
          int i = 0 ;
          do {
            fcontrol::mapping_info minfo ;
            for(int j=0;j!=num_maps;++j) {
              FATAL(!facts.is_a_Map(*(miv[j]))) ;
              MapRepP m = MapRepP(facts.get_variable(*(miv[j]))->getRep()) ;  
              minfo.mapvec.push_back(m) ;
              minfo.mapvar.push_back(*(miv[j])) ;
            }
            minfo.v = *vi ;
            
            for(i=0;(i!=num_maps) &&(++(miv[i]) == si->mapping[i].end());++i) {
              miv[i] = si->mapping[i].begin() ;
            }
            fctrl.recursion_maps.push_back(minfo) ;
          } while(i!=num_maps) ;
        }
      }
    } 
    for(si=finfo.targets.begin();si!=finfo.targets.end();++si) {
      int num_maps = si->mapping.size() ;
      vector<variableSet::const_iterator> miv(num_maps) ;
      for(int i=0;i<num_maps;++i)
        miv[i] = si->mapping[i].begin() ;
      
      for(vi=si->var.begin();vi!=si->var.end();++vi) {
        int i = 0 ;
        do {
          fcontrol::mapping_info minfo ;
          for(int j=0;j!=num_maps;++j) {
            FATAL(!facts.is_a_Map(*(miv[j]))) ;
            MapRepP m = MapRepP(facts.get_variable(*(miv[j]))->getRep()) ;
            minfo.mapvec.push_back(m) ;
            minfo.mapvar.push_back(*(miv[j])) ;
          }
          minfo.v = *vi ;
          
          for(i=0;(i!=num_maps) &&(++(miv[i]) == si->mapping[i].end());++i) {
            miv[i] = si->mapping[i].begin() ;
          }
          fctrl.target_maps.push_back(minfo) ;
        } while(i!=num_maps) ;
      }
    }
    fctrl.use_constraints = false ;
    for(si=finfo.constraints.begin();si!=finfo.constraints.end();++si) {
      warn((si->var & tvars) != EMPTY) ;
      constraints &= vmap_source_exist(*si,facts) ;
      fctrl.use_constraints = true ;
    }

    fctrl.nr_sources = sources ;

    fctrl.constraints = constraints ;
    warn(fctrl.recursion_maps.size() == 0) ;

    fatal(tvars.size() != 1 ||
          control_set.recursion_maps.size() != 1 ||
          control_set.target_maps.size() != 1 ||
          fctrl.use_constraints) ;
    
    fcontrol::mapping_info &rmap = fctrl.recursion_maps[0] ;
    fcontrol::mapping_info &tmap = fctrl.target_maps[0] ;
    vector<multiMap> read_maps,read_map_inv ;

    for(int j=0;j<rmap.mapvec.size();++j) {
      read_maps.push_back(rmap.mapvec[j]->get_map()) ;
    }
    for(int j=0;j<read_maps.size();++j) {
      read_map_inv.push_back(multiMap()) ;
    }
    variable rvar = *(tvars.begin()) ;
    entitySet sdelta = facts.variable_existence(rvar) ;
    entitySet domain = fctrl.nr_sources + sdelta ;

    for(int j=read_maps.size()-1;j>=0;--j) {
      entitySet newdomain = facts.preimage(rmap.mapvar[j],domain).first ;
#ifdef VERBOSE
      cout << "j = " << j << ", domain = " << domain << ", newdomain = "
           << newdomain << endl ;
#endif
      if(domain == ~EMPTY) {
        cerr << "problem in recurse compiler for rule = "<< impl << endl ;
        domain = read_maps[j].domain() ;
      }
      inverseMap(read_map_inv[j],read_maps[j],domain,newdomain) ;
      domain = newdomain ;
    }
  
    for(int j=rmap.mapvec.size()-1;j>=0;--j)
      sdelta = rmap.mapvec[j]->preimage(sdelta).first ;
    sdelta &= fctrl.nr_sources ;
    entitySet tdelta = sdelta ;
    for(int j=0;j<tmap.mapvec.size();++j)
      tdelta = tmap.mapvec[j]->image(tdelta) ;
#ifdef VERBOSE
    cout << "sdelta_init = " << sdelta << ", tdelta = " << tdelta << endl ;
#endif

    if(num_threads > 1)
      par_schedule.push_back(sdelta) ;
    else
      fastseq += sequence(sdelta) ;
    
    entitySet generated = tdelta ;
    const entitySet nr_sources = fctrl.nr_sources ;
    store<bool> exists ;
    exists.allocate(nr_sources) ;
    for(entitySet::const_iterator
          ei=nr_sources.begin();ei!=nr_sources.end();++ei)
      exists[*ei] = false ;
    for(entitySet::const_iterator
          ei=tdelta.begin();ei!=tdelta.end();++ei)
      exists[*ei] = true ;
    vector<const int *> mi(read_maps.size()) ;
    vector<const int *> me(read_maps.size()) ;

    do {
      for(int j=read_map_inv.size()-1;j>=0;--j) {
        entitySet candidates ;
        const int *mi ;
        for(entitySet::const_iterator di=sdelta.begin();di!=sdelta.end();++di) {
#ifdef DEBUG
          if(!read_map_inv[j].domain().inSet(*di)) {
            cerr << " di = " << *di << ", j = " << j << endl ;
          } else
#endif
            //        fatal(!read_map_inv[j].domain().inSet(*di)) ;
            for(mi=read_map_inv[j].begin(*di); mi!=read_map_inv[j].end(*di);++mi) 
              candidates += *mi ;
        }
        sdelta = candidates ;
      }
      // we now have candidate sdeltas, check them to see if they
      // are satisfied.
      entitySet satisfied ;
      for(entitySet::const_iterator di=sdelta.begin();di!=sdelta.end();++di) {
        int c = *di ;
        mi[0] = read_maps[0].begin(c) ;
        me[0] = read_maps[0].end(c) ;
        int j=0 ;
        const int je=read_maps.size();
        bool chk = true ;
        // Check that all paths from candidate go to generated cells
        while(chk && j>=0) {
          c = *mi[j] ;
          j++ ;
          if(j==je) {
            chk = chk && exists[c] && nr_sources.inSet(c) ;
            do {
              j-- ;
            } while(j>=0 && (++mi[j]) == me[j]) ;
          } else {
            mi[j] = read_maps[j].begin(c) ;
            me[j] = read_maps[j].end(c) ;
          }
        
        }
        if(chk)
          satisfied += *di ;
      }
      sdelta = satisfied ;
      entitySet tdelta = sdelta ;
      for(int j=0;j<tmap.mapvec.size();++j)
        tdelta = tmap.mapvec[j]->image(tdelta) ;
#ifdef VERBOSE
      cout << "sdelta = " << sdelta << ", tdelta = " << tdelta << endl ;
#endif
      if(num_threads>1)
        par_schedule.push_back(sdelta) ;
      else
        fastseq += sequence(sdelta) ;

      for(entitySet::const_iterator
            ei=tdelta.begin();ei!=tdelta.end();++ei) 
        exists[*ei] = true ;
    
    } while(sdelta != EMPTY) ;

    for(entitySet::const_iterator
          ei=nr_sources.begin();ei!=nr_sources.end();++ei)
      if(exists[*ei]) {
        const int start = *ei ;
        const int end = nr_sources.Max() ;
        int finish = start ;
        for(;*ei!=end && exists[*ei];++ei)
          finish = *ei ;
        if(*ei == end && exists[end])
          finish = end ;
        generated += interval(start,finish) ;
      }
  
  
    fctrl.generated[rvar] = generated ;
#ifdef VERBOSE
    cout << "recursive rule " << impl << " generating " << generated << endl ;
#endif
    
    for(map<variable,entitySet>::const_iterator mi=fctrl.generated.begin();
        mi!=fctrl.generated.end();++mi)
      facts.set_existential_info(mi->first,impl,mi->second) ;
    
  }
  
  void impl_recurse_compiler::process_var_requests(fact_db &facts) {
    process_rule_requests(impl,facts) ;
  }

  executeP impl_recurse_compiler::create_execution_schedule(fact_db &facts) {

    if(num_threads > 1) {
      CPTR<execute_list> el = new execute_list ;
      sequence seq ;
      for(int i=0;i<par_schedule.size();++i) {
        if(par_schedule[i].size() < num_threads*4) {
          seq += sequence(par_schedule[i]) ;
        } else {
          if(seq.size() > 1) {
            el->append_list(new execute_rule(impl,seq,facts)) ;
            el->append_list(new execute_thread_sync) ;
            seq = sequence() ;
          }
          execute_par *ep = new execute_par ;
          parallel_schedule(ep,par_schedule[i],impl,facts) ;
          el->append_list(ep) ;
          el->append_list(new execute_thread_sync) ;
        }
      }
      if(seq.size() > 1) {
        el->append_list(new execute_rule(impl,seq,facts)) ;
      }
      return executeP(el) ;
    }

    return new execute_rule(impl,fastseq,facts) ;

  }

  void recurse_compiler::set_var_existence(fact_db &facts) {
    control_set.clear() ;
  
    ruleSet::const_iterator fi ;
    variableSet::const_iterator vi ;
    variableSet tvars ;
    const ruleSet &fset = recurse_rules ;
    for(fi=fset.begin();fi!=fset.end();++fi) 
      tvars += fi->targets() ;
    for(fi=fset.begin();fi!=fset.end();++fi) {
      fcontrol &fctrl = control_set[*fi] ;
      const rule_impl::info &finfo = fi->get_info().desc ;
      warn(fi->type() == rule::INTERNAL) ;
      entitySet sources = ~EMPTY ;
      entitySet constraints = ~EMPTY ;
      set<vmap_info>::const_iterator si ;
      for(si=finfo.sources.begin();si!=finfo.sources.end();++si) {
        if((si->var & tvars) == EMPTY)
          sources &= vmap_source_exist(*si,facts) ;
        else {
          tvars += si->var ;
          int num_maps = si->mapping.size() ;
          vector<variableSet::const_iterator> miv(num_maps) ;
          for(int i=0;i<num_maps;++i)
            miv[i] = si->mapping[i].begin() ;
        
          for(vi=si->var.begin();vi!=si->var.end();++vi) {
            int i = 0 ;
            do {
              fcontrol::mapping_info minfo ;
              for(int j=0;j!=num_maps;++j) {
                FATAL(!facts.is_a_Map(*(miv[j]))) ;
                minfo.mapvec.push_back(
                                       MapRepP(facts.get_variable(*(miv[j]))->getRep())) ;
              }
              minfo.v = *vi ;
              
              for(i=0;(i!=num_maps)&&(++(miv[i]) == si->mapping[i].end());++i)
                miv[i] = si->mapping[i].begin() ;

              fctrl.recursion_maps.push_back(minfo) ;
            } while(i!=num_maps) ;
          }
        }
      } 
      for(si=finfo.targets.begin();si!=finfo.targets.end();++si) {
        int num_maps = si->mapping.size() ;
        vector<variableSet::const_iterator> miv(num_maps) ;
        for(int i=0;i<num_maps;++i)
          miv[i] = si->mapping[i].begin() ;

        for(vi=si->var.begin();vi!=si->var.end();++vi) {
          int i = 0 ;
          do {
            fcontrol::mapping_info minfo ;
            for(int j=0;j!=num_maps;++j) {
              FATAL(!facts.is_a_Map(*(miv[j]))) ;
              minfo.mapvec.push_back(
                                     MapRepP(facts.get_variable(*(miv[j]))->getRep())) ;
            }
            minfo.v = *vi ;
          
            for(i=0;(i!=num_maps) &&(++(miv[i]) == si->mapping[i].end());++i) {
              miv[i] = si->mapping[i].begin() ;
            }
            fctrl.target_maps.push_back(minfo) ;
          } while(i!=num_maps) ;
        }
      }
      fctrl.use_constraints = false ;
      for(si=finfo.constraints.begin();si!=finfo.constraints.end();++si) {
        warn((si->var & tvars) != EMPTY) ;
        constraints &= vmap_source_exist(*si,facts) ;
        fctrl.use_constraints = true ;
      }
      fctrl.nr_sources = sources ;
      fctrl.constraints = constraints ;
      warn(fctrl.recursion_maps.size() == 0) ;
    }
  

    map<variable,entitySet> tvar_computed,tvar_update ;
    for(vi=tvars.begin();vi!=tvars.end();++vi) {
      tvar_computed[*vi] = facts.variable_existence(*vi) ;
      tvar_update[*vi] = tvar_computed[*vi] ;
    }
    map<variable,entitySet>::iterator mi ;
    for(mi=tvar_update.begin();mi!=tvar_update.end();++mi)
      mi->second = EMPTY ;

    bool finished = false;

    while(!finished) {
      for(fi=fset.begin();fi!=fset.end();++fi) {
        fcontrol &fctrl = control_set[*fi] ;
        // partial vmap_souce_exist
        entitySet srcs = fctrl.nr_sources ;
        for(int i=0;i<fctrl.recursion_maps.size();++i) {
          entitySet sdelta = tvar_computed[fctrl.recursion_maps[i].v] ;
          for(int j=fctrl.recursion_maps[i].mapvec.size()-1;j>=0;--j)
            sdelta = fctrl.recursion_maps[i].mapvec[j]->preimage(sdelta).first ;
          fctrl.recursion_maps[i].total += sdelta ;
          srcs &= fctrl.recursion_maps[i].total ;
        }
        if(fctrl.use_constraints) {
          if((srcs & fctrl.constraints) != fctrl.constraints) {
            cerr << "recursive rule: " << *fi
                 << " cannot supply all entitites in constraint" << endl ;
            cerr << "constraints = " << fctrl.constraints << endl ;
            entitySet sac = srcs & fctrl.constraints ;
            cerr << "srcs & constraints = " << sac << endl ;
            //exit(-1) ;
          }
          srcs &= fctrl.constraints ;
        }

        fctrl.control_list.push_back(srcs) ;

        for(int i=0;i<fctrl.target_maps.size();++i) {
          entitySet trgts = srcs ;
          for(int j=0;j<fctrl.target_maps[i].mapvec.size();++j)
            trgts = fctrl.target_maps[i].mapvec[j]->image(trgts) ;
          tvar_update[fctrl.target_maps[i].v] += trgts ;
          fctrl.generated[fctrl.target_maps[i].v] += trgts ;
        }
#ifdef VERBOSE
        cout << "recursive rule " << *fi << " generating " << srcs << endl ;
#endif
      }
      finished = true ;
      for(vi=tvars.begin();vi!=tvars.end();++vi) {
        entitySet tvar_delta = tvar_update[*vi] - tvar_computed[*vi] ;
        if(tvar_delta!=EMPTY)
          finished = false ;
        tvar_computed[*vi] = tvar_update[*vi] ;
      }
    }

    for(fi=fset.begin();fi!=fset.end();++fi) {
      fcontrol &fctrl = control_set[*fi] ;
      for(map<variable,entitySet>::const_iterator mi=fctrl.generated.begin();
          mi!=fctrl.generated.end();++mi)
        facts.set_existential_info(mi->first,*fi,mi->second) ;
    }
  }

  void recurse_compiler::process_var_requests(fact_db &facts) {
    ruleSet::const_iterator fi ;
    for(fi=recurse_rules.begin();fi!=recurse_rules.end();++fi) {
      fcontrol &fctrl = control_set[*fi] ;
      entitySet control = process_rule_requests(*fi,facts) ;
      list<entitySet>::iterator ci ;
      entitySet total ;
      for(ci=fctrl.control_list.begin();ci!=fctrl.control_list.end();++ci) {
        *ci -= total ;
        total += *ci ;
      }
      do {
        if(fctrl.control_list.back() == EMPTY)
          fctrl.control_list.pop_back() ;
        if(!fctrl.control_list.empty())
          fctrl.control_list.back() &= control ;
      } while(!fctrl.control_list.empty() && fctrl.control_list.back() == EMPTY) ;
    }
  }

  executeP recurse_compiler::create_execution_schedule(fact_db &facts) {
    CPTR<execute_list> el = new execute_list ;
    
    map<rule, list<entitySet>::const_iterator> rpos ;
    ruleSet::const_iterator ri ;
    
    for(ri=recurse_rules.begin();ri!=recurse_rules.end();++ri)
      rpos[*ri] = control_set[*ri].control_list.begin() ;

    bool finished ;
    do {
      finished = true ;
      
      for(ri=recurse_rules.begin();ri!=recurse_rules.end();++ri) {
        const fcontrol &fctrl = control_set[*ri] ;
        list<entitySet>::const_iterator &li = rpos[*ri] ;
        if(li != fctrl.control_list.end()) {
          finished = false ;
          if(li->size() != 0) {
            const entitySet &exec_seq = *li ;

            if(num_threads > 1 && exec_seq.size() > 1 &&
               (*ri).get_info().rule_impl->thread_rule()) {
              execute_par *ep = new execute_par ;
              parallel_schedule(ep,*li,*ri,facts) ;
              el->append_list(ep) ;
            } else {
              el->append_list(new execute_rule(*ri,sequence(*li),facts)) ;
            }
          }
          li++ ;
        }
      }
      if(num_threads > 1)
        el->append_list(new execute_thread_sync) ;
    } while(!finished) ;
    
    if(el->size() == 0)
      return 0 ;
    else
      return executeP(el) ;
  }

  
}
