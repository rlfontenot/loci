#include "comp_tools.h"
#include "distribute.h"

#include <vector>
using std::vector ;
#include <map>
using std::map ;
#include <list>
using std::list ;
#include <set>
using std::set ;

//#define VERBOSE

namespace Loci {
  void impl_recurse_compiler::set_var_existence(fact_db &facts) {
    warn(facts.isDistributed()) ;

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
      debugout[MPI_rank] << "j = " << j << ", domain = " << domain
                         << ", newdomain = " << newdomain << endl ;
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
    debugout[MPI_rank] << "sdelta_init = " << sdelta << ", tdelta = " << tdelta << endl ;
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
      debugout[MPI_rank] << "sdelta = " << sdelta << ", tdelta = " << tdelta << endl ;
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
    debugout[MPI_rank] << "recursive rule " << impl << " generating " << generated << endl ;
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

    entitySet my_entities = ~EMPTY ;
    if(facts.isDistributed()) {
      fact_db::distribute_infoP d = facts.get_distribute_info() ;
      pre_send_entities = barrier_existential_rule_analysis(recurse_vars,facts) ;
      my_entities = d->my_entities ;
    }

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
      entitySet sources = my_entities ;
      entitySet constraints = my_entities ;
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
                MapRepP mp = MapRepP(facts.get_variable(*(miv[j]))->getRep()) ;
                minfo.mapvec.push_back(mp) ;
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
              MapRepP mp = MapRepP(facts.get_variable(*(miv[j]))->getRep()) ;
              minfo.mapvec.push_back(mp) ;
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

            const rule_impl::info &rinfo = fi->get_info().desc ;
            set<vmap_info>::const_iterator si ;

            cerr << "detailed report:" << endl ;
            cerr<< "fctrl.nr_sources = " << fctrl.nr_sources << endl ;
            entitySet constraints = fctrl.constraints ;
            for(si=rinfo.sources.begin();si!=rinfo.sources.end();++si) {
              entitySet sources = vmap_source_exist(*si,facts) ;
              sources &= my_entities ;
              if((sources & constraints) != constraints) {
                cerr << "sources & constraints != constraints for input"
                     << endl
                     << sources  << " -- " << *si << endl ;
                
                if(si->mapping.size() > 0) {
                  entitySet working = constraints ;
                  for(int i=0;i<si->mapping.size();++i) {
                    entitySet images ;
                    variableSet::const_iterator vi ;
                    for(vi=si->mapping[i].begin();vi!=si->mapping[i].end();++vi)
                      images |= facts.image(*vi,working) ;
                    working = images ;
                  }
                  variableSet::const_iterator vi ;
                  for(vi=si->var.begin();vi!=si->var.end();++vi) {
                    entitySet exist = facts.variable_existence(*vi) ;
                    entitySet fails = working & ~exist ;
                    if(fails != EMPTY) {
                      cerr << "expecting to find variable " << *vi
                           << " at entities " << fails << endl
                           << *vi << " exists at entities " << exist << endl ;
                    }
                  }
                }
              }
            }
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
        debugout[MPI_rank] << "recursive rule " << *fi << " generating "
                           << srcs << endl ;
#endif
      }

      recurse_send_entities.push_back(vector<pair<variable,entitySet> >()) ;

      int deltas = 0 ;
      finished = true ;
      for(vi=tvars.begin();vi!=tvars.end();++vi) {
        entitySet tvar_delta = tvar_update[*vi] - tvar_computed[*vi] ;

        entitySet tvar_other = tvar_delta - my_entities ;
        recurse_send_entities.back().push_back(make_pair(*vi,tvar_other)) ;

        entitySet recvset = send_entitySet(tvar_other,facts) ;
        tvar_delta += recvset ;
        tvar_delta += fill_entitySet(tvar_delta,facts) ;
        if(tvar_delta!=EMPTY)
          deltas++ ;
        tvar_update[*vi] += tvar_delta ;
        tvar_computed[*vi] = tvar_update[*vi] ;
      }
      int deltar = 0 ;
      if(!facts.isDistributed())
        deltar = deltas ;
      else
        MPI_Allreduce(&deltas,&deltar, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD) ;
      if(deltar != 0)
        finished = false ;
    }

    for(fi=fset.begin();fi!=fset.end();++fi) {
      fcontrol &fctrl = control_set[*fi] ;
      for(map<variable,entitySet>::const_iterator mi=fctrl.generated.begin();
          mi!=fctrl.generated.end();++mi) {

        entitySet create = mi->second ;
        create += send_entitySet(create,facts) ;
        create += fill_entitySet(create,facts) ;
        
        facts.set_existential_info(mi->first,*fi,create) ;
      }
    }
  }

  void recurse_compiler::process_var_requests(fact_db &facts) {
    entitySet my_entities = ~EMPTY ;
    if(facts.isDistributed()) {
      fact_db::distribute_infoP d = facts.get_distribute_info() ;
      my_entities = d->my_entities ;
    }

    ruleSet::const_iterator fi ;
    if(facts.isDistributed()) {
      list<comm_info> request_comm ;
      request_comm = barrier_process_rule_requests(recurse_vars, facts) ;

      vector<pair<variable,entitySet> >::const_iterator vi ;
      vector<pair<variable,entitySet> > send_requested ;

      for(vi=pre_send_entities.begin();vi!=pre_send_entities.end();++vi) {
        variable v = vi->first ;
        entitySet send_set = vi->second ;
        send_requested.push_back(make_pair(v,send_set &
                                           facts.get_variable_requests(v))) ;
      }
      pre_plist = put_precomm_info(send_requested, facts) ;
      pre_clist = sort_comm(request_comm,facts) ;
    }

    if(facts.isDistributed()) {
      for(fi=recurse_rules.begin();fi!=recurse_rules.end();++fi) {
        fcontrol &fctrl = control_set[*fi] ;
        entitySet control = process_rule_requests(*fi,facts) ;
        list<entitySet>::iterator ci ;
        entitySet total ;
        for(ci=fctrl.control_list.begin();ci!=fctrl.control_list.end();++ci) {
          *ci -= total ;
          total += *ci ;
        }
      }
      map<rule, list<entitySet>::reverse_iterator> rpos ;
      ruleSet::const_iterator ri ;

      map<variable,entitySet> var_requests, recurse_entities ;
      for(variableSet::const_iterator vi=recurse_vars.begin() ;
          vi!=recurse_vars.end();
          ++vi) {
        var_requests[*vi] = facts.get_variable_requests(*vi) ;
        for(ri=recurse_rules.begin();ri!=recurse_rules.end();++ri) {
          recurse_entities[*vi] += facts.get_existential_info(*vi,*ri) ;
        }
      }
      
                       
      for(ri=recurse_rules.begin();ri!=recurse_rules.end();++ri)
        rpos[*ri] = control_set[*ri].control_list.rbegin() ;

      map<variable,vector<entitySet> > recurse_send_req ;
      map<variable,entitySet> all_requests ;
      bool finished = false ;
      do {
        map<variable,entitySet> vreq_map ;
        for(ri=recurse_rules.begin();ri!=recurse_rules.end();++ri) {
          entitySet &context = *rpos[*ri] ;

          fcontrol &fctrl = control_set[*ri] ;
          for(int i=0;i<fctrl.target_maps.size();++i) {
            entitySet ct = var_requests[fctrl.target_maps[i].v] ;
            for(int j=fctrl.target_maps[i].mapvec.size()-1;j>=0;--j)
              ct = fctrl.target_maps[i].mapvec[j]->preimage(ct).first ;
            context &= ct ;
          }
          for(int i=0;i<fctrl.recursion_maps.size();++i) {
            entitySet rq = context ;
            for(int j=0;j<fctrl.recursion_maps[i].mapvec.size();j++)
              rq = fctrl.recursion_maps[i].mapvec[j]->image(rq) ;
            vreq_map[fctrl.recursion_maps[i].v] += rq  ;
           }
          
          
          rpos[*ri]++ ;
          if(rpos[*ri] == control_set[*ri].control_list.rend())
            finished = true ;

          for(variableSet::const_iterator vi = recurse_vars.begin();
              vi!=recurse_vars.end();
              ++vi) {
            all_requests[*vi] += vreq_map[*vi] ;
            entitySet remain = vreq_map[*vi] & recurse_entities[*vi] ;
            remain -= my_entities ;
            recurse_send_req[*vi].push_back(remain) ;
          }
       }
      } while(!finished) ;

      for(variableSet::const_iterator vi = recurse_vars.begin();
          vi!=recurse_vars.end();
          ++vi) {
        std::reverse(recurse_send_req[*vi].begin(),recurse_send_req[*vi].end()) ;
        entitySet req_loc ;
        for(vector<entitySet>::iterator vei=recurse_send_req[*vi].begin();
            vei!=recurse_send_req[*vi].end();
            ++vei) {
          *vei -= req_loc ;
          list<comm_info> req_comm ;
          send_requests(*vei,*vi,facts,req_comm) ;
          send_req_var[*vi].push_back(sort_comm(req_comm,facts)) ;
          req_loc += *vei ;
        }

      }
    } else {
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
    
  }

  executeP recurse_compiler::create_execution_schedule(fact_db &facts) {
    CPTR<execute_sequence> el = new execute_sequence ;
    if(facts.isDistributed()) {
      el->append_list(new execute_thread_sync) ;
      el->append_list(new execute_comm(pre_plist, facts) ) ; 
      el->append_list(new execute_comm(pre_clist, facts)) ;
    }
    
    map<rule, list<entitySet>::const_iterator> rpos ;
    list<vector<pair<variable,entitySet> > >::const_iterator
      sei = recurse_send_entities.begin() ;
    ruleSet::const_iterator ri ;

    int idx = 0 ;
    for(ri=recurse_rules.begin();ri!=recurse_rules.end();++ri)
      rpos[*ri] = control_set[*ri].control_list.begin() ;

    bool finished = false ;
    do {
      for(variableSet::const_iterator vi=recurse_vars.begin();
          vi!=recurse_vars.end();
          ++vi) {
        vector<list<comm_info> > &commv = send_req_var[*vi] ;
        if(idx<commv.size() && commv[idx].size() != 0)
          el->append_list(new execute_comm(commv[idx],facts)) ;
      }
      idx++ ;
      for(ri=recurse_rules.begin();ri!=recurse_rules.end();++ri) {
        const fcontrol &fctrl = control_set[*ri] ;
        list<entitySet>::const_iterator &li = rpos[*ri] ;
        
        if(li==fctrl.control_list.end()) 
          finished = true ;
        else {
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
      if(!finished) {
        if(num_threads > 1)
          el->append_list(new execute_thread_sync) ;
        if(facts.isDistributed()) {
          list<comm_info> plist = put_precomm_info(*sei, facts) ;
          el->append_list(new execute_comm(plist,facts)) ;
        }
        sei++ ;
      }
    } while(!finished) ;

    if(facts.isDistributed()) 
      el->append_list(new execute_comm(pre_clist, facts)) ;
    
    if(el->size() == 0)
      return 0 ;
    else
      return executeP(el) ;
  }

  
}
