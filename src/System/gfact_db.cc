//#############################################################################
//#
//# Copyright 2008, Mississippi State University
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
#include <typeinfo>
#include <Tools/except.h>

#include <gfact_db.h>
//#include <rule.h>
#include <gconstraint.h>
#include <Tools/debugger.h>
#include <gstore.h>
#include "dist_tools.h"

#include <Loci_Datatypes.h>

using std::string ; 
using std::map ;
using std::make_pair ;
using std::vector ;
using std::list ;
using std::sort ;
using std::pair ;


using std::istream ;
using std::ostream ;
using std::ostringstream ;
using std::endl ;
using std::ios ;

#include <Tools/parse.h>

namespace Loci {
  extern int MPI_processes ;
  extern int MPI_rank ;
  ////////////////////////////////////////////////////////////////////////////////////////
  //the following part is temporary
  ////////////////////////////////////////////////////////////////////////////////////////
  extern gfact_db *exec_current_fact_db;
  int factdb_allocated_base = 0 ;
  //////////////////////////////////////////////////////////////////////////////////////
    
  gfact_db::gfact_db() {
    //set default key manager as distKeyManager 
    gkey_manager = CPTR<gKeyManager>(new distKeyManager());
    /////////////////////////////////////////////////////////////////////////////////////////
    //the following part is temporary
    /////////////////////////////////////////////////////////////////////////////////////
    distributed_info = 0 ;
    maximum_allocated = factdb_allocated_base+0 ;
    minimum_allocated = factdb_allocated_base-1 ;
    for(int i = 0; i < MPI_processes; ++i) {
      init_ptn.push_back(EMPTY) ;
    }
    //    exec_current_fact_db = this;
    dist_from_start = 0 ;
    //why need EMPTY and UNIVERSE?
    gConstraint EMPTY_constraint ;
    EMPTY_constraint = EMPTY ;
    create_gfact("EMPTY",EMPTY_constraint) ;
    gConstraint UNIVERSE_constraint ;
    UNIVERSE_constraint = ~EMPTY ;
    create_gfact("UNIVERSE",UNIVERSE_constraint) ;
    
  }

  gfact_db::~gfact_db() {}

  //this function is not tested
  void
  gfact_db::copy_all_from(const gfact_db& f) {
    synonyms = f.synonyms ;
    gfmap = f.gfmap ;
    fmap = f.fmap ;
    tmap = f.tmap ;
    nspace_vec = f.nspace_vec ;
    extensional_facts = f.extensional_facts ;
    // deep copy the key manager
    if(f.gkey_manager != 0) {
      gkey_manager = f.gkey_manager->clone() ;
    }
    /////////////////////////////////////////////////////////////////////////////////////////
    //the following part is temporary
    /////////////////////////////////////////////////////////////////////////////////////
    init_ptn = f.init_ptn ;
    global_comp_entities = f.global_comp_entities ;
    maximum_allocated = f.maximum_allocated ;
    minimum_allocated = f.minimum_allocated ;
    dist_from_start = f.dist_from_start ;

    // create new keyspaces
    keyspace.clear() ;
    const std::map<std::string,KeySpaceP>& fks = f.keyspace ;
    for(std::map<std::string,KeySpaceP>::const_iterator
          mi=fks.begin();mi!=fks.end();++mi) {
      keyspace[mi->first] = mi->second->new_keyspace() ;
    }      
    for(std::map<std::string,KeySpaceP>::iterator
          mi=keyspace.begin();mi!=keyspace.end();++mi) {
      (mi->second)->set_synonyms(&synonyms) ;
    }
    // deep copy the key manager
    if(f.key_manager != 0) {
      key_manager = f.key_manager->clone() ;
    }
    /* we cannot use the following direct assignment
       to copy the distributed_info from f since
       distributed_info is a NPTR pointer and is
       reference counted this would not be a true copy
    */
    // distributed_info = f.distributed_info ;
    distribute_infoP& df = distributed_info ;
    const distribute_infoP& fdf = f.distributed_info ;
    if(fdf == 0) {
      df = 0 ;
      return ;
    }
    df = new distribute_info ;
    df->myid = fdf->myid ;
    df->isDistributed = fdf->isDistributed ;
    // we make a deep copy of the maps
    entitySet l2g_alloc = fdf->l2g.domain() ;
    entitySet g2l_alloc = fdf->g2l.domain() ;
    df->l2g.allocate(l2g_alloc) ;
    df->g2l.allocate(g2l_alloc) ;
    for(entitySet::const_iterator ei=l2g_alloc.begin();
        ei!=l2g_alloc.end();++ei) {
      df->l2g[*ei] = fdf->l2g[*ei] ;
    }
    for(entitySet::const_iterator ei=g2l_alloc.begin();
        ei!=g2l_alloc.end();++ei)
      df->g2l[*ei] = fdf->g2l[*ei] ;
    df->my_entities = fdf->my_entities ;
    df->comp_entities = fdf->comp_entities ;
    df->copy = fdf->copy ;
    df->xmit = fdf->xmit ;
    df->copy_total_size = fdf->copy_total_size ;
    df->xmit_total_size = fdf->xmit_total_size ;
    //      df->remap = fdf->remap ;
    df->g2f = fdf->g2f ; 
  }

  // void gfact_db::make_extensional_fact(const variable& v) {
  //   // if it is already extensional_fact, then we do nothing
  //   if(extensional_facts.inSet(v))
  //     return ;
  //   std::map<variable, gStoreRepP>::iterator mi = gfmap.find(v) ;
  //   // the fact must exist in the fact_db, otherwise we do nothing
  //   if(mi != gfmap.end())
  //     extensional_facts += v ;
  // }
  // void gfact_db::make_intensional_fact(const variable& v) {
  //   // We perform action only if it is already an extensional_fact
  //   // which implies that the fact exists in the fact_db
  //   if(extensional_facts.inSet(v))
  //     extensional_facts -= v ;
  // }
  
  void gfact_db::set_variable_type(variable v, storeRepP st) {
    v = add_namespace(v) ;
    tmap[v] = storeRepP(st->new_store(EMPTY)) ;
  }

  storeRepP gfact_db::get_variable_type(variable v) const {
    v = add_namespace(v) ;
    map<variable,storeRepP>::const_iterator mi ;
    if((mi=tmap.find(v)) != tmap.end())
      return storeRepP(mi->second->new_store(EMPTY)) ;
    else
      return storeRepP(0) ;
  }
  
  variable gfact_db::remove_synonym(variable v) const {
    std::map<variable,variable>::const_iterator mi ;
    while((mi=synonyms.find(v)) != synonyms.end())
      v = mi->second ;
    return v ;
  }
  
  variable gfact_db::get_real_var(variable v) const {
    //  variable v = add_namespace(variable(vname)) ;
    variable real_v = v;
    if(synonyms.find(v) != synonyms.end()) {
      real_v = remove_synonym(v) ;
    }
    return real_v;
  }
    

  //assume when synonym_variable happens, v is in fmap
  //i.e. local numbering is used ???
  void gfact_db::synonym_variable(variable v, variable synonym) {
    // Find all variables that should be synonymous with v
    variableSet synonym_set ;
    std::map<variable,variable>::const_iterator mi ;
    while((mi=synonyms.find(v)) != synonyms.end()) {
      synonym_set += v ;
      v = mi->second ;
    }
    variable s = synonym;
    while((mi=synonyms.find(s)) != synonyms.end()) {
      synonym_set += s ;
      s = mi->second ;
    }
    synonym_set += s ;
 
    // If the two are already synonymous, we are done
    if(s == v)
      return ;

    // Sanity check, make sure v exists
    std::map<variable,fact_info>::iterator vmi, vmj ;
    if((vmi = fmap.find(v)) == fmap.end()) {
      cerr << "WARNING: synonym_variable("<<v<<","<<synonym<<")"<<endl ;
      cerr << "WARNING: type not known for target of synonym, ignoring operation" << endl ;
      //      abort() ;
      return ;
    }
    // If the synonym already points to a different variable instance,
    // remove it
    if((vmj = fmap.find(s)) != fmap.end()) {
      fact_info &finfo = vmj->second ;
      if((finfo.data_rep->domain() != EMPTY &&
          finfo.data_rep->RepType() != PARAMETER)) {
        cerr << "unable to define synonym variable " << synonym
             << " when variable already created in db. "  << endl ;
        cerr << "variable v = " << v << endl ;
        abort() ;
      }
      remove_variable(synonym) ;
    }
    
    // Add new synonyms so that they point to v
    for(variableSet::const_iterator vi = synonym_set.begin();
        vi!=synonym_set.end();
        ++vi) {
      synonyms[*vi] = v ;
    }
    
  }
  
  void gfact_db::copy_facts(){
    fmap.clear();
    for(std::map<variable,gStoreRepP>::const_iterator itr = gfmap.begin(); itr != gfmap.end(); itr++){
      variable v = itr->first;
      create_fact(v, itr->second->copy2store());
    }
    
  }
 
 
  variable gfact_db::add_namespace(variable v) const {
    variable tmp_v ;
    if((nspace_vec.size() > 0) && (v.get_info().namespac.size() == 0)) {
      tmp_v = v ;
      for(size_t i = 0; i < nspace_vec.size(); ++i)
        tmp_v = tmp_v.add_namespace(nspace_vec[i]) ;
    } else
      tmp_v = v ;
    return tmp_v ;
  }
  // this is the basic method that creats a fact
  // in the gfact_db. It is served as the basis
  // for create_fact methods
  //this method adds st to gfmap with duplication check
  //this method will not process keyspace info
  void gfact_db::create_pure_gfact(const variable& v, gStoreRepP st) {
    //get real variable
    variable real_v = get_real_var(v);
    std::map<variable, gStoreRepP>::iterator mi = gfmap.find(real_v) ;
    if(mi==gfmap.end()) {
      gfmap[real_v]= st ;
    } else {
      if(typeid(st->getRep()) != typeid(mi->second->getRep())) {
        cerr << "set_variable_type() method of gfact_db changing type for variable " << v << endl ;
      }
      mi->second = st ;
    }
    return ;
  }

  void gfact_db::set_variable_domain_space(const variable& v, gStoreRepP st, gKeySpaceP space){
    if(space != gKeySpaceP(0)){
      //get real variable   
      variable real_v = get_real_var(v) ;
      st->set_domain_space(&(*space));
      space->add_out_var(real_v);
    }
    return;
  }
  void gfact_db::set_variable_image_space(const variable& v, gMapRepP st, gKeySpaceP space){
    if(space != 0){
      //get real variable   
      variable real_v = get_real_var(v) ;
      st->set_image_space(&(*space));
      space->add_in_var(real_v); 
      return;
    }
  }
  void gfact_db::remove_variable_domain_space(const variable& v){
     
    gStoreRepP st = get_gvariable(v);
    gKeySpaceP space = st->get_domain_space();
    if(space != 0 ){
      variable real_v = get_real_var(v);
      space->remove_out_var(real_v);
      st->set_domain_space(0);
    }
  }
  
  void gfact_db::remove_variable_image_space(const variable& v){
    //get real varaible 
    gStoreRepP st = get_gvariable(v);
    gKeySpaceP space = gMapRepP(st)->get_image_space();
    if(space != 0){
      variable real_v = get_real_var(v);
      space->remove_in_var(real_v);
      gMapRepP(st)->set_image_space(0);
    }
  }
    

  
  void gfact_db::create_gfact(const variable& v, gStoreRepP st,
                              gKeySpaceP domain_space,
                              gKeySpaceP image_space) {
    variable tmp_v =add_namespace(v) ;
    create_pure_gfact(tmp_v,st) ;
    extensional_facts += tmp_v ;
    if(domain_space != 0) set_variable_domain_space(v, st, domain_space);
    if(image_space != 0)set_variable_image_space(v, gMapRepP(st), image_space);
  }
  void gfact_db::create_gfact(const std::string& vname, gStoreRepP st,
                              gKeySpaceP domain_space ,
                              gKeySpaceP image_space ) {
    variable v = variable(vname) ;
    create_gfact(v,st, domain_space, image_space) ;
  }
  
  void gfact_db::create_gfact(const variable& v, gstore_instance &si,
                              gKeySpaceP domain_space,
                              gKeySpaceP image_space) {
    variable tmp_v = add_namespace(v) ;
    gStoreRepP st = si.Rep(); 
    create_pure_gfact(tmp_v,st) ;
    if(domain_space != 0)set_variable_domain_space(v, st, domain_space);
    if(image_space != 0)set_variable_image_space(v, gMapRepP(st), image_space);
    si.setRep(get_gvariable(v)) ;
    extensional_facts += tmp_v ; 
  }

  void gfact_db::create_gfact(const std::string& vname, gstore_instance &si,
                              gKeySpaceP domain_space ,
                              gKeySpaceP image_space ) {
    variable v = variable(vname) ;
    create_gfact(v,si, domain_space, image_space) ;
  }


  //update gfmap,
  //if new keyspaces are provided, remove v from old spaces, and connect v with  new ones
  //otherwise, assume st already has the domain space and image space set up
  void gfact_db::update_gfact(variable v, gStoreRepP st,
                              gKeySpaceP domain_space ,
                              gKeySpaceP image_space ) {
    //get real variable
    variable tmp_v = add_namespace(v) ;
      
    //tmp_v should not have synonyms
    warn(synonyms.find(tmp_v) != synonyms.end()) ;

    //update gfmap and keyspaces
    std::map<variable, gStoreRepP>::iterator mi = gfmap.find(tmp_v) ;
     
    //tmp_v should be in gfmap
    if(mi != gfmap.end()) {
      if(domain_space != 0){
        remove_variable_domain_space(v);
      }
      if(image_space != 0){
        remove_variable_image_space(v);
      }
         
      mi->second = st ;

      if(domain_space != 0)set_variable_domain_space(v, st, domain_space);
      if(image_space != 0)set_variable_image_space(v, gMapRepP(st), image_space);
    } else{
      cerr << "warning: update_gfact: fact does not exist for variable " << tmp_v
           << endl ;
    }
  }

  
  /*! remove from synonym and gfmap */
  //also remove from keyspaces
  void gfact_db::remove_gvariable(variable v) {
    std::map<variable, variable>::iterator si ;
    std::map<variable, gStoreRepP>::iterator mi ;
    if((si=synonyms.find(v)) != synonyms.end()) {
      variable real_var = remove_synonym(v) ;
      synonyms.erase(si) ;
      remove_gvariable(real_var) ;
    } else if((mi=gfmap.find(v)) != gfmap.end()) {
      // First remove any synonyms to this variable.
      variableSet syn_vars ;
      vector<map<variable,variable>::iterator > lrm ;
      for(si=synonyms.begin();si!=synonyms.end();++si)
        if(si->second == v)
          syn_vars += si->first ;
      for(variableSet::const_iterator vi=syn_vars.begin();
          vi!=syn_vars.end();++vi) {
        if((si=synonyms.find(*vi)) != synonyms.end()) {
          synonyms.erase(si) ;
        }
      }
      //before erase the variable, remove it from keyspaces
      remove_variable_domain_space(v);
      remove_variable_image_space(v);
      // Now erase the variable
      gfmap.erase(mi) ;
     
    }
  }

  //should get from gfmap or fmap? assume it is fmap now
  variableSet gfact_db::get_typed_variables() const {
    //  std::map<variable, gStoreRepP>::const_iterator mi ;
    std::map<variable, fact_info>::const_iterator mi ;
    std::map<variable, variable>::const_iterator si ;
    variableSet all_vars ;
    // for(mi=gfmap.begin();mi!=gfmap.end();++mi)
    //   all_vars += mi->first ;
    for(mi=fmap.begin();mi!=fmap.end();++mi)
      all_vars += mi->first ;
    for(si=synonyms.begin();si!=synonyms.end();++si)
      all_vars += si->first ;
    return all_vars ;
  }

  gStoreRepP gfact_db::get_gvariable(variable v) {
    variable tmp_v = add_namespace(v) ;
    tmp_v = remove_synonym(tmp_v) ;
    std::map<variable, gStoreRepP>::iterator mi =
      gfmap.find(remove_synonym(tmp_v)) ;
    if(mi == gfmap.end()) {
      // if(Loci::MPI_rank == 0)
      // 	cerr << " returning null  gStoreRep for variable " << tmp_v<< endl ;
      return gStoreRepP(0) ;
    }
    else
      return mi->second ;
  }


 
  
  // this read in routine reads the .vars file and set
  // the facts according to the default and optional rules
  std::istream& gfact_db::read_vars(std::istream& s, const rule_db& rdb) {
    gKeySpaceP universe_space = gKeySpace::get_space("UniverseSpace", "");
    bool syntax_error = false ;
    try {
      // first of all, we need to process the default and optional rules
      ruleSet special_rules = rdb.get_default_rules() ;
      // first we process the default rules
      variableSet working_vars ;
      for(ruleSet::const_iterator ri=special_rules.begin();
          ri!=special_rules.end();++ri) {
        // first we need to create the facts in the fact_db
        variableSet targets = ri->targets() ;
        rule_implP rp = ri->get_rule_implP() ;
        bool UseRule = true ;
        for(variableSet::const_iterator vi=targets.begin();
            vi!=targets.end();++vi) 
          if(vi->get_info().namespac != nspace_vec) {
            UseRule = false ;
          }
        if(!UseRule)
          continue ;
        for(variableSet::const_iterator vi=targets.begin();
            vi!=targets.end();++vi) {
          working_vars += *vi;
          // we need to get the storeRep for this variable
          storeRepP srp = rp->get_store(*vi) ;
          if(srp == 0) {
            ostringstream oss ;
            oss << "rule " << *ri << " unable to provide type for " << *vi
                << endl ;
            throw StringError(oss.str()) ;
          }
          create_fact(*vi, srp);
          
        }
        // then we need to call the compute method to set
        // the default value for this variable
        rp->initialize(*this) ; //connect var_table to frozen_variable
        rp->compute(sequence(EMPTY)) ; //assign the default value to frozen_variable
      }
      //thaw frozen variable so that the default value is in gfmap
      for(variableSet::const_iterator vi = working_vars.begin();
          vi != working_vars.end();++vi){
        storeRepP srp = get_variable(*vi);
        if(srp==0) std::cerr<< " default variable " << *vi << "is not in gfact_db" << endl; 
        gStoreRepP gsrp = srp->copy2gstore();
        if(gsrp== 0){
          std::cerr<< " default variable " << *vi << " failed in copy2gstore " << endl; 
        }
        create_gfact(*vi,gsrp, universe_space) ;//assume it is a parameter on universe_space
      }
      
      // then we process the optional rules
      special_rules = rdb.get_optional_rules() ;
      for(ruleSet::const_iterator ri=special_rules.begin();
          ri!=special_rules.end();++ri) {
        // first we need to create the facts in the fact_db
        variableSet targets = ri->targets() ;
        bool UseRule = true ;
        for(variableSet::const_iterator vi=targets.begin();
            vi!=targets.end();++vi) 
          if(vi->get_info().namespac != nspace_vec) {
            UseRule = false ;
          }
        if(!UseRule)
          continue ;


        rule_implP rp = ri->get_rule_implP() ;
        for(variableSet::const_iterator vi=targets.begin();
            vi!=targets.end();++vi) {
          // we need to get the storeRep for this variable
          storeRepP srp = rp->get_store(*vi) ;
          if(srp == 0) {
            ostringstream oss ;
            oss << "rule " << *ri << " unable to provide type for " << *vi
                << endl ;
            throw StringError(oss.str()) ;
          }
          // here we only need to set up the variable type in
          // the fact_db
          set_variable_type(*vi,srp) ;
        }
      }
    
      string vname ;
      parse::kill_white_space(s) ;
      if(s.peek()!='{') {
        throw StringError("format error in fact_db::read") ;
        return s ;
      }
      s.get() ;
      
      variableSet read_vars ;
      for(;;) {
        parse::kill_white_space(s) ;
        if(s.peek() == '}') {
          s.get() ;
          break ;
        }
        if(s.peek() == std::char_traits<char>::eof()) {
          throw StringError("unexpected EOF in fact_db::read") ;
        }
        parse::kill_white_space(s) ;
        if(parse::is_name(s)) 
          vname = parse::get_name(s) ;
        else {
          throw StringError("syntax error in fact_db::read") ;
        }
        try {
          parse::kill_white_space(s) ;
          if(!parse::get_token(s,":")) {
            throw StringError("syntax error in fact_db::read, no ':' separator") ;
          }
          
          variable var(vname) ;
          if(read_vars.inSet(var)) {
            cerr << "WARNING: Redefining variable '" << var << "' while reading in fact_db!!!!!" << endl ;
          }
          read_vars += var ;
          gStoreRepP gvp = get_gvariable(var) ;
          if(gvp == 0) {//this variable is not in gfmap
            storeRepP vp = get_variable_type(var) ;
            if(vp != 0) { //it is typed
              //transform it to gStoreRep
              gStoreRepP gsrp = vp->copy2gstore();//also assume it is a parameter in universe_space
              if(gsrp==0){
                std::cerr<< " Typed variable " << var << " failed in copy2gstore " << endl; 
              }
              create_gfact(var,gsrp) ;
            }
            gvp = get_gvariable(var) ;
          }
          if(gvp == 0) {//the variable is neither in gfmap nor tmap
            ostringstream oss ;
            oss << "variable named '" << vname
                << "' not found in database in fact_db::read." << endl ;
            throw StringError(oss.str()) ;
          }
          gvp->Input(s) ; //read in the value of variable
          parse::kill_white_space(s) ;
          if(s.peek() != '}' && !parse::is_name(s)) {
            ostringstream oss ;
            oss << "syntax error in fact_db::read while reading variable '" << vname << "'" ;
            throw StringError(oss.str()) ;
          }
        } catch(const BasicException &err) {
          err.Print(cerr) ;
          cerr << "input failed for variable " << vname << endl ;
          cerr << "near the following text: \"" ;
          while(s.peek()!=EOF&&s.peek()!='\n'&&s.peek()!='\r')
            cerr << char(s.get()) ;
          cerr << "\"" << endl ;
         
 
          syntax_error = true ;
        }
      }
    } catch(const BasicException &err) {
      err.Print(cerr) ;
      throw StringError("read fact_db failed") ;
    }
    if(syntax_error) {
      cerr << "syntax error reading fact db"  << endl ;
      throw StringError("read fact_db failed") ;
    }
    return s ;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  //the following part is temporary
  ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  void gfact_db::set_maximum_allocated(int i) {
    maximum_allocated = i ;
  }
  void gfact_db::set_minimum_allocated(int i) {
    minimum_allocated = i ;
  } 
  void gfact_db::create_pure_fact(const variable& v, storeRepP st) {
    
    if(st->RepType() == Loci::MAP || st->RepType() == Loci::STORE) {
      int max_val = st->domain().Max() ;
      maximum_allocated = max(maximum_allocated,max_val+1) ;
    }
    variable tmp_v ;
    tmp_v = v ;
    
    if(synonyms.find(tmp_v) != synonyms.end()) {
      tmp_v = remove_synonym(tmp_v) ;
      std::map<variable, fact_info>::iterator mi = fmap.find(tmp_v) ;
      if(mi==fmap.end()) {
        fmap[tmp_v].data_rep = new store_ref ;
        fmap[tmp_v].data_rep->setRep(st->getRep()) ;
      } else {
        if(typeid(st->getRep()) != typeid(mi->second.data_rep->getRep())) {
          cerr << "set_variable_type() method of fact_db changing type for variable " << tmp_v << endl ;
        }
        mi->second.data_rep->setRep(st->getRep()) ;
      }
      return ;
    }
    
    std::map<variable, fact_info>::iterator mi = fmap.find(tmp_v) ;
    if(mi != fmap.end()) {
      cerr << "WARNING: fact_db::set_variable_type retyping variable "
	   << tmp_v << endl ;
      mi->second.data_rep->setRep(st->getRep()) ;
    } else {
      fmap[tmp_v].data_rep = new store_ref ;
      fmap[tmp_v].data_rep->setRep(st->getRep()) ;
    }
    // cout << " tmp_v = " << tmp_v << endl ;
  } 
  
  std::pair<entitySet, entitySet> gfact_db::get_dist_alloc(int size) {

    if(MPI_processes > 1) {
      if(!dist_from_start) {
	dist_from_start = 1 ;
	distributed_info = new distribute_info;
      }

      int* send_buf = new int[MPI_processes] ;
      int* size_send = new int[MPI_processes] ;
      int* size_recv = new int[MPI_processes] ;
      int* recv_buf = new int[MPI_processes] ;
      for(int i = 0; i < MPI_processes; ++i) {
	send_buf[i] = maximum_allocated ;
	size_send[i] = size ;
      } 
      MPI_Alltoall(send_buf, 1, MPI_INT, recv_buf, 1, MPI_INT, MPI_COMM_WORLD) ;
      MPI_Alltoall(size_send, 1, MPI_INT, size_recv, 1, MPI_INT, MPI_COMM_WORLD) ;
      std::sort(recv_buf, recv_buf+MPI_processes) ;
      maximum_allocated = recv_buf[MPI_processes-1] ;
      int local_max = maximum_allocated ;
      int global_max = 0 ;
      for(int i = 0; i < MPI_rank; ++i)
	local_max += size_recv[i] ;
      for(int i = 0; i < MPI_processes; ++i) 
	global_max += size_recv[i] ;
      
      for(int i = 0 ; i < MPI_processes; ++i) {
	int local = maximum_allocated ;
	for(int j = 0; j < i; ++j)
	  local += size_recv[j] ;
	if(size_recv[i] > 0 )
	  init_ptn[i] += interval(local, local+size_recv[i]-1) ;
      }
      entitySet local_ivl, global_ivl ;
      if(size > 0 )
	local_ivl = entitySet(interval(local_max, local_max + size - 1)) ;
      else {
	local_ivl = EMPTY ;
      }
      
      global_ivl = entitySet(interval(maximum_allocated, maximum_allocated+global_max-1)) ;

      global_comp_entities += local_ivl;
      
      delete [] send_buf ;
      delete [] recv_buf ;
      delete [] size_send ;
      delete [] size_recv ;
      maximum_allocated = max(maximum_allocated,global_ivl.Max()+1) ;
      return(make_pair(local_ivl, global_ivl)) ;
    }
    entitySet alloc = entitySet(interval(maximum_allocated,maximum_allocated+size-1)) ;
    maximum_allocated += size ;
    
    init_ptn[0] += alloc ;
    global_comp_entities += alloc;
    return (make_pair(alloc, alloc)) ;
  }
    
  void gfact_db::update_remap(const std::vector<std::pair<int, int> > &remap_update) {
    if(Loci::MPI_processes > 1) {
      warn(!dist_from_start);
      fatal(distributed_info == NULL);
      
      for(std::vector<std::pair<int, int> >::const_iterator vi = remap_update.begin(); vi != remap_update.end(); vi++) {
        //	distributed_info->remap[vi->first] = vi->second;
        distributed_info->g2f[vi->second] = vi->first ;
      }
    }
  }

  std::pair<entitySet, entitySet> gfact_db::get_distributed_alloc(int size) {
    pair<entitySet, entitySet> allocation = get_dist_alloc(size);
    vector<pair<int, int> > remap_update;
   
    FORALL(allocation.first, ai) {
      remap_update.push_back(make_pair(ai, ai));
    }ENDFORALL;

    update_remap(remap_update);    
    return allocation;
  }

  std::pair<entitySet, entitySet> gfact_db::get_distributed_alloc(const std::vector<int> &remap_entities) {
    pair<entitySet, entitySet> allocation = get_dist_alloc(remap_entities.size());
    vector<pair<int, int> > remap_update;
    
    int i = 0;   
    FORALL(allocation.first, ai) {
      remap_update.push_back(make_pair(remap_entities[i], ai));
      i++;
    }ENDFORALL;

    update_remap(remap_update);
    return allocation;
  }

  std::pair<entitySet, entitySet> gfact_db::get_distributed_alloc(int size, int offset) {  
    pair<entitySet, entitySet> allocation = get_dist_alloc(size);
    vector<pair<int, int> > remap_update;
    
    FORALL(allocation.first, ai) {
      remap_update.push_back(make_pair(ai+offset, ai));
    }ENDFORALL;
    
    update_remap(remap_update);
    return allocation;
  }

  /*! remove from synonym and fmap */
  void gfact_db::remove_variable(variable v) {
    std::map<variable, variable>::iterator si ;
    std::map<variable, fact_info>::iterator mi ;
    if((si=synonyms.find(v)) != synonyms.end()) {
      variable real_var = remove_synonym(v) ;
      synonyms.erase(si) ;
      remove_variable(real_var) ;
    } else if((mi=fmap.find(v)) != fmap.end()) {
      // First remove any synonyms to this variable.
      variableSet syn_vars ;
      vector<map<variable,variable>::iterator > lrm ;
      for(si=synonyms.begin();si!=synonyms.end();++si)
        if(si->second == v)
          syn_vars += si->first ;
      for(variableSet::const_iterator vi=syn_vars.begin();
          vi!=syn_vars.end();++vi) {
        if((si=synonyms.find(*vi)) != synonyms.end()) {
          synonyms.erase(si) ;
        }
      }

      // Now erase the variable
      fmap.erase(mi) ;
    }
  }

  storeRepP gfact_db::get_variable(variable v) {
    variable tmp_v = add_namespace(v) ;
    tmp_v = remove_synonym(tmp_v) ;
    std::map<variable, fact_info>::iterator mi =
      fmap.find(remove_synonym(tmp_v)) ;
    if(mi == fmap.end()) {
      //      if(Loci::MPI_rank == 0)
      //	cout << " returning null  storeRep for variable " << tmp_v<< endl ;
      return storeRepP(0) ;
    }
    else
      return storeRepP(mi->second.data_rep) ;
  }
  
  gfact_db::distribute_infoP gfact_db::get_distribute_info() {
    return(distributed_info);
  }
  
  void gfact_db::put_distribute_info(distribute_infoP dp) {
    distributed_info = dp ;
  }
 
  bool gfact_db::isDistributed() {
    if(distributed_info == 0)
      return 0 ;
    else 
      return 1 ;
  }
  
  void gfact_db::rotate_vars(const std::list<variable> &lvars) {
    list<variable>::const_iterator jj ;
    jj = lvars.begin() ;
    storeRepP cp = fmap[remove_synonym(*jj)].data_rep->getRep() ;
    ++jj ;
    if(jj != lvars.end()) {
      for(;jj!=lvars.end();++jj) {
        fact_info &fd = fmap[remove_synonym(*jj)] ;
        storeRepP tmp = fd.data_rep->getRep() ;
        fd.data_rep->setRep(cp) ;
        cp = tmp ;
      }
    }
    fmap[remove_synonym(lvars.front())].data_rep->setRep(cp) ;
  }

  void gfact_db::adjust_rotation_vars(const std::list<variable>& lvars) {
    list<variable>::const_iterator jj ;
    jj = lvars.begin() ;
    storeRepP cur = fmap[remove_synonym(*jj)].data_rep->getRep() ;
    entitySet cur_dom = cur->domain() ;

    ++jj ;
    if(jj != lvars.end()) {
      for(;jj!=lvars.end();++jj) {
        fact_info& fd = fmap[remove_synonym(*jj)] ;
        storeRepP his = fd.data_rep->getRep() ;
        // first we will erase the domains outside of cur_dom
        entitySet his_dom = his->domain() ;
        entitySet out = his_dom - cur_dom ;
        if(out != EMPTY) {
          his->erase(out) ;
        }
        // then copy those missing
        entitySet missing = cur_dom - his_dom ;
        if(missing != EMPTY) {
          his->copy(cur, missing) ;
        }
      }
    }
  }

  
  void gfact_db::update_fact(variable v, storeRepP st) {
    //if st is STORE or MAP, update maximum_allocated
    if(st->RepType() == Loci::MAP || st->RepType() == Loci::STORE) {
      int max_val = st->domain().Max() ;
      maximum_allocated = max(maximum_allocated,max_val+1) ;
    }
    //add namespace
    variable tmp_v = add_namespace(v) ;

    //tmp_v should not have synonyms
    warn(synonyms.find(tmp_v) != synonyms.end()) ;
    std::map<variable, fact_info>::iterator mi = fmap.find(tmp_v) ;
    //tmp_v should be in fmap
    if(mi != fmap.end()) {
      mi->second.data_rep->setRep(st->getRep()) ;
    } else
      cerr << "warning: update_fact: fact does not exist for variable " << tmp_v
	   << endl ;
  } 
  void reorder_facts(gfact_db &facts, dMap &remap) {
    variableSet vars = facts.get_typed_variables() ;
    for(variableSet::const_iterator vi=vars.begin();vi!=vars.end();++vi) {
      storeRepP  p = facts.get_variable(*vi) ;
      if(facts.is_distributed_start())
        facts.replace_fact(*vi,(p->remap(remap))->freeze()) ;
      else
        facts.update_fact(*vi,(p->remap(remap))->freeze()) ;
    }
  }

  
  
  void serial_freeze(gfact_db &facts) {
    variableSet vars = facts.get_typed_variables() ;

    for(variableSet::const_iterator vi=vars.begin();vi!=vars.end();++vi) {
      storeRepP p = facts.get_variable(*vi) ;
      facts.replace_fact(*vi, p->freeze()) ;
    }
  }

  // experimental code to create keyspace from the
  // global registered keyspace list
  // returns "true" to indicate the methods succeeded,
  // "false" to indicate an error.
  bool
  gfact_db::create_keyspace(KeySpaceList& global_list) {
    for(KeySpaceList::Iterator ki=global_list.begin();
        ki!=global_list.end();++ki) {
      KeySpaceP kp = ki.get_p()->rr->get_space() ;
      // first get the space name
      if(!kp->named_space()) {
        if(Loci::MPI_rank == 0)
          cerr << "fact_db Error: Initializing Unnamed Keyspace!"
               << " typeid = " << typeid(*kp).name() << endl ;
        return false ;
      }
      string name = kp->get_name() ;
      map<string,KeySpaceP>::const_iterator mi = keyspace.find(name) ;
      if(mi!=keyspace.end()) {
        if(Loci::MPI_rank == 0)
          cerr << "fact_db Error: Duplicated Keyspace: "
               << name << endl ;
        return false ;
      }
      kp->set_synonyms(&synonyms) ;
      keyspace[name] = kp ;
    }
    return true ;
  }

  KeySpaceP
  gfact_db::get_keyspace(const string& kname) {
    map<string,KeySpaceP>::const_iterator mi = keyspace.find(kname) ;
    if(mi == keyspace.end())
      return KeySpaceP(0) ;
    else
      return mi->second ;
  }
  
  void
  gfact_db::init_key_manager() {
    int max_alloc = get_max_alloc() ;
    int global_max = 0 ;
    MPI_Allreduce(&max_alloc, &global_max, 1,
                  MPI_INT, MPI_MAX, MPI_COMM_WORLD) ;
    key_manager = new KeyManager(global_max+1) ;
  }
  
  void gfact_db::setupDefaults(const rule_db &rdb) {
    ruleSet special_rules = rdb.get_default_rules() ;
    // first we process the default rules
    variableSet working_vars ;
    for(ruleSet::const_iterator ri=special_rules.begin();
        ri!=special_rules.end();++ri) {
      // first we need to create the facts in the fact_db
      variableSet targets = ri->targets() ;
      rule_implP rp = ri->get_rule_implP() ;
      if(targets.size() > 1) // ignore all default rules with more than 1 arg
        continue ;
      variable v = *targets.begin() ;
      // we need to get the storeRep for this variable
      storeRepP srp = rp->get_store(v) ;
      if(srp == 0) {
        ostringstream oss ;
        oss << "default rule " << *ri << " unable to provide type for " << v
            << endl ;
        throw StringError(oss.str()) ;
      }
      std::map<variable, fact_info>::iterator mi = fmap.find(v) ;
      if(mi == fmap.end()) { // does not exist, so install
        create_fact(v,srp) ;        
        rp->initialize(*this) ;
        rp->compute(sequence(EMPTY)) ;
      }
    }
  }

  void gfact_db::write_all_hdf5(const char *filename) {
    variableSet vars = get_typed_variables() ;
    write_hdf5(filename, vars) ;
  }
  void gfact_db::read_all_hdf5(const char *filename) {
    variableSet vars = get_typed_variables() ;
    read_hdf5(filename, vars) ; 
  }
  void gfact_db::write_hdf5(const char *filename, variableSet &vars) {
    hid_t  file_id=0 ;
    if(Loci::MPI_rank == 0) 
      file_id =  H5Fcreate(filename, H5F_ACC_TRUNC,
			   H5P_DEFAULT, H5P_DEFAULT) ;
    
    
    for(variableSet::const_iterator vi = vars.begin(); vi != vars.end(); ++vi) {
      storeRepP  p = get_variable(*vi) ;
      if(p->RepType() == STORE) {
        writeContainer(file_id,variable(*vi).get_info().name,p,*this) ;
      }
    }
    if(Loci::MPI_rank == 0) 
      H5Fclose(file_id) ;
  }
  
  
  void gfact_db::read_hdf5(const char *filename, variableSet &vars) {
    hid_t  file_id=0 ;
    if(Loci::MPI_rank == 0) 
      file_id =  H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT) ;
    for(variableSet::const_iterator vi = vars.begin(); vi != vars.end(); ++vi) {
      storeRepP  p = get_variable(*vi) ;
      if(p->RepType() == STORE) {
        readContainer(file_id,variable(*vi).get_info().name,p,EMPTY,*this) ;
      }
    }
    if(Loci::MPI_rank == 0) 
      H5Fclose(file_id) ;
    
  }
  void gfact_db::Print_diagnostics() {
    std::map<variable, fact_info>::iterator mi ;
    ostringstream oss ;
    oss << "memory." ;
    oss << MPI_rank ;
    string file_name = oss.str() ;
    std::ofstream ofile(file_name.c_str(), ios::out) ;
    double total_size = 0 ;
    double total_wasted = 0 ;
    entitySet dom, total, unused ;
    for(mi = fmap.begin(); mi != fmap.end(); ++mi) {
      fact_info &finfo = mi->second ;
      if(finfo.data_rep->RepType() == STORE) {
	dom = finfo.data_rep->domain() ;
	total = interval(dom.Min(), dom.Max()) ;
	unused = total - dom ;
	total_size += finfo.data_rep->pack_size(dom) ; 
	total_wasted += finfo.data_rep->pack_size(unused) ;
      }
    }
    for(mi = fmap.begin(); mi != fmap.end(); ++mi) {
      fact_info &finfo = mi->second ;
      if(finfo.data_rep->RepType() == STORE) {
	dom = finfo.data_rep->domain() ;
	double size = finfo.data_rep->pack_size(dom) ;
	total = interval(dom.Min(), dom.Max()) ;
	unused = total - dom ;
	double wasted_space = finfo.data_rep->pack_size(unused) ;
	ofile << " ****************************************************" << endl ;
	ofile << " Total_size = " << total_size << endl ;
	ofile << "Variable = "  << mi->first << endl ;
	ofile << "Domain = " << dom << endl ;
	ofile << "Size allocated = " << size << endl ; 
	if( isDistributed() )  {
	  Loci::gfact_db::distribute_infoP d ;
	  d   = Loci::exec_current_fact_db->get_distribute_info() ;
	  entitySet my_entities = d->my_entities ; 
	  entitySet clone = dom - my_entities ;
	  double clone_size = finfo.data_rep->pack_size(clone) ;
	  ofile << "----------------------------------------------------" << endl;
	  ofile << " My_entities = " << my_entities << endl ;
	  ofile << " Clone entities = " << clone << endl ;
	  ofile << "Memory required for the  clone region  = " << clone_size << endl ;
	  ofile << "Percentage of clone memory required (of size allocated)  = " << double(double(100*clone_size) / size)<< endl ;
	  ofile << "Percentage of clone memory required (of total size allocated)  = " << double(double(100*clone_size) / total_size)<< endl ;
	  ofile << "----------------------------------------------------" << endl;
	}
	ofile << "Percentage of total memory allocated  = " << double(double(100*size) / (total_size+total_wasted)) << endl ;
	ofile << "----------------------------------------------------" << endl;
	ofile << "Total wasted size = " << total_wasted << endl ;
	ofile << "Unused entities = " << unused << endl ;
	ofile << "Wasted space = " << wasted_space << endl ;
	ofile << "Percentage of total memory wasted  = " << double(double(100*wasted_space) / (total_size + total_wasted)) << endl ;
	ofile << " ***************************************************" << endl << endl << endl ;
      }
    }
  }
  
}

  
