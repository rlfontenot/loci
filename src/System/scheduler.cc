#include "loci_globs.h"
#include "sched_tools.h"
#include "dist_tools.h"
#include "param_rule.h"

using std::map ;
using std::vector ;
using std::set ;
using std::list ;

using std::pair ;
using std::make_pair ;

using std::ostringstream ;
using std::string ;
using std::endl ;
using std::cout ;
using std::ios ;
using std::ofstream ;

///////////////////////////////////
#include <sstream>
#include <algorithm>
#include <sys/time.h> // for gettimeofday function
///////////////////////////////////

#define PROFILE_CODE

namespace Loci {
  double LociAppPeakMemory = 0 ;
  double LociAppAllocRequestBeanCounting = 0 ;
  double LociAppFreeRequestBeanCounting = 0 ;
  double LociAppPeakMemoryBeanCounting = 0 ;
  double LociAppLargestAlloc = 0 ;
  variable LociAppLargestAllocVar("EMPTY") ;
  double LociAppLargestFree = 0 ;
  variable LociAppLargestFreeVar("EMPTY") ;
  double LociAppPMTemp = 0 ;
  ////////////////////////////
  extern bool profile_memory_usage ;
  extern bool show_graphs ;
  extern void deco_depend_gr(digraph& gr,const variableSet& given) ;

  ////////////////////////////
  namespace {
    double difftime(timeval t1, timeval t2) {
      double dt1 = t1.tv_sec + t1.tv_usec*1e-6 ;
      double dt2 = t2.tv_sec + t2.tv_usec*1e-6 ;
      return dt2-dt1 ;
    }
  }

  double get_timer() {
#ifdef PROFILE_CODE
    clock_t tc ;
    static double to = 0;
    double tn,t ;
    tc = clock() ;
    tn = tc/1000000.0 ;
    t = tn - to ;
    to = tn ;
    return t ;
#else
    return -1.0 ;
#endif
  }  

  namespace {
    // pretty printing of a rule's signature
    // i.e. remove namespace info, if any
    inline std::string pretty_sig(const rule& r) {
      //the following line is the simplest, but it does not
      //include the qualify, such like "SN1:" in a super node
      //return r.get_info().desc.rule_identifier() ;
      ostringstream oss ;
      oss << "R" << r.ident() << endl ;
      return oss.str() ;
      std::string name = r.get_info().name() ;
      if(r.type() == rule::INTERNAL) {
        return name ;
      }else {
        std::string::iterator pos ;
        pos = std::find(name.begin(),name.end(),'#') ;
        return std::string( (pos==name.end()?name.begin():pos+1),name.end()) ;
      }
    }
  } // end of unnamed namespace
  
  ////////////////////////////////////////////////////////////////////
  // the following part are functions to visualize loci internal graphs
  // they create files to be used by "dot", "lefty" and "dotty"
  // programs from AT&T research
  void create_digraph_dot_file(const digraph &dg, const char* fname)
  {
    digraph dgt = dg.transpose() ;
    digraph::vertexSet allvertices = dg.get_all_vertices() ;
    digraph::vertexSet::const_iterator ri ;
    ofstream outf(fname) ;

    outf<<"digraph G {\n" ;
    //outf<<"size = \"8.5,11\";\n" ;
    
    for(ri=allvertices.begin();ri!=allvertices.end();++ri) {

      digraph::vertexSet outvertices = dg[*ri] ;
      digraph::vertexSet incomevertices = dgt[*ri] ;
      
      if(*ri < 0) {
        rule r(*ri) ;
        if(r.type() == rule::INTERNAL)
          outf<<"\""<<pretty_sig(r)<<"\""
              <<"[shape=doubleoctagon,style=filled,color=gold];\n" ;
        else
          outf<<"\""<<pretty_sig(r)<<"\""
              <<"[shape=box,style=filled,color=gold];\n" ;
      }
      else {
        variable v(*ri) ;
        outf<<"\""<<v<<"\""<<"[style=filled,color=green];\n" ;
      }
      if(incomevertices == EMPTY) {
        if(*ri >= 0) {
          variable v(*ri) ;
          outf<<"\""<<v<<"\""<<"[style=filled,color=red];\n" ;
        }else {
          rule r(*ri) ;
          outf<<"\""<<pretty_sig(r)<<"\""<<"[style=filled,color=red];\n" ;
        }
      }
      if(outvertices == EMPTY) {
        if(*ri >= 0) {
          variable v(*ri) ;
          outf<<"\""<<v<<"\""<<"[style=filled,color=blue];\n" ;
        }else {
          rule r(*ri) ;
          outf<<"\""<<pretty_sig(r)<<"\""<<"[style=filled,color=blue];\n" ;
        }
        continue ;
      }
        
      if(*ri < 0) {
        rule r(*ri) ;
        digraph::vertexSet::const_iterator ii ;
        for(ii=outvertices.begin();ii!=outvertices.end();++ii) {
          if(*ii < 0) {
            rule r2(*ii) ;
            outf<<"\""<<pretty_sig(r)<<"\""<<" -> "
                <<"\""<<pretty_sig(r2)<<"\""<<";\n" ;
          }else {
            variable v(*ii) ;
            outf<<"\""<<pretty_sig(r)<<"\""<<" -> "<<"\""<<v<<"\""<<";\n" ;
          }
        }
      }else {
        variable v(*ri) ;
        digraph::vertexSet::const_iterator ii ;
        for(ii=outvertices.begin();ii!=outvertices.end();++ii) {
          if(*ii < 0) {
            rule r(*ii) ;
            outf<<"\""<<v<<"\""<<" -> "<<"\""<<pretty_sig(r)<<"\""<<";\n" ;
          }else {
            variable v2(*ii) ;
            outf<<"\""<<v<<"\""<<" -> "<<"\""<<v2<<"\""<<";\n" ;
          }
        }        
      }
    }
    outf<<"}\n" ;
    outf.close() ;
  }

  bool is_super_rule(int rid)
  {
    if(rid >= 0) // a variable
      return false ;
    rule r(rid) ;
    std::string rqualifier = r.get_info().qualifier() ;

    return (rqualifier.substr(0,2) == "SN") ;
  }

  std::string get_assigned_cluster_name(int rid)
  {
    rule r(rid) ;
    std::string rqualifier = r.get_info().qualifier() ;
    // reads until ':'
    // e.g. SN12: .....
    std::string number = rqualifier.substr(2,rqualifier.size()-2) ;
    
    // get the cluster name
    std::string clustername = "cluster" ;
    clustername += number ;

    return clustername ;
  }

  // given rid is a supernode, connect it to its targets
  void writeout_super_rule(multiLevelGraph &mlg, int rid,
                           const digraph::vertexSet &targets,
                           std::ofstream &outf)
  {
    int repnode = *(mlg.find(rid)->graph_v-mlg.subgraphs).begin() ;
    // if the picked node is a rule
    if(repnode < 0) {
      rule r(repnode) ;
      digraph::vertexSet::const_iterator ii ;
      for(ii=targets.begin();ii!=targets.end();++ii) {
        if(*ii < 0) {// if the target is a rule
          if(is_super_rule(*ii)) {// if the target rule is a super node
            //then we pick a node in that super node
            int repnode2 = *(mlg.find(*ii)->graph_v-mlg.subgraphs).begin() ;
            if(repnode2 < 0) { // if the picked node is a rule
              rule r2(repnode2) ;
              outf<<"\""<<pretty_sig(r)<<"\""<<" -> "<<"\""
                  <<pretty_sig(r2)<<"\""
                  <<"[style=dotted,color=red,ltail="
                  << get_assigned_cluster_name(rid) << ",lhead="
                  << get_assigned_cluster_name(*ii) << "]" << ";\n" ;
            }else { // the picked node is a variable
              variable v2(repnode2) ;
              outf<<"\""<<pretty_sig(r)<<"\""<<" -> "<<"\""<<v2<<"\""
                  <<"[style=dotted,color=red,ltail="
                  << get_assigned_cluster_name(rid) << ",lhead="
                  << get_assigned_cluster_name(*ii) << "]" << ";\n" ;
            }
          }else { // if the target rule is a common single rule
            rule r2(*ii) ;
            outf<<"\""<<pretty_sig(r)<<"\""<<" -> "<<"\""<<pretty_sig(r2)<<"\""
                <<"[style=dotted,color=red,ltail="
                << get_assigned_cluster_name(rid) << "]" << ";\n" ;
          }
        }else { // if the target is variable
          variable v2(*ii) ;
          outf<<"\""<<pretty_sig(r)<<"\""<<" -> "<<"\""<<v2<<"\""
              <<"[style=dotted,color=red,ltail="
              << get_assigned_cluster_name(rid) << "]" << ";\n" ;
        }
      }
    }else { // the picked node is a variable
      variable v(repnode) ;
      digraph::vertexSet::const_iterator ii ;
      for(ii=targets.begin();ii!=targets.end();++ii) {
        if(*ii < 0) {// if the target is a rule
          if(is_super_rule(*ii)) {// if the target rule is a super node
            //then we pick a node in that super node
            int repnode2 = *(mlg.find(*ii)->graph_v-mlg.subgraphs).begin() ;
            if(repnode2 < 0) { // if the picked node is a rule
              rule r2(repnode2) ;
              outf<<"\""<<v<<"\""<<" -> "<<"\""<<pretty_sig(r2)<<"\""
                  <<"[style=dotted,color=red,ltail="
                  << get_assigned_cluster_name(rid) << ",lhead="
                  << get_assigned_cluster_name(*ii) << "]" << ";\n" ;
            }else { // the picked node is a variable
              variable v2(repnode2) ;
              outf<<"\""<<v<<"\""<<" -> "<<"\""<<v2<<"\""
                  <<"[style=dotted,color=red,ltail="
                  << get_assigned_cluster_name(rid) << ",lhead="
                  << get_assigned_cluster_name(*ii) << "]" << ";\n" ;
            }
          }else { // if the target rule is a common single rule
            rule r2(*ii) ;
            outf<<"\""<<v<<"\""<<" -> "<<"\""<<pretty_sig(r2)<<"\""
                <<"[style=dotted,color=red,ltail="
                << get_assigned_cluster_name(rid) << "]" << ";\n" ;
          }
        }else { // if the target is variable
          variable v2(*ii) ;
          outf<<"\""<<v<<"\""<<" -> "<<"\""<<v2<<"\""
              <<"[style=dotted,color=red,ltail="
              << get_assigned_cluster_name(rid) << "]" << ";\n" ;
        }
      }
    }
  }

  // given sid is a NON-super node, tid IS a super node
  // connect sid to tid
  void writeout_super_rule2(multiLevelGraph &mlg, int sid, int tid,
                            std::ofstream &outf)
  {
    if(sid < 0) { // sid is a rule
      rule srule(sid) ;

      // pick a node in tid
      int repnode = *(mlg.find(tid)->graph_v-mlg.subgraphs).begin() ;
      if(repnode < 0) { // if the picked node is a rule
        rule r2(repnode) ;
        outf<<"\""<<srule<<"\""<<" -> "<<"\""<<pretty_sig(r2)<<"\""
            <<"[style=dotted,color=red,lhead="
            << get_assigned_cluster_name(tid) << "]" << ";\n" ;
      }else { // the picked node is a variable
        variable v2(repnode) ;
        outf<<"\""<<srule<<"\""<<" -> "<<"\""<<v2<<"\""
            <<"[style=dotted,color=red,lhead="
            << get_assigned_cluster_name(tid) << "]" << ";\n" ;
      }
    }else { // sid is a variable
      variable sv(sid) ;
      
      // pick a node in tid
      int repnode = *(mlg.find(tid)->graph_v-mlg.subgraphs).begin() ;
      if(repnode < 0) { // if the picked node is a rule
        rule r2(repnode) ;
        outf<<"\""<<sv<<"\""<<" -> "<<"\""<<pretty_sig(r2)<<"\""
            <<"[style=dotted,color=red,lhead="
            << get_assigned_cluster_name(tid) << "]" << ";\n" ;
      }else { // the picked node is a variable
        variable v2(repnode) ;
        outf<<"\""<<sv<<"\""<<" -> "<<"\""<<v2<<"\""
            <<"[style=dotted,color=red,lhead="
            << get_assigned_cluster_name(tid) << "]" << ";\n" ;
      }
    }
  }

  // recursively layout from the toplevel to the bottom
  // this function just lists which nodes belongs to which cluster(supernode)
  // and does not do connections. 
  void layout_super_node(multiLevelGraph &mlg, int sid, std::ofstream &outf)
  {
    static int baseclustercounter = 0 ;
    // colors that fill the clusters
    static vector<string> cluster_color ;
    /*
    cluster_color.push_back("thistle") ;
    cluster_color.push_back("wheat") ;
    cluster_color.push_back("seashell") ;
    cluster_color.push_back("honeydew") ;
    */
    cluster_color.push_back("blueviolet") ;
    cluster_color.push_back("darkorange") ;
    cluster_color.push_back("chartreuse") ;
    cluster_color.push_back("aquamarine") ;
    static int color_counter = 0 ;
    // output cluster header
    std::string clustername = "subgraph cluster" ;
    std::stringstream ss ;
    ss << baseclustercounter++ ;
    clustername += ss.str() ;
      
    outf << clustername << " {" << '\n' ;
    // label the cluster
    outf << "label = " << "\"" << rule(sid) << "\"" << ";" << '\n' ;
    //    outf << "style=filled;" << '\n' ;
    outf << "color=" << cluster_color[(color_counter++)%4] << ";" << '\n' ;
    outf << "fontcolor=slateblue;" << '\n' ;
    outf << "fontsize=10;" << '\n' ;
    outf << "fontname=Helvetica;" << '\n' ;
    
    multiLevelGraph::subGraph *p = mlg.find(sid) ;
    ruleSet level_rules = extract_rules(p->graph_v-mlg.subgraphs) ;
    digraph gr = p->gr ;

    // Remove any rules in that are not in graph_v but are in gr
    ruleSet errs = extract_rules(gr.get_all_vertices()-p->graph_v) ;
    gr.remove_vertices(errs) ;

    digraph grt = gr.transpose() ;

    // do the layout of this subgraph
    digraph::vertexSet allvertices = gr.get_all_vertices() ;
    digraph::vertexSet::const_iterator ri ;

    // vertex set that holds all the previously placed vertices
    static digraph::vertexSet placed ;
    
    for(ri=allvertices.begin();ri!=allvertices.end();++ri) {
      digraph::vertexSet outvertices = gr[*ri] ;
      digraph::vertexSet incomevertices = grt[*ri] ;

      // only list the nodes that belongs to this super node
      // exclude any super nodes
      if( (!placed.inSet(*ri)) && (!is_super_rule(*ri))) {
        placed += *ri ;
        
        if(*ri < 0) {
          rule r(*ri) ;
          if(r.type() == rule::INTERNAL)
            outf<<"\""<<pretty_sig(r)<<"\""
                <<"[shape=doubleoctagon,style=filled,color=gold];\n" ;
          else
            outf<<"\""<<pretty_sig(r)<<"\""
                <<"[shape=box,style=filled,color=gold];\n" ;
        }
        else {
          variable v(*ri) ;
          outf<<"\""<<v<<"\""<<"[style=filled,color=green];\n" ;
        }
        if(incomevertices == EMPTY) {
          variable v(*ri) ;
          outf<<"\""<<v<<"\""<<"[style=filled,color=red];\n" ;
        }
        if(outvertices == EMPTY) {
          variable v(*ri) ;
          outf<<"\""<<v<<"\""<<"[style=filled,color=blue];\n" ;
          continue ;
        }
      }
    }
    
    for(ri=allvertices.begin();ri!=allvertices.end();++ri) {

      if(*ri < 0) {
        // recursively layout all the super nodes in this node
        if(is_super_rule(*ri)) {
          layout_super_node(mlg,*ri,outf) ;
        }
      }
    }
    // end of the subgraph cluster
    outf << "}" << '\n' ;
  }

  // connect all the nodes together
  void fill_mlg_edges(multiLevelGraph &mlg, int sid, std::ofstream &outf)
  {
    // get all the levels of the multi level graph
    vector<int> levels ;
    digraph::vertexSet working ;
    working += mlg.toplevel ;
    digraph::vertexSet::const_iterator vi ;
    while(working != EMPTY) {
      digraph::vertexSet next_level ;
      for(vi=working.begin();vi!=working.end();++vi) {
        levels.push_back(*vi) ;
        next_level += mlg.find(*vi)->graph_v & mlg.subgraphs ;
      }
      working = next_level ;
    }

    // define the usual edge color
    std::string edge_color = "slategray" ;
    
    for(vector<int>::const_iterator levelItr=levels.begin();
        levelItr!=levels.end();++levelItr) {
      multiLevelGraph::subGraph *p = mlg.find(*levelItr) ;
      digraph gr = p->gr ;

      // Remove any rules in that are not in graph_v but are in gr
      ruleSet errs = extract_rules(gr.get_all_vertices()-p->graph_v) ;
      gr.remove_vertices(errs) ;      
      digraph grt = gr.transpose() ;
      
      digraph::vertexSet allvertices = gr.get_all_vertices() ;
      digraph::vertexSet::const_iterator ri ;
    
      for(ri=allvertices.begin();ri!=allvertices.end();++ri) {
        digraph::vertexSet outvertices = gr[*ri] ;
        digraph::vertexSet incomevertices = grt[*ri] ;
      
        if(*ri < 0) {// a rule
          //if it is a super rule
          if(is_super_rule(*ri)) {
            writeout_super_rule(mlg,*ri,outvertices,outf) ;
          }else {
            rule r(*ri) ;
            digraph::vertexSet::const_iterator ii ;
            for(ii=outvertices.begin();ii!=outvertices.end();++ii) {
              if(*ii < 0) {// the target is a rule
                if(is_super_rule(*ii)) { // if the rule is a super rule
                  writeout_super_rule2(mlg,*ri,*ii,outf) ;
                }else {// the rule is a common single rule
                  rule r2(*ii) ;
                  outf<<"\""<<pretty_sig(r)<<"\""<<" -> "
                      <<"\""<<pretty_sig(r2)<<"\""
                      <<"[style=bold,color="<<edge_color<<"]"<<";\n" ;
                }
              }else {// the target is a variable
                variable v(*ii) ;
                outf<<"\""<<pretty_sig(r)<<"\""<<" -> "<<"\""<<v<<"\""
                    <<"[style=bold,color="<<edge_color<<"]"<<";\n" ;
              }
            }
          }
        }else {// a variable
          variable v(*ri) ;
          digraph::vertexSet::const_iterator ii ;
          for(ii=outvertices.begin();ii!=outvertices.end();++ii) {
            if(*ii < 0) {//the target is a rule
              if(is_super_rule(*ii)) { // the rule is a super rule
                writeout_super_rule2(mlg,*ri,*ii,outf) ;
              }else {// the rule is a common single rule
                rule r(*ii) ;
                outf<<"\""<<v<<"\""<<" -> "<<"\""<<pretty_sig(r)<<"\""
                    <<"[style=bold,color="<<edge_color<<"]"<<";\n" ;
              }
            }else {// the target is a variable
              variable v2(*ii) ;
              outf<<"\""<<v<<"\""<<" -> "<<"\""<<v2<<"\""
                  <<"[style=bold,color="<<edge_color<<"]"<<";\n" ;
            }
          }        
        }
      }
    }
  }
  
  void create_mlg_dot_file(decomposed_graph &deco, const char* fname)
  {
    multiLevelGraph &mlg = deco.mlg ;

    // out put dot file
    ofstream outf(fname) ;
    outf << "digraph G {\n" ;
    // allow edges between clusters
    outf << "compound=true;\n" ;

    // beginning by layout the toplevel of the multi-level-graph
    layout_super_node(mlg,mlg.toplevel,outf) ;
    fill_mlg_edges(mlg,mlg.toplevel,outf) ;
    // end of digraph
    outf << "}" << '\n' ;
    outf.close() ;
  }

  // visualize the content of each super node
  // from the toplevel down to the bottom
  void visualize_mlg(decomposed_graph &deco)
  {
    multiLevelGraph &mlg = deco.mlg ;
    // get all the levels of the multi level graph
    vector<int> levels ;
    digraph::vertexSet working ;
    working += mlg.toplevel ;
    digraph::vertexSet::const_iterator vi ;
    while(working != EMPTY) {
      digraph::vertexSet next_level ;
      for(vi=working.begin();vi!=working.end();++vi) {
        levels.push_back(*vi) ;
        next_level += mlg.find(*vi)->graph_v & mlg.subgraphs ;
      }
      working = next_level ;
    }

    int levelcounter = 1 ;
    for(vector<int>::const_iterator levelItr=levels.begin();
        levelItr!=levels.end();++levelItr) {
      multiLevelGraph::subGraph *p = mlg.find(*levelItr) ;
      digraph gr = p->gr ;

      // Remove any rules in that are not in graph_v but are in gr
      ruleSet errs = extract_rules(gr.get_all_vertices()-p->graph_v) ;
      gr.remove_vertices(errs) ;
      
      std::string sname = "levels" ;
      std::stringstream ss ;
      ss << levelcounter++ ;
      sname.append(ss.str()) ;
      sname.append(".dot") ;
      
      create_digraph_dot_file(gr,sname.c_str()) ;
      std::string cmd = "dotty " ;
      cmd += sname ;
      //cmd += "&" ;
      cout << "This is super node: " << "SN" << levelcounter-2 << '\n' ;
      system(cmd.c_str()) ;
      // rm the generated "dot" file
      cmd = "rm -fr " ;
      cmd += sname ;
      system(cmd.c_str()) ;
    }
  }
  /////////////////////////////////////////////////////////////////////

  executeP create_execution_schedule(const rule_db &rdb,
                                     fact_db &facts,
                                     std::string target_string,
                                     int nth) {
    num_threads = min(nth,max_threads) ;

    if(facts.is_distributed_start()) {
      if((MPI_processes > 1)) 
        get_clone(facts, rdb) ;
      else
        Loci::serial_freeze(facts) ; 
    } else {
      Loci::serial_freeze(facts) ;
    }
    
    //double timer = get_timer() ;
    sched_db scheds(facts) ;
    variableSet given = facts.get_typed_variables() ;
    variableSet target(expression::create(target_string)) ;
    if(Loci::MPI_rank==0)
      cout << "generating dependency graph..." << endl ;
    double start_time = MPI_Wtime() ;
    rule_db par_rdb ;
    par_rdb = parametric_rdb(rdb,target) ;
    digraph gr = dependency_graph(par_rdb,given,target).get_graph() ;
    // If graph is empty, return a null schedule 
    if(gr.get_target_vertices() == EMPTY)
      return executeP(0) ;

    ////////////////decorate the dependency graph/////////////////////
    //if(Loci::MPI_rank==0)
    //cout << "decorating dependency graph to include allocation..." << endl ;
    //deco_depend_gr(gr,given) ;
    //////////////////////////////////////////////////////////////////
    
    ////////////////////////////////////////////////////////////////////////
    std::string dottycmd = "dotty " ;
    if(Loci::MPI_rank==0) {
      if(show_graphs) {
        cout << "creating visualization file for dependency graph..." << endl ;
        create_digraph_dot_file(gr,"dependgr.dot") ;
        std::string cmd = dottycmd + "dependgr.dot" ;
        system(cmd.c_str()) ;
      }
    }
    ////////////////////////////////////////////////////////////////////////
    
    if(Loci::MPI_rank==0)
      cout << "setting up variable types..." << endl ;
    set_var_types(facts,gr, scheds) ;
    if(Loci::MPI_rank==0)
      cout << "decomposing graph..." << endl ;
    decomposed_graph decomp(gr,given,target) ;

    //////////////////////////////////////////////////////////////////
    if(Loci::MPI_rank==0) {
      if(show_graphs) {
        cout << "creating visualization file for decomposed graph..." << '\n' ;
        cout << "visualizing decomposed graph..." << '\n' ;
        create_mlg_dot_file(decomp,"decogr.dot") ;
        std::string cmd = dottycmd + "decogr.dot&" ;
        system(cmd.c_str()) ;
        visualize_mlg(decomp) ;
      }
    }
    //////////////////////////////////////////////////////////////////
    variableSet fact_vars, initial_vars ;
    fact_vars = facts.get_typed_variables() ;
    for(variableSet::const_iterator vi=fact_vars.begin();vi!=fact_vars.end();++vi) {
      if(variable(*vi).time().level_name() == "*" ) {
	storeRepP vp = facts.get_variable(*vi) ;
	if(vp->RepType() == STORE) {
	  ostringstream oss ;
	  oss << "source(" <<"EMPTY"<<')' ;
	  oss << ",target(" << *vi << ')' ;
	  string sig = oss.str() ;
	  rule r(sig) ;
	  entitySet t ;
	  if(par_rdb.rules_by_target(*vi) == EMPTY) {
	    if(facts.isDistributed()) {
	      fact_db::distribute_infoP d = facts.get_distribute_info() ;
	      for(size_t i = 0; i < d->copy.size(); ++i)
		t += d->copy[i].entities ;
	      initial_vars += *vi ;
	      scheds.set_existential_info(*vi, r, t) ;
	    }
	  }
	}
      }
    }
    if(Loci::MPI_rank==0)
      cout << " initial_vars = " << initial_vars << endl ;
    
    graph_compiler compile_graph(decomp, initial_vars) ;
    compile_graph.compile(facts,scheds,given,target) ;
    
    double end_time = MPI_Wtime() ;
    Loci::debugout << "Time taken for graph processing  = "
                   << end_time  - start_time << "  seconds " << endl ;
#ifdef PROFILE_CODE
    //timer = get_timer() ;
    //cout << "Graph Processing Time: "<<timer << " seconds" << endl ;
#endif
    if(Loci::MPI_rank==0)
      cout << "existential analysis..." << endl ;
    start_time = MPI_Wtime() ;
    compile_graph.existential_analysis(facts, scheds) ;
    end_time = MPI_Wtime() ;
    Loci::debugout << "Time taken for existential_analysis  = "
                   << end_time  - start_time << "  seconds " << endl ;
    ///////////////////////////////////
    /*
    if(Loci::MPI_rank==0) {
      ofstream exinfo(".exis") ;
      scheds.print_summary(facts,exinfo) ;
      exinfo.close() ;
    }
    */
    ///////////////////////////////////
    if(Loci::MPI_rank==0)
      cout << "creating execution schedule..." << endl;
    executeP sched =  compile_graph.execution_schedule
      (facts,scheds,initial_vars,num_threads) ;
    
    if(GLOBAL_OR(scheds.errors_found())) {
      if(MPI_rank == 0) {
        cerr << "error in generating schedule, dumping schedule files" << endl ;
        if(MPI_processes != 1)
          cerr << "see debug files for more information" << endl ;

      }
      ostringstream oss ;
      oss << ".schedule" ;

      if(MPI_processes > 1) {
        oss << "-" << MPI_rank ;
      }
    
      string sched_filename = oss.str() ;
      ofstream sched_file(sched_filename.c_str(),ios::out) ;
      sched->Print(sched_file) ;
      sched_file.close() ;
    
      
      Loci::Abort() ;
    }
    //scheds.print_summary(facts,Loci::debugout) ;
#ifdef PROFILE_CODE    
    //timer = get_timer() ;
    //cout << "Schedule Generation Time: " << timer << " seconds" << endl ;
#endif
    return sched ;
  }

  bool makeQuery(const rule_db &rdb, fact_db &facts, std::string query) {
    executeP schedule =  create_execution_schedule(rdb, facts, query) ;
    if(schedule == 0)
      return false ;

    // If a schedule was generated, execute it
    if(MPI_rank == 0)
      cout << "begin execution" << endl ;

    if(schedule_output) {
      // Save the schedule in the file .schedule for reference
      ostringstream oss ;
      oss << ".schedule" ;
      
      if(MPI_processes > 1) {
        oss << "-" << MPI_rank ;
      }
      
      string sched_filename = oss.str() ;
      ofstream sched_file(sched_filename.c_str(),ios::out) ;
      schedule->Print(sched_file) ;
      sched_file.close() ;
    }
    // execute schedule
    //timeval t1,t2 ;
    double st = MPI_Wtime() ;
    //gettimeofday(&t1,NULL) ;
    schedule->execute(facts) ;
    //gettimeofday(&t2,NULL) ;
    double et = MPI_Wtime() ;
    Loci::debugout << " Time taken for exectution of the schedule = " << et-st << " seconds " << endl ;
    //Loci::debugout << " Time taken for exectution of the schedule = "
    //             << difftime(t1,t2) << " seconds " << endl ;

    if(profile_memory_usage) {
      Loci::debugout << "Peak Memory used: " << LociAppPeakMemory
                     << " bytes ("
                     << LociAppPeakMemory/(1024*1024)
                     << "MB)" << endl ;
      
      Loci::debugout << "Total allocation requests: "
                     << LociAppAllocRequestBeanCounting
                     << " bytes ("
                     << LociAppAllocRequestBeanCounting/(1024*1024)
                     << "MB)" << endl ;
      
      Loci::debugout << "Total recycle requests: "
                     << LociAppFreeRequestBeanCounting
                     << " bytes ("
                     << LociAppFreeRequestBeanCounting/(1024*1024)
                     << "MB)" << endl ;
      
      Loci::debugout << "Peak Memory in bean counting: "
                     << LociAppPeakMemoryBeanCounting << " bytes ("
                     << LociAppPeakMemoryBeanCounting/(1024*1024)
                     << "MB)" << endl ;
      Loci::debugout << "The largest allocation was: "
                     << LociAppLargestAlloc << " bytes ("
                     << LociAppLargestAlloc/(1024*1024) << "MB)"
                     << " for variable: " << LociAppLargestAllocVar
                     << endl ;
      Loci::debugout << "The largest recycle was: "
                     << LociAppLargestFree << " bytes ("
                     << LociAppLargestFree/(1024*1024) << "MB)"
                     << " for variable: " << LociAppLargestFreeVar
                     << endl ;
      if(MPI_processes > 1) {
        // code to find out the largest memory bounds on all processes
        double LargestPeakMemory = 0 ;
        MPI_Allreduce(&LociAppPeakMemory,
                      &LargestPeakMemory,
                      1, MPI_DOUBLE,
                      MPI_MAX,MPI_COMM_WORLD) ;
        
        double LargestPeakMemoryBeanCounting = 0 ;
        MPI_Allreduce(&LociAppPeakMemoryBeanCounting,
                      &LargestPeakMemoryBeanCounting,
                      1, MPI_DOUBLE,
                      MPI_MAX,MPI_COMM_WORLD) ;

        Loci::debugout << endl ;
        Loci::debugout << "The global largest Peak Memory: "
                       << LargestPeakMemory << " bytes ("
                       << LargestPeakMemory/(1024*1024) << "MB)"
                       << endl ;
        Loci::debugout << "The global largest Peak Memory in bean counting: "
                       << LargestPeakMemoryBeanCounting << " bytes ("
                       << LargestPeakMemoryBeanCounting/(1024*1024) << "MB)"
                       << endl ;
      }
    }

    // communicate the execution time
    if(MPI_processes > 1) {
      double mytime = et-st ;
      double maxtime = 0 ;
      MPI_Allreduce(&mytime,&maxtime,1,
                    MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD) ;
      Loci::debugout << "Global max time taken for exectution"
                     << " of the schedule = "
                     << maxtime << " seconds " << endl ;
    }
    
    return true ;
  }

} // end of namespace Loci
