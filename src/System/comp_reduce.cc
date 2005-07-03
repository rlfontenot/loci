#include "comp_tools.h"
#include "dist_tools.h"

#include <vector>
using std::vector ;
#include <set>
using std::set ;

using std::ostream ;
using std::endl ;
using std::ostringstream ;

#include <Tools/hash_map.h>
#include "loci_globs.h"
using std::list ; 

//#define VERBOSE

namespace Loci {

  class joiner_oper : public execute_modules {
    variable joiner_var ;
    storeRepP joiner_store ;
    vector<entitySet> partition ;
    vector<storeRepP> var_vec ;
    vector<CPTR<joiner> > join_ops ;
    CPTR<joiner> joiner_op ;
  public:
    joiner_oper(variable jv, storeRepP &js, vector<entitySet> &ptn,
                vector<storeRepP> &vv, CPTR<joiner> &jo) :
      joiner_var(jv), joiner_store(js), partition(ptn),var_vec(vv),
      joiner_op(jo) {
      for(size_t i=0;i<var_vec.size();++i) {
        join_ops.push_back(jo->clone()) ;
        join_ops[i]->SetArgs(joiner_store,var_vec[i]) ;
      }
    }
    virtual void execute(fact_db &facts) ;
    virtual void Print(ostream &s) const ;
  } ;

  void joiner_oper::execute(fact_db &facts) {
    for(size_t i=0;i<var_vec.size();++i) 
      join_ops[i]->Join(sequence(partition[i])) ;
    //joiner_op->Join(joiner_store,var_vec[i],sequence(partition[i])) ;
  }

  void joiner_oper::Print(ostream &s) const {
    s << "reducing thread results for variable " << joiner_var << endl ;
    s << "reducing partitions = " << endl ;
    for(size_t i=0;i<var_vec.size();++i)
      s << "p["<<i<< "]="<<partition[i]<<endl ;
  }
  
  entitySet vmap_source_exist_apply(const vmap_info &vmi, fact_db &facts,
                                    variable reduce_var, sched_db &scheds) {
    variableSet::const_iterator vi ;
    entitySet sources = ~EMPTY ;
    for(vi=vmi.var.begin();vi!=vmi.var.end();++vi)
      if(*vi != reduce_var)
        sources &= scheds.variable_existence(*vi) ;
    vector<variableSet>::const_reverse_iterator mi ;
    for(mi=vmi.mapping.rbegin();mi!=vmi.mapping.rend();++mi) {
      entitySet working = ~EMPTY ;
      for(vi=mi->begin();vi!=mi->end();++vi) {
        FATAL(!scheds.is_a_Map(*vi)) ;
        working &= scheds.preimage(*vi,sources).first ;
      }
      sources = working ;
    }
    return sources ;
  }
  
  void apply_compiler::set_var_existence(fact_db &facts, sched_db &scheds) {
    existential_applyrule_analysis(apply, facts, scheds);
  }

  void apply_compiler::process_var_requests(fact_db &facts, sched_db &scheds) {
    exec_seq = process_applyrule_requests(apply, unit_tag, output_mapping, facts, scheds);
  }
  
  executeP apply_compiler::create_execution_schedule(fact_db &facts, sched_db &scheds) {
    CPTR<execute_list> el = new execute_list ;
    if(num_threads == 1 || !apply.get_info().rule_impl->thread_rule() ||
       exec_seq.size() < num_threads*30 ) {
      execution_factory ef(apply,sequence(exec_seq),facts, scheds) ;
      el->append_list(ef.create_product());
    } else if(!apply.get_info().output_is_parameter &&!output_mapping) {
      execute_par *ep = new execute_par ;
      parallel_schedule(ep,exec_seq,apply,facts, scheds) ;
      el->append_list(ep)  ;
    } else if(apply.get_info().output_is_parameter) {
      variableSet target = apply.targets() ;
      fatal(target.size() != 1) ;
      variable v = *(target.begin()) ;
      storeRepP sp = facts.get_variable(v) ;
      vector<entitySet> partition = partition_set(exec_seq,num_threads) ;
      vector<storeRepP> var_vec ;

      execute_par *ep = new execute_par ;
      for(size_t i=0;i<partition.size();++i) {
        storeRepP rp = sp->new_store(EMPTY) ;
        rp->allocate(partition[i]) ;
        var_vec.push_back(rp) ;
        execute_sequence *es = new execute_sequence ;

	execution_factory ef_unit(unit_tag,sequence(partition[i]), facts, scheds);
        es->append_list(ef_unit.create_product(v, rp));

	execution_factory ef_apply(apply,sequence(partition[i]), facts, scheds);
        es->append_list(ef_apply.create_product(v, rp));

        ep->append_list(es) ;
      }
      el->append_list(ep) ;
      el->append_list(new execute_thread_sync) ;
      rule_implP arule = (apply.get_info().rule_impl);
      fatal(arule == 0) ;
      CPTR<joiner> j_op = (arule)->get_joiner() ;
      fatal(j_op == 0) ;
      el->append_list(new joiner_oper(v,sp,partition,var_vec,j_op)) ;
    } else {
      variableSet target = apply.targets() ;
      fatal(target.size() != 1) ;
      variable v = *(target.begin()) ;
      storeRepP sp = facts.get_variable(v) ;
      vector<entitySet> partition = partition_set(exec_seq,num_threads) ;

      const rule_impl::info &rinfo = apply.get_info().desc ;
      execute_par *ep = new execute_par ;
      entitySet apply_domain,all_contexts ;
      vector<entitySet> shards, shard_domains ;
      for(size_t i=0;i<partition.size();++i) {
        fatal(rinfo.targets.size() != 1) ;
        entitySet context = partition[i] ;
        entitySet pdom = vmap_target_exist(*rinfo.targets.begin(),facts,context, scheds) ;
        entitySet rem = pdom & apply_domain ;
        if(rem != EMPTY) {
          entitySet compute = rem ;
          const vmap_info &vmi = *rinfo.targets.begin() ;
          vector<variableSet>::const_reverse_iterator mi ;
          for(mi=vmi.mapping.rbegin();mi!=vmi.mapping.rend();++mi) {
            entitySet working = EMPTY ;
            variableSet::const_iterator vi;
            for(vi=mi->begin();vi!=mi->end();++vi) {
              FATAL(!scheds.is_a_Map(*vi)) ;
              working |= scheds.preimage(*vi,compute).second ;
            }
            compute = working ;
          }
          compute &= partition[i] ;
          shards.push_back(compute) ;
          entitySet sdom = vmap_target_exist(vmi,facts,compute, scheds) ;
          shard_domains.push_back(sdom) ;
          context &= ~compute ;
        }
        apply_domain |= pdom ;
        all_contexts |= partition[i] ;
	execution_factory ef(apply,sequence(context),facts, scheds);
	ep->append_list(ef.create_product());
      }
      if(shards.size() == 0) {
        el->append_list(ep) ;
        el->append_list(new execute_thread_sync) ;
      } else {
        ep->append_list(new execute_sequence) ;
        bool disjoint = true ;
        entitySet dom_tot ;
        for(size_t i=0;i<shards.size();++i) {
          if((shard_domains[i] & dom_tot) != EMPTY)
            disjoint = false ;
          dom_tot += shard_domains[i] ;
        }
        variableSet target = apply.targets() ;
        fatal(target.size() != 1) ;
        variable v = *(target.begin()) ;
        storeRepP sp = facts.get_variable(v) ;
        fatal(sp == 0) ;
      
        vector<storeRepP> var_vec ;
      
        for(size_t i=0;i<shards.size();++i) {
          storeRepP rp = sp->new_store(EMPTY) ;
          rp->allocate(shard_domains[i]) ;

          var_vec.push_back(rp) ;
          execute_sequence *es = new execute_sequence ;
	  execution_factory ef_unit(unit_tag,sequence(shard_domains[i]),
				    facts,scheds);
          es->append_list(ef_unit.create_product(v, rp));

	  execution_factory ef_apply(apply,sequence(shards[i]),
				    facts, scheds);
          es->append_list(ef_apply.create_product(v, rp));
          ep->append_list(es) ;
        }
      
        el->append_list(ep) ;
        el->append_list(new execute_thread_sync) ;
        //-----Now join the partial results

        rule_implP arule = (apply.get_info().rule_impl);
        fatal(arule == 0) ;

        if(disjoint) {
          execute_par *epj = new execute_par ;
          epj->append_list(new execute_sequence) ;
          for(size_t i=0;i<shard_domains.size();++i) {
            vector<entitySet> ve ;
            vector<storeRepP> vv ;
            ve.push_back(shard_domains[i]) ;
            vv.push_back(var_vec[i]) ;
            CPTR<joiner> j_op = (arule)->get_joiner() ;
            fatal(j_op == 0) ;
            epj->append_list(new joiner_oper(v,sp,ve,vv,j_op)) ;
          }
          el->append_list(epj) ;
          el->append_list(new execute_thread_sync) ;
        } else {
          for(size_t i=0;i<shard_domains.size();++i) {
            execute_par *epj = new execute_par ;
            vector<entitySet> decompose = partition_set(shard_domains[i],num_threads) ;
            vector<storeRepP> vv ;
            vv.push_back(var_vec[i]) ;
            for(size_t j=0;j<decompose.size();++j) {
              vector<entitySet> ve ;
              ve.push_back(decompose[j]) ;
              CPTR<joiner> j_op = (arule)->get_joiner() ;
              fatal(j_op == 0) ;
              epj->append_list(new joiner_oper(v,sp,ve,vv,j_op)) ;
            }
            el->append_list(epj) ;
            el->append_list(new execute_thread_sync) ;
          }
        }
      }
    }
    el->append_list(new execute_thread_sync) ;
    return executeP(el) ;
  }

  //targetsize >= tcount, is real size(bytes) of target buffer while tcount is valid bytes in target buffer
  void myJoin(void * source, int scount, unsigned char *&target, int &tcount,  int &targetSize, vector<CPTR<joiner> > join_op) { 
    int targetUnpackPosition = 0; //position pointer for unpacking target to tp
    int sourceUnpackPosition = 0; //position pointer for unpacking source to sp
    int targetPackPosition = 0; //position pointer for packing target from tp
    int totalCountLeft = tcount; //total number of bytes still to be processesed in next loops
    entitySet e;
    sequence seq;  
    tcount = 0;
    int num_reduce = join_op.size();
    for(int i = 0; i < num_reduce; i++) {
      storeRepP tp = join_op[i]->getTargetRep();
      storeRepP sp = join_op[i]->getTargetRep();
      targetUnpackPosition = targetPackPosition;
      tp->unpack(target, targetUnpackPosition, targetSize, e);
      sp->unpack(source, sourceUnpackPosition, scount, e);
      join_op[i]->SetArgs(tp, sp);
      int sizeBeforeJoin = tp->pack_size(e);
      totalCountLeft -= sizeBeforeJoin;
      join_op[i]->Join(seq);
      int sizeAfterJoin = tp->pack_size(e);
      if(targetSize < targetPackPosition + sizeAfterJoin + totalCountLeft ) {
	while(targetSize < targetPackPosition + sizeAfterJoin + totalCountLeft)
	  targetSize *= 2;
	unsigned char * temp_target = new unsigned char[targetSize];
	for(int j = 0; j < targetPackPosition ; j++) {
	  temp_target[j] = target[j];
	}
	tp->pack(temp_target, targetPackPosition, targetSize, e);
	for(int j = 0; j < totalCountLeft; j++) {
	  temp_target[targetPackPosition+j] = target[targetPackPosition-sizeAfterJoin+sizeBeforeJoin+j];
	}
	delete [] target;
	target = temp_target;
      }
      else {
	if(sizeAfterJoin != sizeBeforeJoin) {
	  if(sizeBeforeJoin < sizeAfterJoin) {
	    for(int j = 0; j < totalCountLeft; j++) {
	      target[targetPackPosition+sizeAfterJoin+totalCountLeft-1-j] = target[targetPackPosition+sizeBeforeJoin+totalCountLeft-1-j];
	    }
	  }
	  else {
	    for(int j = 0; j < totalCountLeft; j++) {
	      target[targetPackPosition+sizeAfterJoin+j] = target[targetPackPosition+sizeBeforeJoin+j];
	    }
	  }
	}
	tp->pack(target, targetPackPosition, targetSize, e);
      }
    
      tcount += sizeAfterJoin;
    }
  }

  unsigned char * groupAllReduce(void *sbuf, int scount, int &rcount, 
				 vector<CPTR<joiner> > join_op, MPI_Comm comm) {
    //get number of processors and my rank
    int myid, numprocs;
    MPI_Comm_size(comm,&numprocs) ;
    MPI_Comm_rank(comm,&myid) ;
    rcount = scount;
    int validCount = scount;
    unsigned char * rbuf = new unsigned char[rcount];
    memcpy(rbuf, sbuf, rcount);
    if(numprocs == 1) 
      return rbuf;

    MPI_Status status;

    int dimension = (int)(ceil(log10((double)numprocs)/log10((double)2))); //dimension of hypercube
    int extraProcessors = (1<<dimension) - numprocs; //extra virtual processors needed
    bool meVirtual = false;
    int myVirtualId = myid ;
    unsigned char *rVirtualBuf; //buffer for virtual processor
    int virtualRCount = 0; //number of elements in rVirtualBuf
    int virtualValidCount = 0;
    //if I am also acting as virtual processor
    if((myid ^ (1<<(dimension-1))) >= numprocs) {
      meVirtual = true;
      virtualRCount = rcount;
      rVirtualBuf = new unsigned char[virtualRCount];
      myVirtualId = myid ^ (1<<(dimension-1));
    }

    int countTag = 0; //tag to send send and recvcount between processors
    int tag1 = 1; //tag for sendrecv call    
    int tag2 = 2; //tag for send and recv calls
    unsigned char *tempRecvBuf; //temperory receive buffer
    int tempCount;
    int partner; //communication partner's id
  
    //for keeping data about virtual processor's buffer has been initialized or not
    //assigned[i] keeps the data of the virtual processor whose id is i+numprocs
    bool *assigned = 0;
    if(extraProcessors > 0) 
      assigned = new bool[extraProcessors]; 
  
  //assign false value for all virtual processors since their buffer is not initalized
    for(int i = 0; i < extraProcessors; i++)
      assigned[i] = false;
    storeRepP tp, sp;
  
    for(int i = 0; i < dimension ; i++) {
      partner = (myid ^ (1<<i));
    
      if(partner >= numprocs) { //if partner is virtual
	if(assigned[partner-numprocs]) { //if virtual partner's buffer has been initialized
	  partner = partner ^ (1<<(dimension-1)); //find the real processor who own the virtual id
	  if(partner == myid) { //if I own the virtual id
	  
	    myJoin(rVirtualBuf, virtualValidCount, rbuf, validCount, rcount, join_op);
	    if(virtualRCount < validCount) {
	      while(virtualRCount < validCount)
		virtualRCount *= 2;
	      delete [] rVirtualBuf;
	      rVirtualBuf = new unsigned char[virtualRCount];
	    }
	    virtualValidCount = validCount;
	    memcpy(rVirtualBuf, rbuf, validCount);
	  }
	  else {
	    MPI_Sendrecv(&validCount, 1, MPI_INT, partner, countTag, 
			 &tempCount, 1, MPI_INT, partner, countTag, comm, &status);
	    tempRecvBuf = new unsigned char[tempCount];
	    MPI_Sendrecv(rbuf, validCount, MPI_PACKED, partner, tag1, 
			 tempRecvBuf, tempCount, MPI_PACKED, partner, tag1, comm, &status);
	    myJoin(tempRecvBuf, tempCount, rbuf, validCount, rcount, join_op);
	    delete [] tempRecvBuf;
	  }
	}
	else { //if virtual partner's buffer has not been initialized
	  partner = partner ^ (1<<(dimension-1));
	  if(partner == myid) { // if I own the virtual id
	    if(virtualRCount < validCount) {
	      while(virtualRCount < validCount)
		virtualRCount *= 2;
	      delete [] rVirtualBuf;
	      rVirtualBuf = new unsigned char[virtualRCount];
	    }
	    memcpy(rVirtualBuf, rbuf, validCount);
	    virtualValidCount = validCount;
	  }
	  else {
	    MPI_Send(&validCount, 1, MPI_INT, partner, countTag, comm);
	 
	    MPI_Send(rbuf, validCount, MPI_PACKED, partner, tag2, comm);
	  }
	}
      }
   
      else { //if partner is not virtual
	MPI_Sendrecv(&validCount, 1, MPI_INT, partner, countTag, 
		     &tempCount, 1, MPI_INT, partner, countTag, comm, &status);
      
	tempRecvBuf = new unsigned char[tempCount];
      
	MPI_Sendrecv(rbuf, validCount, MPI_PACKED, partner, tag1, 
		     tempRecvBuf, tempCount, MPI_PACKED, partner, tag1, comm, &status);
	myJoin(tempRecvBuf, tempCount, rbuf, validCount, rcount, join_op);
	delete [] tempRecvBuf;
      }
    
      if(meVirtual) { //if i act also as a virtual processor
	partner = (myVirtualId ^ (1<<i));
	if(assigned[myVirtualId-numprocs]) { //if my virtual buffer has been initialized
	  if(partner >= numprocs) { //if partner is virtual
	    if(assigned[partner-numprocs]) { //if virtual partner's buffer has been initialized
	      partner = partner ^ (1<<(dimension-1)); //find the real processor who own the virtual id
	      MPI_Sendrecv(&virtualValidCount, 1, MPI_INT, partner, countTag, 
			   &tempCount, 1, MPI_INT, partner, countTag, comm, &status);
	      tempRecvBuf = new unsigned char[tempCount];
	      MPI_Sendrecv(rVirtualBuf, virtualValidCount, MPI_PACKED, partner, tag1, 
			   tempRecvBuf, tempCount, MPI_PACKED, partner, tag1, comm, &status);

	      myJoin(tempRecvBuf, tempCount, rVirtualBuf, virtualValidCount, virtualRCount, join_op);
	      delete [] tempRecvBuf;
	    }
	    else { //if virtual partner's buffer has not been initialized yet
	      partner = partner ^ (1<<(dimension-1)); //find real processor who own the virtual id
	      MPI_Send(&virtualValidCount, 1, MPI_INT, partner, countTag, comm);
	      MPI_Send(rVirtualBuf, virtualValidCount, MPI_PACKED, partner, tag2, comm);
	    }
	  }
	  else { //if partner is not virtual
	    //if virtual partner is myid then don't do anything because it has been taken care of 
	    //in the loops before this loop. if virtual partner is not myid then
	    if(myid != partner) { 
	      MPI_Sendrecv(&virtualValidCount, 1, MPI_INT, partner, countTag, 
			   &tempCount, 1, MPI_INT, partner, countTag, comm, &status);
	      tempRecvBuf = new unsigned char[tempCount];
	      MPI_Sendrecv(rVirtualBuf, virtualValidCount, MPI_PACKED, partner, tag1,
			   tempRecvBuf, tempCount, MPI_PACKED, partner, tag1, comm, &status);
	      myJoin(tempRecvBuf, tempCount, rVirtualBuf, virtualValidCount, virtualRCount, join_op);
	      delete [] tempRecvBuf;
	    }
	  }
	}
	else { //if my virtual buffer has not been initialized
	  if(partner >= numprocs) { //if partner is virtual
	    if(assigned[partner - numprocs]) { //if virtual partner's buffer has been initialized
	      partner = partner ^ (1<<(dimension-1));
	      MPI_Recv(&tempCount, 1 , MPI_INT, partner, countTag, comm, &status);
	      if(virtualRCount < tempCount) {
		while(virtualRCount < tempCount)
		  virtualRCount *= 2;
		delete [] rVirtualBuf;
		rVirtualBuf = new unsigned char[virtualRCount];
	      }
	      MPI_Recv(rVirtualBuf, virtualRCount, MPI_PACKED, partner,tag2, comm, &status);
	      virtualValidCount = tempCount;
	    }
	  }
	  else { 
	    //if partner is not virtual
	    //if virtual partner is myid then don't do anything because it has been taken care of 
	    //in the loops before this loop. if virtual partner is not myid then
	    if(partner != myid) { 
	      MPI_Recv(&tempCount, 1, MPI_INT, partner, countTag, comm, &status);
	      if(virtualRCount < tempCount) {
		while(virtualRCount < tempCount)
		  virtualRCount *= 2;
		delete [] rVirtualBuf;
		rVirtualBuf = new unsigned char[virtualRCount];
	      }
	      MPI_Recv(rVirtualBuf, virtualRCount, MPI_PACKED, partner,tag2, comm, &status);
	      virtualValidCount = tempCount;
	    }
	  }
	}
      }
      //find which virtual processors are initialized in this loop
      for(int j = 0; j < extraProcessors; j++) {
	//if virtual processor's partner is virtual
	if(((j+numprocs) ^ (1<<i)) >= numprocs) {
	  //if virtual partner is initialized then now virtual processor is initialized
	  if(assigned[((j+numprocs)^(1<<i))-numprocs])
	    assigned[j] = true;
	}
	else //if virtual processor's partner is not virtual then virtual processor is initialized now
	  assigned[j] = true;
      }
    }
    //rcount = count;
    if(extraProcessors > 0)
      delete [] assigned;
    if(meVirtual) {
      delete [] rVirtualBuf;
    }
    return rbuf;
  }

  vector<CPTR<joiner> > global_join_ops ;
  void create_user_function(void *send_ptr, void *result_ptr, int *size, MPI_Datatype* dptr) {
    entitySet e ;
    sequence seq ;
    int unpack_send_position = 0 ;
    int unpack_result_position = 0;
    int pack_result_position = 0;
    for(size_t i = 0; i < global_join_ops.size(); i++) {
      storeRepP sp, tp ;
      sp = global_join_ops[i]->getTargetRep() ;
      tp = global_join_ops[i]->getTargetRep() ;
      sp->unpack(send_ptr, unpack_send_position, *size, seq) ;
      tp->unpack(result_ptr, unpack_result_position, *size, seq) ;
      global_join_ops[i]->SetArgs(tp, sp) ;
      global_join_ops[i]->Join(seq) ;
      tp->pack(result_ptr, pack_result_position, *size, e) ;
    }
  } 
  
  execute_param_red::execute_param_red(vector<variable> red, vector<rule> unit, vector<CPTR<joiner> > j_op) {
    reduce_vars = red ;
    unit_rules = unit ;
    join_ops = j_op ;
    MPI_Op_create(&create_user_function, 0, &create_join_op) ;
  }
  execute_param_red::~execute_param_red() {
    MPI_Op_free(&create_join_op) ;
  }
  
  void execute_param_red::execute(fact_db &facts) {
    unsigned char *send_ptr, *result_ptr;
    int size = 0;
    entitySet e ;
    sequence seq ;
    vector<storeRepP> sp;
    for(size_t i = 0; i < reduce_vars.size(); i++) {
      sp.push_back(facts.get_variable(reduce_vars[i])) ;
      size += sp[i]->pack_size(e);
    }
    send_ptr = new unsigned char[size] ;
    int position = 0;
    for(size_t i = 0; i < sp.size(); i++) {
      sp[i]->pack(send_ptr, position, size, e) ;
    }
#ifndef GROUP_ALLREDUCE
    result_ptr = new unsigned char[size] ;
   
    global_join_ops = join_ops;
    MPI_Allreduce(send_ptr, result_ptr, size, MPI_PACKED, create_join_op, MPI_COMM_WORLD) ;
    position = 0;
    for(size_t i = 0; i < sp.size(); i++) {
      sp[i]->unpack(result_ptr, position, size, seq) ;
    }
#else
    int recv_buf_size;
    result_ptr = groupAllReduce(send_ptr, size, recv_buf_size, join_ops, MPI_COMM_WORLD);
    position = 0;
    for(size_t i = 0; i < sp.size(); i++) {
      sp[i]->unpack(result_ptr, position, recv_buf_size, seq) ;
    }
#endif
    delete [] send_ptr ;
    delete [] result_ptr ;
  }
  void execute_param_red::Print(ostream &s) const {
    for(size_t i = 0 ; i < reduce_vars.size(); i++) 
      s << "param reduction on " << reduce_vars[i] << endl ;
  }
  void reduce_param_compiler::set_var_existence(fact_db &facts, sched_db &scheds)  {
    if(facts.isDistributed()) {
      	fact_db::distribute_infoP d = facts.get_distribute_info() ;
      for(size_t i = 0; i < unit_rules.size(); i++) {
	entitySet targets ;
	targets = scheds.get_existential_info(reduce_vars[i], unit_rules[i]) ;
	targets += send_entitySet(targets, facts) ;
	targets &= d->my_entities ;
	targets += fill_entitySet(targets, facts) ;
	scheds.set_existential_info(reduce_vars[i],unit_rules[i],targets) ;
      }
    }
  }
  
  void reduce_param_compiler::process_var_requests(fact_db &facts, sched_db &scheds) {
    if(facts.isDistributed()) {
      for(size_t i = 0; i < unit_rules.size(); i++) {
	entitySet requests = scheds.get_variable_requests(reduce_vars[i]) ;
	requests += send_entitySet(requests, facts) ;
	scheds.variable_request(reduce_vars[i],requests) ;
      }
    }
  } 
  
  executeP reduce_param_compiler::create_execution_schedule(fact_db &facts, sched_db &scheds) {
   
    if(facts.isDistributed()) {
      CPTR<execute_sequence> el = new execute_sequence ; 
      //el->append_list(new execute_thread_sync) ;
 
      el->append_list(new execute_param_red(reduce_vars, unit_rules, join_ops)) ;
      return executeP(el) ;
    }
    ostringstream oss ;
    for(size_t i = 0; i < reduce_vars.size(); i++) 
      oss << "reduce param " << reduce_vars[i] << std::endl;
    return executeP(new execute_msg(oss.str())) ;
  }
  
  void reduce_store_compiler::set_var_existence(fact_db &facts, sched_db &scheds)  {
  }
  
  void swap_send_recv(list<comm_info> &cl) {
    list<comm_info>::iterator li ;
    for(li=cl.begin();li!=cl.end();++li) {
      entitySet tmp  = li->send_set ;
      li->send_set = entitySet(li->recv_set) ;
      li->recv_set = sequence(tmp) ;
    }
  }

  void reduce_store_compiler::process_var_requests(fact_db &facts, sched_db &scheds) {
    if(facts.isDistributed()) {
      fact_db::distribute_infoP d = facts.get_distribute_info() ;
      variableSet vars ;
      vars += reduce_var ;

      //Find out duplication of variables that are associated with rules
      //that compute reduce variables
      if(duplicate_work)
	set_duplication_of_variables(vars, scheds);
      list<comm_info> request_comm = barrier_process_rule_requests(vars,facts, scheds) ;
      entitySet requests = scheds.get_variable_requests(reduce_var) ;
      entitySet shadow;

      //shadow should be empty for the variable which is going to be duplicated
      //otherwise find out shadow that is going to be entities need to be sent
      //to the owner processor 
      if(!duplicate_work || !scheds.is_duplicate_variable(reduce_var))
	shadow = scheds.get_variable_shadow(reduce_var) ;

      shadow &= requests ;
      list<comm_info> slist ;

      //If variable is duplicate variable, we don't need to send request to the 
      //other procesors because owner processor can ablways compute requests
      if(!duplicate_work || !scheds.is_duplicate_variable(reduce_var)) {
	entitySet response = send_requests(shadow, reduce_var,facts,slist) ;
	swap_send_recv(slist) ;
	rlist = sort_comm(slist,facts) ;
      }
      clist = sort_comm(request_comm,facts) ;
      
#ifdef VERBOSE
      if(shadow != EMPTY) {
	debugout << "shadow = " << shadow << endl ;
	shadow -= d->my_entities ;
	debugout << "shadow/my_entites = " << shadow << endl ;
      }
#endif
    }
  }

class execute_comm_reduce : public execute_modules {
  vector<pair<int,vector<send_var_info> > > send_info ;
  vector<pair<int,vector<recv_var_info> > > recv_info ;
  std::vector<std::vector<storeRepP> > send_vars ; 
  std::vector<std::vector<storeRepP> > recv_vars ; 
  CPTR<joiner> join_op ;
  int *maxr_size, *maxs_size, *r_size, *s_size, *recv_sizes ;
  unsigned char **recv_ptr , **send_ptr ;
  MPI_Request *request;
  MPI_Status *status ;
public:
  execute_comm_reduce(list<comm_info> &plist, fact_db &facts,
		      CPTR<joiner> jop) ;
  ~execute_comm_reduce() ;
  virtual void execute(fact_db &facts) ;
  virtual void Print(std::ostream &s) const ;
} ; 

execute_comm_reduce::execute_comm_reduce(list<comm_info> &plist,
					 fact_db &facts,
					 CPTR<joiner> jop) {
  join_op = jop ;
  HASH_MAP(int,vector<send_var_info> ) send_data ;
  HASH_MAP(int,vector<recv_var_info> ) recv_data ;
  list<comm_info>::const_iterator cli ;
  intervalSet send_procs, recv_procs ;
  for(cli=plist.begin();cli!=plist.end();++cli) {
    variable v = cli->v ;
    if(cli->send_set.size() > 0) {
      int send_proc = cli->processor ;
      send_procs += send_proc ;
      entitySet send_set = cli->send_set ;
      send_data[send_proc].push_back(send_var_info(v,send_set)) ;
    }
    if(cli->recv_set.size() > 0) {
      int recv_proc = cli->processor ;
      sequence recv_seq = cli->recv_set ;
      recv_procs += recv_proc ;
      recv_data[recv_proc].push_back(recv_var_info(v,recv_seq)) ;
    }
  }
    
    for(intervalSet::const_iterator ii=send_procs.begin();
        ii!=send_procs.end();
        ++ii) {
      send_info.push_back(make_pair(*ii,send_data[*ii])) ;
      send_vars.push_back(std::vector<storeRepP>()) ; 
      for(size_t i=0;i<send_data[*ii].size();++i) 
        send_vars.back().push_back(facts.get_variable(send_data[*ii][i].v)) ; 
    }
    
    for(intervalSet::const_iterator ii=recv_procs.begin();
        ii!=recv_procs.end();
        ++ii) {
      recv_info.push_back(make_pair(*ii,recv_data[*ii])) ;
      recv_vars.push_back(std::vector<storeRepP>()) ; 
      for(size_t i=0;i<recv_data[*ii].size();++i) 
        recv_vars.back().push_back(facts.get_variable(recv_data[*ii][i].v)) ; 
    }
    
    int nsend = send_info.size() ;
    int nrecv = recv_info.size() ;
    r_size = new int[nrecv] ;
    recv_sizes = new int[nrecv] ; 
    maxr_size = new int[nrecv] ; 
    maxs_size = new int[nsend] ;
    s_size = new int[nsend] ;
    for(int i = 0; i < nrecv; ++i) {
      r_size[i] = 0 ;
      maxr_size[i] = sizeof(int) ;
    }
    for(int i = 0; i < nsend; ++i) {
      maxs_size[i] = sizeof(int) ;
      s_size[i] = 0 ;
    }
    
    recv_ptr = new unsigned char*[max(nrecv,1)] ;
    send_ptr = new unsigned char*[max(nsend,1)] ;
    request =  new MPI_Request[nrecv] ;
    status =  new MPI_Status[nrecv] ;
}
  execute_comm_reduce::~execute_comm_reduce() {
    delete [] maxr_size ;
    delete [] maxs_size ;
    delete [] recv_sizes ; 
    delete [] r_size ;
    delete [] s_size ;
    delete [] recv_ptr ;
    delete [] send_ptr ;
    delete [] request ;
    delete [] status ;
  }
  static unsigned char *recv_ptr_buf = 0; 
  static int recv_ptr_buf_size = 0; 
  static unsigned char *send_ptr_buf = 0 ; 
  static int send_ptr_buf_size = 0 ; 
  
  
  void execute_comm_reduce::execute(fact_db  &facts) {
    const int nrecv = recv_info.size() ;
    int resend_size = 0, rerecv_size = 0 ;
    std::vector<int> send_index ;
    std::vector<int> recv_index ; 
    int total_size = 0 ;
    MPI_Request *re_request = 0;
    MPI_Status *re_status = 0 ;
    for(int i=0;i<nrecv;++i) {
      r_size[i] = maxr_size[i] ;
      total_size += maxr_size[i] ;
    }

    if(recv_ptr_buf_size < total_size) { 
      if(recv_ptr_buf) 
        delete[] recv_ptr_buf ; 
      recv_ptr_buf = new unsigned char[total_size] ; 
      recv_ptr_buf_size = total_size ; 
    } 
    recv_ptr[0] = recv_ptr_buf ; 
    for(int i=1;i<nrecv;++i)
      recv_ptr[i] = recv_ptr[i-1] + r_size[i-1] ;
    
    for(int i=0;i<nrecv;++i) {
      int proc = recv_info[i].first ;
      MPI_Irecv(recv_ptr[i], r_size[i], MPI_PACKED, proc, 1,
                MPI_COMM_WORLD, &request[i]) ;
    }
    total_size = 0 ;
    const int nsend = send_info.size() ;
    entitySet resend_procs, rerecv_procs ;
    for(int i=0;i<nsend;++i) {
      s_size[i] = 0 ;
      for(size_t j=0;j<send_info[i].second.size();++j) {
        //facts.get_variable(send_info[i].second[j].v) ;
	storeRepP sp = send_vars[i][j] ;

        s_size[i] += sp->pack_size(send_info[i].second[j].set) ;
      }
      if((s_size[i] > maxs_size[i]) || ( s_size[i] == sizeof(int))) {
	if(s_size[i] > maxs_size[i])
	  maxs_size[i] = s_size[i] ;
	int proc = send_info[i].first ;
	s_size[i] = sizeof(int) ;
	resend_procs += proc ;
	send_index.push_back(i) ;
      }
      total_size += maxs_size[i] ;
    }
    if(send_ptr_buf_size < total_size) { 
      if(send_ptr_buf) 
        delete[] send_ptr_buf ; 
      send_ptr_buf = new unsigned char[total_size] ; 
      send_ptr_buf_size = total_size ; 
    } 
    send_ptr[0] = send_ptr_buf ; 
    for(int i = 1; i < nsend; i++)
      send_ptr[i] = send_ptr[i-1] + s_size[i-1] ;
    // Pack the buffer for sending 
    for(int i=0;i<nsend;++i) {
      int loc_pack = 0 ;
      if(!resend_procs.inSet(send_info[i].first)) {
	for(size_t j=0;j<send_info[i].second.size();++j) {
	  storeRepP sp = send_vars[i][j] ;//facts.get_variable(send_info[i].second[j].v) ;
	  sp->pack(send_ptr[i], loc_pack,s_size[i],send_info[i].second[j].set);
	}
      }
      else
	MPI_Pack(&maxs_size[i], sizeof(int), MPI_BYTE, send_ptr[i], s_size[i], &loc_pack, MPI_COMM_WORLD) ; 
      
    }
    // Send Buffer
    for(int i=0;i<nsend;++i) {
      int proc = send_info[i].first ;
      MPI_Send(send_ptr[i],s_size[i],MPI_PACKED,proc,1,MPI_COMM_WORLD) ;
    }
    if(nrecv > 0) {
#ifdef DEBUG
      int err =
#endif
        MPI_Waitall(nrecv, request, status) ;
      FATAL(err != MPI_SUCCESS) ;
      for(int i = 0 ; i < nrecv; i++) {
        int rcv_sizes ;
	MPI_Get_count(&status[i], MPI_BYTE, &rcv_sizes) ;  
	if(rcv_sizes == sizeof(int)) {
	  rerecv_procs += recv_info[i].first ;
	  recv_index.push_back(i) ;
	}
      }
    }
    for(int i=0;i<nrecv;++i) {
      int loc_unpack = 0;
      if(rerecv_procs.inSet(recv_info[i].first)) {
	int temp ;
	MPI_Unpack(recv_ptr[i], r_size[i], &loc_unpack, &temp, sizeof(int), MPI_BYTE, MPI_COMM_WORLD) ;
	if(temp > maxr_size[i])
	  maxr_size[i] = temp ;
      }
      else
	for(size_t j=0;j<recv_info[i].second.size();++j) {
	  storeRepP sp = recv_vars[i][j] ;//facts.get_variable(recv_info[i].second[j].v) ;
	  storeRepP sr = sp->new_store(EMPTY) ;
	  sr->allocate(entitySet(recv_info[i].second[j].seq)) ;
          
	  sr->unpack(recv_ptr[i], loc_unpack, r_size[i],
		     recv_info[i].second[j].seq) ;
	  CPTR<joiner> op = join_op->clone() ;
	  op->SetArgs(sp,sr) ;
	  op->Join(recv_info[i].second[j].seq) ;
	}
    }
    rerecv_size = rerecv_procs.size() ;
    resend_size = resend_procs.size() ;
    
    if(rerecv_size > 0) {
      re_request =  new MPI_Request[rerecv_size] ;
      re_status =  new MPI_Status[rerecv_size] ;
    }
    for(int i = 0; i < rerecv_size; i++) {
      int proc = recv_info[recv_index[i]].first ;
      recv_ptr[recv_index[i]] = new unsigned char[maxr_size[recv_index[i]]] ;
      MPI_Irecv(recv_ptr[recv_index[i]], maxr_size[recv_index[i]], MPI_PACKED, proc, 2, MPI_COMM_WORLD, &re_request[i]) ;
    }
    
    for(int i=0;i<resend_size;++i) {
      int loc_pack = 0 ;
      send_ptr[send_index[i]] = new unsigned char[maxs_size[send_index[i]]] ;
      for(size_t j=0;j<send_info[send_index[i]].second.size();++j) {
	storeRepP sp = send_vars[i][j] ;//facts.get_variable(send_info[send_index[i]].second[j].v) ;
	sp->pack(send_ptr[send_index[i]], loc_pack,maxs_size[send_index[i]],send_info[send_index[i]].second[j].set);
      }
    }
    
    // Send Buffer
    for(int i=0;i<resend_size;++i) {
      int proc = send_info[send_index[i]].first ;
      MPI_Send(send_ptr[send_index[i]],maxs_size[send_index[i]],MPI_PACKED,proc,2,MPI_COMM_WORLD) ;
      delete [] send_ptr[send_index[i]] ;
    }
    if(rerecv_size > 0) {
#ifdef DEBUG
      int err =
#endif
        MPI_Waitall(rerecv_size, re_request, re_status) ;
      FATAL(err != MPI_SUCCESS) ;
    }
    for(int i=0;i<rerecv_size;++i) {
      int loc_unpack = 0;
      for(size_t j=0;j<recv_info[recv_index[i]].second.size();++j) {
	storeRepP sp = recv_vars[i][j] ;//facts.get_variable(recv_info[recv_index[i]].second[j].v) ;
	storeRepP sr = sp->new_store(EMPTY) ;
	sr->allocate(entitySet(recv_info[recv_index[i]].second[j].seq)) ;
	
	sr->unpack(recv_ptr[recv_index[i]], loc_unpack, maxr_size[recv_index[i]],
		   recv_info[recv_index[i]].second[j].seq) ;
	
	CPTR<joiner> op = join_op->clone() ;
	op->SetArgs(sp,sr) ;
	op->Join(recv_info[recv_index[i]].second[j].seq) ;
      }
      delete [] recv_ptr[recv_index[i]] ;
    }
    if(rerecv_size > 0) {
      delete [] re_status ;
      delete [] re_request ;
    }
  }
  
  void execute_comm_reduce::Print(ostream &s) const {
    int sz = 0 ;
    if(send_info.size()+recv_info.size() > 0) {
      s << "reduction block {" << endl ;
      if(send_info.size() > 0) {
        s << "Send:" << endl ;
        for(size_t i=0;i<send_info.size();++i) {
          for(size_t j=0;j<send_info[i].second.size();++j) {
            s << send_info[i].second[j].v << "  " ;
	    sz += (send_info[i].second[j].set).size() ;
	  }
	  s << " to " << send_info[i].first << endl ;
        }
	s << " Total entities sent = " << sz << endl ;	
      }
      sz = 0 ;
      if(recv_info.size() > 0) {
        s << "Recv:" << endl ;
        for(size_t i=0;i<recv_info.size();++i) {
          for(size_t j=0;j<recv_info[i].second.size();++j) {
            s << recv_info[i].second[j].v << "  " ;
	    sz += (recv_info[i].second[j].seq).size() ;
	  }
          s << " from " << recv_info[i].first << endl ;
        }
	s << " Total entities recieved = " << sz << endl ;	
      }
      s << "}" << endl ;
    }
  }
  
  executeP reduce_store_compiler::create_execution_schedule(fact_db &facts, sched_db &scheds) {
    if(facts.isDistributed()) {
      CPTR<execute_sequence> el = new execute_sequence ;
      
      el->append_list(new execute_comm_reduce(rlist, facts, join_op)) ;
      el->append_list(new execute_comm(clist, facts)) ;
      ostringstream oss ;
      oss << "reduce store " << reduce_var ;
      el->append_list(new execute_msg(oss.str())) ;

      return executeP(el) ;
    }

    ostringstream oss ;
    oss << "reduce store " << reduce_var ;
    return executeP(new execute_msg(oss.str())) ;
  }
}
