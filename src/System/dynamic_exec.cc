//
// Iterative weighted static load balancing and dynamic loop
//   scheduling rule for Loci
// by RL Carino (RLC@HPC.MsState.Edu)
//
// Ordinarily, Loci executes 
//      rp->compute (sequence (exec_set));
// on each processor; however, the exec_sets may have different 
// sizes, or the amount of computations required by the items 
// may be non-uniform. Therefore, the processor loads may be 
// imbalanced, leading to application performance degradation.
//
// This rule adds a load balancing functionality into Loci, 
// hopefully to increase performance when the distribution 
// of items among processors result in load imbalance. Two
// strategies are implemented: an iterative weighted static
// load balancing strategy, and a dynamic loop scheduling 
// strategy.  Descriptions of these are given in the "NOTES"
// below. 
//
// This version assumes that exec_set.num_intervals()==1
//

#include "sched_tools.h"
#include "distribute.h"

namespace Loci
{

// max no. of procs
#define MAXPROCS 16

// set to 1 to activate debugging statements for:
#define DEBUGOUT 0		// load balancing messages
#define SHOWTIMES 1		// summary timings

// static methods
#define NLB		0	// no load balancing
#define IWS		1	// iterative weighted static load balancing

// loop scheduling methods
#define FSC		2	// fixed size chunking (same no. of chunks as FAC)
#define FAC		3	// factoring
#define GSS		4	// guided self scheduling

//=======================================================
// NOTES for iterative weighted static load balancing (IWS)
//=======================================================
//
//  static load balancing - proc loads are assigned BEFORE any load is executed
//  weighted - proc work times are used as WEIGHTs to predetermine loads
//  iterative - work times are GEOMETRIC-AVERAGED over calls to the module
//
//  items = elements of exec_set
//  chunk = a group of contiguous items
//  chunk time = execution time of chunk (from MPI_Wtime())
//  work time = total execution time of items owned by a proc,
//      some of which may have been executed by other procs
//
// ------------------------------------------
// IWS algorithm overview
// ------------------------------------------
// Step 0. Initializations
//      Set nCalls = 0
//      Compute a common size for chunking items (chunkSize)
// Step 1. Load balance
//      If (nCalls <= SKIPSTEPS) then
// Step 2. Set this proc to neither send nor receive work
//      Else
// Step 3. Based on work times, compute an ideal time (optTime)
// Step 4. Determine amounts of work transfers from work times & optTime
// Step 5. If this proc is a source of work, identify chunks to send
//      End if
//
//      If this proc is a receiver of work, 
// Step 6. Receive chunks of items from others into local_facts
// Step 7. Execute chunks of items in local_facts; collect chunk times
// Step 8. Execute items in facts; collect chunk times & proc work time
// Step 9. Return outputs of items in local_facts, chunk times
//
//      Else if this proc is a sender of work
// Step 9. Send chunks
// Step 10. Execute remaining items; collect chunk times & proc work time
// Step 11. Receive outputs for items sent, chunk times
//
//      Else // this proc is neither sender nor receiver of work
// Step 12. Execute items in facts; collect chunk times & proc work time
//      End if
//
// Step 13. Exchange proc work times; increment nCalls
// Step 14. Compute averages of chunk times and work times
//
// On the next call, start at Step 1.
//
// ------------------------------------------

//=======================================================
// #define's for IWS
//=======================================================

// set to 1 to activate debugging statements for:
#define SHOWLB1 0		// calculations for ideal loads
#define SHOWLB2 0		// estimated loads
#define SHOWCHUNKS 0		// timings for transferred chunks

// skip load balancing for the first few time steps
// use timings collected here to initialize average time
#define SKIPSTEPS 5

// weight/contribution of latest time to average time
#define CONTRIB 0.5

// average no. of chunks per proc
#define AVECHUNKS 50
#define MAXCHUNKS 2*AVECHUNKS

//  skip the JSTART most expensive chunks when deciding transfers
#define JSTART 1

// "negligible" time relative to optimal time
#define tolLB  0.075

// message tags
#define TAG_ACK      1		// handshake acknowledgement
#define TAG_INFO     2		// size of succeeding message
#define TAG_INPUT    3		// input data
#define TAG_OUTPUT   4		// output data

// info[i][] index; info regarding chunk of items to be migrated
#define MSG_LEN      0		// packing size
#define ITEM_COUNT   1		// number of items
#define LOCAL_START  2		// index in facts/local_facts


//=======================================================
// NOTES for dynamic loop scheduling (DLS)
//=======================================================
// 
// Load balancing by DLS is based on a foreman-worker strategy.
// The workers execute loads specified by the foreman, and the
// foreman, in addition to executing loads, is responsible for
// determining the loads and for detecting termination. A loop
// scheduling technique is utilized by the foreman to compute
// loads. A load is basically a chunk of items, which may be
// local to a worker, or migrated from another worker. The
// following "#define"s symbolize the actions that can be taken
// by a processor.

//=======================================================
// #define's for DLS
//=======================================================

// actions
#define WAIT4_MSG          1
#define TEST4_MSG          2
#define RETRIEVE_MSG       3
#define WORK_LOCAL         4
#define WORK_REMOTE        5
#define SEND_INPUT         6
#define RECV_OUTPUT        7
#define FILL_REQUEST       8
#define QUIT               9
#define EXEC_LOCAL_PART   10
#define EXEC_REMOTE_PART  11
#define EXEC_REMOTE_WHOLE 12

// The actions FILL_REQUEST and EXEC_REMOTE_PART are unique to
// the foreman; EXEC_REMOTE_WHOLE is unique to workers; and
// the rest of the actions are common to the foreman and to
// workers.

// =======================================================
// Outline of algorithm
// =======================================================

// =======================================================
// Initial loads
// =======================================================
// if (i_am_foreMan)
//   Send initial local work info to workers
//   Set own local work info, nextAction = WORK_LOCAL
// else // i_am_worker
//   Set nextAction = WAIT4_MSG
// 
// =======================================================
// Main loop
// =======================================================
// while (all termination conditions not met)
// 
//   switch (nextAction)
// 
// =======================================================
//     case WAIT4_MSG
//       MPI_Probe (MPI_ANY_SOURCE, MPI_ANY_TAG, ...)
//       Set nextAction = RETRIEVE_MSG, lastAction = WAIT4_MSG
// =======================================================
// 
// =======================================================
//     case TEST4_MSG
//       MPI_Iprobe (MPI_ANY_SOURCE, MPI_ANY_TAG, ...)
//       if (there was a message) Set nextAction = RETRIEVE_MSG
//       else Set nextAction = lastAction
// =======================================================
// 
// =======================================================
//     case RETRIEVE_MSG
//       Set nextAction from message status.MPI_TAG (one of the #define)
//       Retrieve content of message, which may be
//         - signal to quit
//         - info on local work to execute
//         - info on input items to migrate
//         - a request for load & timing data (for foreman)
//         - input items
//         - output items & timing data
// =======================================================
// 
// =======================================================
//     case QUIT
//       Set one of the termination flags, nextAction = WAIT4_MSG
// =======================================================
// 
// =======================================================
//     case WORK_LOCAL
//       Update local work info
//       Compute probing frequency, when to request next work info
//       Set nextAction = EXEC_LOCAL_PART
// =======================================================
// 
// =======================================================
//     case WORK_REMOTE
//       Update remote work info
//       if (i_am_foreman) Set nextAction = EXEC_REMOTE_PART
//       else Set nextAction = EXEC_REMOTE_WHOLE
// =======================================================
// 
// =======================================================
//     case SEND_INPUT
//       Send input according to info retrieved from message
//       Set nextAction = lastAction
// =======================================================
// 
// =======================================================
//     case RECV_OUTPUT
//       Decrement counter for outputs to receive
//       Set nextAction = lastAction
// =======================================================
// 
// =======================================================
//     case EXEC_LOCAL_PART
//       If this is the last subchunk, then send request for next work info 
//       Execute a subchunk from current local work
//       Update counters and accumulators
//       if (current local chunk not finished)
//         Set nextAction = TEST4_MSG, lastAction = EXEC_LOCAL_PART
//       else if (i_am_foreman) // foreman finished current local chunk 
//         // pretend foreman sent request (to itself)
//         Set lastAction = WAIT4_MSG, nextAction = FILL_REQUEST
//       else Set nextAction = WAIT4_MSG
// =======================================================
// 
// =======================================================
//     case EXEC_REMOTE_PART
//       Foreman executes subchunk of remote items
//       if (chunk not finished)
//         Set nextAction = TEST4_MSG, lastAction = EXEC_REMOTE_PART
//       else // foreman finished current remote chunk 
//         Return outputs, Set lastAction = WAIT4_MSG
//         // pretend foreman sent a request to self
//         Set nextAction = FILL_REQUEST 
// =======================================================
// 
// =======================================================
//     case EXEC_REMOTE_WHOLE
//       Worker executes remote chunk, return the outputs, 
//         requests for new work info
//       Set nextAction = WAIT4_MSG
// =======================================================
// 
// =======================================================
//     case FILL_REQUEST
//       Foreman takes one of the following actions,
//         - inform requesting proc that there are no more items to execute
//         - send new local work info to requesting proc
//         - inform a heavily loaded worker to send items to requesting proc
//       Foreman checks for termination
// =======================================================
// 
//   end switch 
// end while 
// 
// MPI_Barrier (MPI_COMM_WORLD)
// 


//=======================================================
// additional #define's for IWS
//=======================================================

// send request for next chunk when the
// remaining iterations in current chunk reaches PCT_SEND_REQUEST % 
// subject to the minimum of MIN_SEND_REQUEST iterations
#define PCT_SEND_REQUEST 0.10
#define MIN_SEND_REQUEST 50

// IF THE REQUEST FOR NEXT CHUNK HAS NOT BEEN SENT,
// probe for messages after executing
// PCT_PROBE_FREQ % iterations in current chunk 
// subject to the minimum MIN_PROBE_FREQ iterations
#define PCT_PROBE_FREQ 0.05
#define MIN_PROBE_FREQ 10

  //=======================================================
  // common global variables
  //=======================================================

  rule_implP rp1;		// copies, to facilitate shared code
  rule_implP local_comp1;	//  between IWS & loop scheduling methods
  fact_db *facts1;
  fact_db *local_facts1;
  variableSet inputs1;
  variableSet outputs1;

  unsigned char *buf = 0;	// message buffer
  int buf_size = 0;		// size of message buffer

  int nCalls;			// no. of calls to execute()
  int allItems;			// total  no. of items
  int nProcs;			// no. of procs
  int myRank;			// process rank
  int myItemCount;		// local item count
  int myFirstItem;		// first local item

  double ideal1=0.0;		// perfect balance time
  double local1=0.0;		// time spent by this proc on own items
  double remote1=0.0;		// execution time of migrated items
  double wall1=0.0;		// local1 + remote1 + load balancing overheads
  double wall2=0.0;		// wall1 + any post-processing (i.e., MPI_Allreduce() )

  double workTime[MAXPROCS];	// execution time of items belonging to proc i
  double aveWorkTime[MAXPROCS];	// average work time in proc i
  double ewt[MAXPROCS];		// expected work time (EWT) of proc i

  //=======================================================
  // used by IWS
  //=======================================================

  int iwsSize;			// IWS chunk size
  int myChunks;			// local chunk count
  int doneBy[MAXCHUNKS];	// proc which will execute chunk
  double chunkTime[MAXPROCS][MAXCHUNKS];	// latest chunk execution time
  double aveChunkTime[MAXCHUNKS];	// average chunk execution time

  //=======================================================
  // used by DLS
  //=======================================================

  int foreMan;			// rank of loop scheduler 
  int yMapSave[2 * MAXPROCS];	// copy of item count, first item

//=======================================================
// common routines
//=======================================================

  void GetBuffer (int size)
  {
    if (size > buf_size)
      {
	if (buf != 0)
	  delete[]buf;
	buf_size = size + size / 5;
	buf = new unsigned char[buf_size];
      }
  }

  void AllocateLBspace (int nItems)
  {
    for (variableSet::const_iterator vi = inputs1.begin ();
	 vi != inputs1.end (); ++vi)
      {
	storeRepP sp = local_facts1->get_variable (*vi);
	sp->allocate (interval (0, nItems - 1));
      }
    for (variableSet::const_iterator vi = outputs1.begin ();
	 vi != outputs1.end (); ++vi)
      {
	storeRepP sp = local_facts1->get_variable (*vi);
	sp->allocate (interval (0, nItems - 1));
      }
  }

  void FreeLBspace ()
  {
    for (variableSet::const_iterator vi = inputs1.begin ();
	 vi != inputs1.end (); ++vi)
      {
	storeRepP sp = local_facts1->get_variable (*vi);
	sp->allocate (EMPTY);
      }
    for (variableSet::const_iterator vi = outputs1.begin ();
	 vi != outputs1.end (); ++vi)
      {
	storeRepP sp = local_facts1->get_variable (*vi);
	sp->allocate (EMPTY);
      }
  }


  // packing size of a chunk of inputs in facts
  int inputPackSize (int tStart, int tSize) {
    int size;
    size = 0;
    for (variableSet::const_iterator vi = inputs1.begin ();
	 vi != inputs1.end (); ++vi)
      {
	storeRepP s_ptr = facts1->get_variable (*vi);
	size += s_ptr->pack_size (interval (tStart, tStart + tSize - 1));
      }
    return size;
  }


  // packing size of a chunk of outputs in local_facts
  int outputPackSize (int tStart, int tSize) {
    int size;
    size = 0;
    for (variableSet::const_iterator vi = outputs1.begin ();
	 vi != outputs1.end (); ++vi)
      {
	storeRepP s_ptr = local_facts1->get_variable (*vi);
	size += s_ptr->pack_size (interval (tStart, tStart+tSize - 1));
      }
    return size;
  }

  //
  // rule initializer
  //
  dynamic_schedule_rule::dynamic_schedule_rule (rule fi, entitySet eset,
						fact_db & facts,
						sched_db & scheds)
  {
    int i2[2];

    // Loci initializations
    rp = fi.get_rule_implP ();
    rule_tag = fi;
    exec_set = eset;
    local_compute1 = rp->new_rule_impl ();
    entitySet in = rule_tag.sources ();
    outputs = rule_tag.targets ();

    //Setup local facts input variables(types only no allocation)
    for (variableSet::const_iterator vi = in.begin (); vi != in.end (); ++vi)
      {
	storeRepP store_ptr = rp->get_store (*vi);
	if ((store_ptr != 0) && store_ptr->RepType () == Loci::STORE)
	  {
	    inputs += *vi;
	    local_facts.create_fact (*vi, store_ptr->new_store (EMPTY));
	  }
	else
	  {
	    local_facts.create_fact (*vi, facts.get_variable (*vi));
	  }
      }

    //Setup local facts output variables
    for (variableSet::const_iterator vi = outputs.begin ();
	 vi != outputs.end (); ++vi)
      {
	storeRepP store_ptr = rp->get_store (*vi);
	local_facts.create_fact (*vi, store_ptr->new_store (EMPTY));
      }

    //Initialize both functions for remote and local execution.
    local_compute1->initialize (local_facts);
    rp->initialize (facts);

// ========================================================================
// Step 0. Initializations
// ========================================================================

    // proc count, rank 
    nProcs = Loci::MPI_processes;
    myRank = Loci::MPI_rank;
    foreMan = 0;

    // first items, item counts
    myFirstItem = exec_set.Min ();
    myItemCount = exec_set.size ();

    //Total work
    i2[0] = myFirstItem;
    i2[1] = myItemCount;
    MPI_Allgather (i2, 2, MPI_INT, yMapSave, 2, MPI_INT, MPI_COMM_WORLD);
    allItems = 0;
    for (int proc = 0; proc < nProcs; proc++)
      allItems = allItems + yMapSave[2 * proc + 1];

    // IWS chunk size, local chunk count
    iwsSize = (allItems + nProcs * AVECHUNKS - 1) / (nProcs * AVECHUNKS);
    myChunks = (myItemCount + iwsSize - 1) / iwsSize;
    for (int chunkIdx = 0; chunkIdx < MAXCHUNKS; chunkIdx++)
      {
	aveChunkTime[chunkIdx] = 0.0;
	doneBy[chunkIdx] = myRank;
	for (int proc = 0; proc < MAXPROCS; proc++)
	  chunkTime[proc][chunkIdx] = 0.0;
      }
    // initialize work times
    for (int proc = 0; proc < nProcs; proc++)
      {
	workTime[proc] = 0.0;
	aveWorkTime[proc] = 0.0;
      }

    // number of calls to execute()
    nCalls = 0;

#if DEBUGOUT
    Loci::debugout << myRank << " has " << myFirstItem << "," << myItemCount
      << "," << myChunks << "," << iwsSize << endl;
#endif

  }

  dynamic_schedule_rule::~dynamic_schedule_rule ()
  {
    if (buf != 0)
      {
	delete[]buf;
	buf = 0;
	buf_size = 0;
      }
  }

  void iterative_weighted_static ()
  {

    int sortIdx[MAXCHUNKS];	// sort index
    int info[MAXPROCS][3];	// information about items received from other procs

    double optTime = 0.0;	// "ideal" work time (perfect load balance)
    double local0 = 0.0;	// predicted work time for (remaining) local items
    double remote0 = 0.0;	// predicted work time for (sent or received) items
    double toGive[MAXPROCS];	// work to transfer to proc i
    double x[MAXPROCS][MAXPROCS];	// total chunk time to move from proc i to proc j

    int recvFrom, sendTo;
    int nChunks, nItems, msgLen, maxLen;
    int tSize, tStart, position, tIdx;
    int tBuf[2];		// message buffer for 2 ints
    MPI_Status tStatus;

    int i, j, k, dest, src, proc;
    int chunkIdx, jmax, jmin;
    double toFill, tTime;
    double timerStart, timerEnd = 0.0;

    // BARRIER
    // MPI_Barrier (MPI_COMM_WORLD);

// ========================================================================
// Step 1. Load balance
// ========================================================================

// ========================================================================
// Step 2. Set this proc to neither send nor receive work
// ========================================================================

    sendTo = 0;			// how many procs to send work to
    recvFrom = 0;		// how many procs to recv work from
    local1 = 0.0;
    remote1 = 0.0;
    for (proc = 0; proc < nProcs; proc++)
      {
	workTime[proc] = 0.0;	// work times
	info[proc][MSG_LEN] = 0;	// buffer size
	info[proc][ITEM_COUNT] = 0;	// no. of items
	info[proc][LOCAL_START] = 0;	// position of items in local_facts
      }

    // do load balancing only after a few steps
    if (nCalls > SKIPSTEPS)
      {

// ========================================================================
// Step 3. Based on work times, compute an ideal time (optTime)
// ========================================================================

	// set chunk ownership
	for (k = 0; k < myChunks; k++)
	  {
#if SHOWCHUNKS
	    if (doneBy[k] != myRank)
	      {
		i = myFirstItem + k * iwsSize;
		j = min (iwsSize, myFirstItem + myItemCount - i);
		Loci::
		  debugout << nCalls << "," << myRank << "," << k << "," << i
		  << "," << j << "," << chunkTime[myRank][k] << ", done by "
		  << doneBy[k] << endl;
	      }
#endif
	    doneBy[k] = myRank;
	  }

	// ideal time, copy of work times (ewt[])
	optTime = 0.0;
	for (i = 0; i < nProcs; i++)
	  {
	    ewt[i] = aveWorkTime[i];
	    optTime += aveWorkTime[i];
	  }
	optTime /= nProcs;
	local0 = aveWorkTime[myRank];	// predicted work time from local items
	remote0 = 0.0;		// predicted work time for/from remote items

// ========================================================================
// Step 4. Determine amounts of work transfers from work times & optTime
// ========================================================================

	// work time to give away
	for (i = 0; i < nProcs; i++)
	  if (abs (ewt[i] - optTime) / optTime >= tolLB)
	    toGive[i] = ewt[i] - optTime;
	  else
	    toGive[i] = 0.0;
#if SHOWLB1
	Loci::
	  debugout << "ST5: rank=" << i << ", expected=" << ewt[myRank]
	  << ", ideal=" << optTime << ", diff=" << toGive[myRank] << endl;
#endif

	// sort work to give away (increasing)
	for (i = 0; i < nProcs; i++)
	  sortIdx[i] = i;
	for (i = 0; i < (nProcs - 1); i++)
	  for (j = i + 1; j < nProcs; j++)
	    if (toGive[sortIdx[i]] > toGive[sortIdx[j]])
	      {
		k = sortIdx[i];
		sortIdx[i] = sortIdx[j];
		sortIdx[j] = k;
	      }
	// determine communications pattern
	// x[src][dest] is the amount of work from src to dest
	for (i = 0; i < nProcs; i++)
	  for (j = 0; j < nProcs; j++)
	    x[i][j] = 0.0;
	for (k = 0; k < nProcs; k++)
	  {
	    src = sortIdx[k];
	    if (toGive[src] > 0.0)
	      {
		// reduction in src is significant?
		while (toGive[src] / optTime > tolLB)
		  {
		    // find destination
		    dest = k;
		    for (j = 0; j < k; j++)
		      if (toGive[sortIdx[dest]] > toGive[sortIdx[j]])
			dest = j;
		    dest = sortIdx[dest];
		    // destination is already full ?
		    if (abs (toGive[dest]) / optTime <= tolLB)
		      break;
		    toFill = -toGive[dest];
		    // not enough from src ?
		    if (toGive[src] < toFill)
		      toFill = toGive[src];
		    // update accumulators
		    x[src][dest] = toFill;
		    ewt[dest] += toFill;
		    toGive[dest] += toFill;
		    ewt[src] -= toFill;
		    toGive[src] -= toFill;
#if SHOWLB1
		    Loci::
		      debugout << "ST5 xfer : " << src << " -> " << dest <<
		      " : " << toFill << endl;
#endif
		  }
	      }
	  }

	// how many procs to receive from?
	recvFrom = 0;
	for (src = 0; src < nProcs; src++)
	  if (x[src][myRank] > 0.0)
	    {
#if SHOWLB2
	      // receive estimated work from src
	      MPI_Recv (&tTime, 1, MPI_DOUBLE, src, TAG_INPUT,
			MPI_COMM_WORLD, &tStatus);
	      Loci::
		debugout << "ST5: work by " << myRank << " from " << src
		<< " : est. by self=" << x[src][myRank]
		<< ", est. from src=" << tTime
		<< ", %diff=" << 100.0 * (x[src][myRank] -
					  tTime) / x[src][myRank] << endl;
	      remote0 += tTime;
#endif
	      recvFrom++;
	    }

	// how many procs to send to?
	sendTo = 0;
	for (dest = 0; dest < nProcs; dest++)
	  if (x[myRank][dest] > 0.0)
	    sendTo++;

      }				// end if (nCalls > SKIPSTEPS)

// ========================================================================
// Step 5. If this proc is a source of work, identify chunks to send
// ========================================================================

    if (sendTo)
      {

	// reuse EWT for computed amount to transfer to others
	for (i = 0; i < nProcs; i++)
	  ewt[i] = 0.0;
	local0 = aveWorkTime[myRank];

	for (dest = 0; dest < nProcs; dest++)
	  // send something to dest?
	  if (x[myRank][dest] > 0.0)
	    {
	      // sort chunks accd to chunk time
	      jmax = -1;
	      for (i = 0; i < (myChunks - 1); i++)
		if (doneBy[i] == myRank)
		  {
		    jmax++;
		    sortIdx[jmax] = i;
		  }
	      // exclude last chunk (might be smaller than iwsSize)
	      for (i = 0; i < (jmax - 1); i++)
		for (j = i + 1; j < jmax; j++)
		  if (aveChunkTime[sortIdx[i]] < aveChunkTime[sortIdx[j]])
		    {
		      k = sortIdx[i];
		      sortIdx[i] = sortIdx[j];
		      sortIdx[j] = k;
		    }

	      // collect chunks for dest
	      toFill = x[myRank][dest];
	      jmin = JSTART;
	      tTime = aveChunkTime[sortIdx[jmin]];

	      while ((jmin < (jmax - 1)) && (toFill / optTime > tolLB))
		{
		  // toFill needs more than 1 chunk?
		  while ((jmin < (jmax - 1)) && (toFill >= tTime))
		    {
#if SHOWLB2
		      Loci::
			debugout << "ST6: " << myRank << " -> " << dest <<
			" : " << toFill << " - " << tTime <<
			" from chunk " << jmin << " (=" << sortIdx[jmin]
			<< "), rem is " << toFill - tTime << endl;
#endif
		      toFill -= tTime;
		      local0 -= tTime;
		      ewt[dest] += tTime;
		      doneBy[sortIdx[jmin]] = dest;	// mark this chunk as for dest
		      jmin++;
		      tTime = aveChunkTime[sortIdx[jmin]];
		    }
		  if (toFill / optTime <= tolLB)
		    break;
		  // tTime greater than toFill
		  while ((jmin < (jmax - 1)) && (tTime > toFill))
		    {
		      jmin++;
		      tTime = aveChunkTime[sortIdx[jmin]];
		    }
		}
#if SHOWLB2
	      // send estimate of work to dest
	      MPI_Send (&ewt[dest], 1, MPI_DOUBLE, dest, TAG_INPUT,
			MPI_COMM_WORLD);
	      remote0 += ewt[dest];
#endif
	    }
      }

#if SHOWLB2
    Loci::
      debugout << "step=" << nCalls << ",rank=" << myRank
      << "  Predictions: ideal=" << optTime
      << ", local=" << local0
      << ", remote=" << remote0 << ", total=" << local0 + remote0 << endl;
#endif

    // Any proc to receive from ? (Is this rank lightly loaded?)
    if (recvFrom)
      {

// ========================================================================
// Step 6. Receive chunks of items from others into local_facts
// ========================================================================

	// receive info on incoming items
	maxLen = 0;
	nItems = 0;		// total items to receive
	for (proc = 0; proc < recvFrom; proc++)
	  {
	    MPI_Recv (tBuf, 2, MPI_INT, MPI_ANY_SOURCE, TAG_INFO,
		      MPI_COMM_WORLD, &tStatus);
	    src = tStatus.MPI_SOURCE;
	    tSize = tBuf[0];	// number of items
	    msgLen = tBuf[1];	// packed size
#if DEBUGOUT
	    Loci::
	      debugout << "DEST1: input info from " << src << " : " <<
	      tSize << "," << msgLen << endl;
#endif
	    // store into info[][]
	    info[src][MSG_LEN] = msgLen;
	    info[src][ITEM_COUNT] = tSize;
	    if (maxLen < msgLen)
	      maxLen = msgLen;
	    nItems += tSize;
	  }
	// reset recvFrom in case tSize=0 was specified
	recvFrom = 0;
	for (src = 0; src < nProcs; src++)
	  if (info[src][ITEM_COUNT])
	    recvFrom++;

	// receiving from any other proc still?
	if (recvFrom)
	  {
#if DEBUGOUT
	    Loci::
	      debugout << "DEST1: local_facts to accommodate " << nItems <<
	      endl;
#endif
	    // allocate space for items to receive
	    AllocateLBspace (nItems);

	    // prepare buffer, ACKnowledge handshake, receive & unpack inputs
	    GetBuffer (maxLen);
	    tStart = 0;		// item numbers in local_facts start with 0
	    for (src = 0; src < nProcs; src++)
	      {
		tSize = info[src][ITEM_COUNT];
		if (tSize == 0)
		  continue;

		// send acknowledgement
		MPI_Send (NULL, 0, MPI_INT, src, TAG_ACK, MPI_COMM_WORLD);
		msgLen = info[src][MSG_LEN];
		info[src][LOCAL_START] = tStart;
		// receive items
		MPI_Recv (buf, msgLen, MPI_PACKED, src, TAG_INPUT,
			  MPI_COMM_WORLD, &tStatus);
		nChunks = tSize / iwsSize;
		// unpack items into local_facts
		position = 0;
		for (chunkIdx = 0; chunkIdx < nChunks; chunkIdx++)
		  {
		    for (variableSet::const_iterator vi = inputs1.begin ();
			 vi != inputs1.end (); ++vi)
		      {
			storeRepP s_ptr = local_facts1->get_variable (*vi);
			s_ptr->unpack (buf, position, msgLen,
				       sequence (interval
						 (tStart,
						  tStart + iwsSize - 1)));
		      }
#if SHOWCHUNKS
		    Loci::
		      debugout << "DEST2: unpacked input chunk " << chunkIdx
		      << " from " << src << ", start at " << tStart << endl;
#endif
		    tStart += iwsSize;
		  }
	      }

// ========================================================================
// Step 7. Execute chunks of items in local_facts; collect chunk times
// ========================================================================

	    // execute items in local_facts
	    for (src = 0; src < nProcs; src++)
	      {
		tSize = info[src][ITEM_COUNT];
		if (tSize == 0)
		  continue;

		// number of chunks from proc src in local_facts
		nChunks = tSize / iwsSize;
		timerStart = MPI_Wtime ();
		workTime[src] = timerStart;
		timerEnd = timerStart;
		// execute chunks
		for (chunkIdx = 0; chunkIdx < nChunks; chunkIdx++)
		  {
		    tStart = info[src][LOCAL_START] + chunkIdx * iwsSize;
		    local_comp1->
		      compute (sequence
			       (interval (tStart, tStart + iwsSize - 1)));
		    timerEnd = MPI_Wtime ();
		    chunkTime[src][chunkIdx] = (timerEnd - timerStart);
#if DEBUGOUT+SHOWCHUNKS
		    Loci::
		      debugout << "DEST3: executed chunk " << chunkIdx <<
		      " from " << src << ", start at " << tStart << " in "
		      << chunkTime[src][chunkIdx] << endl;
#endif
		    timerStart = timerEnd;
		  }
		workTime[src] = timerEnd - workTime[src];
		remote1 += workTime[src];
#if DEBUGOUT
		Loci::
		  debugout << "DEST3: local_compute->compute() time for " <<
		  src << " is " << workTime[src] << endl;
#endif
	      }

	  }			// end if (recvFrom)

// ========================================================================
// Step 8. Execute items in facts; collect chunk times & proc work time
// ========================================================================

	// execute items in facts, by chunks
	timerStart = MPI_Wtime ();
	local1 = timerStart;
	timerEnd = timerStart;
	for (chunkIdx = 0; chunkIdx < myChunks; chunkIdx++)
	  {
	    tStart = myFirstItem + chunkIdx * iwsSize;
	    tSize = min (iwsSize, myFirstItem + myItemCount - tStart);
	    rp1->compute (sequence (interval (tStart, tStart + tSize - 1)));
	    timerEnd = MPI_Wtime ();
	    chunkTime[myRank][chunkIdx] = (timerEnd - timerStart);
	    timerStart = timerEnd;
	  }
	local1 = timerEnd - local1;
	workTime[myRank] = local1;
#if DEBUGOUT
	Loci::
	  debugout << "DEST4: rp->compute() time=" << local1 << endl;
#endif

	if (recvFrom)
	  {

// ========================================================================
// Step 9. Return outputs of items in local_facts, chunk times
// ========================================================================

	    // calculate message sizes, maximum size, initiate handshake by sending chunk times
	    maxLen = 0;
	    for (dest = 0; dest < nProcs; dest++)
	      {
		tSize = info[dest][ITEM_COUNT];
		if (tSize == 0)
		  continue;

		// calculate pack size
		nChunks = tSize / iwsSize;
		msgLen = 0;
		for (chunkIdx = 0; chunkIdx < nChunks; chunkIdx++)
		  {
		    tStart = info[dest][LOCAL_START] + chunkIdx * iwsSize;
		    msgLen += outputPackSize (tStart, iwsSize);
		//    for (variableSet::const_iterator vi = outputs1.begin ();
		//	 vi != outputs1.end (); ++vi)
		//      {
		//	storeRepP s_ptr = local_facts1->get_variable (*vi);
		//	msgLen +=
		//	  s_ptr->
		//	  pack_size (interval (tStart, tStart + iwsSize - 1));
		//      }
		  }
		if (maxLen < msgLen)
		  maxLen = msgLen;
		info[dest][MSG_LEN] = msgLen;
		// place pack size at end of message
		chunkTime[dest][nChunks] = double (msgLen);
		MPI_Send (chunkTime[dest], nChunks + 1, MPI_DOUBLE,
			  dest, TAG_INFO, MPI_COMM_WORLD);
#if DEBUGOUT
		Loci::
		  debugout << "DEST5: output info to " << dest << " : " <<
		  tSize << "," << msgLen << ", first at " <<
		  info[dest][LOCAL_START] << endl;
#endif
	      }
	    // prepare sufficient buffer
	    GetBuffer (maxLen);
	    // receive ACK, pack and send output
	    for (proc = 0; proc < recvFrom; proc++)
	      {
		MPI_Recv (NULL, 0, MPI_INT, MPI_ANY_SOURCE, TAG_ACK,
			  MPI_COMM_WORLD, &tStatus);
		dest = tStatus.MPI_SOURCE;
		msgLen = info[dest][MSG_LEN];
		nChunks = info[dest][ITEM_COUNT] / iwsSize;
		position = 0;
		for (chunkIdx = 0; chunkIdx < nChunks; chunkIdx++)
		  {
		    tStart = info[dest][LOCAL_START] + chunkIdx * iwsSize;
		    for (variableSet::const_iterator vi = outputs1.begin ();
			 vi != outputs1.end (); ++vi)
		      {
			storeRepP s_ptr = local_facts1->get_variable (*vi);
			s_ptr->pack (buf, position, msgLen,
				     interval (tStart, tStart + iwsSize - 1));
		      }
		  }
		MPI_Send (buf, msgLen, MPI_PACKED, dest, TAG_OUTPUT,
			  MPI_COMM_WORLD);
#if DEBUGOUT
		Loci::
		  debugout << "DEST6: sent outputs to " << dest << " : " <<
		  msgLen << "," << info[dest][ITEM_COUNT] << ", start=" <<
		  info[dest][LOCAL_START] << endl;
#endif
	      }

	    //Free space after LB
	    FreeLBspace ();

	  }			// end if (recvFrom)

      }

    // Send to other procs ? (Is this proc heavily loaded?)
    else if (sendTo)
      {

// ========================================================================
// Step 9. Send chunks
// ========================================================================

	// calculate maximum message size and initiate handshake
	maxLen = 0;
	for (dest = 0; dest < nProcs; dest++)
	  {
	    if (x[myRank][dest] == 0.0)
	      continue;		// nothing for dest
	    msgLen = 0;
	    nItems = 0;
	    for (chunkIdx = 0; chunkIdx < myChunks; chunkIdx++)
	      {
		if (doneBy[chunkIdx] != dest)
		  continue;	// not for dest
		tStart = myFirstItem + chunkIdx * iwsSize;
		nItems += iwsSize;
		msgLen += inputPackSize (tStart, iwsSize);

		//for (variableSet::const_iterator vi = inputs1.begin ();
		//     vi != inputs1.end (); ++vi)
		//  {
		//    storeRepP s_ptr = facts1->get_variable (*vi);
		//    msgLen +=
		//      s_ptr->
		//      pack_size (interval (tStart, tStart + iwsSize - 1));
		//  }
	      }
	    tBuf[0] = nItems;
	    tBuf[1] = msgLen;
	    MPI_Send (tBuf, 2, MPI_INT, dest, TAG_INFO, MPI_COMM_WORLD);
#if DEBUGOUT
	    Loci::
	      debugout << "SRC1: sent input info to " << dest << " : " <<
	      nItems << "," << msgLen << endl;
#endif
	    if (maxLen < msgLen)
	      maxLen = msgLen;
	    info[dest][MSG_LEN] = msgLen;
	    info[dest][ITEM_COUNT] = nItems;

	    // exclude dest in case a chunk was not specified for it 
	    if (nItems == 0)
	      sendTo--;
	  }

	// Send to other procs still?
	if (sendTo)
	  {

	    // prepare sufficient buffer
	    GetBuffer (maxLen);
	    // receive ACK, pack and send input
	    for (proc = 0; proc < sendTo; proc++)
	      {
		MPI_Recv (NULL, 0, MPI_INT, MPI_ANY_SOURCE, TAG_ACK,
			  MPI_COMM_WORLD, &tStatus);
		dest = tStatus.MPI_SOURCE;
		tSize = info[dest][ITEM_COUNT];
		msgLen = info[dest][MSG_LEN];
		//Pack inputs from facts
		position = 0;
		for (chunkIdx = 0; chunkIdx < myChunks; chunkIdx++)
		  {
		    if (doneBy[chunkIdx] != dest)
		      continue;	// not for dest
		    tStart = myFirstItem + chunkIdx * iwsSize;
		    for (variableSet::const_iterator vi = inputs1.begin ();
			 vi != inputs1.end (); ++vi)
		      {
			storeRepP s_ptr = facts1->get_variable (*vi);
			s_ptr->pack (buf, position, msgLen,
				     interval (tStart, tStart + iwsSize - 1));
		      }
		  }
		MPI_Rsend (buf, msgLen, MPI_PACKED, dest, TAG_INPUT,
			   MPI_COMM_WORLD);
#if DEBUGOUT
		Loci::debugout << "SRC2: sent inputs to " << dest << " : "
		  << info[dest][ITEM_COUNT] << "," << msgLen << endl;
#endif
	      }

	  }			//end if (sendTo)

// ========================================================================
// Step 10. Execute remaining items; collect chunk times & proc work time
// ========================================================================

	// execute remaining inputs
	timerStart = MPI_Wtime ();
	local1 = timerStart;
	timerEnd = timerStart;
	for (chunkIdx = 0; chunkIdx < myChunks; chunkIdx++)
	  {
	    if (doneBy[chunkIdx] != myRank)
	      continue;		// sent to another proc
	    tStart = myFirstItem + chunkIdx * iwsSize;
	    tSize = min (iwsSize, myFirstItem + myItemCount - tStart);
	    rp1->compute (sequence (interval (tStart, tStart + tSize - 1)));
	    timerEnd = MPI_Wtime ();
	    chunkTime[myRank][chunkIdx] = (timerEnd - timerStart);
	    timerStart = timerEnd;
	  }
	local1 = timerEnd - local1;
	workTime[myRank] = local1;
#if DEBUGOUT
	Loci::
	  debugout << "SRC4: rp->compute() time=" << local1 << endl;
#endif

	if (sendTo)
	  {

// ========================================================================
// Step 11. Receive outputs for items sent, chunk times
// ========================================================================

	    // receive info on returning outputs
	    maxLen = 0;
	    for (proc = 0; proc < sendTo; proc++)
	      {
		MPI_Probe (MPI_ANY_SOURCE, TAG_INFO, MPI_COMM_WORLD,
			   &tStatus);
		src = tStatus.MPI_SOURCE;
		MPI_Get_count (&tStatus, MPI_DOUBLE, &tSize);
		MPI_Recv (chunkTime[src], tSize, MPI_DOUBLE, src,
			  TAG_INFO, MPI_COMM_WORLD, &tStatus);
		// unpack msgLen from last position
		msgLen = int (chunkTime[src][tSize - 1]);
		info[src][MSG_LEN] = msgLen;
		if (maxLen < msgLen)
		  maxLen = msgLen;
#if DEBUGOUT
		Loci::
		  debugout << "SRC5: recv output info from " << src <<
		  " : " << iwsSize * (tSize -
				      1) << "," << msgLen << " (expecting "
		  << info[src][ITEM_COUNT] << " items)" << endl;
#endif
	      }
	    // prepare buffer, ACKnowledge handshake, receive & unpack inputs
	    GetBuffer (maxLen);
	    for (src = 0; src < nProcs; src++)
	      {
		tSize = info[src][ITEM_COUNT];
		if (tSize == 0)
		  continue;

		MPI_Send (NULL, 0, MPI_INT, src, TAG_ACK, MPI_COMM_WORLD);
		msgLen = info[src][MSG_LEN];
		MPI_Recv (buf, msgLen, MPI_PACKED, src, TAG_OUTPUT,
			  MPI_COMM_WORLD, &tStatus);
		//unpack outputs into facts
		tIdx = 0;	// index to chunk times
		position = 0;
		for (chunkIdx = 0; chunkIdx < myChunks; chunkIdx++)
		  {
		    if (doneBy[chunkIdx] != src)
		      continue;	// not for src
		    tStart = myFirstItem + chunkIdx * iwsSize;
		    for (variableSet::const_iterator vi = outputs1.begin ();
			 vi != outputs1.end (); ++vi)
		      {
			storeRepP s_ptr = facts1->get_variable (*vi);
			s_ptr->unpack (buf, position, msgLen,
				       interval (tStart,
						 tStart + iwsSize - 1));
		      }
		    chunkTime[myRank][chunkIdx] = chunkTime[src][tIdx];
		    remote1 += chunkTime[src][tIdx];
#if DEBUGOUT+SHOWCHUNKS
		    Loci::
		      debugout << "SRC6: unpacked output for chunk " <<
		      chunkIdx << ", first at " << tStart << ", completed in "
		      << chunkTime[src][tIdx] << endl;
#endif
		    tIdx++;
		  }
#if DEBUGOUT
		Loci::
		  debugout << "SRC6: unpacked output from " << src << " : "
		  << tSize << "," << msgLen << endl;
#endif
	      }

	  }			//end if (sendTo)
      }


// ========================================================================
// Step 12. Execute items in facts; collect chunk times & proc work time
// ========================================================================

    // no items to send or receive; simply execute items in facts
    else
      {
	timerStart = MPI_Wtime ();
	local1 = timerStart;
	timerEnd = timerStart;
	for (chunkIdx = 0; chunkIdx < myChunks; chunkIdx++)
	  {
	    tStart = myFirstItem + chunkIdx * iwsSize;
	    tSize = min (iwsSize, myFirstItem + myItemCount - tStart);
	    rp1->compute (sequence (interval (tStart, tStart + tSize - 1)));
	    timerEnd = MPI_Wtime ();
	    chunkTime[myRank][chunkIdx] = (timerEnd - timerStart);
	    timerStart = timerEnd;
	  }
	local1 = timerEnd - local1;
	workTime[myRank] = local1;
#if DEBUGOUT
	Loci::debugout << "NLB4: rp->compute() time=" << local1 << endl;
#endif
      }

#if SHOWTIMES
    wall1 = MPI_Wtime () - wall1;
#endif

// ========================================================================
// Step 13. Exchange proc work times
// ========================================================================

    // BARRIER: gather work times on all procs
    MPI_Allreduce (workTime, ewt, nProcs, MPI_DOUBLE, MPI_SUM,
		   MPI_COMM_WORLD);
    for (proc = 0; proc < nProcs; proc++)
      workTime[proc] = ewt[proc];

// ========================================================================
// Step 14. Compute averages of chunk times and work times
// ========================================================================

    // compute geometric average work time for each proc
    if (nCalls <= SKIPSTEPS)
      {
	optTime = workTime[myRank];	// give optTime a value
	for (proc = 0; proc < nProcs; proc++)
	  aveWorkTime[proc] += workTime[proc];
	for (chunkIdx = 0; chunkIdx < myChunks; chunkIdx++)
	  aveChunkTime[chunkIdx] += chunkTime[myRank][chunkIdx];
	if (nCalls == SKIPSTEPS)
	  {
	    for (proc = 0; proc < nProcs; proc++)
	      aveWorkTime[proc] /= (nCalls + 1);
	    for (chunkIdx = 0; chunkIdx < myChunks; chunkIdx++)
	      aveChunkTime[chunkIdx] /= (nCalls + 1);
	  }
#if DEBUGOUT
	Loci::
	  debugout << "SKIP: average worktime=" << aveWorkTime[myRank] <<
	  endl;
#endif
      }
    else
      {
	for (proc = 0; proc < nProcs; proc++)
	  aveWorkTime[proc] =
	    (1.0 - CONTRIB) * aveWorkTime[proc] + CONTRIB * workTime[proc];
	for (chunkIdx = 0; chunkIdx < myChunks; chunkIdx++)
	  aveChunkTime[chunkIdx] = (1.0 - CONTRIB) * aveChunkTime[chunkIdx] +
	    CONTRIB * chunkTime[myRank][chunkIdx];
#if DEBUGOUT
	Loci::
	  debugout << "CONT: average worktime=" << aveWorkTime[myRank] <<
	  endl;
#endif
      }
#if SHOWLB2
    if (nCalls > SKIPSTEPS)
      Loci::
	debugout << "step=" << nCalls << ",rank=" << myRank
	<< "  Errors (%) : ideal="
	<< 100.0 * (optTime - wall1) / optTime
	<< ", local=" << 100.0 * (local0 - local1) /
	((local0 > 0.0) ? local0 : 1.0)
	<< ", remote=" << 100.0 * (remote0 - remote1) /
	((remote0 > 0.0) ? remote0 : 1.0) << endl;
#endif

    ideal1 = optTime;
  }



//=======================================================
// routines called by DLS
//=======================================================


  void
    SendInfo (int dest, int action, int chunkStart, int chunkSize,
	      int chunkDest, double mu)
  {
    int bufInfo[3], pSize, tPos;
    pSize = 3 * sizeof (int) + sizeof (double);
    GetBuffer (pSize);
    bufInfo[0] = chunkStart;
    bufInfo[1] = chunkSize;
    bufInfo[2] = chunkDest;
    tPos = 0;
    MPI_Pack (bufInfo, 3, MPI_INT, buf, pSize, &tPos, MPI_COMM_WORLD);
    MPI_Pack (&mu, 1, MPI_DOUBLE, buf, pSize, &tPos, MPI_COMM_WORLD);
    MPI_Send (buf, pSize, MPI_PACKED, dest, action, MPI_COMM_WORLD);
#if DEBUGOUT
    //if (myRank == foreMan)
      Loci::
	debugout << "Info: action=" << action << ", " << myRank << "->" << dest 
	<< "->" << chunkDest
	<< " : " << chunkStart << "," << chunkSize 
	<< endl;
#endif
  }

  void
    RecvInfo (int src, int action, int *chunkStart, int *chunkSize,
	      int *chunkDest, double *mu)
  {
    int bufInfo[3], pSize, tPos;
    MPI_Status tStatus;

    pSize = 3 * sizeof (int) + sizeof (double);
    GetBuffer (pSize);
    MPI_Recv (buf, pSize, MPI_PACKED, src, action, MPI_COMM_WORLD, &tStatus);
    tPos = 0;
    MPI_Unpack (buf, pSize, &tPos, bufInfo, 3, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack (buf, pSize, &tPos, mu, 1, MPI_DOUBLE, MPI_COMM_WORLD);
    *chunkStart = bufInfo[0];
    *chunkSize = bufInfo[1];
    *chunkDest = bufInfo[2];
  }

  void SendInput (int dest, int tStart, int tSize)
  {
    int position, size=0, pSize, tPos;
    int bufInfo[3];

    //compute buffer size
//    size = 0;
//    for (variableSet::const_iterator vi = inputs1.begin ();
//	 vi != inputs1.end (); ++vi)
//      {
//	storeRepP s_ptr = facts1->get_variable (*vi);
//	size += s_ptr->pack_size (interval (tStart, tStart + tSize - 1));
//      }
//    pSize = size + 12;		// 3 ints = 12 bytes
    pSize = 12 + inputPackSize (tStart, tSize);
    GetBuffer (pSize);
    //Pack chunk info
    bufInfo[0] = tStart;
    bufInfo[1] = tSize;
    bufInfo[2] = size;
    tPos = 0;
    MPI_Pack (bufInfo, 3, MPI_INT, buf, pSize, &tPos, MPI_COMM_WORLD);
    //Pack inputs from facts
    position = 0;
    entitySet myent = interval (tStart, tStart + tSize - 1);
    for (variableSet::const_iterator vi = inputs1.begin ();
	 vi != inputs1.end (); ++vi)
      {
	storeRepP s_ptr = facts1->get_variable (*vi);
	s_ptr->pack (&buf[tPos], position, size, myent);
      }
    MPI_Send (buf, pSize, MPI_PACKED, dest, WORK_REMOTE, MPI_COMM_WORLD);

#if DEBUGOUT
    Loci::
      debugout << "SendInput " << myRank << "->" << dest << ", chunk=" <<
      tStart << "," << tSize << "," << size << endl;
#endif
  }

  void ReceiveInput (int src, int msgSize, int *tStart, int *tSize)
  {
    int position, tPos, size;
    MPI_Status tStatus;
    int bufInfo[3];

    GetBuffer (msgSize);
    MPI_Recv (buf, msgSize, MPI_PACKED, src, WORK_REMOTE, MPI_COMM_WORLD,
	      &tStatus);
    //unpack chunk info
    tPos = 0;
    MPI_Unpack (buf, msgSize, &tPos, bufInfo, 3, MPI_INT, MPI_COMM_WORLD);
    *tStart = bufInfo[0];
    *tSize = bufInfo[1];
    size = bufInfo[2];
    //unpack inputs into local facts
    position = 0;
    for (variableSet::const_iterator vi = inputs1.begin ();
	 vi != inputs1.end (); ++vi)
      {
	storeRepP s_ptr = local_facts1->get_variable (*vi);
	s_ptr->unpack (&buf[tPos], position, size,
		       sequence (interval (0, *tSize - 1)));
      }

#if DEBUGOUT
    Loci::
      debugout << "ReceiveInput " << myRank << "<-" << src << ", chunk=" <<
      *tStart << "," << *tSize << "," << size << endl;
#endif
  }

  void SendOutput (int dest, int tStart, int tSize, double *tTime)
  {
    int position, size=0, tPos, pSize;
    int bufInfo[3];

    //compute buffer size
//    size = 0;
//    for (variableSet::const_iterator vi = outputs1.begin ();
//	 vi != outputs1.end (); ++vi)
//      {
//	storeRepP s_ptr = local_facts1->get_variable (*vi);
//	size += s_ptr->pack_size (interval (0, tSize - 1));
//      }
//    pSize = size + 3 * sizeof (int) + sizeof (double);
    pSize = 3 * sizeof (int) + sizeof (double) + outputPackSize (0, tSize);
    GetBuffer (pSize);
    //Pack chunk info
    bufInfo[0] = tStart;
    bufInfo[1] = tSize;
    bufInfo[2] = size;
    tPos = 0;
    MPI_Pack (bufInfo, 3, MPI_INT, buf, pSize, &tPos, MPI_COMM_WORLD);
    MPI_Pack (tTime, 1, MPI_DOUBLE, buf, pSize, &tPos, MPI_COMM_WORLD);
    //Pack outputs
    position = 0;
    entitySet myent2 = interval (0, tSize - 1);
    for (variableSet::const_iterator vi = outputs1.begin ();
	 vi != outputs1.end (); ++vi)
      {
	storeRepP s_ptr = local_facts1->get_variable (*vi);
	s_ptr->pack (&buf[tPos], position, size, myent2);
      }
    //Send outputs
    MPI_Send (buf, pSize, MPI_PACKED, dest, RECV_OUTPUT, MPI_COMM_WORLD);

#if DEBUGOUT
    Loci::
      debugout << "SendOutput " << myRank << "->" << dest << ", chunk=" <<
      tStart << "," << tSize << "," << size << endl;
#endif
  }

  void ReceiveOutput (int src, int msgSize, int *iters, double *tTime)
  {
    MPI_Status tStatus;
    int position, tPos, tStart, tSize, size;
    int bufInfo[3];

    GetBuffer (msgSize);
    MPI_Recv (buf, msgSize, MPI_PACKED, src, RECV_OUTPUT, MPI_COMM_WORLD,
	      &tStatus);
    //unpack chunk info
    tPos = 0;
    MPI_Unpack (buf, msgSize, &tPos, bufInfo, 3, MPI_INT, MPI_COMM_WORLD);
    MPI_Unpack (buf, msgSize, &tPos, tTime, 1, MPI_DOUBLE, MPI_COMM_WORLD);
    tStart = bufInfo[0];
    tSize = bufInfo[1];
    size = bufInfo[2];
    *iters = tSize;
    //unpack outputs into facts
    position = 0;
    for (variableSet::const_iterator vi = outputs1.begin ();
	 vi != outputs1.end (); ++vi)
      {
	storeRepP s_ptr = facts1->get_variable (*vi);
	s_ptr->unpack (&buf[tPos], position, size,
		       sequence (interval (tStart, tStart + tSize - 1)));
      }

#if DEBUGOUT
    Loci::
      debugout << "ReceiveOutput " << myRank << "<-" << src << ", chunk=" <<
      tStart << "," << tSize << "," << size << endl;
#endif
  }


  void
    GetChunkSize (int method, int minChunkSize, int source,
		  int *yMap, int *chunkSize, int *batchSize, int *batchRem)
  {
    int i, tChunk, rem;

    rem = 0;
    for (i = 0; i < nProcs; i++)
      rem += yMap[2 * i + 1];

    switch (method)
      {

      case FSC:		//fixed size scheduling, P * log(N / P) chunks
	i = (allItems + nProcs - 1) / nProcs;
	tChunk = int ((0.55 + i * log (2.0e0) / log ((double) i)));
	tChunk = min (tChunk, rem);
	*batchSize = tChunk;
	*batchRem = min (*batchSize, rem);
	break;

      case GSS:		//guided self scheduling
	tChunk = max (minChunkSize, (rem + nProcs - 1) / nProcs);
	tChunk = min (rem, tChunk);
	*batchSize = tChunk;
	*batchRem = min (*batchSize, rem);
	break;

      case FAC:		//factoring
	if (*batchRem == 0)
	  {
	    tChunk =
	      max (minChunkSize, (rem + 2 * nProcs - 1) / (2 * nProcs));
	    *batchSize = nProcs * tChunk;
	    *batchRem = min (*batchSize, rem);
	  }
	// use current batchSize 
	tChunk = max (minChunkSize, *batchSize / nProcs);
	tChunk = min (rem, tChunk);
	break;

      default:
	//no scheduling
	tChunk = yMap[2 * source + 1];
	*batchSize = tChunk;
	*batchRem = min (*batchSize, rem);
	break;

      }

    //adjust according to Remaining[source]
    if (yMap[2 * source + 1] > (tChunk + minChunkSize / 2))
      *chunkSize = tChunk;
    else
      *chunkSize = yMap[2 * source + 1];

    //adjust remaining in batch
    *batchRem -= *chunkSize;
    if (*batchRem < minChunkSize)
      *batchRem = 0;

    //adjust counters for source
    yMap[2 * source] += *chunkSize;
    yMap[2 * source + 1] -= *chunkSize;
  }


  void loop_scheduling (int method)
  {
    int yMap[2 * MAXPROCS];	// working copy of yMapSave[]

    MPI_Status mStatus;
    int action = 0, lastAction = 0;
    int chunkSize = 0, chunkStart = 0, chunkDest = 0;
    int batchSize = 0, batchRem = 0, minChunkSize = 0;
    int inputSent[MAXPROCS];

    int localStart, localSize;
    int remoteStart, remoteSize;
    int saveStart, saveSize, saveSrc;

    int numIdle;
    int returns;		//output data
    int incoming;		//input data
    int gotWork;

    int msgSrc = 0, msgLen = 0;
    int tStart = 0, tSize = 0, tSource = 0;
    int MsgInQueue = 0;
    int sendRequest, probeFreq;

    double tTime, currentChunkTime, latestFinish, execTime;
    double timerStart, timerDiff;
    //double idleTime;          // time in MPI_Probe

    int i;

    local1 = 0.0;
    remote1 = 0.0;
    currentChunkTime = 0.0;
    for (i = 0; i < nProcs; i++)
    {
      workTime[i] = 0.0;
      ewt[i] = aveWorkTime[i];
    }

    // allocate space for LB 
    tSize = (allItems + 4 * nProcs - 1) / (4 * nProcs);
    AllocateLBspace (tSize);

    returns = 0;
    incoming = 0;
    numIdle = 0;
    probeFreq = MIN_PROBE_FREQ;
    sendRequest = MIN_SEND_REQUEST;
    gotWork = 1;

    localSize = 0;
    localStart = 0;

    remoteSize = 0;
    remoteStart = 0;
    saveSrc = 0;
    saveSize = 0;
    saveStart = 0;

// =======================================================
// Initial loads
// =======================================================
// if (i_am_foreMan)
//   Send initial local work info to workers
//   Set own local work info, nextAction = WORK_LOCAL
// else // i_am_worker
//   Set nextAction = WAIT4_MSG

    if (myRank == foreMan)
      {
	//scheduler initializations
	batchSize = 0;
	batchRem = 0;
	minChunkSize = 2 * max (MIN_PROBE_FREQ, MIN_SEND_REQUEST);
	for (i = 0; i < (2 * nProcs); i++)
	  yMap[i] = yMapSave[i];

	// use static scheduling during first call to execute()
	//   to initialize average work time
	if (nCalls)
	  i = method;
	else
	  i = -1;
	//send local work info to others
	for (tSource = 0; tSource < nProcs; tSource++)
	  if (tSource != myRank)
	    {
	      inputSent[tSource] = 0;
	      tStart = yMap[2 * tSource];
	      GetChunkSize (i, minChunkSize, tSource,
			    yMap, &tSize, &batchSize, &batchRem);
	      SendInfo (tSource, WORK_LOCAL, tStart, tSize, -1, 0.0);
	    }
	//foreman collects work for itself
	inputSent[myRank] = 0;
	chunkStart = yMap[2 * myRank];
	GetChunkSize (i, minChunkSize, myRank,
		      yMap, &chunkSize, &batchSize, &batchRem);
	action = WORK_LOCAL;
      }

    else
      {
	//workers block for first message
	action = WAIT4_MSG;
      }
// 
// =======================================================
// Main loop
// =======================================================

    while (gotWork + localSize + remoteSize + incoming + returns)
      {

	switch (action)
	  {

// =======================================================
//     case WAIT4_MSG
//       MPI_Probe (MPI_ANY_SOURCE, MPI_ANY_TAG, ...)
//       Set nextAction = RETRIEVE_MSG, lastAction = WAIT4_MSG
// =======================================================

	  case WAIT4_MSG:	// wait for a message
	    //timerStart = MPI_Wtime ();
#if DEBUGOUT
	    Loci::
	      debugout << "MPI_Probe() ..." << endl; 
#endif
	    MPI_Probe (MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &mStatus);
	    //timerDiff = MPI_Wtime () - timerStart;
	    //idleTime += timerDiff;
	    action = RETRIEVE_MSG;
	    lastAction = WAIT4_MSG;
	    break;

// =======================================================
//     case TEST4_MSG
//       MPI_Iprobe (MPI_ANY_SOURCE, MPI_ANY_TAG, ...)
//       if (there was a message) Set nextAction = RETRIEVE_MSG
//       else Set nextAction = lastAction
// =======================================================

	  case TEST4_MSG:	// test for a message
	    MPI_Iprobe (MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD,
			&MsgInQueue, &mStatus);
	    if (MsgInQueue)
	      action = RETRIEVE_MSG;
	    else
	      action = lastAction;
	    break;

// =======================================================
//     case RETRIEVE_MSG
//       Set nextAction from message status.MPI_TAG (one of the #define)
//       Retrieve content of message, which may be
//         - signal to quit
//         - info on local work to execute
//         - info on input items to migrate
//         - a request for load & timing data (for foreman)
//         - input items
//         - output items & timing data
// =======================================================

	  case RETRIEVE_MSG:	// retrieve message contents according to tag
	    msgSrc = mStatus.MPI_SOURCE;
	    action = mStatus.MPI_TAG;
	    if ((action == QUIT) || (action == WORK_LOCAL) ||
		(action == SEND_INPUT) || (action == FILL_REQUEST))
	      {
		RecvInfo (msgSrc, action, &chunkStart, &chunkSize, &chunkDest,
			  &execTime);
	      }
	    else
	      {			// message contains MPI_PACKED data
		MPI_Get_count (&mStatus, MPI_PACKED, &msgLen);
		if (action == WORK_REMOTE)
		  ReceiveInput (msgSrc, msgLen, &chunkStart, &chunkSize);
		else		// if (action==RECV_OUTPUT)
		  {
		    ReceiveOutput (msgSrc, msgLen, &i, &tTime);
		    remote1 += tTime;
		  }

	      }
	    break;

// =======================================================
//     case QUIT
//       Set one of the termination flags, nextAction = WAIT4_MSG
// =======================================================

	  case QUIT:		// no more work, but there may be messages
	    gotWork = 0;
	    action = WAIT4_MSG;
	    break;

// =======================================================
//     case WORK_LOCAL
//       Update local work info
//       Compute probing frequency, when to request next work info
//       Set nextAction = EXEC_LOCAL_PART
// =======================================================

	  case WORK_LOCAL:	// update local work info
	    localStart = chunkStart;
	    localSize = chunkSize;
	    currentChunkTime = 0.0;
	    if (myRank != foreMan)
	      {
		probeFreq =
		  max (MIN_PROBE_FREQ, int (localSize * PCT_PROBE_FREQ));
		sendRequest =
		  max (MIN_SEND_REQUEST, int (localSize * PCT_SEND_REQUEST));
		if (returns)
		  probeFreq = MIN_PROBE_FREQ;
	      }
	    action = EXEC_LOCAL_PART;	// execute local subchunk
	    break;

// =======================================================
//     case WORK_REMOTE
//       Update remote work info
//       if (i_am_foreman) Set nextAction = EXEC_REMOTE_PART
//       else Set nextAction = EXEC_REMOTE_WHOLE
// =======================================================

	  case WORK_REMOTE:	// receive remote input data
	    saveStart = chunkStart;
	    saveSize = chunkSize;
	    saveSrc = msgSrc;
	    remoteStart = 0;
	    remoteSize = chunkSize;
	    currentChunkTime = 0.0;
	    if (myRank == foreMan)
	      action = EXEC_REMOTE_PART;	// subchunk
	    else
	      action = EXEC_REMOTE_WHOLE;	// whole chunk
	    break;

// =======================================================
//     case SEND_INPUT
//       Send input according to info retrieved from message
//       Set nextAction = lastAction
// =======================================================

	  case SEND_INPUT:	// send input data
	    SendInput (chunkDest, chunkStart, chunkSize);
	    returns++;
	    probeFreq = min (probeFreq, MIN_PROBE_FREQ);
	    action = lastAction;
	    break;

// =======================================================
//     case RECV_OUTPUT
//       Decrement counter for outputs to receive
//       Set nextAction = lastAction
// =======================================================

	  case RECV_OUTPUT:	// receive output
	    returns--;
	    action = lastAction;
	    break;

// =======================================================
//     case EXEC_LOCAL_PART
//       If this is the last subchunk, then send request for next work info 
//       Execute a subchunk from current local work
//       Update counters and accumulators
//       if (current local chunk not finished)
//         Set nextAction = TEST4_MSG, lastAction = EXEC_LOCAL_PART
//       else if (i_am_foreman) // foreman finished current local chunk 
//         // pretend foreman sent request (to itself)
//         Set lastAction = WAIT4_MSG, nextAction = FILL_REQUEST
//       else Set nextAction = WAIT4_MSG
// =======================================================

	  case EXEC_LOCAL_PART:	// execute a local subchunk
	    // Not yet the last subchunk?
	    if (localSize > (sendRequest + probeFreq))
	      tSize = probeFreq;	// no
	    else
	      tSize = localSize;	// yes
	    //Send request before doing last subchunk
	    if ((localSize == tSize) && (myRank != foreMan))
	      {
		SendInfo (foreMan, FILL_REQUEST, 0, 0, myRank, currentChunkTime);
	      }

	    timerStart = MPI_Wtime ();
	    rp1->
	      compute (sequence
		       (interval (localStart, localStart + tSize - 1)));
	    timerDiff = MPI_Wtime () - timerStart;
	    workTime[myRank] += timerDiff;

	    localStart += tSize;
	    localSize -= tSize;
	    currentChunkTime += timerDiff;
	    if (localSize)
	      {			// current local chunk not finished
		action = TEST4_MSG;	// test for a message
		lastAction = EXEC_LOCAL_PART;
	      }
	    else if (myRank == foreMan)
	      {			// foreman finished current local chunk 
		// pretend foreman sent request (to itself)
		lastAction = WAIT4_MSG;	// foreman was waiting for a message
		msgSrc = myRank;	// foreman received a message (from self)
		ewt[myRank] -= currentChunkTime;	// decrement remaining work time
		action = FILL_REQUEST;	// message was a request 
	      }
	    else
	      {
		action = WAIT4_MSG;	// wait for a message
	      }
	    break;

// =======================================================
//     case EXEC_REMOTE_PART
//       Foreman executes subchunk of remote items
//       if (chunk not finished)
//         Set nextAction = TEST4_MSG, lastAction = EXEC_REMOTE_PART
//       else // foreman finished current remote chunk 
//         Return outputs, Set lastAction = WAIT4_MSG
//         // pretend foreman sent a request to self
//         Set nextAction = FILL_REQUEST 
// =======================================================

	  case EXEC_REMOTE_PART:	// foreman executes remote subchunk 
	    tSize = min (remoteSize, MIN_PROBE_FREQ);

	    timerStart = MPI_Wtime ();
	    local_comp1->
	      compute (sequence
		       (interval (remoteStart, remoteStart + tSize - 1)));
	    timerDiff = MPI_Wtime () - timerStart;
	    currentChunkTime += timerDiff;

	    remoteStart += tSize;
	    remoteSize -= tSize;
	    if (remoteSize)
	      {			// current remote chunk not finished
		action = TEST4_MSG;	// test for a message
		lastAction = EXEC_REMOTE_PART;
	      }
	    else
	      {			// foreman finished current remote chunk 
		SendOutput (saveSrc, saveStart, saveSize, &currentChunkTime);
		workTime[saveSrc] += currentChunkTime;
		remote1 += currentChunkTime;
		incoming--;	// decrement counter for incoming data
		// pretend foreman sent request (to itself)
		lastAction = WAIT4_MSG;	// foreman waited for a message
		msgSrc = myRank;	// foreman received a message (from self)
		action = FILL_REQUEST;	// message was a request 
		// RecvInfo (msgSrc, action, &chunkStart, &chunkSize, &chunkDest,
		// 	  &execTime);
		chunkStart = 0;
		chunkSize  = 0;
		chunkDest  = saveSrc;
		execTime = currentChunkTime;
	      }
	    break;

// =======================================================
//     case EXEC_REMOTE_WHOLE
//       Worker executes remote chunk, return the outputs, 
//         requests for new work info
//       Set nextAction = WAIT4_MSG
// =======================================================

	  case EXEC_REMOTE_WHOLE:	// execute remote chunk, send output, request
	    timerStart = MPI_Wtime ();
	    local_comp1->compute (sequence (interval (0, remoteSize - 1)));
	    currentChunkTime = MPI_Wtime () - timerStart;
	    workTime[saveSrc] += currentChunkTime;
	    remote1 += currentChunkTime;

	    SendOutput (saveSrc, saveStart, saveSize, &currentChunkTime);
	    SendInfo (foreMan, FILL_REQUEST, 0, 0, saveSrc, currentChunkTime);
	    action = WAIT4_MSG;	// wait for a message
	    remoteSize = 0;
	    break;

// =======================================================
//     case FILL_REQUEST
//       Foreman takes one of the following actions,
//         - inform requesting proc that there are no more items to execute
//         - send new local work info to requesting proc
//         - inform a heavily loaded worker to send items to requesting proc
// =======================================================

	  case FILL_REQUEST:	// received a request (foreman only)
	    ewt[chunkDest] -= execTime;	// decrement remaining work time of chunk owner
	    // find a data source
	    if (yMap[2 * msgSrc + 1] > 0)
	      tSource = msgSrc;	// requester has unfinished data
	    else if (inputSent[msgSrc])	// requester has none
	      tSource = -1;	// prohibit requester from accepting work since it gave away data
	    else
	      {			// find worker expected to finish last
		tSource = -1;
		latestFinish = 0.0;
		for (i = 0; i < nProcs; i++)
		  if (yMap[2*i+1] && (ewt[i]>latestFinish))
		    {
			tSource = i;
			latestFinish = ewt[i];
		    }
		// Uncomment/comment next line to disable/enable load balancing
		//tSource = -1;
	      }

	    // who's requesting ?
	    if (msgSrc == foreMan)
	      {
		chunkStart = 0;
		chunkSize = 0;
		chunkDest = 0;
		if (tSource == -1)
		  {		// no more unscheduled iterates
		    action = WAIT4_MSG;	// wait for a message
		  }		// end if (tSource == -1)
		else if (tSource == msgSrc)
		  {		// msgSrc has some
		    chunkStart = yMap[2 * msgSrc];
		    GetChunkSize (method, minChunkSize, msgSrc,
				  yMap, &chunkSize, &batchSize, &batchRem);
		    action = WORK_LOCAL;	// update local work info
#if DEBUGOUT
		    Loci::
		      debugout << "Info0:" << myRank << "->" << tSource <<
		      "->-1 : " << chunkStart << "," << chunkSize << endl;
#endif
		  }		//end if (tSource == msgSrc)
		else
		  {		// tSource sends input data to msgSrc
		    tStart = yMap[2 * tSource];
		    GetChunkSize (method, minChunkSize, tSource,
				  yMap, &tSize, &batchSize, &batchRem);
		    SendInfo (tSource, SEND_INPUT, tStart, tSize, msgSrc,
			      0.0);
		    inputSent[tSource]++;	// mark tSource as data giver
		    action = WAIT4_MSG;	// wait for a message
		    incoming++;	// increment counter for incoming work
		  }
	      }			// end if (msgSrc==foreMan)
	    else
	      {			// request was from a worker
		if (tSource == -1)
		  {		// no more unscheduled iterates
		    numIdle++;
		    SendInfo (msgSrc, QUIT, 0, 0, 0, 0.0);
		    gotWork = numIdle != (nProcs - 1);
		    action = TEST4_MSG;	// test for a message
		  }		// end if (tSource == -1)
		else if (tSource == msgSrc)
		  {		// msgSrc has some
		    tStart = yMap[2 * msgSrc];
		    GetChunkSize (method, minChunkSize, msgSrc,
				  yMap, &tSize, &batchSize, &batchRem);
		    SendInfo (msgSrc, WORK_LOCAL, tStart, tSize, -1, 0.0);
		    action = TEST4_MSG;	// test for a message
		  }		//end if (tSource == msgSrc)
		else
		  {		// tSource sends input data to msgSrc
		    tStart = yMap[2 * tSource];
		    GetChunkSize (method, minChunkSize, tSource,
				  yMap, &tSize, &batchSize, &batchRem);
		    inputSent[tSource]++;	// mark tSource as data giver
		    if (tSource != foreMan)
		      {
			SendInfo (tSource, SEND_INPUT, tStart, tSize, msgSrc,
				  0.0);
			action = TEST4_MSG;	// test for a message
		      }
		    else
		      {
			// pretend foreman sent request (to itself)
			// pretend foreman waited for a message
			// pretend foreman received a message (from self)
			action = SEND_INPUT;	// message was to send input
			chunkStart = tStart;
			chunkSize = tSize;
			chunkDest = msgSrc;
#if DEBUGOUT
			Loci::
			  debugout << "Info0:" << myRank << "->" << tSource <<
			  "->" << chunkDest << " : " << chunkStart << "," << chunkSize 
			  << endl;
#endif
		      }
		  }
	      }			// end else (msgSrc!=foreMan)
	    break;

	  }			// end switch (action)
      }				// end while (gotWork + localSize + remoteSize + returns)

// 
// MPI_Barrier (MPI_COMM_WORLD)
// 

#if DEBUGOUT
    Loci::debugout << "Exits: " << myRank << endl;
#endif
    FreeLBspace ();

#if SHOWTIMES
    wall1 = MPI_Wtime () - wall1;
    local1 = workTime[myRank];

    // BARRIER: gather work times on all procs
    for (i = 0; i < nProcs; i++)
      ewt[i] = 0.0;
    MPI_Allreduce (workTime, ewt, nProcs, MPI_DOUBLE, MPI_SUM,
		   MPI_COMM_WORLD);

    // running average of work times
    if (nCalls)
      for (i = 0; i < nProcs; i++)
	aveWorkTime[i] = (1.0-CONTRIB )*aveWorkTime[i] + CONTRIB * ewt[i];
    else
      for (i = 0; i < nProcs; i++)
	aveWorkTime[i] = ewt[i];

    // perfectly balance work time
    ideal1 = 0.0;
    for (i = 0; i < nProcs; i++)
      ideal1 += ewt[i];
    ideal1 /= nProcs;
#else
    // simple BARRIER
    MPI_Barrier (MPI_COMM_WORLD);
#endif
  }


//=======================================================
// The execute() routine
//=======================================================

  void dynamic_schedule_rule::execute (fact_db & facts)
  {

    extern int method;

#if SHOWTIMES
    wall1 = MPI_Wtime ();
    wall2 = wall1;
    if (method == NLB)
      {
	rp->compute (sequence (exec_set));
	wall1 = MPI_Wtime () - wall1;
	local1 = wall1;
	ideal1 = wall1;
	wall2 = wall1;
      }
    else
      {
	// globals
	rp1 = rp;
	local_comp1 = local_compute1;
	facts1 = &facts;
	local_facts1 = &local_facts;
	inputs1 = inputs;
	outputs1 = outputs;

	if (method == IWS)
	  iterative_weighted_static ();
	else
	  loop_scheduling (method);

	// ideal1, local1, remote1 & wall1 are computed 
	//   by the routines
	wall2 = MPI_Wtime () - wall2;

	facts1 = 0;
	local_facts1 = 0;
      }
    Loci::debugout << method << "," << nProcs << "," << myRank << "," << nCalls
      << ", ideal=" << ideal1
      << ", wall1=" << wall1
      << " (" << int (1000.0 * (ideal1 - wall1) / ideal1) / 10.0 << "%)"
      << ", wall2=" << wall2
      << " (" << int (1000.0 * (ideal1 - wall2) / ideal1) / 10.0 << "%)"
      << ", local=" << local1
      << ", remote=" << remote1
      << ", ovrhd1=" << int (1000.0 * (wall1 - local1 - remote1) / wall1) /
      10.0 << "%" << endl;

#else

    if (method == NLB)
      rp->compute (sequence (exec_set));
    else
      {
	rp1 = rp;
	local_comp1 = local_compute1;
	facts1 = &facts;
	local_facts1 = &local_facts;
	inputs1 = inputs;
	outputs1 = outputs;
	if (method == IWS)
	    iterative_weighted_static ();
	else
	    loop_scheduling (method);
	facts1 = 0;
	local_facts1 = 0;
      }
#endif
    nCalls++;
  }


  void dynamic_schedule_rule::Print (std::ostream & s) const
  {
    s << "dynamic schedule " << rule_tag << "  over set " << exec_set << endl;
  }


}
