#include <stdio.h>
#include <math.h>
#include "mpi.h"
#include <Loci.h>
#include <iostream>
using std::cout ;
using std::cerr ; 
using std::endl ;

namespace Loci {
  
#define min(a, b)       ((a) < (b) ? (a) : (b))
#define max(a, b)       ((a) > (b) ? (a) : (b))
#define MAX_PROCESSES 128
#define MAX_CHUNKS 100
  /* message tags */
#define WKP_MSG    9990     /* I'm alive! W->F */
#define WRK_MSG    9980     /* Do chunk with these specs, F->W */
#define REQ_MSG    9970     /* I'm nearly done; send next chunk specs, W->F */
#define GIV_MSG    9960     /* Transfer data for this chunk to ..., F->W */
#define HLP_MSG    9950     /* Here's some chunk data, W->W */
#define RES_MSG    9940     /* Here are the results, W->W */
#define END_MSG    8900     /* We're done! T->F; F->W  */

  /* set to 1 in order to trace messages */

#define TRACEOUT stdout
#define SEND_END_TRACE 0
#define RECV_END_TRACE 0
#define SEND_GIV_TRACE 0
#define RECV_GIV_TRACE 0
#define SEND_HLP_TRACE 0
#define RECV_HLP_TRACE 0
#define SEND_REQ_TRACE 0
#define RECV_REQ_TRACE 0
#define SEND_RES_TRACE 0
#define RECV_RES_TRACE 0
#define SEND_WKP_TRACE 0
#define RECV_WKP_TRACE 0
#define SEND_WRK_TRACE 0
#define RECV_WRK_TRACE 0

#define STATIC  0
#define FSC     1
#define GSS     2
#define FAC     3
#define AF      4
#define AWF_B   5
#define AWF_C   6
#define AWF_D   7
#define AWF_E   8
#define EXPT    9

  
  int numMethods = 10;
  char *methodNameShort[] = {
    "STATIC\0", "FSC\0", "GSS\0", "FAC\0", "AWF-B\0",
    "AWF-C\0", "AF\0", "AWF-D\0", "AWF-E\0", "EXPT\0"
  };

  char *methodNameLong[] = {
    "STATIC SCHEDULING\0",
    "FIXED SIZE SCHEDULING\0",  "GUIDED SELF SCHEDULING\0",
    "FACTORING\0",              "BATCH AWF\0",
    "CHUNK AWF\0",              "ADAPTIVE FACTORING\0",
    "BATCH AWF (chunk times)\0","CHUNK AWF (chunk times)\0",
    "EXPERIMENTAL FSC\0"
  };
  
  int gP;               /* size of group  */
  int gN;               /* work for group  */
  int myRank;           /* MPI rank of this process  */

  int probeFreq;        /* iterates to do before message probe  */
  int sendRequest;      /* iterates left before sending request  */

  int method;

  int itersScheduled;   /* total iterates scheduled */
  int batchSize;        /* iterations in batch */
  int batchRem;         /* remaining in batch */
  int minChunkSize;     /* minimum chunk size */
  int maxChunkSize;     /* maximum chunk size */
  int chunkFSC;         /* assume # of chunks is same as FAC */

  unsigned char *buf;  /*Send and receive data*/
  int position;        
  int size;            /*Size of buffer*/
  int local1=0;        /*Flag for data received from other processor in local_facts1*/
  int local2=0;        /*Flag for data received from other processor in local_facts2*/
  int running1=0;      /*Flag for local_compute1 started*/
  int running2=0;      /*Flag for local_compute2 started*/

void GetChunkSize ( int method, int source, int *yMap, int *chunkSize,
		      double *stats) {
    int i, tChunk, rem;
    double awap, trw, weight;
    double bigD, bigT;
    double tMu, tSigma;

    rem = gN-itersScheduled;
    switch ( method ) {

    case FSC : // fixed size scheduling 
      tChunk = min(minChunkSize, rem);
      batchSize = tChunk; 
      batchRem = min( batchSize, rem);
      break;
    
    case GSS : // guided self scheduling 
      tChunk = max ( minChunkSize, (rem+gP-1)/gP);
      tChunk = min ( rem, tChunk );
      batchSize = tChunk; 
      batchRem = min (batchSize, rem);
      break;
      
    case FAC : // factoring 
      if (batchRem == 0) { 
        tChunk = max ( minChunkSize, (rem+2*gP-1)/(2*gP) );
        batchSize = gP*tChunk;
        batchRem = min (batchSize, rem);
      }
      // else use current batchSize 
      tChunk = max(minChunkSize, batchSize/gP);
      tChunk = min ( rem, tChunk );
      break;

    case AWF_B : case AWF_D : // batch adaptive weighted factoring 
      if (stats[4*source] <= 0.0) {
        tChunk = minChunkSize;
        batchSize = min (tChunk, rem);
        batchRem = batchSize;
      }

      else { 
        tMu = 0.0;  // default mu 
        for (i=0; i<gP; i++)
          if ( stats[4*i+1] > tMu ) tMu = stats[4*i+1];

        awap = 0.0; // average weighted performance 
        for (i=0; i<gP; i++) 
          if (stats[4*i+1] > 0.0) 
            awap = awap + stats[4*i+1];
          else  
            awap = awap + tMu;
        awap = awap/(double) gP;

        trw = 0.0; // total ref weight (refwt[i] = awap/wap[i] 
        for (i=0; i<gP; i++) 
          if (stats[4*i+1] > 0.0) 
            trw = trw + awap/stats[4*i+1];
          else  
            trw = trw + awap/tMu;

        // normalized weight for source 
        weight = ((awap/stats[4*source+1])*(double) gP)/trw;

        if (batchRem == 0) { 
          tChunk = max( minChunkSize, (rem+2*gP-1)/(2*gP) );
          batchSize = gP*tChunk;
          batchRem = min (batchSize, rem);
        }
        // else use current batchSize 
        tChunk = (int) (weight*(batchSize/gP) + 0.55);
        tChunk = max( minChunkSize, tChunk);
        tChunk = min( rem, tChunk);
      }
      break;

    case AWF_C : case AWF_E : // batch adaptive weighted factoring 
      if (stats[4*source] <= 0.0)
        tChunk = minChunkSize;
      else { 
        tMu = 0.0;  // default mu 
        for (i=0; i<gP; i++)
          if (stats[4*i+1] > tMu) tMu = stats[4*i+1];

        awap = 0.0; // average weighted performance 
        for (i=0; i<gP; i++) 
          if (stats[4*i+1] > 0.0) 
            awap = awap + stats[4*i+1];
          else  
            awap = awap + tMu;
        awap = awap/(double) gP;

        trw = 0.0; // total ref weight (refwt[i] = awap/wap[i] 
        for (i=0; i<gP; i++) 
          if (stats[4*i+1] > 0.0) 
            trw = trw + awap/stats[4*i+1];
          else  
            trw = trw + awap/tMu;

        // normalized weight for source 
        weight = ((awap/stats[4*source+1])*(double) gP)/trw;
        tChunk = (int) (weight*((rem+2*gP-1)/(2*gP)) + 0.55);
      }
      tChunk = max( minChunkSize, tChunk);
      batchSize = tChunk;
      batchRem =  min (tChunk, rem);
      break;

    case AF : // adaptive factoring 
      if (stats[4*source] <= 0.0) 
        tChunk = minChunkSize;
      else {
        tMu = 0.0;  // default mu 
        for (i=0; i<gP; i++)
          if (stats[4*i+1] > tMu) {
            tMu = stats[4*i+1];
            tSigma = stats[4*i+2];
          }
        bigD = 0.0;
        bigT = 0.0;
        for (i=0; i<gP; i++) 
          if (stats[4*i+1] > 0.0) {
            bigD = bigD + stats[4*i+2]/stats[4*i+1];
            bigT = bigT + 1.0e+00/stats[4*i+1];
          }
          else  {
            bigD = bigD + tSigma/tMu;
            bigT = bigT + 1.0e+00/tMu;
          }
        bigT = 1.0/bigT;
        // compute chunk size for worker 
        tChunk = (int) (0.5*(bigD + 2.0*bigT*rem - 
			     sqrt(bigD*(bigD + 4.0*rem*bigT)))/stats[4*source+1]);
        tChunk = min( maxChunkSize, tChunk);
      }
      tChunk = max( minChunkSize, tChunk);
      batchSize = tChunk;
      batchRem = min(batchSize,rem);

      break;

    case EXPT : // fixed size scheduling, P*log(N/P) chunks 
      tChunk = min(chunkFSC, rem);
      batchSize = tChunk; 
      batchRem = min( batchSize, rem);
      break;

    default: // fixed size scheduling, P*log(N/P) chunks 
      tChunk = min(chunkFSC, rem);
      batchSize = tChunk; 
      batchRem = min (batchSize, rem);
      break;

    }

    // adjust according to Remaining[source] 
    if (yMap[2*source+1] > (tChunk+minChunkSize/2) ) 
      *chunkSize = tChunk;
    else
      *chunkSize = yMap[2*source+1];

  // adjust remaining in batch 
    batchRem -= *chunkSize;
    if (batchRem < minChunkSize) batchRem = 0;

}


void  SetBreaks ( int breakAfter, int requestWhen, int wSize ) {
    // execute subchunks of size 10% of original chunk 
    if (breakAfter<0) probeFreq = max( 1.0, wSize*0.075);
    // send request before executing last subchunk 
    if (requestWhen<0) sendRequest = probeFreq;
}


void  procPerformance (int method, double *perfInfo, double *stats ) {
    // perfInfo(:) 0=src, 1=chunksize, 2=sumt1, 3=sumt2(AF) 
    // stats()     0=chunks done, 1=mu, 2=sigma(AF)/chunktimes(AWF), 3=chunksize(AWF) 

    int pos;

    pos = 4*(int)perfInfo[0];
    stats[pos] = stats[pos] + 1.0;                 // number of chunks done 
    if (method == AF) { //adaptive factoring 
      stats[pos+1] = perfInfo[2]/perfInfo[1];    // mu 
      if (perfInfo[1] > 1.0) {
	stats[pos+2] = (perfInfo[3] - perfInfo[1]*stats[pos+1]*stats[pos+1]) / 
	  (perfInfo[1]-1.0);        // sigma 
	if (stats[pos+2] < 0.0) stats[pos+2] = 0.0;
      }
      else stats[pos+2] = 0.0;
    }
    else {
      // all other methods use weighted performance 
      stats[pos+2] = stats[pos+2] + stats[pos]*perfInfo[2]; // weighted chunk times 
      stats[pos+3] = stats[pos+3] + stats[pos]*perfInfo[1]; // weighted chunk sizes 
      stats[pos+1] = stats[pos+2]/stats[pos+3];
    }
}

  /*To allocate inputs and outputs over the iterate space*/
void Allocate_func(fact_db &local_facts,entitySet &exec_set,variableSet inputs,variableSet outputs){

    
    for(variableSet::const_iterator vi=inputs.begin();vi!=inputs.end();++vi) {
      storeRepP sp = local_facts.get_variable(*vi) ;
      sp->allocate(interval(0,exec_set.size()-1)) ;
    }
    for(variableSet::const_iterator vi=outputs.begin();vi!=outputs.end();++vi) {
      storeRepP sp = local_facts.get_variable(*vi) ;
      sp->allocate(interval(0,exec_set.size()-1)) ;
    } 
    
}
  /*To deallocate inputs and outputs*/
void Deallocate_func(fact_db &local_facts,variableSet inputs,variableSet outputs){

    //deallocate the temporaries
    for(variableSet::const_iterator vi=inputs.begin();vi!=inputs.end();++vi) {
      storeRepP sp = local_facts.get_variable(*vi) ;
      sp->allocate(EMPTY) ;
    } 
    for(variableSet::const_iterator vi=outputs.begin();vi!=outputs.end();++vi) {
      storeRepP sp = local_facts.get_variable(*vi) ;
      sp->allocate(EMPTY) ;
    }
}

  // Transfer inputs to dest 
void SendInput (int tStart, int tSize, int dest, int tag,MPI_Comm procGrp,entitySet &exec_set,variableSet &inputs,fact_db &facts) {  
 
    // Pack inputs from facts
    position = 0 ;
    size = 0 ;
    for(variableSet::const_iterator vi=inputs.begin();vi!=inputs.end();++vi) {
      storeRepP s_ptr = facts.get_variable(*vi) ;
      size += s_ptr->pack_size(interval(tStart,tStart+tSize-1)) ;
    }
   
    buf = new unsigned char[size] ;
    entitySet myent=interval(tStart,tStart+tSize-1);
    for(variableSet::const_iterator vi=inputs.begin();vi!=inputs.end();++vi) {
      storeRepP s_ptr = facts.get_variable(*vi) ;
      s_ptr->pack(buf,position,size,myent) ;
    }
    //Send inputs 
    int send_array[2];
    send_array[0] = size;
    send_array[1] = tSize;
    MPI_Send(send_array, 2, MPI_INT, dest, tag, procGrp);
    MPI_Send(buf,size,MPI_UNSIGNED_CHAR,dest,tag,procGrp);  
    delete [] buf;
}

  // Receive inputs from sender 
void ReceiveInput (int rcvStart,int &rcvSize,MPI_Status *Status,MPI_Comm procGrp,entitySet &exec_set,variableSet &inputs,fact_db &local_facts) { 

    MPI_Status tStatus;
    int recvArray[2];
    int size;
    MPI_Recv(recvArray, 2, MPI_INT, (*Status).MPI_SOURCE, (*Status).MPI_TAG, procGrp, &tStatus);
    size = recvArray[0];
    rcvSize = recvArray[1] ;
    buf = new unsigned char[size];
    MPI_Recv (buf,size,MPI_UNSIGNED_CHAR, 
	      (*Status).MPI_SOURCE, (*Status).MPI_TAG, procGrp, &tStatus);

   
  
    // unpack inputs into local facts
    position = 0 ;
    for(variableSet::const_iterator vi=inputs.begin();vi!=inputs.end();++vi) {
      storeRepP s_ptr = local_facts.get_variable(*vi) ;
      s_ptr->unpack(buf,position,size,sequence(interval(rcvStart,rcvStart+rcvSize-1))) ;
    }

    // Delete allocated temporaries
    delete[] buf ;
 
}
  //Send outputs to dest
void SendOutput (int tStart, int tSize, int dest,int tag,MPI_Comm procGrp,entitySet &exec_set,variableSet &outputs,fact_db &local_facts) { 
  
    // Pack outputs
    position = 0 ;
    size = 0 ;
    for(variableSet::const_iterator vi=outputs.begin();vi!=outputs.end();++vi) {
      storeRepP s_ptr = local_facts.get_variable(*vi) ;
      size += s_ptr->pack_size(interval(0,tSize-1)) ;
    }
  
    buf = new unsigned char[size] ;
    entitySet myent2=interval(0,tSize-1);
    for(variableSet::const_iterator vi=outputs.begin();vi!=outputs.end();++vi) {
      storeRepP s_ptr = local_facts.get_variable(*vi) ;
      s_ptr->pack(buf,position,size,myent2) ;
    }
    //Send outputs 
    MPI_Send(&size, 1, MPI_INT, dest, tag, procGrp);
    MPI_Send(buf,size,MPI_UNSIGNED_CHAR, dest, tag, procGrp);
    delete [] buf;

  
}
  //Receive outputs from src
void ReceiveOutput (int rcvStart, int rcvSize, int src, int tag,MPI_Comm procGrp,entitySet &exec_set,variableSet &inputs,variableSet &outputs,fact_db &facts) { 

    MPI_Status tStatus;
    int size;
    MPI_Recv(&size, 1, MPI_INT, src, tag, procGrp, &tStatus);
    buf = new unsigned char[size];
    MPI_Recv (buf, size, MPI_UNSIGNED_CHAR, src, tag, procGrp, &tStatus);
    
    // unpack outputs into facts
    position = 0 ;
    for(variableSet::const_iterator vi=outputs.begin();vi!=outputs.end();++vi) {
      storeRepP s_ptr = facts.get_variable(*vi) ;
      s_ptr->unpack(buf,position,size,sequence(interval(rcvStart,rcvStart+rcvSize-1))) ; 
    }
 
    // Delete allocated temporaries
    delete[] buf ;
}




void ExecuteLoop (rule_implP rp1,entitySet &e,int method,int *yMap,fact_db &facts,fact_db &local_facts1,fact_db &local_facts2,rule_implP local_compute1,rule_implP local_compute2,variableSet &inputs,variableSet &outputs) { 
   
 
    int foreMan;      // MPI rank of foreMan 
    int minChunk=-1;     // minimum chunk size 
    int breakAfter=-1;   // iterates to execute before probing for messages 
    int requestWhen=-1;  //iterates remaining in chunk before requesting next chunk 
    int chunkMap[3*MAX_CHUNKS];
    double stats[4*MAX_PROCESSES];
    MPI_Comm procGrp=MPI_COMM_WORLD;
 

/*
   yMap(4*r+0) - block start            
   yMap(4*r+1) - block size             
   chunkMap(3*c+0) - chunk start        
   chunkMap(3*c+1) - chunk size         
   chunkMap(3*c+2) - worker             
   stats(4*r+0) - chunks done            
   stats(4*r+1) - mu(AF)/wap(AWF)        
   stats(4*r+2) - sigma(AF)/wtd chunk times(AWF)
   stats(4*r+3) - wtd chunk sizes(AWF)
   perfInfo() 0=src, 1=chunksize, 2=sumt1, 3=sumt2(AF) 
 */

  // message buffers, MPI_Recv status 
    int chunkInfo[3];             // Chunk info buffer 
    double perfInfo[4];           // PerformancE info buffer 
    MPI_Status mStatus, tStatus;

    // variables used by foreMan 
    int worker;                   // Rank of worker of a chunk 
    int chunkSize;                // size of chunk for a worker 
    int loc;                      // source of chunk to be migrated 
    int numENDed;
  
  // variables used everybody 
    int myChunks;                 // no. of chunks in this process 
    int rStart, rSize, rSource;   // start, size, source of current chunk 
    int wStart, wSize;            // start, size of remaining subchunk 
    int nextStart, nextSize;      // foreMan's response to advance request 
    int nextSource;               // owner of chunk in advance request 
    //  int GIVbuffer1, GIVbuffer2;   // where to store GIV chunks 
    int rcvStart;
    int rcvEnd;
    int GIVpending;               // results of GIV chunks to wait for 
    int i, j, tStart, tSize;

    int gotWork;                  // termination flag 
    int MsgInQueue;               // message came in 
    int req4WRKsent, nextWRKrcvd; // flags for advance request, response 

    double t0, t1, t2, t3, tk, sumt1, sumt2;
    double workTime;              // time spent doing useful work 
    double maxCost, tCost;

    // Initializations 
    gP = Loci::MPI_processes ;
    myRank = Loci::MPI_rank ;   
    //  foreMan=gP-1;
    foreMan=gP-1;
 
  
    chunkMap[0] = yMap[2*myRank]; // start of data 
    chunkMap[1] = yMap[2*myRank+1]; // size of data 
    chunkMap[2] = 0; // chunks in this rank 

    if (method == STATIC) { // no load balancing 
      if (chunkMap[1] > 0) {  
	workTime = MPI_Wtime();
	rp1->compute(sequence(interval(chunkMap[0],chunkMap[0]+chunkMap[1]-1))) ;
	stats[0] = MPI_Wtime() - workTime;
	chunkMap[2] = 1;
	chunkMap[3] = chunkMap[0];
	chunkMap[4] = chunkMap[1];
	chunkMap[5] = myRank;
      } else {
	stats[0] = 0.0;
	chunkMap[2] = 0;
      }
   
      return;
    }

    // total work 
    gN = 0; 
    for (loc=0; loc<gP; loc++) 
      gN = gN +  yMap[2*loc+1]; 
    maxChunkSize = (gN+2*gP-1)/(2*gP); // This is based on their buffers for receiving data...maximum limit is this

   

    if (yMap[2*myRank+1] == 0) {
      // tell foreMan that I'm ready if I don't have data 
      MPI_Send (NULL, 0, MPI_INT, foreMan, WKP_MSG, procGrp);
#if SEND_WKP_TRACE
      fprintf(TRACEOUT, "WKP_SEND  %d->%d\n", myRank, foreMan);
#endif
    }
   
    if (myRank == foreMan) { // scheduler initializations 
    
      itersScheduled = 0; 
      batchSize = 0; 
      batchRem = 0; 
      tSize = (gN+gP-1)/gP;
      chunkFSC = (0.55+tSize*log(2.0e0)/log( (double) tSize ) );
      if (minChunk>0) minChunkSize = minChunk;
      else minChunkSize = max(1,chunkFSC/2);       // default min chunk size ...in case user does not indicate
     
      // initialize performance measures
      for (loc=0; loc<gP; loc++) {
	stats[4*loc  ] = 0.0;  // chunks done
	stats[4*loc+1] = 0.0;  // mu /wap
	stats[4*loc+2] = 0.0;  // sigma/wtd chunk times 
	stats[4*loc+3] = 0.0;  // wtd chunk sizes
      }
   
      // send work to procs with data
      for (worker=0; worker<gP; worker++) {
	if (yMap[2*worker+1] > 0) { // worker has some 

	  GetChunkSize (method, worker, yMap, &chunkSize, stats);

	  chunkInfo[0] = yMap[2*worker]; 
	  chunkInfo[1] = chunkSize;
	  chunkInfo[2] = -1;
#if SEND_WRK_TRACE
	  fprintf(TRACEOUT, 
		  "WRK_MSG   to %d: start=%6d size=%6d, remIters=%6d; batch=%6d, batchrem=%6d\n",
		  worker, chunkInfo[0], chunkInfo[1], 
		  yMap[2*worker+1]-chunkSize,
		  batchSize, batchRem);
#endif
	  MPI_Send (chunkInfo, 3, MPI_INT, worker, WRK_MSG,  procGrp);

	  // update process status 
	  yMap[2*worker] += chunkSize;
	  yMap[2*worker+1] -= chunkSize;
	  itersScheduled += chunkSize;
	
	}
      }
    } // if (myRank == foreMan) 

    
    GIVpending = 0;  // GIV chunk results to receive 
    numENDed = 0;
    probeFreq = max(1, breakAfter); 
    sendRequest = max(1, requestWhen); 
    wSize = 0;  // remaining iterates in current chunk 
    gotWork = 1;
    workTime = 0.0; 
    myChunks = 0; 
    nextWRKrcvd = 0;
    req4WRKsent = 0;
    
    // check for any messages 
    MPI_Iprobe (MPI_ANY_SOURCE, MPI_ANY_TAG, procGrp, &MsgInQueue, &mStatus);
    while (gotWork+wSize+GIVpending+MsgInQueue) {
     
      
      // COMMUNICATIONS 
      while (MsgInQueue) {
     
	switch ( mStatus.MPI_TAG ) {
  
        case WRK_MSG :
          MPI_Recv (chunkInfo, 3, MPI_INT, mStatus.MPI_SOURCE, mStatus.MPI_TAG, 
		    procGrp, &tStatus);
          myChunks++;
          chunkMap[3*myChunks  ] = chunkInfo[0];
          chunkMap[3*myChunks+1] = chunkInfo[1];
          chunkMap[3*myChunks+2] = myRank;
          if (wSize == 0) { // no pending chunk 
            t0 = MPI_Wtime(); // elapsed time for chunk starts here 
            wStart = chunkInfo[0]; 
            wSize = chunkInfo[1]; 
            rStart = wStart; 
            rSize = wSize; 
            rSource = myRank; 
            req4WRKsent = 0; // cancel request for work 
	   
            SetBreaks ( breakAfter, requestWhen, wSize );
#if RECV_WRK_TRACE
            fprintf(TRACEOUT, "WRK_RECV0 %d<-%d: start=%6d size=%6d probe=%6d\n",
		    myRank, mStatus.MPI_SOURCE, wStart, wSize, probeFreq);
#endif
            sumt1 = 0.0;   //for mu/wap 
            sumt2 = 0.0;  // for sigma 
          }
          else {  // current chunk is not finished  save as next chunk 
            nextStart = chunkInfo[0]; 
            nextSize = chunkInfo[1]; 
            nextSource = myRank; 
            nextWRKrcvd = 1;
#if RECV_WRK_TRACE
            fprintf(TRACEOUT, "WRK_RECV1 %d<-%d: nextStart=%6d nextSize=%6d\n",
		    myRank, mStatus.MPI_SOURCE, nextStart, nextSize);
#endif
          }
	  break;

        case GIV_MSG :
          MPI_Recv (chunkInfo, 3, MPI_INT, mStatus.MPI_SOURCE, mStatus.MPI_TAG, 
		    procGrp, &tStatus);
#if RECV_GIV_TRACE
          fprintf(TRACEOUT, "GIV_RECV  %d<-%d: start=%6d, size=%6d\n",
		  myRank, chunkInfo[2], chunkInfo[0], chunkInfo[1]);
#endif
          
	  SendInput(chunkInfo[0], chunkInfo[1], chunkInfo[2], HLP_MSG, procGrp,e,inputs,facts);
          GIVpending++;
	  
#if SEND_HLP_TRACE
          fprintf(TRACEOUT, "HLP_SEND  %d->%d: helpStart=%6d, helpSize=%6d, pending=%6d\n",
		  myRank, chunkInfo[2], chunkInfo[0], chunkInfo[1], GIVpending);
#endif
	  
          myChunks++;
          chunkMap[3*myChunks  ] = chunkInfo[0]; 
          chunkMap[3*myChunks+1] = -chunkInfo[1];  // flag as GIV chunk 
          chunkMap[3*myChunks+2] = chunkInfo[2]; 
	  break;

        case HLP_MSG :
	 
          if (wSize == 0) { // no pending chunk 
            t0 = MPI_Wtime(); // elapsed time for chunk starts here 
            Allocate_func(local_facts1,e,inputs,outputs); 
	    ReceiveInput(0,tSize, &mStatus, procGrp,e,inputs,local_facts1);
#if RECV_HLP_TRACE
            fprintf(TRACEOUT, "HLP_RECV0 %d<-%d: helpStart=%6d, helpSize=%6d\n",
		    myRank, mStatus.MPI_SOURCE, 0,tSize);
#endif
            wStart =0;
            wSize =tSize; 
            rStart =0;
            rSize = tSize;
            rSource = mStatus.MPI_SOURCE;  
            req4WRKsent = 0; // cancel request for work 

            SetBreaks ( breakAfter, requestWhen, wSize );

            sumt1 = 0.0; // for mu/wap 
            sumt2 = 0.0; // for sigma 
          
            local1=1; //set flag for local_facts1
	    running1 = 0;
          } 
          else { // current chunk is not finished; save as next chunk 
	   
	    if(local1==1 && local2 == 0){
	    Allocate_func(local_facts2,e,inputs,outputs);
            ReceiveInput(0,tSize, &mStatus, procGrp,e,inputs,local_facts2); 
#if RECV_HLP_TRACE
            fprintf(TRACEOUT, "HLP_RECV1 %d<-%d: nextHelpStart=%6d, nextHelpSize=%6d\n",
		    myRank, mStatus.MPI_SOURCE, 0,tSize);
#endif
            nextStart = 0; 
            nextSize = tSize; 
            nextSource = mStatus.MPI_SOURCE;  
            nextWRKrcvd = 1;
	    local2=1; //set flag for local_facts2
	    running2 = 0;
            }
            else if(local2==1 && local1 == 0){
	      Allocate_func(local_facts1,e,inputs,outputs);
	      ReceiveInput(0,tSize, &mStatus, procGrp,e,inputs,local_facts1); 
#if RECV_HLP_TRACE
            fprintf(TRACEOUT, "HLP_RECV1 %d<-%d: nextHelpStart=%6d, nextHelpSize=%6d\n",
		    myRank, mStatus.MPI_SOURCE, 0,tSize);
#endif
            nextStart = 0; 
            nextSize = tSize; 
            nextSource = mStatus.MPI_SOURCE;  
            nextWRKrcvd = 1;
	    local1=1; //set flag for local_facts2
	    running1= 0;
	    }
	    else if(local2 == 1 && local1 == 1) {
	      std::cerr<<"Problem in the local protocol" << std::endl;
	      exit(-1);
	    }
	    else {
	      std::cerr<<"Problem in the local protocol" << std::endl;
	      exit(-1);
	    }
          }

          
	  break;

        case RES_MSG :
          // find matching chunk info 
          for (i=myChunks; i>0; i--)
            if ( (chunkMap[3*i+1] < 0) && 
                 (chunkMap[3*i+2] == mStatus.MPI_SOURCE) ) loc = i;
          chunkMap[3*loc+1] = -chunkMap[3*loc+1];
	 
	  ReceiveOutput(chunkMap[3*loc],chunkMap[3*loc+1], mStatus.MPI_SOURCE, mStatus.MPI_TAG, procGrp,e,inputs,outputs,facts);
          GIVpending--;
#if RECV_RES_TRACE
          fprintf(TRACEOUT, "RES_RECV  %d<-%d:  resStart=%6d,  resSize=%6d, pending=%6d\n",
		  myRank, mStatus.MPI_SOURCE, chunkMap[3*loc], chunkMap[3*loc+1], GIVpending-1); 
#endif
        
	  break;

        case REQ_MSG : // received by foreMan only 
        case WKP_MSG : 
          worker = mStatus.MPI_SOURCE;
          if (mStatus.MPI_TAG == WKP_MSG) { // process just woke up 
            MPI_Recv (NULL, 0, MPI_INT, mStatus.MPI_SOURCE, 
		      mStatus.MPI_TAG, procGrp, &tStatus); 
#if RECV_WKP_TRACE
            fprintf(TRACEOUT, "WKP_RECV %d<-%d, share=%d\n",
		    myRank, worker, yMap[2*worker+1]);
#endif
          } 
          else { // REQ_MSG 
            MPI_Recv (perfInfo, 4, MPI_DOUBLE, mStatus.MPI_SOURCE, 
		      mStatus.MPI_TAG, procGrp, &tStatus); 
#if RECV_REQ_TRACE
            fprintf(TRACEOUT, "REQ_RECV %d<-%d, owner=%d, size=%e, time=%e\n",
		    myRank, worker, (int) perfInfo[0], perfInfo[1], perfInfo[2]);
#endif
            procPerformance (method, perfInfo, stats );
          }

          if (yMap[2*worker+1] > 0) {// worker has some 
            GetChunkSize (method, worker, yMap, &chunkSize, stats);
            // send chunk info to worker 
            chunkInfo[0] = yMap[2*worker]; 
            chunkInfo[1] = chunkSize;
            chunkInfo[2] = -1;
            MPI_Send (chunkInfo, 3, MPI_INT, worker, WRK_MSG,  procGrp);
#if SEND_WRK_TRACE
            fprintf(TRACEOUT, 
		    "WRK_SEND %d->%d: start=%d size=%d, rem=%d; bsize=%d, brem=%d\n",
		    myRank, worker, chunkInfo[0], chunkInfo[1], 
		    yMap[2*worker+1]-chunkSize, batchSize, batchRem);
#endif
            // update process status 
            yMap[2*worker] += chunkSize; 
            yMap[2*worker+1] -= chunkSize; 
            itersScheduled += chunkSize; 

          }
	   else { // worker has none; find source and migrate a chunk 
           loc = -1;
           maxCost = 0.0;
           for (i=0; i<gP; i++) {
             tCost = yMap[2*i+1]; // *stats[4*i+1]
             if (tCost > maxCost) {
                loc = i;           //which processor has largest
                maxCost = tCost;   //gives size remaining
             }
           }
           if (loc > -1) { // loc has largest remaining cost 
             GetChunkSize (method, loc, yMap, &chunkSize, stats);
             // send memo to loc to migrate chunk to worker 
             chunkInfo[0] = yMap[2*loc]; 
             chunkInfo[1] = chunkSize; 
             chunkInfo[2] = worker; 
             MPI_Send (chunkInfo, 3, MPI_INT, loc, GIV_MSG, procGrp);
#if SEND_GIV_TRACE
	      fprintf(TRACEOUT, "GIV_SEND %d->%d->%d, start=%6d size=%6d, rem=%6d\n",
		      myRank, loc, worker, chunkInfo[0], chunkInfo[1], yMap[2*loc+1]);
#endif
             // update process status 
             yMap[2*loc] += chunkSize;
             yMap[2*loc+1] -= chunkSize; 
             itersScheduled += chunkSize; 

	     } 
            else if (worker!=myRank) {
              // no chunk to migrate; terminate requesting process 
              numENDed++;
#if SEND_END_TRACE
              fprintf(TRACEOUT, "END_MSG   to  %d; %d active\n", worker, gP-numENDed);
#endif      
              MPI_Send (NULL, 0, MPI_INT, worker, END_MSG, procGrp);
            }

	   }
	    // foreman check for termination 
	    gotWork = numENDed != (gP-1);

	    break;

	  case END_MSG :
	 
	    MPI_Recv (NULL, 0, MPI_INT, mStatus.MPI_SOURCE, 
		      mStatus.MPI_TAG, procGrp, &tStatus); 
#if RECV_END_TRACE
            fprintf(TRACEOUT, "END_RECV%d<-%d; wSize=%d\n",
		    myRank, mStatus.MPI_SOURCE, wSize);
#endif
	    gotWork = 0;
	    break;

	  } // end switch 
	  MPI_Iprobe (MPI_ANY_SOURCE, MPI_ANY_TAG, procGrp, &MsgInQueue, &mStatus);

	} // while (MsgInQueue)  
  
	// COMPUTATIONS  
	tSize = min (wSize, probeFreq);
	if (tSize > 0) {
	  tStart = wStart;
	 
	  if (method == AF) { 
	    // adaptive factoring  
	    t1 = MPI_Wtime();
	    t2 = 0.0;
	    t3 = t1;
	    for (loc=tStart; loc<(tStart+tSize); loc++) {
	      i = loc;

	     if(local1==1 && local2==0) {
	      running1 = 1;   
	      local_compute1->compute(sequence(interval(i,i))) ; 
	     }
	     else if(local2==1 && local1==0){
	      running2 = 1;
	      local_compute2->compute(sequence(interval(i,i))) ; 
	     }
	     else if(local1==1 && local2 == 1) {
	      if(running1 == 1) {
		local_compute1->compute(sequence(interval(i,i))) ; 
	      }
	      else if(running2 == 1) {
		local_compute2->compute(sequence(interval(i,i))) ; 
	      }
	      else {
		std::cerr<<"Problem in the running protocol" << std::endl;
		exit(-1);
	      }
	     }
	     else if(local1==0 && local2==0){
	         
              rp1->compute(sequence(interval(i,i))) ;  
	     }
	      // workProc (i,1,rp1);
	      tk = MPI_Wtime();
	      t2 = t2 + (tk-t3)*(tk-t3);
	      t3 = tk;
	    }
	    t1 = tk - t1;
	    sumt2 = sumt2 + t2;
	  } 
	  else {
       
	    // non-timestepping adaptive weighted factoring, etc  
	    t1 = MPI_Wtime();
	    loc=tStart;
	    i = tSize;
	    

	    if(local1==1 && local2==0) {
	      running1 = 1;   
	      local_compute1->compute(sequence(interval(loc,loc+i-1))) ; 
	    }
	    else if(local2==1 && local1==0){
	      running2 = 1;
	      local_compute2->compute(sequence(interval(loc,loc+i-1))) ; 
	    }
	    else if(local1==1 && local2 == 1) {
	      if(running1 == 1) {
		local_compute1->compute(sequence(interval(loc,loc+i-1))) ; 
	      }
	      else if(running2 == 1) {
		local_compute2->compute(sequence(interval(loc,loc+i-1))) ; 
	      }
	      else {
		std::cerr<<"Problem in the running protocol" << std::endl;
		exit(-1);
	      }
	    }
	    else if(local1==0 && local2==0){
	         
              rp1->compute(sequence(interval(loc,loc+i-1))) ;  
	    }

	    tk = MPI_Wtime();
	    t1 = tk - t1;
	  }
	  wStart += tSize;
	  wSize -= tSize;
	  sumt1 += t1;
	  workTime += t1;

	  if (wSize <= sendRequest) { // time to send request ? 
	    if (req4WRKsent == 0) { // request not yet sent ? 

	      // send REQ_MSG  
	      perfInfo[0] = rSource;
	      perfInfo[1] = rSize-wSize;
	      if ( (method==AWF_D) || (method==AWF_E) )
		perfInfo[2] = tk-t0;    // chunk elapsed time 
	      else perfInfo[2] = sumt1;    // chunk worktime 
	      perfInfo[3] = sumt2;
	      MPI_Send (perfInfo, 4, MPI_DOUBLE, foreMan, REQ_MSG, procGrp);
	      req4WRKsent = 1; // toggle: request sent 
	      nextWRKrcvd = 0; // toggle: work not received 
#if SEND_REQ_TRACE
	      fprintf(TRACEOUT, "REQ_SEND %d->%d, perf[0]=%e, perf[1]=%e, size done=%d\n",
		      myRank, foreMan, perfInfo[0], perfInfo[1], rSize);
#endif
	    }
	  } // if (...sendRequest...)  

	  //Send from local_facts1 or local_facts2
	  if (wSize == 0) { // chunk finished 

	    if (rSource != myRank) { // return results ? 
#if SEND_RES_TRACE
	      fprintf(TRACEOUT, "RES_SEND  %d->%d:  resStart=%6d,  resSize=%6d\n",
		      myRank, rSource, rStart, rSize);
#endif
              if(local1==1 && local2==0){
	      SendOutput(0, rSize, rSource, RES_MSG, procGrp,e,outputs,local_facts1);
	      local1=0;
	      running1=0;    
	      Deallocate_func(local_facts1,inputs,outputs); 
	      }
              if(local2==1 && local1==0){
	      SendOutput(0, rSize, rSource, RES_MSG, procGrp,e,outputs,local_facts2);
	      local2=0;
	      running2=0;  
	      Deallocate_func(local_facts2,inputs,outputs); 
	      }
	      if(local1==1 && local2==1){
		if(running1==1){
		  SendOutput(0, rSize, rSource, RES_MSG, procGrp,e,outputs,local_facts1);
	          local1=0;
		  running1=0;    
		  Deallocate_func(local_facts1,inputs,outputs); 
		}
                else if(running2==1){
		SendOutput(0, rSize, rSource, RES_MSG, procGrp,e,outputs,local_facts2);
		local2=0;
		running2=0;  
		Deallocate_func(local_facts2,inputs,outputs);   
		}
		else {
		  std::cerr<<"Problem in the running protocol" << std::endl;
		  exit(-1);
		}
	      }
	    }
	  
	    

	      // reset accumulators
	      sumt1 = 0.0; // for mu 
	      sumt2 = 0.0; // for sigma 
	      rSize = 0;

	      if (nextWRKrcvd) { // foreMan already responded to advance request ? 
		t0 = MPI_Wtime(); // elapsed time for chunk starts here 
		wStart = nextStart;
		wSize = nextSize;
		rStart = wStart;
		rSize = wSize;
		rSource = nextSource;
		SetBreaks ( breakAfter, requestWhen, wSize );
		nextSize = 0;
		nextWRKrcvd = 0;
		req4WRKsent = 0;
               
	      } 
	    } // if (wSize == 0) 
	  } // if (tSize > 0) 
      
	  MPI_Iprobe (MPI_ANY_SOURCE, MPI_ANY_TAG, procGrp, &MsgInQueue, &mStatus);
	} // while (gotWork+...) 
	chunkMap[2] = myChunks; // chunks in this rank 
	stats[0] = workTime; // useful work time 
      
      }//end of ExecuteLoop

    }//end of namespace Loci
