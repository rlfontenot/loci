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

  /* message tags */
#define WKP_MSG    9990     /* I'm alive! W->F */
#define WRK_MSG    9980     /* Do chunk with these specs, F->W */
#define REQ_MSG    9970     /* I'm nearly done; send next chunk specs, W->F */
#define GIV_MSG    9960     /* Transfer data for this chunk to ..., F->W */
#define HLP_MSG    9950     /* Here's some chunk data, W->W */
#define RES_MSG    9940     /* Here are the results, W->W */
#define END_MSG    8900     /* We're done! T->F; F->W  */
#define RTS_HLP    9930
#define RTR_HLP    9920
#define RTS_RES    9910
#define RTR_RES    9900
#define INI_MSG    9000
  /* set to 1 in order to trace messages */

#define TRACEOUT stderr
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
  /*set to 1 in order to do performance measurements*/

  /*List of load balancing techniques*/
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
  
  int gP=0;               /* size of group  */
  int gN=0;               /* work for group  */
  int myRank=0;           /* MPI rank of this process  */

  int probeFreq=0;        /* iterates to do before message probe  */
  int sendRequest=0;      /* iterates left before sending request  */

  int itersScheduled=0;   /* total iterates scheduled */
  int batchSize=0;        /* iterations in batch */
  int batchRem=0;         /* remaining in batch */
  int minChunkSize=0;     /* minimum chunk size */
  int maxChunkSize=0;     /* maximum chunk size */
  int chunkFSC=0;         /* assume # of chunks is same as FAC */
  int Nchunks=0;
 


void GetChunkSize ( int method, int source, int *yMap, int *chunkSize,
		      double *stats) {
    int i=0, tChunk=0, rem=0;
    double awap=0.0, trw=0.0, weight=0.0;
    double bigD=0.0, bigT=0.0;
    double tMu=0.0, tSigma=0.0;
    Nchunks++;
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
    if (breakAfter<0) probeFreq =int( max( 1.0, wSize*0.075));
    // send request before executing last subchunk 
    if (requestWhen<0) sendRequest = probeFreq;
}


void  procPerformance (int method, double *perfInfo, double *stats ) {
    // perfInfo(:) 0=src, 1=chunksize, 2=sumt1, 3=sumt2(AF) 
    // stats()     0=chunks done, 1=mu, 2=sigma(AF)/chunktimes(AWF), 3=chunksize(AWF) 

    int pos=0;

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


void ExecuteLoop (void (*workCompute) (int,int,int),void (*SendInput) (int,int,int,int,MPI_Comm),void (*ReceiveInput) (int,int *,int,int,MPI_Comm),void (*SendOutput) (int,int,int,int,MPI_Comm),void (*ReceiveOutput) (int,int,int,int,MPI_Comm),void (*Allocate_func) (),void (*Deallocate_func) (),int method,int *yMap,double *stats,int *chunkMap){ 
    
    int foreMan=0;      // MPI rank of foreMan 
    int minChunk=-1;     // minimum chunk size 
    int breakAfter=-1;   // iterates to execute before probing for messages 
    int requestWhen=-1;  //iterates remaining in chunk before requesting next chunk    
    MPI_Comm procGrp;//=MPI_COMM_WORLD;
    MPI_Comm_dup(MPI_COMM_WORLD, &procGrp);

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


    MPI_Status mStatus, tStatus;
    MPI_Request req1;

    // variables used by foreMan 
    int worker=0;                   // Rank of worker of a chunk 
    int chunkSize=0;                // size of chunk for a worker 
    int loc=0;                      // source of chunk to be migrated 
    int numENDed=0;
  
  // variables used everybody 
    int myChunks=0;                 // no. of chunks in this process 
    int rStart=0, rSize=0, rSource=0;   // start, size, source of current chunk 
    int wStart=0, wSize=0;            // start, size of remaining subchunk 
    int nextStart=0, nextSize=0;      // foreMan's response to advance request 
    int nextSource=0;               // owner of chunk in advance request 
    int myRemaining=0;
    int GIVpending=0;               // results of GIV chunks to wait for 
    int i=0,tStart=0, tSize=0, tSource=0;

    int gotWork=0;                  // termination flag 
    int MsgInQueue=0;               // message came in 
    int req4WRKsent=0, nextWRKrcvd=0; // flags for advance request, response 

    double t0=0.0, t1=0.0, t2=0.0, t3=0.0, tk=0.0, sumt1=0.0, sumt2=0.0;
    double workTime=0.0;              // time spent doing useful work 
    double maxCost=0.0, tCost=0.0;
    int HLP_pending=0;
    //    int Signal1=0;
    int Signal2=0;
    // Initializations 
    gP = Loci::MPI_processes ;
    myRank = Loci::MPI_rank ;   
    foreMan=gP-1;
  
    chunkMap[0] = yMap[2*myRank]; // start of data 
    chunkMap[1] = yMap[2*myRank+1]; // size of data 
    chunkMap[2] = 0; // chunks in this rank 
    myRemaining=yMap[2*myRank+1]; //to keep track of its own work

   
    if (method == STATIC) { // no load balancing 
      if (chunkMap[1] > 0) {  
	workTime = MPI_Wtime();
        workCompute(chunkMap[0],chunkMap[1],Signal2);
	stats[0] = MPI_Wtime() - workTime;
	chunkMap[2] = 1;
	chunkMap[3] = chunkMap[0];
	chunkMap[4] = chunkMap[1];
	chunkMap[5] = myRank;
      } 
      else {
	stats[0] = 0.0;
	chunkMap[2] = 0;
      }  
      return;
    }
      // message buffers, MPI_Recv status 
   
    int chunkInfo[4];             // send Chunk info buffer 
    for(int j=0;j<4;j++){
      chunkInfo[j]=0;
    }
    int SchunkInfo[4];             // GIV Chunk info buffer 
    for(int j=0;j<4;j++){
      SchunkInfo[j]=0;
    }
    int RchunkInfo[4];             // receive Chunk info buffer 
    for(int j=0;j<4;j++){
      RchunkInfo[j]=0;
    }
    int FchunkInfo[4];            //for foreman ->foreman WRK_MSG
     for(int j=0;j<4;j++){
      FchunkInfo[j]=0;
    }
    int GchunkInfo[4];            //for foreman ->foreman GIV_MSG Send
     for(int j=0;j<4;j++){
      GchunkInfo[j]=0;
    }
    int GRchunkInfo[4];            //for foreman ->foreman GIV_MSG Recv
     for(int j=0;j<4;j++){
      GRchunkInfo[j]=0;
    }
    double perfInfo[4];           // PerformancE info buffer 
    for(int m=0;m<4;m++){
      perfInfo[m]=0.0;
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
      chunkFSC = int((0.55+tSize*log(2.0e0)/log( (double) tSize ) ));
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
          chunkInfo[3]=Nchunks;
          MPI_Send (chunkInfo, 4, MPI_INT, worker, INI_MSG,  procGrp);
	
#if SEND_WRK_TRACE
	  fprintf(TRACEOUT, 
		  "INI_MSG %d  to %d: start=%6d size=%6d, remIters=%6d; batch=%6d, batchrem=%6d\n",
		  Nchunks,worker, chunkInfo[0], chunkInfo[1], 
		  yMap[2*worker+1]-chunkSize,
		  batchSize, batchRem);
#endif
	  
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
  
     //to receive intial wrkmsg from foreman
    MPI_Irecv (RchunkInfo, 3, MPI_INT, MPI_ANY_SOURCE,INI_MSG, procGrp,&req1);
    int flag_recv=0;
    while(flag_recv==0){
    MPI_Test(&req1,&flag_recv,&mStatus);
    }
   
    myChunks++;
    chunkMap[3*myChunks  ] = RchunkInfo[0];
    chunkMap[3*myChunks+1] = RchunkInfo[1];
    chunkMap[3*myChunks+2] = myRank;
    
    if (wSize == 0) { // no pending chunk 
            
     t0 = MPI_Wtime(); // elapsed time for chunk starts here 
     wStart = RchunkInfo[0]; 
     wSize = RchunkInfo[1]; 
     rStart = wStart; 
     rSize = wSize; 
     rSource = myRank;  
#if RECV_WRK_TRACE
            fprintf(TRACEOUT, "WRK_RECV0 %d<-%d: start=%6d size=%6d probe=%6d\n",
		    myRank, mStatus.MPI_SOURCE, wStart, wSize, probeFreq);
#endif
            sumt1 = 0.0;   //for mu/wap 
            sumt2 = 0.0;  // for sigma 
     }
     else {  
       std::cerr<<"Problem in receiving initial work!"<<std::endl;
     }



    // check for any messages   
    MPI_Iprobe (MPI_ANY_SOURCE, MPI_ANY_TAG, procGrp, &MsgInQueue, &mStatus);  
      if(MsgInQueue) {
	//   cerr << "Iprobe1" << MsgInQueue << " " << mStatus.MPI_TAG << endl ;
      if(mStatus.MPI_TAG == 1) {
	abort() ;
      }
      }
    while (gotWork+wSize+GIVpending+MsgInQueue+HLP_pending) {   
      // COMMUNICATIONS 
      while (MsgInQueue) {
	switch ( mStatus.MPI_TAG ) {
        case WRK_MSG :         
	  MPI_Recv (RchunkInfo, 4, MPI_INT, mStatus.MPI_SOURCE,WRK_MSG,procGrp, &tStatus);
          myChunks++;
          chunkMap[3*myChunks  ] = RchunkInfo[0];
          chunkMap[3*myChunks+1] = RchunkInfo[1];
          chunkMap[3*myChunks+2] = myRank;
        
          if (wSize == 0) { // no pending chunk 
            t0 = MPI_Wtime(); // elapsed time for chunk starts here 
            wStart = RchunkInfo[0]; 
            wSize = RchunkInfo[1]; 
            rStart = wStart; 
            rSize = wSize; 
            rSource = myRank; 
            req4WRKsent = 0; // cancel request for work 
	   
            SetBreaks ( breakAfter, requestWhen, wSize );
	   
#if RECV_WRK_TRACE
            fprintf(TRACEOUT, "WRK_RECV0 %d %d<-%d: start=%6d size=%6d probe=%6d\n",
		    RchunkInfo[3],myRank, mStatus.MPI_SOURCE, wStart, wSize, probeFreq);
#endif
            sumt1 = 0.0;   //for mu/wap 
            sumt2 = 0.0;  // for sigma 
          }
          else {  // current chunk is not finished  save as next chunk 
            nextStart = RchunkInfo[0]; 
            nextSize = RchunkInfo[1]; 
            nextSource = myRank; 
            nextWRKrcvd = 1;
#if RECV_WRK_TRACE
            fprintf(TRACEOUT, "WRK_RECV1 %d %d<-%d: nextStart=%6d nextSize=%6d\n",
		    RchunkInfo[3],myRank, mStatus.MPI_SOURCE, nextStart, nextSize);
#endif
          }
	  break;

        case GIV_MSG :
	
	   MPI_Recv (RchunkInfo, 4, MPI_INT, mStatus.MPI_SOURCE,GIV_MSG,procGrp, &tStatus);
	  
#if RECV_GIV_TRACE
          fprintf(TRACEOUT, "GIV_RECV %d  %d<-%d: start=%6d, size=%6d\n",
		 RchunkInfo[3], myRank, RchunkInfo[2], RchunkInfo[0], RchunkInfo[1]);
#endif    
	  SendInput(RchunkInfo[0], RchunkInfo[1], RchunkInfo[2], HLP_MSG, procGrp);
          GIVpending++;
	  
#if SEND_HLP_TRACE
          fprintf(TRACEOUT, "HLP_SEND  %d->%d: helpStart=%6d, helpSize=%6d, pending=%6d\n",
		  myRank, RchunkInfo[2], RchunkInfo[0], RchunkInfo[1], GIVpending);
#endif
	  myRemaining-=RchunkInfo[1];
          myChunks++;
          chunkMap[3*myChunks  ] = RchunkInfo[0]; 
          chunkMap[3*myChunks+1] = -RchunkInfo[1];  // flag as GIV chunk 
          chunkMap[3*myChunks+2] = RchunkInfo[2]; 
	 
	  break;

        case HLP_MSG :
          //	  Signal1=1;

          if (wSize == 0) { // no pending chunk 
	    Signal2=1;
            t0 = MPI_Wtime(); // elapsed time for chunk starts here 
            Allocate_func(); 
	    ReceiveInput(0,&tSize,mStatus.MPI_SOURCE,HLP_MSG, procGrp);
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
	    sendRequest=0;
            sumt1 = 0.0; // for mu/wap 
            sumt2 = 0.0; // for sigma 
          } 
          else { // current chunk is not finished; save as next chunk 
	    Allocate_func();
            ReceiveInput(0,&tSize,mStatus.MPI_SOURCE,HLP_MSG, procGrp); 
#if RECV_HLP_TRACE
            fprintf(TRACEOUT, "HLP_RECV1 %d<-%d: nextHelpStart=%6d, nextHelpSize=%6d\n",
		    myRank,mStatus.MPI_SOURCE, 0,tSize);
#endif 
            nextStart = 0; 
            nextSize = tSize; 
            nextSource =mStatus.MPI_SOURCE;  
            nextWRKrcvd = 1;
	   
          }  
	  if(myRank==foreMan){
	    HLP_pending--;
	  }
	  break;

        case RTS_RES :
         
          // find matching chunk info 
	  tSource = mStatus.MPI_SOURCE ;
	  loc=0;
          for (i=myChunks; i>0; i--)
            if ( (chunkMap[3*i+1] < 0) && 
                 (chunkMap[3*i+2] == tSource) ) loc = i;
	  if(loc==0) {
	    cerr << "matching chunk not found!" << endl ;
	    cerr << "tSource = " << tSource << endl ;
	    abort() ;
	  }
	  chunkMap[3*loc+1] = -chunkMap[3*loc+1];
	  ReceiveOutput(chunkMap[3*loc],chunkMap[3*loc+1],mStatus.MPI_SOURCE, RES_MSG, procGrp);
	  GIVpending--;
#if RECV_RES_TRACE
          fprintf(TRACEOUT, "RES_RECV  %d<-%d:  resStart=%6d,  resSize=%6d, pending=%6d\n",
		  myRank, mStatus.MPI_SOURCE, chunkMap[3*loc], chunkMap[3*loc+1], GIVpending); 
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
	    MPI_Recv (perfInfo, 4, MPI_DOUBLE, mStatus.MPI_SOURCE,REQ_MSG, procGrp, &tStatus); 	    
#if RECV_REQ_TRACE
            fprintf(TRACEOUT, "REQ_RECV %d<-%d, owner=%d, size=%e, time=%e\n",
		    myRank, worker, (int) perfInfo[0], perfInfo[1], perfInfo[2]);
#endif
            procPerformance (method, perfInfo, stats );
          }
         
          if (yMap[2*worker+1] > 0) {// worker has some 
            GetChunkSize (method, worker, yMap, &chunkSize, stats);
            //use separate buffer for foreman to send to itself
            if(worker==foreMan){
            // send chunk info to worker 
            FchunkInfo[0] = yMap[2*worker]; 
            FchunkInfo[1] = chunkSize;
            FchunkInfo[2] = -1;
            FchunkInfo[3] =Nchunks;
            MPI_Send (FchunkInfo, 4, MPI_INT, worker, WRK_MSG,  procGrp);
	  
#if SEND_WRK_TRACE
            fprintf(TRACEOUT, 
		    "WRK_SEND %d %d->%d: start=%d size=%d, rem=%d; bsize=%d, brem=%d\n",
		    Nchunks,myRank, worker, FchunkInfo[0], FchunkInfo[1], 
		    yMap[2*worker+1]-chunkSize, batchSize, batchRem);
#endif
	   
	    }
	    else{ //send to other processsors

	      // send chunk info to worker 
            chunkInfo[0] = yMap[2*worker]; 
            chunkInfo[1] = chunkSize;
            chunkInfo[2] = -1;
            chunkInfo[3] =Nchunks;
            MPI_Send (chunkInfo, 4, MPI_INT, worker, WRK_MSG,  procGrp);
	    
#if SEND_WRK_TRACE
            fprintf(TRACEOUT, 
		    "WRK_SEND %d %d->%d: start=%d size=%d, rem=%d; bsize=%d, brem=%d\n",
		    Nchunks,myRank, worker, chunkInfo[0], chunkInfo[1], 
		    yMap[2*worker+1]-chunkSize, batchSize, batchRem);
#endif
	    
	    }

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
             if(loc==foreMan){
             // send memo to loc to migrate chunk to worker 
             GchunkInfo[0] = yMap[2*loc]; 
             GchunkInfo[1] = chunkSize; 
             GchunkInfo[2] = worker; 
	     GchunkInfo[3]=Nchunks;
#if SEND_GIV_TRACE
	      fprintf(TRACEOUT, "GIV_SEND %d %d->%d->%d, start=%6d size=%6d, rem=%6d\n",
		      Nchunks,myRank, loc, worker, GchunkInfo[0], GchunkInfo[1], yMap[2*loc+1]);
#endif	   
	     SendInput(GchunkInfo[0], GchunkInfo[1], GchunkInfo[2], HLP_MSG, procGrp);
             GIVpending++;
	  
#if SEND_HLP_TRACE
          fprintf(TRACEOUT, "HLP_SEND  %d->%d: helpStart=%6d, helpSize=%6d, pending=%6d\n",
		  myRank, GchunkInfo[2], GchunkInfo[0], GchunkInfo[1], GIVpending);
#endif
	     myRemaining-=GchunkInfo[1];
             myChunks++;
             chunkMap[3*myChunks  ] = GchunkInfo[0]; 
             chunkMap[3*myChunks+1] = -GchunkInfo[1];  // flag as GIV chunk 
             chunkMap[3*myChunks+2] = GchunkInfo[2]; 
	    
	     }
	     else{
	       // send memo to loc to migrate chunk to worker 
             SchunkInfo[0] = yMap[2*loc]; 
             SchunkInfo[1] = chunkSize; 
             SchunkInfo[2] = worker; 
	     SchunkInfo[3]=Nchunks;
	   
             MPI_Send (SchunkInfo, 4, MPI_INT, loc, GIV_MSG, procGrp);
	   
#if SEND_GIV_TRACE
	      fprintf(TRACEOUT, "GIV_SEND %d %d->%d->%d, start=%6d size=%6d, rem=%6d\n",
		      Nchunks,myRank, loc, worker, SchunkInfo[0], SchunkInfo[1], yMap[2*loc+1]);
#endif
	     }
             // update process status 
             yMap[2*loc] += chunkSize;
             yMap[2*loc+1] -= chunkSize; 
             itersScheduled += chunkSize; 
	     if(worker==myRank){
 	       HLP_pending++;
	     }
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
	    MPI_Recv (NULL, 0, MPI_INT, mStatus.MPI_SOURCE,END_MSG, procGrp,&tStatus); 
           
#if RECV_END_TRACE
            fprintf(TRACEOUT, "END_RECV%d<-%d; wSize=%d\n",
		    myRank, mStatus.MPI_SOURCE, wSize);
#endif
	    gotWork = 0;
	   
	  break;

	  default:
	  cerr << "Protocol error, unknown tag " << mStatus.MPI_TAG
	       << " source " << mStatus.MPI_SOURCE << endl ;
	  abort() ;

	  } // end switch 
	  MPI_Iprobe (MPI_ANY_SOURCE, MPI_ANY_TAG, procGrp, &MsgInQueue, &mStatus);
	  //  if(MsgInQueue) {
	      // cerr << "Iprobe2" << MsgInQueue << " " << mStatus.MPI_TAG << endl ;
	  //   }

	} // while (MsgInQueue)  
  
	// COMPUTATIONS  
     
      //  if(myRemaining==0){
      //	Signal2=1;
      // }
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
	      workCompute(i,1,Signal2);   //i->i
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
	    workCompute(loc,i,Signal2);   //loc->loc+i-1
	    tk = MPI_Wtime();
	    t1 = tk - t1;
	  }
	  wStart += tSize;
	  wSize -= tSize;
	  sumt1 += t1;
	  workTime += t1;
	      //Send Results
	 
	  if (wSize == 0) { // chunk finished 
           	
	    if (rSource != myRank) { // return results ? 

	       SendOutput(0, rSize, rSource, RES_MSG, procGrp);
	       Deallocate_func(); 
#if SEND_RES_TRACE
	      fprintf(TRACEOUT, "RES_SEND  %d->%d:  resStart=%6d,  resSize=%6d\n",
		      myRank, rSource, rStart, rSize);
#endif            
	    }	    
	    else{	   
	     myRemaining-=rSize;	   
	     // if(myRemaining==0){
	     //  Signal2=1;
	     // }
	    }
	  }  
	 
	  if (wSize <= sendRequest) { // time to send request ? 
	    if (req4WRKsent == 0 && nextWRKrcvd == 0) { // request not yet sent ? 

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
	 
	   if (wSize == 0) { 
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
		if(rSource!=myRank){
		  sendRequest=0;
                  Signal2=1;
		}
		nextSize = 0;
		nextWRKrcvd = 0;
		req4WRKsent = 0;
               
	      } 
	    } // if (wSize == 0) 
	  } // if (tSize > 0) 

	  MPI_Iprobe (MPI_ANY_SOURCE, MPI_ANY_TAG, procGrp, &MsgInQueue, &mStatus);
	  //  if(MsgInQueue)
	  //   cerr << "Iprobe3" << MsgInQueue << " " << mStatus.MPI_TAG << endl ;
	} // while (gotWork+...) 
  
	chunkMap[2] = myChunks; // chunks in this rank 
	stats[0] = workTime; // useful work time 
      }//end of ExecuteLoop

    }//end of namespace Loci
