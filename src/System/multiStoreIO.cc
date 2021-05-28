#include "multiStoreIO.h"

namespace Loci {
  /*
   * Create the appropriate File access property list
   */
  namespace hdf5_const {
    extern const int PPN(1);
    extern const int facc_type(FACC_MPIO);		/*Test file access type */
    extern const int dxfer_coll_type(DXFER_COLLECTIVE_IO); /*use collective IO*/
  }


  hid_t
  create_faccess_plist(MPI_Comm comm, MPI_Info info, int l_facc_type)
  {
#ifdef H5_HAVE_PARALLEL
    hid_t ret_pl = -1;
    herr_t ret;                 /* generic return value */
    int mpi_rank;		/* mpi variables */
                                                           
    /* need the rank for error checking macros */
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  
    ret_pl = H5Pcreate (H5P_FILE_ACCESS);
    WARN(ret_pl<0) ;

    if (l_facc_type == FACC_DEFAULT)
      return (ret_pl);
  
    if (l_facc_type == FACC_MPIO){
      /* set Parallel access with communicator */
      ret = H5Pset_fapl_mpio(ret_pl, comm, info);
      WARN(ret<0) ;

      //unleased collective metadata read/write options
      //one process read small data and broadcast it to avoid I/O access
    
      /*ret = H5Pset_all_coll_metadata_ops(ret_pl, TRUE);
	WARN(ret < 0);
	ret = H5Pset_coll_metadata_write(ret_pl, TRUE);
	WARN(ret<0);
      */
      return(ret_pl);
    }
    return(ret_pl);
#else
    return H5P_DEFAULT ;
#endif
  }

  hid_t
  create_xfer_plist(int l_xfer_type){
#ifdef H5_HAVE_PARALLEL
    hid_t xfer_plist = H5Pcreate (H5P_DATASET_XFER);
    if(l_xfer_type == DXFER_COLLECTIVE_IO)
      H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);
    else
      H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_INDEPENDENT);
    return xfer_plist;
#else
    return H5P_DEFAULT ;
#endif
  }

  using std::vector ;
  using std::ostringstream ;
  extern string PFS_Script ;

  hid_t hdf5PCreateFile(const char *name, unsigned flags, hid_t create_id, hid_t access_id, size_t file_size_estimate,MPI_Comm comm) {
#ifndef H5_HAVE_PARALLEL
    int rank = 0 ;
    MPI_Comm_rank(comm,&rank) ;
    if(rank==0)
      return H5Fcreate(name,flags,create_id,access_id) ;
    else
      return 0 ;
#else
    if(MPI_rank == 0 && PFS_Script != "") {
      ostringstream oss ;
      oss << PFS_Script << " " << name << " " << file_size_estimate ;
      string script = oss.str() ;
      int ret =system(script.c_str()) ;
      if(ret !=0)
	cerr << "Error, script '" << script << "' failed!"
	     << endl ;
    }
    return H5Fcreate(name,flags,create_id,access_id) ;
#endif
  }


  void findFileDist(int &fmin, int &fstart, int &fsz, int &delta,
		    entitySet read_set, fact_db &facts) {
    // Now create data structure distributed over read set to
    // contain mapping information from file to local
    fact_db::distribute_infoP df = facts.get_distribute_info() ;
    read_set = read_set & df->my_entities ;
    const_Map l2g ;
    const_dMap g2f ;
    l2g = df->l2g.Rep() ;
    g2f = df->g2f.Rep() ;

    // Find bounds of file number
    int fmin_local = g2f[l2g[read_set.Min()]] ;
    int fmax_local = fmin_local ;
    FORALL(read_set,ii) {
      int fid = g2f[l2g[ii]] ;
      fmin_local = min(fmin_local,fid) ;
      fmax_local = max(fmax_local,fid) ;
    } ENDFORALL ;
    int fmin_global = fmin_local ;
    MPI_Allreduce(&fmin_local,&fmin_global,1,MPI_INT,
		  MPI_MIN,MPI_COMM_WORLD) ;
    fmin = fmin_global ;
    int fmax_global = fmax_local ;
    MPI_Allreduce(&fmax_local,&fmax_global,1,MPI_INT,
		  MPI_MAX,MPI_COMM_WORLD) ;

    // Find distribution of file numbers to processors
    int nfilenums = fmax_global-fmin_global+1 ;
    const int p = MPI_processes ;
    delta = (nfilenums+p-1)/MPI_processes ;
    const int r = MPI_rank ;
    fstart = delta*r ; // start of file number on this processor
    fsz = delta ; // size of file numbers on this processor
    if(delta*(r+1) > nfilenums)
      fsz = max(nfilenums-delta*r,0) ;
  }    

  void gatherFileMapping(vector<int> &map_local,
			 vector<int> &map_proc,
			 int fmin, int fstart, int fsz, int delta,
			 entitySet read_set,
			 fact_db &facts) {
    // Now create data structure distributed over read set to
    // contain mapping information from file to local
    fact_db::distribute_infoP df = facts.get_distribute_info() ;
    read_set = read_set & df->my_entities ;
    const_Map l2g ;
    const_dMap g2f ;
    l2g = df->l2g.Rep() ;
    g2f = df->g2f.Rep() ;

    const int p = MPI_processes ;
    // Collect information about mapping between file number and
    // global-processor numbering
    // Find file number to processor mapping
    vector<int> sendto(p,0) ;
    FORALL(read_set,ii) {
      int fid = g2f[l2g[ii]]-fmin ;
      int ps = fid/delta ;
      sendto[ps]++ ;
    } ENDFORALL ;
    vector<int> recvfrom(p,0) ;
    MPI_Alltoall(&sendto[0],1,MPI_INT,&recvfrom[0],1,MPI_INT,
		 MPI_COMM_WORLD) ;
    const int r = MPI_rank ;

    vector<int> sendbuf(read_set.size()*3) ;
    vector<int> offsets(p,0) ;
    vector<int> send_offsets(p+1,0) ;
    vector<int> recv_offsets(p+1,0) ;
    for(int i=0;i<p;++i) {
      send_offsets[i+1] = send_offsets[i]+sendto[i] ;
      recv_offsets[i+1] = recv_offsets[i]+recvfrom[i] ;
    }
    FORALL(read_set,ii) {
      int fid = g2f[l2g[ii]]-fmin ;
      int p = fid/delta ;
      int offset = offsets[p]+send_offsets[p] ; ;
      offsets[p]++ ;
      sendbuf[offset*3+0] = fid ; // File Id
      sendbuf[offset*3+1] = ii ; // Local Id
      sendbuf[offset*3+2] = r ; // local processor
    } ENDFORALL ;
    int rsz = 0 ;
    for(int i=0;i<p;++i)
      rsz += recvfrom[i] ;
    vector<int> recvbuf(rsz*3) ;
    int nreq = 0 ;
    for(int i=0;i<p;++i) {
      if(sendto[i] > 0)
	nreq++ ;
      if(recvfrom[i] > 0)
	nreq++ ;
    }
    vector<MPI_Request> recv_Requests(nreq) ;
    int req = 0 ;
    for(int i=0;i<p;++i)
      if(recvfrom[i] > 0) {
	MPI_Irecv(&recvbuf[recv_offsets[i]*3],recvfrom[i]*3,MPI_INT,i,3,
		  MPI_COMM_WORLD,&recv_Requests[req]) ;
	req++ ;
      }
    for(int i=0;i<p;++i)
      if(sendto[i] > 0) {
	MPI_Isend(&sendbuf[send_offsets[i]*3],sendto[i]*3,MPI_INT,i,3,
		  MPI_COMM_WORLD,&recv_Requests[req]) ;
	req++ ;
      }
    vector<MPI_Status> statuslist(nreq) ;
    MPI_Waitall(nreq,&recv_Requests[0],&statuslist[0]) ;

    for(int i=0;i<rsz;++i) {
      int lf =recvbuf[i*3+0] - fstart ; // Find file number to local filenumber
      map_local[lf] = recvbuf[i*3+1] ;
      map_proc[lf] = recvbuf[i*3+2] ;
    }
    
  }
  
  // Find mapping from fileID space to local number/procid
  void findMapping(vector<int> &local_num,
		   vector<int> &procID,
		   entitySet read_set,
		   const vector<int> &fileID,
		   fact_db &facts) {

    // Find lowest number in fileID list since we are assuming that
    // the read_sets match
    int sz = fileID.size() ;
    int local_min = fileID[0] ;
    for(int i=0;i<sz;++i) {
      local_min = min(fileID[i],local_min) ;
    }

    // Handle single processor case
    const int p = MPI_processes ;
    if(p == 1) {
      int imin = read_set.Min() ;
      for(int i=0;i<sz;++i) {
	local_num[i] = fileID[i]-local_min+imin ;
	procID[i] = 0 ;
      }
      return ;
    }

    // Find minimum id from fileID input
    int global_min = local_min ;
    MPI_Allreduce(&local_min,&global_min,1,MPI_INT,
		  MPI_MIN,MPI_COMM_WORLD) ;

    int fmin=0, fstart=0, fsz=0,delta=0 ;
    // fmin is minimum file number in read_set
    // fstart is start of file number partition on this processor
    // fsz is the size of the partition on this processor
    // delta is the delta for mapping between file number and processor
    findFileDist(fmin, fstart, fsz, delta, read_set, facts) ;
    
    vector<int> map_local(fsz) ;
    vector<int> map_proc(fsz,-1) ;

    // Gather local number and processor to file distribution
    gatherFileMapping(map_local,map_proc,fmin, fstart, fsz, delta,
		      read_set,facts) ;

    // Now go through fileID and send requests to each processor to fill in
    // mapping information from fileID domain.
    vector<int> send_req(p,0) ;
    for(int i=0;i<sz;++i) {
      int fid = fileID[i]-global_min ;
      int rp = fid/delta ;
      send_req[rp]++ ;
    }
    vector<int> recv_req(p,0) ;
    MPI_Alltoall(&send_req[0],1,MPI_INT,&recv_req[0],1,MPI_INT,
		 MPI_COMM_WORLD) ;
    vector<int> sendbuffer(fileID.size()) ;
    vector<int> soffsets(p+1,0) ;
    vector<int> roffsets(p+1,0) ;
    int rsz = 0 ;
    for(int i=0;i<p;++i) {
      rsz += recv_req[i] ;
      soffsets[i+1] = soffsets[i]+send_req[i] ;
      roffsets[i+1] = roffsets[i]+recv_req[i] ;
    }
    vector<int> scounts(p,0) ;
    for(int i=0;i<sz;++i) {
      int fid = fileID[i]-global_min ;
      int rp = fid/delta ;
      sendbuffer[soffsets[rp]+scounts[rp]] = fid ;
      scounts[rp]++ ;
    }
    vector<int> recvbuffer(rsz) ;
    
    int nreq = 0 ;
    for(int i=0;i<p;++i) {
      if(send_req[i] > 0)
	nreq++ ;
      if(recv_req[i] > 0)
	nreq++ ;
    }
    vector<MPI_Request> recv_Requests(nreq) ;
    int req = 0 ;
    for(int i=0;i<p;++i)
      if(recv_req[i] > 0) {
	MPI_Irecv(&recvbuffer[roffsets[i]],recv_req[i],MPI_INT,i,3,
		  MPI_COMM_WORLD,&recv_Requests[req]) ;
	req++ ;
      }
    for(int i=0;i<p;++i)
      if(send_req[i] > 0) {
	MPI_Isend(&sendbuffer[soffsets[i]],send_req[i],MPI_INT,i,3,
		  MPI_COMM_WORLD,&recv_Requests[req]) ;
	req++ ;
      }
    vector<MPI_Status> statuslist(nreq) ;
    MPI_Waitall(nreq,&recv_Requests[0],&statuslist[0]) ;
    // Now we go the requests, get information to return
    vector<int> datasend(rsz*2) ;
    int map_error = false ;
    for(int i=0;i<rsz;++i) {
      int fid = recvbuffer[i] ;
      int flocal = fid-fstart ;
      if(flocal >= fsz) {
	map_error = true ;
	datasend[i*2+0] = -1 ;
	datasend[i*2+1] = -1 ;
      } else {
	datasend[i*2+0] = map_local[flocal] ;
	datasend[i*2+1] = map_proc[flocal] ;
      }
    }
    if(map_error)
      cerr << "error in mapping, file number out of bounds" << endl ;

    vector<int> datarecv(fileID.size()*2) ;
    
    req = 0 ;
    for(int i=0;i<p;++i)
      if(send_req[i] > 0) {
	MPI_Irecv(&datarecv[soffsets[i]*2],send_req[i]*2,MPI_INT,i,3,
		  MPI_COMM_WORLD,&recv_Requests[req]) ;
	req++ ;
      }
    for(int i=0;i<p;++i)
      if(recv_req[i] > 0) {
	MPI_Isend(&datasend[roffsets[i]*2],recv_req[i]*2,MPI_INT,i,3,
		  MPI_COMM_WORLD,&recv_Requests[req]) ;
	req++ ;
      }
    MPI_Waitall(nreq,&recv_Requests[0],&statuslist[0]) ;
    bool error = false ;

    vector<int> sc(p,0) ;
    for(int i=0;i<sz;++i) {
      int fid = fileID[i]-global_min ;
      int rp = fid/delta ;
      local_num[i] = datarecv[(soffsets[rp]+sc[rp])*2+0] ;
      procID[i] = datarecv[(soffsets[rp]+sc[rp])*2+1] ;
      sc[rp]++ ;
      if(procID[i] < 0) {
	error = true ;
      }
    }

    if(error)
      cerr << "unable to find processor mapping" << endl ;      
  }

  void distributeMapMultiStore(vector<int> &send_sz,
			       vector<int> &recv_sz,
			       vector<int> &recv_local_num,
			       const vector<int> &local_num,
			       const vector<int> &procID) {
    const int p = MPI_processes ;
    for(int i=0;i<p;++i)
      send_sz[i] = 0 ;
    for(size_t i=0;i<procID.size();++i)
      send_sz[procID[i]]++ ;
    MPI_Alltoall(&send_sz[0],1,MPI_INT, &recv_sz[0],1,MPI_INT,
		 MPI_COMM_WORLD) ;
    std::vector<int> send_data(procID.size()) ;
    std::vector<int> soffsets(p+1,0) ;
    for(int i=0;i<p;++i)
      soffsets[i+1] = soffsets[i]+send_sz[i] ;
    std::vector<int> scounts(p,0) ;
    for(size_t i=0;i<procID.size();++i) {
      int pc = procID[i] ;
      send_data[soffsets[pc]+scounts[pc]] = local_num[i] ;
      scounts[pc]++ ;
    }
    std::vector<int> roffsets(p+1,0) ;
    for(int i=0;i<p;++i)
      roffsets[i+1] = roffsets[i]+recv_sz[i] ;
    vector<int> recv_data(roffsets[p]) ;
    int nreq = 0 ;
    for(int i=0;i<p;++i) {
      if(send_sz[i] > 0)
	nreq++ ;
      if(recv_sz[i] > 0)
	nreq++ ;
    }
    vector<MPI_Request> recv_Requests(nreq) ;
    int req = 0 ;
    for(int i=0;i<p;++i)
      if(recv_sz[i] > 0) {
	MPI_Irecv(&recv_data[roffsets[i]],recv_sz[i],MPI_INT,i,3,
		  MPI_COMM_WORLD,&recv_Requests[req]) ;
	req++ ;
      }
    for(int i=0;i<p;++i)
      if(send_sz[i] > 0) {
	MPI_Isend(&send_data[soffsets[i]],send_sz[i],MPI_INT,i,3,
		  MPI_COMM_WORLD,&recv_Requests[req]) ;
	req++ ;
      }
    vector<MPI_Status> statuslist(nreq) ;
    MPI_Waitall(nreq,&recv_Requests[0],&statuslist[0]) ;
    recv_local_num.swap(recv_data) ;
  }    

  void sendCounts(std::vector<int>&recv_count,
		  const std::vector<int> &send_sz,
		  const std::vector<int> &recv_sz,
		  const std::vector<int> &recv_local_num,
		  const std::vector<int> &counts,
		  const std::vector<int> &procID) {
    const int p = MPI_processes ;
    std::vector<int> soffsets(p+1,0) ;
    std::vector<int> roffsets(p+1,0) ;
    for(int i=0;i<p;++i) {
      soffsets[i+1] = soffsets[i] + send_sz[i] ;
      roffsets[i+1] = roffsets[i] + recv_sz[i] ;
    }
    std::vector<int> so(p,0) ;
    std::vector<int> sbuf(counts.size()) ;
    for(size_t i=0;i<counts.size();++i) {
      int ps = procID[i] ;
      sbuf[soffsets[ps]+so[ps]] = counts[i] ;
      so[ps]++ ;
    }

    std::vector<int> rbuf(roffsets[p]) ;
    int nreq = 0 ;
    for(int i=0;i<p;++i) {
      if(send_sz[i] > 0)
	nreq++ ;
      if(recv_sz[i] > 0)
	nreq++ ;
    }
    vector<MPI_Request> recv_Requests(nreq) ;
    int req = 0 ;
    for(int i=0;i<p;++i)
      if(recv_sz[i] > 0) {
	MPI_Irecv(&rbuf[roffsets[i]],recv_sz[i],MPI_INT,i,3,
		  MPI_COMM_WORLD,&recv_Requests[req]) ;
	req++ ;
      }
    for(int i=0;i<p;++i)
      if(send_sz[i] > 0) {
	MPI_Isend(&sbuf[soffsets[i]],send_sz[i],MPI_INT,i,3,
		  MPI_COMM_WORLD,&recv_Requests[req]) ;
	req++ ;
      }
    vector<MPI_Status> statuslist(nreq) ;
    MPI_Waitall(nreq,&recv_Requests[0],&statuslist[0]) ;
    recv_count.swap(rbuf) ;
    
  }
}
