#ifndef MPI_H_STUBB
#define MPI_H_STUBB



#define MPI_STUBB
#include <iostream>

/******************/
/* MPI-1 bindings */
/******************/


#define MPI_BOTTOM		((MPI_Aint)0)

typedef long			MPI_Aint;
typedef unsigned int		MPI_Request;
typedef unsigned int		MPI_Group;
typedef unsigned int		MPI_Comm;
typedef unsigned int		MPI_Errhandler;
typedef unsigned int		MPI_Op;
typedef unsigned int		MPI_Datatype;

typedef struct { 
	int MPI_SOURCE;
	int MPI_TAG;
	int MPI_ERROR;
	int size;
	int reserved[2];
} MPI_Status;

enum {
	MPI_COMM_NULL		= 0,
	MPI_COMM_WORLD		= 1,
	MPI_COMM_SELF		= 2
};

enum {
	MPI_ERRHANDLER_NULL	= 0,
	MPI_ERRORS_ARE_FATAL	= 1,
	MPI_ERRORS_RETURN	= 2
};

enum {
	MPI_GROUP_NULL		= 0,
	MPI_GROUP_EMPTY		= 1
};

enum {
	MPI_REQUEST_NULL	= 0
};

enum {
	MPI_OP_NULL		= 0,
	MPI_MAX			= 1,
	MPI_MIN			= 2,
	MPI_SUM			= 3,
	MPI_PROD		= 4,
	MPI_LAND		= 5,
	MPI_BAND 		= 6,
	MPI_LOR			= 7,
	MPI_BOR			= 8,
	MPI_LXOR		= 9,
	MPI_BXOR		= 10,
	MPI_MAXLOC		= 11,
	MPI_MINLOC		= 12
};

enum {
	MPI_DATATYPE_NULL	= 0,

	MPI_CHAR		= 1,
	MPI_SHORT		= 2,
	MPI_INT			= 3,
	MPI_LONG		= 4,
	MPI_UNSIGNED_CHAR	= 5,
	MPI_UNSIGNED_SHORT	= 6,
	MPI_UNSIGNED		= 7,
	MPI_UNSIGNED_LONG	= 8,
	MPI_FLOAT		= 9,
	MPI_DOUBLE		= 10,
	MPI_LONG_DOUBLE		= 11,
	MPI_LONG_LONG		= 12,

	MPI_INTEGER		= 13,
	MPI_REAL		= 14,
	MPI_DOUBLE_PRECISION	= 15,
	MPI_COMPLEX		= 16,
	MPI_DOUBLE_COMPLEX	= 17,
	MPI_LOGICAL		= 18,
	MPI_CHARACTER		= 19,
	MPI_INTEGER1		= 20,
	MPI_INTEGER2		= 21,
	MPI_INTEGER4		= 22,
	MPI_INTEGER8		= 23,
	MPI_REAL4		= 24,
	MPI_REAL8		= 25,
	MPI_REAL16		= 26,

	MPI_BYTE		= 27,
	MPI_PACKED		= 28,
	MPI_UB			= 29,
	MPI_LB			= 30,

	MPI_FLOAT_INT		= 31,
	MPI_DOUBLE_INT		= 32,
	MPI_LONG_INT		= 33,
	MPI_2INT		= 34,
	MPI_SHORT_INT		= 35,
	MPI_LONG_DOUBLE_INT	= 36,

	MPI_2REAL		= 37,
	MPI_2DOUBLE_PRECISION	= 38,
	MPI_2INTEGER		= 39,
	_MPI_SGI_TYPE_LAST
};
const int MPI_TYPE_SIZE[] =
  {0,
   sizeof(char),sizeof(short),sizeof(int),sizeof(long),
   sizeof(unsigned char), sizeof(unsigned short), sizeof(unsigned int),
   sizeof(unsigned long),
   sizeof(float),sizeof(double),sizeof(long double), sizeof(long long),
   sizeof(int), sizeof(float),sizeof(double), 2*sizeof(float),2*sizeof(double),
   sizeof(bool), sizeof(char), 1,2,4,8,4,8,16,1,1,1,1,
   sizeof(float)+sizeof(int),sizeof(double)+sizeof(int),
   sizeof(long)+sizeof(int),2*sizeof(int)} ;

inline int MPI_GET_TYPE_SIZE(int type) {
  if(size_t(type) > sizeof(MPI_TYPE_SIZE)/sizeof(int))
    return 1 ;
  else
    return MPI_TYPE_SIZE[type] ;
}

  

#define MPI_LONG_LONG_INT	MPI_LONG_LONG

enum {
	MPI_SUCCESS			= 0,

	/* These 19 error codes are specified by the MPI-1 standard */

	MPI_ERR_BUFFER			= 1,
	MPI_ERR_COUNT			= 2,
	MPI_ERR_TYPE			= 3,
	MPI_ERR_TAG			= 4,
	MPI_ERR_COMM			= 5,
	MPI_ERR_RANK			= 6,
	MPI_ERR_REQUEST			= 7,
	MPI_ERR_ROOT			= 8,
	MPI_ERR_GROUP			= 9,
	MPI_ERR_OP			= 10,
	MPI_ERR_TOPOLOGY		= 11,
	MPI_ERR_DIMS			= 12,
	MPI_ERR_ARG			= 13,
	MPI_ERR_UNKNOWN			= 14,
	MPI_ERR_TRUNCATE		= 15,
	MPI_ERR_OTHER			= 16,
	MPI_ERR_INTERN			= 17,
	MPI_ERR_IN_STATUS		= 18,
	MPI_ERR_PENDING			= 19,


	/* Error codes 20-27 used by MPI on T3E systems. */

	/* These 34 error codes are specified by the MPI-2 standard */

	MPI_ERR_ACCESS			= 28,
	MPI_ERR_AMODE			= 29,
	MPI_ERR_ASSERT			= 30,
	MPI_ERR_BAD_FILE		= 31,
	MPI_ERR_BASE			= 32,
	MPI_ERR_CONVERSION		= 33,
	MPI_ERR_DISP			= 34,
	MPI_ERR_DUP_DATAREP		= 35,
	MPI_ERR_FILE_EXISTS		= 36,
	MPI_ERR_FILE_IN_USE		= 37,
	MPI_ERR_FILE			= 38,
	MPI_ERR_INFO_KEY		= 39,
	MPI_ERR_INFO_NOKEY		= 40,
	MPI_ERR_INFO_VALUE		= 41,
	MPI_ERR_INFO			= 42,
	MPI_ERR_IO			= 43,
	MPI_ERR_KEYVAL			= 44,
	MPI_ERR_LOCKTYPE		= 45,
	MPI_ERR_NAME			= 46,
	MPI_ERR_NO_MEM			= 47,
	MPI_ERR_NOT_SAME		= 48,
	MPI_ERR_NO_SPACE		= 49,
	MPI_ERR_NO_SUCH_FILE		= 50,
	MPI_ERR_PORT			= 51,
	MPI_ERR_QUOTA			= 52,
	MPI_ERR_READ_ONLY		= 53,
	MPI_ERR_RMA_CONFLICT		= 54,
	MPI_ERR_RMA_SYNC		= 55,
	MPI_ERR_SERVICE			= 56,
	MPI_ERR_SIZE			= 57,
	MPI_ERR_SPAWN			= 58,
	MPI_ERR_UNSUPPORTED_DATAREP	= 59,
	MPI_ERR_UNSUPPORTED_OPERATION	= 60,
	MPI_ERR_WIN			= 61,
	MPI_ERR_LASTCODE		= 100	/* last built-in error code */
};

enum {
	MPI_KEYVAL_INVALID	= 0,
	MPI_TAG_UB		= 1,
	MPI_HOST		= 2,
	MPI_IO			= 3,
	MPI_WTIME_IS_GLOBAL	= 4,
	MPI_UNIVERSE_SIZE	= 9,
	MPI_APPNUM		= 11
};

enum {
	MPI_IDENT		= 0,
	MPI_CONGRUENT		= 1,
	MPI_SIMILAR		= 2,
	MPI_UNEQUAL		= 3
};

enum {
	MPI_GRAPH		= 1,
	MPI_CART		= 2
};

enum {
	MPI_UNDEFINED		= -3,
	MPI_ANY_SOURCE		= -2,
	MPI_PROC_NULL		= -1
};

enum {
	MPI_ANY_TAG		= -1
};

enum {
	MPI_BSEND_OVERHEAD	= 32
};

enum {
	MPI_MAX_PROCESSOR_NAME	= 256
};

enum {
	MPI_MAX_ERROR_STRING	= 256
};


typedef int MPI_Copy_function(MPI_Comm, int, void *, void *, void *, int *);
typedef int MPI_Delete_function(MPI_Comm, int, void *, void *);
typedef void MPI_Handler_function(MPI_Comm *, int *, ...);
typedef void MPI_User_function(void *, void *, int *, MPI_Datatype *); 

MPI_Copy_function		MPI_NULL_COPY_FN, MPI_DUP_FN;
MPI_Delete_function		MPI_NULL_DELETE_FN;


/*************************************/
/* MPI-1 bindings, sorted by chapter */
/*************************************/


/* 3.2 */

inline void err_report() {
  std::cerr << "MPI Function stub!" << std::endl ;
#ifndef NO_CSTDLIB
  std::abort() ;
#else
  abort() ;
#endif
}
inline  int  MPI_Send(void *, int, MPI_Datatype, int, int, MPI_Comm) {err_report(); return -1;}

inline  int  MPI_Recv(void *, int, MPI_Datatype, int, int, MPI_Comm, MPI_Status *) {err_report(); return -1;}

inline int  MPI_Get_count(MPI_Status *, MPI_Datatype, int *){return 0 ;}

/* 3.4 */

inline int  MPI_Bsend(void *, int, MPI_Datatype, int, int, MPI_Comm){err_report(); return -1;}

inline int  MPI_Ssend(void *, int, MPI_Datatype, int, int, MPI_Comm){err_report(); return -1;}

inline int  MPI_Rsend(void *, int, MPI_Datatype, int, int, MPI_Comm){err_report(); return -1;}

/* 3.6 */

inline int  MPI_Buffer_attach(void *, int){err_report(); return -1;}

inline int  MPI_Buffer_detach(void *, int *){err_report(); return -1;}

/* 3.7 */

inline int  MPI_Isend(void *, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *){err_report(); return -1;}

inline int  MPI_Ibsend(void *, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *){err_report(); return -1;}

inline int  MPI_Issend(void *, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *){err_report(); return -1;}

inline int  MPI_Irsend(void *, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *){err_report(); return -1;}

inline int  MPI_Irecv(void *, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *){err_report(); return -1;}

inline int  MPI_Wait(MPI_Request *, MPI_Status *){err_report(); return -1;}

inline int  MPI_Test(MPI_Request *, int *, MPI_Status *){err_report(); return -1;}

inline int  MPI_Request_free(MPI_Request *){err_report(); return -1;}

inline int  MPI_Waitany(int, MPI_Request *, int *, MPI_Status *){err_report(); return -1;}

inline int  MPI_Testany(int, MPI_Request *, int *, int *, MPI_Status *){err_report(); return -1;}

inline int  MPI_Waitall(int, MPI_Request *, MPI_Status *){err_report(); return -1;}

inline int  MPI_Testall(int, MPI_Request *, int *, MPI_Status *){err_report(); return -1;}

inline int  MPI_Waitsome(int, MPI_Request *, int *, int *, MPI_Status *){err_report(); return -1;}

inline int  MPI_Testsome(int, MPI_Request *, int *, int *, MPI_Status *){err_report(); return -1;}

/* 3.8 */

inline int  MPI_Iprobe(int, int, MPI_Comm, int *, MPI_Status *){err_report(); return -1;}

inline int  MPI_Probe(int, int, MPI_Comm, MPI_Status *){err_report(); return -1;}

inline int  MPI_Cancel(MPI_Request *){err_report(); return -1;}

inline int  MPI_Test_cancelled(MPI_Status *, int *){err_report(); return -1;}

/* 3.9 */

inline int  MPI_Send_init(void *, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *){err_report(); return -1;}

inline int  MPI_Bsend_init(void *, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *){err_report(); return -1;}

inline int  MPI_Ssend_init(void *, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *){err_report(); return -1;}

inline int  MPI_Rsend_init(void *, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *){err_report(); return -1;}

inline int  MPI_Recv_init(void *, int, MPI_Datatype, int, int, MPI_Comm, MPI_Request *){err_report(); return -1;}

inline int  MPI_Start(MPI_Request *){err_report(); return -1;}

inline int  MPI_Startall(int, MPI_Request *){err_report(); return -1;}

/* 3.10 */

inline int  MPI_Sendrecv(void *, int, MPI_Datatype, int, int, void *, int, MPI_Datatype, int, int, MPI_Comm, MPI_Status *){err_report(); return -1;}

inline int  MPI_Sendrecv_replace(void *, int, MPI_Datatype, int, int, int, int, MPI_Comm, MPI_Status *){err_report(); return -1;}

/* 3.12 */

inline int  MPI_Type_contiguous(int, MPI_Datatype, MPI_Datatype *){err_report(); return -1;}

inline int  MPI_Type_vector(int, int, int, MPI_Datatype, MPI_Datatype *){err_report(); return -1;}

inline int  MPI_Type_hvector(int, int, MPI_Aint, MPI_Datatype, MPI_Datatype *){err_report(); return -1;}

inline int  MPI_Type_indexed(int, int *, int *, MPI_Datatype, MPI_Datatype *){err_report(); return -1;}

inline int  MPI_Type_hindexed(int, int *, MPI_Aint *, MPI_Datatype, MPI_Datatype *){err_report(); return -1;}

inline int  MPI_Type_struct(int, int *, MPI_Aint *, MPI_Datatype *, MPI_Datatype *){err_report(); return -1;}

inline int  MPI_Address(void *, MPI_Aint *){err_report(); return -1;}

inline int  MPI_Type_extent(MPI_Datatype, MPI_Aint *){err_report(); return -1;}

inline int  MPI_Type_size(MPI_Datatype, int *){err_report(); return -1;}

inline int  MPI_Type_lb(MPI_Datatype, MPI_Aint *){err_report(); return -1;}

inline int  MPI_Type_ub(MPI_Datatype, MPI_Aint *){err_report(); return -1;}

inline int  MPI_Type_commit(MPI_Datatype *){err_report(); return -1;}

inline int  MPI_Type_free(MPI_Datatype *){err_report(); return -1;}

inline int  MPI_Get_elements(MPI_Status *, MPI_Datatype, int *){err_report(); return -1;}

/* 3.13 */

inline int  MPI_Pack(void *buf, int count, MPI_Datatype datatype,
                     void *packbuf, int packsize, int *packpos,
                     MPI_Comm comm){
  int size = MPI_GET_TYPE_SIZE(datatype) ;
  size *= count ;
  memcpy(((char *)packbuf)+*packpos,buf,size) ;
  *packpos += size;
  return 0 ;
}

inline int  MPI_Unpack(void *packbuf, int packsize, int *packpos,
                       void *buf, int count, MPI_Datatype datatype,
                       MPI_Comm comm) {
  int size = MPI_GET_TYPE_SIZE(datatype) ;
  size *= count ;

  memcpy(buf,((char *)packbuf)+*packpos,size) ;
  *packpos += size;
  return 0 ;
}

inline int  MPI_Pack_size(int incount, MPI_Datatype datatype , MPI_Comm comm, int *size){
  *size = MPI_GET_TYPE_SIZE(datatype)*incount ;
  return 0 ;
}

/* 4.3 */

inline int  MPI_Barrier(MPI_Comm){ return 0;}

/* 4.4 */

inline int  MPI_Bcast(void *, int, MPI_Datatype, int, MPI_Comm){ return 0 ;}

/* 4.5 */

inline int  MPI_Gather(void *sendbuf, int sendcnt, MPI_Datatype sendtype,
                       void *recvbuf, int recvcnt, MPI_Datatype recvtype,
                       int root , MPI_Comm comm)
{  memcpy(recvbuf,sendbuf,MPI_GET_TYPE_SIZE(recvtype)*recvcnt) ; return 0 ; }


inline int  MPI_Gatherv(void *, int, MPI_Datatype, void *, int *, int *, MPI_Datatype, int, MPI_Comm){err_report(); return -1;}

/* 4.6 */

inline int  MPI_Scatter(void *, int, MPI_Datatype, void *, int, MPI_Datatype, int, MPI_Comm){err_report(); return -1;}

inline int  MPI_Scatterv(void *, int *, int *, MPI_Datatype, void *, int, MPI_Datatype, int, MPI_Comm){err_report(); return -1;}

/* 4.7 */

inline int  MPI_Allgather(void *sendbuf, int sendcnt, MPI_Datatype sendtype,
                          void *recvbuf, int recvcnt, MPI_Datatype recvtype,
                          MPI_Comm comm )
{memcpy(recvbuf,sendbuf,MPI_GET_TYPE_SIZE(recvtype)*recvcnt) ; return 0 ; }

inline int  MPI_Allgatherv(void *sendbuf, int sendcnt, MPI_Datatype sendtype,
                           void *recvbuf, int *recvcounts, int *displs,
                           MPI_Datatype recvtype, MPI_Comm comm)
{ const int sz = MPI_GET_TYPE_SIZE(recvtype) ;
 memcpy(((char *)recvbuf)+displs[0]*sz, sendbuf,sz*recvcounts[0]);
 return 0 ; } 

/* 4.8 */

inline int  MPI_Alltoall(void *send, int ssz, MPI_Datatype dts, void *recv, int rsz, MPI_Datatype dtr, MPI_Comm comm){
  memcpy(recv,send,MPI_GET_TYPE_SIZE(dts)*rsz) ;
  return 0;
}

inline int  MPI_Alltoallv(void *sendbuf, int *sendcnts, int *sdispls,
                          MPI_Datatype sendtype,
                          void *recvbuf, int *recvcnts, int *rdispls,
                          MPI_Datatype recvtype, MPI_Comm comm)
{ const int sz = MPI_GET_TYPE_SIZE(recvtype) ;
 memcpy(((char *)recvbuf)+rdispls[0]*sz,((char *)sendbuf)+sdispls[0]*sz,
        recvcnts[0]*sz) ;
 return 0 ;
}

/* 4.9 */

inline int  MPI_Reduce(void *sendbuf, void *recvbuf, int count,
                       MPI_Datatype type, MPI_Op op, int root, MPI_Comm comm)
{memcpy(recvbuf,sendbuf,MPI_GET_TYPE_SIZE(type)*count) ; return 0 ; } 

inline int  MPI_Op_create(MPI_User_function *, int, MPI_Op *){err_report(); return -1;}

inline int  MPI_Op_free(MPI_Op *){err_report(); return -1;}

inline int  MPI_Allreduce(void *sendbuf, void *recvbuf, int count,
                          MPI_Datatype type, MPI_Op op , MPI_Comm comm)
{ memcpy(recvbuf,sendbuf,MPI_GET_TYPE_SIZE(type)*count) ; return 0 ; } 


/* 4.10 */

inline int  MPI_Reduce_scatter(void *, void *, int *, MPI_Datatype, MPI_Op, MPI_Comm){
  err_report() ; return -1 ;
}

/* 4.11 */

inline int  MPI_Scan(void *, void *, int, MPI_Datatype, MPI_Op, MPI_Comm){err_report(); return -1;}

/* 5.3 */

inline int  MPI_Group_size(MPI_Group, int *){err_report(); return -1;}

inline int  MPI_Group_rank(MPI_Group, int *){err_report(); return -1;}

inline int  MPI_Group_translate_ranks(MPI_Group, int, int *, MPI_Group, int *){err_report(); return -1;}

inline int  MPI_Group_compare(MPI_Group, MPI_Group, int *){err_report(); return -1;}

inline int  MPI_Comm_group(MPI_Comm, MPI_Group *){err_report(); return -1;}

inline int  MPI_Group_union(MPI_Group, MPI_Group, MPI_Group *){err_report(); return -1;}

inline int  MPI_Group_intersection(MPI_Group, MPI_Group, MPI_Group *){err_report(); return -1;}

inline int  MPI_Group_difference(MPI_Group, MPI_Group, MPI_Group *){err_report(); return -1;}

inline int  MPI_Group_incl(MPI_Group, int, int *, MPI_Group *){err_report(); return -1;}

inline int  MPI_Group_excl(MPI_Group, int, int *, MPI_Group *){err_report(); return -1;}

inline int  MPI_Group_range_incl(MPI_Group, int, int [][3], MPI_Group *){err_report(); return -1;}

inline int  MPI_Group_range_excl(MPI_Group, int, int [][3], MPI_Group *){err_report(); return -1;}

inline int  MPI_Group_free(MPI_Group *){err_report(); return -1;}

/* 5.4 */

inline int  MPI_Comm_size(MPI_Comm, int *val){*val = 1 ; return 0;}

inline int  MPI_Comm_rank(MPI_Comm, int *val){*val = 0 ; return 0;}

inline int  MPI_Comm_compare(MPI_Comm, MPI_Comm, int *){err_report(); return -1;}

inline int  MPI_Comm_dup(MPI_Comm, MPI_Comm *){err_report(); return -1;}

inline int  MPI_Comm_create(MPI_Comm, MPI_Group, MPI_Comm *){err_report(); return -1;}

inline int  MPI_Comm_split(MPI_Comm, int, int, MPI_Comm *){err_report(); return -1;}

inline int  MPI_Comm_free(MPI_Comm *){err_report(); return -1;}

/* 5.6 */

inline int  MPI_Comm_test_inter(MPI_Comm, int *){err_report(); return -1;}

inline int  MPI_Comm_remote_size(MPI_Comm, int *){err_report(); return -1;}

inline int  MPI_Comm_remote_group(MPI_Comm, MPI_Group *){err_report(); return -1;}

inline int  MPI_Intercomm_create(MPI_Comm, int, MPI_Comm, int, int, MPI_Comm *){err_report(); return -1;}

inline int  MPI_Intercomm_merge(MPI_Comm, int, MPI_Comm *){err_report(); return -1;}

/* 5.7 */

inline int  MPI_Keyval_create(MPI_Copy_function *, MPI_Delete_function *, int *, void *){err_report(); return -1;}

inline int  MPI_Keyval_free(int *){err_report(); return -1;}

inline int  MPI_Attr_put(MPI_Comm, int, void *){err_report(); return -1;}

inline int  MPI_Attr_get(MPI_Comm, int, void *, int *){err_report(); return -1;}

inline int  MPI_Attr_delete(MPI_Comm, int){err_report(); return -1;}

/* 6.5 */

inline int  MPI_Cart_create(MPI_Comm, int, int *, int *, int, MPI_Comm *){err_report(); return -1;}

inline int  MPI_Dims_create(int, int, int *){err_report(); return -1;}

inline int  MPI_Graph_create(MPI_Comm, int, int *, int *, int, MPI_Comm *){err_report(); return -1;}

inline int  MPI_Topo_test(MPI_Comm, int *){err_report(); return -1;}

inline int  MPI_Graphdims_get(MPI_Comm, int *, int *){err_report(); return -1;}

inline int  MPI_Graph_get(MPI_Comm, int, int, int *, int *){err_report(); return -1;}

inline int  MPI_Cartdim_get(MPI_Comm, int *){err_report(); return -1;}

inline int  MPI_Cart_get(MPI_Comm, int, int *, int *, int *){err_report(); return -1;}

inline int  MPI_Cart_rank(MPI_Comm, int *, int *){err_report(); return -1;}

inline int  MPI_Cart_coords(MPI_Comm, int, int, int *){err_report(); return -1;}

inline int  MPI_Graph_neighbors_count(MPI_Comm, int, int *){err_report(); return -1;}

inline int  MPI_Graph_neighbors(MPI_Comm, int, int, int *){err_report(); return -1;}

inline int  MPI_Cart_shift(MPI_Comm, int, int, int *, int *){err_report(); return -1;}

inline int  MPI_Cart_sub(MPI_Comm, int *, MPI_Comm *){err_report(); return -1;}

inline int  MPI_Cart_map(MPI_Comm, int, int *, int *, int *){err_report(); return -1;}

inline int  MPI_Graph_map(MPI_Comm, int, int *, int *, int *){err_report(); return -1;}

/* 7.1 */

inline int  MPI_Get_processor_name(char *, int *){err_report(); return -1;}

/* 7.2 */

inline int  MPI_Errhandler_create(MPI_Handler_function *, MPI_Errhandler *){return 0 ;}

inline int  MPI_Errhandler_set(MPI_Comm, MPI_Errhandler){return 0;}

inline int  MPI_Errhandler_get(MPI_Comm, MPI_Errhandler *){err_report(); return -1;}

inline int  MPI_Errhandler_free(MPI_Errhandler *){err_report(); return -1;}

inline int  MPI_Error_string(int, char *, int *){err_report(); return -1;}

/* 7.3 */

inline int  MPI_Error_class(int, int *){err_report(); return -1;}

/* 7.4 */

inline double  MPI_Wtime(void){return 0; }

inline double  MPI_Wtick(void){return 0; }

/* 7.5 */

inline int  MPI_Init(int *, char ***){ return 0;}

inline int  MPI_Finalize(void){return 0;}

inline int  MPI_Initialized(int *){return 0;}

inline int  MPI_Abort(MPI_Comm, int){
#ifndef NO_CSTDLIB
  std::exit(-1);
#else
  exit(-1) ;
#endif
  return -1;}

/* 8.3 */

inline int  MPI_Pcontrol(int, ...){err_report(); return -1;}

/********************/
/* MPI-1.2 bindings */
/********************/

#define MPI_VERSION		1
#define MPI_SUBVERSION		2

inline int  MPI_Get_version(int *, int *){return 1;}


/*************************************/
/* MPI-2 bindings, sorted by chapter */
/*************************************/

/* 4.10 */

typedef unsigned int          MPI_Info;

enum {
	MPI_INFO_NULL		= 0,
	MPI_MAX_INFO_KEY        = 255,
	MPI_MAX_INFO_VAL        = 1024
};

enum {
	MPI_FUNDAMENTAL		= -1
};

inline int  MPI_Info_create(MPI_Info *){err_report(); return -1;}

inline int MPI_Info_delete(MPI_Info, char *){err_report(); return -1;}

inline int MPI_Info_dup(MPI_Info, MPI_Info *){err_report(); return -1;}

inline int MPI_Info_free(MPI_Info *){err_report(); return -1;}

inline int MPI_Info_get(MPI_Info, char *, int, char *, int *){err_report(); return -1;}

inline int MPI_Info_get_nkeys(MPI_Info, int *){err_report(); return -1;}

inline int MPI_Info_get_nthkey(MPI_Info, int, char *){err_report(); return -1;}

inline int MPI_Info_get_valuelen(MPI_Info, char *, int *, int *){err_report(); return -1;}

inline int MPI_Info_set(MPI_Info, char *, char *){err_report(); return -1;}

/* 4.11 */

inline int MPI_Alloc_mem(MPI_Aint,MPI_Info,void *){err_report(); return -1;}

inline int MPI_Free_mem(void *){err_report(); return -1;}

/* 4.12 */

typedef int MPI_Fint;

inline MPI_Fint MPI_Info_c2f(MPI_Info){err_report(); return -1;}

inline MPI_Info MPI_Info_f2c(MPI_Fint){err_report(); return 0;}

inline MPI_Fint MPI_Comm_c2f(MPI_Comm){err_report(); return 0;}

inline MPI_Comm MPI_Comm_f2c(MPI_Fint){err_report(); return 0;}

inline MPI_Fint MPI_Type_c2f(MPI_Datatype){err_report(); return 0;}

inline MPI_Datatype MPI_Type_f2c(MPI_Fint){err_report(); return 0;}

inline MPI_Fint MPI_Group_c2f(MPI_Group){err_report(); return 0;}

inline MPI_Group MPI_Group_f2c(MPI_Fint){err_report(); return 0;}

inline MPI_Fint MPI_Request_c2f(MPI_Request){err_report(); return 0;}

inline MPI_Request MPI_Request_f2c(MPI_Fint){err_report(); return 0;}

inline MPI_Fint MPI_Op_c2f(MPI_Op){err_report(); return 0;}

inline MPI_Op MPI_Op_f2c(MPI_Fint){err_report(); return 0;}


/* 4.14 */

inline int  MPI_Type_create_hvector(int, int, MPI_Aint, MPI_Datatype, MPI_Datatype *){err_report(); return 0;}

inline int  MPI_Type_create_hindexed(int, int *, MPI_Aint *, MPI_Datatype, MPI_Datatype *){err_report(); return 0;}

inline int  MPI_Type_create_struct(int, int *, MPI_Aint *, MPI_Datatype *, MPI_Datatype *){err_report(); return 0;}

inline int MPI_Get_address(void *, MPI_Aint *){err_report(); return 0;}

/* 5.3 */

#define MPI_ARGV_NULL              ((char **)NULL)
#define MPI_ARGVS_NULL              ((char ***)NULL)
#define MPI_ERRCODES_IGNORE        ((int *)NULL)

inline int MPI_Comm_spawn(char *, char **, int, MPI_Info, int, MPI_Comm, MPI_Comm *, int *){err_report(); return 0;}
inline int MPI_Comm_spawn_multiple(int , char **, char ***, int *, MPI_Info *,
                            int , MPI_Comm, MPI_Comm *, int *){err_report(); return 0;} 
inline int MPI_Comm_get_parent(MPI_Comm *){err_report(); return 0;}

/* 6 */

/* MPI one-sided is supported only under ABI 64.  */
#if 	!_ABIN32

typedef unsigned int MPI_Win;

enum {
	MPI_WIN_NULL 	= 0
};

/* referenced in 4.12 of MPI-2 standard with other transfer of handles functions */

inline MPI_Fint MPI_Win_c2f(MPI_Win){err_report(); return 0;}

inline MPI_Win MPI_Win_f2c(MPI_Fint){err_report(); return 0;}

/* 6.2 */

inline int  MPI_Win_create(void *, MPI_Aint, int, MPI_Info, MPI_Comm, MPI_Win *){err_report(); return 0;}

inline int  MPI_Win_fence(int, MPI_Win){err_report(); return 0;}

inline int  MPI_Win_free(MPI_Win *){err_report(); return 0;}

/* 6.3 */


inline int  MPI_Put(void *, int, MPI_Datatype, int, MPI_Aint, int, MPI_Datatype,
	MPI_Win){err_report(); return 0;}

inline int  MPI_Get(void *, int, MPI_Datatype, int, MPI_Aint, int, MPI_Datatype,
	MPI_Win){err_report(); return 0;}

inline int  MPI_Accumulate(void *, int, MPI_Datatype, int, MPI_Aint, int,
	MPI_Datatype, MPI_Op, MPI_Win){err_report(); return 0;}

/* 6.4 */

enum {
	MPI_MODE_NOCHECK	= 1,
	MPI_MODE_NOSTORE	= 2,
	MPI_MODE_NOPUT		= 4,
	MPI_MODE_NOPRECEDE	= 8,
	MPI_MODE_NOSUCCEED 	= 16
};


#endif	/* !_ABIN32 */


/* 7.5 */

inline int  MPI_Type_get_envelope(MPI_Datatype, int *, int *, int *, int *){err_report(); return 0;}

inline int  MPI_Type_get_contents(MPI_Datatype, int, int, int, int *, MPI_Aint *, MPI_Datatype *){err_report(); return 0;}

/* 7.8 */

inline int  MPI_Type_dup(MPI_Datatype, MPI_Datatype *){err_report(); return 0;}

/* 8.6 */

enum {
	MPI_COMBINER_NAMED      	= (0),
	MPI_COMBINER_CONTIGUOUS 	= 0,
	MPI_COMBINER_VECTOR     	= 1,
	MPI_COMBINER_HVECTOR    	= 2,
	MPI_COMBINER_INDEXED    	= 3,
	MPI_COMBINER_HINDEXED   	= 4,
	MPI_COMBINER_STRUCT     	= 5,
	MPI_COMBINER_DARRAY		= 6,
	MPI_COMBINER_DUP		= 7,
	MPI_COMBINER_F90_COMPLES	= 8,
	MPI_COMBINER_F90_INTEGER	= 9,
	MPI_COMBINER_F90_REAL		= 10,
	MPI_COMBINER_HINDEXED_INTEGER	= 11,
	MPI_COMBINER_HVECTOR_INTEGER	= 12,
	MPI_COMBINER_INDEXED_BLOCK	= 13,
	MPI_COMBINER_RESIZED		= 14,
	MPI_COMBINER_STRUCT_INTEGER	= 15,
	MPI_COMBINER_SUBARRAY		= 16
};

/* 8.7 */

inline int  MPI_Init_thread(int *, char ***, int, int *){err_report(); return 0;}
inline int  MPI_Query_thread(int *){err_report(); return 0;}
inline int  MPI_Is_thread_main(int *){err_report(); return 0;}

enum {
	MPI_THREAD_SINGLE		= 0,
	MPI_THREAD_FUNNELED		= 1,
	MPI_THREAD_SERIALIZED		= 2,
	MPI_THREAD_MULTIPLE		= 3
};

/* 9.6 */

inline int  MPI_Finalized(int *){err_report(); return 0;}

#endif	/* MPI_H_INCLUDED */
