#include <sys/stat.h>
#include <unistd.h>
#include <iostream>
#include <fstream>
#include <Tools/xdr.h>
#include <algorithm>
#include <mpi.h>

int myId = 0;
int numProcessors = 1;

using std::cout ;
using std::endl ;
using std::cerr ;
using std::ios ;

typedef int fourtuple[4] ;
typedef int fivetuple[5] ;
typedef int sixtuple[6] ;
typedef int eighttuple[8] ;
template< class T > class autoArray {
  T *ptr ;
public:
  explicit autoArray(T* p=0) throw() { ptr = p; }
  ~autoArray() throw() { delete[] ptr ;}

  autoArray(autoArray &a) throw() { ptr = a.ptr; a.ptr = 0 ; }
  autoArray &operator=(autoArray &a) throw() {ptr = a.ptr; a.ptr = 0 ; }

  T& operator*() const throw() { return *ptr; }
  T* operator->() const throw() { return ptr; }
  T& operator[](size_t i) throw() { return ptr[i]; } 
} ;

    
//---------------------Array----------------------//
template <class T,size_t n> class Array {
  T x[n] ;
public:
  typedef T * iterator ;
  typedef const T * const_iterator ;
    
  Array() {} ;
  Array(const Array<T,n> &v)
  { for(size_t i=0;i<n;++i) x[i] = v.x[i] ; } 
  Array<T,n> &operator=(const Array<T,n> &v)
  { for(size_t i=0;i<n;++i) x[i] = v.x[i] ; return *this ; } 
    
  Array<T,n> &operator +=(const Array<T,n> &v)
  { for(size_t i=0;i<n;++i) x[i] += v.x[i] ; return *this ; }
  Array<T,n> &operator -=(const Array<T,n> &v)
  { for(size_t i=0;i<n;++i) x[i] -= v.x[i] ; return *this ; }
  Array<T,n> &operator *=(const Array<T,n> &v)
  { for(size_t i=0;i<n;++i) x[i] *= v.x[i] ; return *this ; }
  Array<T,n> &operator /=(const Array<T,n> &v)
  { for(size_t i=0;i<n;++i) x[i] /= v.x[i] ; return *this ; }
    
  T &operator[](size_t indx) { return x[indx]; }
  const T &operator[](size_t indx) const { return x[indx] ; }
    
  iterator begin() { return &x[0] ; }
  iterator end() { return begin()+n ; }
  const_iterator begin() const { return &x[0] ; }
  const_iterator end() const { return begin()+n ; }
    
  size_t size() const  { return n ; }
} ;

template <class T,size_t n> inline std::ostream &
operator<<(std::ostream &s, const Array<T,n> &v) {
  for(int i=0;i<n;++i)
    s << v[i] << ' ' ;
  return s ;
}

template <class T,size_t n> inline std::istream &
operator>>(std::istream &s, Array<T,n> &v) {
  for(int i=0;i<n;++i)
    s >> v[i] ;
  return s ;
}


bool reverse_byteorder = false ;

void check_order() {
  static int test = 15 ;
  char *p = (char *)&test ;
  if(int(*p) == test) {
    reverse_byteorder = true ;
  }
}


void ug_io_reverse_byte_order
(void * Data,
 size_t Size,
 int Number)

{
 
  /*
   * Set file format and host to big or little endian byte ordering.
   * 
   */

  char *Data_Byte_Ptr;
  char Temp_Data_Byte;

  int Byte_Index, Index, Number_of_Bytes, Reverse_Byte_Index;

  Number_of_Bytes = int(Size);

  Data_Byte_Ptr = (char *) Data;

  for (Index = 0; Index < Number; ++Index)
    {
      Reverse_Byte_Index = Number_of_Bytes;

      for (Byte_Index = 0; Byte_Index < Number_of_Bytes/2; ++Byte_Index)
        {
          --Reverse_Byte_Index;

          Temp_Data_Byte = Data_Byte_Ptr[Byte_Index];

          Data_Byte_Ptr[Byte_Index] = Data_Byte_Ptr[Reverse_Byte_Index];

          Data_Byte_Ptr[Reverse_Byte_Index] = Temp_Data_Byte;
        }

      Data_Byte_Ptr += Number_of_Bytes;
    }

  return;
}

typedef Array<int,4> tri_info ;
typedef Array<int,5> quad_info ;

inline bool tri_sort(const tri_info &a1,const tri_info &a2) {
  return a1[0] < a2[0] || (a1[0]==a2[0] && a1[1]<a2[1]) ||
    (a1[0]==a2[0] && a1[1] == a2[1] && a1[2]<a2[2]) ||
    (a1[0]==a2[0] && a1[1] == a2[1] && a1[2] == a2[2] && a1[3] < a2[3]) ;
}

inline bool quad_sort(const quad_info &a1,const quad_info &a2) {
  return a1[0] < a2[0] || (a1[0]==a2[0] && a1[1]<a2[1]) ||
    (a1[0]==a2[0] && a1[1] == a2[1] && a1[2]<a2[2]) ||
    (a1[0]==a2[0] && a1[1] == a2[1] && a1[2] == a2[2] && a1[3] < a2[3]) ||
    (a1[0]==a2[0] && a1[1] == a2[1] && a1[2] == a2[2] && a1[3] == a2[3]
     && a1[4] < a2[4]) ;
}


inline bool bnd_tri_sort(const Array<int,5> &a1,
                         const Array<int,5> &a2) {
  return a1[4] < a2[4] || (a1[4]==a2[4] && a1[3]<a2[3]) ||  (a1[4]==a2[4] && a1[3]==a2[3] && a1[2]< a2[2]) || (a1[4]==a2[4] && a1[3]==a2[3] && a1[2] == a2[2] && a1[1] < a2[1]) || (a1[4]==a2[4] && a1[3]==a2[3] && a1[2] == a2[2] && a1[1] == a2[1] && a1[0] < a2[0]); }


inline bool bnd_quad_sort(const Array<int,6> &a1,
                          const Array<int,6> &a2) {
  return a1[5] < a2[5] || (a1[5] == a2[5] && a1[4] < a2[4]) || (a1[5] == a2[5] && a1[4] == a2[4] && a1[3] < a2[3] ) || (a1[5] == a2[5] && a1[4] == a2[4] && a1[3] ==  a2[3] && a1[2] < a2[2]) ||  (a1[5] == a2[5] && a1[4] == a2[4] && a1[3] ==  a2[3] && a1[2] == a2[2] && a1[1]<a2[1] ) || (a1[5] == a2[5] && a1[4] == a2[4] && a1[3] ==  a2[3] && a1[2] == a2[2] && a1[1]==a2[1] && a1[0]< a2[0]);}


int Server(int ac, char *av[],int , int );
void Client(int , int );

/***************************************************/
int main(int ac, char *av[]) {
 

  MPI_Init(&ac,&av) ;
 

  /* Get the number of processors and my processor identification */
  MPI_Comm_size( MPI_COMM_WORLD, &numProcessors) ;
  MPI_Comm_rank( MPI_COMM_WORLD, &myId) ;

  if(myId == 0) {
    // Processor 0 runs the server code
    Server(ac,av, myId, numProcessors) ;    
  } else {
    // all other processors run the client code.
    Client(myId, numProcessors);
  }

  // All MPI programs must call this before exit  
  MPI_Finalize() ;
  return 0 ; 
}  
int Server(int ac, char *av[], int Id, int numProcessors)
{ 
  int myId, npes; 
  myId = Id;
  npes = numProcessors;
 
  const char *filename ;
  std::string tmp_str ;
  bool binary = false ;
  if(ac == 3) {
    tmp_str.append(av[2]) ;
    if(!strcmp(av[1],"-b")) 
      binary = 1 ;
    else {
      cerr << "Right now the only option supported is '-b' which is the binary ugrid format" << endl ;
      exit(-1) ;
    }
  } else if(ac == 2) {
    tmp_str.append(av[1]) ;
  } else {
    cerr << "solidMesh2xdr requires one argument" << endl 
         << " (the -b flag may be specified for binary files)" << endl;
    exit(-1) ;
  }

  check_order() ;
 
  int loc = 0;
  loc = tmp_str.find('.') ;
  std::string new_str = tmp_str.substr(0, loc) ;
  filename = new_str.c_str() ;

  FILE* IFP ;
  char buf[512] ;
  if(!binary) {
    struct stat fstat ;
    sprintf(buf,"%s.ugrid",filename) ;
    if(stat(buf,&fstat)<0) {
      binary = true ;
    }
  }
    
  if(!binary)
    sprintf(buf,"%s.ugrid",filename) ;
  else
    sprintf(buf,"%s.b8.ugrid",filename) ;
  if(!binary)
    IFP = fopen(buf, "r") ;
  else
    IFP = fopen(buf, "rb") ;
  if(IFP == NULL) {
    cerr << "can't open '" << buf << "'" << endl ;
    exit(-1) ;
  }
  char out_buf[512] ;
  sprintf(out_buf,"%s.xdr",filename) ;
  FILE *FP = fopen(out_buf, "w") ;
  if(FP == NULL) {
    cerr << "can't open " << out_buf <<  endl ;
    return(-1);
  }
  
  MPI_Status stat;
  int num_nodes, num_sf_trias, num_sf_quads ;
  int num_vol_tets, num_vol_pents5, num_vol_pents6, num_vol_hexs ;
          
  if(!binary) {
    fscanf(IFP, "%d%d%d", &num_nodes, & num_sf_trias, & num_sf_quads) ;
    fscanf(IFP, "%d%d%d%d", &num_vol_tets, &num_vol_pents5, &num_vol_pents6, &num_vol_hexs) ;
  }
  else {
    fread(&num_nodes, sizeof(int), 1, IFP) ;
    if(reverse_byteorder)
      ug_io_reverse_byte_order(&num_nodes,sizeof(int),1) ;
    fread(&num_sf_trias, sizeof(int), 1, IFP) ;
    if(reverse_byteorder)
      ug_io_reverse_byte_order(&num_sf_trias,sizeof(int),1) ;
    fread(&num_sf_quads, sizeof(int), 1, IFP) ;
    if(reverse_byteorder)
      ug_io_reverse_byte_order(&num_sf_quads,sizeof(int),1) ;
    fread(&num_vol_tets, sizeof(int), 1, IFP) ;
    if(reverse_byteorder)
      ug_io_reverse_byte_order(&num_vol_tets,sizeof(int),1) ;
    fread(&num_vol_pents5, sizeof(int), 1, IFP) ;
    if(reverse_byteorder)
      ug_io_reverse_byte_order(&num_vol_pents5,sizeof(int),1) ;
    fread(&num_vol_pents6, sizeof(int), 1, IFP) ;
    if(reverse_byteorder)
      ug_io_reverse_byte_order(&num_vol_pents6,sizeof(int),1) ;
    fread(&num_vol_hexs, sizeof(int), 1, IFP) ;
    if(reverse_byteorder)
      ug_io_reverse_byte_order(&num_vol_hexs,sizeof(int),1) ;
  }

  if(myId == 0) {
    cout << " Number of  nodes = " << num_nodes << endl ;
    cout << " Number of surface triangles = " << num_sf_trias << endl ;
    cout << " Number of surface quads = " << num_sf_quads << endl ;
    cout << " Number of volume tetrahedras = " << num_vol_tets << endl ;
    cout << " Number of volume pents_5 = " << num_vol_pents5 << endl ;
    cout << " Number of volume pents6 = " << num_vol_pents6 << endl ;
    cout << " Number of volume hexahedra = " << num_vol_hexs << endl ;
  }
  
  size_t num_quad_faces = num_sf_quads + num_vol_pents5 + num_vol_pents6*3 + num_vol_hexs*6 ;
 
  if((num_quad_faces & 1) == 1) {
    cerr << "not all quad faces can pair!" << endl ;
    exit(-1) ;
  }
  // We counted face pairs, this is twice the final number
  num_quad_faces = num_quad_faces >> 1 ;

  size_t num_tri_faces = num_sf_trias + num_vol_tets*4+num_vol_pents5*4+num_vol_pents6*2 ;

  //  if(myId == 0)
  //    cout<<"num_tri_faces  "<<num_tri_faces<<endl;

  if((num_tri_faces & 1) == 1) {
    cerr << "not all trianglular faces can pair!" << endl ;
    exit(-1) ;
  }

  num_tri_faces = num_tri_faces >> 1 ;
 
  int Info[7] ;
 
  Info[0] = num_nodes ; 
  Info[1] = num_sf_trias ;
  Info[2] = num_sf_quads ; 
  Info[3] = num_vol_tets ; 
  Info[4] = num_vol_pents5 ;
  Info[5] = num_vol_pents6 ;
  Info[6] = num_vol_hexs ; 
  

  // Broadcasting information to all the other processor's
  MPI_Bcast( &Info, 7, MPI_INT, 0, MPI_COMM_WORLD);

   
  int ndim = 3 ;
  int nzones = 1 ;
  int npatch = 0 ;
  int ncells = num_vol_tets + num_vol_pents5 + num_vol_pents6 + num_vol_hexs ;
  int nfaces = num_quad_faces+num_tri_faces ;
  int max_ppf = 3 ;
  int max_fpc = 4 ;
  if(num_vol_pents5>0 || num_vol_pents6>0 || num_vol_hexs>0)
    max_ppf = 4 ;
  if(num_vol_pents5>0 || num_vol_pents6>0)
    max_fpc = 5 ;
  if(num_vol_hexs>0)
    max_fpc = 6 ;
      
  XDR xdr_handle ;
  xdrstdio_create(&xdr_handle, FP, XDR_ENCODE) ;

  xdr_int(&xdr_handle, &ndim) ;
  xdr_int(&xdr_handle, &nzones) ;
  xdr_int(&xdr_handle, &npatch) ;
  xdr_int(&xdr_handle, &num_nodes) ;
  xdr_int(&xdr_handle, &nfaces) ;
  xdr_int(&xdr_handle, &ncells) ;
  xdr_int(&xdr_handle, &max_ppf) ;
  xdr_int(&xdr_handle, &max_fpc) ;

  if(myId == 0) {
    cout << "XDR grid contains " << num_nodes << " nodes, "
         << nfaces << " faces, and " << ncells << " cells" << endl ;

    cout << "copying node information..." << endl ;
  }
  
  if(!binary) {
    double ptmp ;
    for(int i = 0; i < 3 * num_nodes; ++i) {
      fscanf(IFP, "%lf", &ptmp) ;
      xdr_double(&xdr_handle, &ptmp) ;
    }
  } else {
    double ptmp[3] ;
    if(reverse_byteorder)
      for(int i = 0; i < num_nodes; ++i) {
        fread(ptmp,sizeof(double),3,IFP) ;
        ug_io_reverse_byte_order(ptmp,sizeof(double),3) ;
        xdr_double(&xdr_handle,&ptmp[0]) ;
        xdr_double(&xdr_handle,&ptmp[1]) ;
        xdr_double(&xdr_handle,&ptmp[2]) ;
      }
    else
      for(int i = 0; i < num_nodes; ++i) {
        fread(ptmp,sizeof(double),3,IFP) ;
        xdr_double(&xdr_handle,&ptmp[0]) ;
        xdr_double(&xdr_handle,&ptmp[1]) ;
        xdr_double(&xdr_handle,&ptmp[2]) ;
      }
  }
 

  /*************************************************************/
  // CALCULATING THE NUMBER OF ELEMENTS TO SEND TO EACH PROCESSOR

 
  // Calculating number of surface trias for each processor

  int num_sf_tri=0 , num_sf_tri_lastpe=0, num=0;

  if(npes==1)
    num_sf_tri=Info[1];
  else {
    num = Info[1]/npes;   
  
    if( num % 2 == 1)  {
      num_sf_tri = num - 1;
      num_sf_tri_lastpe = Info[1] - (num-1)*(npes-1);			  }
    else {
      num_sf_tri = num ;
      num_sf_tri_lastpe=num;
      if(Info[1] % npes != 0 )		   
        num_sf_tri_lastpe = Info[1] - num*(npes-1);
    }		   
  
  }
  // Calculating number of boundary quads for each processor    
  int num_sf_quad=0, num_sf_quad_lastpe=0;
 
  if(npes==1)
    num_sf_quad=Info[2];
  else {
    num = Info[2]/npes ;  
    if(num % 2 == 1) {	
      num_sf_quad = num - 1; 
      num_sf_quad_lastpe = Info[2] - (num - 1)*(npes-1); 
    } else { 
      num_sf_quad = num ;
      num_sf_quad_lastpe=num;
      if(Info[2] % npes != 0)	    
        num_sf_quad_lastpe = Info[2] - num*(npes-1) ;
    }    
  }

  // Calculating number of tets to give to each processor

  int num_tets=0, num_tets_lastpe=0;   

  if(npes==1)
    num_tets=Info[3];
  else {
    if(Info[3]==0) {
      num_tets=0;
      num_tets_lastpe=0;
    } else {
      num = Info[3]/npes;
	  
      if(Info[3] % npes != 0)  {
        num_tets = num ;
        num_tets_lastpe = Info[3] - num*(npes-1);
      } else { 
        num_tets = num; 
        num_tets_lastpe = num; 
      } 
    }
  }
  //  cout<<"num of tets at the server "<<num_tets<<endl;
 
  // Calculating number of pent5 to give to each processor for trifaces
  // This should be allocated evenly to each processor  
 
  int num_pent5=0, num_pent5_lastpe=0;
  
  if(npes==1)
    num_pent5=Info[4];
  else {
    if(Info[4]==0) {
      num_pent5 =0;
      num_pent5_lastpe = 0;
    } else {
      num = Info[4]/npes;
	 
      if(num % 2 == 1)  {
        num_pent5 = num - 1;
        num_pent5_lastpe = Info[4] - (num - 1)*(npes-1);  }
      else {
        if(Info[4] % npes != 0) {
          num_pent5 = num ;
          num_pent5_lastpe = Info[4] - num*(npes-1);
        }
      }    
    }  
  }

  // Calculating number of Pent6 to give to each processor
  int num_pent6=0, num_pent6_lastpe=0 ;
  
  if(npes==1)
    num_pent6=Info[5];
  else {      
    if(Info[5]==0) {
      num_pent6=0;
      num_pent6_lastpe=0;
    } else {     
      num = Info[5]/npes;
	  
      if(num % 2 == 1)  {
        num_pent6 = num - 1 ; 
        num_pent6_lastpe = Info[5] - (num-1)*(npes-1); 
      } else {
        num_pent6 = num ; 	
        num_pent6_lastpe = num;
        if(Info[5] % npes != 0) {     
          num_pent6 = num ; 	
          num_pent6_lastpe = Info[5] - num*(npes-1);
        } 
      }
    }
  }

  // Calculating the number of Hex to give to each processor
  int num_hex=0,num_hex_lastpe=0 ;
  
  if(npes==1)
    num_hex=Info[6];
  else {
    if(Info[6]==0) {
      num_hex=0;
      num_hex_lastpe=0;
    } else {

      num = Info[6]/npes;
	  
      if(Info[6] % npes != 0)  {
        num_hex = num ;	
        num_hex_lastpe = Info[6] - num*(npes-1); 
      } else {  
        num_hex = num;
        num_hex_lastpe = num; 
      }  
    }
  }

  // calculating the total number of volume elements recv by last PE
  // And broadcasting them to each processor for the purpose of giving proper cellnum   
 
  int total_vol_elmn = num_hex_lastpe + num_tets_lastpe + num_pent5_lastpe + num_pent6_lastpe;

  // Array used to store the volume elements used by the last processor 
  int info_cellnum[4];

  info_cellnum[0]=num_tets_lastpe;
  info_cellnum[1]=num_pent5_lastpe;
  info_cellnum[2]=num_pent6_lastpe;
  info_cellnum[3]=num_hex_lastpe;

  MPI_Bcast(&info_cellnum[0], 4, MPI_INT, 0, MPI_COMM_WORLD); 

  MPI_Bcast(&total_vol_elmn, 1, MPI_INT, 0, MPI_COMM_WORLD);  
  
  
  // Calculating the total number of Tri_faces & quad faces generated at each PE // Initializing mem for each  
  
  size_t tri_total = num_sf_tri + num_tets*4 + num_pent5*4 + num_pent6*2 ;
  size_t quad_total = num_sf_quad + num_pent5 + num_pent6*3 + num_hex*6;
 
  autoArray<tri_info> tri_faces(new tri_info[tri_total]) ;
  autoArray<quad_info> quad_faces(new quad_info[quad_total]);

  if(myId == 0) {
//    cout<<"tri_total "<<tri_total<<endl;

    cout << "reading in boundary information..." << endl ;
  }
  
  size_t tf = 0 ; // triangle and quad face pointers
  size_t qf = 0 ; // size_t is unsigned long 
  
  // Read in boundary triangles
  int ctr=0, count=0,c=0,temp=0,tag=0; 
  

  /************************************************************/
  // SENDING THE ELEMENTS TO OTHER PROCESSORS

  

  // Sending Tri element information to each processor starting with the
  // last processor first to the first processor
   
  autoArray<tri_info> tri_buffer(new tri_info[num_sf_tri_lastpe]) ;

  while(ctr < npes-1)  {
    c=0;
    if(ctr==0) {	    
      temp = num_sf_tri_lastpe ;
      //	cout<<"temp  0 "<<temp<<endl;
    } else {	 	
      temp = num_sf_tri ;
      //	cout<<"temp   "<<ctr+1<<temp<<endl;
    }		

    for(int i = 0; i < temp; ++i) {
      tri_info tmp_tria ;
      tmp_tria[3] = -100000 ;
      if(!binary)
	fscanf(IFP, "%d%d%d", &tmp_tria[0], &tmp_tria[1], &tmp_tria[2]);
      else { 
        fread(&tmp_tria[0], sizeof(int), 3, IFP) ;
        if(reverse_byteorder)
          ug_io_reverse_byte_order(&tmp_tria[0],sizeof(int),3) ;
      }        
      tri_buffer[c++] = tmp_tria ; 
    } 

    count=c*4;
    //   int d=0;

    MPI_Recv(&tag, 1, MPI_INT, npes-ctr-1, 1, MPI_COMM_WORLD, &stat); 

    if(tag==1)
      MPI_Send(&tri_buffer[0], count, MPI_INT, npes-ctr-1, 0, MPI_COMM_WORLD );  

    
    ctr++;     
  }  
  
  
  // Copying the Tri element information for processor 0 
  temp = num_sf_tri; 

//  if(myId == 0)
//    cout<<"num_sf_tri  "<<num_sf_tri<<endl;

  for(int i = 0; i < temp; ++i) {
    tri_info tmp_tria ;
    tmp_tria[3] = -100000 ;
    if(!binary)
      fscanf(IFP, "%d%d%d", &tmp_tria[0], &tmp_tria[1], &tmp_tria[2]) ;  
    else { 
      fread(&tmp_tria[0], sizeof(int), 3, IFP) ;
      if(reverse_byteorder)
        ug_io_reverse_byte_order(&tmp_tria[0],sizeof(int),3) ;
    }
    tri_faces[tf++] = tmp_tria  ;
  }    

  
  // creating buffer for sending the surface quad 

  autoArray<quad_info> quad_buffer(new quad_info[num_sf_quad_lastpe]) ;
  ctr=0;
  temp=0;  
  tag=0;
  
  // Read in boundary quads
  while(ctr < npes-1)  {
    c=0; 
    if(ctr==0) {	
      temp = num_sf_quad_lastpe ;
    } else {	
      temp = num_sf_quad ;
    }
	
    for(int i = 0; i < temp; ++i) {
      quad_info tmp_quad ;
      tmp_quad[4] = -100000 ;
      
      if(!binary)
	fscanf(IFP, "%d%d%d%d", &tmp_quad[0], &tmp_quad[1], &tmp_quad[2], &tmp_quad[3]);  
      else {
        fread(&tmp_quad[0], sizeof(int), 4, IFP) ;
        if(reverse_byteorder)
          ug_io_reverse_byte_order(&tmp_quad[0],sizeof(int),4) ;
      }
      quad_buffer[c++] = tmp_quad ;
    }  
    count = 5*c;
    
    MPI_Recv(&tag, 1, MPI_INT,  npes-ctr-1, 1, MPI_COMM_WORLD, &stat); 
    
    if(tag==1)
      MPI_Send(&quad_buffer[0], count, MPI_INT, npes-ctr-1, 0, MPI_COMM_WORLD );    
    ctr++;
    
  }
  
  // Copying the element information of quad for processor 0 

  temp = num_sf_quad;
//  if(myId == 0)
//    cout<<"num_sf_quad "<<num_sf_quad<<endl;
   
  for(int i = 0; i < temp; ++i) {
    quad_info tmp_quad ;
    tmp_quad[4] = -100000 ;
    if(!binary)
      fscanf(IFP, "%d%d%d%d", &tmp_quad[0], &tmp_quad[1], &tmp_quad[2], &tmp_quad[3]) ;  
    else {
      fread(&tmp_quad[0], sizeof(int), 4, IFP) ;
      if(reverse_byteorder)
	ug_io_reverse_byte_order(&tmp_quad[0],sizeof(int),4) ;
    }
    quad_faces[qf++] = tmp_quad ;
  }
  
  // Read in boundary flags for surface trias
  // creating buffer to hold boundary information for flag and sending to each pe     
  tag=0;
  ctr=0;  
  autoArray<int> buffer(new int[num_sf_tri_lastpe]) ;
  
  while(ctr < npes-1)  { 
    c=0;
    if(ctr==0) {	
      temp = num_sf_tri_lastpe ;
	
    } else {		
      temp = num_sf_tri ;
    }	
    
    for(int i = 0; i < temp; ++i) {	
      if(!binary)
        fscanf(IFP, "%d", &buffer[i]) ;
      else {
        fread(&buffer[i], sizeof(int), 1, IFP) ;
        if(reverse_byteorder)
          ug_io_reverse_byte_order(&buffer[i],sizeof(int),1) ;
      }
      buffer[c++] = -buffer[i] ;
    }  
    count=c;
    MPI_Recv(&tag, 1, MPI_INT, npes-ctr-1 , 1, MPI_COMM_WORLD, &stat); 

    if(tag==1)
      MPI_Send(&buffer[0], count, MPI_INT, npes-ctr-1, 0, MPI_COMM_WORLD);
    ctr++;
   
  } 
  
  // reading in boundary flag for Processor 0 
  temp = num_sf_tri; 
  for(int i = 0; i < temp; ++i) {	
    if(!binary)
      fscanf(IFP, "%d", &tri_faces[i][3]) ;
    else {
      fread(&tri_faces[i][3], sizeof(int), 1, IFP) ;
      if(reverse_byteorder)
        ug_io_reverse_byte_order(&tri_faces[i][3],sizeof(int),1) ;
    }
    tri_faces[i][3] = -tri_faces[i][3];
  }
      
  // Read in boundary flags for surface quad
  // creating buffer to hold boundary information for flag and sending to each pe     

  autoArray<int> qbuffer(new int[num_sf_quad_lastpe]) ;
  ctr=0; 
  tag=0;
  // Read in boundary flags for surface quads
      
  while(ctr < npes-1)  {
    c=0;
    if(ctr==0) {
      temp = num_sf_quad_lastpe ;
    } else {
      temp = num_sf_quad ;
    }
	
    for(int i = 0; i < temp; ++i) {
      if(!binary)
        fscanf(IFP, "%d", &qbuffer[i]) ;
      else {
        fread(&qbuffer[i], sizeof(int), 1, IFP) ;
        if(reverse_byteorder)
          ug_io_reverse_byte_order(&qbuffer[i],sizeof(int),1) ;      
      }
	  
      qbuffer[c++] = -qbuffer[i] ;
    }
    count = c;

    MPI_Recv(&tag, 1, MPI_INT, npes-ctr-1 , 1, MPI_COMM_WORLD, &stat); 	
    if(tag==1)
      MPI_Send(&qbuffer[0], count, MPI_INT, npes-ctr-1, 0, MPI_COMM_WORLD);
    ctr++;
	
  }
      
  // reading in boundary flag for Processor 0 
  temp = num_sf_quad; 
  for(int i = 0; i < temp; ++i) {
    if(!binary)
      fscanf(IFP, "%d", &quad_faces[i][4]);
    else {
      fread(&quad_faces[i][4], sizeof(int), 1, IFP) ;
      if(reverse_byteorder)
        ug_io_reverse_byte_order(&quad_faces[i][4],sizeof(int),1);
    } 
    quad_faces[i][4] = -quad_faces[i][4] ;
  }
      

  if(myId == 0) 
    cout << "reading volume elements..." << endl ;
   
  // Intializing buffer for sending the tet elements to processors 
  // Read in volume tets and sending to the Processors
  tag=0;
  ctr=0;
  autoArray<tri_info> tri_buffer1(new tri_info[num_tets_lastpe]) ;
      
      
  while(ctr < npes-1) { 
    c=0;
    if(ctr==0) {	
      temp = num_tets_lastpe ;	 
    } else {	
      temp = num_tets ;
    }    	    
	  
    for(int i = 0; i < temp; ++i) {
      Array<int,4> tmp_quad ;
      if(!binary)
        fscanf(IFP, "%d%d%d%d", &tmp_quad[0], &tmp_quad[1], &tmp_quad[2], &tmp_quad[3]) ;
      else {
        fread(&tmp_quad[0], sizeof(int), 4, IFP) ;
        if(reverse_byteorder)
          ug_io_reverse_byte_order(&tmp_quad[0],sizeof(int),4) ;
      }
      tri_buffer1[c++] = tmp_quad ;
    } 
    count = 4*c;
    MPI_Recv(&tag, 1, MPI_INT, npes-ctr-1 , 1, MPI_COMM_WORLD, &stat); 

    if(tag==1)
      MPI_Send(&tri_buffer1[0], count, MPI_INT, npes-ctr-1, 0, MPI_COMM_WORLD);  
    ctr++;
	  
  }   

  /*****************************************************/
  // READING VOLUME ELEMENTS , MAKING FACES & GIVING CELLNUM   
      

  // Calculating total number of vol elmn that would be received by PE 0 
  int cellnum, vol_elmn;

  vol_elmn = num_hex+num_tets+num_pent5+num_pent6;       
    
  if(npes>1) 
    cellnum= num_tets_lastpe + num_tets*(npes-2)+1; 
  else
    cellnum=1;
      
  // Reading in volume tets for processor 0 and converting it into faces  
      
  for(int i = 0; i < num_tets; ++i) {
    int tmp_quad[4] ;
    if(!binary)
      fscanf(IFP, "%d%d%d%d", &tmp_quad[0], &tmp_quad[1], &tmp_quad[2], &tmp_quad[3]) ;
    else {
      fread(&tmp_quad[0], sizeof(int), 4, IFP) ;
      if(reverse_byteorder)
        ug_io_reverse_byte_order(&tmp_quad[0],sizeof(int),4) ;
    }
	
    tri_faces[tf][0] = tmp_quad[0] ;
    tri_faces[tf][1] = tmp_quad[1] ;
    tri_faces[tf][2] = tmp_quad[3] ;
    tri_faces[tf++][3] = cellnum ;
	
    tri_faces[tf][0] = tmp_quad[1] ;
    tri_faces[tf][1] = tmp_quad[2] ;
    tri_faces[tf][2] = tmp_quad[3] ; 
    tri_faces[tf++][3] = cellnum ;
	
    
    tri_faces[tf][0] = tmp_quad[3] ;
    tri_faces[tf][1] = tmp_quad[2] ;
    tri_faces[tf][2] = tmp_quad[0] ;
    tri_faces[tf++][3] = cellnum ;
  
    tri_faces[tf][0] = tmp_quad[0] ;
    tri_faces[tf][1] = tmp_quad[2] ;
    tri_faces[tf][2] = tmp_quad[1] ;
    tri_faces[tf++][3] = cellnum ;
	
    cellnum++ ;
  }
     
  // Intializing buffer for sending the pent5 elements to processors     
  // Read in volume pent5 and sending to the Processors  
  tag=0;
  ctr=0;
  autoArray<quad_info> pent5_buffer(new quad_info[num_pent5_lastpe]) ;
   
  while(ctr < npes-1) {
    c=0;
    if(ctr==0) {	
      temp = num_pent5_lastpe ;
    } else {
      temp = num_pent5 ;
    }		
      
    for(int i = 0; i < temp; ++i) {
      quad_info tmp_pents5 ;
      if(!binary)
        fscanf(IFP, "%d%d%d%d%d", &tmp_pents5[0], &tmp_pents5[1], &tmp_pents5[2], &tmp_pents5[3], &tmp_pents5[4]) ;
      else {
        fread(&tmp_pents5[0], sizeof(int), 5, IFP) ;
        if(reverse_byteorder)
          ug_io_reverse_byte_order(&tmp_pents5[0],sizeof(int),5) ;
      }
      pent5_buffer[c++] = tmp_pents5 ;
    }
      
    count = c*5;
     
    MPI_Recv(&tag, 1, MPI_INT, npes-ctr-1 , 1, MPI_COMM_WORLD, &stat); 
      
    if(tag==1)
      MPI_Send(&pent5_buffer[0], count, MPI_INT, npes-ctr-1, 0, MPI_COMM_WORLD);
      
    ctr++;
  } 

  // calculating cellnum for pents5 elements
  if(npes>1)
    cellnum = num_vol_tets +  num_pent5_lastpe + num_pent5*(npes-2)+1; 
  
  // read in volume pent 5's and convert to faces for pe 0


  for(int i = 0; i < num_pent5; ++i) {
    int tmp_pents5[5] ;
    if(!binary)
      fscanf(IFP, "%d%d%d%d%d", &tmp_pents5[0], &tmp_pents5[1], &tmp_pents5[2], &tmp_pents5[3], &tmp_pents5[4]) ;
    else {
      fread(&tmp_pents5[0], sizeof(int), 5, IFP) ;
      if(reverse_byteorder)
        ug_io_reverse_byte_order(&tmp_pents5[0],sizeof(int),5) ;
    }
    tri_faces[tf][0] = tmp_pents5[4] ;
    tri_faces[tf][1] = tmp_pents5[1] ;
    tri_faces[tf][2] = tmp_pents5[2] ;
    tri_faces[tf++][3] = cellnum ;
    
    tri_faces[tf][0] = tmp_pents5[4] ;
    tri_faces[tf][1] = tmp_pents5[2] ;
    tri_faces[tf][2] = tmp_pents5[3] ;
    tri_faces[tf++][3] = cellnum ;
    
    tri_faces[tf][0] = tmp_pents5[3] ;
    tri_faces[tf][1] = tmp_pents5[2] ;
    tri_faces[tf][2] = tmp_pents5[0] ;
    tri_faces[tf++][3] = cellnum ;
    
    tri_faces[tf][0] = tmp_pents5[0] ;
    tri_faces[tf][1] = tmp_pents5[2] ;
    tri_faces[tf][2] = tmp_pents5[1] ;
    tri_faces[tf++][3] = cellnum ;
    
    quad_faces[qf][0] = tmp_pents5[0] ;
    quad_faces[qf][1] = tmp_pents5[1] ;
    quad_faces[qf][2] = tmp_pents5[4] ;
    quad_faces[qf][3] = tmp_pents5[3] ;
    quad_faces[qf++][4] = cellnum ;
    cellnum++ ;
  } 
  
  // Intializing buffer for sending the pent6 elements to processors
  // Read in volume pent6 and sending to the Processors  
  
  Array<int,6> tmp_pents6 ;
  autoArray<Array<int,6> > pent6_buffer(new Array<int,6>[num_pent6_lastpe]) ;
  tag=0;
  ctr=0;
  while(ctr < npes-1) {
    c=0;
    if(ctr==0) {		
      temp = num_pent6_lastpe ;
    } else {
      temp = num_pent6 ;
    }
      
    for(int i = 0; i < temp; ++i) {
      if(!binary)
        fscanf(IFP, "%d%d%d%d%d%d", &tmp_pents6[0], &tmp_pents6[1], &tmp_pents6[2], &tmp_pents6[3], &tmp_pents6[4], &tmp_pents6[5]) ;
      else {
        fread(&tmp_pents6[0], sizeof(int), 6, IFP) ;
        if(reverse_byteorder)
          ug_io_reverse_byte_order(&tmp_pents6[0],sizeof(int),6) ;
      }
	
      pent6_buffer[c++]=tmp_pents6; 
    }  
    count = 6*c;
 
    MPI_Recv(&tag, 1, MPI_INT, npes-ctr-1 , 1, MPI_COMM_WORLD, &stat); 
      
    if(tag==1)
      MPI_Send(&pent6_buffer[0], count, MPI_INT, npes-ctr-1, 0, MPI_COMM_WORLD);
      
    ctr++;
  }  
  
  // Calculating cellnum for pent6 elements
  if(npes>1)
    cellnum = num_vol_pents5 +  num_pent6_lastpe + num_pent6*(npes-2)+1; 
  

  for(int i = 0; i < num_pent6; ++i) {
    
    if(!binary)
      fscanf(IFP, "%d%d%d%d%d%d", &tmp_pents6[0], &tmp_pents6[1], &tmp_pents6[2], &tmp_pents6[3], &tmp_pents6[4], &tmp_pents6[5]) ;
    else {
      fread(&tmp_pents6[0], sizeof(int), 6, IFP) ;
      if(reverse_byteorder)
	ug_io_reverse_byte_order(&tmp_pents6[0],sizeof(int),6) ;
    } 
  
    
    tri_faces[tf][0] = tmp_pents6[3] ;
    tri_faces[tf][1] = tmp_pents6[4] ;
    tri_faces[tf][2] = tmp_pents6[5] ;
    tri_faces[tf++][3] = cellnum ;
    
    tri_faces[tf][0] = tmp_pents6[0] ;
    tri_faces[tf][1] = tmp_pents6[2] ;
    tri_faces[tf][2] = tmp_pents6[1] ;
    tri_faces[tf++][3] = cellnum ;
    
    quad_faces[qf][0] = tmp_pents6[0] ;
    quad_faces[qf][1] = tmp_pents6[1] ;
    quad_faces[qf][2] = tmp_pents6[4] ;
    quad_faces[qf][3] = tmp_pents6[3] ;
    quad_faces[qf++][4] = cellnum ;
    
    quad_faces[qf][0] = tmp_pents6[1] ;
    quad_faces[qf][1] = tmp_pents6[2] ;
    quad_faces[qf][2] = tmp_pents6[5] ;
    quad_faces[qf][3] = tmp_pents6[4] ;
    quad_faces[qf++][4] = cellnum ;
    
    quad_faces[qf][0] = tmp_pents6[3] ;
    quad_faces[qf][1] = tmp_pents6[5] ;
    quad_faces[qf][2] = tmp_pents6[2] ;
    quad_faces[qf][3] = tmp_pents6[0] ;
    quad_faces[qf++][4] = cellnum ;
    
    cellnum++ ;
  }
 
  // Intializing buffer for sending the Hex elements to processors
  // Read in volume hex and sending to the Processors
  tag=0;
  ctr=0;
  autoArray<Array<int,8> > hex_buffer(new Array<int,8>[num_hex_lastpe]) ;
  
  while(ctr < npes-1) {
    c=0;
    if(ctr==0) {	
      temp = num_hex_lastpe ;
    } else {
      temp = num_hex ;
    }	
	
    for(int i = 0; i < temp; ++i) {
      Array<int,8> tmp_hexs ;
      if(!binary)
        fscanf(IFP, "%d%d%d%d%d%d%d%d", &tmp_hexs[0], &tmp_hexs[1], &tmp_hexs[2], &tmp_hexs[3], &tmp_hexs[4], &tmp_hexs[5], &tmp_hexs[6], &tmp_hexs[7] ) ;
      else {
        fread(&tmp_hexs[0], sizeof(int), 8, IFP) ;
        if(reverse_byteorder)
          ug_io_reverse_byte_order(&tmp_hexs[0],sizeof(int),8) ;
      }
      hex_buffer[c++]=tmp_hexs;  
    }

    count = 8*c;
    
    MPI_Recv(&tag, 1, MPI_INT, npes-ctr-1 , 1, MPI_COMM_WORLD, &stat); 

    if(tag==1)
      MPI_Send(&hex_buffer[0], count, MPI_INT, npes-ctr-1, 0, MPI_COMM_WORLD);  
    ctr++;
  }  

  // calculating cellnum for the hex elements  
  if(npes>1)
    cellnum = num_vol_pents6 +  num_hex_lastpe + num_hex*(npes-2)+1; 

  // Reading in volume Hex element to Processor 0 and converting to faces

  for(int i = 0; i < num_hex; ++i) {
    int tmp_hexs[8] ;
    if(!binary)
      fscanf(IFP, "%d%d%d%d%d%d%d%d", &tmp_hexs[0], &tmp_hexs[1], &tmp_hexs[2], &tmp_hexs[3], &tmp_hexs[4], &tmp_hexs[5], &tmp_hexs[6], &tmp_hexs[7] ) ;
    else {
      fread(&tmp_hexs[0], sizeof(int), 8, IFP) ;
      if(reverse_byteorder)
	ug_io_reverse_byte_order(&tmp_hexs[0],sizeof(int),8) ;
    }
    
    quad_faces[qf][0] = tmp_hexs[0] ;
    quad_faces[qf][1] = tmp_hexs[1] ;
    quad_faces[qf][2] = tmp_hexs[5] ;
    quad_faces[qf][3] = tmp_hexs[4] ;
    quad_faces[qf++][4] = cellnum ;
    
    quad_faces[qf][0] = tmp_hexs[1] ;
    quad_faces[qf][1] = tmp_hexs[2] ;
    quad_faces[qf][2] = tmp_hexs[6] ;
    quad_faces[qf][3] = tmp_hexs[5] ;
    quad_faces[qf++][4] = cellnum ;
    
    quad_faces[qf][0] = tmp_hexs[2] ;
    quad_faces[qf][1] = tmp_hexs[3] ;
    quad_faces[qf][2] = tmp_hexs[7] ;
    quad_faces[qf][3] = tmp_hexs[6] ;
    quad_faces[qf++][4] = cellnum ;
    
    quad_faces[qf][0] = tmp_hexs[4] ;
    quad_faces[qf][1] = tmp_hexs[7] ;
    quad_faces[qf][2] = tmp_hexs[3] ;
    quad_faces[qf][3] = tmp_hexs[0] ;
    quad_faces[qf++][4] = cellnum ;
    
    quad_faces[qf][0] = tmp_hexs[4] ;
    quad_faces[qf][1] = tmp_hexs[5] ;
    quad_faces[qf][2] = tmp_hexs[6] ;
    quad_faces[qf][3] = tmp_hexs[7] ;
    quad_faces[qf++][4] = cellnum ;
    
    quad_faces[qf][0] = tmp_hexs[3] ;
    quad_faces[qf][1] = tmp_hexs[2] ;
    quad_faces[qf][2] = tmp_hexs[1] ;
    quad_faces[qf][3] = tmp_hexs[0] ;
    quad_faces[qf++][4] = cellnum ;
    
    cellnum++ ;
  }


  if(myId == 0)
    cout << "finished reading ugrid file, matching faces..." << endl ;
  fclose(IFP) ;

  

  // prepare triangle faces (sort them)  
  for(size_t i=0;i<tf;++i) {
    // xdr numbers nodes from zero
    tri_faces[i][0] -= 1 ;
    tri_faces[i][1] -= 1 ;
    tri_faces[i][2] -= 1 ;
    
    if(tri_faces[i][0] > tri_faces[i][1])
      std::swap(tri_faces[i][0],tri_faces[i][1]) ;
    if(tri_faces[i][0] > tri_faces[i][2])
      std::swap(tri_faces[i][0],tri_faces[i][2]) ;
    if(tri_faces[i][1] > tri_faces[i][2])
      std::swap(tri_faces[i][1],tri_faces[i][2]) ;
  }

  // prepare quad faces (sort them, but be careful)
  for(size_t i=0;i<qf;++i) {
    // xdr numbers nodes from zero
    quad_faces[i][0] -=1 ;
    quad_faces[i][1] -=1 ;
    quad_faces[i][2] -=1 ;
    quad_faces[i][3] -=1 ;
    // First make sure first entry is lowest number
    int tmp_face[4] ;
    int vs = quad_faces[i][0] ;
    size_t nv = 0 ;
    for(size_t j=1;j<4;++j)
      if(vs > quad_faces[i][j]) {
        vs = quad_faces[i][j] ;
        nv = j ;
      }
    for(size_t j=0;j<4;++j)
      tmp_face[j] = quad_faces[i][(j+nv)&0x3] ;
    // next make orientation so that it will match other face 
    if(tmp_face[1] < tmp_face[3])
      for(int j=0;j<4;++j)
        quad_faces[i][j] = tmp_face[j] ;
    else
      for(size_t j=0;j<4;++j)
        quad_faces[i][j] = tmp_face[(4 - j) &0x3 ] ;
  }


  if(myId == 0) 
    cout << "sorting faces..." << endl ;

  /*******************************************************/
  //       STARTS BUCKET SORT FOR TRIAS
  

  // Local sort 
  std::sort(&tri_faces[0],&tri_faces[tf],tri_sort);     
  
  num = Info[0]/npes ;
  int rem = Info[0] % npes ;
  int cnt=0;  
   
  autoArray<int> count1(new int[npes*2]) ;
  for(int j=0; j<npes*2; j++)
    {
      count1[j]=0;
    }     
   
  // Count array determines which node belongs to which PE'S
  for(int i = 0; i<npes; i++ )	{	
    if(i==npes-1 && rem>0)
      {		
	count1[cnt]=i*num;			
	cnt++;		
	count1[cnt]=Info[0];			
      }
    else										      {	
      count1[cnt]=i*num;			
      cnt++;
      count1[cnt]=(i+1)*num;
      cnt++;
    }	
  }	
   
  // Determines the number of elements that would be broadcasted
  autoArray<int> splitter(new int[npes]) ;
  autoArray<int> recv_num_elmn(new int[npes]);
  autoArray<int> sdisp_tri(new int[npes]) ;
  autoArray<int> rdisp_tri(new int[npes]);
   
  for(int j=0; j<npes; j++) {
    splitter[j] = 0;
  }  
   
  for(int i=0; i < (int)tf ; i++) {
    for(int j=0; j< npes; j++) {   
	   
      if((tri_faces[i][0]>count1[j*2] || tri_faces[i][0]==count1[j*2])&& tri_faces[i][0]< count1[j*2+1]) {
        splitter[j]++; 
	       
      } 
    }             
  }    
  // Calculating displacement for sending the elements of the array
   
  for(int i=0; i<npes ; i++)
    splitter[i]=splitter[i]*4;
   
  sdisp_tri[0]=0;
  for(int i=1; i<npes; i++) {
    sdisp_tri[i] = splitter[i-1]+sdisp_tri[i-1];
  } 
  
  
  // Using All to All personalized communication 
  int sum_tri=0;
  
  MPI_Alltoall(&splitter[0], 1, MPI_INT, &recv_num_elmn[0], 1, MPI_INT, MPI_COMM_WORLD);  
  
  // calculating the size of the receiving array
  
  for(int i=0;i<npes;i++) 
    sum_tri = sum_tri + recv_num_elmn[i]/4;  
  
 
  // Calculating the displacement for storing the array
 
  rdisp_tri[0] = 0;
  for(int i=1; i<npes; i++)  {
    rdisp_tri[i] = recv_num_elmn[i-1] + rdisp_tri[i-1] ;
  }	
 
  // initialize array  
  autoArray<Array<int,4> > new_tri_faces(new Array<int,4>[sum_tri]) ;
  
  MPI_Alltoallv(&tri_faces[0], &splitter[0], &sdisp_tri[0], MPI_INT, &new_tri_faces[0], &recv_num_elmn[0] , &rdisp_tri[0] , MPI_INT, MPI_COMM_WORLD);  

  // final LOCAL sort  
  std::sort(&new_tri_faces[0],&new_tri_faces[sum_tri],tri_sort);        



  /********************************************************/
  //          BUCKET SORT FOR QUADS 
    
  std::sort(&quad_faces[0], &quad_faces[qf], quad_sort);  
    
  // similar procedure for quad array   
  // Determines the number of elements that would be broadcasted
  autoArray<int> splitter_quad(new int[npes]); 
  autoArray<int> sdisp_quad(new int[npes]);
  autoArray<int> rdisp_quad(new int[npes]);
  autoArray<int> recv_elmn_quad(new int[npes]);

  for(int j=0; j<npes; j++) {
    splitter_quad[j] = 0;
  }  
  for(int i=0; i < (int)qf ; i++) {	
    for(int j=0; j< npes; j++) { 
	    
      if((quad_faces[i][0]>count1[j*2] || quad_faces[i][0]==count1[j*2])&& quad_faces[i][0]< count1[j*2+1]) {
        splitter_quad[j]++;
		
      } 
    }             
  }   
  // Calculating displacement for sending the elements of the array

  for(int i=0;i<npes;i++)
    splitter_quad[i]=splitter_quad[i]*5;    
    

  sdisp_quad[0]=0;
  for(int i=1; i<npes; i++) {
    sdisp_quad[i] = splitter_quad[i-1]+sdisp_quad[i-1];
  }     
    
  // Using All to All personalized communication 
    
    
  MPI_Alltoall(&splitter_quad[0], 1, MPI_INT, &recv_elmn_quad[0], 1, MPI_INT, MPI_COMM_WORLD); 
   
    
  // calculating the size of the receiving array
  int sum_quad=0;
    
  for(int i=0;i<npes;i++) 
    sum_quad = sum_quad + recv_elmn_quad[i]/5;
  
   
  // Calculating the displacement for storing in the array
  rdisp_quad[0] = 0;
  for(int i=1; i<npes; i++)  {
    rdisp_quad[i] = recv_elmn_quad[i-1] + rdisp_quad[i-1] ;
  }	 
   
  // initialize array
  autoArray<Array<int,5> > new_quad_faces(new Array<int,5>[sum_quad]) ;
    
  MPI_Alltoallv(&quad_faces[0], &splitter_quad[0], &sdisp_quad[0], MPI_INT, &new_quad_faces[0], &recv_elmn_quad[0], &rdisp_quad[0], MPI_INT, MPI_COMM_WORLD);  
  // final sort 
  std::sort(&new_quad_faces[0], &new_quad_faces[sum_quad], quad_sort);      
     
    

  if(myId == 0)
    cout << "writing face information..." << endl ;
    
  /*********************************************************/
  // WRITING FACE INFORMATION TO THE XDR FILE 

          
  size_t btf = 0 ;
  size_t bqf = 0 ;
  
  int off = 0 ;

  for(int i=0; i<sum_tri/2 ; ++i) {
    if(new_tri_faces[i*2][0] != new_tri_faces[i*2+1][0] ||
       new_tri_faces[i*2][1] != new_tri_faces[i*2+1][1] ||
       new_tri_faces[i*2][2] != new_tri_faces[i*2+1][2] ||
       (new_tri_faces[i*2][3] < 0 && new_tri_faces[i*2+1][3] < 0)) {
      cerr << "trouble matching triangle faces! " << endl ;
      cerr << "perhaps an embedded surface remains in the grid?" << endl ;
      exit(-1) ;
    }
	
    if(new_tri_faces[i*2][3] < 0) 
      btf++;  
  }
      
      
  // MAKING BUFFER FOR BOUDARY TRIAS 

  autoArray<Array<int,5> > bnd_tri(new Array<int,5>[btf]) ;
      
  btf=0;
 
  for(int i=0; i<sum_tri/2 ; ++i) {
    if(new_tri_faces[i*2][3] < 0) {
      bnd_tri[btf][0] = new_tri_faces[i*2][0] ;
      bnd_tri[btf][1] = new_tri_faces[i*2][1] ;
      bnd_tri[btf][2] = new_tri_faces[i*2][2] ;
      bnd_tri[btf][3] = new_tri_faces[i*2+1][3] ;
      bnd_tri[btf++][4] = new_tri_faces[i*2][3] ;
    } else {
      xdr_int(&xdr_handle, &off) ;
      xdr_int(&xdr_handle, &new_tri_faces[i*2+1][3]) ;
      xdr_int(&xdr_handle, &new_tri_faces[i*2][3]) ;
	  
      off += 3 ;
    }  
  }
      
      
  // for recv the sorted internal tri cellnum frm the pe's   
            
  ctr = 1;
  int dum = 0, flag = 1;
              

  while(ctr < npes) {
	  
    MPI_Recv(&dum, 1, MPI_INT, ctr, 1, MPI_COMM_WORLD, &stat);
	  	  
    autoArray<Array<int,2> > buff_tri_faces(new Array<int,2>[dum/2]) ;
	  
    MPI_Send(&flag, 1, MPI_INT, ctr, 0, MPI_COMM_WORLD);

    MPI_Recv( &buff_tri_faces[0], dum, MPI_INT, ctr, 1, MPI_COMM_WORLD, &stat);
    
    for(int i=0; i<dum/2 ; ++i) {
      xdr_int(&xdr_handle, &off) ;
      xdr_int(&xdr_handle, &buff_tri_faces[i][0]);
      xdr_int(&xdr_handle, &buff_tri_faces[i][1]);
	    
      off += 3 ;
    }
	  
    ctr++;
  }
      
      

  for(int i=0;i<sum_quad/2;++i) {
    if(new_quad_faces[i*2][0] != new_quad_faces[i*2+1][0] ||
       new_quad_faces[i*2][1] != new_quad_faces[i*2+1][1] ||
       new_quad_faces[i*2][2] != new_quad_faces[i*2+1][2] ||
       new_quad_faces[i*2][3] != new_quad_faces[i*2+1][3] ||
       (new_quad_faces[i*2][4] < 0 && new_quad_faces[i*2+1][4] < 0)) {
      cerr << "trouble matching quad faces!" << endl ;
      cerr << "perhaps an embedded surface remains in the grid?" << endl ;
      exit(-1) ;
    }
	
    if(new_quad_faces[i*2][4] < 0) 
      bqf++;
	
  }
      
      
  // MAKING BUFFER FOR BOUNDARY QUADS

  autoArray<Array<int,6> > bnd_quad(new Array<int,6>[bqf]) ;
      
  bqf=0;
      
  for(int i=0;i<sum_quad/2;++i) {
    if(new_quad_faces[i*2][4] < 0) {
      bnd_quad[bqf][0] = new_quad_faces[i*2][0] ;
      bnd_quad[bqf][1] = new_quad_faces[i*2][1] ;
      bnd_quad[bqf][2] = new_quad_faces[i*2][2] ;
      bnd_quad[bqf][3] = new_quad_faces[i*2][3] ;
      bnd_quad[bqf][4] = new_quad_faces[i*2+1][4] ;
      bnd_quad[bqf++][5] = new_quad_faces[i*2][4] ;
    } else {
      xdr_int(&xdr_handle, &off) ;
      xdr_int(&xdr_handle, &new_quad_faces[i*2+1][4]) ;
      xdr_int(&xdr_handle, &new_quad_faces[i*2][4]) ;

      off += 4 ;
    }
  }
      
      
      
  // for recv the sorted internal quad faces frm the pe's
      
  ctr = 0;
  int sum = 0;
  while(ctr < npes-1) {	
	  
    MPI_Recv( &sum, 1, MPI_INT, ctr+1, 1, MPI_COMM_WORLD, &stat);
    
    autoArray<Array<int,2> > buff_quad_faces(new Array<int,2>[sum/2]) ;

    // Sending the acknowledgement to send the buffer
    MPI_Send(&flag, 1, MPI_INT, ctr+1, 0, MPI_COMM_WORLD);
    
    MPI_Recv( &buff_quad_faces[0], sum,  MPI_INT, ctr+1, 1, MPI_COMM_WORLD, &stat);
    
    for(int i=0; i<sum/2 ; ++i) {
      xdr_int(&xdr_handle, &off) ;
      xdr_int(&xdr_handle, &buff_quad_faces[i][0]) ;
      xdr_int(&xdr_handle, &buff_quad_faces[i][1]) ;
	    
      off += 4 ;
    }
	  
    ctr++;
  }
           

  // Getting the number of boundary trias from each pe 

  ctr = 0;
  sum = 0;
  autoArray<int> b_num(new int[npes-1]);
  int a;

  size_t total = btf*5;
  autoArray<int> b_total(new int[npes]);
      
  b_total[0] = btf*5;
  while(ctr < npes-1) {
    MPI_Recv( &a, 1, MPI_INT, ctr+1, 3, MPI_COMM_WORLD, &stat);
	  
    b_num[ctr]=a;
	  
    b_total[ctr+1] = b_total[ctr] + b_num[ctr];
       
    ctr++; 
  }
      
      
  // total number of boundary trias
  int t_btf = b_total[npes-1]/5;
   
  autoArray<Array<int,5> > buff_btri(new Array<int,5>[t_btf]) ;


  for(size_t i=0 ; i<btf ;i++) {
    buff_btri[i][0]= bnd_tri[i][0] ;
    buff_btri[i][1] = bnd_tri[i][1];
    buff_btri[i][2] = bnd_tri[i][2] ;
    buff_btri[i][3] = bnd_tri[i][3]; 
    buff_btri[i][4] = bnd_tri[i][4]; 
  }
   
  // Getting the boundary trias from each pe

  ctr = 0;
  sum = 0;
  int t = b_total[ctr]/5;      
      
  while(ctr < npes-1) {

    autoArray<Array<int,5> > temp_bnd_tri(new Array<int,5>[b_num[ctr]/5]) ;
       
    // Sending the acknowledgement to send the buffer
    MPI_Send(&flag, 1, MPI_INT, ctr+1, 0, MPI_COMM_WORLD);
    

    MPI_Recv( &temp_bnd_tri[0], b_num[ctr], MPI_INT, ctr+1, 1, MPI_COMM_WORLD, &stat);
       
    for(int i=0;i<b_num[ctr]/5;i++)
      {
        buff_btri[t][0]= temp_bnd_tri[i][0];
        buff_btri[t][1]= temp_bnd_tri[i][1];
        buff_btri[t][2]= temp_bnd_tri[i][2];
        buff_btri[t][3]= temp_bnd_tri[i][3];
        buff_btri[t++][4]= temp_bnd_tri[i][4];
      }
  
    ctr++;
  
  }
      
     
  std::sort(&buff_btri[0], &buff_btri[t_btf], bnd_tri_sort) ;

  //    WRITING THE BND FACES 

  for(int i=0;i<t_btf;++i) {
    xdr_int(&xdr_handle, &off) ;
    xdr_int(&xdr_handle, &buff_btri[i][3]) ;
    xdr_int(&xdr_handle, &buff_btri[i][4]) ;

    off += 3 ;
 
  }
 

  // Getting the number of boundary quads from each pe 

  ctr = 0;
  sum = 0;
      
  total = bqf*6;
      
  b_total[0] = bqf*6;
      
  while(ctr < npes-1) {
       
    MPI_Recv( &a, 1, MPI_INT, ctr+1, 1, MPI_COMM_WORLD, &stat);

    b_num[ctr] = a;
    b_total[ctr+1] = b_total[ctr] + b_num[ctr];

    ctr++; 
  }


//  cout<<"after recieving numbers "<<endl;

  // total number of boundary quads
  int q_bqf = b_total[npes-1]/6;

  autoArray<Array<int,6> > buff_bquad(new Array<int,6>[q_bqf]) ;

  for(size_t i=0 ; i<bqf ;i++) {
    buff_bquad[i][0]= bnd_quad[i][0] ;
    buff_bquad[i][1] = bnd_quad[i][1];
    buff_bquad[i][2] = bnd_quad[i][2] ;
    buff_bquad[i][3] = bnd_quad[i][3] ;
    buff_bquad[i][4] = bnd_quad[i][4];
    buff_bquad[i][5] = bnd_quad[i][5];
  }


  // Getting the boundary trias from each pe

  ctr = 0;
  sum = 0;
  t = b_total[0]/6;
      
  while(ctr < npes-1) {
    autoArray<Array<int,6> > temp_bnd_quad(new Array<int,6>[b_num[ctr]/6]) ;

    // Sending the acknowledgement to send the buffer
    MPI_Send(&flag, 1, MPI_INT, ctr+1, 0, MPI_COMM_WORLD);

    MPI_Recv( &temp_bnd_quad[0], b_num[ctr], MPI_INT, ctr+1, 1, MPI_COMM_WORLD, &stat);

    for(int i=0;i<b_num[ctr]/6;i++) {
      buff_bquad[t][0]=temp_bnd_quad[i][0];
      buff_bquad[t][1]=temp_bnd_quad[i][1];
      buff_bquad[t][2]=temp_bnd_quad[i][2];
      buff_bquad[t][3]=temp_bnd_quad[i][3];
      buff_bquad[t][4]=temp_bnd_quad[i][4];
      buff_bquad[t++][5]=temp_bnd_quad[i][5];
    }
	  
    ctr++;
  }


     

  std::sort(&buff_bquad[0], &buff_bquad[q_bqf], bnd_quad_sort) ;

  // Writing the bnd quad faces 

  for(int i=0;i<q_bqf;++i) {
    xdr_int(&xdr_handle, &off) ;
    xdr_int(&xdr_handle, &buff_bquad[i][4]) ;
    xdr_int(&xdr_handle, &buff_bquad[i][5]) ;
    
    off += 4 ;
  }
  
  xdr_int(&xdr_handle, &off) ;



   
  for(int i=0; i<sum_tri/2; ++i) {
    if(new_tri_faces[i*2][3] >=0) {
      xdr_int(&xdr_handle, &new_tri_faces[i*2][0]) ;
      xdr_int(&xdr_handle, &new_tri_faces[i*2][1]) ;
      xdr_int(&xdr_handle, &new_tri_faces[i*2][2]) ;
    }
  }
  
  
  // for recv the sorted internal tri nodes frm the pe's     
  
  ctr=0;
  sum=0;
  while(ctr < npes-1) {


    MPI_Recv( &sum, 1, MPI_INT, ctr+1, 1, MPI_COMM_WORLD, &stat);
       
    autoArray<Array<int,3> > buff_tri_faces(new Array<int,3>[sum/3]) ;
    // Sending the acknowledgement to send the buffer
    MPI_Send(&flag, 1, MPI_INT, ctr+1, 0, MPI_COMM_WORLD);
  
    MPI_Recv( &buff_tri_faces[0], sum, MPI_INT, ctr+1, 1, MPI_COMM_WORLD, &stat);
 
    for(int i=0; i<sum/3 ; ++i) {
      xdr_int(&xdr_handle, &buff_tri_faces[i][0]);
      xdr_int(&xdr_handle, &buff_tri_faces[i][1]);
      xdr_int(&xdr_handle, &buff_tri_faces[i][2]);
	 
    }
       
    ctr++;
  }
  
  // writng the  quad nodes to the file

  for(int i=0;i<sum_quad/2;++i) {
    if(new_quad_faces[i*2][4] >=0) {
      xdr_int(&xdr_handle, &new_quad_faces[i*2][0]) ;
      xdr_int(&xdr_handle, &new_quad_faces[i*2][1]) ;
      xdr_int(&xdr_handle, &new_quad_faces[i*2][2]) ;
      xdr_int(&xdr_handle, &new_quad_faces[i*2][3]) ;
      
    }
  }
  
  
  //for recv the sorted internal quad nodes frm the pe's     
  
  ctr=0;
  sum=0;
  
  while(ctr < npes-1) {

    //    if(myId == 0)
    //      cout<<endl<<"buff_quad_nodes  from "<<ctr<<endl;

    MPI_Recv( &sum, 1, MPI_INT, ctr+1, 1, MPI_COMM_WORLD, &stat);
      
    autoArray<Array<int,4> > buff_quad_faces(new Array<int,4>[sum/4]) ; 
    // Sending the acknowledgement to send the buffer
    MPI_Send(&flag, 1, MPI_INT, ctr+1, 0, MPI_COMM_WORLD);

         
    MPI_Recv( &buff_quad_faces[0], sum, MPI_INT, ctr+1, 1, MPI_COMM_WORLD, &stat);
    
    for(int i=0; i<sum/4 ; ++i) {
      xdr_int(&xdr_handle, &buff_quad_faces[i][0]) ;
      xdr_int(&xdr_handle, &buff_quad_faces[i][1]) ;
      xdr_int(&xdr_handle, &buff_quad_faces[i][2]) ;
      xdr_int(&xdr_handle, &buff_quad_faces[i][3]) ;

    }
      
    ctr++;
  }
  
  // writing the bnd nodes to the file

 
  for(int i=0;i<t_btf;++i) {
    xdr_int(&xdr_handle, &buff_btri[i][0]) ;
    xdr_int(&xdr_handle, &buff_btri[i][1]) ;
    xdr_int(&xdr_handle, &buff_btri[i][2]) ;
  }
  
  
  
  
  // writing the bnd quad nodes info
  
  
  for(int i=0; i<q_bqf; ++i) {
    xdr_int(&xdr_handle, &buff_bquad[i][0]) ;
    xdr_int(&xdr_handle, &buff_bquad[i][1]) ;
    xdr_int(&xdr_handle, &buff_bquad[i][2]) ;
    xdr_int(&xdr_handle, &buff_bquad[i][3]) ;

  }
  
  
  xdr_destroy(&xdr_handle) ;
  fclose(FP) ;
  return 0 ;  
}



/*************************************************************************/
// Start of client
void Client(int id, int numProcessors) {
  
  int myrank, npes, vol_elmn_lastpe; 
  MPI_Status stat;
    
  // Info array contains the grid information like num of nodes, trias etc..
  int Info[7] ;
  int info_cellnum[4] ;
  npes = numProcessors;
  myrank = id;

 
  MPI_Bcast( Info, 7, MPI_INT, 0, MPI_COMM_WORLD);  


  MPI_Bcast( info_cellnum, 4, MPI_INT, 0, MPI_COMM_WORLD);  

  MPI_Bcast(&vol_elmn_lastpe, 1, MPI_INT, 0, MPI_COMM_WORLD);  

  /***********************************************************/
  // CALCULATING THE NUMBER OF ELEMENTS EACH PE WILL RECIEVE

  // Calculating number of boundary trias for each processor and allocating mem 
 
  int num = Info[1]/npes ; 
  int num_sf_tri;   
  if(num % 2 == 1)  {
    num_sf_tri = num - 1;  
    if(myrank==npes-1)
      num_sf_tri = Info[1] - (num-1)*(npes-1);
  } else {
    num_sf_tri = num ;
    if(Info[1] % npes != 0 && myrank==npes-1)		
      num_sf_tri = Info[1] - num*(npes-1);
  }		   
		      
  // Calculating number of boundary quads for each processor and allocating mem    
  num = Info[2]/npes ;
  int num_sf_quad;
  if(num % 2 == 1) {	
    num_sf_quad = num - 1; 
    if(myrank==npes-1)
      num_sf_quad = Info[2] - (num - 1)*(npes-1);
  } else {    
    num_sf_quad = num ;
    if(Info[2] % npes != 0 && myrank == npes-1)	 
      num_sf_quad = Info[2] - num*(npes-1);
  }

      
  // Calculating number of tets to give to each processor
  num = Info[3]/npes; 
  int num_tets; 
  if(Info[3] % npes != 0)  {
    num_tets = num ;
    if(myrank == npes-1)
      num_tets = Info[3] - num*(npes-1);
  } else {
    num_tets = num;
    if(myrank == npes-1) 
      num_tets = num;
  } 

  //Calculating number of pent5 to give to each processor for trifaces
  // This should be allocated evenly to each processor  

  num = Info[4]/npes; 
  int num_pent5;   
  if(num % 2 == 1)  {
    num_pent5 = num - 1; 
    if(myrank==npes-1)
      num_pent5 = Info[4] - (num - 1)*(npes-1);
  } else { 
    num_pent5 = num ;
    if(Info[4] % npes != 0 && myrank==npes-1)
      num_pent5 = Info[4] - num*(npes-1);
  }    
  
 
  // Calculating number of Pent6 to give to each processor
  num = Info[5]/npes;  
  int num_pent6 ;
  
  if(num % 2 == 1)  {
    num_pent6 = num - 1 ;
    if(myrank==npes-1)
      num_pent6 = Info[5] - (num-1)*(npes-1);
  } else {
    num_pent6 = num ; 
    if(Info[5] % npes != 0 && myrank == npes-1) 
      num_pent6 = Info[5] - num*(npes-1);
  }
  
  // Calculating the number of Hex to give to each processor
  num = Info[6]/npes;
  int num_hex;

  if(Info[6] % npes != 0)  {
    num_hex = num ;
    if(myrank==npes-1) 	
      num_hex = Info[6] - num*(npes-1);
  } else { 
    num_hex = num;
  } 
 
  // Calculating the total number of Tri_faces & quad faces generated at each PE // Initializing mem for each

  size_t tri_total = num_sf_tri + num_tets*4 + num_pent5*4 + num_pent6*2 ;
  size_t quad_total = num_sf_quad + num_pent5 + num_pent6*3 + num_hex*6;

  autoArray<tri_info> tri_faces(new tri_info[tri_total]) ;
  autoArray<quad_info> quad_faces(new quad_info[quad_total]) ;
  
  size_t tf = 0 ;// triangle and quad face pointers
  size_t qf = 0 ;
  

  /**************************************************************/
  //         RECEIVING ELEMENTS FROM SOURCE PE 

 
  // receiving the set of sf_trias
  int tp = num_sf_tri*4,tag=1; 

  MPI_Send(&tag, 1, MPI_INT, 0, 1, MPI_COMM_WORLD); 

  MPI_Recv(&tri_faces[0],tp, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat);
    
  MPI_Send(&tag, 1, MPI_INT, 0, 1, MPI_COMM_WORLD); 
    
    
  // receiving the set of sf_quads
  MPI_Recv(&quad_faces[0], (num_sf_quad*5), MPI_INT, 0, 0, MPI_COMM_WORLD, &stat);
    
    
    
  // Creating buffer to recieve the boundary flags for sf_trias
  autoArray<int> tri_bf(new int[num_sf_tri]);
    
  MPI_Send(&tag, 1, MPI_INT, 0, 1, MPI_COMM_WORLD); 
    
  MPI_Recv(&tri_bf[0], num_sf_tri, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat);
    
    
  for(int i=0; i<num_sf_tri; i++) {
    tri_faces[i][3] = tri_bf[i];
  }
    
  // Creating buffer to recieve the boundary flags for sf_quads
    
  autoArray<int> quad_bf(new int[num_sf_quad]);
    
  MPI_Send(&tag, 1, MPI_INT, 0, 1, MPI_COMM_WORLD); 
    
  MPI_Recv(&quad_bf[0], num_sf_quad, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat);
    
  for(int i=0; i<num_sf_quad; i++) {
    quad_faces[i][4] = quad_bf[i];
  }
    
    
  // Initializing cell number according to the number of PE
    
  int cellnum ;
    
    
  if(myrank==npes-1)
    cellnum = 1;
  else 
    cellnum = info_cellnum[0]+ num_tets*(npes - myrank -2)+1; 
  
    
  /***************************************************/
  // READING IN VOLUME ELEMENTS AND CONVERTING THEM INTO FACES
    
    
  // Receive volume tets from PE 0
  // creating buffer for the receiving tet elements

  autoArray<fourtuple> tmp_tet(new fourtuple[num_tets]) ;
  
  MPI_Send(&tag, 1, MPI_INT, 0, 1, MPI_COMM_WORLD); 
  
  MPI_Recv(&tmp_tet[0], num_tets*4, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat);
  tf=num_sf_tri;
  
  
  for(int i=0; i<num_tets; i++) {
    tri_faces[tf][0] = tmp_tet[i][0] ;
    tri_faces[tf][1] = tmp_tet[i][1] ;
    tri_faces[tf][2] = tmp_tet[i][3] ;
    tri_faces[tf++][3] = cellnum ;
      
    tri_faces[tf][0] = tmp_tet[i][1] ;
    tri_faces[tf][1] = tmp_tet[i][2] ;
    tri_faces[tf][2] = tmp_tet[i][3] ; 
    tri_faces[tf++][3] = cellnum ;    
   
    tri_faces[tf][0] = tmp_tet[i][3] ;
    tri_faces[tf][1] = tmp_tet[i][2] ;
    tri_faces[tf][2] = tmp_tet[i][0] ;
    tri_faces[tf++][3] = cellnum ;
  
    tri_faces[tf][0] = tmp_tet[i][0] ;
    tri_faces[tf][1] = tmp_tet[i][2] ;
    tri_faces[tf][2] = tmp_tet[i][1] ;
    tri_faces[tf++][3] = cellnum ;

    cellnum++ ;
  }
 
 
 
  if(myrank==npes-1)
    cellnum = Info[3]+1;
  else 
    cellnum = Info[3] + info_cellnum[1] + num_pent5*(npes - myrank -2)+1; 
  
  // Receive volume pent5 from PE 0
  // creating buffer for the receiving pent5 elements

  autoArray<fivetuple> tmp_pents5(new fivetuple[num_pent5]) ;
  
  MPI_Send(&tag, 1, MPI_INT, 0, 1, MPI_COMM_WORLD); 

  MPI_Recv(&tmp_pents5[0], num_pent5*5, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat);
  qf = num_sf_quad;

  for(int i=0; i<num_pent5; i++) { 
  
    tri_faces[tf][0] = tmp_pents5[i][4] ;
    tri_faces[tf][1] = tmp_pents5[i][1] ;
    tri_faces[tf][2] = tmp_pents5[i][2] ;
    tri_faces[tf++][3] = cellnum ;
   
    tri_faces[tf][0] = tmp_pents5[i][4] ;
    tri_faces[tf][1] = tmp_pents5[i][2] ;
    tri_faces[tf][2] = tmp_pents5[i][3] ;
    tri_faces[tf++][3] = cellnum ;
   
    tri_faces[tf][0] = tmp_pents5[i][3] ;
    tri_faces[tf][1] = tmp_pents5[i][2] ;
    tri_faces[tf][2] = tmp_pents5[i][0] ;
    tri_faces[tf++][3] = cellnum ;
   
    tri_faces[tf][0] = tmp_pents5[i][0] ;
    tri_faces[tf][1] = tmp_pents5[i][2] ;
    tri_faces[tf][2] = tmp_pents5[i][1] ;
    tri_faces[tf++][3] = cellnum ;
    
    quad_faces[qf][0] = tmp_pents5[i][0] ;
    quad_faces[qf][1] = tmp_pents5[i][1] ;
    quad_faces[qf][2] = tmp_pents5[i][4] ;
    quad_faces[qf][3] = tmp_pents5[i][3] ;
    quad_faces[qf++][4] = cellnum ;
    cellnum++ ;
  } 
   



  if(myrank==npes-1)
    cellnum = Info[4]+1;
  else 
    cellnum = Info[4] + info_cellnum[2] + num_pent6*(npes - myrank -2)+1; 

 
  // Receive volume pent6 from PE 0
  // creating buffer for the receiving pent6 elements

  autoArray<sixtuple> tmp_pents6(new sixtuple[num_pent6]) ;
  
  MPI_Send(&tag, 1, MPI_INT, 0, 1, MPI_COMM_WORLD); 

  MPI_Recv(&tmp_pents6[0], num_pent6*6, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat);
   
  for(int i=0; i<num_pent6; i++) {   
    tri_faces[tf][0] = tmp_pents6[i][3] ;
    tri_faces[tf][1] = tmp_pents6[i][4] ;
    tri_faces[tf][2] = tmp_pents6[i][5] ;
    tri_faces[tf++][3] = cellnum ;
   
    tri_faces[tf][0] = tmp_pents6[i][0] ;
    tri_faces[tf][1] = tmp_pents6[i][2] ;
    tri_faces[tf][2] = tmp_pents6[i][1] ;
    tri_faces[tf++][3] = cellnum ;
   
    quad_faces[qf][0] = tmp_pents6[i][0] ;
    quad_faces[qf][1] = tmp_pents6[i][1] ;
    quad_faces[qf][2] = tmp_pents6[i][4] ;
    quad_faces[qf][3] = tmp_pents6[i][3] ;
    quad_faces[qf++][4] = cellnum ;
    
    quad_faces[qf][0] = tmp_pents6[i][1] ;
    quad_faces[qf][1] = tmp_pents6[i][2] ;
    quad_faces[qf][2] = tmp_pents6[i][5] ;
    quad_faces[qf][3] = tmp_pents6[i][4] ;
    quad_faces[qf++][4] = cellnum ;
    
    quad_faces[qf][0] = tmp_pents6[i][3] ;
    quad_faces[qf][1] = tmp_pents6[i][5] ;
    quad_faces[qf][2] = tmp_pents6[i][2] ;
    quad_faces[qf][3] = tmp_pents6[i][0] ;
    quad_faces[qf++][4] = cellnum ;

    cellnum++ ;
  }


  if(myrank==npes-1)
    cellnum = Info[5]+1;
  else 
    cellnum = Info[5] + info_cellnum[3] + num_hex*(npes - myrank -2)+1; 

 
  // Receive volume Hex from PE 0
  // creating buffer for the receiving hex elements

  autoArray<eighttuple> tmp_hexs(new eighttuple[num_hex]) ;

  MPI_Send(&tag, 1, MPI_INT, 0, 1, MPI_COMM_WORLD); 
  
  MPI_Recv(&tmp_hexs[0], num_hex*8, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat);
   
  for(int i=0; i<num_hex; i++) {   
   
    quad_faces[qf][0] = tmp_hexs[i][0] ;
    quad_faces[qf][1] = tmp_hexs[i][1] ;
    quad_faces[qf][2] = tmp_hexs[i][5] ;
    quad_faces[qf][3] = tmp_hexs[i][4] ;
    quad_faces[qf++][4] = cellnum ;
    
    quad_faces[qf][0] = tmp_hexs[i][1] ;
    quad_faces[qf][1] = tmp_hexs[i][2] ;
    quad_faces[qf][2] = tmp_hexs[i][6] ;
    quad_faces[qf][3] = tmp_hexs[i][5] ;
    quad_faces[qf++][4] = cellnum ;
    
    quad_faces[qf][0] = tmp_hexs[i][2] ;
    quad_faces[qf][1] = tmp_hexs[i][3] ;
    quad_faces[qf][2] = tmp_hexs[i][7] ;
    quad_faces[qf][3] = tmp_hexs[i][6] ;
    quad_faces[qf++][4] = cellnum ;
    
    quad_faces[qf][0] = tmp_hexs[i][4] ;
    quad_faces[qf][1] = tmp_hexs[i][7] ;
    quad_faces[qf][2] = tmp_hexs[i][3] ;
    quad_faces[qf][3] = tmp_hexs[i][0] ;
    quad_faces[qf++][4] = cellnum ;
    
    quad_faces[qf][0] = tmp_hexs[i][4] ;
    quad_faces[qf][1] = tmp_hexs[i][5] ;
    quad_faces[qf][2] = tmp_hexs[i][6] ;
    quad_faces[qf][3] = tmp_hexs[i][7] ;
    quad_faces[qf++][4] = cellnum ;
    
    quad_faces[qf][0] = tmp_hexs[i][3] ;
    quad_faces[qf][1] = tmp_hexs[i][2] ;
    quad_faces[qf][2] = tmp_hexs[i][1] ;
    quad_faces[qf][3] = tmp_hexs[i][0] ;
    quad_faces[qf++][4] = cellnum ;
    
    cellnum++ ;
  }
  

  //    if((qf != size_t(num_quad_faces)*2) || tf != size_t(num_tri_faces)*2) {
  //   cerr << "face numbers not consistent!" ;
  //  exit(-1) ;
  // }  


  // prepare triangle faces (sort them)  
  for(size_t i=0;i<tf;++i) {
    // xdr numbers nodes from zero
    tri_faces[i][0] -= 1 ;
    tri_faces[i][1] -= 1 ;
    tri_faces[i][2] -= 1 ;
    
    if(tri_faces[i][0] > tri_faces[i][1])
      std::swap(tri_faces[i][0],tri_faces[i][1]) ;
    if(tri_faces[i][0] > tri_faces[i][2])
      std::swap(tri_faces[i][0],tri_faces[i][2]) ;
    if(tri_faces[i][1] > tri_faces[i][2])
      std::swap(tri_faces[i][1],tri_faces[i][2]) ;
  }

  // prepare quad faces (sort them, but be careful)
  for(size_t i=0;i<qf;++i) {
    // xdr numbers nodes from zero
    quad_faces[i][0] -=1 ;
    quad_faces[i][1] -=1 ;
    quad_faces[i][2] -=1 ;
    quad_faces[i][3] -=1 ;
    // First make sure first entry is lowest number
    int tmp_face[4] ;
    int vs = quad_faces[i][0] ;
    size_t nv = 0 ;
    for(size_t j=1;j<4;++j)
      if(vs > quad_faces[i][j]) {
        vs = quad_faces[i][j] ;
        nv = j ;
      }
    for(size_t j=0;j<4;++j)
      tmp_face[j] = quad_faces[i][(j+nv)&0x3] ;
    // next make orientation so that it will match other face 
    if(tmp_face[1] < tmp_face[3])
      for(int j=0;j<4;++j)
        quad_faces[i][j] = tmp_face[j] ;
    else
      for(size_t j=0;j<4;++j)
        quad_faces[i][j] = tmp_face[(4 - j) &0x3 ] ;
  }


  if(myId == 0) 
    cout << "sorting faces..." << endl ;

  /*************************************/
  // BUCKET SORT FOR THE TRIAS 


  std::sort(&tri_faces[0],&tri_faces[tf],tri_sort) ;
  std::sort(&quad_faces[0],&quad_faces[qf],quad_sort) ;


  num = Info[0]/npes ;
  int rem = Info[0] % npes ;
  int cnt=0;  

  // Count array determines which node belongs to which PE'S
  // For the purpose of determining the splitters in the array     
 
  autoArray<int> count(new int[npes*2]) ;
  for(int j=0; j<npes*2; j++) {
    count[j]=0;
  } 
      
  for(int i = 0; i<npes; i++ ) {	
    if(i==npes-1 && rem>0) {
      count[cnt]=i*num;	
      cnt++;	
      count[cnt]=Info[0];
    } else {		
      count[cnt]=i*num;
      cnt++;
      count[cnt]=(i+1)*num;
      cnt++;
    }	
  }	
 

 
  // Determines the number of elements that would be broadcasted
  // Array sdisp keep track where in the tri_faces array chunk changes i.e splitters are    


  autoArray<int> splitter(new int[npes]);
  autoArray<int> recv_num_elmn(new int[npes]);  
  autoArray<int> sdisp_tri(new int[npes]);
  autoArray<int> rdisp_tri(new int[npes]);

  for(int j=0; j<npes; j++) {
    splitter[j] = 0;
  }  
  
  for(int i=0; i < (int)tf ; i++) {	
    for(int j=0; j< npes; j++) {
      
      if((tri_faces[i][0]>count[j*2] || tri_faces[i][0]==count[j*2])&&( tri_faces[i][0]<count[j*2+1])) {
        splitter[j]++;       
        
      }
    }             
  }    


  // Calculating displacement for sending the elements of the array
  for(int i=0;i<npes;i++)
    splitter[i]=splitter[i]*4;	   

  sdisp_tri[0]=0;
  for(int i=1; i<npes; i++) {
    sdisp_tri[i] = splitter[i-1]+sdisp_tri[i-1];
  } 


  // Using All to All personalized communication 
  MPI_Alltoall(&splitter[0], 1, MPI_INT, &recv_num_elmn[0], 1, MPI_INT, MPI_COMM_WORLD);  
 

  int sum_tri=0;
  for(int i=0;i<npes;i++) 
    sum_tri = sum_tri + recv_num_elmn[i]/4;  
 
 
  rdisp_tri[0] = 0;
 
  // Calculating the displacement for storing the array

  for(int i=1; i<npes; i++) {
    rdisp_tri[i] = recv_num_elmn[i-1] + rdisp_tri[i-1] ; 
  }	
 

  // initialize array 
  autoArray<Array<int,4> > new_tri_faces(new Array<int,4>[sum_tri]) ;
 
 
  MPI_Alltoallv(&tri_faces[0], &splitter[0], &sdisp_tri[0], MPI_INT, &new_tri_faces[0], &recv_num_elmn[0], &rdisp_tri[0], MPI_INT, MPI_COMM_WORLD);  
 

  // final sort 
  std::sort(&new_tri_faces[0],&new_tri_faces[sum_tri], tri_sort);        

  /****************************************************************/
  //         BUCKET SORT FOR THE QUADS


  // similar procedure for quad array 
  // Determines the number of elements that would be broadcasted
 
  autoArray<int> splitter_quad(new int[npes]); 
  autoArray<int> sdisp_quad(new int[npes]);
  autoArray<int> rdisp_quad(new int[npes]);
  autoArray<int> recv_elmn_quad(new int[npes]);
 
  for(int j=0; j<npes; j++) {
    splitter_quad[j] = 0;
  }
  // Calculating displacement for sending the elements of the array 
  // Calculating the number of elements to be sent to each PE  
  for(int i=0; i < (int)qf ; i++) {
    for(int j=0; j< npes; j++) {
	
      if((quad_faces[i][0]>count[j*2] || quad_faces[i][0]==count[j*2])&& quad_faces[i][0]< count[j*2+1]) {
        splitter_quad[j]++; 
	   
      }
    }             
  }        

   
  // Calculating displacement for sending the elements of the array
  for(int i=0;i<npes;i++)
    splitter_quad[i]=splitter_quad[i]*5;	   
   
  sdisp_quad[0]=0;
  for(int i=1; i<npes; i++) {
    sdisp_quad[i] = splitter_quad[i-1]+sdisp_quad[i-1];
  } 
    

  // Using All to All personalized communication 
  MPI_Alltoall(&splitter_quad[0], 1, MPI_INT, &recv_elmn_quad[0], 1, MPI_INT, MPI_COMM_WORLD);  
 

  // calculating the size of the receiving array  
  int sum_quad=0;
  for(int i=0;i<npes;i++) 
    sum_quad = sum_quad + recv_elmn_quad[i]/5;  
 
  // Calculating the displacement for storing in the array
  
  rdisp_quad[0] = 0;
  for(int i=1; i<npes; i++) {
    rdisp_quad[i] = recv_elmn_quad[i-1] + rdisp_quad[i-1] ;
  }	 
 

  // initialize array 
  autoArray<Array<int,5> > new_quad_faces(new Array<int,5>[sum_quad]) ;
  
  MPI_Alltoallv(&quad_faces[0], &splitter_quad[0], &sdisp_quad[0], MPI_INT, &new_quad_faces[0], &recv_elmn_quad[0], &rdisp_quad[0], MPI_INT, MPI_COMM_WORLD);  
 
  // final sort  
  std::sort(&new_quad_faces[0], &new_quad_faces[sum_quad], quad_sort);    
 
  
  size_t btf = 0 ;
  size_t bqf = 0 ;
  size_t itf = 0 ;
  size_t iqf = 0 ;
    
  for(int i=0; i<sum_tri/2 ; ++i) {
    if(new_tri_faces[i*2][0] != new_tri_faces[i*2+1][0] ||
       new_tri_faces[i*2][1] != new_tri_faces[i*2+1][1] ||
       new_tri_faces[i*2][2] != new_tri_faces[i*2+1][2] ||
       (new_tri_faces[i*2][3] < 0 && new_tri_faces[i*2+1][3] < 0)) {
      cerr << "trouble matching triangle faces! " << endl ;
      cerr << "perhaps an embedded surface remains in the grid?" << endl ;
      exit(-1) ;
    }
    
    if(new_tri_faces[i*2][3] < 0)   
      btf++;
  }

 
  int ntri;
  if(sum_tri==0)
    ntri=0;
  else
    ntri = sum_tri/2 - btf;
  
  autoArray<Array<int,5> > bnd_tri(new Array<int,5>[btf]) ;
 
  autoArray<Array<int,2> > int_tri_faces(new Array<int,2>[ntri]) ;
  autoArray<Array<int,3> > int_tri_nodes(new Array<int,3>[ntri]) ;
     
  btf=0;  

  for(int i=0; i<sum_tri/2 ; ++i) {
    if(new_tri_faces[i*2][3] < 0) {
      bnd_tri[btf][0] = new_tri_faces[i*2][0] ;
      bnd_tri[btf][1] = new_tri_faces[i*2][1] ;
      bnd_tri[btf][2] = new_tri_faces[i*2][2] ;
      bnd_tri[btf][3] = new_tri_faces[i*2+1][3] ;
      bnd_tri[btf++][4] = new_tri_faces[i*2][3] ;
    } else {
      int_tri_nodes[itf][0] = new_tri_faces[i*2][0] ;
      int_tri_nodes[itf][1] = new_tri_faces[i*2][1] ;
      int_tri_nodes[itf][2] = new_tri_faces[i*2][2] ;
      int_tri_faces[itf][0] = new_tri_faces[i*2+1][3] ;
      int_tri_faces[itf++][1] = new_tri_faces[i*2][3] ;
    }  
  }   
 
 
  int tum = (int)itf*2,flag=0; 

  // Sending size of the array to PE 0
  MPI_Send( &tum, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);

  // Receiving acknowledgement from pe0 to send the data
  MPI_Recv(&flag, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat); 

  // Sending the sorted internal tri faces

  if(flag == 1) 
    MPI_Send(&int_tri_faces[0], tum, MPI_INT, 0, 1, MPI_COMM_WORLD);
 
   
  for(int i=0;i<sum_quad/2;++i) {
    if(new_quad_faces[i*2][0] != new_quad_faces[i*2+1][0] ||
       new_quad_faces[i*2][1] != new_quad_faces[i*2+1][1] ||
       new_quad_faces[i*2][2] != new_quad_faces[i*2+1][2] ||
       new_quad_faces[i*2][3] != new_quad_faces[i*2+1][3] ||
       (new_quad_faces[i*2][4] < 0 && new_quad_faces[i*2+1][4] < 0)) {
      cerr << "trouble matching quad faces!" << endl ;
      cerr << "perhaps an embedded surface remains in the grid?" << endl ;
      exit(-1) ;
    }
    
    if(new_quad_faces[i*2][4] < 0)  
      bqf++   ;
  }  

  

  int nquad=sum_quad/2 - bqf;
 
  autoArray<Array<int,6> > bnd_quad(new Array<int,6>[bqf]) ;
 
  autoArray<Array<int,2> > int_quad_faces(new Array<int,2>[nquad]) ;
  autoArray<Array<int,4> > int_quad_nodes(new Array<int,4>[nquad]) ;
  
  int q=0;


  for(int i=0; i<sum_quad/2 ; ++i) {
    if(new_quad_faces[i*2][4] < 0) {
      bnd_quad[q][0] = new_quad_faces[i*2][0] ;
      bnd_quad[q][1] = new_quad_faces[i*2][1] ;
      bnd_quad[q][2] = new_quad_faces[i*2][2] ;   
      bnd_quad[q][3] = new_quad_faces[i*2][3] ;
      bnd_quad[q][4] = new_quad_faces[i*2+1][4] ;
      bnd_quad[q++][5] = new_quad_faces[i*2][4] ;
    } else {
      int_quad_nodes[iqf][0] = new_quad_faces[i*2][0] ;
      int_quad_nodes[iqf][1] = new_quad_faces[i*2][1] ;
      int_quad_nodes[iqf][2] = new_quad_faces[i*2][2] ;
      int_quad_nodes[iqf][3] = new_quad_faces[i*2][3] ;
      int_quad_faces[iqf][0] = new_quad_faces[i*2+1][4] ;
      int_quad_faces[iqf++][1] = new_quad_faces[i*2][4] ;
    }  
  }
 

  int sum = iqf*2;
  flag=0;    

  // Sending size of the array to PE 0
  MPI_Send(&sum, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);

 
  // Receiving acknowledgement from pe0 to send the data
  MPI_Recv(&flag, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat); 

  // Sending the sorted internal quad faces

  if(flag == 1)
    MPI_Send(&int_quad_faces[0], sum, MPI_INT, 0, 1, MPI_COMM_WORLD);
  
 
  sum = btf*5;
  
  // Sending size of the array to PE 0
  MPI_Send(&sum, 1, MPI_INT, 0, 3, MPI_COMM_WORLD);
  
  flag=0;
  
  // Receiving acknowledgement from pe to send the data
  MPI_Recv(&flag, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat); 
  
  // Sending the sorted  bnd trias
  if(flag==1)
    MPI_Send(&bnd_tri[0], sum, MPI_INT, 0, 1, MPI_COMM_WORLD);
  
  sum = bqf*6 ;
  flag=0;
  
  // Sending size of the array to pe 0
  MPI_Send(&sum, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
  
  // Receiving acknowledgement from pe to send the data
  MPI_Recv(&flag, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat); 
   
  // Sending the sorted bnd quads
  if(flag==1)
    MPI_Send( &bnd_quad[0], sum, MPI_INT, 0, 1, MPI_COMM_WORLD);
   
   
  sum = itf*3,flag=0;
   
  // Sending size of the array to PE 0
   
  MPI_Send(&sum, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
  
  // Receiving acknowledgement from pe to send the data
  MPI_Recv(&flag, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat); 
   

  // Sending the sorted internal tri nodes
  if(flag==1)
    MPI_Send(&int_tri_nodes[0], sum, MPI_INT, 0, 1, MPI_COMM_WORLD);
  
 

  sum = iqf*4,flag=0;
  
  // Sending size of the array to PE 0
  MPI_Send(&sum, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
   
  // Receiving acknowledgement from pe to send the data
  MPI_Recv(&flag, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat); 
   
  // Sending the sorted internal quad nodes
  if(flag==1)
    MPI_Send(&int_quad_nodes[0], sum, MPI_INT, 0, 1, MPI_COMM_WORLD);
  
}
