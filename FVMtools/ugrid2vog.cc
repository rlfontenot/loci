#include <Loci.h>
#include "vogtools.h"
#include <sys/stat.h>
#include <unistd.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <string>
#include <vector>

using std::cout ;
using std::endl ;
using std::cerr ;
using std::ios ;
using std::string ;
using std::vector ;
using std::istringstream ;
using Loci::Array ;
using Loci::MPI_rank ;
using Loci::MPI_processes ;

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

void readUGRID(string filename,bool binary, store<vector3d<double> > &pos,
               vector<Array<int,5> > &qfaces, vector<Array<int,4> > &tfaces,
               vector<Array<int,4> > &tets, vector<Array<int,5> > &pyramids,
               vector<Array<int,6> > &prisms, vector<Array<int,8> > &hexs) {

  pos.allocate(EMPTY) ;
  qfaces.clear() ;
  tfaces.clear() ;
  tets.clear() ;
  pyramids.clear() ;
  prisms.clear() ;
  hexs.clear() ;

  int num_nodes, num_sf_trias, num_sf_quads ;
  int num_vol_tets, num_vol_pents5, num_vol_pents6, num_vol_hexs ;

  FILE* IFP = NULL ;

  const int P = MPI_processes ;
  const int R = MPI_rank ;
  if(R == 0) {
    if(!binary)
      IFP = fopen(filename.c_str(), "r") ;
    else
      IFP = fopen(filename.c_str(), "rb") ;
    if(IFP == NULL) {
      cerr << "can't open '" << filename << "'" << endl ;
      exit(-1) ;
    }


    if(!binary) {
      fscanf(IFP, "%d%d%d", &num_nodes, & num_sf_trias, & num_sf_quads) ;
      fscanf(IFP, "%d%d%d%d", &num_vol_tets, &num_vol_pents5, &num_vol_pents6, &num_vol_hexs) ;
    } else {
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

    cout << "nnodes=" << num_nodes <<",ntria="<<num_sf_trias
         << ",nquad="<<num_sf_quads<<",ntets="<<num_vol_tets
         << ",npyrm="<<num_vol_pents5<<",nprsm="<<num_vol_pents6
         << ",nhex="<<num_vol_hexs << endl ;
  }

  Array<int,7> data ;
  data[0] = num_nodes ;
  data[1] = num_sf_trias ;
  data[2] = num_sf_quads ;
  data[3] = num_vol_tets ;
  data[4] = num_vol_pents5 ;
  data[5] = num_vol_pents6 ;
  data[6] = num_vol_hexs ;
  MPI_Bcast(&data[0],7,MPI_INT,0,MPI_COMM_WORLD) ;
  num_nodes      = data[0] ;
  num_sf_trias   = data[1] ;
  num_sf_quads   = data[2] ;
  num_vol_tets   = data[3] ;
  num_vol_pents5 = data[4] ;
  num_vol_pents6 = data[5] ;
  num_vol_hexs   = data[6] ;

  vector<int> node_ptns = VOG::simplePartitionVec(0,num_nodes-1,P) ;
  vector<entitySet> local_nodes(P) ;
  for(int i=0;i<P;++i)
    local_nodes[i] = interval(node_ptns[i],node_ptns[i+1]-1) ;

  pos.allocate(local_nodes[R]) ;

  if(R == 0) { // Read in positions

    FORALL(local_nodes[0],nd) {
      double ptmp[3] ;
      if(!binary) {
        for(int i = 0; i < 3; ++i) {
          fscanf(IFP, "%lf", &ptmp[i]) ;
        }
      } else {
        fread(ptmp,sizeof(double),3,IFP) ;
        if(reverse_byteorder)
          ug_io_reverse_byte_order(ptmp,sizeof(double),3) ;
      }
      pos[nd] = vector3d<double>(ptmp[0],ptmp[1],ptmp[2]) ;
    } ENDFORALL ;

    int mxsize = max(local_nodes[0].size(),
                     local_nodes[P-1].size()) ;
    vector<double> buf(mxsize*3) ;

    for(int i=1;i<P;++i) {
      int cnt = 0 ;
      FORALL(local_nodes[i],nd) {
        double ptmp[3] ;
        if(!binary) {
          for(int i = 0; i < 3; ++i) {
            fscanf(IFP, "%lf", &ptmp[i]) ;
          }
        } else {
          fread(ptmp,sizeof(double),3,IFP) ;
          if(reverse_byteorder)
            ug_io_reverse_byte_order(ptmp,sizeof(double),3) ;
        }
        buf[cnt++] = ptmp[0] ;
        buf[cnt++] = ptmp[1] ;
        buf[cnt++] = ptmp[2] ;
      } ENDFORALL ;
      MPI_Send(&buf[0],local_nodes[i].size()*3,MPI_DOUBLE,i,9,MPI_COMM_WORLD) ;
    }
  } else {
    MPI_Status status ;
    entitySet nodeSet = local_nodes[R] ;
    int recv_count = nodeSet.size()*3 ;
    vector<double> tmp_pos(recv_count) ;

    MPI_Recv(&tmp_pos[0],recv_count,MPI_DOUBLE,0,9,MPI_COMM_WORLD,&status) ;

    int tmp = 0 ;
    FORALL(nodeSet,nd) {
      vector3d<double> t(tmp_pos[tmp], tmp_pos[tmp+1], tmp_pos[tmp+2]) ;
      tmp += 3 ;
      pos[nd] = t ;
    } ENDFORALL ;

  }

  vector<int> triadist = VOG::simplePartitionVec(0,num_sf_trias-1,P) ;
  vector<int> quaddist = VOG::simplePartitionVec(0,num_sf_quads-1,P) ;

  qfaces = vector<Array<int,5> >(quaddist[R+1]-quaddist[R]) ;
  tfaces = vector<Array<int,4> >(triadist[R+1]-triadist[R]) ;

  // Read in boundary information
  if(R == 0) {
    { // triangles
      for(int i=triadist[R];i<triadist[R+1];++i) {
        tfaces[i][3] = 0 ;
        if(!binary)
          fscanf(IFP, "%d%d%d", &tfaces[i][0], &tfaces[i][1],
                 &tfaces[i][2]) ;
        else {
          fread(&tfaces[i][0], sizeof(int), 3, IFP) ;
          if(reverse_byteorder)
            ug_io_reverse_byte_order(&tfaces[i][0],sizeof(int),3) ;
        }
      }

      int tsz = max(triadist[1]-triadist[0],triadist[P]-triadist[P-1]) ;
      vector<Array<int,4> > ttmp(tsz) ;
      for(int p=1;p<P;++p) {
        int ltsz = triadist[p+1]-triadist[p] ;
        for(int i=0;i<ltsz;++i) {
          ttmp[i][3] = 0 ;
          if(!binary)
            fscanf(IFP, "%d%d%d", &ttmp[i][0], &ttmp[i][1], &ttmp[i][2]) ;
          else {
            fread(&ttmp[i][0], sizeof(int), 3, IFP) ;
            if(reverse_byteorder)
              ug_io_reverse_byte_order(&ttmp[i][0],sizeof(int),3) ;
          }
        }
        if(ltsz != 0)
          MPI_Send(&ttmp[0][0],ltsz*4,MPI_INT,p,8,MPI_COMM_WORLD) ;
      }
    }

    { // Quads
      for(int i=quaddist[R];i<quaddist[R+1];++i) {
        qfaces[i][4] = 0 ;
        if(!binary)
          fscanf(IFP, "%d%d%d%d", &qfaces[i][0], &qfaces[i][1], &qfaces[i][2],
                 &qfaces[i][3]) ;
        else {
          fread(&qfaces[i][0], sizeof(int), 4, IFP) ;
          if(reverse_byteorder)
            ug_io_reverse_byte_order(&qfaces[i][0],sizeof(int),4) ;
        }
      }

      int qsz = max(quaddist[1]-quaddist[0],quaddist[P]-quaddist[P-1]) ;
      vector<Array<int,5> > qtmp(qsz) ;
      for(int p=1;p<P;++p) {
        int lqsz = quaddist[p+1]-quaddist[p] ;
        for(int i=0;i<lqsz;++i) {
          qtmp[i][4] = 0 ;
          if(!binary)
            fscanf(IFP, "%d%d%d%d", &qtmp[i][0], &qtmp[i][1], &qtmp[i][2],
                   &qtmp[i][3]) ;
          else {
            fread(&qtmp[i][0], sizeof(int), 4, IFP) ;
            if(reverse_byteorder)
              ug_io_reverse_byte_order(&qtmp[i][0],sizeof(int),4) ;
          }
        }
        if(lqsz !=0)
          MPI_Send(&qtmp[0][0],lqsz*5,MPI_INT,p,8,MPI_COMM_WORLD) ;
      }
    }

    // Boundary tags

    { // triangles
      for(int i=triadist[R];i<triadist[R+1];++i) {
        if(!binary)
          fscanf(IFP, "%d", &tfaces[i][3]) ;
        else {
          fread(&tfaces[i][3], sizeof(int), 1, IFP) ;
          if(reverse_byteorder)
            ug_io_reverse_byte_order(&tfaces[i][3],sizeof(int),1) ;
        }
        tfaces[i][3] = -tfaces[i][3] ; // We use negative bc tags
      }

      int tsz = max(triadist[1]-triadist[0],triadist[P]-triadist[P-1]) ;
      vector<int> ttmp(tsz) ;
      for(int p=1;p<P;++p) {
        int ltsz = triadist[p+1]-triadist[p] ;
        for(int i=0;i<ltsz;++i) {
          if(!binary)
            fscanf(IFP, "%d", &ttmp[i]) ;
          else {
            fread(&ttmp[i], sizeof(int), 1, IFP) ;
            if(reverse_byteorder)
              ug_io_reverse_byte_order(&ttmp[i],sizeof(int),1) ;
          }
        }
        if(ltsz != 0)
          MPI_Send(&ttmp[0],ltsz,MPI_INT,p,8,MPI_COMM_WORLD) ;
      }
    }

    { // Quads
      for(int i=quaddist[R];i<quaddist[R+1];++i) {
        if(!binary)
          fscanf(IFP, "%d", &qfaces[i][4]) ;
        else {
          fread(&qfaces[i][4], sizeof(int), 1, IFP) ;
          if(reverse_byteorder)
            ug_io_reverse_byte_order(&qfaces[i][4],sizeof(int),1) ;
        }
        qfaces[i][4] = -qfaces[i][4] ;
      }

      int qsz = max(quaddist[1]-quaddist[0],quaddist[P]-quaddist[P-1]) ;
      vector<int> qtmp(qsz) ;
      for(int p=1;p<P;++p) {
        int lqsz = quaddist[p+1]-quaddist[p] ;
        for(int i=0;i<lqsz;++i) {
          if(!binary)
            fscanf(IFP, "%d", &qtmp[i]) ;
          else {
            fread(&qtmp[i], sizeof(int), 1, IFP) ;
            if(reverse_byteorder)
              ug_io_reverse_byte_order(&qtmp[i],sizeof(int),1) ;
          }
        }
        if(lqsz != 0)
          MPI_Send(&qtmp[0],lqsz,MPI_INT,p,8,MPI_COMM_WORLD) ;
      }
    }

  } else {
    MPI_Status status ;
    int tsz = tfaces.size() ;
    if(tsz != 0)
      MPI_Recv(&tfaces[0][0],tsz*4,MPI_INT,0,8,MPI_COMM_WORLD,&status) ;

    int qsz = qfaces.size();
    if(qsz != 0)
      MPI_Recv(&qfaces[0][0],qsz*5,MPI_INT,0,8,MPI_COMM_WORLD,&status) ;

    vector<int> tbuf(tsz) ;
    if(tsz !=0)
      MPI_Recv(&tbuf[0],tsz,MPI_INT,0,8,MPI_COMM_WORLD,&status) ;
    for(int i=0;i<tsz;++i)
      tfaces[i][3] = -tbuf[i] ;

    vector<int> qbuf(qsz) ;
    if(qsz != 0)
      MPI_Recv(&qbuf[0],qsz,MPI_INT,0,8,MPI_COMM_WORLD,&status) ;
    for(int i=0;i<qsz;++i)
      qfaces[i][4] = -qbuf[i] ;
  }

  // Read in volume elements

  // Read in tetrahedra
  vector<int> tetsdist = VOG::simplePartitionVec(0,num_vol_tets-1,P) ;
  tets = vector<Array<int,4> >(tetsdist[R+1]-tetsdist[R]) ;

  if(R == 0) {
    for(int i=tetsdist[R];i<tetsdist[R+1];++i) {
      if(!binary)
        fscanf(IFP, "%d%d%d%d", &tets[i][0], &tets[i][1],
               &tets[i][2], &tets[i][3]) ;
      else {
        fread(&tets[i][0], sizeof(int), 4, IFP) ;
        if(reverse_byteorder)
          ug_io_reverse_byte_order(&tets[i][0],sizeof(int),4) ;
      }
    }

    int tsz = max(tetsdist[1]-tetsdist[0],tetsdist[P]-tetsdist[P-1]) ;
    vector<Array<int,4> > ttmp(tsz) ;
    for(int p=1;p<P;++p) {
      int ltsz = tetsdist[p+1]-tetsdist[p] ;
      for(int i=0;i<ltsz;++i) {
        if(!binary)
          fscanf(IFP, "%d%d%d%d", &ttmp[i][0], &ttmp[i][1], &ttmp[i][2],
                 &ttmp[i][3]) ;
        else {
          fread(&ttmp[i][0], sizeof(int), 4, IFP) ;
          if(reverse_byteorder)
            ug_io_reverse_byte_order(&ttmp[i][0],sizeof(int),4) ;
        }
      }
      if(ltsz != 0)
        MPI_Send(&ttmp[0][0],ltsz*4,MPI_INT,p,7,MPI_COMM_WORLD) ;
    }
  } else {
    MPI_Status status ;
    int tsz = tets.size() ;
    if(tsz != 0)
      MPI_Recv(&tets[0][0],tsz*4,MPI_INT,0,7,MPI_COMM_WORLD,&status) ;
  }

  // Read in pyramids
  vector<int> pyrmdist = VOG::simplePartitionVec(0,num_vol_pents5-1,P) ;
  pyramids = vector<Array<int,5> >(pyrmdist[R+1]-pyrmdist[R]) ;

  if(R == 0) {
    for(int i=pyrmdist[R];i<pyrmdist[R+1];++i) {
      if(!binary)
        fscanf(IFP, "%d%d%d%d%d", &pyramids[i][0], &pyramids[i][1],
               &pyramids[i][2], &pyramids[i][3], &pyramids[i][4]) ;
      else {
        fread(&pyramids[i][0], sizeof(int), 5, IFP) ;
        if(reverse_byteorder)
          ug_io_reverse_byte_order(&pyramids[i][0],sizeof(int),5) ;
      }
    }

    int tsz = max(pyrmdist[1]-pyrmdist[0],pyrmdist[P]-pyrmdist[P-1]) ;
    vector<Array<int,5> > ttmp(tsz) ;
    for(int p=1;p<P;++p) {
      int ltsz = pyrmdist[p+1]-pyrmdist[p] ;
      for(int i=0;i<ltsz;++i) {
        if(!binary)
          fscanf(IFP, "%d%d%d%d%d", &ttmp[i][0], &ttmp[i][1], &ttmp[i][2],
                 &ttmp[i][3], &ttmp[i][4]) ;
        else {
          fread(&ttmp[i][0], sizeof(int), 5, IFP) ;
          if(reverse_byteorder)
            ug_io_reverse_byte_order(&ttmp[i][0],sizeof(int),5) ;
        }
      }
      if(ltsz != 0)
        MPI_Send(&ttmp[0][0],ltsz*5,MPI_INT,p,6,MPI_COMM_WORLD) ;
    }
  } else {
    MPI_Status status ;
    int tsz = pyramids.size() ;
    if(tsz != 0)
      MPI_Recv(&pyramids[0][0],tsz*5,MPI_INT,0,6,MPI_COMM_WORLD,&status) ;
  }

  // Read in prisms
  vector<int> prsmdist = VOG::simplePartitionVec(0,num_vol_pents6-1,P) ;
  prisms = vector<Array<int,6> >(prsmdist[R+1]-prsmdist[R]) ;

  if(R == 0) {
    for(int i=prsmdist[R];i<prsmdist[R+1];++i) {
      if(!binary)
        fscanf(IFP, "%d%d%d%d%d%d", &prisms[i][0], &prisms[i][1],
               &prisms[i][2], &prisms[i][3], &prisms[i][4],
               &prisms[i][5]) ;
      else {
        fread(&prisms[i][0], sizeof(int), 6, IFP) ;
        if(reverse_byteorder)
          ug_io_reverse_byte_order(&prisms[i][0],sizeof(int),6) ;
      }
    }

    int tsz = max(prsmdist[1]-prsmdist[0],prsmdist[P]-prsmdist[P-1]) ;
    vector<Array<int,6> > ttmp(tsz) ;
    for(int p=1;p<P;++p) {
      int ltsz = prsmdist[p+1]-prsmdist[p] ;
      for(int i=0;i<ltsz;++i) {
        if(!binary)
          fscanf(IFP, "%d%d%d%d%d%d", &ttmp[i][0], &ttmp[i][1], &ttmp[i][2],
                 &ttmp[i][3], &ttmp[i][4], &ttmp[i][5]) ;
        else {
          fread(&ttmp[i][0], sizeof(int), 6, IFP) ;
          if(reverse_byteorder)
            ug_io_reverse_byte_order(&ttmp[i][0],sizeof(int),6) ;
        }
      }
      if(ltsz != 0)
        MPI_Send(&ttmp[0][0],ltsz*6,MPI_INT,p,5,MPI_COMM_WORLD) ;
    }
  } else {
    MPI_Status status ;
    int tsz = prisms.size() ;
    if(tsz != 0)
      MPI_Recv(&prisms[0][0],tsz*6,MPI_INT,0,5,MPI_COMM_WORLD,&status) ;
  }

  // Read in hexahdra
  vector<int> hexsdist = VOG::simplePartitionVec(0,num_vol_hexs-1,P) ;
  hexs = vector<Array<int,8> >(hexsdist[R+1]-hexsdist[R]) ;

  if(R == 0) {
    for(int i=hexsdist[R];i<hexsdist[R+1];++i) {
      if(!binary)
        fscanf(IFP, "%d%d%d%d%d%d%d%d", &hexs[i][0], &hexs[i][1],
               &hexs[i][2], &hexs[i][3], &hexs[i][4],
               &hexs[i][5], &hexs[i][6], &hexs[i][7]) ;
      else {
        fread(&hexs[i][0], sizeof(int), 8, IFP) ;
        if(reverse_byteorder)
          ug_io_reverse_byte_order(&hexs[i][0],sizeof(int),8) ;
      }
    }

    int tsz = max(hexsdist[1]-hexsdist[0],hexsdist[P]-hexsdist[P-1]) ;
    vector<Array<int,8> > ttmp(tsz) ;
    for(int p=1;p<P;++p) {
      int ltsz = hexsdist[p+1]-hexsdist[p] ;
      for(int i=0;i<ltsz;++i) {
        if(!binary)
          fscanf(IFP, "%d%d%d%d%d%d%d%d", &ttmp[i][0], &ttmp[i][1],
                 &ttmp[i][2], &ttmp[i][3], &ttmp[i][4], &ttmp[i][5],
                 &ttmp[i][6], &ttmp[i][7]) ;
        else {
          fread(&ttmp[i][0], sizeof(int), 8, IFP) ;
          if(reverse_byteorder)
            ug_io_reverse_byte_order(&ttmp[i][0],sizeof(int),8) ;
        }
      }
      if(ltsz != 0)
        MPI_Send(&ttmp[0][0],ltsz*8,MPI_INT,p,4,MPI_COMM_WORLD) ;
    }
  } else {
    MPI_Status status ;
    int tsz = hexs.size() ;
    if(tsz != 0)
      MPI_Recv(&hexs[0][0],tsz*8,MPI_INT,0,4,MPI_COMM_WORLD,&status) ;
  }
}

struct quadFace {
  Array<int,4> nodes ;
  int cell ;
  bool left ;
} ;

struct triaFace {
  Array<int,3> nodes ;
  int cell ;
  bool left ;
} ;

inline bool quadCompare(const quadFace &f1, const quadFace &f2) {
  return ((f1.nodes[0] < f2.nodes[0]) ||
          (f1.nodes[0] == f2.nodes[0] && f1.nodes[1] < f2.nodes[1]) ||
          (f1.nodes[0] == f2.nodes[0] && f1.nodes[1] == f2.nodes[1] &&
           f1.nodes[2] < f2.nodes[2]) ||
          (f1.nodes[0] == f2.nodes[0] && f1.nodes[1] == f2.nodes[1] &&
           f1.nodes[2] == f2.nodes[2] && f1.nodes[3] < f2.nodes[3])) ;
}

inline bool triaCompare(const triaFace &f1, const triaFace &f2) {
  return ((f1.nodes[0] < f2.nodes[0]) ||
          (f1.nodes[0] == f2.nodes[0] && f1.nodes[1] < f2.nodes[1]) ||
          (f1.nodes[0] == f2.nodes[0] && f1.nodes[1] == f2.nodes[1] &&
           f1.nodes[2] < f2.nodes[2])) ;
}

void convert2face(store<vector3d<double> > &pos,
                  vector<Array<int,5> > &qfaces, vector<Array<int,4> > &tfaces,
                  vector<Array<int,4> > &tets, vector<Array<int,5> > &pyramids,
                  vector<Array<int,6> > &prisms, vector<Array<int,8> > &hexs,
                  multiMap &face2node,Map &cl, Map &cr) {
  int maxid = 0 ;
  entitySet posDom = pos.domain() ;
  if(posDom != EMPTY)
    maxid = posDom.Max()+1 ;
  int cellid ;

  MPI_Allreduce(&maxid,&cellid,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD) ;
  int cellbase = cellid ;
  int ncells = tets.size()+pyramids.size()+prisms.size()+hexs.size() ;
  vector<int> cellsizes(MPI_processes) ;
  MPI_Allgather(&ncells,1,MPI_INT,&cellsizes[0],1,MPI_INT,MPI_COMM_WORLD) ;
  for(int i=0;i<MPI_rank;++i)
    cellid += cellsizes[i] ;

  int num_quad_faces =
    qfaces.size() + pyramids.size()+ prisms.size()*3 + hexs.size()*6 ;

  int num_tria_faces =
    tfaces.size() + tets.size()*4 + pyramids.size()*4 + prisms.size()*2 ;

  vector<triaFace> tria(num_tria_faces) ;
  vector<quadFace> quad(num_quad_faces) ;
  int tf = 0 ;
  int qf = 0 ;
  for(size_t i=0;i<tfaces.size();++i) {
    for(int n=0;n<3;++n)
      tria[tf].nodes[n] = tfaces[i][n] ;
    tria[tf].left = true ;
    tria[tf++].cell = tfaces[i][3] ;
  }
  for(size_t i=0;i<qfaces.size();++i) {
    for(int n=0;n<4;++n)
      quad[qf].nodes[n] = qfaces[i][n] ;
    quad[qf].left = true ;
    quad[qf++].cell = qfaces[i][4] ;
  }

  // Create faces generated by tetrahedra
  for(size_t i=0;i<tets.size();++i) {
    tria[tf].nodes[0] = tets[i][0] ;
    tria[tf].nodes[1] = tets[i][1] ;
    tria[tf].nodes[2] = tets[i][3] ;
    tria[tf].left = true ;
    tria[tf++].cell = cellid ;

    tria[tf].nodes[0] = tets[i][1] ;
    tria[tf].nodes[1] = tets[i][2] ;
    tria[tf].nodes[2] = tets[i][3] ;
    tria[tf].left = true ;
    tria[tf++].cell = cellid ;
    

    tria[tf].nodes[0] = tets[i][3] ;
    tria[tf].nodes[1] = tets[i][2] ;
    tria[tf].nodes[2] = tets[i][0] ;
    tria[tf].left = true ;
    tria[tf++].cell = cellid ;

    tria[tf].nodes[0] = tets[i][0] ;
    tria[tf].nodes[1] = tets[i][2] ;
    tria[tf].nodes[2] = tets[i][1] ;
    tria[tf].left = true ;
    tria[tf++].cell = cellid ;

    cellid++ ;
  }

  // create faces generated by pyramids
  for(size_t i=0;i<pyramids.size();++i) {
    tria[tf].nodes[0] = pyramids[i][4] ;
    tria[tf].nodes[1] = pyramids[i][1] ;
    tria[tf].nodes[2] = pyramids[i][2] ;
    tria[tf].left = true ;
    tria[tf++].cell = cellid ;

    tria[tf].nodes[0] = pyramids[i][4] ;
    tria[tf].nodes[1] = pyramids[i][2] ;
    tria[tf].nodes[2] = pyramids[i][3] ;
    tria[tf].left = true ;
    tria[tf++].cell = cellid ;

    tria[tf].nodes[0] = pyramids[i][3] ;
    tria[tf].nodes[1] = pyramids[i][2] ;
    tria[tf].nodes[2] = pyramids[i][0] ;
    tria[tf].left = true ;
    tria[tf++].cell = cellid ;

    tria[tf].nodes[0] = pyramids[i][0] ;
    tria[tf].nodes[1] = pyramids[i][2] ;
    tria[tf].nodes[2] = pyramids[i][1] ;
    tria[tf].left = true ;
    tria[tf++].cell = cellid ;

    quad[qf].nodes[0] = pyramids[i][0] ;
    quad[qf].nodes[1] = pyramids[i][1] ;
    quad[qf].nodes[2] = pyramids[i][4] ;
    quad[qf].nodes[3] = pyramids[i][3] ;
    quad[qf].left = true ;
    quad[qf++].cell = cellid ;
    cellid++ ;
  }

  // create faces generated by prisms
  for(size_t i=0;i<prisms.size();++i) {
    tria[tf].nodes[0] = prisms[i][3] ;
    tria[tf].nodes[1] = prisms[i][4] ;
    tria[tf].nodes[2] = prisms[i][5] ;
    tria[tf].left = true ;
    tria[tf++].cell = cellid ;

    tria[tf].nodes[0] = prisms[i][0] ;
    tria[tf].nodes[1] = prisms[i][2] ;
    tria[tf].nodes[2] = prisms[i][1] ;
    tria[tf].left = true ;
    tria[tf++].cell = cellid ;

    quad[qf].nodes[0] = prisms[i][0] ;
    quad[qf].nodes[1] = prisms[i][1] ;
    quad[qf].nodes[2] = prisms[i][4] ;
    quad[qf].nodes[3] = prisms[i][3] ;
    quad[qf].left = true ;
    quad[qf++].cell = cellid ;

    quad[qf].nodes[0] = prisms[i][1] ;
    quad[qf].nodes[1] = prisms[i][2] ;
    quad[qf].nodes[2] = prisms[i][5] ;
    quad[qf].nodes[3] = prisms[i][4] ;
    quad[qf].left = true ;
    quad[qf++].cell = cellid ;

    quad[qf].nodes[0] = prisms[i][3] ;
    quad[qf].nodes[1] = prisms[i][5] ;
    quad[qf].nodes[2] = prisms[i][2] ;
    quad[qf].nodes[3] = prisms[i][0] ;
    quad[qf].left = true ;
    quad[qf++].cell = cellid ;

    cellid++ ;
  }

  // create faces generated by hexahedra
  for(size_t i=0;i<hexs.size();++i) {
    quad[qf].nodes[0] = hexs[i][0] ;
    quad[qf].nodes[1] = hexs[i][1] ;
    quad[qf].nodes[2] = hexs[i][5] ;
    quad[qf].nodes[3] = hexs[i][4] ;
    quad[qf].left = true ;
    quad[qf++].cell = cellid ;

    quad[qf].nodes[0] = hexs[i][1] ;
    quad[qf].nodes[1] = hexs[i][2] ;
    quad[qf].nodes[2] = hexs[i][6] ;
    quad[qf].nodes[3] = hexs[i][5] ;
    quad[qf].left = true ;
    quad[qf++].cell = cellid ;

    quad[qf].nodes[0] = hexs[i][2] ;
    quad[qf].nodes[1] = hexs[i][3] ;
    quad[qf].nodes[2] = hexs[i][7] ;
    quad[qf].nodes[3] = hexs[i][6] ;
    quad[qf].left = true ;
    quad[qf++].cell = cellid ;

    quad[qf].nodes[0] = hexs[i][4] ;
    quad[qf].nodes[1] = hexs[i][7] ;
    quad[qf].nodes[2] = hexs[i][3] ;
    quad[qf].nodes[3] = hexs[i][0] ;
    quad[qf].left = true ;
    quad[qf++].cell = cellid ;

    quad[qf].nodes[0] = hexs[i][4] ;
    quad[qf].nodes[1] = hexs[i][5] ;
    quad[qf].nodes[2] = hexs[i][6] ;
    quad[qf].nodes[3] = hexs[i][7] ;
    quad[qf].left = true ;
    quad[qf++].cell = cellid ;

    quad[qf].nodes[0] = hexs[i][3] ;
    quad[qf].nodes[1] = hexs[i][2] ;
    quad[qf].nodes[2] = hexs[i][1] ;
    quad[qf].nodes[3] = hexs[i][0] ;
    quad[qf].left = true ;
    quad[qf++].cell = cellid ;

    cellid++ ;
  }

  if(qf != num_quad_faces) {
    cerr << "internal consistency error on quad faces" << endl ;
    Loci::Abort() ;
  }
  if(tf != num_tria_faces) {
    cerr << "internal consistency error on triangle faces" << endl ;
    Loci::Abort() ;
  }

  // prepare triangle faces (sort them)
  for(size_t i=0;i<tria.size();++i) {
    // pos numbers nodes from zero
    tria[i].nodes[0] -= 1 ;
    tria[i].nodes[1] -= 1 ;
    tria[i].nodes[2] -= 1 ;

    if(tria[i].nodes[0] > tria[i].nodes[1]) {
      std::swap(tria[i].nodes[0],tria[i].nodes[1]) ;
      tria[i].left = !tria[i].left ;
    }
    if(tria[i].nodes[0] > tria[i].nodes[2]) {
      std::swap(tria[i].nodes[0],tria[i].nodes[2]) ;
      tria[i].left = !tria[i].left ;
    }
    if(tria[i].nodes[1] > tria[i].nodes[2]) {
      std::swap(tria[i].nodes[1],tria[i].nodes[2]) ;
      tria[i].left = !tria[i].left ;
    }
  }

  // prepare quad faces (sort them, but be careful)
  for(size_t i=0;i<quad.size();++i) {
    // pos numbers nodes from zero
    quad[i].nodes[0] -=1 ;
    quad[i].nodes[1] -=1 ;
    quad[i].nodes[2] -=1 ;
    quad[i].nodes[3] -=1 ;
    // First make sure first entry is lowest number
    int tmp_face[4] ;
    int vs = quad[i].nodes[0] ;
    size_t nv = 0 ;
    for(size_t j=1;j<4;++j)
      if(vs > quad[i].nodes[j]) {
        vs = quad[i].nodes[j] ;
        nv = j ;
      }
    for(size_t j=0;j<4;++j)
      tmp_face[j] = quad[i].nodes[(j+nv)&0x3] ;
    // next make orientation so that it will match other face
    if(tmp_face[1] < tmp_face[3])
      for(int j=0;j<4;++j)
        quad[i].nodes[j] = tmp_face[j] ;
    else {
      for(size_t j=0;j<4;++j)
        quad[i].nodes[j] = tmp_face[(4 - j) &0x3 ] ;
      quad[i].left = !quad[i].left ;
    }
  }

  int tsz = tria.size() ;
  int mtsz ;
  MPI_Allreduce(&tsz,&mtsz,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD) ;
  if(mtsz != 0) {
    Loci::parSampleSort(tria,triaCompare,MPI_COMM_WORLD) ;
  }

  int qsz = quad.size() ;
  int mqsz ;
  MPI_Allreduce(&qsz,&mqsz,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD) ;
  if(mqsz != 0) {
    Loci::parSampleSort(quad,quadCompare,MPI_COMM_WORLD) ;
  }
  if((tria.size() & 1) == 1) {
    cerr << "non-even number of triangle faces! inconsistent!" << endl ;
    Loci::Abort() ;
  }
  if((quad.size() & 1 ) == 1) {
    cerr << "non-even number of quad faces! inconcistent!" << endl ;
    Loci::Abort() ;
  }

  int ntria = tria.size()/2 ;
  int nquad = quad.size()/2 ;

  int nfaces = ntria+nquad ;

  ncells = 0 ;
  for(int i=0;i<MPI_processes;++i)
    ncells += cellsizes[i] ;
  int facebase = cellbase + ncells ;
  vector<int> facesizes(MPI_processes) ;
  MPI_Allgather(&nfaces,1,MPI_INT,&facesizes[0],1,MPI_INT,MPI_COMM_WORLD) ;
  for(int i=0;i<MPI_rank;++i)
    facebase += facesizes[i] ;
  entitySet faces = interval(facebase,facebase+nfaces-1) ;
  store<int> count ;
  cl.allocate(faces) ;
  cr.allocate(faces) ;
  count.allocate(faces) ;
  int fc = facebase ;
  for(int i=0;i<ntria;i++)
    count[fc++] = 3 ;
  for(int i=0;i<nquad;i++)
    count[fc++] = 4 ;
  face2node.allocate(count) ;
  fc = facebase ;
  for(int i=0;i<ntria;++i) {
    for(int j=0;j<3;++j) {
      face2node[fc][j] = tria[i*2].nodes[j] ;
      FATAL(tria[i*2].nodes[j] != tria[i*2+1].nodes[j]) ;
    }
    int c1 = tria[i*2].cell ;
    int c2 = tria[i*2+1].cell ;
    if(c1 < 0 && c2 < 0) {
      cerr << "two boundary faces glued together, probably a tranparent surface is causing the problem!"<< endl ;
      Loci::Abort() ;
    }
    if(c1 < 0) {
      cl[fc] = c2 ;
      cr[fc] = c1 ;
      if(tria[i*2+1].left)
        for(int j=0;j<3;++j)
          face2node[fc][j] = tria[i*2+1].nodes[j] ;
      else
        for(int j=0;j<3;++j)
          face2node[fc][j] = tria[i*2+1].nodes[2-j] ;
    } else if(c2 < 0) {
      cl[fc] = c1 ;
      cr[fc] = c2 ;
      if(tria[i*2].left)
        for(int j=0;j<3;++j)
          face2node[fc][j] = tria[i*2].nodes[j] ;
      else
        for(int j=0;j<3;++j)
          face2node[fc][j] = tria[i*2].nodes[2-j] ;
    } else {
      if(!tria[i*2].left) 
        std::swap(c1,c2) ;
      cl[fc] = c1 ;
      cr[fc] = c2 ;
      if((tria[i*2].left && tria[i*2+1].left) ||
         (!tria[i*2].left && !tria[i*2+1].left)) {
        cerr << "consistency error" << endl ;
      }
    }
    fc++ ;
  }
  for(int i=0;i<nquad;++i) {
    for(int j=0;j<4;++j) {
      face2node[fc][j] = quad[i*2].nodes[j] ;
      FATAL(quad[i*2].nodes[j] != quad[i*2+1].nodes[j]) ;
    }
    int c1 = quad[i*2].cell ;
    int c2 = quad[i*2+1].cell ;
    if(c1 < 0 && c2 < 0) {
      cerr << "two boundary faces glued together, probably a tranparent surface is causing the problem!"<< endl ;
      Loci::Abort() ;
    }
    
    if(c1 < 0) {
      cl[fc] = c2 ;
      cr[fc] = c1 ;
      if(quad[i*2+1].left)
        for(int j=0;j<4;++j)
          face2node[fc][j] = quad[i*2+1].nodes[j] ;
      else
        for(int j=0;j<4;++j)
          face2node[fc][j] = quad[i*2+1].nodes[3-j] ;
    } else if(c2 < 0) {
      cl[fc] = c1 ;
      cr[fc] = c2 ;
      if(quad[i*2].left)
        for(int j=0;j<4;++j)
          face2node[fc][j] = quad[i*2].nodes[j] ;
      else
        for(int j=0;j<4;++j)
          face2node[fc][j] = quad[i*2].nodes[3-j] ;
    } else {
      if(!quad[i*2].left) 
        std::swap(c1,c2) ;
      cl[fc] = c1 ;
      cr[fc] = c2 ;
      if((quad[i*2].left && quad[i*2+1].left) ||
         (!quad[i*2].left && !quad[i*2+1].left)) {
        cerr << "consistency error" << endl ;
      }
    }
    fc++ ;
  }

}

int main(int ac, char* av[]) {
  using namespace VOG ;

  bool optimize = true ;

  Loci::Init(&ac,&av) ;
  const char *filename ;
  std::string tmp_str ;
  bool binary = 0;
  string Lref = "1 meter" ;
  while(ac>=2 && av[1][0] == '-') {
    // If user specifies an alternate query, extract it from the
    // command line.
    if(ac >= 3 && !strcmp(av[1],"-Lref")) {
      Lref = av[2] ;
      ac -= 2 ;
      av += 2 ;
    } else if(ac >= 2 && !strcmp(av[1],"-v")) {
      cout << "Loci version: " << Loci::version() << endl ;
      if(ac == 2) {
        Loci::Finalize() ;
        exit(0) ;
      }
      ac-- ;
      av++ ;
    } else if(ac >= 2 && !strcmp(av[1],"-b")) {
      binary = 1 ;
      ac-- ;
      av++ ;
    } else if(ac >= 2 && !strcmp(av[1],"-o")) {
      optimize = false ;
      ac-- ;
      av++ ;
    } else if(ac >= 2 && !strcmp(av[1],"-in")) {
      Lref = "1 inch" ;
      ac-- ;
      av++ ;
    } else if(ac >= 2 && !strcmp(av[1],"-ft")) {
      Lref = "1 foot" ;
      ac-- ;
      av++ ;
    } else if(ac >= 2 && !strcmp(av[1],"-cm")) {
      Lref = "1 centimeter" ;
      ac-- ;
      av++ ;
    } else {
      cerr << "argument " << av[1] << " is not understood." << endl ;
      ac-- ;
      av++ ;
    }
  }


  if(ac == 2) {
    tmp_str.append(av[1]) ;
  } else {
    cerr << "Usage: ugrid2vog <options> <file>" << endl
         << "Where options are listed below and <file> is the filename sans postfix" << endl
         << "flags:" << endl
         << "  -b  : assume input file is binary" << endl
         << "  -o  : disable optimization that reorders nodes and faces" << endl
         << "  -v  : display version" << endl
         << "  -in : input grid is in inches" << endl
         << "  -ft : input grid is in feet" << endl
         << "  -cm : input grid is in centimeters" << endl
         << "  -Lref <units> : 1 unit in input grid is <units> long" << endl
         << endl ;
      exit(-1) ;
  }

  if(Lref == "")
    Lref = "1 meter" ;
  
  if(!isdigit(Lref[0])) {
    Lref = string("1") + Lref ;
  }

  Loci::UNIT_type tp ;
  istringstream iss(Lref) ;
  iss >> tp ;
  double posScale = tp.get_value_in("meter") ;

  if(Loci::MPI_rank == 0) {
    cout << "input grid file units = " << tp ;
    if(posScale != 1.0) 
      cout << " = " << posScale << " meters " ;
    cout << endl ;
  }
  
  // Check machines internal byte order
  check_order() ;

  int loc = 0;
  loc = tmp_str.find('.') ;
  std::string new_str = tmp_str.substr(0, loc) ;
  filename = new_str.c_str() ;
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
  string infile = buf ;

  string outfile = string(filename) + string(".vog") ;

  store<vector3d<double> > pos ;
  vector<Array<int,5> > qfaces ;
  vector<Array<int,4> > tfaces ;
  vector<Array<int,4> > tets ;
  vector<Array<int,5> > pyramids ;
  vector<Array<int,6> > prisms ;
  vector<Array<int,8> > hexs ;

  readUGRID(infile, binary, pos,qfaces,tfaces,tets,pyramids,prisms,hexs) ;

  if(posScale != 1.0) {
    FORALL(pos.domain(),nd) {
      pos[nd] *= posScale ;
    } ENDFORALL ;
  }
  
  vector<BC_descriptor> bcs ;
  vector<int> transsurf ;

  if(MPI_rank == 0) {
    string tagsfile = string(filename) + ".tags" ;
    bcs = readTags(tagsfile) ;
    if(bcs.size() == 0) {
      cerr << "unable to read '" << tagsfile << "'" << endl ;
    }

    for(size_t i=0;i<bcs.size();++i)
      if(bcs[i].Trans)
        transsurf.push_back(bcs[i].id) ;

    cout << "boundary faces:" ;
    for(size_t i=0;i<bcs.size();++i)
      cout << ' ' << bcs[i].name ;
    cout << endl ;
  }
  int trans_size = transsurf.size() ;
  MPI_Bcast(&trans_size,1,MPI_INT,0,MPI_COMM_WORLD) ;
  if(trans_size != 0) {
    if(MPI_rank != 0)
      transsurf = vector<int>(trans_size) ;
    MPI_Bcast(&transsurf[0],trans_size,MPI_INT,0,MPI_COMM_WORLD) ;

    if(MPI_rank == 0)
      cout << "removing transparent surfaces" << endl ;

    vector<Array<int,5> > qtfaces ;
    for(size_t i=0;i<qfaces.size();++i) {
      bool trans = false ;
      for(int j=0;j<trans_size;++j)
        if(qfaces[i][4] == -transsurf[j])
          trans = true ;
      if(!trans)
        qtfaces.push_back(qfaces[i]) ;
    }
    qfaces.swap(qtfaces) ;

    vector<Array<int,4> > ttfaces ;
    for(size_t i=0;i<tfaces.size();++i) {
      bool trans = false ;
      for(int j=0;j<trans_size;++j)
        if(tfaces[i][3] == -transsurf[j])
          trans = true ;
      if(!trans)
        ttfaces.push_back(tfaces[i]) ;
    }
    tfaces.swap(ttfaces) ;
  }
  vector<pair<int,string> > surf_ids ;
  for(size_t i=0;i<bcs.size();++i)
    surf_ids.push_back(pair<int,string>(bcs[i].id,bcs[i].name)) ;

  multiMap face2node ;
  Map cl,cr ;

  convert2face(pos,qfaces,tfaces,tets,pyramids,prisms,hexs,
               face2node,cl,cr) ;

  // This code is no longer needed, the face orientation is established
  // in convert2face
  
  // establish face left-right orientation
  //  if(MPI_rank == 0)
  //    cerr << "orienting faces" << endl ;
  //  VOG::orientFaces(pos,cl,cr,face2node) ;

  if(MPI_rank == 0)
    cerr << "coloring matrix" << endl ;
  VOG::colorMatrix(pos,cl,cr,face2node) ;

  if(optimize) {
    if(MPI_rank == 0)
      cerr << "optimizing mesh layout" << endl ;
    VOG::optimizeMesh(pos,cl,cr,face2node) ;
  }

  if(MPI_rank == 0)
    cerr << "writing VOG file" << endl ;
  Loci::writeVOG(outfile, pos, cl, cr, face2node,surf_ids) ;

  Loci::Finalize() ;
}
