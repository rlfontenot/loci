#include <Loci.h>
//#include "read_grid.h"
//#include "sciTypes.h"
#include <iostream>
#include <fstream>

#include <utility>
#include <algorithm>
#include <string>
#include <vector>
#include <list>
#include <string>

#include <rpc/rpc.h>
#include <rpc/xdr.h>

using std::string ;

//using namespace std ;
using std::ifstream ;
using std::cout ;
using std::cin ;
using std::endl ;
using std::cerr ;
using std::vector ;
using std::list ;
using std::ios ;

//using chem::vect3d ;

typedef vector3d<double> vect3d ;


int twoDExtract(int ac, char *av[]) {
  bool found_flag = false ;
  vector<int> bc ;
  enum coord_view { XR,XY,XZ,YZ } ;
  coord_view view_type = XY ;
  
  do {
    found_flag = true ;
    if(ac > 2 && !strcmp(av[1],"-bc")) {
      int ibc = -atoi(av[2]) ;
      bc.push_back(ibc) ;
      av+=2 ;
      ac-=2 ;
    } else if(ac > 1 && !strcmp(av[1],"-xr")) {
      view_type = XR ;
      av++ ;
      ac-- ;
    } else if(ac > 1 && !strcmp(av[1],"-xy")) {
      view_type = XY ;
      av++ ;
      ac-- ;
    } else if(ac > 1 && !strcmp(av[1],"-xz")) {
      view_type = XZ ;
      av++ ;
      ac-- ;
    } else if(ac > 1 && !strcmp(av[1],"-yz")) {
      view_type = YZ ;
      av++ ;
      ac-- ;
    } else {
      found_flag = false ;
    }
  } while(found_flag) ;
  
  char *problem_name = av[1] ;
  ac-- ;
  av++ ;
  bool data = false ;
  char bstr[] = "0" ;
  char *ncyc = bstr ;
  int type = 0 ;
  int number_type= 0 ;
  char *varfile = 0 ;
  bool vectorval = false ;
  if(ac >= 2) {
    data = true ;
    ncyc = av[1] ;
  }
  if(ac == 3) {
    type = *av[2] ;
    if(type >= '0' && type <= '9') {
      number_type = atoi(av[2]) ;
    }
    if(type == 'n') {
      varfile = av[2]+1 ;
    }
    if(type == 'v') {
        varfile = av[2]+1 ;
        vectorval = true ;
    }
  }

  fact_db facts ;
  string gridfile = string(problem_name) + string(".xdr") ;
  if(!Loci::readFVMGrid(facts,gridfile)) {
    cerr << "Unable to read file '" << gridfile << "'" << endl
         << "unable to continue." << endl ;
    return -1 ;
  }

  Map cr ;
  cr = facts.get_fact("cr") ;
  Map cl ;
  cl = facts.get_fact("cl") ;
  multiMap face2node ;
  face2node = facts.get_fact("face2node") ;
  entitySet img = cr.image(cr.domain()) ;
  img &= interval(Loci::UNIVERSE_MIN,-1) ;

  if(bc.size() == 0)
    bc.push_back(-1) ;
  
  entitySet negative_numbers  ;
  for(size_t i=0;i<bc.size();++i)
    negative_numbers += bc[i] ;

  entitySet boundary_faces = cr.preimage(negative_numbers).first ;
  cout << "boundary_faces = " << boundary_faces.size() << endl ;

  store<vect3d> pos ;
  pos = facts.get_fact("pos") ;

  
  multiMap node2face ;
  Loci::inverseMap(node2face,face2node,pos.domain(),face2node.domain()) ;

  cerr << "node2face.domain()" << node2face.domain() << endl ;
  dMap El,Er,N1,N2 ;
  int edge_count = 0 ;
  int count2 = 0 ;
  store<bool> visited ;
  visited.allocate(face2node.domain()) ;
  FORALL(face2node.domain(),i) {
    visited[i] = false ;
  } ENDFORALL ;
  FORALL(boundary_faces,fc) {
    visited[fc] = true ;
    int nf = face2node.end(fc) - face2node.begin(fc) ;
    for(int i=0;i<nf;++i) {
      int n1 = face2node[fc][(i==0)?(nf-1):(i-1)] ;
      int n2 = face2node[fc][i] ;
      bool found_match = false ;
      for(Entity *e = node2face.begin(n1);e!=node2face.end(n1);++e) {
        if(cr[*e] >= 0 || *e == fc) 
          continue ; // not a boundary face, continue
        int nf2 = face2node.end(*e)-face2node.begin(*e) ;
        int k ;
        for(k=0;k<nf2;++k)
          if(face2node[*e][k] == n1)
            break ;
        if(k==nf2) {
          cerr << "didn't find n1!" << endl ;
        }
        if(face2node[*e][(k==0)?(nf2-1):(k-1)] == n2 ||
           face2node[*e][(k+1==nf2)?0:(k+1)] == n2) {
          // found match
          if(found_match) {
            cerr << "found 2 matches for face " << endl ;
          }
          if(!visited[*e]) { // Only add an edge if it hasn't been seen before
            El[edge_count] = fc ;
            Er[edge_count] = *e ;
            if(!negative_numbers.inSet(cr[*e])) 
              Er[edge_count] = cr[*e] ;
            N1[edge_count] = n1 ;
            N2[edge_count] = n2 ;
            edge_count++ ;
          } else
            count2++ ;
          found_match = true ;
          break ;
        }
      }
      if(!found_match) {
        cerr << "didn't find match for face" << endl ;
      }
    }
  } ENDFORALL ;

  entitySet subnodes = N1.image(N1.domain()) + N2.image(N2.domain()) ;
  cout << "nodes count = " << subnodes.size() << endl ;
  cout << "edge_count = " << edge_count << endl ;
  store<int> node_numbers ;
  node_numbers.allocate(subnodes) ;
  int cnt = 1 ;
  FORALL(subnodes,nd) {
    node_numbers[nd] = cnt ;
    cnt++ ;
  } ENDFORALL ;
  store<int> cell_numbers ;
  cell_numbers.allocate(boundary_faces) ;
  cnt = 1 ;
  FORALL(boundary_faces,cl) {
    cell_numbers[cl] = cnt ;
    cnt++ ;
  } ENDFORALL ;

  string *species=0 ;
  int squery=-1 ;
  char filename[512] ;

  store<float> value ;
  entitySet valdom = pos.domain() ;
  value.allocate(valdom) ;
  for(entitySet::const_iterator ei=valdom.begin();ei!=valdom.end();++ei)
    value[*ei] = 0 ;
  
  if(data) {
    entitySet::const_iterator ptr = valdom.begin() ;

    if(varfile == 0) {
      sprintf(filename,"output/qn_xdr.%s",ncyc) ;
      FILE *FP = fopen(filename, "r") ;
      if(FP == NULL) {
        cerr << "open failed on " << filename << endl ;
        exit(EXIT_FAILURE) ;
      }
      XDR xdr_handle ;
      xdrstdio_create(&xdr_handle, FP, XDR_DECODE) ;
      int ns ;
      char tmp_char[256] ;
      int tmp_size = 0 ;
      xdr_int(&xdr_handle, &ns) ;
      species = new string[ns];
      int numvars = 7 + ns ;
      int nv = numvars ;
      for(int i=0;i<ns;++i) {
        xdr_int(&xdr_handle, &tmp_size) ;
        for(int j = 0; j < tmp_size; ++j)
          xdr_char(&xdr_handle, &tmp_char[j]) ;
        for(int k = 0; k < tmp_size; ++k)
          species[i] += tmp_char[k] ;
      }
      int nnodes  =0;
      xdr_int(&xdr_handle, &nnodes) ;
      
      if(type == 'f') {
        string s = string(av[2]+1) ;
        for(int i=0;i<ns;++i)
          if(species[i] == s)
            squery = i ;
        if(squery == -1) {
          cerr << "species " << s << " not in data " << endl ;
          exit(EXIT_FAILURE) ;
        }
      }
      
      double *qb = new double[nv] ;
      for(int i=0;i<nnodes;++i) {
        for(int ii=0;ii<nv;++ii) 
          xdr_double(&xdr_handle, &qb[ii]) ;
        double r =  qb[0] ;
        double u =  qb[1] ;
        double v =  qb[2] ;
        double w =  qb[3] ;
        double a =  qb[4] ;
        double T =  qb[5] ;
        double P =  qb[6] ;
        double U = sqrt(u*u+v*v+w*w) ;
        double M = U/a ;
        double Gamma = a*a*r/P ;
        double Rt = P/(r*T) ;
        int t ;
        double val = 0 ;
        switch(type) {
        case 'r':
          val = r ;
          break ;
        case 'P':
          val = P ;
          break ;
        case 'p':
          val = log10(P)  ;
          break ;
        case 'u':
          val = U ;
          break ;
        case 'm':
          val = M ;
          break ;
        case 't':
          val = T ;
          break ;
        case 'a':
          val = a ;
          break ;
        case 'g':
          val = Gamma ;
          break ;
        case 'R':
          val = Rt ;
          break ;
        case '0':
          val = u ;
          break ;
        case '1':
          val = v ;
          break ;
        case '2':
          val = w ;
          break ;
        case 'f':
          t = squery+7 ;
          val =  qb[t] ;
          break ;
        default:
          cerr << "wrong type"<< endl ;
          exit(EXIT_FAILURE) ;
        }
        if(ptr == valdom.end()) {
          cerr << "node_domain smaller than nodal data file!" << endl ;
          break ;
        }
        value[*ptr]  = val ;
        ptr++ ;
        
      }
    } else {
      char tmp_name[512] ;
      sprintf(tmp_name, "output/%s_hdf5.%s",varfile, ncyc) ;
      cout << "opening file " << tmp_name << endl ;
      hid_t file_id, group_id ;
      file_id = H5Fopen(tmp_name,H5F_ACC_RDONLY, H5P_DEFAULT) ;
      group_id = H5Gopen(file_id,varfile) ;
      store<float> var ;
      entitySet tmpdom = EMPTY ;
      Loci::read_container(group_id, var.Rep(), tmpdom) ;
      H5Gclose(group_id) ;
      H5Fclose(file_id) ;
      if(tmpdom.size() != valdom.size()) {
        cerr << " domains don't match in size for variable " << tmp_name
             << endl ;
      } else {
        for(entitySet::const_iterator ei=valdom.begin(),ti=tmpdom.begin();
            ei!=valdom.end();++ei,++ti) {
          value[*ei] = var[*ti] ;
        }
      }
    }
  }
  char out_name[512] ;
  if(!data) {
    sprintf(out_name,"%s.2dgv",problem_name) ;
  } else if(varfile != 0) {
    sprintf(out_name,"%s.%s",varfile,ncyc) ;
  } else if((type >='0' && type <='9'))
    sprintf(out_name,"%d.%s",number_type,ncyc) ;
  else if(type != 'f')
    sprintf(out_name,"%c.%s",type,ncyc) ;
  else
    sprintf(out_name,"%s.%s",species[squery].c_str(),ncyc) ;
  
  std::ofstream out(out_name,ios::out) ;
  out.precision(16) ;
  out << "general" << endl ;
  out << subnodes.size() << ' ' << '1' << endl ;
  FORALL(subnodes,nd) {
    double x = pos[nd].x ;
    double y = pos[nd].y ;
    double z = pos[nd].z ;

    double r = sqrt(y*y+z*z) ;
    switch(view_type) {
    case XR:
      out << x << ' ' << r << endl ;
      break ;
    case XY:
      out << x << ' ' << y << endl ;
      break ;
    case XZ:
      out << x << ' ' << z << endl ;
      break ;
    case YZ:
      out << y << ' ' << z << endl ;
      break ;
    default:
      cerr << "wrong output coordinate system!" << endl ;
      exit(EXIT_FAILURE) ;
    }
  } ENDFORALL ;
  out << edge_count << " 1 " << cnt << " 1" << endl ;
  
  for(int i=0;i<edge_count;++i) {
    out << node_numbers[N1[i]] << ' '
        << node_numbers[N2[i]] << ' '
        << cell_numbers[El[i]] << ' ' ;
    if(Er[i] < 0)
      out << Er[i] << endl ;
    else
      out << cell_numbers[Er[i]] << endl ;
  }
  
  
  FORALL(subnodes,nd) {
    out << value[nd] << endl ;
  } ENDFORALL ;

  return 0 ;
}

int asciiExtract(int ac, char *av[]) {
  char *problem_name = av[1] ;
  ac-- ;
  av++ ;
  char bstr[] = "0" ;
  char *ncyc = bstr ;
  int type = 0 ;
  //  int number_type= 0 ;
  if(ac >= 2)
    ncyc = av[1] ;

  fact_db facts ;
  string gridfile = string(problem_name) + string(".xdr") ;
  if(!Loci::readFVMGrid(facts,gridfile)) {
    cerr << "Unable to read file '" << gridfile << "'" << endl
         << "unable to continue." << endl ;
    exit(EXIT_FAILURE) ;
  }


  store<vect3d> pos ;
  pos = facts.get_fact("pos") ;

  

  string *species=0 ;
  int squery=-1 ;
  char filename[512] ;

  list<Loci::storeRepP> value_list ;

  while(ac >= 3) {
    store<float> value ;
    entitySet valdom = pos.domain() ;
    value.allocate(valdom) ;
    entitySet::const_iterator ptr = valdom.begin() ;

    char *varfile = 0 ;
    bool vectorvar = false ;
    char *speciesname =0;
    type = *av[2] ;
    //    if(type >= '0' && type <= '9') {
    //      number_type = atoi(av[2]) ;
    //    }
    if(type == 'n') {
      varfile = av[2]+1 ;
    }
    if(type == 'v') {
      varfile = av[2]+1 ;
      vectorvar = true ;
    }
    if(type == 'f') {
      speciesname = av[2]+1 ;
    }
    ac-- ;
    av++ ;
    if(type == 'x' || type == 'y' || type == 'z') {
      for(entitySet::const_iterator ii=valdom.begin();ii != valdom.end();++ii) {
        if(type == 'x')
          value[*ii] = pos[*ii].x ;
        if(type == 'y')
          value[*ii] = pos[*ii].y ;
        if(type == 'z')
          value[*ii] = pos[*ii].z ;
      }
    } else if(varfile == 0) {
      sprintf(filename,"output/qn_xdr.%s",ncyc) ;
      FILE *FP = fopen(filename, "r") ;
      if(FP == NULL) {
        cerr << "open failed on " << filename << endl ;
        exit(EXIT_FAILURE) ;
      }
      XDR xdr_handle ;
      xdrstdio_create(&xdr_handle, FP, XDR_DECODE) ;
      int ns ;
      char tmp_char[256] ;
      int tmp_size = 0 ;
      xdr_int(&xdr_handle, &ns) ;
    
      species = new string[ns];
      int numvars = 7 + ns ;
      int nv = numvars ;
      for(int i=0;i<ns;++i) {
        xdr_int(&xdr_handle, &tmp_size) ;
        for(int j = 0; j < tmp_size; ++j)
          xdr_char(&xdr_handle, &tmp_char[j]) ;
        for(int k = 0; k < tmp_size; ++k)
          species[i] += tmp_char[k] ;
      }
      int nnodes  =0;
      xdr_int(&xdr_handle, &nnodes) ;
    
      if(type == 'f') {
        string s = string(speciesname) ;
        for(int i=0;i<ns;++i)
          if(species[i] == s)
            squery = i ;
        if(squery == -1) {
          cerr << "species " << s << " not in data " << endl ;
          exit(EXIT_FAILURE) ;
        }
      }
    
      double *qb = new double[nv] ;
      for(int i=0;i<nnodes;++i) {
        for(int ii=0;ii<nv;++ii) 
          xdr_double(&xdr_handle, &qb[ii]) ;
        double r =  qb[0] ;
        double u =  qb[1] ;
        double v =  qb[2] ;
        double w =  qb[3] ;
        double a =  qb[4] ;
        double T =  qb[5] ;
        double P =  qb[6] ;
        double U = sqrt(u*u+v*v+w*w) ;
        double M = U/a ;
        double Gamma = a*a*r/P ;
        double Rt = P/(r*T) ;
        int t ;
        double val = 0 ;
        switch(type) {
        case 'r':
          val = r ;
          break ;
        case 'P':
          val = P ;
          break ;
        case 'p':
          val = log10(P)  ;
          break ;
        case 'u':
          val = U ;
          break ;
        case 'm':
          val = M ;
          break ;
        case 't':
          val = T ;
          break ;
        case 'a':
          val = a ;
          break ;
        case 'g':
          val = Gamma ;
          break ;
        case 'R':
          val = Rt ;
          break ;
        case '0':
          val = u ;
          break ;
        case '1':
          val = v ;
          break ;
        case '2':
          val = w ;
          break ;
        case 'f':
          t = squery+7 ;
          val =  qb[t] ;
          break ;
        default:
          cerr << "wrong type"<< endl ;
          exit(EXIT_FAILURE) ;
        }
        if(ptr == valdom.end()) {
          cerr << "node_domain smaller than nodal data file!" << endl ;
          break ;
        }
        value[*ptr]  = val ;
        ptr++ ;
      
      }
    } else {
        if(vectorvar) {
            char tmp_name[512] ;
            sprintf(tmp_name, "output/%s_hdf5.%s",varfile, ncyc) ;
            hid_t file_id, group_id ;
            file_id = H5Fopen(tmp_name,H5F_ACC_RDONLY, H5P_DEFAULT) ;
            group_id = H5Gopen(file_id,varfile) ;
            store<vector3d<float> > var ;
            cerr << "reading " << tmp_name << endl ;
            Loci::read_container(group_id, var.Rep(), valdom) ;
            H5Gclose(group_id) ;
            H5Fclose(file_id) ;


            store<float> valuex,valuey ;
            valuex.allocate(valdom) ;
            valuey.allocate(valdom) ;
            for(entitySet::const_iterator ei=valdom.begin();
                ei!=valdom.end();++ei) {
                valuex[*ei] = var[*ei].x ;
                valuey[*ei] = var[*ei].y ;
                value[*ei] = var[*ei].z ;
            }
            value_list.push_back(valuex.Rep()) ;
            value_list.push_back(valuey.Rep()) ;
            
        } else {
            char tmp_name[512] ;
            sprintf(tmp_name, "output/%s_hdf5.%s",varfile, ncyc) ;
            hid_t file_id, group_id ;
            file_id = H5Fopen(tmp_name,H5F_ACC_RDONLY, H5P_DEFAULT) ;
            group_id = H5Gopen(file_id,varfile) ;
            store<double> var ;
            Loci::read_container(group_id, var.Rep(), valdom) ;
            H5Gclose(group_id) ;
            H5Fclose(file_id) ;
            for(entitySet::const_iterator ei=valdom.begin();
                ei!=valdom.end();++ei) {
                value[*ei] = var[*ei] ;
            }
        }
    }
    value_list.push_back(value.Rep()) ;
  }
  entitySet dom = pos.domain() ;
  for(entitySet::const_iterator ei=dom.begin();ei!=dom.end();++ei) {
    list<Loci::storeRepP>::const_iterator li ;
    for(li=value_list.begin();li!=value_list.end();++li) {
      store<float> v ;
      v = *li ;
      cout << v[*ei] << ' ' ;
    }
    cout << endl ;
  }
  return 0 ;
}
