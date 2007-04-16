#include <Loci.h> 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string>
using std::string ;
#include <iostream> 
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>
using std::vector ;
using std::string ;
using std::cerr ;
using std::endl ;
using std::cout ;
using std::map ;
using std::ofstream ;
using std::ios ;
using std::sort ;
using std::unique ;

#include <sys/types.h>
#include <sys/stat.h>

#include "extract.h"

    

void Usage(int ac, char *av[]) {
  cerr << av[0] << ": Incorrect Usage" << endl ;
  cout << "Usage:" << endl;
  cout << av[0] << " <package> [package options] <case_name> <time_step> <variable(s)>" << endl ;
  cout << endl ;
  cout << "where <package> may be:" << endl
       << "-2d :  extract for the 2dgv plotting package" << endl
       << "-fv :  extract for the FieldView post-processing package" << endl
       << "-en :  extract for the Ensight post-processing package" << endl
       << "-tec:  extract for the TecPlot post-procesing package" << endl
       << "-ascii: extract to an ascii file" << endl
       << endl ;
  cout << "Variables are defined by the solver, but typically include: " << endl
       << "r     - nodal density" << endl 
       << "p     - nodal log10 pressure" << endl
       << "P     - nodal absolute pressure" << endl
       << "pg    - nodal gage pressure" << endl 
       << "u     - nodal velocity magnitude" << endl
       << "m     - nodal mach number" << endl
       << "t     - nodal temperature" << endl
       << "a     - nodal soundspeed" << endl
    //       << "G - extract Gamma" << endl
    //       << "R - extract Gas Constant R~" << endl
       << "f<species> - nodal mass fractions for <species>" << endl
       << "v     - nodal velocity vector"<<endl
       << "x     - x coordinate" << endl 
       << "y     - y coordinate" << endl 
       << "z     - z coordinate" << endl
       << "Boundary Variables:" << endl 
       << "qdot  - wall boundary heat flux" << endl
       << "yplus - wall boundary y plus" << endl
       << "tau   - viscous wall shear stress vector" << endl
       << "tw    - viscous wall boundary temperature" << endl
       << "pw    - viscous wall boundary pressure" << endl
       << "n     - boundary normal (-ascii only)" << endl
       << "area  - boundary area   (-ascii only)" << endl
       << endl ;

  cout << "extra options for the 2dgv postprocessing package" << endl 
       << "  -bc <boundary_tag>  : specify which boundary to extract for (multiple -bc's"
       << endl 
       << "             are merged)" << endl 
       << "  -xy : project boundary onto z=0 plane (default)" << endl
       << "  -yz : project boundary onto x=0 plane" << endl
       << "  -xz : project boundary onto y=0 plane" << endl
       << "  -xr : project boundary onto x,radius plane (for axisymmetric grids)"
       << endl << endl ;
  cout << "example:  to extract OH species from time step 50 of ssme simulation for visualization with 2dgv use:" << endl
       << av[0] << " -2d -bc 1 -xr ssme 50 fOH" << endl ;
  cout << "example: to extract an ascii table of boundary heat flux and x locations:"<<endl
       << av[0] << " -ascii -bc 4 nozzle 0 x qdot" << endl ;

  exit(-1) ;
}

int  sizeElementType(hid_t group_id, const char *element_name) {
  hid_t dataset = H5Dopen(group_id,element_name) ;
  if(dataset < 0) {
    H5Eclear() ;
    return 0 ;
  }
  hid_t dspace = H5Dget_space(dataset) ;

  hsize_t size = 0 ;
  H5Sget_simple_extent_dims(dspace,&size,NULL) ;
  
  
  H5Dclose(dataset) ;
  return int(size) ;
  
}


void getDerivedVar(vector<float> &dval, string var_name,
                   string casename, string iteration) {
  if(var_name == "m") {
    string filename = "output/a_sca."+iteration + "_" + casename ;
    
    hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                       H5F_ACC_RDONLY,
                                       H5P_DEFAULT) ;
    if(file_id < 0) {
      cerr << "unable to open file '" << filename << "'!" << endl ;
      return ;
    }

    fact_db facts ;
    store<float> soundSpeed ;
    Loci::readContainer(file_id,"a",soundSpeed.Rep(),EMPTY,facts) ;
    Loci::hdf5CloseFile(file_id) ;

    filename = "output/v_vec." + iteration +"_" + casename ;
    file_id = Loci::hdf5OpenFile(filename.c_str(),
                                       H5F_ACC_RDONLY,
                                       H5P_DEFAULT) ;
    if(file_id < 0) {
      cerr << "unable to open file '" << filename << "'!" << endl ;
      return ;
    }

    store<vector3d<float> > u ;
    Loci::readContainer(file_id,"v",u.Rep(),EMPTY,facts) ;
    Loci::hdf5CloseFile(file_id) ;

    entitySet dom = u.domain() ;
    int c = 0 ;
    FORALL(dom,nd) {
      float m = norm(u[nd])/soundSpeed[nd] ;
      dval[c++] = m ;
    } ENDFORALL ;
  } else if(var_name == "p" || var_name == "P") {
    string filename = "output/pg_sca."+iteration + "_" + casename ;
    
    hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                       H5F_ACC_RDONLY,
                                       H5P_DEFAULT) ;
    if(file_id < 0) {
      cerr << "unable to open file '" << filename << "'!" << endl ;
      return ;
    }

    fact_db facts ;
    store<float> pg ;
    Loci::readContainer(file_id,"pg",pg.Rep(),EMPTY,facts) ;
    Loci::hdf5CloseFile(file_id) ;

    filename = "output/Pambient_par." + iteration +"_" + casename ;
    file_id = Loci::hdf5OpenFile(filename.c_str(),
                                       H5F_ACC_RDONLY,
                                       H5P_DEFAULT) ;
    if(file_id < 0) {
      cerr << "unable to open file '" << filename << "'!" << endl ;
      return ;
    }

    param<float> Pambient ;
    Loci::readContainer(file_id,"Pambient",Pambient.Rep(),EMPTY,facts) ;
    Loci::hdf5CloseFile(file_id) ;

    entitySet dom = pg.domain() ;
    bool log = (var_name == "p") ;
    int c = 0 ;
    FORALL(dom,nd) {
      float p = pg[nd]+Pambient[nd] ;
      if(log)
        dval[c++] = log10(p) ;
      else
        dval[c++] = p ;
    } ENDFORALL ;
  } else if (var_name == "u") {
    fact_db facts ;
    string filename = "output/v_vec." + iteration +"_" + casename ;
    hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                       H5F_ACC_RDONLY,
                                       H5P_DEFAULT) ;
    if(file_id < 0) {
      cerr << "unable to open file '" << filename << "'!" << endl ;
      return ;
    }

    store<vector3d<float> > u ;
    Loci::readContainer(file_id,"v",u.Rep(),EMPTY,facts) ;
    Loci::hdf5CloseFile(file_id) ;

    entitySet dom = u.domain() ;
    int c = 0 ;
    FORALL(dom,nd) {
      float m = norm(u[nd]) ;
      dval[c++] = m ;
    } ENDFORALL ;
  } else if(var_name == "x" || var_name =="y" || var_name == "z") {
    store<vector3d<float> > pos ;
    string posname = "output/grid_pos." + iteration + "_" + casename ;
    hid_t file_id = Loci::hdf5OpenFile(posname.c_str(),
                                       H5F_ACC_RDONLY,
                                       H5P_DEFAULT) ;
    if(file_id < 0) {
      cerr << "unable to get grid positions for iteration " << iteration
           << endl ;
      cerr << "does file '" << posname << "' exist?" << endl ;
      exit(-1) ;
    }

    fact_db facts ;
    Loci::readContainer(file_id,"pos",pos.Rep(),EMPTY,facts) ;
    Loci::hdf5CloseFile(file_id) ;
    entitySet dom = pos.domain() ;
    int c = 0 ;
    if(var_name == "x") {
      FORALL(dom,nd) {
        dval[c++] = pos[nd].x ;
      } ENDFORALL ;
    }
    if(var_name == "y") {
      FORALL(dom,nd) {
        dval[c++] = pos[nd].x ;
      } ENDFORALL ;
    }
    if(var_name == "z") {
      FORALL(dom,nd) {
        dval[c++] = pos[nd].x ;
      } ENDFORALL ;
    }
  } else if(var_name == "0" || var_name =="1" || var_name == "2") {
    fact_db facts ;
    string filename = "output/v_vec." + iteration +"_" + casename ;
    hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                       H5F_ACC_RDONLY,
                                       H5P_DEFAULT) ;
    if(file_id < 0) {
      cerr << "unable to open file '" << filename << "'!" << endl ;
      return ;
    }

    store<vector3d<float> > u ;
    Loci::readContainer(file_id,"v",u.Rep(),EMPTY,facts) ;
    Loci::hdf5CloseFile(file_id) ;

    entitySet dom = u.domain() ;
    int c = 0 ;
    if(var_name == "0") {
      FORALL(dom,nd) {
        dval[c++] = u[nd].x ;
      } ENDFORALL ;
    }
    if(var_name == "1") {
      FORALL(dom,nd) {
        dval[c++] = u[nd].y ;
      } ENDFORALL ;
    }
    if(var_name == "2") {
      FORALL(dom,nd) {
        dval[c++] = u[nd].z ;
      } ENDFORALL ;
    }
  } else {
    cerr << "don't know how to get derived variable " << var_name << endl ;
  }
}

void setup_grid_topology(string casename, string iteration) {
  fact_db facts ;
  string file = casename + ".xdr" ;
  if(!Loci::setupFVMGrid(facts,file)) {
    cerr << "unable to read grid " << file << endl ;
  }
  createLowerUpper(facts) ;
  string filename = "output/"+casename+".topo" ;
  multiMap upper,lower,boundary_map,face2node ;
  Map ref ;
  store<string> boundary_names ;
  constraint geom_cells ;
  upper = facts.get_variable("upper") ;
  lower = facts.get_variable("lower") ;
  boundary_map = facts.get_variable("boundary_map") ;
  face2node = facts.get_variable("face2node") ;
  ref = facts.get_variable("ref") ;
  boundary_names = facts.get_variable("boundary_names") ;
  geom_cells = facts.get_variable("geom_cells") ;
  store<vector3d<double> > pos ;
  pos = facts.get_variable("pos") ;
  
  Loci::parallelWriteGridTopology(filename.c_str(),
                                  upper.Rep(),lower.Rep(),boundary_map.Rep(),
                                  face2node.Rep(),
                                  ref.Rep(),
                                  boundary_names.Rep(),
                                  pos.Rep(),
                                  *geom_cells,
                                  facts) ;

  filename = "output/grid_pos." + iteration + "_" + casename ;
  hid_t file_id = Loci::hdf5CreateFile(filename.c_str(),H5F_ACC_TRUNC,
                                       H5P_DEFAULT, H5P_DEFAULT) ;
  
  Loci::writeContainer(file_id,"pos",pos.Rep(),facts) ;
  
  Loci::hdf5CloseFile(file_id) ;

}


void extract_grid(string casename, string iteration,
                  grid_topo_handler *topo,
                  vector<string> variables,
                  vector<int> variable_types,
                  vector<string> variable_filenames) {
  
  vector<string> bnd_scalar_vars,bnd_scalar_filenames ;
  vector<string> bnd_vector_vars,bnd_vector_filenames ;
  vector<string> mfvars ;

  {
    vector<string> vname ;
    vector<int> vtype;
    vector<string> vfile ;
  
    for(size_t i=0;i<variables.size();++i) {
      const string var_name(variables[i]) ;
      const int vt = variable_types[i] ;
      const string filename(variable_filenames[i]) ;
      switch(vt) {
      case NODAL_SCALAR:
      case NODAL_DERIVED:
      case NODAL_VECTOR:
        vname.push_back(var_name) ;
        vfile.push_back(filename) ;
        vtype.push_back(vt) ;
        break ;
      case NODAL_MASSFRACTION:
        mfvars.push_back(var_name) ;
        break ;
      case BOUNDARY_SCALAR:
        bnd_scalar_vars.push_back(var_name) ;
        bnd_scalar_filenames.push_back(filename) ;
        break ;
      case BOUNDARY_VECTOR:
        bnd_vector_vars.push_back(var_name) ;
        bnd_vector_filenames.push_back(filename) ;
        break ;
      default:
        cerr << "unable to process variable " << var_name << endl ;
        break ;
      }
    }
    for(size_t i=0;i<mfvars.size();++i) {
      vname.push_back(mfvars[i]) ;
      vfile.push_back(string("")) ;
      vtype.push_back(NODAL_MASSFRACTION) ;
    }
    for(size_t i=0;i<bnd_scalar_vars.size();++i) {
      vname.push_back(bnd_scalar_vars[i]) ;
      vfile.push_back(bnd_scalar_filenames[i]) ;
      vtype.push_back(BOUNDARY_SCALAR) ;
    }
    for(size_t i=0;i<bnd_vector_vars.size();++i) {
      vname.push_back(bnd_vector_vars[i]) ;
      vfile.push_back(bnd_vector_filenames[i]) ;
      vtype.push_back(BOUNDARY_VECTOR) ;
    }

    variables = vname ;
    variable_filenames = vfile ;
    variable_types = vtype ;
  }

  Array<int,5> events ;
  topo->fileWritingSequence(events) ;
  FATAL(Loci::MPI_processes != 1) ;
  store<vector3d<float> > pos ;
  string posname = "output/grid_pos." + iteration + "_" + casename ;
  hid_t file_id = Loci::hdf5OpenFile(posname.c_str(),
                                     H5F_ACC_RDONLY,
                                     H5P_DEFAULT) ;
  if(file_id < 0) {
    cerr << "unable to get grid positions for iteration " << iteration
         << endl ;
    cerr << "does file '" << posname << "' exist?" << endl ;
    exit(-1) ;
  }

  fact_db facts ;
  Loci::readContainer(file_id,"pos",pos.Rep(),EMPTY,facts) ;
  Loci::hdf5CloseFile(file_id) ;
  int npnts = pos.domain().size() ;


  Loci::hdf5CloseFile(file_id) ;
  
  string gridtopo = "output/" + casename +".topo" ;


  file_id = H5Fopen(gridtopo.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT) ;

  hid_t elg = H5Gopen(file_id,"elements") ;

  int ntets = sizeElementType(elg,"tetrahedra") ;
  int nhexs = sizeElementType(elg,"hexahedra") ;
  int nprsm = sizeElementType(elg,"prism") ;
  int npyrm = sizeElementType(elg,"pyramid") ;
  int ngenc = sizeElementType(elg,"GeneralCellNfaces") ;

  hid_t bndg = H5Gopen(file_id,"boundaries") ;
  hsize_t num_bcs = 0 ;
  H5Gget_num_objs(bndg,&num_bcs) ;
  vector<string>  bc_names ;
  for(hsize_t bc=0;bc<num_bcs;++bc) {
    char buf[1024] ;
    memset(buf, '\0', 1024) ;
    H5Gget_objname_by_idx(bndg,bc,buf,sizeof(buf)) ;
    buf[1023]='\0' ;
    bc_names.push_back(string(buf)) ;
  }
  
  cout << "npnts = " << npnts << ' '
       << "ntets = " << ntets << ' '
       << "npyrm = " << npyrm << ' '
       << "nprsm = " << nprsm << ' '
       << "nhexs = " << nhexs << ' '
       << "ngenc = " << ngenc << endl ;

  cout << "bcs = " ;
  for(size_t i=0;i<bc_names.size();++i)
    cout << bc_names[i] << ' ' ;
  cout << endl ;
  
  topo->open(casename,iteration,npnts,ntets,nprsm,npyrm,nhexs,ngenc,
             bc_names, variables,variable_types) ;


  for(int i=0;i<5;++i) {
    switch(events[i]) {
    case GRID_POSITIONS:
      {
        int minpos = pos.domain().Min() ;
        topo->create_mesh_positions(&pos[minpos],npnts) ;
        pos.allocate(EMPTY) ;
      }
      break ;
    case GRID_VOLUME_ELEMENTS:
      topo->create_mesh_elements() ;
      if(ntets > 0) {
        vector<Array<int,4> > tets(ntets) ;
        readElementType(elg,"tetrahedra",tets) ;
        topo->write_tets(&tets[0],ntets) ;
      }
      if(npyrm > 0) {
        vector<Array<int,5> > pyrm(npyrm) ;
        readElementType(elg,"pyramid",pyrm) ;
        topo->write_pyrm(&pyrm[0],npyrm) ;
      }
      if(nprsm > 0) {
        vector<Array<int,6> > prsm(nprsm) ;
        readElementType(elg,"prism",prsm) ;
        topo->write_prsm(&prsm[0],nprsm) ;
      }
      if(nhexs > 0) {
        vector<Array<int,8> > hexs(nhexs) ;
        readElementType(elg,"hexahedra",hexs) ;
        topo->write_hexs(&hexs[0],nhexs) ;
      }
      if(ngenc > 0) {
        vector<int> GeneralCellNfaces(ngenc) ;
        readElementType(elg,"GeneralCellNfaces",GeneralCellNfaces) ;
        int nside = sizeElementType(elg,"GeneralCellNsides") ;
        vector<int> GeneralCellNsides(nside) ;
        readElementType(elg,"GeneralCellNsides",GeneralCellNsides) ;
        int nnodes = sizeElementType(elg,"GeneralCellNodes") ;
        vector<int> GeneralCellNodes(nnodes) ;
        readElementType(elg,"GeneralCellNodes",GeneralCellNodes) ;
        topo->write_general_cell(&GeneralCellNfaces[0],ngenc,
                                 &GeneralCellNsides[0],nside,
                                 &GeneralCellNodes[0],nnodes) ;
      }
      topo->close_mesh_elements() ;
      H5Gclose(elg) ;
      break ;
    case GRID_BOUNDARY_ELEMENTS:
      for(hsize_t bc=0;bc<num_bcs;++bc) {
        hid_t bcg = H5Gopen(bndg,bc_names[bc].c_str()) ;
        
        int nquads = sizeElementType(bcg,"quads") ;
        int ntrias = sizeElementType(bcg,"triangles") ;
        int ngeneral = sizeElementType(bcg,"nside_sizes") ;
        
        vector<Array<int,3> > trias(ntrias) ;
        readElementType(bcg,"triangles",trias) ;
        vector<Array<int,4> > quads(nquads) ;
        readElementType(bcg,"quads",quads) ;
        vector<int> nside_sizes(ngeneral) ;
        readElementType(bcg,"nside_sizes",nside_sizes) ;
        int nside_nodes_size = sizeElementType(bcg,"nside_nodes") ;
        vector<int> nside_nodes(nside_nodes_size) ;
        readElementType(bcg,"nside_nodes",nside_nodes) ;
        
        vector<int> node_set ;
        for(int i=0;i<nside_nodes_size;++i)
          node_set.push_back(nside_nodes[i]) ;
        
        for(int i=0;i<ntrias;++i) {
          node_set.push_back(trias[i][0]) ;
          node_set.push_back(trias[i][1]) ;
          node_set.push_back(trias[i][2]) ;
        }
        for(int i=0;i<nquads;++i) {
          node_set.push_back(quads[i][0]) ;
          node_set.push_back(quads[i][1]) ;
          node_set.push_back(quads[i][2]) ;
          node_set.push_back(quads[i][3]) ;
        }
        sort(node_set.begin(),node_set.end()) ;
        node_set.erase(unique(node_set.begin(),node_set.end()),node_set.end()) ;
        
        map<int,int> nmap ;
        for(size_t i=0;i<node_set.size();++i) {
          nmap[node_set[i]] = i+1 ;
        }
        
        
        topo->create_boundary_part(bc_names[bc],&node_set[0],node_set.size()) ;
        
        for(int i=0;i<ntrias;++i) 
          for(int j=0;j<3;++j) 
            trias[i][j] = nmap[trias[i][j]] ;
        for(int i=0;i<nquads;++i) 
          for(int j=0;j<4;++j) 
            quads[i][j] = nmap[quads[i][j]] ;
        for(int i=0;i<nside_nodes_size;++i)
          nside_nodes[i] = nmap[nside_nodes[i]] ;
        
        vector<int > trias_id(ntrias) ;
        readElementType(bcg,"triangles_id",trias_id) ;
        vector<int > quads_id(nquads) ;
        readElementType(bcg,"quads_id",quads_id) ;
        vector<int > nside_id(ngeneral) ;
        readElementType(bcg,"nside_id",nside_id) ;
        
        topo->write_trias(&trias[0],&trias_id[0],ntrias) ;
        topo->write_quads(&quads[0],&quads_id[0],nquads) ;
        topo->write_general_face(&nside_sizes[0], &nside_id[0], ngeneral,
                                 &nside_nodes[0], nside_nodes_size) ;
        H5Gclose(bcg) ;
        topo->close_boundary_part() ;
        
      }
      break ;
    case NODAL_VARIABLES:
      {
        
        for(size_t i=0;i<variables.size();++i) {
          const string var_name(variables[i]) ;
          const int vt = variable_types[i] ;
          const string filename(variable_filenames[i]) ;
          switch(vt) {
          case NODAL_SCALAR:
            {
              hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                                 H5F_ACC_RDONLY,
                                                 H5P_DEFAULT) ;
              if(file_id < 0) {
                cerr << "unable to open file '" << filename << "'!" << endl ;
                continue ;
              }

              fact_db facts ;
              store<float> scalar ;
              Loci::readContainer(file_id,var_name,scalar.Rep(),EMPTY,facts) ;
              entitySet dom = scalar.domain() ;
              int sz = dom.size() ;
              int min_val= dom.Min() ;
              topo->output_nodal_scalar(&scalar[min_val],sz,var_name) ;
              Loci::hdf5CloseFile(file_id) ;
            }
            break;
          case NODAL_DERIVED:
            {
              vector<float> dval(npnts) ;
              getDerivedVar(dval,var_name,casename,iteration) ;
              int sz = dval.size() ;
              topo->output_nodal_scalar(&dval[0],sz,var_name) ;
            }
            break;
          case NODAL_VECTOR:
            {
              hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                                 H5F_ACC_RDONLY,
                                                 H5P_DEFAULT) ;
              if(file_id < 0) {
                cerr << "unable to open file '" << filename << "'!" << endl ;
                continue ;
              }
              
              fact_db facts ;
              store<vector3d<float> > vec ;
              Loci::readContainer(file_id,var_name,vec.Rep(),EMPTY,facts) ;
              entitySet dom = vec.domain() ;
              int sz = dom.size() ;
              int min_val= dom.Min() ;
              topo->output_nodal_vector(&vec[min_val],sz,var_name) ;
              Loci::hdf5CloseFile(file_id) ;
            }
            break;
          case NODAL_MASSFRACTION:
            break ;
          case BOUNDARY_SCALAR:
          case BOUNDARY_VECTOR:
            break ;
          default:
            cerr << "unable to process variable " << var_name << endl ;
            break ;
          }
        }

        if(mfvars.size() > 0) {
          string filename = "output/mix." + iteration + "_" + casename ;
          
          hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                             H5F_ACC_RDONLY,
                                             H5P_DEFAULT) ;
          if(file_id < 0) {
            cerr << "unable to open file '" << filename << "'!" << endl ;
            exit(-1) ;
          }
          
          fact_db facts ;
          storeVec<float> mix ;
          Loci::readContainer(file_id,"mixture",mix.Rep(),EMPTY,facts) ;
          param<string> species_names ;
          Loci::readContainer(file_id,"species_names",species_names.Rep(),EMPTY,facts) ;
          Loci::hdf5CloseFile(file_id) ;
          
          map<string,int> smap ;
          std::istringstream iss(*species_names) ;
          for(int i=0;i<mix.vecSize();++i) {
            string s ;
            iss >> s ;
            smap[s] = i ;
          }
      
          entitySet dom = mix.domain() ;
          
          vector<float> vec(npnts) ;
          
          for(size_t i=0;i<mfvars.size();++i) {
            const string var_name(mfvars[i]) ;
            string sp = string(mfvars[i].c_str()+1) ;
            map<string,int>::const_iterator mi = smap.find(sp) ;
            if(mi == smap.end()) {
              cerr << "warning, species " << sp << " does not exist in dataset!"
                   << endl ;
            } else {
              const int ind = mi->second ;
              int c = 0 ;
              FORALL(dom,nd) {
                vec[c++] = mix[nd][ind] ;
              } ENDFORALL ;
              topo->output_nodal_scalar(&vec[0],npnts,var_name) ;
            }
          }

        }
      }
      break;
    case BOUNDARY_VARIABLES:
      if(bnd_scalar_vars.size() > 0) {
        for(size_t b=0;b<bnd_scalar_vars.size();++b) {
          const string filename(bnd_scalar_filenames[b]) ;
          hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                             H5F_ACC_RDONLY,
                                             H5P_DEFAULT) ;
          if(file_id < 0) {
            cerr << "unable to open file '" << filename << "'!" << endl ;
            continue ;
          }
          
          hid_t di = H5Gopen(file_id,"dataInfo") ;
          int nbel = sizeElementType(di,"entityIds") ;
          
          vector<int> elemIds(nbel) ;
          readElementType(di,"entityIds",elemIds) ;
          
          H5Gclose(di) ;
          vector<float> var(nbel) ;
          readElementType(file_id,bnd_scalar_vars[b].c_str(),var) ;
          
          topo->output_boundary_scalar(&var[0],&elemIds[0],nbel,
                                       bnd_scalar_vars[b]) ;
        }
      }
      
      if(bnd_vector_vars.size() > 0) {
        for(size_t b=0;b<bnd_vector_vars.size();++b) {
          const string filename(bnd_vector_filenames[b]) ;
          hid_t file_id = Loci::hdf5OpenFile(filename.c_str(),
                                             H5F_ACC_RDONLY,
                                             H5P_DEFAULT) ;
          if(file_id < 0) {
            cerr << "unable to open file '" << filename << "'!" << endl ;
            continue ;
          }
          
          hid_t di = H5Gopen(file_id,"dataInfo") ;
          int nbel = sizeElementType(di,"entityIds") ;
          
          vector<int> elemIds(nbel) ;
          readElementType(di,"entityIds",elemIds) ;
          
          H5Gclose(di) ;
          vector<vector3d<float> > var(nbel) ;
          readElementType(file_id,bnd_vector_vars[b].c_str(),var) ;
          
          topo->output_boundary_vector(&var[0],&elemIds[0],nbel,
                                       bnd_vector_vars[b]) ;
        }
      }
      break;
    default:
      cerr << "internal error, problem with topo sequence" << endl ;
    }
  }
  topo->close() ;
  H5Gclose(bndg) ;
  H5Fclose(file_id) ;
}

int main(int ac, char *av[]) {
  Loci::Init(&ac,&av) ;

  enum {ASCII,TWODGV,ENSIGHT,FIELDVIEW,TECPLOT, NONE} plot_type = NONE ;

  string casename ;
  bool found_casename = false ;
  bool found_iteration = false ;
  string iteration ;
  vector<string> variables ;
  vector<string> boundaries ;

  int view = VIEWXY ;
  
  for(int i=1;i<ac;++i) {
    if(av[i][0] == '-') {
      if(!strcmp(av[i],"-ascii"))
        plot_type = ASCII ;
      else if(!strcmp(av[i],"-2d"))
        plot_type = TWODGV ;
      else if(!strcmp(av[i],"-en"))
        plot_type = ENSIGHT ;
      else if(!strcmp(av[i],"-fv"))
        plot_type = FIELDVIEW ;
      else if(!strcmp(av[i],"-tec"))
        plot_type = TECPLOT ;
      else if(!strcmp(av[i],"-xy")) 
        view=VIEWXY ;
      else if(!strcmp(av[i],"-yz")) 
        view=VIEWYZ ;
      else if(!strcmp(av[i],"-xz")) 
        view=VIEWXZ ;
      else if(!strcmp(av[i],"-xr")) 
        view=VIEWXR ;
      else if(!strcmp(av[i],"-bc")) {
        i++ ;
        string v(av[i]) ;
        if(av[i][0] >='0' && av[i][0] <= '9')
          v = "BC_"+v ;
        boundaries.push_back(v) ;
      } else {
        cerr << "unknown option " << av[i] << endl ;
        Usage(ac,av) ;
      }
      
    } else
      if(found_iteration)
        variables.push_back(string(av[i])) ;
      else if(found_casename) {
        iteration = string(av[i]) ;
        found_iteration = true ;
      } else {
        casename = string(av[i]) ;
        found_casename = true ;
      }
  }
  if(plot_type == NONE) {
    Usage(ac,av) ;
  }

  vector<int> variable_type(variables.size()) ;
  vector<string> variable_file(variables.size()) ;
  for(size_t i=0;i<variables.size();++i) {
    const string var(variables[i]) ;
    string filename = "output/" + var + "_hdf5." + iteration ;
    struct stat tmpstat ;
    if(stat(filename.c_str(),&tmpstat)== 0) {
      variable_type[i] = NODAL_SCALAR ;
      variable_file[i] = filename ;
      continue ;
    }
    // Handle old nodal scalar format
    if(var[0] == 'n') {
      filename = "output/" + string(var.c_str()+1) + "_hdf5." + iteration ;
      if(stat(filename.c_str(),&tmpstat)== 0) {
        variables[i] = string(var.c_str()+1) ;
        variable_type[i] = NODAL_SCALAR ;
        variable_file[i] = filename ;
        continue ;
      }
    }
      
    filename = "output/" + var + "_sca." + iteration + "_" + casename ;
    if(stat(filename.c_str(),&tmpstat)== 0) {
      variable_type[i] = NODAL_SCALAR ;
      variable_file[i] = filename ;
      continue ;
    }

    filename = "output/" + var + "_vec." + iteration + "_" + casename ;
    if(stat(filename.c_str(),&tmpstat)== 0) {
      variable_type[i] = NODAL_VECTOR ;
      variable_file[i] = filename ;
      continue ;
    }
    filename = "output/" + var + "_bnd." + iteration + "_" + casename ;
    if(stat(filename.c_str(),&tmpstat)== 0) {
      variable_type[i] = BOUNDARY_SCALAR ;
      variable_file[i] = filename ;
      continue ;
    }
    filename = "output/" + var + "_bndvec." + iteration + "_" + casename ;
    if(stat(filename.c_str(),&tmpstat)== 0) {
      variable_type[i] = BOUNDARY_VECTOR ;
      variable_file[i] = filename ;
      continue ;
    }

    if(plot_type == ASCII && boundaries.size() > 0) {
      if(var == "n") {
        variable_type[i] = BOUNDARY_DERIVED_VECTOR ;
        continue ;
      }
      if(var == "area" || var == "x" || var == "y" || var == "z") {
        variable_type[i] = BOUNDARY_DERIVED_SCALAR ;
        continue ;
      }
    }
    if(var == "m") {
      variable_type[i] = NODAL_DERIVED ;
      continue ;
    }
    if(var == "p") {
      variable_type[i] = NODAL_DERIVED ;
      continue ;
    }
    if(var == "P") {
      variable_type[i] = NODAL_DERIVED ;
      continue ;
    }
    if(var == "u") {
      variable_type[i] = NODAL_DERIVED ;
      continue ;
    }
    if(var == "0" || var == "1" || var == "2") {
      variable_type[i] = NODAL_DERIVED ;
      continue ;
    }
    if(var == "x" || var == "y" || var == "z") {
      variable_type[i] = NODAL_DERIVED ;
      continue ;
    }
      
    if(var.size()>1 && var[0] == 'f') {
      variable_type[i] = NODAL_MASSFRACTION ;
      continue ;
    }

    // ... other derived variables here

    cerr << "Warning, variable '" << var << "' is unknown and will not be processed." << endl ;
    variable_type[i] = UNDEFINED ;
  }

  H5Eset_auto(NULL,NULL) ;

  string filename = "output/" +  casename + ".topo" ;
  struct stat tmpstat ;
  if(stat(filename.c_str(),&tmpstat)!= 0) {
    cerr << "Warning, no grid topology file.  Will attempt to generate!"
         << endl ;
    setup_grid_topology(casename,iteration) ;
  }
      

  
  if(plot_type == ASCII) {
    if(boundaries.size() == 0) {
      process_ascii_nodal(casename,iteration,
                          variables,variable_type,variable_file) ;
      exit(0) ;
    } else {
      process_ascii_bndry(casename,iteration,
                          variables,variable_type,variable_file,
                          boundaries) ;
      exit(0) ;
    }
  }
  
  if(plot_type == TWODGV) {
    if(variables.size() != 1) {
      cerr << "2dgv extract can only extract one variable at a time."
           << endl ;
      Usage(ac,av) ;
    }
    if(boundaries.size() == 0) {
      cerr << "2dgv extract must have the projected boundaries identified using the '-bc' flag" << endl ;
      Usage(ac,av) ;
    }
    get_2dgv(casename,iteration,variables,variable_type,variable_file,
             boundaries,view) ;
    exit(0) ;
  }
  
  // process grid topology
  grid_topo_handler *topo_out = 0 ;
  switch(plot_type) {
  case ENSIGHT:
    topo_out = new ensight_topo_handler ;
    break ;
  case FIELDVIEW:
    topo_out = new fv_topo_handler ;
    break ;
  case TECPLOT:
    topo_out = new tecplot_topo_handler ;
    break ;
  default:
    cerr << "Unknown export method!" << endl ;
    break ;
  }

  if(topo_out != 0) {
    extract_grid(casename,iteration,topo_out,variables,variable_type,variable_file) ;
  }
}

