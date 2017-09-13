//#############################################################################
//#
//# Copyright 2008, 2015, Mississippi State University
//#
//# This file is part of the Loci Framework.
//#
//# The Loci Framework is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The Loci Framework is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# Lesser GNU General Public License for more details.
//#
//# You should have received a copy of the Lesser GNU General Public License
//# along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################
#include <Loci>
#include <LociGridReaders.h>
#include <Tools/tools.h>
#include <map>
#include "pnn.h"
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
#include <list>
using std::list ;
#include <string>
using std::string ;
#include <vector>
using std::vector ;
using std::pair ;
#include <algorithm>
using std::sort ;
using std::unique ;

#include "dist_tools.h"
using std::cout ;

#define vect3d vector3d<double>

namespace Loci{
  // Convert container from local numbering to file numbering
  // pass in store rep pointer: sp
  // entitySet to write: dom
  // return offset in file numbering (each processor will allocate from zero,
  // add offset to domain to get actual file numbering)
  // distribution info pointer (dist)
  // MPI Communicator
  storeRepP Local2FileOrder(storeRepP sp, entitySet dom, int &offset,
                            fact_db::distribute_infoP dist, MPI_Comm comm);
  // Convert container from local numbering to output file numbering
  // pass in store rep pointer: sp
  // entitySet to write: dom
  // fact_db pointer  (facts)
  // MPI Communicator
  storeRepP Local2FileOrder_output(storeRepP sp, entitySet dom,
                                   fact_db& facts, MPI_Comm comm);
  
  
  //map from local numbering to input file numbering
  //assume the union of nodes on all processors will be either all the nodes,
  //all the faces, or all the cells. i.e., one interval in file numbering
  storeRepP get_node_remap(fact_db &facts,entitySet nodes) {

    if(MPI_processes == 1) {
      int minNode = nodes.Min() ;
      
      Map nm ;
      nm.allocate(nodes) ;
      FORALL(nodes,nd) {
        nm[nd] = nd - minNode + 1 ;
      } ENDFORALL ;
      return nm.Rep() ;
    }
    
    vector<entitySet> init_ptn = facts.get_init_ptn() ;
    fact_db::distribute_infoP df = facts.get_distribute_info() ;
    Map l2g ;
    l2g = df->l2g.Rep() ;
    dMap g2f ;
    g2f = df->g2f.Rep() ;

    entitySet gnodes = l2g.image(nodes&l2g.domain()) ;
    entitySet gset = findBoundingSet(gnodes) ;

    int minNode = gset.Min() ;

    Map newnum ;
    newnum.allocate(nodes) ;

    // Expand g2f to include clone regions
    entitySet out_of_dom = gnodes - init_ptn[MPI_rank] ;
    g2f.setRep(MapRepP(g2f.Rep())->expand(out_of_dom, init_ptn)) ;

    FORALL(nodes,i) {
      newnum[i] = g2f[l2g[i]]-minNode+1 ;
    } ENDFORALL ;
    return newnum.Rep() ;
  }


  //return a map from local numbering to output file numbering
  //the file number starts with 1
  //input: nodes : the local store domain need to be output,
  //assume nodes will be written out in the global ordering.   
  storeRepP get_output_node_remap(fact_db &facts,entitySet nodes) {
   
    if(MPI_processes == 1) {
      int index = 1;
      
      Map nm ;
      nm.allocate(nodes) ;
      FORALL(nodes,nd) {
        nm[nd] = index++ ;
      } ENDFORALL ;
      return nm.Rep() ;
    }
    
    vector<entitySet> init_ptn = facts.get_init_ptn() ;
    fact_db::distribute_infoP df = facts.get_distribute_info() ;
    Map l2g ;
    l2g = df->l2g.Rep() ;
    
    entitySet gnodes = l2g.image(nodes&l2g.domain()) ;
    if(nodes.size() != gnodes.size()){
      debugout<<"ERROR: l2g.domain is smaller than  dom in get_output_node_remap " << endl;
    }
    
    entitySet gset = all_collect_entitySet(gnodes) ;
    Map g2f;
    g2f.allocate(gset);
    
    
    int index = 1; //node index starts with 1
    FORALL(gset, n){
      g2f[n] = index++;
    }ENDFORALL;
    
    Map newnum ;
    newnum.allocate(nodes) ;

    FORALL(nodes,i) {
      newnum[i] = g2f[l2g[i]];
    } ENDFORALL ;
   
    return newnum.Rep() ;
    
  }
  
  int classify_cell(Entity *faces,int nfaces,const_multiMap &face2node) {
    int num_triangles = 0 ;
    int num_quads = 0 ;
    int num_others = 0 ;
    int triangle_nodes[3][2] ;
    for(int f=0;f<nfaces;++f) {
      Entity fc = faces[f] ;
      int count = face2node[fc].size() ;
      if(count == 3) {
        if(num_triangles < 2) {
          triangle_nodes[0][num_triangles] = face2node[fc][0] ;
          triangle_nodes[1][num_triangles] = face2node[fc][1] ;
          triangle_nodes[2][num_triangles] = face2node[fc][2] ;
        }
        num_triangles++ ;
      } else if(count == 4)
        num_quads++ ;
      else
        num_others++ ;
    }
    bool prism_test = false ;

    if((num_triangles == 2) && (num_quads == 3) && (num_others == 0)) {
      prism_test = true ;
      for(int i=0;i<3;++i)
        for(int j=0;j<3;++j)
          if(triangle_nodes[i][0] == triangle_nodes[j][1])
            prism_test = false ;
    }

    bool hex_test = false ;
    if( (num_triangles == 0) && (num_quads == 6) && (num_others == 0)) {
      const Entity ef = faces[0] ;
      int count = 0 ;
      for(int fj = 1;fj<nfaces;++fj) {
        Entity ej = faces[fj] ;
        bool find = false ;
        for(int i=0;i<4;++i)
          for(int j=0;j<4;++j)
            if(face2node[ef][i] == face2node[ej][j])
              find = true ;
        if(find)
          count++ ;
      }
      if(count == 4)
        hex_test = true ;
    }


    // new classification code
    if( (num_triangles == 4) && (num_quads == 0) && (num_others == 0)) {
      return 0 ;
    } else if( hex_test ) {
      return 1 ;
    } else if( prism_test ) {
      return 2 ;
    } else if( (num_triangles == 4) && (num_quads == 1) && (num_others == 0)) {
      return 3 ;
    }
    return 4 ;
    
  }

  void fillTet(Array<int,4> &tet, Array<int,3> *tri_faces) {
    tet[0] = tri_faces[0][2] ;
    tet[1] = tri_faces[0][1] ;
    tet[2] = tri_faces[0][0] ;
    for(int i=0;i<3;++i) {
      tet[3] = tri_faces[1][i] ;
      if(tet[3] != tet[0] && tet[3] != tet[1] && tet[3] != tet[2])
        return ;
    }
    cerr << "unable to form valid tet!" << endl ;
  }

  void fillHex(Array<int,8> &hex, Array<int,4> *quad_faces) {
    int quad_id[6] ;
    for(int i=0;i<6;++i) 
      quad_id[i] = i ;
    bool degenerate = quad_faces[quad_id[0]][0] == quad_faces[quad_id[0]][3];
    for(int j=0;j<3;++j) 
      if(quad_faces[quad_id[0]][j] == quad_faces[quad_id[0]][j+1])
        degenerate = true ;
    if(degenerate) {
      for(int i=1;i<6;++i) {
        degenerate = quad_faces[quad_id[i]][0] == quad_faces[quad_id[i]][3];
        for(int j=0;j<3;++j) 
          if(quad_faces[quad_id[i]][j] == quad_faces[quad_id[i]][j+1])
            degenerate = true ;
        if(!degenerate) {
          std::swap(quad_id[i],quad_id[0]) ;
          break ;
        }
      }
    }
    hex[0] = quad_faces[quad_id[0]][3] ;
    hex[1] = quad_faces[quad_id[0]][2] ;
    hex[2] = quad_faces[quad_id[0]][1] ;
    hex[3] = quad_faces[quad_id[0]][0] ;
    hex[4] = hex[0] ;
    hex[5] = hex[1] ;
    hex[6] = hex[2] ;
    hex[7] = hex[3] ;
    for(int i = 0; i < 4; i+=2) {
      int n1 = hex[i] ;
      int n2 = hex[i+1] ;
        
      int cnt = 0 ;
      for(int j=1; j<6; ++j) {
        for(int k=0;k<4;++k) {
          int f1 = quad_faces[quad_id[j]][k] ;
          int f2 = quad_faces[quad_id[j]][(k+1)%4] ;
          if((f1 == n1 && f2 == n2)) {
            hex[i+4] = quad_faces[quad_id[j]][(k-1+4)%4] ;
            hex[i+1+4] = quad_faces[quad_id[j]][(k+2)%4] ;
            cnt++ ;
          }
        }
      }
      if(cnt != 1) {
        cerr << "Error: Hex elem ordering screwed up " <<  endl ;
      }
    }
  }
  void fillPrism(Array<int,6> &prism, Array<int,3> *tri_faces,
                 Array<int,4> *quad_faces) {
    prism[0] = tri_faces[0][2] ;
    prism[1] = tri_faces[0][1] ;
    prism[2] = tri_faces[0][0] ;
    prism[3] = prism[0] ;
    prism[4] = prism[1] ;
    prism[5] = prism[2] ;
    
    int n1 = prism[0] ;
    int n2 = prism[1] ;
    int n3 = prism[2] ;
    int cnt = 0 ;
    for(int j=0;j<3;++j) {
      for(int k=0;k<4;++k) {
        int f1 = quad_faces[j][k] ;
        int f2 = quad_faces[j][(k+1)%4] ;
        if((f1 == n1 && f2 == n2)) {
          prism[3] = quad_faces[j][(k+3)%4] ;
          cnt++ ;
        }
          
        if((f1 == n2 && f2 == n3)) {
          prism[4] = quad_faces[j][(k+3)%4] ;
          prism[5] = quad_faces[j][(k+2)%4] ;
          cnt++ ;
        }
      }
    }
    if(cnt != 2) {
      cerr << "prism ordering screwed up" << endl ;
    }
  }
  void fillPyramid(Array<int,5> &pyramid, Array<int,3> *tri_faces,
                   Array<int,4> *quad_faces) {
    pyramid[0] = quad_faces[0][3] ;
    pyramid[1] = quad_faces[0][2] ;
    pyramid[2] = quad_faces[0][1] ;
    pyramid[3] = quad_faces[0][0] ;
    pyramid[4] = pyramid[0] ;
    for(int i=0;i<3;++i) {
      int nd = tri_faces[0][i] ;
      if(nd != pyramid[0] && nd != pyramid[1] &&
         nd != pyramid[2] && nd != pyramid[3]) {
        pyramid[4] = nd ;
        return ;
      }
    }
    cerr << "pyramid ordering screwed up!" << endl ;
  }


  string MPIConcatStrings(string input, MPI_Comm comm) {
    int sz = input.size()+1 ;
    int p = 1 ;
    MPI_Comm_size(comm,&p) ;
    int *sizes = new int[p] ;
    MPI_Allgather(&sz,1,MPI_INT,sizes,1,MPI_INT,comm) ;
    int tot = 0 ;
    for(int i=0;i<p;++i)
      tot += sizes[i] ;
    char *buf = new char[tot] ;
    int *displ = new int[p] ;
    displ[0] = 0 ;
    for(int i=1;i<p;++i)
      displ[i] = displ[i-1]+sizes[i-1] ;
    char *ibuf = new char[sz] ;
    strcpy(ibuf,input.c_str()) ;
    MPI_Allgatherv(ibuf,sz,MPI_CHAR,
                   buf, sizes, displ,MPI_CHAR,comm) ;
    string retval ;
    for(int i=0;i<p;++i)
      retval += string(&(buf[displ[i]])) ;

    delete[] sizes ;
    delete[] buf ;
    delete[] ibuf ;
    delete[] displ ;
    
    return retval ;
  }
  //node_id range is 1~numNodes
  //face_id range is numNodes~numNodes+numFaces
  //cell_id range is 1~numCells
  //Reason: in FVMstuff.cc, parallelWriteGridTopo(),
  //get_node_remap is used to get the file number of nodes and cells, which set the start index to 1
  //however, g2f(l2g(ff)) is used to get the file number of faces.
  //Since boundary variables also use the face_id, to avoid changing too many places, we keep the face_id this way.
  void parallelWriteGridTopology(const char *filename,
                                 storeRepP upperRep,
                                 storeRepP lowerRep,
                                 storeRepP boundary_mapRep,
                                 storeRepP face2nodeRep,
                                 storeRepP refRep,
                                 storeRepP bnamesRep,
                                 storeRepP posRep,
                                 entitySet localCells,
                                 fact_db &facts) {
    const_multiMap upper(upperRep),lower(lowerRep),
      boundary_map(boundary_mapRep),face2node(face2nodeRep) ;
    const_Map ref(refRep) ;

    store<int> elem_type ;
    elem_type.allocate(localCells) ;

    int ntets = 0 ;
    int nhexs = 0 ;
    int nprsm = 0 ;
    int npyrm = 0 ;
    int ngnrl = 0 ;
    
    // Classify Cells
    FORALL(localCells,cc) {
      int nfaces = upper[cc].size()+lower[cc].size()+boundary_map[cc].size() ;
      tmp_array<Entity> faces(nfaces) ;
      int cnt = 0 ;
      for(int i=0;i<upper[cc].size();++i) {
        faces[cnt++] = upper[cc][i] ;
      }
      for(int i=0;i<lower[cc].size();++i) {
        faces[cnt++] = lower[cc][i] ;
      }
      for(int i=0;i<boundary_map[cc].size();++i) {
        faces[cnt++] = boundary_map[cc][i] ;
      }
      
      elem_type[cc] = classify_cell(faces,nfaces,face2node) ;
      switch(elem_type[cc]) {
      case 0:
        ntets++; break ;
      case 1:
        nhexs++ ; break ;
      case 2:
        nprsm++ ; break ;
      case 3:
        npyrm++ ; break ;
      default:
        ngnrl++ ;
      }
    } ENDFORALL ;

    // Collect 4 cell type info
    vector<Array<int,4> > tets(ntets) ;
    vector<Array<int,5> > pyrm(npyrm) ;
    vector<Array<int,6> > prsm(nprsm) ;
    vector<Array<int,8> > hexs(nhexs) ;

    //for ids
    vector<int> tets_ids(ntets) ;
    vector<int> pyrm_ids(npyrm) ;
    vector<int> prsm_ids(nprsm) ;
    vector<int> hexs_ids(nhexs) ;
    vector<int> generalCell_ids;
    int tet_no = 0 ;
    int hex_no = 0 ;
    int pyramid_no = 0 ;
    int prism_no = 0 ;

    Map node_remap ;
    node_remap = get_node_remap(facts,posRep->domain()) ;

    Map cell_remap;
    cell_remap = get_node_remap(facts,localCells) ;

    vector<int> generalCellNfaces ;
    vector<int> generalCellNsides ;
    vector<int> generalCellNodes ;
    
    // Generate Cells
    FORALL(localCells,cc) {
      int nfaces = upper[cc].size()+lower[cc].size()+boundary_map[cc].size() ;
      tmp_array<int> faces(nfaces) ;
      tmp_array<int> swapface(nfaces) ;
      tmp_array<Array<int,3> > tri_faces(nfaces) ;
      tmp_array<Array<int,4> > quad_faces(nfaces) ;
      
      int tcnt = 0 ;
      int qcnt = 0 ;
      int nf = 0 ;
      for(int i=0;i<upper[cc].size();++i) {
        int fc = upper[cc][i] ;
        swapface[nf] = 0 ;
        faces[nf] = fc ;
        nf++ ;
        int fsz = face2node[fc].size() ;
        if(fsz == 3) {
          tri_faces[tcnt][0] = node_remap[face2node[fc][0]] ;
          tri_faces[tcnt][1] = node_remap[face2node[fc][1]] ;
          tri_faces[tcnt][2] = node_remap[face2node[fc][2]] ;
          tcnt++ ;
        }
        if(fsz == 4) {
          quad_faces[qcnt][0] = node_remap[face2node[fc][0]] ;
          quad_faces[qcnt][1] = node_remap[face2node[fc][1]] ;
          quad_faces[qcnt][2] = node_remap[face2node[fc][2]] ;
          quad_faces[qcnt][3] = node_remap[face2node[fc][3]] ;
          qcnt++ ;
        }
      }

      for(int i=0;i<lower[cc].size();++i) {
        int fc = lower[cc][i] ;
        swapface[nf] = 1 ;
        faces[nf] = fc ;
        nf++ ;
        int fsz = face2node[fc].size() ;
        if(fsz == 3) {
          tri_faces[tcnt][0] = node_remap[face2node[fc][2]] ;
          tri_faces[tcnt][1] = node_remap[face2node[fc][1]] ;
          tri_faces[tcnt][2] = node_remap[face2node[fc][0]] ;
          tcnt++ ;
        }
        if(fsz == 4) {
          quad_faces[qcnt][0] = node_remap[face2node[fc][3]] ;
          quad_faces[qcnt][1] = node_remap[face2node[fc][2]] ;
          quad_faces[qcnt][2] = node_remap[face2node[fc][1]] ;
          quad_faces[qcnt][3] = node_remap[face2node[fc][0]] ;
          qcnt++ ;
        }
      }

      for(int i=0;i<boundary_map[cc].size();++i) {
        int fc = boundary_map[cc][i] ;
        swapface[nf] = 0 ;
        faces[nf] = fc ;
        nf++ ;
        int fsz = face2node[fc].size() ;
        if(fsz == 3) {
          tri_faces[tcnt][0] = node_remap[face2node[fc][0]] ;
          tri_faces[tcnt][1] = node_remap[face2node[fc][1]] ;
          tri_faces[tcnt][2] = node_remap[face2node[fc][2]] ;
          tcnt++ ;
        }
        if(fsz == 4) {
          quad_faces[qcnt][0] = node_remap[face2node[fc][0]] ;
          quad_faces[qcnt][1] = node_remap[face2node[fc][1]] ;
          quad_faces[qcnt][2] = node_remap[face2node[fc][2]] ;
          quad_faces[qcnt][3] = node_remap[face2node[fc][3]] ;
          qcnt++ ;
        }
      }

      switch(elem_type[cc]) {
      case 0:
        tets_ids[tet_no] = cell_remap[cc];
        fillTet(tets[tet_no++],tri_faces) ;
        break ;
      case 1:
        hexs_ids[hex_no] = cell_remap[cc];
        fillHex(hexs[hex_no++],quad_faces) ;
        break ;
      case 2:
        prsm_ids[prism_no] = cell_remap[cc];
        fillPrism(prsm[prism_no++],tri_faces,quad_faces) ;
        break ;
      case 3:
        pyrm_ids[pyramid_no] = cell_remap[cc];
        fillPyramid(pyrm[pyramid_no++],tri_faces,quad_faces) ;
        break ;
      default:
        generalCellNfaces.push_back(nfaces) ;
        generalCell_ids.push_back(cell_remap[cc]);
        
        for(int i =0;i<nfaces;++i) {
          int fc = faces[i] ;
          int fsz = face2node[fc].size() ;
          generalCellNsides.push_back(fsz) ;
          if(swapface[i] == 1) {
            for(int j=0;j<fsz;++j)
              generalCellNodes.push_back(node_remap[face2node[fc][fsz-j-1]]) ;
          } else { 
            for(int j=0;j<fsz;++j)
              generalCellNodes.push_back(node_remap[face2node[fc][j]]) ;
          }
        }
      }
    } ENDFORALL ;


    // write grid topology file
    hid_t file_id = 0, group_id = 0 ;
    if(MPI_rank == 0) {
      file_id = H5Fcreate(filename,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT) ;
      group_id = H5Gcreate(file_id,"elements",
			   H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;
    }

    writeUnorderedVector(group_id, "tetrahedra",tets) ;
    writeUnorderedVector(group_id, "tetrahedra_ids",tets_ids) ;
    
    writeUnorderedVector(group_id, "hexahedra",hexs) ;
    writeUnorderedVector(group_id, "hexahedra_ids",hexs_ids) ;
    
    writeUnorderedVector(group_id, "prism",prsm) ;
    writeUnorderedVector(group_id, "prism_ids",prsm_ids) ;
    
    writeUnorderedVector(group_id, "pyramid",pyrm) ;
    writeUnorderedVector(group_id, "pyramid_ids",pyrm_ids) ;
    
    writeUnorderedVector(group_id, "GeneralCellNfaces",generalCellNfaces) ;
    writeUnorderedVector(group_id, "GeneralCellNsides",generalCellNsides) ;
    writeUnorderedVector(group_id, "GeneralCellNodes", generalCellNodes) ;
    writeUnorderedVector(group_id, "GeneralCell_ids", generalCell_ids) ;
    
    if(MPI_rank == 0) {
      H5Gclose(group_id) ;
      group_id = H5Gcreate(file_id,"boundaries",
			   H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;
    }

    const_store<string> boundary_names(bnamesRep) ;

    entitySet boundaries = boundary_names.domain() ;
    if(MPI_processes > 1) {
      entitySet local_boundaries ;
      std::vector<entitySet> init_ptn = facts.get_init_ptn() ;
      Map l2g ;
      fact_db::distribute_infoP df = facts.get_distribute_info() ;
      l2g = df->l2g.Rep() ;
      FORALL(boundaries,bb) {
        if(init_ptn[MPI_rank].inSet(l2g[bb]))
          local_boundaries += bb ;
      } ENDFORALL ;
      boundaries = local_boundaries ;
    }
    string bnames ;
    FORALL(boundaries,bb) {
      bnames += '"' + boundary_names[bb] + '"' ;
    } ENDFORALL ;


    bnames = MPIConcatStrings(bnames,MPI_COMM_WORLD) ;

    vector<string> bnamelist ;
    for(size_t i=0;i<bnames.size();++i) {
      if(bnames[i] != '"') {
        cerr << "confused in boundary name extraction" << endl ;
        break ;
      }
      ++i ;
      string tmp ;
      while(i < bnames.size() && bnames[i] != '"')
        tmp += bnames[i++] ;
      bnamelist.push_back(tmp) ;
    }


    // Identify Boundary faces to be written by this processor
    entitySet  fset = (MapRepP(boundary_mapRep)->image(localCells)+
                       MapRepP(upperRep)->image(localCells)) & ref.domain() ;

    for(size_t i=0;i<bnamelist.size();++i) {
      hid_t bc_id = 0 ;
      string current_bc = bnamelist[i] ;
      if(MPI_rank==0) {
        bc_id = H5Gcreate(group_id,current_bc.c_str(),
			  H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT) ;
      }

      bool found_ref = false ;
      Entity id = 0 ;
      FORALL(boundary_names.domain(),bb) {
        if(boundary_names[bb] == current_bc) {
          found_ref = true ;
          id = bb ;
        }
      } ENDFORALL ;
      entitySet bfaces ;
      if(found_ref) {
        FORALL(fset,fc) {
          if(ref[fc] == id)
            bfaces+= fc ;
        }ENDFORALL ;
      }
      Map l2g ;
      dMap g2f ;
      if(MPI_processes > 1) {
        fact_db::distribute_infoP df = facts.get_distribute_info() ;
        l2g = df->l2g.Rep() ;
        g2f = df->g2f.Rep() ;//why no map expanding? 
      } else {
        l2g.allocate(bfaces) ;
        FORALL(bfaces,fc) {
          l2g[fc] = fc ;
          g2f[fc] = fc ;
        } ENDFORALL ;
      }
      
      int ntria=0, nquad=0, nsided =0;
      FORALL(bfaces,fc) {
        if(face2node[fc].size() == 3)
          ntria++ ;
        else if(face2node[fc].size() == 4)
          nquad++ ;
        else
          nsided++ ;
      } ENDFORALL ;
      vector<Array<int,3> > Trias(ntria) ;
      vector<Array<int,4> > Quads(nquad) ;
      vector<int> tria_ids(ntria) ;
      vector<int> quad_ids(nquad) ;
      vector<int> genc_ids(nsided) ;
      int nt = 0 ;
      int nq = 0 ;
      int ng = 0 ;

      vector<int> nsizes(nsided) ;
      vector<int> nsidenodes ;

      FORALL(bfaces,fc) {
        if(face2node[fc].size() == 3) {
          Trias[nt][0] = node_remap[face2node[fc][0]] ;
          Trias[nt][1] = node_remap[face2node[fc][1]] ;
          Trias[nt][2] = node_remap[face2node[fc][2]] ;
          tria_ids[nt] = g2f[l2g[fc]] ;
          nt++ ;
        } else if(face2node[fc].size() == 4) {
          Quads[nq][0] = node_remap[face2node[fc][0]] ;
          Quads[nq][1] = node_remap[face2node[fc][1]] ;
          Quads[nq][2] = node_remap[face2node[fc][2]] ;
          Quads[nq][3] = node_remap[face2node[fc][3]] ;
          quad_ids[nq] = g2f[l2g[fc]] ;
          nq++ ;
        } else {
          nsizes[ng] = face2node[fc].size() ;
          for(int i=0;i<nsizes[ng];++i)
            nsidenodes.push_back(node_remap[face2node[fc][i]]) ;
          genc_ids[ng] = g2f[l2g[fc]] ;
          ng++ ;
        }
      } ENDFORALL ;
          
      writeUnorderedVector(bc_id,"triangles",Trias) ;
      writeUnorderedVector(bc_id,"triangles_id",tria_ids) ;

      writeUnorderedVector(bc_id,"quads",Quads) ;
      writeUnorderedVector(bc_id,"quads_id",quad_ids) ;

      writeUnorderedVector(bc_id,"nside_sizes",nsizes) ;
      writeUnorderedVector(bc_id,"nside_nodes",nsidenodes) ;
      writeUnorderedVector(bc_id,"nside_id",genc_ids) ;
      
      if(MPI_rank == 0) {
        H5Gclose(bc_id) ;
      }
      
    }
      
    if(MPI_rank == 0) {
      H5Gclose(group_id) ;
      H5Fclose(file_id) ;
    } 



  }
  
  //collect all boundary names
  std::vector<string> get_boundary_names(storeRepP bnamesRep, fact_db &facts){

   
    const_store<string> boundary_names(bnamesRep) ;
    entitySet boundaries = boundary_names.domain() ;
    if(MPI_processes > 1) {
      entitySet local_boundaries ;
      std::vector<entitySet> init_ptn = facts.get_init_ptn() ;
      Map l2g ;
      fact_db::distribute_infoP df = facts.get_distribute_info() ;
      l2g = df->l2g.Rep() ;
      FORALL(boundaries,bb) {
        if(init_ptn[MPI_rank].inSet(l2g[bb]))
          local_boundaries += bb ;
      } ENDFORALL ;
      boundaries = local_boundaries ;
    }
    
    string bnames ;
    FORALL(boundaries,bb) {
      bnames += '"' + boundary_names[bb] + '"' ;
    } ENDFORALL ;


    bnames = MPIConcatStrings(bnames,MPI_COMM_WORLD) ;

    vector<string> bnamelist ;
    for(size_t i=0;i<bnames.size();++i) {
      if(bnames[i] != '"') {
        cerr << "confused in boundary name extraction" << endl ;
        break ;
      }
      ++i ;
      string tmp ;
      while(i < bnames.size() && bnames[i] != '"')
        tmp += bnames[i++] ;
      bnamelist.push_back(tmp) ;
    }
    std::sort(bnamelist.begin(), bnamelist.end()); 
    return bnamelist;
  }
  
  //mkdir /output/$bc_name
  void get_bc_directory(string bc_name){
    string directory_name = "output/"+bc_name;
    struct stat statbuf ;
    int fid = open(directory_name.c_str(), O_RDONLY) ;
    if(fid < 0) {
      mkdir(directory_name.c_str(),0755) ;
    } else {
      fstat(fid,&statbuf) ;
      if(!S_ISDIR(statbuf.st_mode)) {
        cerr << "file " << directory_name <<" should be a directory!, rename it and start again."
             << endl ;
        Loci::Abort() ;
      }
      close(fid) ;
    }
  }
  //open /output/$bc_name/$file_name
  hid_t open_boundary_file(string bc_name,
                           string file_name
                           ){
    // open the file
    hid_t file_id = 0;
    if(MPI_rank == 0) {
      get_bc_directory(bc_name);
      string dirname = "output/"+bc_name+"/";
      string filename = dirname+file_name;
      file_id = H5Fcreate(filename.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT) ;
    
      if(file_id == 0){
        cerr << "ERROR: can not file " << filename << endl ;
        Loci::Abort() ;
      }
    }
    return file_id;
  }
  
 
  
  //get boundary faces that belong to a boundary surface current_bc
  entitySet get_boundary_faces(string current_bc,//boundary name
                               storeRepP refRep, // ref map
                               storeRepP bnamesRep,//bounadry name store
                               entitySet fset //all boundary faces 
                               ){
    
    const_store<string> boundary_names(bnamesRep) ;
    const_Map ref(refRep) ;
    
    //find its ref id
    bool found_ref = false ;
    Entity id = 0 ;
    FORALL(boundary_names.domain(),bb) {
      if(boundary_names[bb] == current_bc) {
        found_ref = true ;
        id = bb ;
      }
    } ENDFORALL ;
    entitySet bfaces ;
    if(found_ref) {
      FORALL(fset,fc) {
        if(ref[fc] == id)
          bfaces+= fc ;
      }ENDFORALL ;
    }
    return bfaces;
  }

  //get boundary nodes that belong to a boundary surface current_bc
  entitySet get_boundary_nodes(string current_bc,//boundary name
                               storeRepP face2nodeRep,
                               storeRepP refRep, // ref map
                               storeRepP bnamesRep,//bounadry name store
                               entitySet fset, //all boundary faces 
                               fact_db &facts ){
    entitySet bfaces = get_boundary_faces(current_bc, refRep, bnamesRep, fset);
    //get containers
    const_multiMap face2node(face2nodeRep) ;
    //get all the local nodes on the boundary 
    entitySet nodes_local = MapRepP(face2node)->image(bfaces);
    // get the nodes that belong to this processor
    if(MPI_processes > 1) {
      entitySet dom = nodes_local;
      fact_db::distribute_infoP dist = facts.get_distribute_info() ;
      entitySet  my_entities = dist->my_entities ;
      nodes_local = my_entities & nodes_local ;
    }
    return nodes_local;
  }
       
  void writeBoundaryTopo(hid_t file_id, //file_id of this boudnary surface
                         storeRepP face2nodeRep,
                         entitySet bfaces,//boundary faces define this surface 
                         fact_db &facts ){ 
    const_multiMap face2node(face2nodeRep) ;

    //collect the local boundary nodes belong to this boundary
    entitySet nodes_local = (MapRepP(face2node)->image(bfaces));
          
    //map each node to its file number in output file
    Map node_remap ;
    node_remap = get_output_node_remap(facts, nodes_local) ;

    // Compute face reordering for topo sorting
    store<int> faceorder ;
    faceorder.allocate(bfaces) ;
    int sz = bfaces.size();
    int off = 0 ;
    MPI_Scan(&sz,&off,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD) ;
    off -= sz ;
    FORALL(bfaces,fc) {
      faceorder[fc] = off ;
      off++ ;
    } ENDFORALL ;
    //get output vectors
    int ntria=0, nquad=0, nsided =0;
    FORALL(bfaces,fc) {
      if(face2node[fc].size() == 3)
        ntria++ ;
      else if(face2node[fc].size() == 4)
        nquad++ ;
      else
        nsided++ ;
    } ENDFORALL ;
    vector<Array<int,3> > Trias(ntria) ;
    vector<Array<int,4> > Quads(nquad) ;
    vector<int> tria_ids(ntria) ;
    vector<int> quad_ids(nquad) ;
    vector<int> genc_ids(nsided) ;
    int nt = 0 ;
    int nq = 0 ;
    int ng = 0 ;
    
    vector<int> nsizes(nsided) ;
    vector<int> nsidenodes ;
   
    
    FORALL(bfaces,fc) {
      if(face2node[fc].size() == 3) {
        Trias[nt][0] = node_remap[face2node[fc][0]] ;
        Trias[nt][1] = node_remap[face2node[fc][1]] ;
        Trias[nt][2] = node_remap[face2node[fc][2]] ;
        tria_ids[nt] = faceorder[fc] ;
        nt++ ;
      } else if(face2node[fc].size() == 4) {
        Quads[nq][0] = node_remap[face2node[fc][0]] ;
        Quads[nq][1] = node_remap[face2node[fc][1]] ;
        Quads[nq][2] = node_remap[face2node[fc][2]] ;
        Quads[nq][3] = node_remap[face2node[fc][3]] ;
          
        quad_ids[nq] = faceorder[fc] ; 
        nq++ ;
      } else {
        nsizes[ng] = face2node[fc].size() ;
        for(int i=0;i<nsizes[ng];++i)
          nsidenodes.push_back(node_remap[face2node[fc][i]]) ;
        
        genc_ids[ng] = faceorder[fc] ; 
        ng++ ;
      }
    } ENDFORALL ;
    
      
    //write out vectors
    writeUnorderedVector(file_id,"triangles",Trias) ;
    writeUnorderedVector(file_id,"triangles_ord",tria_ids) ;
    writeUnorderedVector(file_id,"quads",Quads) ;
    writeUnorderedVector(file_id,"quads_ord",quad_ids) ;
    
    writeUnorderedVector(file_id,"nside_sizes",nsizes) ;
    writeUnorderedVector(file_id,"nside_nodes",nsidenodes) ;
    writeUnorderedVector(file_id,"nside_ord",genc_ids) ;
  }
      
  //this function find the index of an inner edge.
  //an inner edge is the edge from facecenter to one of its nodes during face triangulation
  //if the edge already exists, the index is returned
  //otherwise, the edge is created and registered in inner_edges and edgeIndex
  //the return value is the index of the edge in inner_edges, starts with 1
  int create_inner_edge(int face, //face id
                        int node, //node id
                        std::map<pair<int, int>, int>& edgeIndex, //map of inner edges to its index in inner_edges, the index starts with 1 
                        vector<pair<int, int> >& inner_edges) //the edges from facecenter and one of its node
  {
    pair<int, int> inner_edge = make_pair(face, node);
    //make sure no duplicate in inner_edges
    if(edgeIndex.find(inner_edge) == edgeIndex.end()){
      inner_edges.push_back(inner_edge);
      edgeIndex[inner_edge] = inner_edges.size();
    }
    return edgeIndex[inner_edge];
  }
  
  //If a face has more than 2 edges cut, to disambiguate the connection between the cut points,
  //the face will be triangulated and each tri face will be registered
  void disambiguateFace(Entity f, 
                        const_multiMap& face2node,
                        const_multiMap& face2edge,
                        const_MapVec<2>& edge2node,
                        const_store<real_t>& val,
                        entitySet& edgesCut, //the edges(node 2 node) cut
                        std::map<pair<int, int>, int >& edgeIndex,//map of inner edges(face 2 node) to its index in inner_edges, the index starts with 1 
                        vector<pair<int, int> >& intersects, //pair of edges, if the second number is negative, it is index to inner_edges
                        vector<pair<int, int> >& inner_edges //the edges from facecenter and one of its node 
                        ){
    
    int nNodes = face2node.num_elems(f);
    
    //get the value of face center
    real_t  newNode = 0.0 ;
    for (int i = 0; i < nNodes; ++i) {
      newNode += val[face2node[f][i]];
    }
    newNode /= double(nNodes);
    
    //triangulate the face
    //and register each tri face.
    int faceCut[2];
    for (int i = 0;i < face2edge.num_elems(f); ++i) {
      int cutsFound = 0;
      //check this edge
      Entity edge = face2edge[f][i];
      if (signbit(val[edge2node[edge][0]] * val[edge2node[edge][1]])) {
        edgesCut += edge;
        faceCut[cutsFound] = edge;
        cutsFound++; 
      }
      //check face center to first node
      if (signbit(newNode * val[edge2node[edge][0]])) {
        int node = edge2node[edge][0];
        int edgeId = create_inner_edge(f, node, edgeIndex,inner_edges);
        faceCut[cutsFound]  = -edgeId;
        cutsFound++;                     
      }
  
      //check face center to second node
      if (signbit(newNode * val[edge2node[edge][1]])) {
        if(cutsFound >1){
          cerr<<"ERROR: tri face has more than two cuts" << endl;
          exit(-1);
        }
        int node = edge2node[edge][1];
        int edgeId = create_inner_edge(f, node, edgeIndex,inner_edges);
        faceCut[cutsFound]  = -edgeId;
        cutsFound++;                     
      }
      if (cutsFound == 2) intersects.push_back(make_pair(faceCut[0], faceCut[1]));
      //else ????holes happen?
    }
  }

  // Find edges cut in the face and register them into intersects
  bool registerFace(Entity f,
                    const_multiMap& face2node,
                    const_multiMap& face2edge,
                    const_MapVec<2>& edge2node,
                    const_store<real_t>& val,
                    entitySet& edgesCut,//the edges(node 2 node) cut
                    std::map<pair<int, int>, int >& edgeIndex,//map of inner edges(face center 2 node) to its index in inner_edges, the index starts with 1
                    vector<pair<int, int> >& intersects, //pair of edges, if the value is negative, it is index to inner_edges
                    vector<pair<int, int> >& inner_edges//pair<face, node>, the edges from facecenter and one of the face node 
                    ){
    
    int faceCut[2];
    int cutsFound = 0;
    //check each edge
    for(int i = 0; i < face2edge.num_elems(f); ++i) {
      Entity edge = face2edge[f][i];
      //if it is cut
      if (signbit(val[edge2node[edge][0]] * val[edge2node[edge][1]])) {
        if (cutsFound < 2) {
          //store it in faceCut
          faceCut[cutsFound] = edge;
          edgesCut += edge;
        }
        cutsFound++; 
      }
    }
    //no ambiguation, register faceCut in inertsection
    //otherwise, disambiguate the face
    if (cutsFound == 2)
      intersects.push_back(make_pair(faceCut[0], faceCut[1]));
    else if (cutsFound > 2)
      disambiguateFace( f,
                        face2node,
                        face2edge,
                        edge2node,
                        val,
                        edgesCut,
                        edgeIndex,
                        intersects,
                        inner_edges
                        );
    return (cutsFound > 0);
  }
    
  //after all the faces of a cell is registered,
  //check the intersects to form the faces in cutting plane
  void checkLoop( const vector<pair<int, int> >& intersects,//pair of edges, if the value number is negative, it is index to inner_edges
                  vector<vector<int > >& faceLoops,//the loops(faces in cutting plane) formed 
                  int start,//the start index to intersects 
                  int end) { //the end index to intersects
    bool loopIsGood = true;
    vector<int > loop;
  
    list<int> edges;
    for (int i = start+1; i < end; i++)
      edges.push_back(i);

    int firstNode, nextNode;
    firstNode = intersects[start].first;
    nextNode = intersects[start].second;
    loop.push_back(firstNode);
    loop.push_back(nextNode);

    list<int>::iterator iter;
    while (!edges.empty() && loopIsGood) {
      for (iter = edges.begin(); iter != edges.end(); iter++) {
        if (intersects[*iter].first == nextNode) {
          nextNode = intersects[*iter].second;
          if(firstNode != nextNode)loop.push_back(nextNode);
          edges.erase(iter);
          break;
        } else if (intersects[*iter].second == nextNode) {
          nextNode = intersects[*iter].first;
          if(firstNode != nextNode) loop.push_back(nextNode);
          edges.erase(iter);
          break;
        }
      }
      //can not find the next edge
      if (iter == edges.end())
        loopIsGood = false;

      //a loop is formed
      if (firstNode == nextNode){
        if(!edges.empty()) {
          faceLoops.push_back(loop);
          loop.clear();
        
          firstNode = intersects[edges.front()].first;
          nextNode = intersects[edges.front()].second;
          loop.push_back(firstNode);
          loop.push_back(nextNode);
          edges.erase(edges.begin());
        }else{
          faceLoops.push_back(loop);
        }
      }
    }
      
    if (firstNode != nextNode){
      debugout << "ERROR: ** Problem cell:  ** (failed loop test)" << endl;
    }
  }
  //after the loops are formed, this function will remove all negative value in loops.
  vector<vector<int > > remove_inner_edges(const vector<vector<int > >& faceLoops)//the loops(faces in cutting plane) formed 
  {
    vector<vector<int > >  new_faceLoops;
    for(unsigned int i = 0; i < faceLoops.size(); i++){
      vector<int> loop;
      for(unsigned int j = 0; j < faceLoops[i].size(); j++){
        if(faceLoops[i][j] >=0) loop.push_back(faceLoops[i][j]);
      }
      if(loop.size() >0)new_faceLoops.push_back(loop);
    }
    return new_faceLoops;
  }
    
  //return the cutting position of an edge 
  double get_edge_weight(Entity e,
                         const_MapVec<2>& edge2node,
                         const_store<real_t>& val){
    
    real_t a = val[edge2node[e][0]];
    real_t b = val[edge2node[e][1]];
    return ((b)/(b - a));
  }
 

  //generate Cutplane
  CutPlane getCutPlane(storeRepP upperRep,
                       storeRepP lowerRep,
                       storeRepP boundary_mapRep,
                       storeRepP face2nodeRep,
                       storeRepP face2edgeRep,
                       storeRepP edge2nodeRep,
                       storeRepP valRep,
                       entitySet localCells,//all geom_cells
                       fact_db &facts) {
    
    const_multiMap upper(upperRep),lower(lowerRep),
      boundary_map(boundary_mapRep),face2node(face2nodeRep), face2edge(face2edgeRep) ;
    
    const_MapVec<2> edge2node(edge2nodeRep);
    const_store<real_t> val(valRep);

    // data structure:
    entitySet edgesCut;
    store<double> edgesWeight; //the weight for interpoplation for each edgesCut, allocated on edgesCut
    vector<pair<int, int> > inner_edges; //the inner edges(facecenter to one of its nodes) cut, the values stored are pair<local_faceid, noderank>  
    vector<vector<int > > faceLoops;  //loops formed, the values stored are edge ids, which is either local edge entity or index to inner_edges
   
   
    //extra data structure
    vector<pair<int, int> > intersects; //the pairs of edges  or inner edges cut
    std::map<pair<int, int>, int > edgeIndex; //map of inner edges to their indexes in inner_edges
    // for each cell
    FORALL(localCells,cc) {
      bool  isCut = false;
      int  intersectStart = intersects.size();
      
      //collect all its faces
      entitySet faces;
      for(int i=0;i<upper[cc].size();++i)faces += upper[cc][i];
      for(int i=0;i<lower[cc].size();++i)faces += lower[cc][i];
      for(int i=0;i<boundary_map[cc].size();++i)faces += boundary_map[cc][i];
      
      //for each face, find the edges that is cut and register it
      for(entitySet::const_iterator ei = faces.begin(); ei != faces.end(); ei++){
        if (registerFace(*ei,
                         face2node,
                         face2edge,
                         edge2node,
                         val,
                         edgesCut,
                         edgeIndex,
                         intersects, 
                         inner_edges 
                         )) isCut = true;
      }

      //check the loops formed by this cell and register the loop
      if(isCut){
        checkLoop(intersects, faceLoops, intersectStart, intersects.size());
        intersectStart = intersects.size();
      }
    }ENDFORALL;
    //remove the inner edges cut from loops
    vector<vector<int > > new_faceLoops = remove_inner_edges(faceLoops);
    
    //compute the cutting postions of edges 
    edgesWeight.allocate(edgesCut);
    FORALL(edgesCut, e){
      double t = get_edge_weight(e, edge2node, val);
      edgesWeight[e] = t;
    }ENDFORALL;
    
   
    return CutPlane(edgesWeight.Rep(),
                    new_faceLoops);
                     
   
  }

  
  void writeCutPlaneTopo(hid_t bc_id,
                         const CutPlane& cp,
                         fact_db &facts){ 
    Map node_remap;//map from local numbering to output file numbering for edge entities
    entitySet edgesCut = (cp.edgesWeight)->domain();
    node_remap = get_output_node_remap(facts, edgesCut);
    
   
    //write out the sizes of faceLoops  
    int num_faces = cp.faceLoops.size();
    vector<int> nsizes(num_faces);
    for(int i = 0; i < num_faces; i++){
      nsizes[i] = cp.faceLoops[i].size();
    }
    
    writeUnorderedVector(bc_id,"nside_sizes",nsizes) ;
    
    //write out face nodes
    vector<int> nsidenodes ;
    for(int i = 0; i < num_faces; i++){
      for(int j = 0; j < nsizes[i]; j++){
        int node = cp.faceLoops[i][j];
        if(node >=0){
          nsidenodes.push_back(node_remap[node]); //edge nodes
        }
      }
    }
    writeUnorderedVector(bc_id,"nside_nodes",nsidenodes) ;
  }
  
  // namespace {
  void get_vect3dOption(const options_list &ol,std::string vname,
                        std::string units, vector3d<real_t> &vec, real_t Lref) {
    option_value_type ovt= ol.getOptionValueType(vname) ;
    if(ovt == REAL) {
      double v ;
      ol.getOption(vname,v) ;
      vec = vector3d<real_t>(v*Lref,0,0) ;
    } else if(ol.getOptionValueType(vname) == UNIT_VALUE) {
      UNIT_type vu ;
      ol.getOption(vname,vu) ;
      if(!vu.is_compatible(units)) {
        std::cerr << "wrong type of units for vector " << vname
                  << ": " << vu << std::endl ;
        Abort() ;
      } else {
        double v ;
        v = vu.get_value_in(units) ;
        vec = vector3d<real_t>(v,0,0) ;
      }
    } else if(ovt == LIST) {
      options_list::arg_list value_list ;
      ol.getOption(vname,value_list) ;
      if(value_list.size() != 3) {
        std::cerr << "error on reading '" << vname
                  <<"': vector input must contain 3 terms"
                  << std::endl ;
        Abort() ;
      }
      for(int i=0;i<3;++i)
        if(value_list[i].type_of() != REAL &&
           value_list[i].type_of() != UNIT_VALUE) {
          std::cerr << "improper vector specification for '"
                    << vname << std::endl ;
          Abort() ;
        }
      double vecval[3] ;
      for(int i=0;i<3;++i) {
        if(value_list[i].type_of() == UNIT_VALUE) {
          UNIT_type vu ;
          value_list[i].get_value(vu) ;
          if(!vu.is_compatible(units)) {
            std::cerr << "wrong type of units for vector " << vname
                      << ": " << vu << std::endl ;
            Abort() ;
          }
          vecval[i] = vu.get_value_in(units) ;
        } else {
          value_list[i].get_value(vecval[i]) ;
          vecval[i] *= Lref ;
        }
      }
      vec.x = vecval[0] ;
      vec.y = vecval[1] ;
      vec.z = vecval[2] ;
    } else if(ovt == FUNCTION) {
      string name ;
      options_list::arg_list value_list ;
      ol.getOption(vname,name,value_list) ;
      if(name != "polar") {
        std::cerr << "don't know coordinate function '" << name
                  <<"', defaulting to polar" << std::endl ;
        Abort() ;
      }
      if(value_list.size() != 3) {
        std::cerr << "error on reading '"
                  << vname << "': vector input must contain 3 terms"
                  << std::endl ;
        Abort() ;
      }
      for(int i=0;i<3;++i)
        if(value_list[i].type_of() != REAL &&
           value_list[i].type_of() != UNIT_VALUE) {
          std::cerr << "improper vector specification for '"
                    << vname << std::endl ;
          Abort() ;
        }
      real_t r=1 ,theta=0 ,eta=0 ;
      real_t conv = M_PI/180.0 ;
      if(value_list[0].type_of() == UNIT_VALUE) {
        UNIT_type vu ;
        value_list[0].get_value(vu) ;
        if(!vu.is_compatible(units)) {
          std::cerr << "wrong type of units for vector " << vname
                    << ": " << vu << std::endl ;
          Abort() ;
        }
        r = vu.get_value_in(units) ;
      } else {
        value_list[0].get_value(r) ;
        r *= Lref ;
      }
      if(value_list[1].type_of() == UNIT_VALUE) {
        UNIT_type vu ;
        value_list[1].get_value(vu) ;
        if(!vu.is_compatible("radians")) {
          std::cerr << "wrong type of units for vector " << vname
                    << ": " << vu << std::endl ;
          Abort() ;
        }
        theta = vu.get_value_in("radians") ;
      } else {
        value_list[1].get_value(theta) ;
        theta *= conv  ;
      }
      if(value_list[2].type_of() == UNIT_VALUE) {
        UNIT_type vu ;
        value_list[2].get_value(vu) ;
        if(!vu.is_compatible("radians")) {
          std::cerr << "wrong type of units for vector " << vname
                    << ": " << vu << std::endl ;
          Abort() ;
        }
        eta = vu.get_value_in("radians") ;
      } else {
        value_list[2].get_value(eta) ;
        eta *= conv  ;
      }
      
      vec.x = r*cos(theta)*cos(eta) ;
      vec.y = r*sin(theta)*cos(eta) ;
      vec.z = r*sin(eta) ;
    } else {
      std::cerr << "unable to get vector type!" << std::endl ;
      Abort() ;
    }
  }    

  void get_vect3d(const options_list &ol,std::string vname,
                  vector3d<real_t> &vec) {
    option_value_type ovt= ol.getOptionValueType(vname) ;
    if(ovt == LIST) {
      options_list::arg_list value_list ;
      ol.getOption(vname,value_list) ;
      if(value_list.size() != 3) {
        std::cerr << "error on reading '" << vname
                  <<"': vector input must contain 3 terms"
                  << std::endl ;
        Abort() ;
      }
      for(int i=0;i<3;++i)
        if(value_list[i].type_of() != REAL) {
          std::cerr << "improper vector specification for '"
                    << vname << std::endl ;
          Abort() ;
        }
      double vecval[3] ;
      for(int i=0;i<3;++i) {
        value_list[i].get_value(vecval[i]) ;
      }
      vec.x = vecval[0] ;
      vec.y = vecval[1] ;
      vec.z = vecval[2] ;
    } else if(ovt == FUNCTION) {
      string name ;
      options_list::arg_list value_list ;
      ol.getOption(vname,name,value_list) ;
      if(name != "polar") {
        std::cerr << "don't know coordinate function '" << name
                  <<"', defaulting to polar" << std::endl ;
        Abort() ;
      }
      if(value_list.size() != 3) {
        std::cerr << "error on reading '"
                  << vname << "': vector input must contain 3 terms"
                  << std::endl ;
        Abort() ;
      }
      if(value_list[0].type_of() != REAL) {
        std::cerr << "improper vector specification for '"
                  << vname << std::endl ;
        Abort() ;
      }
      for(int i=1;i<3;++i)
        if(value_list[i].type_of() != REAL &&
           value_list[i].type_of() != UNIT_VALUE) {
          std::cerr << "improper vector specification for '"
                    << vname << std::endl ;
          Abort() ;
        }
      real_t r=1 ,theta=0 ,eta=0 ;
      real_t conv = M_PI/180.0 ;
      value_list[0].get_value(r) ;
      if(value_list[1].type_of() == UNIT_VALUE) {
        UNIT_type vu ;
        value_list[1].get_value(vu) ;
        if(!vu.is_compatible("radians")) {
          std::cerr << "wrong type of units for vector " << vname
                    << ": " << vu << std::endl ;
          Abort() ;
        }
        theta = vu.get_value_in("radians") ;
      } else {
        value_list[1].get_value(theta) ;
        theta *= conv  ;
      }
      if(value_list[2].type_of() == UNIT_VALUE) {
        UNIT_type vu ;
        value_list[2].get_value(vu) ;
        if(!vu.is_compatible("radians")) {
          std::cerr << "wrong type of units for vector " << vname
                    << ": " << vu << std::endl ;
          Abort() ;
        }
        eta = vu.get_value_in("radians") ;
      } else {
        value_list[2].get_value(eta) ;
        eta *= conv  ;
      }

      vec.x = r*cos(theta)*cos(eta) ;
      vec.y = r*sin(theta)*cos(eta) ;
      vec.z = r*sin(eta) ;
    } else {
      std::cerr << "unable to get vector type!" << std::endl ;
      Abort() ;
    }
  }  
}

