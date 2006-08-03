#include <Loci>
#include <LociGridReaders.h>


#include <Tools/tools.h>
#include <map>

#include <list>
using std::list ;
#include <string>
using std::string ;
#include <vector>
using std::vector ;
using std::pair ;

using std::cout ;

namespace Loci {

  // Oct Tree code
  
  void octree::collect_points(vector<int> &v, int node) {
    for(int i=0;i<8;++i) {
      if(otree[node][i]>0)
        collect_points(v,otree[node][i]) ;
      if(otree[node][i]<0)
        v.push_back(otree[node][i]) ;
    }
  }

  void octree::insert_tree(int node, int i,
                           double xl, double xh,
                           double yl, double yh,
                           double zl, double zh) {
    double xm = 0.5*(xl+xh) ;
    double ym = 0.5*(yl+yh) ;
    double zm = 0.5*(zl+zh) ;

    bool tx = vlist[i].x < xm ;
    bool ty = vlist[i].y < ym ;
    bool tz = vlist[i].z < zm ;
    int idx = (tx?0:4)+(ty?0:2)+(tz?0:1) ;

    double xo = tx?xl:xh ;
    double yo = ty?yl:yh ;
    double zo = tz?zl:zh ;
    
    xl = min(xo,xm) ;
    xh = max(xo,xm) ;
    yl = min(yo,ym) ;
    yh = max(yo,ym) ;
    zl = min(zo,zm) ;
    zh = max(zo,zm) ; 
    
    if(otree[node][idx] == 0) {
      otree[node][idx] = -(i+1) ;
      return ;
    }
    if(otree[node][idx] > 0) {
      insert_tree(otree[node][idx],i,xl,xh,yl,yh,zl,zh) ;
    } else {
      int j = (-otree[node][idx])-1 ;
      vector3d<double> dif = vlist[i]-vlist[j] ;
      double n1 = abs(dif.x)+abs(dif.y)+abs(dif.z) ;
      if(n1 < 1e-9 ) {
        cerr << "n1 = " << n1 << endl ;
        cerr << "duplicate points, ignoring second entry" << endl ;
        cerr << "pnt["<< i << "] = " << vlist[i]
             << ", pnt["<< j << "] = " << vlist[j] << endl ;

        return ;
      }
      otree.push_back(Loci::Array<int,8>()) ;
      int id = otree.size() - 1; 
      for(int k=0;k<8;++k) 
        otree[id][k] = 0 ;
      otree[node][idx] = id ;
      insert_tree(id,i,xl,xh,yl,yh,zl,zh) ;
      insert_tree(id,j,xl,xh,yl,yh,zl,zh) ;
    }
  }

  void octree::find_path(vector<unsigned char> &v,
                         int node, vector3d<double> &pt,
                         double xl, double xh,
                         double yl, double yh,
                         double zl, double zh) {

    vector<int> node_path ;
    do {
      double xm = 0.5*(xl+xh) ;
      double ym = 0.5*(yl+yh) ;
      double zm = 0.5*(zl+zh) ;

      bool tx = pt.x < xm ;
      bool ty = pt.y < ym ;
      bool tz = pt.z < zm ;
      int idx = (tx?0:4)+(ty?0:2)+(tz?0:1) ;

      double xo = tx?xl:xh ;
      double yo = ty?yl:yh ;
      double zo = tz?zl:zh ;
    
      xl = min(xo,xm) ;
      xh = max(xo,xm) ;
      yl = min(yo,ym) ;
      yh = max(yo,ym) ;
      zl = min(zo,zm) ;
      zh = max(zo,zm) ; 
      if(xl > pt.x || xh <pt.x )
        cerr << "xerror"<< endl ;
      if(yl > pt.y || yh <pt.y )
        cerr << "yerror"<< endl ;
      if(zl > pt.z || zh <pt.z )
        cerr << "zerror"<< endl ;
  
      v.push_back((unsigned char) idx) ;
      node_path.push_back(node) ;
      node = otree[node][idx] ;
    } while(node > 0) ;
  }
  void octree::find_neighbors(vector<int> &points,const vector<unsigned char> &v) {

    if(v.size() <= 1) {
      collect_points(points,0) ;
      for(unsigned int i=0;i<points.size();++i) {
        points[i] = -points[i]-1 ;
      }
      return ;
    } 
    vector<unsigned char > tmp = v ;
    tmp.pop_back() ;

    for(int i=0;i<3;++i) 
      for(int j=0;j<3;++j)
        for(int k=0;k<3;++k) {
          vector<unsigned char> path = tmp ;
          if(i==1)
            if(inc_path(path,0))
              continue ;
          if(i==2) 
            if(dec_path(path,0))
              continue ;
          if(j==1)
            if(inc_path(path,1))
              continue ;
          if(j==2) 
            if(dec_path(path,1))
              continue ;
          if(k==1)
            if(inc_path(path,2))
              continue ;
          if(k==2) 
            if(dec_path(path,2))
              continue ;
          int nd = get_node_from_path(path) ;
          if(nd > 0) {
            collect_points(points,nd) ;
          } if(nd < 0) {
            points.push_back(nd) ;
          }
        
        }
  
    for(unsigned int i=0;i<points.size();++i) {
      points[i] = -points[i]-1 ;
    }
    return ;
  }

  void octree::print_path(Loci::vector3d<double> &pt) {
    vector<unsigned char> path ;
    find_path(path,0,pt,xmin,xmax,ymin,ymax,zmin,zmax) ;
    std::cout << pt 
         << ' ' << xmin
         << ' ' << xmax
         << ' ' << ymin
         << ' ' << ymax
         << ' ' << zmin
         << ' ' << zmax
         << ' ' << "path = " ;
    for(unsigned int i=0;i<path.size();++i) {
      cout << (int)path[i] ;
    }
    cout << endl ;
  }
  
  namespace {
    void get_vect3dOption(const options_list &ol,std::string vname,
                          std::string units, vector3d<real_t> &vec, real_t Lref) {
      Loci::option_value_type ovt= ol.getOptionValueType(vname) ;
      if(ovt == Loci::REAL) {
        double v ;
        ol.getOption(vname,v) ;
        vec = vector3d<real_t>(v*Lref,0,0) ;
      } else if(ol.getOptionValueType(vname) == Loci::UNIT_VALUE) {
        Loci::UNIT_type vu ;
        ol.getOption(vname,vu) ;
        if(!vu.is_compatible(units)) {
          std::cerr << "wrong type of units for vector " << vname
                    << ": " << vu << std::endl ;
          Loci::Abort() ;
        } else {
          double v ;
          v = vu.get_value_in(units) ;
          vec = vector3d<real_t>(v,0,0) ;
        }
      } else if(ovt == Loci::LIST) {
        Loci::options_list::arg_list value_list ;
        ol.getOption(vname,value_list) ;
        if(value_list.size() != 3) {
          std::cerr << "error on reading '" << vname
                    <<"': vector input must contain 3 terms"
                    << std::endl ;
          Loci::Abort() ;
        }
        for(int i=0;i<3;++i)
          if(value_list[i].type_of() != Loci::REAL &&
             value_list[i].type_of() != Loci::UNIT_VALUE) {
            std::cerr << "improper vector specification for '"
                      << vname << std::endl ;
            Loci::Abort() ;
          }
        double vecval[3] ;
        for(int i=0;i<3;++i) {
          if(value_list[i].type_of() == Loci::UNIT_VALUE) {
            Loci::UNIT_type vu ;
            value_list[i].get_value(vu) ;
            if(!vu.is_compatible(units)) {
              std::cerr << "wrong type of units for vector " << vname
                        << ": " << vu << std::endl ;
              Loci::Abort() ;
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
      } else if(ovt == Loci::FUNCTION) {
        string name ;
        Loci::options_list::arg_list value_list ;
        ol.getOption(vname,name,value_list) ;
        if(name != "polar") {
          std::cerr << "don't know coordinate function '" << name
                    <<"', defaulting to polar" << std::endl ;
          Loci::Abort() ;
        }
        if(value_list.size() != 3) {
          std::cerr << "error on reading '"
                    << vname << "': vector input must contain 3 terms"
                    << std::endl ;
          Loci::Abort() ;
        }
        for(int i=0;i<3;++i)
          if(value_list[i].type_of() != Loci::REAL &&
             value_list[i].type_of() != Loci::UNIT_VALUE) {
            std::cerr << "improper vector specification for '"
                      << vname << std::endl ;
            Loci::Abort() ;
          }
        real_t r=1 ,theta=0 ,eta=0 ;
        real_t conv = M_PI/180.0 ;
        if(value_list[0].type_of() == Loci::UNIT_VALUE) {
          Loci::UNIT_type vu ;
          value_list[0].get_value(vu) ;
          if(!vu.is_compatible(units)) {
            std::cerr << "wrong type of units for vector " << vname
                      << ": " << vu << std::endl ;
            Loci::Abort() ;
          }
          r = vu.get_value_in(units) ;
        } else {
          value_list[0].get_value(r) ;
          r *= Lref ;
        }
        if(value_list[1].type_of() == Loci::UNIT_VALUE) {
          Loci::UNIT_type vu ;
          value_list[1].get_value(vu) ;
          if(!vu.is_compatible("radians")) {
            std::cerr << "wrong type of units for vector " << vname
                      << ": " << vu << std::endl ;
            Loci::Abort() ;
          }
          theta = vu.get_value_in("radians") ;
        } else {
          value_list[1].get_value(theta) ;
          theta *= conv  ;
        }
        if(value_list[2].type_of() == Loci::UNIT_VALUE) {
          Loci::UNIT_type vu ;
          value_list[2].get_value(vu) ;
          if(!vu.is_compatible("radians")) {
            std::cerr << "wrong type of units for vector " << vname
                      << ": " << vu << std::endl ;
            Loci::Abort() ;
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
        Loci::Abort() ;
      }
    }    

    void get_vect3d(const options_list &ol,std::string vname,
                    vector3d<real_t> &vec) {
      Loci::option_value_type ovt= ol.getOptionValueType(vname) ;
      if(ovt == Loci::LIST) {
        Loci::options_list::arg_list value_list ;
        ol.getOption(vname,value_list) ;
        if(value_list.size() != 3) {
          std::cerr << "error on reading '" << vname
                    <<"': vector input must contain 3 terms"
                    << std::endl ;
          Loci::Abort() ;
        }
        for(int i=0;i<3;++i)
          if(value_list[i].type_of() != Loci::REAL) {
            std::cerr << "improper vector specification for '"
                      << vname << std::endl ;
            Loci::Abort() ;
          }
        double vecval[3] ;
        for(int i=0;i<3;++i) {
          value_list[i].get_value(vecval[i]) ;
        }
        vec.x = vecval[0] ;
        vec.y = vecval[1] ;
        vec.z = vecval[2] ;
      } else if(ovt == Loci::FUNCTION) {
        string name ;
        Loci::options_list::arg_list value_list ;
        ol.getOption(vname,name,value_list) ;
        if(name != "polar") {
          std::cerr << "don't know coordinate function '" << name
                    <<"', defaulting to polar" << std::endl ;
          Loci::Abort() ;
        }
        if(value_list.size() != 3) {
          std::cerr << "error on reading '"
                    << vname << "': vector input must contain 3 terms"
                    << std::endl ;
          Loci::Abort() ;
        }
        if(value_list[0].type_of() != Loci::REAL) {
          std::cerr << "improper vector specification for '"
                    << vname << std::endl ;
          Loci::Abort() ;
        }
        for(int i=1;i<3;++i)
          if(value_list[i].type_of() != Loci::REAL &&
             value_list[i].type_of() != Loci::UNIT_VALUE) {
            std::cerr << "improper vector specification for '"
                      << vname << std::endl ;
            Loci::Abort() ;
          }
        real_t r=1 ,theta=0 ,eta=0 ;
        real_t conv = M_PI/180.0 ;
        value_list[0].get_value(r) ;
        if(value_list[1].type_of() == Loci::UNIT_VALUE) {
          Loci::UNIT_type vu ;
          value_list[1].get_value(vu) ;
          if(!vu.is_compatible("radians")) {
            std::cerr << "wrong type of units for vector " << vname
                      << ": " << vu << std::endl ;
            Loci::Abort() ;
          }
          theta = vu.get_value_in("radians") ;
        } else {
          value_list[1].get_value(theta) ;
          theta *= conv  ;
        }
        if(value_list[2].type_of() == Loci::UNIT_VALUE) {
          Loci::UNIT_type vu ;
          value_list[2].get_value(vu) ;
          if(!vu.is_compatible("radians")) {
            std::cerr << "wrong type of units for vector " << vname
                      << ": " << vu << std::endl ;
            Loci::Abort() ;
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
        Loci::Abort() ;
      }
    }  

    struct BCinfo {
      std::string name ;
      int key ;
      entitySet apply_set ;
      options_list bc_options ;
      BCinfo() {}
      BCinfo(const std::string &n, int k,
             const entitySet &a, const options_list &o) :
        name(n),key(k),apply_set(a),bc_options(o) {}
      
    } ; 
  }

  void point_connect_3d(const dstore<vector3d<double> > &p1, entitySet p1set,
                        const dstore<vector3d<double> > &p2, entitySet p2set,
                        dMap &leftc, dMap &rightc, double tol,
                        fact_db &facts) {
    int npnts = p2set.size() ;
    vector<vector3d<double> > pnts(npnts) ;
    vector<Entity> ids(npnts) ;
    int cnt = 0 ;
    vector3d<double>   bmin(1e30,1e30,1e30),bmax(-1e30,-1e30,-1e30) ;
    for(entitySet::const_iterator ei=p2set.begin();ei!=p2set.end();++ei) {
      const vector3d<double>  pt = p2[*ei] ;
      pnts[cnt] = pt ;
      ids[cnt++] = *ei ;
      bmin.x = min(pt.x,bmin.x) ;
      bmax.x = max(pt.x,bmax.x) ;
      bmin.y = min(pt.y,bmin.y) ;
      bmax.y = max(pt.y,bmax.y) ;
      bmin.z = min(pt.z,bmin.z) ;
      bmax.z = max(pt.z,bmax.z) ;
    }
    for(entitySet::const_iterator ei=p1set.begin();ei!=p1set.end();++ei) {
      const vector3d<double>  pt = p1[*ei] ;
      bmin.x = min(pt.x,bmin.x) ;
      bmax.x = max(pt.x,bmax.x) ;
      bmin.y = min(pt.y,bmin.y) ;
      bmax.y = max(pt.y,bmax.y) ;
      bmin.z = min(pt.z,bmin.z) ;
      bmax.z = max(pt.z,bmax.z) ;
    }

    bmax.x += 1e-33 ;
    bmax.y += 1e-33 ;
    bmax.z += 1e-33 ;
    bmin.x -= 1e-32 ;
    bmin.y -= 1e-32 ;
    bmin.z -= 1e-32 ;
    
    double refdist = norm(bmin-bmax)*1e-3 ;
    
    octree search_tree(pnts,bmin,bmax) ;

    vector<Entity> connect1, connect2 ;
    
    for(entitySet::const_iterator ei=p1set.begin();ei!=p1set.end();++ei) {
      const vector3d<double>  pt = p1[*ei] ;
      vector<int> neighbors ;
      search_tree.find_close_points(neighbors,pt) ;
      
      if(neighbors.size() == 0) {
        cerr << "neighbors size is zero in node_connect_3d!" << endl ;
        continue ;
      }

      double closedist = dot((pt-pnts[neighbors[0]]),
                           (pt-pnts[neighbors[0]])) ;
      int close_id = 0 ;
      for(size_t i=1;i<neighbors.size();++i) {
        double dist2 = dot((pt-pnts[neighbors[i]]),
                         (pt-pnts[neighbors[i]])) ;
        if(closedist > dist2) {
          closedist = dist2 ;
          close_id = i ;
        }
      }
      if(closedist <= refdist*refdist) {
        connect1.push_back(*ei) ;
        connect2.push_back(ids[neighbors[close_id]]) ;
      }
    }
    
    for(size_t i=0;i<connect1.size();++i) {
      leftc[connect1[i]] = connect1[i] ;
      rightc[connect1[i]] = connect2[i] ;
    }
    entitySet t1set = create_entitySet(connect1.begin(),connect1.end()) ;
    entitySet t2set = create_entitySet(connect2.begin(),connect2.end()) ;

#ifdef OLD
    if(t1set != p1set || t2set != p2set) {
      std::cerr << "no point matched connectivity in interface!" << std::endl ;
      cerr << "p1set = " << p1set << "t1set = " << t1set << endl ;
      cerr << "p2set = " << p2set << "t2set = " << t2set << endl ;
      cerr << "outputing debug file 'connect.debug'" << endl ;
      std::ofstream outfile("connect.debug",std::ios::out) ;
      outfile.precision(16) ;
      outfile << p1set.size() << endl ;
      for(entitySet::const_iterator ei= p1set.begin();ei!=p1set.end();++ei)
        outfile << *ei << ' '  <<p1[*ei] << endl ;
      outfile << p2set.size() << endl ;
      for(entitySet::const_iterator ei= p2set.begin();ei!=p2set.end();++ei)
        outfile << *ei << ' ' <<  p2[*ei] << endl ;
      
      outfile.close() ;
      Loci::Abort() ;
    }
#endif
  }

  void setup_periodic_bc(list<pair<periodic_info,periodic_info> >
                         &periodic_list,fact_db &facts) {

    dMap pmap ;
    dstore<rigid_transform> periodic_transform ;


    // Compute fluid face centers
    store<vector3d<real_t> > pos ;
    pos = facts.get_variable("pos") ;

    // First fill in tmp_pos so that it is valid for any reference to
    // it from face2node on this processor.
    multiMap face2node ;
    face2node = facts.get_variable("face2node") ;
    entitySet f2n_image = Loci::MapRepP(face2node.Rep())->image(face2node.domain()) ;
    entitySet out_of_dom = f2n_image - pos.domain() ;
    dstore<vector3d<real_t> > tmp_pos ;
    FORALL(pos.domain(), pi) {
      tmp_pos[pi] = pos[pi] ;
    } ENDFORALL ;
    Loci::storeRepP sp = tmp_pos.Rep() ;
    int tmp_out = out_of_dom.size() ;
    std::vector<entitySet> init_ptn ;
    if(facts.is_distributed_start()) {
      init_ptn = facts.get_init_ptn() ;
      if(GLOBAL_OR(tmp_out)) 
	fill_clone(sp, out_of_dom, init_ptn) ;
    }

    list<pair<periodic_info,periodic_info> >::const_iterator ii ;

    for(ii=periodic_list.begin();ii!=periodic_list.end();++ii) {
      int bc1 = ii->first.bc_num ;
      int bc2 = ii->second.bc_num ;
      real_t angle = ii->first.angle ;
      vector3d<real_t> center = ii->first.center ;
      vector3d<real_t> v = ii->first.v ;
      vector3d<real_t> trans = ii->first.translate ;
      
      periodic_transform[bc1] = rigid_transform(center,v,angle,trans) ;
      periodic_transform[bc2] = rigid_transform(center,v,-angle,-1.*trans) ;

      dstore<vector3d<real_t> > p1center ;
      entitySet p1Set = ii->first.bset ;
      rigid_transform tran = periodic_transform[bc1] ;
      for(entitySet::const_iterator ei = p1Set.begin();ei!=p1Set.end();++ei) {
        vector3d<real_t> tot = vector3d<real_t>(0.0,0.0,0.0);
        const int sz = face2node.end(*ei)-face2node.begin(*ei) ;
        for(int i=0; i<sz; ++i) {
          tot += tmp_pos[face2node[*ei][i]] ;
        }
        tot *= real_t(1)/real_t(sz) ;
        p1center[*ei] = tran.transform(tot) ;
      }

      dstore<vector3d<real_t> > p2center ;
      entitySet p2Set = ii->second.bset ;
      for(entitySet::const_iterator ei = p2Set.begin();ei!=p2Set.end();++ei) {
        vector3d<real_t> tot = vector3d<real_t>(0.0,0.0,0.0);
        const int sz = face2node.end(*ei)-face2node.begin(*ei) ;
        for(int i=0; i<sz; ++i) {
          tot += tmp_pos[face2node[*ei][i]] ;
        }
        tot *= real_t(1)/real_t(sz) ;
        p2center[*ei] = tot ;
      }

      if(facts.is_distributed_start()) {
        
        Loci::storeRepP sfp = p1center.Rep() ;
        Loci::storeRepP ffp = p2center.Rep() ;
        dstore<vector3d<real_t> > global_p1center ;
        global_p1center = Loci::collect_global_store(sfp) ;
        dstore<vector3d<real_t> > global_p2center ;
        global_p2center = Loci::collect_global_store(ffp) ;

      
        // Connect solid and fluid faces using centers.
        dMap lc,rc ;

        if(Loci::MPI_rank == 0) {
          point_connect_3d(global_p1center,global_p1center.domain(),
                           global_p2center,global_p2center.domain(),
                           lc,rc,
                           1e-5,facts) ;
          if(lc.domain() != global_p1center.domain() ||
             rc.image(rc.domain()) != global_p2center.domain()) {
            cerr << "periodic boundary not point matched!" << endl ;
            Loci::Abort() ;
          }
            
        }

        dMap lpmap_glob ;
        FORALL(lc.domain(),i) {
          lpmap_glob[lc[i]] = rc[i] ;
          lpmap_glob[rc[i]] = lc[i] ;
        } ENDFORALL ;

        dMap lpmap ;
        lpmap = Loci::distribute_dMap(lpmap_glob,init_ptn) ;

        FORALL(lpmap.domain(),i) {
          pmap[i] = lpmap[i] ;
        } ENDFORALL  ;

      } else {
        dMap lc , rc ;
        point_connect_3d(p1center,p1center.domain(),
                         p2center,p2center.domain(),
                         lc,rc,1e-5,facts) ;
        if(lc.domain() != p1center.domain() ||
           rc.image(rc.domain()) != p2center.domain()) {
          cerr << "periodic boundary not point matched!" << endl ;
          Loci::Abort() ;
        }

        FORALL(lc.domain(),i) {
          pmap[lc[i]] = rc[i] ;
          pmap[rc[i]] = lc[i] ;
        } ENDFORALL ;
      }
    }

    
    facts.create_fact("pmap",pmap) ;
    facts.create_fact("periodicTransform",periodic_transform) ;

    constraint pfaces ;
    Map cl ;
    pfaces = facts.get_variable("periodicFaces") ;
    *pfaces  = Loci::all_collect_entitySet(*pfaces) ;
    
    cl = facts.get_variable("cl") ;
    entitySet pcells = Loci::MapRepP(cl.Rep())->image(*pfaces) ;

    pcells = Loci::all_collect_entitySet(pcells) ;
    constraint periodicCells ;
    *periodicCells = pcells ;

    facts.create_fact("periodicCells",periodicCells) ;
    constraint notPeriodicCells ;
    *notPeriodicCells = ~pcells ;
    facts.create_fact("notPeriodicCells",notPeriodicCells) ;
  }

  void create_ci_map(fact_db &facts) {
    constraint boundary_faces ;
    boundary_faces = facts.get_variable("boundary_faces") ;
    entitySet ci_faces = *boundary_faces ;
    Loci::storeRepP pfacesP = facts.get_variable("periodicFaces") ;
    if(pfacesP != 0) {
      constraint periodicFaces ;
      periodicFaces = pfacesP ;
      Loci::debugout << "periodicFaces = " << periodicFaces << endl ;
      ci_faces -= *periodicFaces ;
    }

    Map cl,ci ;

    cl = facts.get_variable("cl") ;
    ci.allocate(ci_faces) ;
    
    FORALL(ci_faces,fc) {
      ci[fc] = cl[fc] ;
    } ENDFORALL ;
    facts.create_fact("ci",ci) ;
    Loci::debugout << "boundary_faces = " << *boundary_faces << endl ;
    Loci::debugout << "ci_faces = " << ci_faces << endl ;
  }

  void setupBoundaryConditions(fact_db &facts) {
    list<BCinfo> BCinfo_list ;
    std::map<std::string,entitySet> BCsets ;
    
    /*Boundary Conditions*/
    entitySet periodic ;
    constraint periodic_faces;
    constraint no_symmetry_BC ;

    entitySet symmetry ;
    
    store<string> boundary_names ;
    boundary_names = facts.get_variable("boundary_names") ;
    Map ref ;
    ref = facts.get_variable("ref") ;
    entitySet dom = boundary_names.domain() ;
    dom = all_collect_entitySet(dom) ;
    
    param<options_list> bc_info ;
    bc_info = facts.get_variable("boundary_conditions") ;
    param<real_t> Lref ;
    *Lref = 1.0 ;
    storeRepP p = facts.get_variable("Lref") ;
    if(p != 0)
      Lref = p ;
    
    vector<periodic_info> periodic_data ;

    for(entitySet::const_iterator ei=dom.begin();ei!=dom.end();++ei) {
      Entity bc = *ei ;
      entitySet bcset = interval(bc,bc) ;
      entitySet bfaces = ref.preimage(bcset).first ;

      string bname = boundary_names[bc] ;

      //      cout << "boundary_name =" << bname << endl ;
      constraint bconstraint ;
      *bconstraint = bfaces ;

      facts.create_fact(bname,bconstraint) ;
      Loci::debugout << "boundary " << bname << " = " << *bconstraint << endl ;
      
      Loci::option_value_type vt =
        bc_info->getOptionValueType(bname);
      Loci::option_values ov = bc_info->getOption(bname) ;
      options_list::arg_list value_list ;
      string name ;

      switch(vt) {
      case Loci::NAME :
        ov.get_value(name) ;
        bc_info->setOption(bname,name) ;
        {
          BCinfo_list.push_back(BCinfo(name,bc,bfaces,options_list())) ;
          BCsets[name] += bfaces ;
        }
        break ;
      case Loci::FUNCTION:
        ov.get_value(name) ;
        ov.get_value(value_list) ;
        bc_info->setOption(bname,name,value_list) ;
        {
          options_list ol ;
          ol.Input(value_list) ;
          BCinfo_list.push_back(BCinfo(name,bc,bfaces,ol)) ;
          BCsets[name] += bfaces ;
        }
        break ;
      default:
        cerr << "setup_bc can not interpret value assigned to " << bname 
             << " in boundary_conditions" << endl ;
        exit(-1) ;
      }
      if(name == "symmetry") {
        symmetry += bfaces ;
      } else if(name == "reflecting") {
        symmetry += bfaces ;
      } else if(name == "periodic") {
        periodic += bfaces ;
        periodic_info pi ;
        pi.bc_num = bc ;
        pi.bset = bfaces ;
        options_list ol ;
        ol.Input(value_list) ;
        if(ol.optionExists("rotate") || ol.optionExists("translate")) {
          pi.master = true ;
        }
        if(ol.optionExists("name")) {
          ol.getOption("name",pi.name) ;
        }
        if(ol.optionExists("center")) {
          get_vect3dOption(ol,"center","m",pi.center,*Lref) ;
        }
        if(ol.optionExists("vector")) {
          get_vect3d(ol,"vector",pi.v) ;
          pi.v /=  norm(pi.v) ;
        }
        if(ol.optionExists("translate")) {
          get_vect3dOption(ol,"translate","m",pi.translate,*Lref) ;
        }
        if(ol.optionExists("rotate")) {
          ol.getOptionUnits("rotate","radians",pi.angle) ;
        }
        
        periodic_data.push_back(pi) ;
      }
    } 

    
    {
      std::map<std::string, entitySet>::const_iterator mi ;
      for(mi=BCsets.begin();mi!=BCsets.end();++mi) {
        if(GLOBAL_OR(mi->second.size())) {
          constraint bc_constraint ;
          bc_constraint = mi->second ;
          std::string constraint_name = mi->first + std::string("_BC") ;
          facts.create_fact(constraint_name,bc_constraint) ;
          if(Loci::MPI_processes == 1)
            std::cout << constraint_name << ' ' << mi->second << endl ;
          else if(Loci::MPI_rank == 0)
            std::cout << "setting boundary condition " << constraint_name << endl ;
        }
      }
    }

 
    constraint cells ;
    cells = facts.get_variable("cells") ;
    store<options_list> BC_options ;
    BC_options.allocate(dom) ;

    std::map<std::string,entitySet> BC_options_args ;

    for(list<BCinfo>::iterator
          li = BCinfo_list.begin() ; li != BCinfo_list.end() ; ++li) {
      BC_options[li->key] = li->bc_options ;
      options_list::option_namelist onl=li->bc_options.getOptionNameList() ;
      options_list::option_namelist::const_iterator onli  ;

      for(onli=onl.begin();onli!=onl.end();++onli) {
        BC_options_args[*onli] += li->key ;
      }
    }

    {
      std::map<std::string, entitySet>::const_iterator mi ;
      for(mi=BC_options_args.begin();mi!=BC_options_args.end();++mi) {
        constraint bc_constraint ;
        bc_constraint = mi->second ;
        std::string constraint_name = mi->first + std::string("_BCoption") ;
        facts.create_fact(constraint_name,bc_constraint) ;
      }
    }
    
    if(periodic_data.size() != 0) {
      periodic_faces = periodic ;
      facts.create_fact("periodicFaces",periodic_faces) ;

      list<pair<periodic_info,periodic_info> > periodic_list ;
      for(size_t i=0;i<periodic_data.size();++i) {
        if(!periodic_data[i].processed) {
          periodic_data[i].processed = true ;
          periodic_info p1 = periodic_data[i] ;
          periodic_info p2 ;
          p2.name = "," ;
          for(size_t j=i+1;j<periodic_data.size();++j)
            if(periodic_data[i].name == periodic_data[j].name) {
              if(p2.name != ",") {
                cerr << "periodic name appears more than two times!" ;
                Loci::Abort() ;
              }
              p2 = periodic_data[j] ;
              periodic_data[j].processed = true ;
            }
          if(p1.name != p2.name) {
            cerr << "Could not find matching periodic boundary named "
                 << p1.name << endl ;
            Loci::Abort() ;
          }
          int p1inp = p1.bset.size() ;
          int p2inp = p2.bset.size() ;
          int p1size ;
          int p2size ;
          MPI_Allreduce(&p1inp,&p1size,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD) ;
          MPI_Allreduce(&p2inp,&p2size,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD) ;
          if(p1size != p2size) {
            if(Loci::MPI_rank == 0) {
            cerr << "periodic boundaries " << p1.name
                 << " do not match in number of faces" << endl ;
            cerr << "master has " << p1size << " faces and slave has "
                 << p2size << " faces" << endl ;
            }
            Loci::Abort() ;
          }
          if(p1.master & p2.master | (!p1.master & !p2.master)) {
            cerr << "only one master in periodic boundary conditons named "
                 << p1.name << endl ;
            Loci::Abort() ;
          }
          if(p1.master)
            periodic_list.push_back(std::make_pair(p1,p2)) ;
          else
            periodic_list.push_back(std::make_pair(p2,p1)) ;
        
        }
      }

      setup_periodic_bc(periodic_list,facts) ;
    } else {
      constraint notPeriodicCells ;
      *notPeriodicCells = ~EMPTY ;
      facts.create_fact("notPeriodicCells",notPeriodicCells) ;
    }      
    

    entitySet no_symmetry ;
    constraint allfaces ;
    allfaces = facts.get_variable("faces") ;
    no_symmetry  = *allfaces - symmetry ;
    no_symmetry_BC = no_symmetry ;
    facts.create_fact("no_symmetry_BC",no_symmetry_BC) ;

    facts.create_fact("BC_options",BC_options) ;

    create_ci_map(facts) ;
    
  }
  void createLowerUpper(fact_db &facts) {
    constraint geom_cells,interior_faces,boundary_faces ;
    constraint faces = facts.get_variable("faces") ;
    geom_cells = facts.get_variable("geom_cells") ;
    interior_faces = facts.get_variable("interior_faces") ;
    boundary_faces = facts.get_variable("boundary_faces") ;
    entitySet bfaces = *boundary_faces ;
    entitySet ifaces = *interior_faces ;

    Loci::storeRepP pfacesP = facts.get_variable("periodicFaces") ;
    if(pfacesP != 0) {
      constraint periodicFaces ;
      periodicFaces = pfacesP ;
      bfaces -= *periodicFaces ;
      ifaces += *periodicFaces ;
    }
    entitySet global_interior_faces = all_collect_entitySet(ifaces,facts) ;
    entitySet global_boundary_faces = all_collect_entitySet(bfaces,facts) ;
    
    Map cl,cr ;
    cl = facts.get_variable("cl") ;
    cr = facts.get_variable("cr") ;
    entitySet global_geom_cells ; 
    std::vector<entitySet> init_ptn ;

    global_geom_cells = all_collect_entitySet(*geom_cells,facts) ;
    multiMap lower,upper,boundary_map ;
    distributed_inverseMap(upper, cl, global_geom_cells, global_interior_faces,
                           facts) ;
    distributed_inverseMap(lower, cr, global_geom_cells, global_interior_faces,
                           facts) ;
    distributed_inverseMap(boundary_map, cl, global_geom_cells,
                           global_boundary_faces, facts) ;
    
    facts.create_fact("lower",lower) ;
    facts.create_fact("upper",upper) ;
    facts.create_fact("boundary_map",boundary_map) ;
  }

  
}
