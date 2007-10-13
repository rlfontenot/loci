#include <Loci.h>
#include "defines.h"
using std::cout ;
using std::vector;
namespace Loci{
  int classify_cell(Entity *faces,int nfaces,const_multiMap &face2node);

  void parallelClassifyCell(fact_db &facts) {
    //get variables
    const_multiMap face2node ;
    face2node = facts.get_variable("face2node") ;
    const_multiMap upper ;
    upper = facts.get_variable("upper") ;
    const_multiMap lower ;
    lower = facts.get_variable("lower") ;
    const_multiMap boundary_map ;
    boundary_map = facts.get_variable("boundary_map") ;

    constraint geom_cells;
    geom_cells = facts.get_variable("geom_cells");
    constraint faces;
    faces = facts.get_variable("faces");

    //find actually face_dom
    entitySet face_dom;
    
    FORALL(*geom_cells, cc) {
      for(int i=0;i<upper[cc].size();++i) face_dom += upper[cc][i];
      for(int i=0;i<lower[cc].size();++i) face_dom += lower[cc][i];
      for(int i=0;i<boundary_map[cc].size();++i) face_dom += boundary_map[cc][i];
    }ENDFORALL;
    
    //expand map face2node
    if(Loci::MPI_processes > 1){
      vector<entitySet> init_ptn = facts.get_init_ptn() ;
      entitySet out_of_dom = face_dom - face2node.domain();
      face2node.setRep(MapRepP(face2node.Rep())->expand(out_of_dom, init_ptn));
    }

    entitySet hexcell;
    entitySet prism;
    entitySet gnrlcell;
    int elem_type;


  
     // Classify Cells
     FORALL(*geom_cells, cc) {
       int nfaces = upper[cc].size()+lower[cc].size()+boundary_map[cc].size() ;
       tmp_array<Entity> faces(nfaces) ;
       int cnt = 0 ;
       for(int i=0;i<upper[cc].size();++i)
         faces[cnt++] = upper[cc][i] ;
       for(int i=0;i<lower[cc].size();++i)
      faces[cnt++] = lower[cc][i] ;
       for(int i=0;i<boundary_map[cc].size();++i)
         faces[cnt++] = boundary_map[cc][i] ;
       elem_type = classify_cell(faces,nfaces,face2node) ;
       switch(elem_type) {
       case 1:
         hexcell += cc ; break ;
       case 2:
         prism += cc ; break ;
       default:
         gnrlcell += cc ;
       }
     } ENDFORALL ;
     
     constraint hexcells;
     *hexcells = hexcell;
     
     constraint prisms;
     *prisms = prism;
     
     constraint gnrlcells;
     *gnrlcells = gnrlcell;
     
     facts.create_fact("hexcells", hexcells);
     facts.create_fact("prisms", prisms);
     facts.create_fact("gnrlcells", gnrlcells);

     entitySet quadface;
     FORALL(*faces, ff){
       if(face2node.num_elems(ff)==4)quadface += ff;
     }ENDFORALL;

     
     constraint quadrangles;
     *quadrangles = quadface;
     facts.create_fact("quadrangles", quadrangles);
   //   cout <<"myID: " << Loci::MPI_rank <<  " hexcells : " << hexcell.size() << " prisms: " << prism.size()
//           << " gnrlcells : " << gnrlcell.size() << " quadface: " << quadface.size() << endl;
     
  }
}



// class cell_constraint : public constraint_rule {
//   const_multiMap lower;
//   const_multiMap upper;
//   const_multiMap boundary_map;
//   const_multiMap face2node;
//   const_store<vect3d> pos;
//   Constraint hexcells ;
//   Constraint prisms ;
//   Constraint gnrlcells ;
  
// public:
//   cell_constraint() {
//     name_store("lower", lower);
//     name_store("upper", upper);
//     name_store("boundary_map", boundary_map);
//     name_store("face2node", face2node);
//     name_store("pos", pos);
//     name_store("hexcells", hexcells);
//     name_store("prisms", prisms);
//     name_store("gnrlcells", gnrlcells);
//     input("(upper,lower,boundary_map)->face2node->pos");
//     output("(hexcells,prisms,gnrlcells)") ;
//     constraint("geom_cells");
//   }
//   void compute(const sequence& seq) {
//     // do_loop(seq, this);
//     // }
//     // void calculate(Entity cc) {

//     hexcells = EMPTY;
//     prisms = EMPTY;
//     gnrlcells = EMPTY;
//     int elem_type;
//     for(sequence::const_iterator cc = seq.begin(); cc != seq.end(); cc++){
//       int nfaces = upper[*cc].size()+lower[*cc].size()+boundary_map[*cc].size() ;
//       Loci::tmp_array<Entity> faces(nfaces) ;
//       int cnt = 0 ;
//       for(int i=0;i<upper[*cc].size();++i)
//         faces[cnt++] = upper[*cc][i] ;
//       for(int i=0;i<lower[*cc].size();++i)
//         faces[cnt++] = lower[*cc][i] ;
//       for(int i=0;i<boundary_map[*cc].size();++i)
//         faces[cnt++] = boundary_map[*cc][i] ;
//       elem_type = classify_cell(faces,nfaces,face2node) ;
//       switch(elem_type) {
//       case 1:
//         *hexcells += *cc ; break ;
//       case 2:
//         *prisms += *cc ; break ;
//       default:
//          *gnrlcells += *cc ;
//       }
//     }
//   }
//   } ;
// register_rule<cell_constraint> register_cell_constraint ;


// class face_constraint : public constraint_rule {
//   const_multiMap face2node;
//   const_store<vect3d> pos;
//   Constraint quadrangles ;
 
  
// public:
//   face_constraint() {
//     name_store("face2node", face2node);
//     name_store("pos", pos);
//     name_store("quadrangles", quadrangles);
//     input("face2node->pos");
//     output("quadrangles") ;
//     constraint("faces");
//   }
//   void compute(const sequence& seq) {
//     quadrangles = EMPTY;
//     do_loop(seq, this);
//   }
//   void calculate(Entity cc) {
//     if(face2node.num_elems(cc)==4)*quadrangles += cc;
//   }
// } ;
// register_rule<face_constraint> register_face_constraint ;
















//dummy rules to create constraint quadface

class quad_face1:public pointwise_rule{
 
  store<bool > quadfaces;

public:
  quad_face1(){
   
    name_store("quadfaces", quadfaces);
    output("quadfaces");
    constraint("boundary_faces, cl->hexcells");
  }
  virtual void compute(const sequence &seq){
    do_loop(seq, this);
   
  }
  void calculate(Entity f){}
};

register_rule<quad_face1> register_quad_face1;

class quad_face2:public pointwise_rule{
 
  store<bool > quadfaces;

public:
  quad_face2(){
   
    name_store("quadfaces", quadfaces);
    output("quadfaces");
    constraint("boundary_faces, quadrangles, cl->prisms");
  }
  virtual void compute(const sequence &seq){
    do_loop(seq, this);
   
  }
  void calculate(Entity f){}
};

register_rule<quad_face2> register_quad_face2;

class quad_face3:public pointwise_rule{
 
  store<bool > quadfaces;

public:
  quad_face3(){
   
    name_store("priority::cr::quadfaces", quadfaces);
    output("priority::cr::quadfaces");
    constraint("interior_faces, quadrangles,  cl->prisms");
  }
  virtual void compute(const sequence &seq){
    do_loop(seq, this);
   
  }
  void calculate(Entity f){}
};

register_rule<quad_face3> register_quad_face3;

class quad_face4:public pointwise_rule{
 
  store<bool > quadfaces;

public:
  quad_face4(){
   
    name_store("priority::cr::quadfaces", quadfaces);
    output("priority::cr::quadfaces");
    constraint("interior_faces, cl->hexcells");
  }
  virtual void compute(const sequence &seq){
    do_loop(seq, this);
   
  }
  void calculate(Entity f){}
};

register_rule<quad_face4> register_quad_face4;

class quad_face5:public pointwise_rule{
 
  store<bool > quadfaces;

public:
  quad_face5(){
   
    name_store("cr::quadfaces", quadfaces);
    output("cr::quadfaces");
    constraint("interior_faces, quadrangles, cr->prisms");
  }
  virtual void compute(const sequence &seq){
    do_loop(seq, this);
   
  }
  void calculate(Entity f){}
};

register_rule<quad_face5> register_quad_face5;

class quad_face6:public pointwise_rule{
 
  store<bool > quadfaces;

public:
  quad_face6(){
   
    name_store("cr::quadfaces", quadfaces);
    output("cr::quadfaces");
    constraint("interior_faces, cr->hexcells");
  }
  virtual void compute(const sequence &seq){
    do_loop(seq, this);
   
  }
  void calculate(Entity f){}
};

register_rule<quad_face6> register_quad_face6;


//class to check if a face is quadface
class identify_face:public pointwise_rule{
 
  store<bool > is_quadface;

public:
  identify_face(){
   
    name_store("quadfaces::is_quadface", is_quadface);
    output("quadfaces::is_quadface");
    constraint("faces");
  }
  virtual void compute(const sequence &seq){
    do_loop(seq, this);
   
  }
  void calculate(Entity f){
    is_quadface[f] = false;
  }
};

register_rule<identify_face> register_identify_face;


class identify_quadface:public pointwise_rule{
 
  store<bool > is_quadface;

public:
  identify_quadface(){
   
    name_store("priority::quadfaces::is_quadface", is_quadface);
    output("priority::quadfaces::is_quadface");
    constraint("quadfaces");
  }
  virtual void compute(const sequence &seq){
    do_loop(seq, this);
    //std::cout << "quadface in a rule: " << seq.size() << " myID " << Loci::MPI_rank << std::endl;
   
  }
  void calculate(Entity f){
    is_quadface[f] = true;
  }
};

register_rule<identify_quadface> register_identify_quadface;
