#include <Loci.h>
#include <map>
#include "defines.h"








namespace Loci{
  entitySet findBoundingSet(entitySet dom);
  
  storeRepP  my_get_node_remap(entitySet nodes, entitySet faces, entitySet edges) {
    
    entitySet dom = nodes + faces + edges;
    if(MPI_processes == 1) {
      int minNode = nodes.Min() ;
      Map nm ;
     
      nm.allocate(dom) ;
      
      FORALL(nodes,nd) {
        nm[nd] = nd - minNode +1;;
       } ENDFORALL ;
      
      FORALL(faces,nd) {
        nm[nd] = nd;
      } ENDFORALL ;
      
      FORALL(edges,nd) {
        nm[nd] = nd;
      } ENDFORALL ;
      
      return nm.Rep() ;
    }
    
    std::vector<entitySet> init_ptn = Loci::exec_current_fact_db->get_init_ptn() ;
    fact_db::distribute_infoP df = Loci::exec_current_fact_db->get_distribute_info() ;
    Map l2g ;
    l2g = df->l2g.Rep() ;
    dMap g2f ;
    g2f = df->g2f.Rep() ;
    
    entitySet gnodes = l2g.image(nodes&l2g.domain()) ;
    entitySet gset = findBoundingSet(gnodes) ;
    int minNode = gset.Min();
    
    entitySet gelements =l2g.image(dom&l2g.domain()) ; 
    Map newnum ;
    newnum.allocate(dom) ;
    
    // Expand g2f to include clone regions
    entitySet out_of_dom = gelements - init_ptn[MPI_rank] ;
    g2f.setRep(MapRepP(g2f.Rep())->expand(out_of_dom, init_ptn)) ;
    
    FORALL(nodes,i) {
      newnum[i] = g2f[l2g[i]] - minNode +1;
    } ENDFORALL ;
    
    FORALL(faces,i) {
      newnum[i] = g2f[l2g[i]];
    } ENDFORALL ;
    FORALL(edges,i) {
      newnum[i] = g2f[l2g[i]];
    } ENDFORALL ;
    return newnum.Rep() ;
  }
}













// Note that this is a unit_rule, since this is the
// only way we can have stores as input to a rule outputting blackboxes.
class set_cell_remap_unit : public unit_rule {
private:
  const_store<vect3d> pos;
  const_multiMap upper;
  const_multiMap lower;
  const_multiMap boundary_map;
  const_MapVec<2> edge2node;
  const_multiMap face2edge;
  const_multiMap face2node;
  blackbox<Loci::storeRepP > node_remap ;
public:

      // Define input and output.
  set_cell_remap_unit() {
    name_store("pos", pos);
    name_store("lower", lower);
    name_store("upper", upper);
    name_store("boundary_map", boundary_map);
    name_store("face2node", face2node);
    name_store("face2edge", face2edge);
    name_store("edge2node", edge2node);
    name_store("node_remap", node_remap);
    
    input("(lower, upper, boundary_map)->face2node->pos");
    input("(lower, upper, boundary_map)->face2edge->edge2node->pos");
    output("node_remap");
    constraint("geom_cells") ;
    disable_threading() ;
  }
  
  // Do the set-up.
  void compute(const sequence & seq) {
   *node_remap = my_get_node_remap(pos.domain(), face2node.domain(), edge2node.domain());
   
  }
  
};
register_rule<set_cell_remap_unit> register_set_cell_remap_unit ;

 // Empty apply rule required by Loci. The data type and operator do not
  // matter since nothing is done by this rule. Keep the same inputs and
  // outputs as the unit rule, even though we don't have to.
class set_cell_remap_apply : public apply_rule<blackbox<Map>,
                             Loci::NullOp<Map> > {
private:
  const_store<vect3d> pos;
  const_multiMap upper;
  const_multiMap lower;
  const_multiMap boundary_map;
  const_MapVec<2> edge2node;
  const_multiMap face2edge;
  const_multiMap face2node;
  blackbox<Loci::storeRepP > node_remap ;
public:
  
  // Define input and output.
  set_cell_remap_apply() {
    name_store("pos", pos);
    name_store("lower", lower);
    name_store("upper", upper);
    name_store("boundary_map", boundary_map);
    name_store("face2node", face2node);
    name_store("face2edge", face2edge);
    name_store("edge2node", edge2node);
    name_store("node_remap", node_remap);
    input("(lower, upper, boundary_map)->face2node->pos");
    input("(lower, upper, boundary_map)->face2edge->edge2node->pos");
    output("node_remap");
    constraint("geom_cells") ;
    disable_threading() ;
  }
  
  // Do nothing.
  void compute(const sequence & seq) {}
  } ;

  register_rule<set_cell_remap_apply> registerset_cell_remap_apply ;


class set_interior_face_remap_unit : public unit_rule {
private:
  const_store<vect3d> pos;
  const_multiMap upper;
  const_multiMap lower;
  const_multiMap boundary_map;
  
  const_MapVec<2> edge2node;
  const_multiMap face2edge;
  const_multiMap face2node;
   const_Map cl;
  const_Map cr;
  blackbox<Loci::storeRepP > node_remap ;
public:

      // Define input and output.
  set_interior_face_remap_unit() {
    name_store("pos", pos);
    name_store("lower", lower);
    name_store("upper", upper);
    name_store("boundary_map", boundary_map);
    name_store("face2node", face2node);
    name_store("face2edge", face2edge);
    name_store("edge2node", edge2node);
    name_store("iface_remap", node_remap);
    name_store("cl", cl);
    name_store("cr", cr);
    input("(cl, cr)->(lower, upper, boundary_map)->face2node->pos");
    input("(cl, cr)->(lower, upper, boundary_map)->face2edge->edge2node->pos");
    input("face2node->pos");
    output("iface_remap");
    constraint("interior_faces") ;
    disable_threading() ;
  }
  
  // Do the set-up.
  void compute(const sequence & seq) {
   *node_remap = my_get_node_remap(pos.domain(), face2node.domain(), edge2node.domain());
   
  }
  
};
register_rule<set_interior_face_remap_unit> register_set_interior_face_remap_unit ;

 // Empty apply rule required by Loci. The data type and operator do not
  // matter since nothing is done by this rule. Keep the same inputs and
  // outputs as the unit rule, even though we don't have to.
class set_interior_face_remap_apply : public apply_rule<blackbox<Map>,
                             Loci::NullOp<Map> > {
private:
  const_store<vect3d> pos;
  const_multiMap upper;
  const_multiMap lower;
  const_multiMap boundary_map;
  
  const_MapVec<2> edge2node;
  const_multiMap face2edge;
  const_multiMap face2node;
   const_Map cl;
  const_Map cr;
  blackbox<Loci::storeRepP > node_remap ;
public:
  
  // Define input and output.
  set_interior_face_remap_apply() {

    name_store("pos", pos);
    name_store("lower", lower);
    name_store("upper", upper);
    name_store("boundary_map", boundary_map);
    name_store("face2node", face2node);
    name_store("face2edge", face2edge);
    name_store("edge2node", edge2node);
    name_store("iface_remap", node_remap);
    name_store("cl", cl);
    name_store("cr", cr);
    input("(cl, cr)->(lower, upper, boundary_map)->face2node->pos");
    input("(cl, cr)->(lower, upper, boundary_map)->face2edge->edge2node->pos");
    input("face2node->pos");
    output("iface_remap");
    constraint("interior_faces") ;
    disable_threading() ;
  }
    
  
  // Do nothing.
  void compute(const sequence & seq) {}
  } ;

register_rule<set_interior_face_remap_apply> register_set_interior_face_remap_apply ;

class set_boundary_face_remap_unit : public unit_rule {
private:
  const_store<vect3d> pos;
  const_multiMap upper;
  const_multiMap lower;
  const_multiMap boundary_map;
  
  const_MapVec<2> edge2node;
  const_multiMap face2edge;
  const_multiMap face2node;
  const_Map cl;
  blackbox<Loci::storeRepP > node_remap ;
public:

      // Define input and output.
  set_boundary_face_remap_unit() {
    name_store("pos", pos);
    name_store("lower", lower);
    name_store("upper", upper);
    name_store("boundary_map", boundary_map);
    name_store("face2node", face2node);
    name_store("face2edge", face2edge);
    name_store("edge2node", edge2node);
    name_store("bface_remap", node_remap);
    name_store("cl", cl);
  
    input("cl->(lower, upper, boundary_map)->face2node->pos");
    input("cl->(lower, upper, boundary_map)->face2edge->edge2node->pos");
    input("face2node->pos");
    output("bface_remap");
    constraint("boundary_faces") ;
    disable_threading() ;
  }
  
  // Do the set-up.
  void compute(const sequence & seq) {
   *node_remap = my_get_node_remap(pos.domain(), face2node.domain(), edge2node.domain());
   
  }
  
};
register_rule<set_boundary_face_remap_unit> register_set_boundary_face_remap_unit ;

 // Empty apply rule required by Loci. The data type and operator do not
  // matter since nothing is done by this rule. Keep the same inputs and
  // outputs as the unit rule, even though we don't have to.
class set_boundary_face_remap_apply : public apply_rule<blackbox<Map>,
                             Loci::NullOp<Map> > {
private:
  const_store<vect3d> pos;
  const_multiMap upper;
  const_multiMap lower;
  const_multiMap boundary_map;
  
  const_MapVec<2> edge2node;
  const_multiMap face2edge;
  const_multiMap face2node;
   const_Map cl;
 
  blackbox<Loci::storeRepP > node_remap ;
public:
  
  // Define input and output.
  set_boundary_face_remap_apply() {

    name_store("pos", pos);
    name_store("lower", lower);
    name_store("upper", upper);
    name_store("boundary_map", boundary_map);
    name_store("face2node", face2node);
    name_store("face2edge", face2edge);
    name_store("edge2node", edge2node);
    name_store("bface_remap", node_remap);
    name_store("cl", cl);
 
    input("cl->(lower, upper, boundary_map)->face2node->pos");
    input("cl->(lower, upper, boundary_map)->face2edge->edge2node->pos");
    input("face2node->pos");
    output("bface_remap");
    constraint("boundary_faces") ;
    disable_threading() ;
  }
    
  
  // Do nothing.
  void compute(const sequence & seq) {}
  } ;

register_rule<set_boundary_face_remap_apply> register_set_boundary_face_remap_apply ;

//for fl
class set_face_remap_unit : public unit_rule {
private:
  const_store<vect3d> pos;
  const_multiMap upper;
  const_multiMap lower;
  const_multiMap boundary_map;
  
  const_MapVec<2> edge2node;
  const_multiMap face2edge;
  const_multiMap face2node;
  const_Map cl;
  blackbox<Loci::storeRepP > node_remap ;
public:

      // Define input and output.
  set_face_remap_unit() {
    name_store("pos", pos);
    name_store("lower", lower);
    name_store("upper", upper);
    name_store("boundary_map", boundary_map);
    name_store("face2node", face2node);
    name_store("face2edge", face2edge);
    name_store("edge2node", edge2node);
    name_store("face_remap", node_remap);
    name_store("cl", cl);
  
    input("cl->(lower, upper, boundary_map)->face2node->pos");
    input("cl->(lower, upper, boundary_map)->face2edge->edge2node->pos");
    input("face2node->pos");
    output("face_remap");
    constraint("faces") ;
    disable_threading() ;
  }
  
  // Do the set-up.
  void compute(const sequence & seq) {
   *node_remap = my_get_node_remap(pos.domain(), face2node.domain(), edge2node.domain());
   
  }
  
};
register_rule<set_face_remap_unit> register_set_face_remap_unit ;

 // Empty apply rule required by Loci. The data type and operator do not
  // matter since nothing is done by this rule. Keep the same inputs and
  // outputs as the unit rule, even though we don't have to.
class set_face_remap_apply : public apply_rule<blackbox<Map>,
                             Loci::NullOp<Map> > {
private:
  const_store<vect3d> pos;
  const_multiMap upper;
  const_multiMap lower;
  const_multiMap boundary_map;
  
  const_MapVec<2> edge2node;
  const_multiMap face2edge;
  const_multiMap face2node;
   const_Map cl;
 
  blackbox<Loci::storeRepP > node_remap ;
public:
  
  // Define input and output.
  set_face_remap_apply() {

    name_store("pos", pos);
    name_store("lower", lower);
    name_store("upper", upper);
    name_store("boundary_map", boundary_map);
    name_store("face2node", face2node);
    name_store("face2edge", face2edge);
    name_store("edge2node", edge2node);
    name_store("face_remap", node_remap);
    name_store("cl", cl);
 
    input("cl->(lower, upper, boundary_map)->face2node->pos");
    input("cl->(lower, upper, boundary_map)->face2edge->edge2node->pos");
    input("face2node->pos");
    output("face_remap");
    constraint("faces") ;
    disable_threading() ;
  }
    
  
  // Do nothing.
  void compute(const sequence & seq) {}
  } ;

register_rule<set_face_remap_apply> register_set_face_remap_apply ;
