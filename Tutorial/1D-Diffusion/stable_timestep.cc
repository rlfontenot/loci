#include <Loci.h>

#include <Tools/tools.h>

//*******************************************************************
// Compute maximum stable timestep for simulation.  Use reduction rule
// and calculate timestep as a function of local conditions.  The global
// timestep is the minimum of these local timesteps
//*******************************************************************

class timestep_unit : public unit_rule {
  param<float> dt ;
public:
  timestep_unit() {
    name_store("dt",dt) ;
    output("dt") ;
    constraint("UNIVERSE") ;  // This property applies to all entities
  }
  virtual void compute(const sequence &seq) {
    *dt = 1e30 ;              // Largest allowable timestep
  }
} ;

register_rule<timestep_unit> register_timestep_unit ;

//*******************************************************************

class timestep_apply : public apply_rule<param<float>,
                       Loci::Minimum<float> > {
  const_store<float> xc ;
  const_Map cl,cr ;
  const_param<float> nu ;
  param<float> dt ;
public:
  timestep_apply() {
    name_store("xc",xc) ;
    name_store("cl",cl) ;
    name_store("cr",cr) ;
    name_store("nu",nu) ;
    name_store("dt",dt) ;
    input("dt,(cl,cr)->xc,nu") ;
    output("dt") ;
  }
  void calculate(Entity e) {
    float dx = abs(xc[cr[e]]-xc[cl[e]]) ;
    // Compute timestep as 1/2 of maximum stable timestep
    float local_dt = dx*dx*nu[e]/4. ;  

    join(dt[e],local_dt) ; // Set dt = min(dt,local_dt)
  }

  virtual void compute(const sequence &seq) {
    do_loop(seq,this,&timestep_apply::calculate) ;
  }
  
} ;

register_rule<timestep_apply> register_timestep_apply ;

//*******************************************************************
