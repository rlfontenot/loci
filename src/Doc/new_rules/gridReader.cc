// This is an example file for the new grid reader when use with
// the new default_rule, optional_rule, and constraint_rule. 
// This is NOT a final
// version. We are still working on this feature. As we make progress,
// I will add corresponding stuff.

#include "sciTypes.h"
#include "gridReader.h"
#include "coeff.h"
#include <string>
using std::string ;
#include <vector>
using std::vector ;

namespace heat {

  // The read_grid function now requires to have the
  // global rule_db as an argument. The rdb is used
  // in the read in process.
  void read_grid(fact_db &facts, const rule_db& rdb,
                 const char *filename) {
    char buf[512] ;
    sprintf(buf,"%s.vars",filename) ;
    ifstream ifile(buf,ios::in) ;
    if(ifile.fail()) {
      cerr<<"can't open " << buf << endl ;
      exit(-1) ;
    }
    
    //************************read vars file
    // NOTE: this is the new interface for read in .vars file.
    // the original one "read" is still available, but it will
    // not recognize the new "default_rule" and "optinoal_rule"
    facts.read_vars(ifile,rdb) ;
    //************************set variables

    param<string> modelName ;
    *modelName = filename ;
    facts.create_fact("modelName",modelName) ;

    // the struct "grid_file_info" and "bc_info" are
    // now moved into the header file "gridReader.h"
    param<grid_options> grid_file_info ;
    grid_file_info = facts.get_fact("grid_file_info") ;

    string file_type ;
    (*grid_file_info).getOption("file_type",file_type) ;

    if(file_type == "NEW_BINARY_FORMATTED" || file_type == "XDR") {
      string file = string(filename) + string(".xdr") ;
      if(Loci::MPI_rank == 0)
        cout << "Portable XDR Grid File input, reading file = "
             << file << endl ;
      if(!Loci::readFVMGrid(facts,file)) {
        cerr << "Yikes! Reading grid file '"
             << file <<"' failed in grid reader!" << endl ;
        cerr << "aborting!" << endl ;
        Loci::Abort() ;
      }
      if(Loci::MPI_rank == 0)
        cout << "Reading XDR Grid File Complete" << endl ;
    } else {
      cerr << "cannot handle grid " << file_type << " type for reading" << endl ;
    }

    double Lref = 1.0 ;
    if((*grid_file_info).optionExists("Lref")) {
      Loci::option_value_type ovt = (*grid_file_info).getOptionValueType("Lref") ;
      if(ovt == Loci::REAL) {
        (*grid_file_info).getOption("Lref",Lref) ;
      }
      if(ovt == Loci::UNIT_VALUE) {
        Loci::UNIT_type vu ;
        (*grid_file_info).getOption("Lref",vu) ;
        if(!vu.is_compatible("m")) {
          cerr << "wrong type of units for reference length Lref" << endl ;
        } else {
          Lref = vu.get_value_in("m") ;
        }
      }
      store<vect3d> pos ;
      pos = facts.get_fact("pos") ;
      entitySet posdom = pos.domain() ;
      if(Loci::MPI_rank == 0)
        cout << "Scaling " << filename << " grid by " << Lref << endl ;
      for(entitySet::const_iterator ei =posdom.begin();ei!=posdom.end();++ei){
        pos[*ei] *= Lref ;
      }
    }

    setup_faces(facts) ;

    setup_boundary_conditions(facts) ;

    color_matrix(facts, COLOR_DFS) ;
    make_faces_consistent(facts) ;
    create_lower_upper(facts) ;
  }
}

// the end
