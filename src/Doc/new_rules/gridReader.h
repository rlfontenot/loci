#ifndef READ_GRID_H
#define READ_GRID_H

#include <Loci.h>

namespace heat {
  void read_grid(fact_db &facts, const rule_db& rdb,
                 const char *filename);

  class grid_options : public options_list {
  public:
    grid_options() :
      options_list("file_type:Lref") {} ;
  } ;
  
  class bc_options : public options_list {
  public:
    bc_options() :
      options_list("") {} ;
  } ;
}

namespace Loci {
  template<> struct data_schema_traits<heat::grid_options> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    
    typedef char Converter_Base_Type ;
    typedef StringStreamConverter<heat::grid_options> Converter_Type ;
  } ;
  
}

namespace heat {
  void setup_boundary_conditions(fact_db &facts) ;

  enum matrix_coloring_type {COLOR_DEFAULT, COLOR_DFS} ;

  void setup_faces(fact_db &facts) ;
  void color_matrix(fact_db &facts, matrix_coloring_type mct) ;

  void create_ghost_cells(fact_db &facts) ;
  void create_ci_map(fact_db &facts) ;
  
  void make_faces_consistent(fact_db &facts) ;

  void create_lower_upper(fact_db &facts) ;
}

namespace Loci {
  template<> struct data_schema_traits<heat::bc_options> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    
    typedef char Converter_Base_Type ;
    typedef StringStreamConverter<heat::bc_options> Converter_Type ;
  } ;
}

#endif
