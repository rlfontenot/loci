#ifndef DATA_TRAITS_H_
#define DATA_TRAITS_H_

#include <vector>

#include <mpi.h>

extern "C" {
#include <hdf5.h>
#ifdef inline
#undef inline
#endif
}

#include <Loci_Datatypes.h>

namespace Loci {

#ifdef ALLOW_DEFAULT_CONVERTER
  struct DEFAULT_CONVERTER{};
  template <class T> 
  class data_schema_traits { 
  public:
    typedef DEFAULT_CONVERTER Schema_Converter;
    static MPI_Datatype  get_mpi_type()  {return MPI_CHAR;}
    static hid_t         get_hdf5_type() {return H5T_NATIVE_CHAR;}
  };
#else
  template <class T> 
  class data_schema_traits { };
#endif

  struct IDENTITY_CONVERTER{};
  struct USER_DEFINED_CONVERTER{};

  //**************************************************************************

  template <class T> 
  class data_schema_converter_traits { };

  template <>
  struct data_schema_traits<char> {
    typedef IDENTITY_CONVERTER Schema_Converter;

    static AtomicType*   instance()        
                         { AtomicType *atom = new AtomicType(CHAR);
                           return atom;
                         }
    static hid_t         get_hdf5_type()   {return H5T_NATIVE_CHAR;}
    static Datatype      get_loci_type()   {return ATOMIC;}
    static MPI_Datatype  get_mpi_type()    {return MPI_CHAR;}
  };

  template <>
  struct data_schema_traits<short> {
    typedef IDENTITY_CONVERTER Schema_Converter;

    static AtomicType*   instance()        
                         { AtomicType *atom = new AtomicType(SHORT);
                           return atom;
                         }
    static hid_t         get_hdf5_type()   {return H5T_NATIVE_SHORT;}
    static Datatype      get_loci_type()   {return ATOMIC;}
    static MPI_Datatype  get_mpi_type()    {return MPI_SHORT;}
  };

  template <>
  struct data_schema_traits<int> {
    typedef IDENTITY_CONVERTER Schema_Converter;

    static AtomicType*   instance()  
                         { AtomicType *atom = new AtomicType(INT);
                           return atom;
                         }
    static hid_t         get_hdf5_type()   {return H5T_NATIVE_INT;}
    static Datatype      get_loci_type()   {return ATOMIC;}
    static MPI_Datatype  get_mpi_type()    {return MPI_INT;}
  };

  template <>
  struct data_schema_traits<unsigned> {
    typedef IDENTITY_CONVERTER Schema_Converter;

    static AtomicType*   instance()        
                         { AtomicType *atom = new AtomicType(UNSIGNED);
                           return atom;
                         }
    static hid_t         get_hdf5_type()   {return H5T_NATIVE_UINT;}
    static Datatype      get_loci_type()   {return ATOMIC;}
    static MPI_Datatype  get_mpi_type()    {return MPI_UNSIGNED;}
  };

  template <>
  struct data_schema_traits<long> {
    typedef IDENTITY_CONVERTER Schema_Converter;

    static AtomicType*   instance()        
                         { AtomicType *atom = new AtomicType(LONG);
                           return atom;
                         }
    static hid_t         get_hdf5_type()   {return H5T_NATIVE_LONG;}
    static Datatype      get_loci_type()   {return ATOMIC;}
    static MPI_Datatype  get_mpi_type()    {return MPI_LONG;}
  };

  template <>
  struct data_schema_traits<unsigned long> {
    typedef IDENTITY_CONVERTER Schema_Converter;

    static AtomicType*   instance()        
                         { AtomicType *atom = new AtomicType(UNSIGNED_LONG);
                           return atom;
                         }
    static hid_t         get_hdf5_type()   {return H5T_NATIVE_ULONG;}
    static Datatype      get_loci_type()   {return ATOMIC;}
    static MPI_Datatype  get_mpi_type()    {return MPI_UNSIGNED_LONG;}
  };

  template <>
  struct data_schema_traits<float> {
    typedef IDENTITY_CONVERTER Schema_Converter;

    static AtomicType*   instance()        
                         { AtomicType *atom = new AtomicType(FLOAT);
                           return atom;
                         }
    static hid_t         get_hdf5_type()   {return H5T_NATIVE_FLOAT;}
    static Datatype      get_loci_type()   {return ATOMIC;}
    static MPI_Datatype  get_mpi_type()    {return MPI_FLOAT;}
  };

  template <>
  struct data_schema_traits<double> {
    typedef IDENTITY_CONVERTER Schema_Converter;

    static AtomicType*   instance()        
                         { AtomicType *atom = new AtomicType(DOUBLE);
                           return atom;
                         }
    static hid_t         get_hdf5_type()   {return H5T_NATIVE_DOUBLE;}
    static Datatype      get_loci_type()   {return ATOMIC;}
    static MPI_Datatype  get_mpi_type()    {return MPI_DOUBLE;}
  };

  template <>
  struct data_schema_traits<long double> {
    typedef IDENTITY_CONVERTER Schema_Converter;

    static AtomicType*   instance()        
                         { AtomicType *atom = new AtomicType(LONG_DOUBLE);
                           return atom;
                         }
    static hid_t         get_hdf5_type()   {return H5T_NATIVE_LDOUBLE;}
    static Datatype      get_loci_type()   {return ATOMIC;}
    static MPI_Datatype  get_mpi_type()    {return MPI_LONG_DOUBLE;}
  };

  template <>
  struct data_schema_traits<bool> {
    typedef IDENTITY_CONVERTER Schema_Converter;

    static AtomicType*   instance()
                         { AtomicType *atom = new AtomicType(BOOL);
                           return atom;
                         }
    static hid_t         get_hdf5_type()   {return H5T_NATIVE_HBOOL;}
    static Datatype      get_loci_type()   {return ATOMIC;}
  };

  template <>
  struct data_schema_traits<entitySet> {
    typedef IDENTITY_CONVERTER Schema_Converter;

    static AtomicType* instance() 
                       { AtomicType *atom = new AtomicType(BOOL);
                           return atom;
                       }
    static hid_t       get_hdf5_type() {return H5T_NATIVE_INT;}
    static Datatype    get_loci_type() {return ATOMIC;}
  };

}

#endif
