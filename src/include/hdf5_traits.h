#ifndef HDF5_TRAITS_H_
#define HDF5_TRAITS_H_

#include <hdf5CC/H5cpp.h>
#include <vector>
#include <complex>

#include <mpi.h>

namespace Loci {

  //list of selectors
  struct DEFAULT_CONVERTER{};
  struct IDENTITY_CONVERTER{};
  struct USER_DEFINED_CONVERTER{};

  template <class T> 
    class hdf5_schema_traits {
    public:
    typedef DEFAULT_CONVERTER Schema_Converter;
    static H5::DataType get_type() {return H5::PredType::NATIVE_CHAR;}
    static MPI_Datatype get_mpi_type() {return MPI_CHAR;}
  };
 
  template <>
    class hdf5_schema_traits<short> {
    public:
    typedef IDENTITY_CONVERTER Schema_Converter;
    static H5::DataType get_type() {return H5::PredType::NATIVE_SHORT;}
    static MPI_Datatype get_mpi_type() {return MPI_SHORT;}
  };

  template <>
    class hdf5_schema_traits<unsigned short> {
    public:
    typedef IDENTITY_CONVERTER Schema_Converter;
    static H5::DataType get_type() {return H5::PredType::NATIVE_USHORT;}
    static MPI_Datatype get_mpi_type() {return MPI_UNSIGNED_CHAR;}
  };

  template <>
    class hdf5_schema_traits<int> {
    public:
    typedef IDENTITY_CONVERTER Schema_Converter;
    static H5::DataType get_type() {return H5::PredType::NATIVE_INT;}
    static MPI_Datatype get_mpi_type() {return MPI_INT;}
  };
  
  template <>
    class hdf5_schema_traits<unsigned> {
    public:
    typedef IDENTITY_CONVERTER Schema_Converter;
    static H5::DataType get_type() {return H5::PredType::NATIVE_UINT;}
    static MPI_Datatype get_mpi_type() {return MPI_UNSIGNED;}
  };

  template <>
    class hdf5_schema_traits<long> {
    public:
    typedef IDENTITY_CONVERTER Schema_Converter;
    static H5::DataType get_type() {return H5::PredType::NATIVE_LONG;}
    static MPI_Datatype get_mpi_type() {return MPI_LONG;}
  };

  template <>
    class hdf5_schema_traits<unsigned long> {
    public:
    typedef IDENTITY_CONVERTER Schema_Converter;
    static H5::DataType get_type() {return H5::PredType::NATIVE_ULONG;}
    static MPI_Datatype get_mpi_type() {return MPI_UNSIGNED_LONG;}
  };

  template <>
    class hdf5_schema_traits<long long> {
    public:
    typedef IDENTITY_CONVERTER Schema_Converter;
    static H5::DataType get_type() {return H5::PredType::NATIVE_LLONG;}
  };

   template <>
    class hdf5_schema_traits<unsigned long long> {
    public:
    typedef IDENTITY_CONVERTER Schema_Converter;
    static H5::DataType get_type() {return H5::PredType::NATIVE_ULLONG;}
  };

  template <>
  class hdf5_schema_traits<double> {
  public:
    typedef IDENTITY_CONVERTER Schema_Converter;
    static H5::DataType get_type() {
      return H5::PredType::NATIVE_DOUBLE;
    }
    static MPI_Datatype get_mpi_type() {return MPI_DOUBLE;}
  };

  template <>
  class hdf5_schema_traits<float> {
  public:
    typedef IDENTITY_CONVERTER Schema_Converter;
    static H5::DataType get_type() {return H5::PredType::NATIVE_FLOAT;}
    static MPI_Datatype get_mpi_type() {return MPI_FLOAT;}
  };

  template <>
  class hdf5_schema_traits<long double> {
  public:
    typedef IDENTITY_CONVERTER Schema_Converter;
    static H5::DataType get_type() {return H5::PredType::NATIVE_LDOUBLE;}
    static MPI_Datatype get_mpi_type() {return MPI_LONG_DOUBLE;}
  };

  template <>
  class hdf5_schema_traits<bool> {
  public:
    typedef IDENTITY_CONVERTER Schema_Converter;
    static H5::DataType get_type() {return H5::PredType::NATIVE_HBOOL;}
  };

  template <class T1,class T2>
  class hdf5_schema_traits< std::pair<T1,T2> > {
  public:
    typedef IDENTITY_CONVERTER Schema_Converter;
    static H5::DataType get_type() {
      typedef std::pair<T1,T2> par;
      hdf5_schema_traits<T1> hdfT1;
      hdf5_schema_traits<T2> hdfT2;

      H5::CompType ctype(sizeof(par));

      ctype.insertMember("first", HOFFSET(par,first), hdfT1.get_type());
      ctype.insertMember("second", HOFFSET(par,second), hdfT2.get_type());

      return ctype;
      }
  };


  //***************************************************************************
  // Traits for STL Vectors....
  //***************************************************************************

  //selector which tell us if we need variable data
  struct Need_Variable_data{};
  struct No_Need_Variable_data{};

  //***************************************************************************
  // STL Vector...
  //***************************************************************************

  template <class T>
  class hdf5_schema_traits< std::vector<T> > {
  public:
    typedef USER_DEFINED_CONVERTER Schema_Converter;
  };

  template<class T>
  class hdf5_schema_converter_traits 
  { };

  template <class T> 
  class hdf5_schema_converter_traits< std::vector<T> > 
  {
    public:
      typedef T memento_type;
      static const int var_length = 1;
      static const int var_size   = 0;

      static H5::DataType get_variable_HDF5_type() 
      {
          hdf5_schema_traits<T> hdfT;
          return hdfT.get_type();
      };

      static MPI_Datatype get_variable_MPI_type() 
      {
          hdf5_schema_traits<T>  traits;
          return traits.get_mpi_type();
      };
  };

  //***************************************************************************
  // Traits for complex numbers ...
  //***************************************************************************

  template <class T>
  class hdf5_schema_traits< std::complex<T> > {
  public:
    typedef USER_DEFINED_CONVERTER Schema_Converter;
  };

  template <class T> 
     class hdf5_schema_converter_traits< std::complex<T>  > {
     public:
     typedef Need_Variable_data need_variable_selector;
     typedef T memento_type;

     static H5::DataType get_variable_HDF5_type() {
       hdf5_schema_traits<T> hdfT;
       return hdfT.get_type();
     };

     static MPI_Datatype get_variable_MPI_type() {
       hdf5_schema_traits<T>  traits;
       return traits.get_mpi_type();
     };

  };

  //***************************************************************************
}

template<class T>
inline std::ostream& operator << (std::ostream & s, const std::vector<T> &vec) 
{
    for( int i = 0; i < vec.size(); i++)
         s << vec[i] << " ";
    return s;
}

template<class T>
inline std:: istream& operator >> (std::istream & s, std::vector<T> &vec) {
    for( int i = 0; i < vec.size(); i++)
         s >> vec[i];
    return s;
 }

#endif
