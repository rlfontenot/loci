#ifndef DATA_TRAITS_H_
#define DATA_TRAITS_H_

#include <vector>
#include <utility>
#include <entitySet.h>

extern "C" {
#include <hdf5.h>
#ifdef inline
#undef inline
#endif
}

#include <Loci_Datatypes.h>

namespace Loci {

  template <class T> 
  class data_schema_traits { };

  template<class T> class StringStreamConverter {
    T &ref ;
  public:
    StringStreamConverter(T &iref) : ref(iref) {}
    int getSize() const {
      std::ostringstream oss ;
      oss << ref ;
      return oss.str().size() ;
    }
    void getState(char *buf, int &size) {
      std::ostringstream oss ;
      oss << ref ;
      const std::string &s = oss.str() ;
      size = s.size() ;
      for(int i=0;i<size;++i)
        buf[i] = s[i] ;
    }
    void setState(char *buf, int size) {
      std::string s ;
      for(int i=0;i<size;++i)
        s += buf[i] ;
      std::istringstream iss(s) ;
      iss >> ref ;
    }
  } ;
    
  template<class T> class UndefinedConverter {
    T &ref ;
  public:
    UndefinedConverter(T &iref) : ref(iref) {}
    int getSize() const {
      std::cerr << "undefined converter" << std::endl ;
      return 0 ;
    }
    void getState(char *buf, int &size) {
      std::cerr << "undefined converter" << std::endl ;
    }
    void setState(char *buf, int size) {
      std::cerr << "undefined converter" << std::endl ;
    }
  } ;
    
    
  struct IDENTITY_CONVERTER{};
  struct USER_DEFINED_CONVERTER{};

  //**************************************************************************

  template <class T> 
  class data_schema_converter_traits { };

  template <>
  struct data_schema_traits<char> {
    typedef IDENTITY_CONVERTER Schema_Converter;

    static DatatypeP get_type() {
      return DatatypeP(new AtomicType(CHAR)) ;
    }
  };

  template <>
  struct data_schema_traits<unsigned char> {
    typedef IDENTITY_CONVERTER Schema_Converter;

    static DatatypeP get_type() {
      return DatatypeP(new AtomicType(UNSIGNED_CHAR)) ;
    }
  };

  template <>
  struct data_schema_traits<short> {
    typedef IDENTITY_CONVERTER Schema_Converter;

    static DatatypeP get_type() {
      return DatatypeP(new AtomicType(SHORT)) ;
    }
  };

  template <>
  struct data_schema_traits<int> {
    typedef IDENTITY_CONVERTER Schema_Converter;

    static DatatypeP get_type() {
      return DatatypeP(new AtomicType(INT)) ;
    }
  } ;
    
  template <>
  struct data_schema_traits<unsigned> {
    typedef IDENTITY_CONVERTER Schema_Converter;

    static DatatypeP get_type() {
      return DatatypeP(new AtomicType(UNSIGNED)) ;
    }
  };

  template <>
  struct data_schema_traits<unsigned short> {
    typedef IDENTITY_CONVERTER Schema_Converter;

    static DatatypeP get_type() {
      return DatatypeP(new AtomicType(UNSIGNED_SHORT)) ;
    }
  };

  template <>
  struct data_schema_traits<long> {
    typedef IDENTITY_CONVERTER Schema_Converter;

    static DatatypeP get_type() {
      return DatatypeP(new AtomicType(UNSIGNED)) ;
    }      
  };

  template <>
  struct data_schema_traits<unsigned long> {
    typedef IDENTITY_CONVERTER Schema_Converter;

    static DatatypeP get_type() {
      return DatatypeP(new AtomicType(UNSIGNED_LONG)) ;
    }      
  } ;

  template <>
  struct data_schema_traits<float> {
    typedef IDENTITY_CONVERTER Schema_Converter;

    static DatatypeP get_type() {
      return DatatypeP(new AtomicType(FLOAT)) ;
    }      
  };

  template <>
  struct data_schema_traits<double> {
    typedef IDENTITY_CONVERTER Schema_Converter;

    static DatatypeP get_type() {
      return DatatypeP(new AtomicType(DOUBLE)) ;
    }      
  };

  template <>
  struct data_schema_traits<long double> {
    typedef IDENTITY_CONVERTER Schema_Converter;

    static DatatypeP get_type() {
      return DatatypeP(new AtomicType(LONG_DOUBLE)) ;
    }      
  };

  template <>
  struct data_schema_traits<bool> {
    typedef IDENTITY_CONVERTER Schema_Converter;

    static DatatypeP get_type() {
      return DatatypeP(new AtomicType(BOOL)) ;
    }      
  };


  class entitySetSchemaConverter {
    // For the schema converter, we always store a reference to the object
    // we are converting schmata for.
    entitySet &eref ;
  public:
    explicit entitySetSchemaConverter(entitySet &iset): eref(iset) {}
    int getSize() const {
      return eref.size() ;
    }
    void getState(int *buf, int &size) {
      size = getSize() ;
      int ii=0; 
      for(entitySet::const_iterator i=eref.begin();i!=eref.end();++i)
        buf[ii++] = *i ;
    }
    void setState(int *buf, int size) {
      eref = EMPTY ;
      for(int i=0;i<size;++i)
        eref += buf[i] ;
    }
  } ;
  
  template <>
  struct data_schema_traits<entitySet> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;

    typedef int Converter_Base_Type ;
    typedef entitySetSchemaConverter Converter_Type ;
  };

  template<class T> class vectorSchemaConverter {
    // For the schema converter, we always store a reference to the object
    // we are converting schmata for.
    std::vector<T> &eref ;
  public:
    explicit vectorSchemaConverter(std::vector<T> &iset): eref(iset) {}
    int getSize() const {
      return eref.size() ;
    }
    void getState(T *buf, int &size) {
      size = getSize() ;
      for(int i=0;i<size;++i)
        buf[i] = eref[i] ;
    }
    void setState(T *buf, int size) {
      eref.clear() ;
      eref.reserve(size) ;
      for(int i=0;i<size;++i)
        eref.push_back(buf[i]) ;
    }
  } ;

  template<class T> struct data_schema_traits<std::vector<T> > {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;

    typedef T Converter_Base_Type ;
    typedef vectorSchemaConverter<T> Converter_Type ;
  } ;

  template<class T> inline std::ostream &operator << (std::ostream &s,
                                                      const std::vector<T> &v){
    s << v.size() ;
    for(std::vector<T>::const_iterator i=v.begin();i!=v.end(); ++i) 
      s << ' ' << *i ;
    s << std::endl ;
    return s ;
  }

  template<class T> inline std::istream &operator >> (std::istream &s,
                                                      std::vector<T> &v){
    v.clear() ;
    int sz ;
    s >> sz ;
    for(int i=0;i<sz;++i) {
      T val ;
      s >> val ;
      v.push_back(val) ;
    }
    return s ;
  }
  
  template<class T, class S> struct data_schema_traits<std::pair<T,S> > {
    typedef IDENTITY_CONVERTER Schema_Converter ;
    static DatatypeP get_type() {
      typedef std::pair<T,S> pair ;
      pair p ;
      CompoundDatatypeP ct = CompoundFactory(p) ;
      ct->insert("first",offsetof(pair,first), getLociType(p.first)) ;
      ct->insert("second",offsetof(pair,second), getLociType(p.second)) ;
      return DatatypeP(ct) ;
    }
  } ;

  template<class T> inline DatatypeP getLociType(T in) {
    return data_schema_traits<T>::get_type() ;
  }
    
}

#define LOCI_INSERT_TYPE(ct,type,variable) \
   ct->insert(# variable, offsetof(type,variable),getLociType(type().variable))
#endif
