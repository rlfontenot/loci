#ifndef LOCI_DATATYPE_H
#define LOCI_DATATYPE_H

#include <set>
#include <string>
#include <vector>
#include <sstream>

extern "C" {
#include <hdf5.h>
#ifdef inline
#undef inline
#endif
}

namespace Loci
{
  enum Datatype { ATOMIC, ARRAY, COMPOUND };

  enum AtomType {
    BOOL=1, CHAR, UNSIGNED_CHAR, BYTE,
    SHORT, UNSIGNED_SHORT, INT, UNSIGNED, LONG,
    UNSIGNED_LONG, FLOAT, DOUBLE, LONG_DOUBLE,
  };

  class CompoundType;
  class ArrayType;

  struct AbstractDatatype
  { 
     virtual int    getSize() const = 0;
     virtual hid_t  get_hdf5_type() const = 0;
  };

  class AtomicType : public AbstractDatatype
  {
    public:
     AtomicType(){}
     AtomicType(const AtomType &a) {atom = a;}

     AtomicType& operator = (const AtomType &a) 
                            { atom = a; return *this;}
     int          getSize() const ;
     hid_t        get_hdf5_type() const;
    private:
     AtomType     atom;
  };

  struct ArrayType 
  {
       Datatype    type;
       AtomicType  atom;
       int         rank; 
       int         dimension[10];
       int         getSize() const;
  };

  struct field
  {
      Datatype       type;
      std::string    name;
      size_t         offset;
      size_t         ncount;
      AtomicType     atom;
      ArrayType      array;
      CompoundType  *compound;
  };    


  class CompoundType : public AbstractDatatype
  {
  public:
    // Constructors ...
    CompoundType(){}

    // Assignment  ...
    CompoundType &operator = (const CompoundType &c);

    field &operator[] (int i) {return newtype[i];}

    // Member function ...
    void    insert( const std::string &name, size_t  offset, AtomType d);
    void    insert( const std::string &name, size_t  offset, ArrayType  a, AtomType d);
    void    insert( const std::string &name, size_t  offset, ArrayType  a, CompoundType &c);
    void    insert( const std::string &name, size_t  offset, CompoundType &c);

    std::string name( int i ) const;
    size_t      offset(int i ) const;
    AtomicType  type( int i ) const;

    void    setSize(const size_t &buf);
    int     getNumFields();

    int     getSize() const ;
    hid_t   get_hdf5_type() const;

    void    clear();
  private:
    size_t                 numBytes;
    std::vector<field>     newtype;      
    std::set<std::string>  name_used;
  };

  inline CompoundType & CompoundType :: operator = (const CompoundType &rhs)
  {
    if( this == &rhs) return *this;

    clear();
    name_used = rhs.name_used;

    std::vector<field> :: const_iterator iter;

    field  newfield;
    for( iter = rhs.newtype.begin(); iter != rhs.newtype.end(); ++iter){
      numBytes          = rhs.numBytes;
      newfield.name     = iter->name;
      newfield.atom     = iter->atom;
      newfield.array    = iter->array;
      newfield.compound = iter->compound;
      newfield.offset   = iter->offset;
      newtype.push_back(newfield);
    }

    return (*this);
  }

  inline void CompoundType :: clear()
  {
    name_used.clear();
  }

  inline void CompoundType :: insert( const std::string &name, size_t offset, AtomType d)
  {
    if( name_used.find(name) != name_used.end() ) return;

    field   newfield;
    name_used.insert(name);

    newfield.type     = ATOMIC;
    newfield.name     = name;
    newfield.atom     = d;
    newfield.offset   = offset;
    newtype.push_back(newfield);
  }

  inline void CompoundType :: insert( const std::string &name, size_t offset, ArrayType a,
                                      AtomType atom)
  {
    if( name_used.find(name) != name_used.end() ) return;

    field   newfield;
    name_used.insert(name);

    newfield.type        = ARRAY;
    newfield.name        = name;
    newfield.offset      = offset;
    newfield.array       = a;
    newfield.array.atom  = atom;
    newfield.array.type  = ATOMIC;
    newtype.push_back(newfield);

  }

  //***************************************************************************

  inline void CompoundType :: insert( const std::string &name, size_t offset, 
                                      ArrayType a, CompoundType &cmpd)
  {
    if( name_used.find(name) != name_used.end() ) return;

    field   newfield;
    name_used.insert(name);

    newfield.type        = ARRAY;
    newfield.name        = name;
    newfield.offset      = offset;
    newfield.array       = a;
    newfield.array.type  = COMPOUND;
    newfield.compound    = new CompoundType(cmpd);
    newtype.push_back(newfield);

  }

  //***************************************************************************
  inline void CompoundType :: insert( const std::string &name, size_t offset, 
                                      CompoundType &cmpd)
  {
    if( name_used.find(name) != name_used.end() ) return;
    
    field   newfield;
    name_used.insert(name);

    newfield.type        = COMPOUND;
    newfield.name        = name;
    newfield.offset      = offset;
    newfield.compound    = new CompoundType(cmpd);
    newtype.push_back(newfield);

  }
  //***************************************************************************

  inline void CompoundType :: setSize(const size_t &buf)
  {
     numBytes = buf;
  }
  //**************************************************************************

  inline int CompoundType :: getSize() const
  {
     return numBytes;
  }
  //**************************************************************************
  inline int CompoundType :: getNumFields()
  {
      return newtype.size();
  }
  //**************************************************************************
  inline std::string  CompoundType ::  name( int i ) const
  {
      return newtype[i].name;

  }
   //**************************************************************************
   inline size_t CompoundType ::  offset(int i ) const
   {
      return newtype[i].offset;
   }
   //**************************************************************************
   inline AtomicType  CompoundType :: type( int i ) const
   {
      return newtype[i].atom;
   }
  //**************************************************************************
  inline hid_t AtomicType :: getSize() const
  {
     std::cout << " NOT IMPLEMENTED " << std::endl;
     exit(0);

     return 0;
  }
  //**************************************************************************

  inline hid_t AtomicType :: get_hdf5_type() const
  {
     switch( atom )
     {
      case FLOAT:
           return (H5Tcopy(H5T_NATIVE_FLOAT));
      case DOUBLE:
           return (H5Tcopy(H5T_NATIVE_DOUBLE));
      case INT:
           return (H5Tcopy(H5T_NATIVE_INT));
      case CHAR:
           return (H5Tcopy(H5T_NATIVE_CHAR));
     }
     return 0;
  }
  //**************************************************************************

  inline hid_t CompoundType :: get_hdf5_type() const
  {
    hid_t vDatatype  = H5Tcreate( H5T_COMPOUND, numBytes);

    int rank = 1;
    hid_t hdf5T, array_type;
    hsize_t array_dims[10];

    for( int i = 0; i < newtype.size(); i++) {
      switch( newtype[i].type )
        {
        case ATOMIC:
          hdf5T = newtype[i].atom.get_hdf5_type();
          H5Tinsert(vDatatype, newtype[i].name.c_str(), newtype[i].offset, hdf5T);
          break;
        case ARRAY:
          if( newtype[i].array.type == ATOMIC) {
            hdf5T          = newtype[i].array.atom.get_hdf5_type();
            rank           = newtype[i].array.rank;
            for(int k = 0; k < rank; k++)
              array_dims[k]  = newtype[i].array.dimension[k];
#ifdef HDF5V1_2
            H5Tinsert_array(vDatatype, newtype[i].name.c_str(), newtype[i].offset, 
                            rank, array_dims, NULL, hdf5T);
#else
            array_type = H5Tarray_create( hdf5T,rank, array_dims,NULL); 
            H5Tinsert(vDatatype, newtype[i].name.c_str(), newtype[i].offset, 
	                   array_type);
            H5Tclose(array_type);
#endif

          } else  {
            hid_t vDatatype2 = H5Tcreate( H5T_COMPOUND, newtype[i].compound->getSize());
            hdf5T = newtype[i].compound->get_hdf5_type();
            rank  = newtype[i].array.rank;
            for(int k = 0; k < rank; k++)
              array_dims[k]  = newtype[i].array.dimension[k];
#ifdef HDF5V1_2
            H5Tinsert_array(vDatatype, newtype[i].name.c_str(), newtype[i].offset, 
                            rank, array_dims, NULL, hdf5T);
#else
            array_type = H5Tarray_create(hdf5T, rank, array_dims,NULL); 
            H5Tinsert(vDatatype, newtype[i].name.c_str(), newtype[i].offset, 
	                   array_type);
            H5Tclose(array_type);
#endif

          }
          break;
        }
    }
    return vDatatype;
  }
  //*****************************************************************************
  inline hid_t get_hdf5_type(float f)
  { return H5Tcopy(H5T_NATIVE_FLOAT);}

  inline hid_t get_hdf5_type(double d)
  { return H5Tcopy(H5T_NATIVE_DOUBLE);}

  inline hid_t get_hdf5_type(char c)
  { return H5Tcopy(H5T_NATIVE_CHAR);}

  inline hid_t get_hdf5_type( int i)
  { return H5Tcopy(H5T_NATIVE_INT);}

  inline hid_t get_hdf5_type( bool b)
  { return H5Tcopy(H5T_NATIVE_HBOOL);}

  //*****************************************************************************
}

#endif
