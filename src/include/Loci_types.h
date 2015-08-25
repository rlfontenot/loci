//#############################################################################
//#
//# Copyright 2008, 2015, Mississippi State University
//#
//# This file is part of the Loci Framework.
//#
//# The Loci Framework is free software: you can redistribute it and/or modify
//# it under the terms of the Lesser GNU General Public License as published by
//# the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The Loci Framework is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# Lesser GNU General Public License for more details.
//#
//# You should have received a copy of the Lesser GNU General Public License
//# along with the Loci Framework.  If not, see <http://www.gnu.org/licenses>
//#
//#############################################################################
#ifndef LOCI_TYPES_H
#define LOCI_TYPES_H

#ifdef HAVE_CONFIG_H
#include <config.h> // This must be the first file included
#endif
#include <Config/conf.h>


#include <data_traits.h>

#include <Tools/basic_types.h>
#include <Tools/options_list.h>

namespace Loci {
  template <class T>
    struct data_schema_traits< vector3d<T> > {
      typedef IDENTITY_CONVERTER Schema_Converter;
      static DatatypeP get_type() {
        CompoundDatatypeP ct = CompoundFactory(vector3d<T>() ) ;

        LOCI_INSERT_TYPE(ct,vector3d<T>,x) ;
        LOCI_INSERT_TYPE(ct,vector3d<T>,y) ;
        LOCI_INSERT_TYPE(ct,vector3d<T>,z) ;
        return DatatypeP(ct) ;
      }
    };


  template <class T>
    struct data_schema_traits< tensor3d<T> > {
      typedef IDENTITY_CONVERTER Schema_Converter;
      static DatatypeP get_type() {
        CompoundDatatypeP ct = CompoundFactory(tensor3d<T>()) ;

        LOCI_INSERT_TYPE(ct,tensor3d<T>,x) ;
        LOCI_INSERT_TYPE(ct,tensor3d<T>,y) ;
        LOCI_INSERT_TYPE(ct,tensor3d<T>,z) ;
        return DatatypeP(ct) ;
      }
    };
  
    
  
  template <class T>
    struct  data_schema_traits< vector2d<T> > {
    public:
      typedef IDENTITY_CONVERTER Schema_Converter;
      static DatatypeP get_type() {
        vector2d<T> t ;
        CompoundDatatypeP ct = CompoundFactory(t) ;
        LOCI_INSERT_TYPE(ct,vector2d<T>,x) ;
        LOCI_INSERT_TYPE(ct,vector2d<T>,y) ;
        return DatatypeP(ct) ;
      }
    };
  
  
  

  template <class T,size_t n> inline std::ostream &
    operator<<(std::ostream &s, const Array<T,n> &v) {
    for(size_t i=0;i<n;++i)
      s << v[i] << ' ' ;
    return s ;
  }

  template <class T,size_t n> inline std::istream &
    operator>>(std::istream &s, Array<T,n> &v) {
    for(size_t i=0;i<n;++i)
      s >> v[i] ;
    return s ;
  }

  
  template <class T,size_t n> 
  class data_schema_traits< Array<T,n> > {
  public:
    typedef IDENTITY_CONVERTER Schema_Converter;
    static DatatypeP get_type() {
      int dim = n ;
      return new ArrayType(getLociType(T()),sizeof(Array<T,n>),1,&dim) ;
    }
  };

  template <> struct data_schema_traits<options_list> {
    typedef USER_DEFINED_CONVERTER Schema_Converter ;
    typedef char Converter_Base_Type ;
    typedef StringStreamConverter<options_list> Converter_Type ;
  } ;

      
}

#endif
