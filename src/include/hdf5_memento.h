#ifndef HDF5_MEMENTO_H_
#define HDF5_MEMENTO_H_

#include <typeinfo>
#include <data_traits.h>

#include <vector>
#include <algorithm>

namespace Loci {

template<class T>
class Memento
{ };


template< class T>
class Memento< std::vector<T> >
{
   typedef data_schema_converter_traits<std::vector<T> > converter_traits; 

public:
  explicit Memento( const std::vector<T> &buf)
  { obj = buf; };

  int  getSize()  const
  { return obj.size(); };

  void getState(typename converter_traits::memento_type *buf, int &size)
  {
    size = obj.size();
    for( int i = 0; i < size; i++) 
         buf[i] = obj[i];
  }
  std::vector<T> setState(typename converter_traits::memento_type *buf, int size)
  {
    obj.clear();
 
    for( int i = 0; i < size; i++)
         obj.push_back( buf[i] );

    return( obj );
  }

private:
   std::vector<T>  obj;
};

}
#endif
