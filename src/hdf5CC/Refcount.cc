#include "Refcount.h"

namespace H5 {

// Creates a reference counter to be used by an HDF5 object
ReferenceCounter::ReferenceCounter() : counter(1) {}; 
 
// Returns the current value of the reference counter
int ReferenceCounter::getCounter () const { return counter; }

// Increments the reference counter as a copy of the object that uses
// this counter is created.
void ReferenceCounter::increment() { counter++; }

// Decrements the reference counter as a copy of the object that uses
// this counter is destroyed.
void ReferenceCounter::decrement() { counter--; }

// Decrements the reference counter then determines if there are no more 
// reference to the object that uses this counter
bool ReferenceCounter::noReference()
{
   counter--;
   return( counter == 0 ? true:false );
}
} // end namespace
