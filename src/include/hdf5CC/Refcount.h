#ifndef _H5ReferenceCounter_H
#define _H5ReferenceCounter_H
;
namespace H5 {

class ReferenceCounter {
   private:
	int counter; // keeps track of number of copies of an object

   public:
	// Creates a reference counter to be used by an HDF5 object
	ReferenceCounter();

        int getCounter () const;
        void increment();
        void decrement();

	// this bool function is used to determine whether to close an
	// HDF5 object when there are no more reference to that object
	bool noReference();

	~ReferenceCounter() {};
};

}
#endif
