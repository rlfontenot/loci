// Example program illustrating how to manipluate Entity, entitySet,
// and sequence data structures in Loci

// Every Loci program includes the Loci header file
#include <Loci.h>

// includes for standard C++ functions
#include <iostream>
using std::cout ;
using std::endl ;
#include <vector>
using std::vector ;

int main()
{
  // The Entity 
  Entity e10(10) ;  // First, the entity labeled by integer 10

  // The interval (A collection of consecutively labeled entities)
  // For example, onedigit is an interval of entity labels consisting of
  // only one decimal digit.
  interval onedigit = interval(0,9) ;   

  // entitySet provides abilities to hold general sets of entities
  // They can be initialized to other intervals
  entitySet A = onedigit ;
  entitySet B = interval(14,100) ;
  entitySet C = interval(5,15) ;

  // We can also add an individual entity to any existing entitySet using the
  // += operator, for example we can include the entitiy e10 in set B:
  B += e10 ;

  //  A = ([0,9])
  //  B = ([10,10][14,100])
  //  C = ([5,15])
  cout << "A = " << A << endl ;
  cout << "B = " << B << endl ;
  cout << "C = " << C << endl ;

  // Or they can be formed by performing set operations on other entitySet's
  // For example, entitySet D becomse the union of A and B
  // That is D = ([0-9],[14-100])
  entitySet D = A + B ;

  //  D = ([0,10][14,100])
  cout << "D = " << D << endl ;

  // Note that entitySet's remain ordered and do not contain duplicates, that
  // is they are true sets.   For example
  entitySet E = B + C ;
  //  E = B union C =  ([5,100])
  // [gives the set ([5-100]) without duplicating 14 and 15]
  cout << "E = B union C =  " << E << endl ;

  // Other operations are also supported, for example
  // A & C  gives the intersection of A and C (the interval [5-9])
  cout << "A intersect C = " << (A & C) << endl ;
  // A - C gives the set difference (A take away C) (the interval [0-4])
  cout << "A take away C = " << (A - C) << endl ;

  // We can also create entitySet's from arbitrary lists of entity identifiers
  // For example,
  int values[] = {10,15,12,1,14,16,17} ;
  entitySet vset = create_entitySet(values,values+7) ;

  // vset = ([1,1][10,10][12,12][14,17])
  cout << "vset = " << vset << endl ;

  // Note, create_entitySet also works with std::vector and begin() and end()
  // See a standard C++ book for more details on using vector<> and other
  // STL containers.
  vector<int> vec ;
  for(int i=10;i>0;--i)
    vec.push_back(i) ;

  entitySet vecset = create_entitySet(vec.begin(),vec.end()) ;

  // vecset = ([1,10])
  cout << "vecset = " <<vecset << endl ;

  // We can also iterate (loop) over an entitySet in a fashion similar to
  // how we iterate over standard C++ containers.  For example, to iterate
  // over all of the entities in vset we would write a loop as follows

  // First we create iterator ei for entity sets.
  entitySet::const_iterator ei ;

  // Then we use the iterator to loop over a given entitySet using the
  // begin() and end() methods.  For example to loop over the entities
  // contained in vset we write:
  cout << "looping over vset:" ;
  for(ei = vset.begin(); ei != vset.end(); ++ei) 
    cout << " " << *ei ;
  cout << endl ;
  // above outputs:
  // looping over vset: 1 10 12 14 15 16 17

  // Other useful methods include Min() and Max() which can be used to find
  // the largest and smallest integer labels of entities contained in a given
  // entitySet.  For example, vset.Max() == 17 and vset.Min() == 1

  // vset.Min() == 1, vset.Max() == 17
  cout << "vset.Min() = " << vset.Min() << endl ;
  cout << "vset.Max() = " << vset.Max() << endl ;

  // Similarly, we can check the number of entities contained within an
  // entitySet by using the size() method.  For example

  // vset.size() == 7
  cout << "vset contains " << vset.size() << " entities" << endl ;

  // We can also check to see if a particular entity label is in a
  // entitySet using the inSet() method.  For example

  if(A.inSet(5))
    cout << "entity labeled 5 is in entitySet A" << endl ;

  // Sequences provide a way of storing ordered lists of entities.  Usually,
  // users don't need to worry about creating sequences directly in Loci,
  // but instead Loci creates sequences to describe the order in which
  // calculations will be carried out.  However, for completeness we
  // will give some examples here on how to create and iterate of sequences.

  // We can create an arbitrary sequence from a list of integers
  // in a fashion similar to the create_entitySet function.  For
  // example:
  int listvals[] = {10,15,12,1,14,16,17} ;
  
  sequence vseq = create_sequence(listvals,listvals+7) ;

  // vseq = ([10,10][15,15][12,12][1,1][14,14][16,17])
  cout << "vseq = " << vseq << endl ;

  // Note that we can also create a sequence from an entitySet, however
  // in this case the sequence will contain the contents of the entitySet
  // in increasing order.  For example
  sequence Aseq = A ;

  // Aseq = ([0-9])
  cout << "Aseq = " << Aseq << endl ;

  // Also we can append to sequences, for example:
  sequence Cseq = Aseq + vseq;

  // Cseq = ([0,10][15,15][12,12][1,1][14,14][16,17])
  cout << "Cseq = " << Cseq << endl ;

  // Similarly we can reverse the order of a sequence using the Reverse()
  // method.  For example
  Cseq.Reverse() ;

  // Cseq = ([17,16][14,14][1,1][12,12][15,15][10,0])
  cout << "reversed Cseq = " << Cseq << endl ;

  return 0 ;
}
