#include <Loci.h>
#include <iostream.h>

using namespace std ;


int main(int ac, char *av[])
{
/*
  entitySet A,B,C,D,E ;

  A = interval(1,10) ;
  B = interval(7,20) ;

  E = A - B ;  // Replace A with itself but remove elements in C

  cout << "E = " << E << endl ;
*/
  sequence    sset;
  entitySet   eset;
  store<int>  astore;

  eset  +=  7; sset  +=  7;
  eset  +=  6; sset  +=  6;
  eset  +=  5; sset  +=  5;
  eset  +=  4; sset  +=  4;
  eset  +=  3; sset  +=  3;
  eset  +=  2; sset  +=  2;
  eset  +=  1; sset  +=  1;
  eset  +=  0; sset  +=  0;

  astore.allocate( eset );

  astore[0] = 0;
  astore[1] = 1;
  astore[2] = 2;
  astore[3] = 3;
  astore[4] = 4;
  astore[5] = 5;
  astore[6] = 6;
  astore[7] = 7;

  entitySet :: const_iterator  ei;

  for( ei = eset.begin(); ei != eset.end(); ++ei)
       cout << *ei;
  cout << endl;

  sequence  :: const_iterator  si;
  for( si = sset.begin(); si != sset.end(); ++si)
       cout << *si;
  cout << endl;
}
