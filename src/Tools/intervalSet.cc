#include <Tools/debug.h>
#include <Tools/tools.h>
#include <Tools/intervalSet.h>
#include <Tools/stream.h>

#include <algorithm>
#include <functional>



//using namespace std::rel_ops ;
using std::pair ;

using std::binary_search ;
using std::equal_range ;
using std::reverse ;
using std::swap ;


namespace Loci {
    intervalSet::rep_holder *intervalSet::rhp = 0 ;
    sequence::rep_holder *sequence::rhp = 0 ;
    const intervalSet EMPTY ;

    /* The intervalSetRep class does all of the work for the intervalSet class.
     * It is referenced through a counted pointer (Handle) template.
     */
    ostream &operator<<(ostream &s, const interval &i) {
	s << '[' ;
	if(i.first==UNIVERSE_MIN)
	    s << '#' ;
	else 
	    s << i.first ;
	s << ',' ;
	if(i.second==UNIVERSE_MAX)
	    s << '#' ;
	else
	    s << i.second ;
	s << ']' ;
	return s;
    }


  istream &operator>>(istream &s, interval &i) {
    char ch ;
    do{
      ch = s.get() ;
    } while(ch==' ' || ch=='\n') ;
    if(ch!='[') {
      cerr << "Incorrect format when reading interval" << endl ;
      cerr << "expected a '[' but got a '" << ch << "'" << endl ;
      s.putback(ch) ;
      return s ;
    }
    if(s.peek() == '#') {
      ch = s.get() ;
      i.first = UNIVERSE_MIN ;
    } else
      s >> i.first ;
    do{
      ch = s.get() ;
    } while(ch==' ' || ch=='\n') ;
    if(ch!=',') {
      cerr << "Incorrect format when reading interval" << endl ;
      cerr << "expected a ',' but got a '" << ch << "'" << endl ;
      s.putback(ch) ;
      return s ;
    }
    if(s.peek() == '#') {
      ch = s.get() ;
      i.second = UNIVERSE_MAX ;
    } else
      s >> i.second ;
    
    do{
      ch = s.get() ;
    } while(ch==' ' || ch=='\n') ;
    if(ch!=']') {
      cerr << "Incorrect format when reading interval" << endl ;
      cerr << "expected a ']' but got a '" << ch << "'" << endl ;
      s.putback(ch) ;
      return s ;
    }
    return s ;
  }
    
  /*

    IntervalSet Class Code starts Here
   
    */

  ostream &intervalSet::Print(ostream &s) const {
    s << "(" ;
    for(int i=0;i<num_intervals();++i)
      s << (*Rep)[i] ;
    s << ")" ;
    return s ;
  }

  istream &intervalSet::Input(istream &s) {
    Rep.New() ;

    char ch ;
    // Startup ;

    // Kill white space
    do{
      ch = s.get() ;
    } while(isspace(ch)) ;

    // Check for correct startup character
    if(ch != '(') {
      cerr << "Incorrect format when reading intervalSet" << endl ;
      s.putback(ch) ;
      return s ;
    }

    // get intervals
    for(;;) {
      // Skip whitespace
      while(isspace(s.peek()))
        ch = s.get() ;
      // If next input is not an interval, then we are finished.
      if(s.peek()!='[')
        break ;
      // grab interval from input
      interval ivl  ;
      s >> ivl ;
      // Add it up
      Union(ivl) ;
    }
    // Check that everything finished ok, then finish
    ch = s.get() ;
    if(ch!=')') {
      cerr << "Incorrect format when reading intervalSet" << endl ;
      s.putback(ch) ;
    }
      
    return s ;
  }

  // This partial ordering function defines an implicit equivalence
  // relation between intervals if they intersect one another
  namespace {
    inline bool interval_porder_intersect(const interval &i1, 
					  const interval &i2) {
      return (i1.first < i2.first) && (i1.second < i2.first)  ;
    }
  }

  bool intervalSet::inSet(int_type indx) const {
    return binary_search(Rep->begin(),Rep->end(),interval(indx,indx),
                         interval_porder_intersect) ;
  }

  // This partial ordering function defines an implicit equivalence
  // relation between intervals that are ordered lowvalue then high value
  // so that if they intersect or they touch they are in the same equivalence
  // relation.
  namespace {
    inline bool interval_porder_union(const interval &i1, 
				      const interval &i2) {
      return i1.second+1 < i2.first ;
    }
  }

    void Union(Handle<pair_vector> &Rep, const interval &ivl) {
    Rep.MakeUnique() ;
    interval i(min(ivl.first,ivl.second),
               max(ivl.first,ivl.second)) ;
    if(num_intervals(Rep) == 0) {
      Rep->push_back(i) ;
      return ;
    }
    typedef pair_vector intervalSetRep;
    pair<intervalSetRep::iterator,intervalSetRep::iterator> range ;
    range = equal_range(Rep->begin(),Rep->end(),ivl,
                        interval_porder_union) ;
    int range_size = range.second - range.first ;
    fatal(range_size<0) ;
    switch(range_size) {
    case 0:
      Rep->insert(range.first,i) ;
      break ;
    case 1:
      (*range.first).first = min(i.first,(*range.first).first) ;
      (*range.first).second = max(i.second,(*range.first).second) ;
      break ;
    default:
      int_type lowval = (*range.first).first ;
      int_type hival  = (*(range.second-1)).second ;
      (*range.first).first = min(lowval,i.first) ;
      (*range.first).second = max(hival,i.second) ;
      Rep->erase(range.first+1,range.second) ;
      break ;
    }
  }





    Handle<pair_vector> Union(const Handle<pair_vector> &Rep1,
			      const Handle<pair_vector> &Rep2) {
	pair_vector::const_iterator i1 = Rep1->begin() ;
	pair_vector::const_iterator i2 = Rep2->begin() ;
	const pair_vector::const_iterator e1 = Rep1->end() ;
	const pair_vector::const_iterator e2 = Rep2->end() ;
	const unsigned int size = Rep1->size() + Rep2->size() ;
	Handle<pair_vector> Rep ;
	Rep->reserve(size) ;
	while(i1!=e1 && i2!=e2) {
	    if(interval_porder_union(*i1,*i2)) {
		Rep->push_back(*i1) ;
		++i1 ;
	    } else if(interval_porder_union(*i2,*i1)) {
		Rep->push_back(*i2) ;
		++i2 ;
	    } else {
		interval ivl(min((*i1).first,(*i2).first),
			     max((*i1).second,(*i2).second)) ;
		++i1 ;
		++i2 ;
		bool flag = true ;
		while(flag) {
		    flag = false ;
		    if(i1!=e1 && !interval_porder_union(ivl,*i1)
		       && !interval_porder_union(*i1,ivl)) {
			ivl.second = max(ivl.second,(*i1).second) ;
			++i1 ;
			flag = true ;
		    }
		    if(i2!=e2 && !interval_porder_union(ivl,*i2)
		       && !interval_porder_union(*i2,ivl)) {
			ivl.second = max(ivl.second,(*i2).second) ;
			++i2 ;
			flag = true ;
		    }
		} 
		Rep->push_back(ivl) ;
	    }
	}
	while(i1!=e1) {
	    Rep->push_back(*i1) ;
	    ++i1 ;
	}
	while(i2!=e2) {
	    Rep->push_back(*i2) ;
	    ++i2 ;
	}
	return Rep ;
    }


  void intervalSet::Complement() {
    Rep.MakeUnique() ;
    if(Rep->size() == 0) {
      Rep->push_back(interval(UNIVERSE_MIN,UNIVERSE_MAX)) ;
      return ;
    }

    if((*Rep)[0].first == UNIVERSE_MIN) {
      intervalSetRep::iterator stop = Rep->end()-1 ;
      for(intervalSetRep::iterator i=Rep->begin();i!=stop;) {
        i->first = i->second+1 ;
        intervalSetRep::iterator next = i+1 ;
        i->second = next->first-1 ;
        i = next ;
      }
      if((*stop).second == UNIVERSE_MAX)
        Rep->erase(stop) ;
      else {
        stop->first = stop->second+1 ;
        stop->second = UNIVERSE_MAX ;
      }
      return ;
    }

    intervalSetRep::iterator i = Rep->begin() ;
    int_type sec = i->second ;
    int_type fir = i->first ;
    (*i).first = UNIVERSE_MIN ;
    (*i).second = fir-1 ;
    for(++i;i!=Rep->end();++i) {
      fir = i->first ;
      i->first = sec+1 ;
      sec = i->second ;
      i->second = fir-1 ;
    }
    if(sec != UNIVERSE_MAX)
      Rep->push_back(interval(sec+1,UNIVERSE_MAX)) ;
    return ;
            

  }

  intervalSet intervalSet::Complement(const intervalSet &set) {
    const int size = set.Rep->size() ;
    Handle<intervalSetRep> Rep ;
    Rep->reserve(size+1) ;
    if(size == 0)
      Rep->push_back(interval(UNIVERSE_MIN,UNIVERSE_MAX)) ;
    else {
      if(set.Rep->front().first != UNIVERSE_MIN)
        Rep->push_back(interval(UNIVERSE_MIN,set.Rep->front().first-1)) ;

      int sec = set.Rep->front().second ;
      for(int i=1;i<size;++i) {
        const interval &ci = (*(set.Rep))[i] ;
        Rep->push_back(interval(sec+1,ci.first-1)) ;
        sec = ci.second ;
      }
      if(sec!=UNIVERSE_MAX) 
        Rep->push_back(interval(sec+1,UNIVERSE_MAX)) ;
    }
        
    intervalSet result ;
    if(Rep->size() != 0)
      result.Rep = Rep ;
    return result ;
  }

  intervalSet &intervalSet::operator>>=(int rotval) {
    Rep.MakeUnique() ;
    for(intervalSetRep::iterator i = Rep->begin();i!=Rep->end();++i) {
      if(i->first != UNIVERSE_MIN)
        i->first += rotval ;
      if(i->second != UNIVERSE_MAX)
        i->second += rotval ;
    }
    return *this ;
  }

  intervalSet &intervalSet::operator<<=(int rotval) {
    Rep.MakeUnique() ;
    for(intervalSetRep::iterator i = Rep->begin();i!=Rep->end();++i) {
      if(i->first != UNIVERSE_MIN)
        i->first -= rotval ;
      if(i->second != UNIVERSE_MAX)
        i->second -= rotval ;
    }
    return *this ;
  }


  //  Finding an intersection to an interval is a simple matter of extracting
  //  the equivalence set for under the partial ordering for overlaping
  //  intervals and then truncating the range to the interval bounds.
  void intervalSet::Intersection(const interval &ivl) {
    pair<intervalSetRep::iterator,intervalSetRep::iterator> range ;
    range = equal_range(Rep->begin(),Rep->end(),ivl,
                        interval_porder_intersect) ;
    int range_size = range.second-range.first ;
    fatal(range_size<0) ;
    if(range_size == 0) {
      Rep = EMPTY.Rep ;
      return ;
    }
    Rep.MakeUnique() ;
    Rep->erase(range.second,Rep->end()) ;
    Rep->erase(Rep->begin(),range.first) ;
    range.first = Rep->begin() ;
    range.second = Rep->end()-1 ;
    (*range.first).first = max(ivl.first,(*range.first).first) ;
    (*range.second).second = min(ivl.second,(*range.second).second) ;
    return ;
  }


  void intervalSet::Intersection(const intervalSet &ptn) {
    intervalSet tmp = intervalSet::Intersection(*this,ptn) ;
    Rep = tmp.Rep ; 
  }

  intervalSet intervalSet::Intersection(const intervalSet &set1,
                                        const intervalSet &set2) {
    if(set1 == EMPTY || set2 == EMPTY)
      return EMPTY ;
      
    intervalSetRep::const_iterator i1 = set1.Rep->begin() ;
    intervalSetRep::const_iterator i2 = set2.Rep->begin() ;
    const intervalSetRep::const_iterator e1 = set1.Rep->end() ;
    const intervalSetRep::const_iterator e2 = set2.Rep->end() ;
    const unsigned int size = set1.Rep->size() + set2.Rep->size() ;
    Handle<intervalSetRep> Rep ;
    Rep->reserve(size) ;
    while(i1!=e1 && i2!=e2) {
      if(interval_porder_intersect(*i1,*i2))
        ++i1 ;
      else if(interval_porder_intersect(*i2,*i1))
        ++i2 ;
      else {
        interval ivl(max((*i1).first,(*i2).first),
                     min((*i1).second,(*i2).second)) ;
        Rep->push_back(ivl) ;
        if((*i1).second > (*i2).second)
          ++i2 ;
        else
          ++i1 ;
      }
    }
    intervalSet result ;
    if(Rep->size() != 0)
      result.Rep = Rep ;
    return result ;
  }

  ostream &sequence::Print(ostream &s) const {
    s << "(" ;
    for(int i=0;i<num_intervals();++i)
      s << (*Rep)[i] ;
    s << ")" ;
    return s ;
  }

  istream &sequence::Input(istream &s) {
    Rep.New() ;
     
    char ch ;
    // Startup ;

    // Kill white space
    do{
      ch = s.get() ;
    } while(isspace(ch)) ;

    // Check for correct startup character
    if(ch != '(') {
      cerr << "Incorrect format when reading intervalSet" << endl ;
      s.putback(ch) ;
      return s ;
    }

    // get intervals
    for(;;) {
      // Skip whitespace
      while(isspace(s.peek()))
        ch = s.get() ;
      // If next input is not an interval, then were finished
      if(s.peek()!='[')
        break ;
      // grab interval from input
      interval ivl  ;
      s >> ivl ;
      // Add it up
      Append(ivl) ;
    }
    // Check that everything finished ok, then finish
    ch = s.get() ;
    if(ch!=')') {
      cerr << "Incorrect format when reading intervalSet" << endl ;
      s.putback(ch) ;
    }
      
    return s ;
  }

  void sequence::Append(const interval &ivl) {
    Rep.MakeUnique() ;
    if(Rep->size() == 0) {
      Rep->push_back(ivl) ;
      return ;
    }
            
    interval &back = Rep->back() ;
    int dir = back.second-back.first ;
    int dirivl = ivl.second-ivl.first ;
    if((ivl.first == (back.second+1)) && dirivl >= 0 && dir >= 0) {
      (Rep->back()).second = ivl.second ;
      return ;
    }
    if((ivl.first == (back.second-1)) && dirivl <= 0 && dir <= 0) {
      (Rep->back()).second = ivl.second ;
      return ;
    }            
    Rep->push_back(ivl) ;
  }

  void sequence::Append(const sequence &seq) {
    if(seq.Rep->size() == 0)
      return ;
    Rep.MakeUnique() ;
    sequenceRep::const_iterator i = seq.Rep->begin() ;
    Append(*i) ;
    ++i ;
    Rep->insert(Rep->end(),i,seq.Rep->end()) ;
  }

  void sequence::Append(const intervalSet &ptn) {
    sequence seq(ptn) ;
    Append(seq) ;
  }

  sequence & sequence::Reverse() {
    Rep.MakeUnique() ;
    reverse(Rep->begin(),Rep->end()) ;
    for(sequenceRep::iterator i=Rep->begin();i!=Rep->end();++i)
      swap(i->first,i->second) ;
    return *this ;
  }

}
