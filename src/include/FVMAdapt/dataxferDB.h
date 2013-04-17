#ifndef DATAXFERDB_H
#define DATAXFERDB_H

#include <Loci>
#include <map>
#include <vector>
#include <string>

namespace Loci {

  class DataXFERDatabase {
    std::map<std::string,Loci::storeRepP> dbmap ;
  public:
    Loci::storeRepP getItem(const string &s) {
      std::map<std::string,Loci::storeRepP>::iterator mi ;
      mi = dbmap.find(s) ;
      if(mi == dbmap.end()) 
	return 0 ;
      return mi->second ;
    }
    void insertItem(const string &s, Loci::storeRepP item) {
      dbmap[s] = item ;
    }
    void deleteItem(const string &s) {
      std::map<std::string,Loci::storeRepP>::iterator mi ;
      mi = dbmap.find(s) ;
      if(mi != dbmap.end()) 
	dbmap.erase(mi) ;
    }
    void clearDB() {
      dbmap = std::map<std::string,Loci::storeRepP>() ;
    }
    std::vector<std::string> nameList() const {
      std::vector<std::string> nlist ;
      std::map<std::string,Loci::storeRepP>::const_iterator mi ;
      for(mi=dbmap.begin();mi!=dbmap.end();++mi) {
	nlist.push_back(mi->first) ;
      }
      return nlist ;
    }
  } ;

  extern DataXFERDatabase DataXFER_DB ;
}
#endif
