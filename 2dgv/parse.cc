#include "parse.h"
#include <cstdio>
#include "debug.h"

namespace parse {
  
    using namespace std ;
    void kill_white_space(istream &s) {

      while(!s.eof() && isspace(s.peek()))
        s.get() ;

      return ;
    }
    
    bool is_name(istream &s) {
        kill_white_space(s) ;
        int ch = s.peek() ;
        return isalpha(ch) || ch == '_' ;
    }
    
    string get_name(istream &s) {
        if(!is_name(s))
          return "" ;
        string str ;
        while(!s.eof() && (s.peek() != EOF) &&
              (isalnum(s.peek()) || (s.peek() == '_')) )
          str += s.get() ;

        return str ;
    }

    bool is_int(istream &s) {
        kill_white_space(s) ;
        return isdigit(s.peek()) || s.peek()=='-' || s.peek()=='+' ;
    }
    
    long get_int(istream &s) {
        if(!is_int(s))
          return 0 ;
        long l = 0 ;
        s >> l ;
        return l ;
    }

    bool is_real(istream &s) {
        kill_white_space(s) ;
        const char ch = s.peek() ;
        return isdigit(ch) || ch=='-' || ch=='+' || ch =='.' ;
    }
    
    double get_real(istream &s) {
        if(!is_real(s))
          return 0.0 ;
        double r = 0.0 ;
        s >> r ;
        return r ;
    }

    bool is_string(istream &s) {
        kill_white_space(s) ;
        return s.peek() == '\"' ;
    }
    
    string get_string(istream &s) {
        if(!is_string(s))
          return "" ;
        string str ;
#ifdef DEBUG        
        if(s.eof())
          cerr << "s.eof() true in parse::get_string" << endl  ;
#endif
        s.get() ;
        int ch = s.get() ;
        while(ch != '\"' &&!s.eof()) {
            str += ch ;
            ch = s.get() ;
        }
#ifdef DEBUG
        if(ch!='"')
          cerr << "no closing \" in parse::get_string" << endl ;
#endif
        return str ;
    }

    bool  is_token(istream &s, const string &token) {
        kill_white_space(s) ;
        for(int i=0;i<int(token.size());++i) {
            if(s.peek() != token[i]) {
                for(--i;i>=0;--i)
                  s.putback(token[i]) ;
                return false ;
            }
            s.get() ;
        }
        for(int i=token.size()-1;i>=0;--i) 
          s.putback(token[i]) ;
        return true ;
    }
    
    bool get_token(istream &s, const string &token) {
        kill_white_space(s) ;
        for(int i=0;i<int(token.size());++i) {
            if(s.peek() != token[i]) {
                for(--i;i>=0;--i)
                  s.putback(token[i]) ;
                return false ;
            }
            s.get() ;
        }
        return true ;
    }
}    

