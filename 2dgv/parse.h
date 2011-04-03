#ifndef PARSE_H
#define PARSE_H

#include <iostream>
#include <string>

namespace parse {
    
    void kill_white_space(std::istream &s) ;
    
    bool is_name(std::istream &s) ;
    std::string get_name(std::istream &s) ;

    bool is_int(std::istream &s) ;
    long get_int(std::istream &s) ;

    bool is_real(std::istream &s) ;
    double get_real(std::istream &s) ;

    bool is_string(std::istream &s) ;
    std::string get_string(std::istream &s) ;

    bool is_token(std::istream &s, const std::string &token) ;
    bool get_token(std::istream &s, const std::string &token) ;
}

#endif
