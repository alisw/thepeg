// -*- C++ -*-
//
// StringUtils.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the StringUtils class.
//

#include "StringUtils.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;

StringUtils::StringVector StringUtils::split(string s, string ws) {
  StringVector r;
  string::size_type i = s.find_first_not_of(ws);
  while ( i != string::npos ) {
    s = s.substr(i);
    i = s.find_first_of(ws);
    r.push_back(s.substr(0, i));
    if ( i == string::npos ) return r;
    s = s.substr(i);
    i = s.find_first_not_of(ws);
  }
  return r;
}

string StringUtils::car(string s, string ws) {
  string::size_type p = s.find_first_not_of(ws);
  if ( p == string::npos ) return "";
  s = s.substr(p);
  return s.substr(0, s.find_first_of(ws));
}

string StringUtils::cdr(string s, string ws) {
  string::size_type p = s.find_first_not_of(ws);
  if ( p == string::npos ) return "";
  s = s.substr(p);
  p = s.find_first_of(ws);
  if ( p == string::npos ) return "";
  s = s.substr(p);
  p = s.find_first_not_of(ws);
  if ( p == string::npos ) return "";
  return s.substr(p);
}

string StringUtils::stripws(string str) {
  string::size_type i = str.find_first_not_of(" \t\r\n");
  if ( i != string::npos ) str = str.substr(i);
  i = str.find_last_not_of(" \t\r\n");
  return str.substr(0, i + 1);
}

map<string,string> StringUtils::
xmlAttributes(string tag, string line, string::size_type curr) {

  map<string,string> ret;

  // First find the beginning and end of the tag
  pos_t begin = line.find("<" + tag, curr);
  if ( begin == end ) return ret;
  curr = begin + tag.length() + 2;

  while ( true ) {

    pos_t close = line.find(">", curr);

    // Now skip some white space to see if we can find an attribute.
    curr = line.find_first_not_of(" \t\r\n", curr);
    if ( curr == end || curr >= close ) return ret;

    pos_t tend = line.find_first_of("= \t\r\n", curr);
    if ( tend == end || tend >= close ) return ret;

    string name = line.substr(curr, tend - curr);
    curr = line.find("=", curr) + 1;

    // OK now find the beginning and end of the atribute.
    curr = line.find("\"", curr);
    if ( curr == end || curr >= close ) return ret;
    pos_t bega = ++curr;
    curr = line.find("\"", curr);
    while ( curr != end && line[curr - 1] == '\\' )
      curr = line.find("\"", curr + 1);

    string value = line.substr(bega, curr == end? end: curr - bega);

    ret[name] = value;

    ++curr;

  }
}


string StringUtils::dirname(string file) {
  string::size_type pos = file.rfind('/');
  if ( pos == string::npos ) return "";
  return file.substr(0,pos);
}

string StringUtils::basename(string file) {
  string::size_type pos = file.rfind('/');
  if ( pos == string::npos ) return file;
  return file.substr(pos + 1);
}

string StringUtils::remsuf(string file) {
  string::size_type pos = file.rfind('.');
  if ( pos == string::npos ) return "";
  return file.substr(0,pos);
}

string StringUtils::suffix(string file) {
  string::size_type pos = file.rfind('.');
  if ( pos == string::npos ) return file;
  return file.substr(pos + 1);
}

string StringUtils::typeName(const type_info & t) {
  string ret;
  string name = t.name();
  if ( name[0] != 'N' ) return name;

  string::size_type pos = 1;
  while ( pos < name.length() && name[pos] != 'E' ) {
    int size = 0;
    while ( pos < name.length() && isdigit(name[pos]) )
      size = size*10 + int(name[pos++] - '0');
    if ( size == 0 ) return name;
    if ( !ret.empty() ) ret += "::";
    ret += name.substr(pos, size);
    pos += size;
  }
  return ret;
}


string StringUtils::replace(string original, string from, string to){

  if (from == "" || from == to) return(original);

  string result = original;
  size_t n = result.length();
  size_t f;

  for (size_t i=0; i<n;){
    f = result.find(from,i);
    if (f == string::npos) return(result);
    result.replace(f, from.length(), to);
    i = f + to.length();
    n += to.length();
  }
  return(result);

}

