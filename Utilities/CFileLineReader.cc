// -*- C++ -*-
//
// CFileLineReader.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the CFileLineReader class.
//

#include "CFileLineReader.h"
#include "config.h"
#include <cstdlib>

using namespace ThePEG;

CFileLineReader::CFileLineReader()
  : bufflen(defsize), buff(new char[defsize]), pos(buff), bad(false) {}

CFileLineReader::CFileLineReader(string filename, int len)
  : bufflen(len), buff(new char[len]), pos(buff), bad(false) {
  open(filename);
}

CFileLineReader::~CFileLineReader() {
  close();
}

void CFileLineReader::open(string filename) {
  close();
  file = CFile(filename, "r");
  pos = buff;
  bad = !file;
}

void CFileLineReader::close() {
  file.close();
  bad = true;
}

bool CFileLineReader::readline() {
  bad =  ( file.gets(buff, bufflen) == 0 );
  pos = buff;
  return !bad;
}


CFileLineReader::operator void *() {
  // the returned address is not meant to be dereferenced,
  // it only indicates success.
  return bad ? 0 : &bad;
}

bool CFileLineReader::operator!() {
  return bad;
}

 char CFileLineReader::getc() {
  char c = *pos++;
  bad = ( *pos == 0 );
  return c;
}

bool CFileLineReader::skip(char c) {
  while ( *pos != 0 && *pos++ != c ) {}
  bad = ( *pos == 0 );
  return !bad;
}

bool CFileLineReader::find(string str) const {
  return ( string(pos).find(str) != string::npos );
}

std::string CFileLineReader::getline() const {
  return std::string(pos);
}

CFileLineReader & CFileLineReader::operator>>(long & l) {
  char * next;
  l = std::strtol(pos, &next, 0);
  bad = ( next == pos );
  pos = next;
  return *this;
}

 CFileLineReader & CFileLineReader::operator>>(int & i) {
  char * next;
  i = int(std::strtol(pos, &next, 0));
  bad = ( next == pos );
  pos = next;
  return *this;
}

CFileLineReader & CFileLineReader::operator>>(unsigned long & l) {
  char * next;
  l = std::strtoul(pos, &next, 0);
  bad = ( next == pos );
  pos = next;
  return *this;
}

CFileLineReader & CFileLineReader::operator>>(unsigned int & i) {
  char * next;
  i = static_cast<unsigned int>(std::strtoul(pos, &next, 0));
  bad = ( next == pos );
  pos = next;
  return *this;
}

CFileLineReader & CFileLineReader::operator>>(double & d) {
  char * next;
  d = std::strtod(pos, &next);
  bad = ( next == pos );
  pos = next;
  // fortran formatted doubles
  if(!bad && ( *pos == 'd' || *pos == 'D')) {
    ++pos;
    int exp;
    *this >> exp;
    d *= pow(10.,exp);
  }
  return *this;
}

CFileLineReader & CFileLineReader::operator>>(float & f) {
  char * next;
  f = float(std::strtod(pos, &next));
  bad = ( next == pos );
  pos = next;
  // fortran formatted doubles
  if(!bad && ( *pos == 'd' || *pos == 'D')) {
    ++pos;
    int exp;
    *this >> exp;
    f *= pow(10.,exp);
  }
  return *this;
}

CFileLineReader & CFileLineReader::operator>>(std::string & s) {
  while ( *pos != 0 && std::isspace(*pos) ) ++pos;
  bad = ( *pos == 0 );
  char * first = pos;
  while ( *pos != 0 && !std::isspace(*pos) ) ++pos;
  s.assign(first, pos);
  return *this;
}

CFile CFileLineReader::cFile() const {
  return file;
}

void CFileLineReader::resetline() {
  bad = false;
  pos = buff;
}
