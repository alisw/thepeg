// -*- C++ -*-
//
// CFileLineReader.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_CFileLineReader_H
#define THEPEG_CFileLineReader_H
//
// This is the declaration of the CFileLineReader class.
//

#include "ThePEG/Config/ThePEG.h"
#include "CFileLineReader.fh"
#include "CFile.h"
#include <cstdio>
#include <string>

namespace ThePEG {

/**
 * CFileLineReader is a wrapper around a standard C FILE stream. With
 * it one reads one line at the time (with readline()) into an
 * internal buffer from which one can then read as from a standard
 * std::istream with a limited set of operator>> functions. It can be
 * thought of as an std::ifstream where the internal buffer must be
 * filled by hand one line at the time.
 *
 * Contrary to std::ifstream the CFileLineReader can also handle
 * gipped files and pipes. Gzipped files are automatically handles by
 * a pipe using the <code>zcat</code> command if the file name ends
 * with <code>.gz</code>. Also if a file name ends with a
 * <code>|</code> sign, the preceding string is interpreted as a
 * command defining a pipe from which to read.
 *
 * Since CFileLineReader is very close to the standard C FILE stream
 * it is in many cases much faster than eg. reading from lines via
 * std::istringstream.
 */
class CFileLineReader {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  CFileLineReader();

  /**
   * Constructor taking a \a filename as argument. Optionally the size
   * \a len of the line buffer can be specified. If \a filename ends
   * with <code>.gz</code> a pipe is opened where the file is read by
   * <code>zcat</code>. If \a filename ends with a <code>|</code>
   * sign, the preceding string is interpreted as a command defining a
   * pipe from which to read.
   */
  CFileLineReader(string filename, int len = defsize);

  /**
   * The destructor.
   */
  ~CFileLineReader();
  //@}

  /** @name Initialization functions. */
  //@{
  /**
   * Initialize with a \a filename. If \a filename ends with
   * <code>.gz</code> a pipe is opened where the file is read by
   * <code>zcat</code>. If \a filename ends with a <code>|</code>
   * sign, the preceding string is interpreted as a command defining a
   * pipe from which to read.
   */
  void open(string filename);

  /**
   * If the file was opened from within this object, close it.
   */
  void close();
  //@}

  /**
   * Read a line from the underlying c-file into the line buffer.
   */
  bool readline();

  /**
   * Undo reading from the current line, ie. the next read will be
   * from the beginning of the current line. Afterwards the state will
   * be not bad.
   */
  void resetline();

  /**
   * Return a string containing what is left of the line buffer.
   */
  string getline() const;

  /**
   * Return the underlying c-file.
   */
  CFile cFile() const;

  /**
   * Return null if a previous read failed.
   */
  operator void *();

  /**
   * Return true if a previous read failed.
   */
  bool operator!();

  /**
   * Scan forward up and until the first occurrence of the given
   * character.
   * @return true if the given character was found.
   */
  bool skip(char c);

  /**
   * Check if a given string is present in the current line buffer.
   */
  bool find(string str) const;

  /** @name Operators to read from the line buffer. */
  //@{
  /**
   * Return the next character of the line-buffer.
   */
  char getc();

  /**
   * Read a long from the line buffer.
   */
  CFileLineReader & operator>>(long & l);

  /**
   * Read an int from the line buffer.
   */
  CFileLineReader & operator>>(int & i);

  /**
   * Read an unsigned long from the line buffer.
   */
  CFileLineReader & operator>>(unsigned long & l);

  /**
   * Read an unsigned int from the line buffer.
   */
  CFileLineReader & operator>>(unsigned int & i);

  /**
   * Read a double from the line buffer.
   */
  CFileLineReader & operator>>(double & d);

  /**
   * Read a float from the line buffer.
   */
  CFileLineReader & operator>>(float & f);

  /**
   * Read a (whitespace delimited) string from the line buffer.
   */
  CFileLineReader & operator>>(std::string & s);
  //@}

private:

  /**
   * The c-file to be read from.
   */
  CFile file;

  /**
   * The length of the line buffer.
   */
  int bufflen;

  /**
   * The line buffer.
   */
  char * buff;

  /**
   * The current position in the line buffer.
   */
  char * pos;

  /**
   * The current state is bad if a read has failed.
   */
  bool bad;

  /**
   * The default size of the buffer.
   */
  static const int defsize = 1024;

private:

  /**
   * The copy constructor is private and not implemented.
   */
  CFileLineReader(const CFileLineReader &);

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  CFileLineReader & operator=(const CFileLineReader &);

};

}

#endif /* THEPEG_CFileLineReader_H */
