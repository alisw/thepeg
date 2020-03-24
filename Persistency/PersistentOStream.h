// -*- C++ -*-
//
// PersistentOStream.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_PersistentOStream_H
#define ThePEG_PersistentOStream_H
// This is the declaration of the PersistentOStream class.

#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/Utilities/ClassDescription.h"
#include "ThePEG/Utilities/Exception.h"
#include "ThePEG/Utilities/Debug.h"
#include "PersistentOStream.fh"
#include "PersistentOStream.xh"
#include <valarray>

namespace ThePEG {

/** @ingroup Persistency
 * PersistentOStream is used to write objects persistently
 * to a stream from which they can be read in again with a
 * PersistentIStream. Pointers to objects of classes
 * derived from <code>PersistentBase</code> may be written out if a
 * static ClassDescription object is present for the
 * class. Also basic types may be written to the stream, as well as
 * containers of pointers to persistent objects and basic types.
 *
 * The <code>PersistentOStream</code> keeps a list of all pointers to
 * written persistent objects, so that if several pointers to the
 * smame object is written, the object will only be written once.
 *
 * Each base class of a given object will be asked to write its
 * members to the stream starting from the least derived class going to
 * the most derived one. Members may be pointers to other persistent
 * objects or basic types or containers of these. The output for each
 * object part should be implemented by specializing the
 * ClassTraits<T>::output method, which otherwise
 * will call the non-virtual <code>persistentOutput</code> function of
 * the class. Note that for diamond-shaped multiple inheritance
 * structures, the virtual base classes will be written out several
 * times for the same object.
 *
 * @see PersistentIStream
 * @see ClassDescription
 * @see ClassTraits
 */
class PersistentOStream {

public:

  ThePEG_DECLARE_POINTERS(PersistentBase,BPtr);

  /** A map of objects indexed by integers */
  ThePEG_DECLARE_MAP(cBPtr, int, ObjectMap);

  /** A map relating class descriptions to integers. */
  ThePEG_DECLARE_MAP(const ClassDescriptionBase *, int, ClassMap);

  /** A vector of bare pointers to InputDescription objects. */
  typedef ClassDescriptionBase::DescriptionVector DescriptionVector;

public:

  /**
   * Constuctor giving an output stream. Optionally a vector of
   * libraries to be loaded before the resulting file can be read in
   * again can be given in \a libs.
   */
  PersistentOStream(ostream &, const vector<string> & libs = vector<string>());

  /**
   * Constuctor giving a file name to read. If the first
   * character in the string is a '|', the corresponding program is
   * run and its standard input is used instead. If the filename ends
   * in ".gz" the file is compressed with gzip. Optionally a vector of
   * libraries to be loaded before the resulting file can be read in
   * again can be given in \a libs.
   */
  PersistentOStream(string, const vector<string> & libs = vector<string>());

  /**
   * The destructor
   */
  ~PersistentOStream();

  /**
   * Operator for writing persistent objects to the stream.
   * @param p a pointer to the object to be written.
   * @return a reference to the stream.
   */
  template <typename T>
  PersistentOStream & operator<<(const RCPtr<T> & p) { 
    return outputPointer(p); 
  }

  /**
   * Operator for writing persistent objects to the stream.
   * @param p a pointer to the object to be written.
   * @return a reference to the stream.
   */
  template <typename T>
  PersistentOStream & operator<<(const ConstRCPtr<T> & p) { 
    return outputPointer(p); 
  }

  /**
   * Operator for writing persistent objects to the stream.
   * @param p a pointer to the object to be written.
   * @return a reference to the stream.
   */
  template <typename T>
  PersistentOStream & operator<<(const TransientRCPtr<T> & p) { 
    return outputPointer(p); 
  }

  /**
   * Operator for writing persistent objects to the stream.
   * @param p a pointer to the object to be written.
   * @return a reference to the stream.
   */
  template <typename T>
  PersistentOStream & operator<<(const TransientConstRCPtr<T> & p) {
    return outputPointer(p);
  }


  /** @name Operators for extracting built-in types from the stream. */
  //@{
  /**
   * Write a character string.
   */
  PersistentOStream & operator<<(string s) {
    for ( string::const_iterator i = s.begin(); i < s.end(); ++i ) escape(*i);
    put(tSep);
    return *this;
  }

  /**
   * Write a character.
   */
  PersistentOStream & operator<<(char c) {
    escape(c);
    put(tSep);
    return *this;
  }

  /**
   * Write a signed character.
   */
  PersistentOStream & operator<<(signed char c) {
    return (*this) << static_cast<char>(c);
  }
  
  /**
   * Write an unsigned character.
   */
  PersistentOStream & operator<<(unsigned char c) {
    return (*this) << static_cast<char>(c);
  }
  
  /**
   * Write an integer.
   */
  PersistentOStream & operator<<(int i) {
    os() << i;
    put(tSep);
    return *this;
  }

  /**
   * Write an unsigned integer.
   */
  PersistentOStream & operator<<(unsigned int i) {
    os() << i;
    put(tSep);
    return *this;
  }

  /**
   * Write a long integer.
   */
  PersistentOStream & operator<<(long i) {
    os() << i;
    put(tSep);
    return *this;
  }

  /**
   * Write an unsigned long integer.
   */
  PersistentOStream & operator<<(unsigned long i) {
    os() << i;
    put(tSep);
    return *this;
  }

  /**
   * Write a short integer.
   */
  PersistentOStream & operator<<(short i) {
    os() << i;
    put(tSep);
    return *this;
  }

  /**
   * Write an unsigned short integer.
   */
  PersistentOStream & operator<<(unsigned short i) {
    os() << i;
    put(tSep);
    return *this;
  }

  /**
   * Write a double.
   */
  PersistentOStream & operator<<(double d) {
    if ( ! isfinite(d) )
      throw WriteError()
	<< "Tried to write a NaN or Inf double to a persistent stream."
	<< Exception::runerror;
    os() << setprecision(18) << d;
    put(tSep);
    return *this;
  }

  /**
   * Write a float.
   */
  PersistentOStream & operator<<(float f) {
    if ( ! isfinite(f) )
      throw WriteError()
	<< "Tried to write a NaN or Inf float to a persistent stream."
	<< Exception::runerror;
    os() << setprecision(9) << f;
    put(tSep);
    return *this;
  }

  /**
   * Write a boolean.
   */
  PersistentOStream & operator<<(bool t) {
    if (t) put(tYes);
    else put(tNo);
    // This is a workaround for a possible bug in gcc 4.0.0
    // which inserts tYes and tNo as global symbols although
    // they are private
    //  put(t? tYes: tNo);
    put(tSep);
    return *this;
  }

  /**
   * Write a c-style character string (to be read in as a std::string).
   */
  PersistentOStream & operator<<(const char * s) {
    *this << string(s);
    return *this;
  }

  /**
   * Write a Complex.
   */
  PersistentOStream & operator<<(Complex z) {
    *this << z.real() << z.imag();
    return *this;
  }
  //@}

  /**
   * Output of containers of streamable objects.
   */
  template <typename Container>
  void putContainer(const Container & c) {
    *this << c.size();
    for ( typename Container::const_iterator it = c.begin();
	  it != c.end() && good() ; ++it )
      *this << *it;
  }

  /**
   * Write out a persistent object given a pointer to it.
   */
  PersistentOStream & outputPointer(tcBPtr);

  /**
   * For a given object, write the member variables corresponding to a
   * given ClassDescriptionBase object.
   * @param obj the object to be written.
   * @param cd a pointer to a ClassDescriptionBase describing the
   * (sub)class to written.
   */
  void putObjectPart(tcBPtr obj, const ClassDescriptionBase * cd);

  /**
   * Remove all objects that have been written, except those which are
   * to be saved, from the list of written objects.
   */
  PersistentOStream & flush();
  
  /**
   * Instuct the stream to save the following objects (protecting them from
   * being flushed).
   */
  PersistentOStream & push() {
    lastSavedObject.push(writtenObjects.size() - 1);
    return *this;
  }

  /**
   * Instuct the stream not to save the following objects.
   */
  PersistentOStream & pop() {
    lastSavedObject.pop();
    return *this;
  }

  /**
   * Check the state of the stream.
   */
  bool good() const { return !badState && os(); }

  /**
   * Check the state of the stream.
   */
  operator bool() const { return good(); }

  /**
   * Check the state of the stream.
   */
  bool operator!() const { return !good(); }

private:

  /** @cond EXCEPTIONCLASSES */
  /** @ingroup Persistency
   * Internal exception class.
   */
  struct MissingClass: public Exception {};
  /** @ingroup Persistency
   * Internal exception class.
   */
  struct WriteError: public Exception {};
  /** @endcond */

  /**
   * The version of this PersistentOStream implementation.
   */
  static const int version = 0;

  /**
   * The subversion of this PersistentOStream implementation.
   */
  static const int subVersion = 3;

  /** @name Special marker characters */
  //@{
  /**
   * The special marker character indicating the beginning of an object.
   */
  static const char tBegin = '{';

  /**
   * The special marker character indicating the end of an object.
   */
  static const char tEnd = '}';

  /**
   * The marker character indicating the beginning of the next base
   * class in case of multiple inheritance.
   */
  /**
   * The special marker character indicating an escaped marker character.
   */
  static const char tNext = '|';

  /**
   * The special marker character indicating an escaped marker character.
   */
  static const char tNull = '\\';

  /**
   * The special marker character indicating the end of a value.
   */
  static const char tSep = '\n';

  /**
   * The special marker character used to avoid confusion with escaped
   * tSep markers.
   */
  static const char tNoSep = 'n';

  /**
   * The special marker character indicating a true boolean value.
   */
  static const char tYes = 'y';

  /**
   * The special marker character indicating a false boolean value.
   */
  static const char tNo = 'n';
  //@}

  /**
   * Return true if the given character is aspecial marker character.
   */
  bool isToken(char c) const {
    return c == tBegin || c == tEnd || c == tNext || c == tSep || c == tNull;
  }

  /**
   * Set the stream in a bad state.
   */
  void setBadState() {
    breakThePEG();
    badState = true;
  }

  /**
   * Check if the state is ok.
   */
  void checkState() { if ( ! os() ) badState = true; }

  /**
   * Write out class information to the associated ostream.
   */
  const ClassDescriptionBase * writeClassId(tcBPtr);

  /**
   * write out class information to the associated ostream.
   */
  void writeClassDescription(const ClassDescriptionBase *);

  /**
   * Put a "begin object" marker on the associated ostream
   */
  void beginObject() { put(tBegin); }

  /**
   * Put a "end of object" marker on the associated ostream
   */
  void endObject() { put(tEnd); }

  /**
   * Put an "next base class" marker on the associated ostream
   */
  void endBase() { put(tNext); }

  /**
   * Put a character on the associated ostream
   */
  void put(char c) { os().put(c); }

  /**
   * Put a character on the associated ostream but escape it if it is
   * a token.
   */
  void escape(char c) { 
    if ( isToken(c) ) {
      put(tNull);
      put( c == tSep? tNoSep: c );
    } else
      put(c);
  }

  /**
   * Return a reference to the associated ostream.
   */
  ostream & os() { return *theOStream; }

  /**
   * Return a const reference to the associated ostream.
   */
  const ostream & os() const { return *theOStream; }

  /**
   * Write out initial metainfo on the stream.
   */
  void init(const vector<string> & libs);

  /**
   * List of written objects.
   */
  ObjectMap writtenObjects;

  /**
   * List of written objects that are to be saved.
   */
  stack<int> lastSavedObject;

  /**
   * List of written classes.
   */
  ClassMap writtenClasses;

  /**
   * A pointer to the associated ostream.
   */
  ostream * theOStream;

  /**
   * True if no errors has occurred.
   */
  bool badState;

  /**
   * True if the associated ostream should be deleted in the destructor.
   */
  bool allocStream;

private:

  /**
   * Standard ctors and assignment are private and not implemented.
   */
  PersistentOStream();

  /**
   * Standard ctors and assignment are private and not implemented.
   */
  PersistentOStream(const PersistentOStream &);

  /**
   * Standard ctors and assignment are private and not implemented.
   */
  PersistentOStream & operator=(const PersistentOStream &) = delete;

};

/**
 * Operator for applying manipulators to the stream.
 */
inline PersistentOStream &
operator<<(PersistentOStream & os, PersistentOManip func) {
  return (*func)(os);
}


/**
 * The manipulator for calling PersistentOStream::flush().
 */
inline PersistentOStream & flush(PersistentOStream & os) { return os.flush(); }

/**
 * The manipulator for calling PersistentOStream::push().
 */
inline PersistentOStream & push(PersistentOStream & os) { return os.push(); }

/**
 * The manipulator for calling PersistentOStream::pop().
 */
inline PersistentOStream & pop(PersistentOStream & os) { return os.pop(); }


/**
 * @name Partial specializations of operator<< for output of
 * std::containers.
 */
//@{
/** Output a pair of objects. */
template <typename T1, typename T2>
inline PersistentOStream & operator<<(PersistentOStream & os,
				      const pair<T1,T2> & p) {
  return os << p.first << p.second;
}

/**
 * Output a multimap of key/object pairs.
 */
template <typename Key, typename T, typename Cmp, typename A>
inline PersistentOStream & operator<<(PersistentOStream & os,
				      const multimap<Key,T,Cmp,A> & m) {
  os.putContainer(m);
  return os;
}

/**
 * Output a map of key/object pairs.
 */
template <typename Key, typename T, typename Cmp, typename A>
inline PersistentOStream & operator<<(PersistentOStream & os,
				      const map<Key,T,Cmp,A> & m) {
  os.putContainer(m);
  return os;
}

/**
 * Output a set of objects.
 */
template <typename Key, typename Cmp, typename A>
inline PersistentOStream & operator<<(PersistentOStream & os,
				      const set<Key,Cmp,A> & s) {
  os.putContainer(s);
  return os;
}


/**
 * Output a multiset of objects.
 */
template <typename Key, typename Cmp, typename A>
inline PersistentOStream & operator<<(PersistentOStream & os,
				      const multiset<Key,Cmp,A> & s) {
  os.putContainer(s);
  return os;
}


/**
 * Output a list of objects.
 */
template <typename T, typename A>
inline PersistentOStream & operator<<(PersistentOStream & os,
				      const list<T,A> & l) {
  os.putContainer(l);
  return os;
}


/**
 * Output a vector of objects.
 */
template <typename T, typename A>
inline PersistentOStream & operator<<(PersistentOStream & os,
              const vector<T,A> & v) {
  os.putContainer(v);
  return os;
}

/**
 * Output an array of objects.
 */
template <typename T, size_t N>
inline PersistentOStream & operator<<(PersistentOStream & os,
              const array<T,N> & a) {
  for ( auto it = a.cbegin(); it != a.cend() && os.good() ; ++it )
      os << *it;
  return os;
}

/**
 * Output a deque of objects.
 */
template <typename T, typename A>
inline PersistentOStream & operator<<(PersistentOStream & os,
				      const deque<T,A> & d) {
  os.putContainer(d);
  return os;
}

/**
 * Output a valarray.
 */
template <typename T>
inline PersistentOStream & operator<<(PersistentOStream & os,
				      const std::valarray<T> & v) {
  os << v.size();
  for ( int i = 0, N = v.size(); i < N; ++i ) os << v[i];
  return os;
}
//@}

}

#endif /* ThePEG_PersistentOStream_H */
