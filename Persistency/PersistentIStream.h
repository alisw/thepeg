// -*- C++ -*-
//
// PersistentIStream.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_PersistentIStream_H
#define ThePEG_PersistentIStream_H
// This is the declaration of the PersistentIStream class.

#include "ThePEG/Config/ThePEG.h"
#include "InputDescription.h"
#include "PersistentIStream.fh"
#include "ThePEG/Utilities/Exception.h"
#include <climits>
#include <valarray>

namespace ThePEG {

/** @ingroup Persistency
 * PersistentIStream is used to read persistent objects from a stream
 * where they were previously written using PersistentOStream. Basic
 * types and pointers to objects derived from
 * <code>PersistentBase</code> should be read in the same order they
 * were written out. If <code>pedantic()</code> is true the same
 * classes that were written out must be present in the current
 * program. If <code>pedantic()</code> is false and if an object is
 * read for which only a base class is present in the current program,
 * only the parts corresponding to the base class will be read, and
 * the rest will be gracefully skipped.
 *
 * Each base class of a given object will be asked to read its
 * members from the stream starting from the least derived class going to
 * the most derived one. Members may be pointers to other persistent
 * objects or basic types or containers of these. The output for each
 * object part should be implemented by specializing the
 * ClassTraits<T>::input method, which otherwise
 * will call the non-virtual <code>persistentInput</code> function of
 * the class. Note that for diamond-shaped multiple inheritance
 * structures, the virtual base classes will be written out several
 * times for the same object.
 *
 * @see PersistentOStream
 * @see ClassTraits
 */
class PersistentIStream {

public:

  ThePEG_DECLARE_POINTERS(PersistentBase,BPtr);

  /** A vector of pointers to persistent objects */
  typedef vector<BPtr> ObjectVector;

  /** A vector of bare pointers to InputDescription objects. */
  typedef InputDescription::DescriptionVector DescriptionVector;

public:

  /**
   * Constuctor giving an input stream to be used as an underlying
   * istream.
   */
  PersistentIStream(istream & is) 
    : theIStream(&is), isPedantic(true), 
      allocStream(false), badState(false)
  {
    init();
  }



  /**
   * Constuctor giving a file name to read from. If the first
   * character in the string is a '|', the corresponding program is
   * run and its standard output is used instead. If the filename ends
   * in ".gz" the file is uncompressed with gzip.
   */
  PersistentIStream(string);

  /**
   * The destructor.
   */
  ~PersistentIStream();

  /**
   * Operator for extracting persistent objects from the stream.
   * @param ptr this pointer will refer to the extracted object.
   * @return a reference to the stream.
   */
  template <typename T>
  PersistentIStream & operator>>(RCPtr<T> & ptr) {
    BPtr b = getObject();
    ptr = dynamic_ptr_cast< RCPtr<T> >(b);
    if ( b && !ptr ) setBadState();
    return *this;
  }

  /**
   * Operator for extracting persistent objects from the stream.
   * @param ptr this pointer will refer to the extracted object.
   * @return a reference to the stream.
   */
  template <typename T>
  PersistentIStream & operator>>(ConstRCPtr<T> & ptr) {
    BPtr b = getObject();
    ptr = dynamic_ptr_cast< ConstRCPtr<T> >(b);
    if ( b && !ptr ) setBadState();
    return *this;
  }

  /**
   * Operator for extracting persistent objects from the stream.
   * @param ptr this pointer will refer to the extracted object.
   * @return a reference to the stream.
   */
  template <typename T>
  PersistentIStream & operator>>(TransientRCPtr<T> & ptr) {
    BPtr b = getObject();
    ptr = dynamic_ptr_cast< TransientRCPtr<T> >(b);
    if ( b && !ptr ) setBadState();
    return *this;
  }

  /**
   * Operator for extracting persistent objects from the stream.
   * @param ptr this pointer will refer to the extracted object.
   * @return a reference to the stream.
   */
  template <typename T>
  PersistentIStream & operator>>(TransientConstRCPtr<T> & ptr) {
    BPtr b = getObject();
    ptr = dynamic_ptr_cast< TransientConstRCPtr<T> >(b);
    if ( b && !ptr ) setBadState();
    return *this;
  }


  /** @name Operators for extracting built-in types from the stream. */
  //@{
  /**
   * Read a character string.
   */
  PersistentIStream & operator>>(string &);

  /**
   * Read a character.
   */
  PersistentIStream & operator>>(char &);

  /**
   * Read a signed character.
   */
  PersistentIStream & operator>>(signed char &);

  /**
   * Read an unsigned character.
   */
  PersistentIStream & operator>>(unsigned char &);

  /**
   * Read an integer.
   */
  PersistentIStream & operator>>(int & i) {
    is() >> i;
    getSep();
    return *this;
  }

  /**
   * Read an unsigned integer.
   */
  PersistentIStream & operator>>(unsigned int & i) {
    is() >> i;
    getSep();
    return *this;
  }

  /**
   * Read a long integer.
   */
  PersistentIStream & operator>>(long & i) {
    is() >> i;
    getSep();
    return *this;
  }

  /**
   * Read an unsigned long integer.
   */
  PersistentIStream & operator>>(unsigned long & i) {
    is() >> i;
    getSep();
    return *this;
  }

  /**
   * Read a short integer.
   */
  PersistentIStream & operator>>(short & i) {
    is() >> i;
    getSep();
    return *this;
  }

  /**
   * Read an unsigned short integer.
   */
  PersistentIStream & operator>>(unsigned short & i) {
    is() >> i;
    getSep();
    return *this;
  }

  /**
   * Read a double.
   */
  PersistentIStream & operator>>(double & d) {
    is() >> d;
    getSep();
    return *this;
  }

  /**
   * Read a float.
   */
  PersistentIStream & operator>>(float & f) {
    is() >> f;
    getSep();
    return *this;
  }

  /**
   * Read a bool.
   */
  PersistentIStream & operator>>(bool &);

  /**
   * Read a Complex.
   */
  PersistentIStream & operator>>(Complex &);
  //@}

  /**
   * Intput of containers streamable objects.
   * @param c the container into which objects are added.
   */
  template <typename Container> void getContainer(Container & c) {
    long size;
    typename Container::value_type val;
    c.clear();
    *this >> size;
    while ( size-- && good() ) {
      *this >> val;
      c.insert(c.end(), val);
    }
  }

  /**
   * Read in an object. Create an object and read its data from the
   * stream.
   * @return a pointer to the read object.
   */
  BPtr getObject();

  /**
   * For a given object, read the member variables corresponding to a
   * given InputDescription object.
   * @param obj the object to be read into.
   * @param pid a pointer to an InputDescription describing the
   * (sub)class to be read.
   */
  void getObjectPart(tBPtr obj, const InputDescription * pid);

  /**
   * Read a class description from the underlying stream and return a
   * corresponding InputDescription object
   */
  const InputDescription * getClass();
  
  /**
   * Set pedantic mode.  If the stream is set to be tolerant it is
   * allowed to read an object from the stream even if the
   * corresponding class is not known to the running executable, under
   * the condition that a public base class of the unknown class is
   * known. If the stream is set to be pedantic this is not allowed.
   * By default, the stream is pedantic.
   */
  PersistentIStream & setPedantic() {
    isPedantic = true;
    return *this;
  }

  /**
   * Set tolerant mode.  If the stream is set to be tolerant it is
   * allowed to read an object from the stream even if the
   * corresponding class is not known to the running executable, under
   * the condition that a public base class of the unknown class is
   * known. If the stream is set to be pedantic this is not allowed.
   * By default, the stream is pedantic.
   */
  PersistentIStream & setTolerant() {
    isPedantic = false;
    return *this;
  }

  /**
   * Check the state of the stream.
   */
  bool good() const { return !badState && is(); }
  
  /**
   * Check the state of the stream.
   */
  bool operator!() const { return !good(); }

  /**
   * Check the state of the stream.
   */
  operator bool() const { return good(); }

  /**
   * Check the tolerance. Returns true if setPedantic() has been
   * called or if not setTolerant() has been called.
   */
  bool pedantic() const { return isPedantic; }

  /**
   * The global libraries loaded on initialization.
   */
  const vector<string> & globalLibraries() const {
    return theGlobalLibraries;
  }

private:

  /** @cond EXCEPTIONCLASSES */
  /** @ingroup Persistency
      Thrown if a class is missing */
  struct MissingClass: public Exception {};

  /** @ingroup Persistency Thrown if an object which should have been
      read in is missing. */
  struct MissingObject: public Exception {};

  /** @ingroup Persistency Thrown if reading from the stream failed
      for some reason. */
  struct ReadFailure: public Exception {};
  /** @endcond */

  /**
   * Internal initialization.
   */
  void init();

  /**
   * Get the next character from the associated istream.
   */
  char get() { return is().get(); }

  /**
   * Get the next character from the associated istream and decode it
   * if it is escaped.
   */
  char escaped() {
    char c = get();
    return c == tNoSep? tSep: c;
  }

  /**
   * Set the stream in a bad state
   */
  void setBadState() {
    breakThePEG();
    badState = true;
  }

  /**
   * Read a field separator from the stream.
   */
  void getSep() {
    if ( !pedantic() ) skipField();
    else if ( get() != tSep ) setBadState();
  }

  /**
   * Scan the stream for the next field separator.
   */
  void skipField() {
    is().ignore(INT_MAX, tSep);
    if ( !is() ) setBadState();
  }


  /**
   * Check if the next char to be read is a tBegin marker.
   */
  bool beginObject() { return is().peek() == tBegin; }

  /**
   * Scan the stream to the end of the current object. If any new object are
   * found these are read prom the stream to ensure that the pointer structure
   * is preserved.
   */
  void endObject();

  /**
   * Scan stream for "end base class" marker. The \a classname is the
   * name of the class currently being read and is only used for
   * documenting exceptions.
   */
  void endBase(string classname);

  /**
   * Return a reference to the associated stream.
   */
  istream & is() { return *theIStream; }

  /**
   * Return a const reference to the associated stream.
   */
  const istream & is() const { return *theIStream; }

  /**
   * Lists of objects that have been read.
   */
  ObjectVector readObjects;

  /**
   * Lists of classes and corresponding version strings that have been read.
   */
  DescriptionVector readClasses;

  /**
   * A pointer to the associated istream.
   */
  istream * theIStream;

  /**
   * Pedantic or tolerant. See description of the setPedantic() and
   * setTolerant() methods.
   */
  bool isPedantic;

  /**
   * True if the associated istream should be deleted when the PersistentIStream
   * is destroyed.
   */
  bool allocStream;

  /**
   * False if no errors has occurred.
   */
  bool badState;

  /** Version number of the PersistentOStream which has written the
   *  file being read. */
  int version;

  /** Subversion number of the PersistentOStream which has written the
   *  file being read. */
  int subVersion;

  /**
   * Global libraries loaded in the initialization.
   */
  vector<string> theGlobalLibraries;

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

private:

  /**
   * Standard ctors and assignment are private and not implemented.
   */
  PersistentIStream();

  /**
   * Standard ctors and assignment are private and not implemented.
   */
  PersistentIStream(const PersistentIStream &);

  /**
   * Standard ctors and assignment are private and not implemented.
   */
  PersistentIStream & operator=(const PersistentIStream &) = delete;

};


/**
 * Operator for applying manipulators to the stream.
 */
inline PersistentIStream & operator>>(PersistentIStream & is, 
				      PersistentIManip func) {
  return (*func)(is);
}
  
/**
 * The manipulator for setting pedantic mode.
 */
inline PersistentIStream & pedantic(PersistentIStream & is) {
  return is.setPedantic();
}


/**
 * The manipulator for setting tolerant mode.
 */
inline PersistentIStream & tolerant(PersistentIStream & is) {
  return is.setTolerant();
}


/**
 * @name Partial specializations of operator>> for input of
 * std::containers.
 */
//@{
/** Input a pair of objects. */
template <typename T1, typename T2>
inline PersistentIStream & operator>>(PersistentIStream & is, pair<T1,T2> & p) {
  return is >> p.first >> p.second;
}

/** Input a map of key/objects pairs. */
template <typename Key, typename T, typename Cmp, typename A>
inline PersistentIStream & operator>>(PersistentIStream & is, map<Key,T,Cmp,A> & m) {
  m.clear();
  long size;
  Key k;
  is >> size;
  while ( size-- && is ) {
    is >> k;
    is >> m[k];
  }
  return is;
}

/** Input a multimap of key/objects pairs. */
template <typename Key, typename T, typename Cmp, typename A>
inline PersistentIStream & operator>>(PersistentIStream & is,
				      multimap<Key,T,Cmp,A> & m) {
  m.clear();
  long size;
  Key k;
  T t;
  is >> size;
  while ( size-- && is ) {
    is >> k;
    is >> t;
    m.insert(make_pair(k, t));
  }
  return is;
}


/** Input a set of objects. */
template <typename Key, typename Cmp, typename A>
inline PersistentIStream & operator>>(PersistentIStream & is, 
				      set<Key,Cmp,A> & s) {
  is.getContainer(s);
  return is;
}

/** Input a multoset of objects. */
template <typename Key, typename Cmp, typename A>
inline PersistentIStream & operator>>(PersistentIStream & is,
				      multiset<Key,Cmp,A> & s) {
  is.getContainer(s);
  return is;
}


/** Input a list of objects. */
template <typename T, typename A>
inline PersistentIStream & operator>>(PersistentIStream & is, 
				      list<T,A> & l) {
  is.getContainer(l);
  return is;
}


/** Input a vector of objects. */
template <typename T, typename A>
inline PersistentIStream & operator>>(PersistentIStream & is, 
              vector<T,A> & v) {
  is.getContainer(v);
  return is;
}

/** Input an array of objects. */
template <typename T, size_t N>
inline PersistentIStream & operator>>(PersistentIStream & is, 
              array<T,N> & a) {
  for ( size_t i = 0; i < N && is.good(); ++i ) 
    is >> a[i];
  return is;
}

/** Input a deque of objects. */
template <typename T, typename A>
inline PersistentIStream & operator>>(PersistentIStream & is, 
				      deque<T,A> & d) {
  is.getContainer(d);
  return is;
}

/** Input a valarray. */
template <typename T>
inline PersistentIStream & operator>>(PersistentIStream & is, 
				      std::valarray<T> & v) {
  long size;
  is >> size;
  v = std::valarray<T>(size);
  for ( int i = 0; i < size && is.good(); ++i ) is >> v[i];
  return is;
}

}

#endif /* ThePEG_PersistentIStream_H */
