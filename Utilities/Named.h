// -*- C++ -*-
//
// Named.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_Named_H
#define ThePEG_Named_H
// This is the declaration of the Named class.


#include <string>

namespace ThePEG {

/**
 * The <code>Named</code> class is a simple concrete base class to
 * used by classes of objects with a name. It just defines a string
 * member variable with corresponding (protected) set and get
 * functions.
 */
class Named {

public:

  /**
   * Constructor with name.
   */
  Named(const string & newName = string()) 
    : theName(newName) {}

  /**
   *  Explicit default copy-constructor (too avoid compiler warnings)
   */
  Named(const Named & ) = default;
  
  /**
   * Return name.
   */
  const string & name() const { return theName; }

  /**
   * Test for equality.
   */
  bool operator == (const Named & other) const { 
    return theName == other.name(); 
  }

  /**
   * Lexicographical comparison.
   */
  bool operator < (const Named & other) const { 
    return theName < other.name(); 
  }

protected:

  /**
   * Assignment.
   */
  const Named & operator = (const Named & other) {
    if (this != &other)
      theName = other.name(); 
    return *this; 
  }

  /**
   * Set new name.
   */
  const string & name(const string & newName) { 
    return theName = newName; 
  } 

private:

  /**
   * The string containing the name.
   */
  string theName;

};

}

#endif /* ThePEG_Named_H */
