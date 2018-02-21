// -*- C++ -*-
//
// InputDescription.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_InputDescription_H
#define ThePEG_InputDescription_H
// This is the declaration of the InputDescription class.

#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/Utilities/ClassDescription.h"

namespace ThePEG {

/** @ingroup Persistency
 * InputDescription objects are used by the PersistentIStream class to
 * keep track of all classes it has read from a stream. It keeps a
 * pointer to the corresponding ClassDescription in case the class
 * read in was actually present in the current program, a version
 * number of the read class which may be different from the class
 * present in the current program and a list of base class
 * <code>InputDescription</code>s.
 *
 * @see PersistentIStream
 * @see Named
 * @see ClassDescription
 */
class InputDescription: public Named {

public:

  /** A vector of pointers to InputDescription objects. */
  typedef vector<const InputDescription *> DescriptionVector;

  ThePEG_DECLARE_POINTERS(PersistentBase,BPtr);

public:

  /**
   * The standard constructor.
   * @param newName the name of the class being read.
   * @param newVersion the version number of the class when written.
   */
  InputDescription(string newName, int newVersion) 
    : Named(newName), theDescription(0), theVersion(newVersion) {}

  /**
   * Set the ClassDescriptionBase object of the class being read.
   */
  void setDescription(const ClassDescriptionBase * cd) {
    theDescription = cd;
  }

  /**
   * Add a base class description.
   */
  void addBaseClass(const InputDescription * newBase) {
    theBaseClasses.push_back(newBase);
  }

  /**
   * Return the list of base class descriptions.
   */
  const DescriptionVector & descriptions() const {
    return theBaseClasses;
  }

  /**
   * Create an object of the corresponding class.
   */
  BPtr create() const {
    if ( theDescription ) return theDescription->create();
    DescriptionVector::const_iterator dit = theBaseClasses.begin();
    while ( dit != theBaseClasses.end() ) {
      BPtr obj = (*dit++)->create();
      if ( obj ) return obj;
    }
    return BPtr();
  }

  /**
   * Read an object part of the corresponding class from a stream.
   * Will only read the part of the object corresponding to the
   * members of the class represented by this object.
   */
  void input(tBPtr b, PersistentIStream & is) const {
    if ( theDescription ) theDescription->input(b, is, theVersion);
  }

private:

  /**
   * The list of base class descriptions.
   */
  DescriptionVector theBaseClasses;

  /**
   * The description of the corresponding class in the current
   * program.
   */
  const ClassDescriptionBase * theDescription;

  /**
   * The version of the class to be read.
   */
  int theVersion;

};

}

#endif /* ThePEG_InputDescription_H */
