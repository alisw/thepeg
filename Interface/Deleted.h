// -*- C++ -*-
//
// Deleted.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_Deleted_H
#define ThePEG_Deleted_H
// This is the declaration of the Deleted and DeletedBase classes.

#include "ThePEG/Config/ThePEG.h"
#include "InterfaceBase.h"

namespace ThePEG {

/**
 * The DeletedBase and its templated Deleted sub-class defines an
 * interface to a class derived from the InterfacedBase. It should be
 * used when an interface is removed to provide a user-friendly
 * message indicating why it was removed and possibly which interface
 * should be used instead.
 *
 * For each deleted interface to be defined for a class
 * <code>T</code>, exactly one static object of the Deleted<T> must be
 * created and initialized as follows:
 *
 * <code>Deleted<T> delint(name, description);</code>
 *
 * Where <code>name</code> is an identifier <code>std::string</code> which
 * should only contain letters [a-zA-z0-9_] and <code>description</code> is
 * an arbitrary <code>std::string</code>
 *
 * The <code>Deleted</code> class, as all other InterfaceBase classes
 * are mainly used in the BaseRepository class.
 *
 *
 * @see InterfacedBase
 * @see InterfaceBase
 * @see BaseRepository
 * 
 */
class DeletedBase: public InterfaceBase {

public:

  /**
   * Standard constructor.
   *
   * @param newName the name of the interface, may only contain
   * letters [a-zA-z0-9_].
   *
   * @param newDescription a brief description of the interface.
   *
   * @param newClassName the name of the corresponding class.
   *
   * @param newTypeInfo the type_info object of the corresponding
   * class.
   *
   */
  DeletedBase(string newName, string newDescription, string newClassName,
		     const type_info & newTypeInfo)
    : InterfaceBase(newName, newDescription, newClassName, 
		    newTypeInfo, true, false) {
    rank(-1.0e10);
  }

  /**
   * The general interface method overriding the one in
   * InterfaceBase. For this class, an exception will be thrown with a
   * message given by the description string provided in the
   * constructor.
   */
  virtual string
  exec(InterfacedBase &ib, string action, string arguments) const
   ;

  /**
   * Return a string describing the type of interface to be included
   * in the Doxygen documentation.
   */
  virtual string doxygenType() const;

  /**
   * Return a code for the type of this interface.
   */
  virtual string type() const;

};

/**
 * The DeletedBase and its templated Deleted sub-class defines an
 * interface to a class derived from the InterfacedBase. It should be
 * used when an interface is removed to provide a user-friendly
 * message indicating why it was removed and possibly which interface
 * should be used instead.
 *
 * For each deleted interface to be defined for a class
 * <code>T</code>, exactly one static object of the Deleted<T> must be
 * created and initialized as follows:
 *
 * <code>Deleted<T> delint(name, description);</code>
 *
 * Where <code>name</code> is an identifier <code>std::string</code> which
 * should only contain letters [a-zA-z0-9_] and <code>description</code> is
 * an arbitrary <code>std::string</code>
 *
 * The <code>Deleted</code> class, as all other InterfaceBase classes
 * are mainly used in the BaseRepository class.
 *
 *
 * @see InterfacedBase
 * @see InterfaceBase
 * @see BaseRepository
 * 
 */
template <class T>
class Deleted: public DeletedBase {

public:

  /**
   * Standard constructor.
   *
   * @param newName the name of the interface, may only contain
   * letters [a-zA-z0-9_].
   *
   * @param newDescription a brief description of the interface.
   *
   */
  Deleted(string newName, string newDescription)
    : DeletedBase(newName, newDescription, 
		  ClassTraits<T>::className(), typeid(T)) {}

};

}

#endif /* ThePEG_Deleted_H */
