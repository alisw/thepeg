// -*- C++ -*-
//
// Command.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_Command_H
#define ThePEG_Command_H
// This is the declaration of the Command and CommandBase classes.

#include "ThePEG/Config/ThePEG.h"
#include "InterfaceBase.h"
#include "Command.fh"
#include "Command.xh"

namespace ThePEG {

/**
 * The CommandBase and its templated Command sub-class defines an
 * interface to a class derived from the InterfacedBase, through which
 * arbitratry command strings can be sent and
 * received. <code>Command</code> is templated and is derived from the
 * InterfaceBase class via <code>CommandBase</code>.
 *
 * For each command interface to be defined for a class
 * <code>T</code>, exactly one static object of the Command<T> must be
 * created and initialized as follows:
 *
 * <code>Command<T> comint(name, description, &T::memberfn,
 * depsafe);</code>
 *
 * Where <code>name</code> is an identifier <code>std::string</code> which
 * should only contain letters [a-zA-z0-9_], <code>description</code> is
 * an arbitrary <code>std::string</code>, <code>memberfn</code> should be
 * a non-static member function of <code>T</code> and defined as
 * <code>std::string T::memberfn(std::string)</code>. Finally if
 * <code>depsafe</code> is true it can be assumed that a call to the
 * <code>memberfn</code> for an object does not influence other objects
 * which may depend on the first.
 *
 * The <code>Command</code> class, as all other
 * InterfaceBase classes are mainly used in the
 * BaseRepository class.
 *
 *
 * @see InterfacedBase
 * @see InterfaceBase
 * @see BaseRepository
 * 
 */
class CommandBase: public InterfaceBase {

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
   * @param depSafe set to true if calls to this interface for one
   * object does not influence other objects.
   */
  CommandBase(string newName, string newDescription, string newClassName,
		     const type_info & newTypeInfo, bool depSafe)
    : InterfaceBase(newName, newDescription, newClassName, 
		    newTypeInfo, depSafe, false) {
    hasDefault = false;
  }

  /**
   * The general interface method overriding the one in
   * InterfaceBase. For this class, the \a action and \a argument
   * arguments are concatenated (with a space character inserted) and
   * sent to the cmd() function.
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

  /**
   * Execute the member function. For the object \a ib execute the
   * memberfunction (defined in the derived class) with \a c as
   * argument and return the return value.
   */
  virtual string cmd(InterfacedBase & ib, string c) const
    = 0;

};

/**
 * The CommandBase and its templated Command sub-class defines an
 * interface to a class derived from the InterfacedBase, through which
 * arbitratry command strings can be sent and
 * received. <code>Command</code> is templated and is derived from the
 * InterfaceBase class via <code>CommandBase</code>.
 *
 * For each command interface to be defined for a class
 * <code>T</code>, exactly one static object of the Command<T> must be
 * created and initialized as follows:
 *
 * <code>Command<T> comint(name, description, &T::memberfn,
 * depsafe);</code>
 *
 * Where <code>name</code> is an identifier <code>std::string</code> which
 * should only contain letters [a-zA-z0-9_], <code>description</code> is
 * an arbitrary <code>std::string</code>, <code>memberfn</code> should be
 * a non-static member function of <code>T</code> and defined as
 * <code>std::string T::memberfn(std::string)</code>. Finally if
 * <code>depsafe</code> is true it can be assumed that a call to the
 * <code>memberfn</code> for an object does not influence other objects
 * which may depend on the first.
 *
 * The <code>Command</code> class, as all other
 * InterfaceBase classes are mainly used in the
 * BaseRepository class.
 *
 *
 * @see InterfacedBase
 * @see InterfaceBase
 * @see BaseRepository
 * 
 */
template <class T>
class Command: public CommandBase {

public:

  /**
   * The declaration of member functions which can be used by this
   * Command interface.
   */
  typedef string (T::*ExeFn)(string);

public:

  /**
   * Standard constructor.
   *
   * @param newName the name of the interface, may only contain
   * letters [a-zA-z0-9_].
   *
   * @param newDescription a brief description of the interface.
   *
   * @param newExeFn pointer to the function to be called in the
   * corresponding class.
   *
   * @param depSafe set to true if calls to this interface for one
   * object does not influence other objects.
   */
  Command(string newName, string newDescription,
	  ExeFn newExeFn, bool depSafe = false)
    : CommandBase(newName, newDescription, 
		  ClassTraits<T>::className(), typeid(T), depSafe), 
      theExeFn(newExeFn) {}

  /**
   * Execute the member function. For the object \a ib execute the
   * memberfunction with \a c as argument and return the return value.
   */
  virtual string cmd(InterfacedBase & ib, string)
    const;


private:

  /**
   * The pointer to the member function.
   */
  ExeFn theExeFn;

};

}

#ifndef ThePEG_TEMPLATES_IN_CC_FILE
#include "Command.tcc"
#endif

#endif /* ThePEG_Command_H */
