// -*- C++ -*-
//
// ParVector.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_ParVector_H
#define ThePEG_ParVector_H
// This is the declaration of the ParVector, ParVectorTBase and
// ParVectorBase classes.

#include "ThePEG/Config/ThePEG.h"
#include "InterfaceBase.h"
#include "ParVector.fh"
#include <limits>

namespace ThePEG {

/// Helper functions for putUnit
namespace {
  template <typename T>
  /**
   * Helper function for putUnit
   */
  inline void putUnitImpl2(ostream & os, T v, T u, DimensionT) {
    os << v/u;
  }
  
  template <typename T>
  /**
   * Helper function for putUnit
   */
  inline void putUnitImpl2(ostream & os, T v, T u, StandardT) {
    if ( u > T() )
      os << v/u;
    else
      os << v;
  }
}

/**
 * The ParVector and its base classes ParVectorTBase and ParVectorBase
 * defines an interface to a class derived from the InterfacedBase,
 * through which vectors (or any other container) of simple member
 * variables can be manuipulated.  ParVector is templated on the type
 * of the member variable and the type of the InterfacedBase class,
 * and is derived from the InterfaceBase class via
 * ParVectorTBase (which is templated only on the class
 * of the member variable) and ParVectorBase.
 *
 * For each InterfacedBase class exactly one static ParVector object
 * should created for each member variable of container type which
 * should be interfaced. This object will automatically register
 * itself with the BaseRepository class.
 *
 * @see Parameter
 * @see InterfacedBase
 * @see InterfaceBase
 * 
 */
class ParVectorBase: public InterfaceBase {

public:

  /** A vector of strings. */
  typedef vector<string> StringVector;

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
   * @param newSize the size of the container or -1 if varying.
   *
   * @param depSafe set to true if calls to this interface for one
   * object does not influence other objects.
   *
   * @param readonly if this is set true the interface will not be
   * able to manipulate objects of the corresponding class, but will
   * still be able to access information.
   *
   * @param limits determines if the values of the parameters are
   * limited from above and/or below. The possible values are given by
   * Interface::Limits.
   */
  ParVectorBase(string newName, string newDescription,
		string newClassName,
		const type_info & newTypeInfo, int newSize,
		bool depSafe, bool readonly, int limits) 
    : InterfaceBase(newName, newDescription, newClassName, 
		    newTypeInfo, depSafe,
		    readonly), limit(limits), theSize(newSize) {
    hasDefault = false;
  }

  /**
   * Destructor.
   */
  virtual ~ParVectorBase() {}

  /**
   * The general interface method overriding the one in
   * InterfaceBase. For this class, \a action can be any of "set",
   * "insert", "erase", "get", "min", "max", "def" and "setdef" and \a
   * argument should be a something which can be read into an integer
   * and into a variable through a stringstream with the standard '>>'
   * operator.
   */
  virtual string exec(InterfacedBase &, string action,
		      string arguments) const;

  /**
   * Return a complete description of this parameter vector.
   */
  virtual string fullDescription(const InterfacedBase & ib) const;

  /**
   * Set the \a i'th element of a container of member variables of \a
   * ib to \a val. \a val should be able to be read into the
   * corresponding variable type through a stringstream with the
   * standard '>>' operator.
   */
  virtual void set(InterfacedBase & ib, string val, int i)
    const = 0;

  /**
   * Insert a new object before the \a i'th element of a container of
   * member variables of \a ib and set it to \a val. \a val should be
   * able to be read into the corresponding variable type through a
   * stringstream with the standard '>>' operator.
   */
  virtual void insert(InterfacedBase & ib, string val, int i)
    const = 0;

  /**
   * Remove the \a i'th element of a container of member variables of
   * \a ib.
   */
  virtual void erase(InterfacedBase & ib, int i)
    const = 0;

  /**
   * Clear the container of pointers of \a ib.
   */
  virtual void clear(InterfacedBase & ib)
    const = 0;

  /**
   * Return the values of a container of member variables of \a ib.
   */
  virtual StringVector get(const InterfacedBase & ib) const
    = 0;

  /**
   * Return the minimum value allowed for the \a i'th element of a
   * container of member variables of \a ib.
   */
  virtual string minimum(const InterfacedBase & ib, int i) const
    = 0;

  /**
   * Return the maximum value allowed for the \a i'th element of a
   * container of member variables of \a ib.
   */
  virtual string maximum(const InterfacedBase & ib, int i) const
    = 0;

  /**
   * Return the default value for the \a i'th element of a container
   * of member variables of \a ib.
   */
  virtual string def(const InterfacedBase & ib, int i) const
    = 0;

  /**
   * Return the general default value for this parameter vector.
   */
  virtual string def() const = 0;

  /**
   * Set the \a i'th element of a container of member variables of \a
   * ib to its default value.
   */
  virtual void setDef(InterfacedBase & ib, int i) const
    = 0;

  /**
   * True if there the variable is limited from above and below.
   */
  bool limited() const { return limit != Interface::nolimits; }

  /**
   * True if there the variable is limited from abovew.
   */
  bool upperLimit() const { 
    return limit == Interface::limited || limit == Interface::upperlim;
  }

  /**
   * True if there the variable is limited from  below.
   */
  bool lowerLimit() const { 
    return limit == Interface::limited || limit == Interface::lowerlim;
  }

  /**
   * Set a flag indicating that there are limits associated with the
   * variables.
   */
  void setLimited() { limit = Interface::limited; }

  /**
   * Set a flag indicating if there are no limits associated with the
   * variables.
   */
  void setUnlimited() { limit = Interface::nolimits; }

  /**
   * Get the size of the container being interfaced. If the size() is
   * less than 0, the size is allowed to vary.
   */
  int size() const { return theSize; }

  /**
   * Set the size of the container being interfaced. If the size is
   * less than 0, the size is allowed to vary.
   */
  void setSize(int sz) { theSize = sz; }

  /**
   * Set the size of the container being interfaced to -1, i.e. the
   * size is allowed to vary.
   */
  void setVariableSize() { theSize = 0; }

private:

  /**
   * True if there are limits associated with the
   * variables.
   */
  int limit;

  /**
   * The size of the container being interfaced.
   */
  int theSize;

};

/**
 * The ParVector and its base classes ParVectorTBase and ParVectorBase
 * defines an interface to a class derived from the InterfacedBase,
 * through which vectors (or any other container) of simple member
 * variables can be manuipulated.  ParVector is templated on the type
 * of the member variable and the type of the InterfacedBase class,
 * and is derived from the InterfaceBase class via
 * ParVectorTBase (which is templated only on the class
 * of the member variable) and ParVectorBase.
 *
 * For each InterfacedBase class exactly one static ParVector object
 * should created for each member variable of container type which
 * should be interfaced. This object will automatically register
 * itself with the BaseRepository class.
 *
 * @see Parameter
 * @see InterfacedBase
 * @see InterfaceBase
 * 
 */
template <typename Type>
class ParVectorTBase: public ParVectorBase {

public:

  /** A vector of objects of the template argument type. */
  typedef vector<Type> TypeVector;

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
   * @param newUnit the unit assumed when a number is read or written
   * to a stream.
   *
   * @param newSize the size of the container or -1 if varying.
   *
   * @param depSafe set to true if calls to this interface for one
   * object does not influence other objects.
   *
   * @param readonly if this is set true the interface will not be
   * able to manipulate objects of the corresponding class, but will
   * still be able to access information.
   *
   * @param limits determines if the values of the parameters are
   * limited from above and/or below. The possible values are given by
   * Interface::Limits.
   */
  ParVectorTBase(string newName, string newDescription,
		 string newClassName, const type_info & newTypeInfo,
		 Type newUnit, int newSize, bool depSafe,
		 bool readonly, int limits)
    : ParVectorBase(newName, newDescription, newClassName, 
		    newTypeInfo, newSize,
		    depSafe, readonly, limits), theUnit(newUnit) {}


  /**
   * Destructor.
   */
  virtual ~ParVectorTBase() {}

  /**
   * Return a code for the type of this parameter.
   */
  virtual string type() const;

  /**
   * Return a string describing the type of interface to be included
   * in the Doxygen documentation.
   */
  virtual string doxygenType() const;

  /**
   * Return a complete description of this parameter vector.
   */
  virtual string fullDescription(const InterfacedBase & ib) const;

  /**
   * Set the \a i'th element of a container of member variables of \a
   * ib to val. Uses a stringstream to read the \a val into a Type
   * object and then calls tset(InterfacedBase &, Type, int).
   */
  virtual void set(InterfacedBase & ib, string val, int i) const
   ;

  /**
   * Set the \a i'th element of a container of member variables of \a
   * ib to \a val.
   */
  virtual void tset(InterfacedBase & ib, Type val, int i)
    const = 0;

  /**
   * Insert a new object before the \a i'th element of a container of
   * member variables of \a ib and set it to \a val. Uses a
   * stringstream to read the \a val into a Type object and then calls
   * tinsert(InterfacedBase &, Type, int).
   */
  virtual void insert(InterfacedBase & ib, string val, int i) const
   ;

private:
  /// Implementation of set() for standard types.
  void setImpl(InterfacedBase & ib, string val, int i, StandardT)
    const;

  /// Implementation of set() for dimensioned types.
  void setImpl(InterfacedBase & ib, string val, int i, DimensionT) 
    const;

  /// Implementation of insert() for standard types.
  void insertImpl(InterfacedBase & ib, string val, int i, StandardT)
    const;

  /// Implementation of insert() for dimensioned types.
  void insertImpl(InterfacedBase & ib, string val, int i, DimensionT) 
    const;

public:
  
  /**
   * Insert a new object before the \a i'th element of a container of
   * member variables of \a ib and set it to \a val.
   */
  virtual void tinsert(InterfacedBase & ib, Type val, int i)
    const = 0;

  /**
   * Return the values of a container of member variables of \a ib in
   * a vector of strings. Calls the The tget(const InterfacedBase &)
   * and returns a vector of strings converted with ostringstreams.
   */
  virtual StringVector get(const InterfacedBase & ib) const
   ;

  /**
   * Return the values of a container of member variables of \a ib in a
   * vector of Type.
   */
  virtual TypeVector tget(const InterfacedBase & ib) const
    = 0;

  /**
   * Return the minimum value allowed for the \a i'th element of a
   * container of member variables of \a ib. Calls tminimum(const
   * InterfacedBase &, int) and converts the returned value with an
   * ostringstream.
   * 
   */
  virtual string minimum(const InterfacedBase & ib, int i) const
   ;

  /**
   * Return the minimum value allowed for the \a i'th element of a
   * container of member variables of \a ib.
   */
  virtual Type tminimum(const InterfacedBase & ib, int i) const
    = 0;

  /**
   * Return the maximum value allowed for the \a i'th element of a
   * container of member variables of \a ib. Calls tmaximum(const
   * InterfacedBase &, int) and converts the returned value with an
   * ostringstream.
   */
  virtual string maximum(const InterfacedBase & ib, int i) const
   ;

  /**
   * Return the maximum value allowed for the \a i'th element of a
   * container of member variables of \a ib.
   */
  virtual Type tmaximum(const InterfacedBase & ib, int i) const
    = 0;

  /**
   * Return the default value for the \a i'th element of a container
   * of member variables of \a ib. Calls tdef(const InterfacedBase &,
   * int) and converts the returned value with an ostringstream.
   */
  virtual string def(const InterfacedBase & ib, int i) const
   ;

  /**
   * Return the default value for the \a i'th element of a container of
   * member variables of \a ib.
   */
  virtual Type tdef(const InterfacedBase & ib, int i) const
    = 0;

  /**
   * Return the general default value for this parameter vector. Calls
   * tdef() and converts the returned value with an ostringstream.
   */
  virtual string def() const;

  /**
   * Return the general default value for this parameter vector.
   */
  virtual Type tdef() const = 0;

  /**
   * set the \a i'th element of a container of member variables of \a ib to
   * its default value.
   */
  virtual void setDef(InterfacedBase & ib, int i) const
   ;

  /**
   * Get the unit which an Type object is divided (multiplied) by when
   * written to (read from) a stream via a double. If unit() is zero,
   * the Type object is written/read directly.
   */
  Type unit() const { return theUnit; }

  /**
   * Set the unit which an Type object is divided (multiplied) by when
   * written to (read from) a stream via a double. If unit() is zero,
   * the Type object is written/read directly.
   */
  void unit(Type u) { theUnit = u; }

protected:

  /**
   * Write a numer to a stream with the unit specified with unit().
   */
  void putUnit(ostream & os, Type val) const { 
    putUnitImpl2(os, val, unit(), typename TypeTraits<Type>::DimType());
  }

private:

  /**
   * The unit which an Type object is divided (multiplied) by
   * when written to (read from) a stream via a double. If unit() is
   * zero, the Type object is written/read directly.
   */
  Type theUnit;

};

/**
 * The ParVector and its base classes ParVectorTBase and ParVectorBase
 * defines an interface to a class derived from the InterfacedBase,
 * through which vectors (or any other container) of simple member
 * variables can be manuipulated.  ParVector is templated on the type
 * of the member variable and the type of the InterfacedBase class,
 * and is derived from the InterfaceBase class via
 * ParVectorTBase (which is templated only on the class
 * of the member variable) and ParVectorBase.
 *
 * For each InterfacedBase class exactly one static ParVector object
 * should created for each member variable of container type which
 * should be interfaced. This object will automatically register
 * itself with the BaseRepository class.
 *
 * @see Parameter
 * @see InterfacedBase
 * @see InterfaceBase
 * 
 */
template <typename T, typename Type>
class ParVector: public ParVectorTBase<Type> {

public:

  /**
   * The declaration of member functions which can be used by this
   * ParVector interface for 'set' actions.
   */
  typedef void (T::*SetFn)(Type, int);

  /**
   * The declaration of member functions which can be used by this
   * ParVector interface for 'insert' actions.
   */
  typedef void (T::*InsFn)(Type, int);

  /**
   * The declaration of member functions which can be used by this
   * ParVector interface for 'erase' actions.
   */
  typedef void (T::*DelFn)(int);

  /**
   * A vector of objects of the template parameter Type.
   */
  typedef vector<Type> TypeVector;

  /**
   * The declaration of member functions which can be used by this
   * ParVector interface for 'get' actions.
   */
  typedef TypeVector (T::*GetFn)() const;

  /**
   * A vector of strings.
   */
  typedef vector<string> StringVector;

  /**
   * The declaration of member functions which can be used by this
   * ParVector interface for 'get' actions.
   */
  typedef StringVector (T::*StringGetFn)() const;

  /**
   * The declaration of member functions which can be used by this
   * ParVector interface for 'erase' actions.
   */
  typedef Type (T::*DefFn)(int) const;

  /**
   * Declaration of a direct pointer to the member variable in case it
   * is a vector.
   */
  typedef TypeVector T::* Member;

public:

  /**
   * Standard constructor.
   *
   * @param newName the name of the interface, may only contain
   * letters [a-zA-z0-9_].
   *
   * @param newDescription a brief description of the interface.
   *
   * @param newMember a pointer to a Member which is a TypeVector. May
   * be null, in which case the pointers to member functions must be
   * specified.
   *
   * @param newSize the size of the container or -1 if varying.
   *
   * @param newDef the default value of the corresponding parameters.
   *
   * @param newMin the minimum value of the corresponding parameters.
   *
   * @param newMax the maximum value of the corresponding parameters.
   *
   * @param depSafe set to true if calls to this interface for one
   * object does not influence other objects.
   *
   * @param readonly if this is set true the interface will not be
   * able to manipulate objects of the corresponding class, but will
   * still be able to access information.
   *
   * @param limits determines if the values of the parameters are
   * limited from above and below.
   *
   * @param newSetFn optional pointer to member function for the 'set'
   * action.
   *
   * @param newInsFn optional pointer to member function for the
   * 'insert' action.
   *
   * @param newDelFn optional pointer to member function for the
   * 'erase' action.
   *
   * @param newDefFn optional pointer to member function for the
   * 'default' action.
   *
   * @param newGetFn optional pointer to member function for the
   * 'get' action.
   *
   * @param newMinFn optional pointer to member function for the
   * 'minimum' action.
   *
   * @param newMaxFn optional pointer to member function for the
   * 'maximum' action.
   *
   * @param newStringGetFn optional pointer to member function for the
   * 'get' action.
   */
  ParVector(string newName, string newDescription,
		   Member newMember, int newSize, Type newDef, Type newMin,
		   Type newMax, bool depSafe = false, bool readonly = false,
		   bool limits = true, SetFn newSetFn = 0,
		   InsFn newInsFn = 0, DelFn newDelFn = 0, GetFn newGetFn = 0,
		   DefFn newDefFn = 0, DefFn newMinFn = 0, DefFn newMaxFn = 0,
		   StringGetFn newStringGetFn = 0)
    : ParVectorTBase<Type>(newName, newDescription, ClassTraits<T>::className(),
			   typeid(T), Type(), newSize, depSafe, readonly,
			   limits),
      theMember(newMember), theDef(newDef), theMin(newMin), theMax(newMax),
      theSetFn(newSetFn), theInsFn(newInsFn), theDelFn(newDelFn),
      theGetFn(newGetFn), theDefFn(newDefFn), theMinFn(newMinFn),
      theMaxFn(newMaxFn), theStringGetFn(newStringGetFn) {}


  /**
   * Standard constructor.
   *
   * @param newName the name of the interface, may only contain
   * letters [a-zA-z0-9_].
   *
   * @param newDescription a brief description of the interface.
   *
   * @param newMember a pointer to a Member which is a TypeVector. May
   * be null, in which case the pointers to member functions must be
   * specified.
   *
   * @param newUnit the unit assumed when a number is read or written
   * to a stream.
   *
   * @param newSize the size of the container or -1 if varying.
   *
   * @param newDef the default value of the corresponding parameters.
   *
   * @param newGetFn optional pointer to member function for the
   * 'get' action.
   *
   * @param newMin the minimum value of the corresponding parameters.
   *
   * @param newMax the maximum value of the corresponding parameters.
   *
   * @param depSafe set to true if calls to this interface for one
   * object does not influence other objects.
   *
   * @param readonly if this is set true the interface will not be
   * able to manipulate objects of the corresponding class, but will
   * still be able to access information.
   *
   * @param limits determines if the values of the parameters are
   * limited from above and below.
   *
   * @param newSetFn optional pointer to member function for the 'set'
   * action.
   *
   * @param newInsFn optional pointer to member function for the
   * 'insert' action.
   *
   * @param newDelFn optional pointer to member function for the
   * 'erase' action.
   *
   * @param newDefFn optional pointer to member function for the
   * 'default' action.
   *
   * @param newMinFn optional pointer to member function for the
   * 'minimum' action.
   *
   * @param newMaxFn optional pointer to member function for the
   * 'maximum' action.
   *
   * @param newStringGetFn optional pointer to member function for the
   * 'get' action.
   */
  ParVector(string newName, string newDescription, Member newMember,
		   Type newUnit, int newSize, Type newDef, Type newMin,
		   Type newMax, bool depSafe = false, bool readonly = false,
		   bool limits = true, SetFn newSetFn = 0,
		   InsFn newInsFn = 0, DelFn newDelFn = 0, GetFn newGetFn = 0,
		   DefFn newDefFn = 0, DefFn newMinFn = 0, DefFn newMaxFn = 0,
		   StringGetFn newStringGetFn = 0)
    : ParVectorTBase<Type>(newName, newDescription, ClassTraits<T>::className(),
			   typeid(T), newUnit, newSize, depSafe, readonly,
			   limits),
      theMember(newMember), theDef(newDef), theMin(newMin), theMax(newMax),
      theSetFn(newSetFn), theInsFn(newInsFn), theDelFn(newDelFn),
      theGetFn(newGetFn), theDefFn(newDefFn), theMinFn(newMinFn),
      theMaxFn(newMaxFn), theStringGetFn(newStringGetFn) {}

  /**
   * Standard constructor.
   *
   * @param newName the name of the interface, may only contain
   * letters [a-zA-z0-9_].
   *
   * @param newDescription a brief description of the interface.
   *
   * @param newMember a pointer to a Member which is a TypeVector. May
   * be null, in which case the pointers to member functions must be
   * specified.
   *
   * @param newSize the size of the container or -1 if varying.
   *
   * @param newDef the default value of the corresponding parameters.
   *
   * @param newMin the minimum value of the corresponding parameters.
   *
   * @param newMax the maximum value of the corresponding parameters.
   *
   * @param depSafe set to true if calls to this interface for one
   * object does not influence other objects.
   *
   * @param readonly if this is set true the interface will not be
   * able to manipulate objects of the corresponding class, but will
   * still be able to access information.
   *
   * @param limits determines if the values of the parameters are
   * limited from above and/or below. The possible values are given by
   * Interface::Limits.
   *
   * @param newSetFn optional pointer to member function for the 'set'
   * action.
   *
   * @param newInsFn optional pointer to member function for the
   * 'insert' action.
   *
   * @param newDelFn optional pointer to member function for the
   * 'erase' action.
   *
   * @param newDefFn optional pointer to member function for the
   * 'default' action.
   *
   * @param newGetFn optional pointer to member function for the
   * 'get' action.
   *
   * @param newMinFn optional pointer to member function for the
   * 'minimum' action.
   *
   * @param newMaxFn optional pointer to member function for the
   * 'maximum' action.
   *
   * @param newStringGetFn optional pointer to member function for the
   * 'get' action.
   */
  ParVector(string newName, string newDescription,
		   Member newMember, int newSize, Type newDef, Type newMin,
		   Type newMax, bool depSafe = false, bool readonly = false,
		   int limits = Interface::limited, SetFn newSetFn = 0,
		   InsFn newInsFn = 0, DelFn newDelFn = 0, GetFn newGetFn = 0,
		   DefFn newDefFn = 0, DefFn newMinFn = 0, DefFn newMaxFn = 0,
		   StringGetFn newStringGetFn = 0)
    : ParVectorTBase<Type>(newName, newDescription, ClassTraits<T>::className(),
			   typeid(T), Type(), newSize, depSafe, readonly,
			   limits),
      theMember(newMember), theDef(newDef), theMin(newMin), theMax(newMax),
      theSetFn(newSetFn), theInsFn(newInsFn), theDelFn(newDelFn),
      theGetFn(newGetFn), theDefFn(newDefFn), theMinFn(newMinFn),
      theMaxFn(newMaxFn), theStringGetFn(newStringGetFn) {}

  /**
   * Standard constructor.
   *
   * @param newName the name of the interface, may only contain
   * letters [a-zA-z0-9_].
   *
   * @param newDescription a brief description of the interface.
   *
   * @param newMember a pointer to a Member which is a TypeVector. May
   * be null, in which case the pointers to member functions must be
   * specified.
   *
   * @param newUnit the unit assumed when a number is read or written
   * to a stream.
   *
   * @param newSize the size of the container or -1 if varying.
   *
   * @param newDef the default value of the corresponding parameters.
   *
   * @param newGetFn optional pointer to member function for the
   * 'get' action.
   *
   * @param newMin the minimum value of the corresponding parameters.
   *
   * @param newMax the maximum value of the corresponding parameters.
   *
   * @param depSafe set to true if calls to this interface for one
   * object does not influence other objects.
   *
   * @param readonly if this is set true the interface will not be
   * able to manipulate objects of the corresponding class, but will
   * still be able to access information.
   *
   * @param limits determines if the values of the parameters are
   * limited from above and/or below. The possible values are given by
   * Interface::Limits.
   *
   * @param newSetFn optional pointer to member function for the 'set'
   * action.
   *
   * @param newInsFn optional pointer to member function for the
   * 'insert' action.
   *
   * @param newDelFn optional pointer to member function for the
   * 'erase' action.
   *
   * @param newDefFn optional pointer to member function for the
   * 'default' action.
   *
   * @param newMinFn optional pointer to member function for the
   * 'minimum' action.
   *
   * @param newMaxFn optional pointer to member function for the
   * 'maximum' action.
   *
   * @param newStringGetFn optional pointer to member function for the
   * 'get' action.
   */
  ParVector(string newName, string newDescription, Member newMember,
		   Type newUnit, int newSize, Type newDef, Type newMin,
		   Type newMax, bool depSafe = false, bool readonly = false,
		   int limits = Interface::limited, SetFn newSetFn = 0,
		   InsFn newInsFn = 0, DelFn newDelFn = 0, GetFn newGetFn = 0,
		   DefFn newDefFn = 0, DefFn newMinFn = 0, DefFn newMaxFn = 0,
		   StringGetFn newStringGetFn = 0)
    : ParVectorTBase<Type>(newName, newDescription, ClassTraits<T>::className(),
			   typeid(T), newUnit, newSize, depSafe, readonly,
			   limits),
      theMember(newMember), theDef(newDef), theMin(newMin), theMax(newMax),
      theSetFn(newSetFn), theInsFn(newInsFn), theDelFn(newDelFn),
      theGetFn(newGetFn), theDefFn(newDefFn), theMinFn(newMinFn),
      theMaxFn(newMaxFn), theStringGetFn(newStringGetFn) {}

  /**
   * Set the \a i'th element of a container of member variables of \a
   * ib to \a val.
   */
  virtual void tset(InterfacedBase & ib, Type val, int i) const
   ;

  /**
   * Insert a new object before the \a i'th element of a container of
   * member variables of \a ib and set it to \a val.
   */
  virtual void tinsert(InterfacedBase & ib, Type val, int i) const
   ;

  /**
   * Remove the \a i'th element of a container of
   * member variables of \a ib.
   */
  virtual void erase(InterfacedBase & ib, int i)
    const;

  /**
   * Clear the container of pointers of \a ib.
   */
  virtual void clear(InterfacedBase & ib)
    const;

  /**
   * Return the values of a container of member variables of \a ib in
   * a vector of strings. Calls the The tget(const InterfacedBase &)
   * and returns a vector of strings converted with ostringstreams.
   */
  virtual StringVector get(const InterfacedBase & ib) const
   ;

  /**
   * Return the values of a container of member variables of \a ib in a
   * vector of Type.
   */
  virtual TypeVector tget(const InterfacedBase & ib) const
   ;

  /**
   * Return the minimum value allowed for the \a i'th element of a
   * container of member variables of \a ib.
   */
  virtual Type tminimum(const InterfacedBase & ib, int i) const
   ;

  /**
   * Return the maximum value allowed for the \a i'th element of a
   * container of member variables of \a ib.
   */
  virtual Type tmaximum(const InterfacedBase & ib, int i) const
   ;

  /**
   * Return the default value for the \a i'th element of a container of
   * member variables of \a ib.
   */
  virtual Type tdef(const InterfacedBase &, int) const
   ;

  /**
   * Return the general default value for this parameter vector.
   */
  virtual Type tdef() const;

  /**
   * Give a pointer to a member function to be used by tset().
   */
  void setSetFunction(SetFn sf) { theSetFn = sf; }

  /**
   * Give a pointer to a member function to be used by tinsert().
   */
  void setInsertFunction(InsFn ifn) { theInsFn = ifn; }

  /**
   * Give a pointer to a member function to be used by tget().
   */
  void setGetFunction(GetFn gf) { theGetFn = gf; }

  /**
   * Give a pointer to a member function to be used by terase().
   */
  void setEraseFunction(DelFn df) { theDelFn = df; }

  /**
   * Give a pointer to a member function to be used by tdef().
   */
  void setDefaultFunction(GetFn df) { theDefFn = df; }

  /**
   * Give a pointer to a member function to be used by tminimum().
   */
  void setMinFunction(GetFn mf) { theMinFn = mf; }

  /**
   * Give a pointer to a member function to be used by tmaximum().
   */
  void setMaxFunction(GetFn mf) { theMaxFn = mf; }

  /**
   * Give a pointer to a member function to be used by get().
   */
  void setStringGetFunction(StringGetFn gf) { theStringGetFn = gf; }

  /**
   * Print a description to be included in the Doxygen documentation
   * to the given \a stream.
   */
  virtual void doxygenDescription(ostream & stream) const;

private:

  /**
   * The pointer to the member variable.
   */
  Member theMember;

  /**
   * Default value to be used if no corresponding member function
   * pointer is given.
   */
  Type theDef;

  /**
   * Minimum value to be used if no corresponding member function
   * pointer is given.
   */
  Type theMin;

  /**
   * Maximum value to be used if no
   * corresponding member function pointer is given.
   */
  Type theMax;

  /**
   * A pointer to a member function to be used by tset().
   */
  SetFn theSetFn;

  /**
   * A pointer to a member function to be used by tinsert().
   */
  InsFn theInsFn;

  /**
   * A pointer to a member function to be used by terase().
   */
  DelFn theDelFn;

  /**
   * A pointer to a member function to be used by tget().
   */
  GetFn theGetFn;

  /**
   * Pointer to member function to be used by tdef().
   */
  DefFn theDefFn;

  /**
   * Pointer to member function to be used by tminimum().
   */
  DefFn theMinFn;

  /**
   * Pointer to member function to be used by tmaximum().
   */
  DefFn theMaxFn;

  /**
   * A pointer to a member function to be used by set().
   */
  StringGetFn theStringGetFn;

};

}

#ifndef ThePEG_TEMPLATES_IN_CC_FILE
#include "ParVector.tcc"
#endif

#endif /* ThePEG_ParVector_H */
