// -*- C++ -*-
//
// Reference.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_Reference_H
#define ThePEG_Reference_H
// This is the declaration of the Reference and ReferenceBase classes.

#include "ThePEG/Config/ThePEG.h"
#include "InterfaceBase.h"
#include "Reference.xh"
#include "Reference.fh"

namespace ThePEG {

/**
 * The Reference class and its base class ReferenceBase defines an
 * interface to a class derived from the InterfacedBase, through which
 * pointers to other InterfacedBase objects may be manipulated.
 * Reference is templated on the type of the class and the class of
 * the objects pointed to, and is derived from the InterfaceBase class
 * via ReferenceBase and RefInterfaceBase.
 *
 * For each InterfacedBase class exactly one static Reference object
 * should created for each member variable which should be
 * interfaced. This object will automatically register itself with the
 * BaseRepository class.
 *
 * @see InterfacedBase
 * @see InterfaceBase
 * 
 */
class ReferenceBase: public RefInterfaceBase {

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
   * @param newRefClassName the name of the class pointed to.
   *
   * @param newRefTypeInfo the type_info object of the class pointed
   * to.
   *
   * @param depSafe set to true if calls to this interface for one
   * object does not influence other objects.
   *
   * @param readonly if this is set true the interface will not be
   * able to manipulate objects of the corresponding class, but will
   * still be able to access information.
   *
   * @param norebind if set to true, this interface is not responsible
   * for the rebinding of corresponding objects.
   *
   * @param nullable if set to true this corresponding references may
   * be null.
   *
   * @param defnull if set to true and a corresponding reference is
   * null it may be given a a default value in the initialization of
   * an EventGenerator.
   */
  ReferenceBase(string newName, string newDescription,
		string newClassName,
		const type_info & newTypeInfo, 
		string newRefClassName,
		const type_info & newRefTypeInfo, bool depSafe,
		bool readonly, bool norebind, bool nullable, bool defnull);

  /**
   * The general interface method overriding the one in
   * InterfaceBase. For this class, \a action can be any of "set" and
   * "get" and \a argument should correspond to the name of an
   * InterfacedBase object in the BaseRepository.
   */
  virtual string exec(InterfacedBase & ib, string action,
		      string arguments) const;

  /**
   * Return a complete description of this reference.
   */
  virtual string fullDescription(const InterfacedBase & ib) const;

  /**
   * Return a code for the type of this reference.
   */
  virtual string type() const;

  /**
   * Return a string describing the type of interface to be included
   * in the Doxygen documentation.
   */
  virtual string doxygenType() const;

  /**
   * Set the pointer of \a ib to \a ip.
   */
  virtual void set(InterfacedBase & ib, IBPtr ip, bool chk = true)
    const = 0;

  /**
   * Return the pointer of \a ib.
   */
  virtual IBPtr get(const InterfacedBase & ib)
    const = 0;

  /**
   * Check if set(ib, ip) will be successfull but do not do
   * anything.
   */
  virtual bool check(const InterfacedBase & ib, cIBPtr ip) const
    = 0;

  /**
   * In the object \a ib, replace the pointer in this interface with one
   * of the translated ones provided by trans. If the pointer is null,
   * and defaultIfNull() is true, replace it with the first allowed
   * object found in \a defs.
   */
  virtual void rebind(InterfacedBase & ib, const TranslationMap & trans,
		      const IVector & defs) const;

  /**
   * Return the pointer to another object in \a ib (in a vector).
   */
  virtual IVector getReferences(const InterfacedBase & ib) const;

};


/**
 * The Reference and its base class ReferenceBase defines an interface
 * to a class derived from the InterfacedBase, through which pointers
 * to other InterfacedBase objects may be manipulated.  Reference is
 * templated on the type of the class and the class of the objects
 * pointed to, and is derived from the InterfaceBase class via
 * ReferenceBase and RefInterfaceBase.
 *
 * For each InterfacedBase class exactly one static Reference object
 * should created for each member variable which should be
 * interfaced. This object will automatically register itself with the
 * BaseRepository class.
 *
 * @see InterfacedBase
 * @see InterfaceBase
 * 
 */
template <class T, class R>
class Reference: public ReferenceBase {

public:

  /** A pointer to the class of objects referred to. */
  typedef typename Ptr<R>::pointer RefPtr;
  /** A const pointer to the class of objects referred to. */
  typedef typename Ptr<R>::const_pointer cRefPtr;
  /** A pointer to a menberfunction to be used for the 'set' action. */
  typedef void (T::*SetFn)(RefPtr);
  /** A pointer to a menberfunction to be used for the 'check' action. */
  typedef bool (T::*CheckFn)(cRefPtr) const;
  /** A pointer to a menberfunction to be used for the 'get' action. */
  typedef RefPtr (T::*GetFn)() const;
  /** Declaration of a direct pointer to the member variable. */
  typedef RefPtr T::* Member;

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
   * @param depSafe set to true if calls to this interface for one
   * object does not influence other objects.
   *
   * @param readonly if this is set true the interface will not be
   * able to manipulate objects of the corresponding class, but will
   * still be able to access information.
   *
   * @param rebind if set to true, this interface is responsible
   * for the rebinding of corresponding objects.
   *
   * @param nullable if set to true this corresponding references may
   * be null.
   *
   * @param newSetFn optional pointer to member function for the 'set'
   * action.
   *
   * @param newGetFn optional pointer to member function for the
   * 'get' action.
   *
   * @param newCheckFn optional pointer to member function for the
   * 'check' action.
   */
  Reference(string newName, string newDescription,
	    Member newMember, bool depSafe = false,
	    bool readonly = false, bool rebind = true, bool nullable = true,
	    SetFn newSetFn = 0, GetFn newGetFn = 0,
	    CheckFn newCheckFn = 0)
    : ReferenceBase(newName, newDescription,
		    ClassTraits<T>::className(), typeid(T),
		    ClassTraits<R>::className(), typeid(R),
		    depSafe, readonly, !rebind, nullable, false),
      theMember(newMember), theSetFn(newSetFn), theGetFn(newGetFn),
      theCheckFn(newCheckFn) {}

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
   * @param depSafe set to true if calls to this interface for one
   * object does not influence other objects.
   *
   * @param readonly if this is set true the interface will not be
   * able to manipulate objects of the corresponding class, but will
   * still be able to access information.
   *
   * @param rebind if set to true, this interface is responsible
   * for the rebinding of corresponding objects.
   *
   * @param nullable if set to true this corresponding references may
   * be null.
   *
   * @param defnull if set to true and a corresponding reference is
   * null it may be given a a default value in the initialization of
   * an EventGenerator.
   *
   * @param newSetFn optional pointer to member function for the 'set'
   * action.
   *
   * @param newGetFn optional pointer to member function for the
   * 'get' action.
   *
   * @param newCheckFn optional pointer to member function for the
   * 'check' action.
   */
  Reference(string newName, string newDescription,
	    Member newMember, bool depSafe, bool readonly, bool rebind,
	    bool nullable, bool defnull, SetFn newSetFn = 0, GetFn newGetFn = 0,
	    CheckFn newCheckFn = 0)
    : ReferenceBase(newName, newDescription,
		    ClassTraits<T>::className(), typeid(T),
		    ClassTraits<R>::className(), typeid(R),
		    depSafe, readonly, !rebind, nullable, defnull),
      theMember(newMember), theSetFn(newSetFn), theGetFn(newGetFn),
      theCheckFn(newCheckFn) {}


  /**
   * Set the pointer of \a ib to \a ip.
   */
  virtual void set(InterfacedBase & ib, IBPtr ip, bool chk = true) const
   ;

  /**
   * Return the pointer of \a ib.
   */
  virtual IBPtr get(const InterfacedBase & ib) const
   ;

  /**
   * Check if set(ib, ip) will be successfull but do not do
   * anything.
   */
  virtual bool check(const InterfacedBase & ib, cIBPtr newRef) const
   ;

  /**
   * Give a pointer to a member function to be used by 'set()'.
   */
  void setSetFunction(SetFn sf) { theSetFn = sf; }

  /**
   * Give a pointer to a member function to be used by 'get()'.
   */
  void setGetFunction(GetFn gf) { theGetFn = gf; }

  /**
   * Give a pointer to a member function to be used by 'check()'.
   */
  void setCheckFunction(CheckFn cf) { theCheckFn = cf; }

private:

  /**
   * The pointer to the member variable.
   */
  Member theMember;

  /**
   * A pointer to a member function to be used by 'set()'.
   */
  SetFn theSetFn;

  /**
   * Give a pointer to a member function to be used by 'get()'.
   */
  GetFn theGetFn;

  /**
   * Give a pointer to a member function to be used by 'check()'.
   */
  CheckFn theCheckFn;


};

}

#include "Reference.tcc"

#endif /* ThePEG_Reference_H */
