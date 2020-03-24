// -*- C++ -*-
//
// RefVector.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_RefVector_H
#define ThePEG_RefVector_H
// This is the declaration of the RefVector and RefVectorBase classes.

#include "ThePEG/Config/ThePEG.h"
#include "InterfaceBase.h"
#include "RefVector.xh"
#include "RefVector.fh"

namespace ThePEG {

/**
 * The RefVector and its base class RefVectorBase defines an interface
 * to a class derived from the InterfacedBase, through which vectors
 * (or any other container) of pointers to other InterfacedBase
 * objects may be manipulated.  RefVector is templated on the type of
 * the class and the class of the objects pointed to, and is derived
 * from the InterfaceBase class via RefVectorBase and
 * RefInterfaceBase.
 *
 * For each InterfacedBase class exactly one static RefVector object
 * should created for each member variable of container type which
 * should be interfaced. This object will automatically register
 * itself with the BaseRepository class.
 * 
 * @see InterfacedBase
 * @see InterfaceBase
 * 
 */
class RefVectorBase: public RefInterfaceBase {

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
   * @param newSize the size of the container or -1 if varying.
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
  RefVectorBase(string newName, string newDescription,
		string newClassName,
		const type_info & newTypeInfo,
		string newRefClassName,
		const type_info & newRefTypeInfo, 
		int newSize, bool depSafe,
		bool readonly, bool norebind, bool nullable, bool defnull);

  /**
   * The general interface method overriding the one in
   * InterfaceBase. For this class, \a action can be any of "set",
   * "get", "erase" and "insert" and \a argument should be a something
   * which can be read into an integer while the rest of \a argument
   * should correspond to the name of an InterfacedBase object in the
   * BaseRepository.
   */
  virtual string exec(InterfacedBase & ib, string action,
		      string arguments) const;

  /**
   * Return a complete description of this reference vector.
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
   * Set the \a i'th element of a container of pointers of \a ib
   * to \a ip.
   */
  virtual void set(InterfacedBase & ib, IBPtr ip, int i, bool chk = true)
    const = 0;

  /**
   * Insert a new pointer before the \a i'th element of a container of
   * pointers of \a ib and set it to \a ip.
   */
  virtual void insert(InterfacedBase & ib, IBPtr ip, int i, bool chk = true)
    const = 0;

  /**
   * Remove the \a i'th element of a container of pointers of \a ib.
   */
  virtual void erase(InterfacedBase & ib, int i)
    const = 0;

  /**
   * Clear the container of pointers of \a ib.
   */
  virtual void clear(InterfacedBase & ib)
    const = 0;

  /**
   * Return a vector of pointers corresponding to the container of
   * pointers of \a ib.
   */
  virtual IVector get(const InterfacedBase & ib) const
    = 0;

  /**
   * Check if set(ib, ip, i) will be successfull but do not do
   * anything.
   */
  virtual bool check(const InterfacedBase & ib, cIBPtr ip, int i) const
    = 0;

  /**
   * In the object \a ib, replace all pointers to objects in this
   * interface with the translated ones provided by \a trans. If a
   * pointer is null, and defaultIfNull() is true, replace it with
   * the first alowed object found in \a defs.
   */
  virtual void rebind(InterfacedBase & ib, const TranslationMap & trans,
		      const IVector & defs) const;

  /**
   * Return pointers to other objects in \a ib.
   */
  virtual IVector getReferences(const InterfacedBase & ib) const;

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
   * The size of the container being interfaced.
   */
  int theSize;

};


/**
 * The RefVector and its base class RefVectorBase defines an interface
 * to a class derived from the InterfacedBase, through which vectors
 * (or any other container) of pointers to other InterfacedBase
 * objects may be manipulated.  RefVector is templated on the type of
 * the class and the class of the objects pointed to, and is derived
 * from the InterfaceBase class via RefVectorBase and
 * RefInterfaceBase.
 *
 * For each InterfacedBase class exactly one static RefVector object
 * should created for each member variable of container type which
 * should be interfaced. This object will automatically register
 * itself with the BaseRepository class.
 * 
 * @see InterfacedBase
 * @see InterfaceBase
 * 
 */
template <class T, class R>
class RefVector: public RefVectorBase {

public:

  /** A pointer to the class of objects referred to. */
  typedef typename Ptr<R>::pointer RefPtr;
  /** A const pointer to the class of objects referred to. */
  typedef typename Ptr<R>::const_pointer cRefPtr;
  /** A pointer to a menberfunction to be used for the 'set' action. */
  typedef void (T::*SetFn)(RefPtr, int);
  /** A pointer to a menberfunction to be used for the 'insert' action. */
  typedef void (T::*InsFn)(RefPtr, int);
  /** A pointer to a menberfunction to be used for the 'erase' action. */
  typedef void (T::*DelFn)(int);
  /** A pointer to a menberfunction to be used for the 'check' action. */
  typedef bool (T::*CheckFn)(cRefPtr, int) const;
  /** A pointer to a menberfunction to be used for the 'get' action. */
  typedef vector<RefPtr> (T::*GetFn)() const;
  /**
   * Declaration of a direct pointer to the member variable in case it
   * is a vector.
   */
  typedef vector<RefPtr> T::* Member;

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
   * @param newInsFn optional pointer to member function for the
   * 'insert' action.
   *
   * @param newDelFn optional pointer to member function for the
   * 'erase' action.
   *
   * @param newGetFn optional pointer to member function for the
   * 'get' action.
   *
   * @param newCheckFn optional pointer to member function for the
   * 'check' action.
   */
  RefVector(string newName, string newDescription,
	    Member newMember, int newSize, bool depSafe = false,
	    bool readonly = false, bool rebind = true, bool nullable = true,
	    SetFn newSetFn = 0, InsFn newInsFn = 0, DelFn newDelFn = 0,
	    GetFn newGetFn = 0, CheckFn newCheckFn = 0);

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
   * @param newInsFn optional pointer to member function for the
   * 'insert' action.
   *
   * @param newDelFn optional pointer to member function for the
   * 'erase' action.
   *
   * @param newGetFn optional pointer to member function for the
   * 'get' action.
   *
   * @param newCheckFn optional pointer to member function for the
   * 'check' action.
   */
  RefVector(string newName, string newDescription,
	    Member newMember, int newSize, bool depSafe,
	    bool readonly, bool rebind, bool nullable, bool defnull,
	    SetFn newSetFn = 0, InsFn newInsFn = 0, DelFn newDelFn = 0,
	    GetFn newGetFn = 0, CheckFn newCheckFn = 0);

  /**
   * Set the \a i'th element of a container of pointers of \a ib
   * to \a ip.
   */
  virtual void set(InterfacedBase & ib, IBPtr ip, int i, bool chk = true)
    const;

  /**
   * Insert a new pointer before the \a i'th element of a container of
   * pointers of \a ib and set it to \a ip.
   */
  virtual void insert(InterfacedBase & ib, IBPtr ip, int i, bool chk = true)
    const;

  /**
   * Remove the \a i'th element of a container of pointers of \a ib.
   */
  virtual void erase(InterfacedBase & ib, int i)
    const;

  /**
   * Clear the container of pointers of \a ib.
   */
  virtual void clear(InterfacedBase & ib)
    const;

  /**
   * Return a vector of pointers corresponding to the container of
   * pointers of \a ib.
   */
  virtual IVector get(const InterfacedBase & ib) const
   ;

  /**
   * Check if set(ib, ip, i) will be successfull but do not do
   * anything.
   */
  virtual bool check(const InterfacedBase & ib, cIBPtr, int i) const
   ;

  /**
   * Give a pointer to a member function to be used by 'set()'.
   */
  void setSetFunction(SetFn sf) { theSetFn = sf; }

  /**
   * Give a pointer to a member function to be used by 'insert()'.
   */
  void setInsertFunction(InsFn ifn) { theInsFn = ifn; }

  /**
   * Give a pointer to a member function to be used by 'get()'.
   */
  void setGetFunction(GetFn gf) { theGetFn = gf; }

  /**
   * Give a pointer to a member function to be used by 'erase()'.
   */
  void setEraseFunction(DelFn df) { theDelFn = df; }

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
   * Give a pointer to a member function to be used by 'insert()'.
   */
  InsFn theInsFn;

  /**
   * Give a pointer to a member function to be used by 'erase()'.
   */
  DelFn theDelFn;

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

#include "RefVector.tcc"

#endif /* ThePEG_RefVector_H */
