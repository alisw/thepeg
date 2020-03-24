// -*- C++ -*-
//
// InterfaceBase.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_InterfaceBase_H
#define ThePEG_InterfaceBase_H
// This is the declaration of the InterfaceBase and RefInterfaceBase classes.

#include "ThePEG/Config/ThePEG.h"
#include "InterfaceBase.fh"
#include "InterfaceBase.xh"
#include "ThePEG/Utilities/Named.h"
#include "ThePEG/Utilities/ClassTraits.h"
#include "Interface.h"

namespace ThePEG {

/**
 * The InterfaceBase class defines a generic interface to any class
 * derived from the InterfacedBase class. Using the pure virtual
 * exec() function, it is possible to manipulate any InterfacedBase
 * object. InterfaceBase is an abstract base class for derived classes
 * such as Command, Parameter and Reference.
 *
 * InterfaceBase objects are managed by the BaseRepository.
 *
 * InterfaceBase is derived from the Named to manage the name of the
 * interface.
 *
 * From the Repository it is possible to generate a file with doxygen
 * comments which can be included in the documentation describing the
 * InterfaceBase objects defined for a class. For each class,
 * <code>ClassName</code>, there will be produced a file called
 * <code>ClassNameInterfaces.html</code> which can be referred to with
 * a standard html <code>href</code> anchor. Also a specific
 * interface, <code>InterfaceName</code> can be referred to with
 * <code>ClassNameInterfaces.html#InterfaceName</code>. The file can
 * also be referred to with the doxygen <code>\\ref</code>
 * command. Inside the description of an interface, other interfaces
 * in the same class can be tagged with
 * \<interface\>InterfaceName\</interface\> or, if the interface
 * belongs to another class,
 * \<interface\>ClassName::InterfaceName\</interface\>. By running the
 * script in <code>ThePEG/Doc/fixinterfaces.pl</code> these tags will
 * be converted to proper <code>href</code> anchors.
 *
 * @see InterfacedBase
 * @see Command
 * @see Parameter
 * @see Reference
 * @see BaseRepository
 * @see Named
 * 
 */
class InterfaceBase: public Named {

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
   *
   * @param readonly if this is set true the interface will not be
   * able to manipulate objects of the corresponding class, but will
   * still be able to access information.
   */
  InterfaceBase(string newName, string newDescription,
		string newClassName,
		const type_info & newTypeInfo, bool depSafe,
		bool readonly);

  /**
   * The destructor.
   */
  virtual ~InterfaceBase() {}

  /**
   * Create a tag for this interface using its name and optional
   * poisitional argument.
   */
  string tag(int pos = -1) const;

  /**
   * Manipulate an object of the corresponding class. Execute the \a
   * action command with the given \a arguments for the given object
   * \a ib.
   */
  virtual string
  exec(InterfacedBase & ib, string action, string arguments) const
    = 0;

  /**
   * Return a code for the type of this interface.
   */
  virtual string type() const = 0;

  /**
   * Returns true if the setting for this interface has been changed
   * from its default value.
   */
  virtual bool notDefault(InterfacedBase &) const;

  /**
   * Returns the map of objectDefaults of the given object.
   */
  map<string,string> & objectDefaults(InterfacedBase &) const;

  /**
   * Rebind all references in ib according to the translation
   * map. Only used by derived classed interfacing references.
   */
  virtual void rebind(InterfacedBase &,
		      const TranslationMap &,
		      const IVector & = IVector()) const {}

  /**
   * For derived classes interfacing references between Interfaced
   * objects, return the references for this interface.
   */
  virtual IVector getReferences(const InterfacedBase &) const {
    return IVector();
  }

  /**
   * Return the description of this interface.
   */
  string description() const { return theDescription; }

  /**
   * Return a complete description of this interface.
   */
  virtual string fullDescription(const InterfacedBase & ib) const;

  /**
   * Print a description to be included in the Doxygen documentation
   * to the given \a stream.
   */
  virtual void doxygenDescription(ostream & stream) const;

  /**
   * Return a string describing the type of interface to be included
   * in the Doxygen documentation.
   */
  virtual string doxygenType() const = 0;

  /**
   * Return the class name for the class this interface is defined
   * for.
   */
  string className() const { return theClassName; }    

  /**
   * Get the flag saying whether changing an object with this
   * interface may change the state of a dependent object .
   */
  bool dependencySafe() const { return isDependencySafe; } 

  /**
   * Set the flag saying whether changing an object with this
   * interface may change the state of a dependent object .
   */
  void setDependencySafe() { isDependencySafe = true; } 

  /**
   * Set the flag saying whether changing an object with this
   * interface may change the state of a dependent object .
   */
  void setDependencySensitive() { isDependencySafe = false; } 

  /**
   * Get the flag saying whether this interface is allowed to change
   * an object.
   */
  bool readOnly() const { return isReadOnly && (!NoReadOnly); } 

  /**
   * Set the flag saying that this interface is allowed to
   * change an object.
   */
  void setReadOnly() { isReadOnly = true; } 

  /**
   * Unset the flag saying that this interface is allowed to change an
   * object.
   */
  void setReadWrite() { isReadOnly = false; } 

  /**
   * Return true if this interface is anonyous, ie. invisible for the
   * user interface.
   */
  bool anonymous() const { return description().empty(); } 

  /**
   * Get the rank for this interface. Used for sorting by user
   * interface.
   */
  double rank() const { return theRank; } 

  /**
   * Set the rank for this interface. Used for sorting by user
   * interface.
   */
  void rank(double r) { theRank = r; } 

  /**
   * Indicate that this interface has a default value.
   */
  void setHasDefault(bool b) {
    hasDefault = b;
  }

  /**
   * If set to true, all read-only interfaces can be changed.
   */
  static bool NoReadOnly;

private:

  /**
   * The description of this interface.
   */
  string theDescription;

  /**
   * The class name and for the class this interface is defined for.
   */
  string theClassName;

  /**
   * A rank assigned to this interface. Used for sorting by user
   * interface.
   */
  double theRank;

protected:

  /**
   * A flag indicating whether this interface has a default setting.
   */
  bool hasDefault;

  /**
   * The flag saying whether changing an object with this interface
   * may change the state of a dependent object .
   */
  mutable bool isDependencySafe;

  /**
   * The flag saying whether this interface is allowed to change an
   * object.
   */
  mutable bool isReadOnly;


};


/**
 * RefInterfaceBase is an abstract base class inheriting from
 * InterfaceBase used for subclasses dealing with interfaces to do
 * with references in one Interfaced object to another.
 */
class RefInterfaceBase: public InterfaceBase {

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
   * @param newRefClassName the name of the class referred to by the
   * corresponding class..
   *
   * @param newRefTypeInfo the type_info object of the class referred
   * to by the corresponding class.
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
   RefInterfaceBase(string newName, string newDescription,
			  string newClassName, const type_info & newTypeInfo,
			  string newRefClassName,
			  const type_info & newRefTypeInfo,
			  bool depSafe, bool readonly,
			  bool norebind, bool nullable, bool defnull);

  /**
   * Return the class name of the class referred to by
   * this interface.
   */
  string refClassName() const { return theRefClassName; } 

  /**
   * Return the type_info object of the class referred to by this
   * interface.
   */
  const type_info & refTypeInfo() const { return theRefTypeInfo; } 

  /**
   * Get the flag saying whether the interface is responsible for
   * rebinding of the corresponding refenerces.
   */
  bool noRebind() const { return dontRebind; } 

  /**
   * Set the flag saying that the interface is not responsible for
   * rebinding refenerces.
   */
  void setNoRebind() { dontRebind = true; } 

  /**
   * Set the flag saying that the interface is responsible for
   * rebinding refenerces.
   */
  void setRebind() { dontRebind = false; } 

  /**
   * Get the flag saying whether the interface is allowed to set the
   * reference to null.
   */
  bool noNull() const { return !isNullable; } 

  /**
   * Set the flag saying that the interface it is allowed to set the
   * reference to null.
   */
  void setNullable() { isNullable = true; } 

  /**
   * Set the flag saying that the interface it is not allowed to set
   * the reference to null.
   */
  void setNotNullable() { isNullable = false; } 

  /**
   * Get the flag saying wether a null pointer should be replaced by a
   * default of suitable class when rebind is called.
   */
  bool defaultIfNull() const { return theDefaultIfNull; } 

  /**
   * Set the flag saying that a null pointer should be replaced by a
   * default of suitable class when rebind is called.
   */
  void setDefaultIfNull() { theDefaultIfNull = true; } 

  /**
   * Set the flag saying that a null pointer should not be replaced by
   * a default of suitable class when rebind is called.
   */
  void setNoDefaultIfNull() { theDefaultIfNull = false; } 

private:

  /**
   * The class name of the class referred to by this
   * interface.
   */
  string theRefClassName;

  /**
   * The type_info object of the class referred to by this interface.
   */
  const type_info & theRefTypeInfo;

  /**
   * The flag saying whether the interface is responsible for
   * rebinding refenerces.
   */
  bool dontRebind;

  /**
   * The flag saying whether the interface is allowed to set a
   * reference to null.
   */
  bool isNullable;

  /**
   * The flag saying wether a null pointer should be replaced
   * by a default of suitable class when rebind is called.
   */
  bool theDefaultIfNull;

};

/** Dummy function to ensure that strings can be used as arguments
 *  also where numbers are assumed. */
inline double operator/(string,string) { return 0.0; }

/** Dummy function to ensure that strings can be used as arguments
 *  also where numbers are assumed. */
inline string operator*(double,string) { return ""; }

}

#endif /* ThePEG_InterfaceBaseH */
