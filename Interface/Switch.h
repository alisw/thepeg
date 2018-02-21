// -*- C++ -*-
//
// Switch.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_Switch_H
#define ThePEG_Switch_H
// This is the declaration of the Switch, SwitchBase and SwitchOption classes.

#include "ThePEG/Config/ThePEG.h"
#include "Switch.fh"
#include "Switch.xh"
#include "InterfaceBase.h"

namespace ThePEG {

/**
 * SwitchOption is used by the Switch class and its base class
 * SwitchBase to define valid options in a switch.
 *
 * For each InterfacedBase class exactly one static Switch object
 * should created for each member variable which should be
 * interfaced. This object will automatically register itself with the
 * BaseRepository class. Also for each Switch object exactly one
 * static SwitchOption object should be created for each valid integer
 * option.
 *
 * @see InterfacedBase
 * @see InterfacedBase
 * @see Named
 * 
 */
class SwitchOption: public Named {

public:

  /**
   * Standard constructor.
   *
   * @param theSwitch the Switch object for which this option is
   * defined. Note thet the static Switch object must be created
   * before this is created.
   *
   * @param newName the name of the option, may only contain
   * letters [a-zA-z0-9_].
   *
   * @param newDescription a brief description of the option.
   *
   * @param newValue the integer value corresponding to this option.
   */
  template<typename EnumT> 
  SwitchOption(SwitchBase & theSwitch, string newName,
	       string newDescription, EnumT newValue);

  /**
   * Default constructor.
   */
  SwitchOption() : theValue(-999) {}

  /**
   * The description of this option
   */
  const string & description() const { return theDescription; }

  /**
   * The value of this option.
   */
  long value() const { return theValue; }

  /**
   * The value of this option.
   */
  operator long () const;

protected:

private:

  /**
   * The description of this option
   */
  string theDescription;

  /**
   * The value of this option.
   */
  long theValue;

};

/**
 * The Switch class and its base class SwitchBase defines an interface
 * to a class derived from the InterfacedBase, through which simple
 * integer member variables can be manuipulated and set to a
 * pre-defined set of values (options). Switch is
 * templated on the type of the integer member variable and the type
 * of the class, and is derived from the InterfaceBase class via
 * <code>SwitchBase</code>.
 *
 * The Switch class has a set of Named SwitchOptions,
 * which limits the values possible to set.
 *
 * For each InterfacedBase class exactly one static Switch object
 * should created for each member variable which should be
 * interfaced. This object will automatically register itself with the
 * BaseRepository class. Also for each Switch object exactly one
 * static SwitchOption object should be created for each valid integer
 * option.
 *
 * @see InterfacedBase
 * @see InterfacedBase
 * @see Named
 * 
 */
class SwitchBase: public InterfaceBase {

public:

  /** A map with SwitchOptions indexed by their values. */
  typedef map<long, SwitchOption> OptionMap;
  /** A map with SwitchOptions indexed by their names. */
  typedef map<string, SwitchOption> StringMap;

  /** SwitchOption is a friend. */
  friend class SwitchOption;

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
  SwitchBase(string newName, string newDescription,
	     string newClassName, const type_info & newTypeInfo,
	     bool depSafe, bool readonly)
    : InterfaceBase(newName, newDescription, newClassName, 
		    newTypeInfo, depSafe, readonly) {}

  /**
   * The general interface method overriding the one in
   * InterfaceBase. For this class, \a action can be any of "set",
   * "get", "def" and "setdef" and \a argument should be a something
   * which can be read into an integer variable through a stringstream
   * with the standard '>>' operator.
   */
  virtual string exec(InterfacedBase & ib, string action,
		    string arguments) const;

  /**
   * Return a complete description of this switch.
   */
  virtual string fullDescription(const InterfacedBase & ib) const;

  /**
   * Return a code for the type of this switch.
   */
  virtual string type() const;

  /**
   * Set the member variable of \a ib to \a val.
   */
  virtual void set(InterfacedBase & ib, long val)
    const = 0;

  /**
   * Return the value of the member variable of \a ib.
   */
  virtual long get(const InterfacedBase & ib)
    const = 0;

  /**
   * Return the default value for the member variable of \a ib.
   */
  virtual long def(const InterfacedBase & ib)
    const = 0;

  /**
   * Set the member variable of \a ib to its default value.
   */
  void setDef(InterfacedBase & i) const {
    set(i, def(i));
  }

  /**
   * Check if \a val is among the listed options.
   */
  bool check(long newValue) const { return member(theOptions, newValue); }

  /**
   * Return the map relating options to their values
   */
  const OptionMap & options() const { return theOptions; }

  /**
   * Return a string describing the type of interface to be included
   * in the Doxygen documentation.
   */
  virtual string doxygenType() const;

  /**
   * Return a string with the option index and its associated tag.
   */
  string opttag(long opt) const;

protected:

  /**
   * Register a new option.
   */
  void registerOption(const SwitchOption & o) {
    theOptions[o.value()] = o;
    theOptionNames[o.name()] = o;
  }

private:

  /**
   * The map relating options to their values
   */
  OptionMap theOptions;

  /**
   * The map relating options to their names
   */
  StringMap theOptionNames;

};

/**
 * The Switch class and its base class SwitchBase defines an interface
 * to a class derived from the InterfacedBase, through which simple
 * integer member variables can be manuipulated and set to a
 * pre-defined set of values (options). Switch is templated on the
 * type of the integer member variable (also enums and bool are
 * allowed) and the type of the class, and is derived from the
 * InterfaceBase class via SwitchBase.
 *
 * The Switch class has a set of Named SwitchOptions,
 * which limits the values possible to set.
 *
 * For each InterfacedBase class exactly one static Switch object
 * should created for each member variable which should be
 * interfaced. This object will automatically register itself with the
 * BaseRepository class. Also for each Switch object exactly one
 * static SwitchOption object should be created for each valid integer
 * option.
 *
 * @see InterfacedBase
 * @see InterfacedBase
 * @see Named
 * 
 */
template <typename T, typename Int>
class Switch: public SwitchBase {

public:

  /**
   * The declaration of member functions which can be used by this
   * Switch interface for the 'set' action.
   */
  typedef void (T::*SetFn)(Int);
  /**
   * The declaration of member functions which can be used by this
   * Switch interface for the 'get' action.
   */
  typedef Int (T::*GetFn)() const;

  /**
   * Declaration of a direct pointer to the member variable.
   */
  typedef Int T::* Member;

public:

  /**
   * Standard constructor.
   *
   * @param newName the name of the interface, may only contain
   * letters [a-zA-z0-9_].
   *
   * @param newDescription a brief description of the interface.
   *
   * @param newMember a pointer to the member variable. May be null if
   * corresponding set/get functions are provided.
   *
   * @param newDef the default value for the member variable.
   *
   * @param depSafe set to true if calls to this interface for one
   * object does not influence other objects.
   *
   * @param readonly if this is set true the interface will not be
   * able to manipulate objects of the corresponding class, but will
   * still be able to access information.
   *
   * @param newSetFn optional pointer to the member function for the
   * 'set' action.
   *
   * @param newGetFn optional pointer to the member function for the
   * 'get' action.
   *
   * @param newDefFn optional pointer to the member function for the
   * 'def' action.
   */
  Switch(string newName, string newDescription,
	 Member newMember, Int newDef, bool depSafe = false,
	 bool readonly = false, SetFn newSetFn = 0, GetFn newGetFn = 0,
	 GetFn newDefFn = 0)
    : SwitchBase(newName, newDescription, ClassTraits<T>::className(),
		 typeid(T), depSafe, readonly),
      theMember(newMember), theDef(newDef), theSetFn(newSetFn),
      theGetFn(newGetFn), theDefFn(newDefFn) {}
  
  /**
   * Set the member variable of \a ib to \a val.
   */
  virtual void set(InterfacedBase & ib, long val) const
   ;

  /**
   * Return the value of the member variable of \a ib.
   */
  virtual long get(const InterfacedBase & ib) const;

  /**
   * Return the default value for the member variable of \a ib.
   */
  virtual long def(const InterfacedBase & ib) const;

  /**
   * Give a pointer to a member function to be used by 'set()'.
   */
  void setSetFunction(SetFn sf) { theSetFn = sf; }

  /**
   * Give a pointer to a member function to be used by 'get()'.
   */
  void setGetFunction(GetFn gf) { theGetFn = gf; }

  /**
   * Give a pointer to a member function to be used by 'def()'.
   */
  void setDefaultFunction(GetFn df) { theDefFn = df; }

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
   * pointers are given.
   */
  Int theDef;

  /**
   * A pointer to a member function to be used by 'set()'.
   */
  SetFn theSetFn;

  /**
   * Pointer to member function to be used by get().
   */
  GetFn theGetFn;  

  /**
   * Pointer to member function to be used by def().
   */
  GetFn theDefFn;  

};

}

#ifndef ThePEG_TEMPLATES_IN_CC_FILE
#include "Switch.tcc"
#endif

#endif /* ThePEG_Switch_H */
