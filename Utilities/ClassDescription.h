// -*- C++ -*-
//
// ClassDescription.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_ClassDescription_H
#define ThePEG_ClassDescription_H

#include "ThePEG/Config/ThePEG.h"
#include "ClassDescription.fh"
#include "ThePEG/Utilities/Named.h"
#include "ThePEG/Persistency/PersistentOStream.fh"
#include "ThePEG/Persistency/PersistentIStream.fh"
#include "ClassTraits.h"
#include "DescriptionList.h"

namespace ThePEG {

/**
 * ClassDescriptionBase is the base class for all class description
 * classes.  ClassDescriptionBase is the non-templated base class for
 * the templated ClassDescriptionTBase, ClassDescription,
 * AbstractClassDescription, NoPIOClassDescription and
 * AbstractNoPIOClassDescription classes. An instantiation of one of
 * these classes represents a meta class encapsulating information
 * about the template argument class. This information is obtained
 * from the templated ClassTraits class which can be specialized for
 * any given class should the default information not be
 * satisfactory. Information is also obtained from the templated
 * BaseClassTraits class which must be specialized for each class
 * declaring a typedef for each of its base classes.
 *
 * The information about a class stored in a
 * <code>ClassDescriptionBase</code> objects is the following:<BR> a
 * platform-independent class name,<BR> a reference to the
 * <code>type_info</code><BR> an integer version number,<BR> a vector
 * of <code>ClassDescriptionBase*</code> corresponding to the base
 * classes,<BR> methods for reading/writing members of a class from/to
 * a PersistentIStream/PersistentOStream and<BR> a method for creating
 * an object of a class.
 *
 * The ClassDescriptionBase objects should be created before main() is
 * executed, and shouold therfore be instantiated as static member
 * variables. Only one ClassDescriptionBase object shall be
 * instantiated for each class to be described. The
 * ClassDescriptionBase objects are automatically stored in a purely
 * static DescriptionList class.
 *
 * AbstractClassDescription and AbstractNoPIOClassDescription should
 * be used for abstract base classes, while NoPIOClassDescription and
 * AbstractNoPIOClassDescription should be used for classes which do
 * not have any members which need to be read and written
 * persistently.
 *
 * @see ClassTraits
 * @see PersistentIStream
 * @see PersistentOStream
 * @see DescriptionList
 */
class ClassDescriptionBase: public Named {

public:

  /** A vector of class descriptions. */
  typedef vector<const ClassDescriptionBase *> DescriptionVector;

protected:

  /**
   * The constructor used by sub-classes.
   * @param newName the platform independent name of the class.
   * @param newInfo the type_info object corresponding to the class.
   * @param newVersion the implementation version of the class.
   * @param newLibrary the name of a file containing the dynamic
   * @param abst true if the class is abstract.
   * library where the class is implemented.
   */
  ClassDescriptionBase(string newName, 
		       const type_info & newInfo,
		       int newVersion, 
		       string newLibrary,
		       bool abst)
    : Named(newName), theVersion(newVersion), theLibrary(newLibrary),
      theInfo(newInfo), isAbstract(abst), done(false) {}

public:

  /**
   * Empty destructor.
   */
  virtual ~ClassDescriptionBase();

  /**
   * The standart RTTI type_info object for the described class.
   */
  const type_info & info() const { return theInfo; }

  /**
   * The version of the described class.
   */
  int version() const { return theVersion; }

  /**
   * The name of a file containing the dynamic
   * library where the class is implemented.
   */
  string library() const { return theLibrary; }

  /**
   * Return true if this object was set up properly.
   */
  bool check() const { return done; }

  /**
   * Return the descriptions of the base classes of the described
   * class.
   */
  const DescriptionVector & descriptions() const { return theBaseClasses; }

  /**
   * Set up the base class information for this object.
   */
  virtual void setup() = 0;

  /**
   * Create an object of the described class.
   */
  virtual BPtr create() const = 0;

  /**
   * Output the members of an object of the described class to a
   * persistent stream.
   * @param b the object to be written.
   * @param os the persistent stream.
   */
  virtual void output(tcBPtr b, PersistentOStream & os) const = 0;

  /**
   * Read the members of an object of the described class from a
   * persistent stream.
   * @param b the object to be read.
   * @param is the persistent stream.
   * @param oldVersion the version number of the object when it was written.
   */
  virtual void input(tBPtr b, PersistentIStream & is, int oldVersion) const = 0;

  /**
   * Return true if the class described by the argument is a base
   * class of the class described by this.
   */
  bool isA(const ClassDescriptionBase & base) const;

  /**
   * Return true if the corresponding class is abstract.
   */
  bool abstract() const { return isAbstract; }

protected:

  /**
   * Set the base classes.
   * @param first an iterator refering to the first base class
   * @param last an iterator giving the end of the range of base class
   * descriptions.
   */
  void baseClasses(DescriptionVector::iterator first,
		   DescriptionVector::iterator last) 
  {
    theBaseClasses = DescriptionVector(first, last);
    done = true;
  }

private:

  /**
   * The version of the described class.
   */
  int theVersion;

  /**
   * The library file where this class may be found.
   */
  string theLibrary;

  /**
   * The type_info object for the described class
   */
  const type_info & theInfo;

  /**
   * The vector of base classes.
   */
  DescriptionVector theBaseClasses;

  /**
   * True if this class is abstract.
   */
  bool isAbstract;

  /**
   * True if this object was set up properly.
   */
  bool done;

};

/**
 * A helper class for tracing the base classes of a class to be
 * described
 */
template <typename T, int IBase,
          typename B = typename BaseClassTrait<T,IBase>::NthBase>
struct ClassDescriptionHelper {
  /** Add base classes */
  static void addBases(vector<const ClassDescriptionBase *> & c){
    const ClassDescriptionBase * b = DescriptionList::find(typeid(B));
    if ( !b ) return;
    c.push_back(b);
    ClassDescriptionHelper<T,IBase+1>::addBases(c);
  }
};

/** @cond TRAITSPECIALIZATIONS */

/**
 * A helper class for tracing the base classes of a class to be
 * described
 */
template <typename T, int IBase>
struct ClassDescriptionHelper<T, IBase, int> {
  /** Add base classes */
  static void addBases(vector<const ClassDescriptionBase *> & ) {}
};

/** @endcond */

/**
 * An intermediate templated base class derived from
 * ClassDescriptionBase.
 */
template <typename T>
class ClassDescriptionTBase: public ClassDescriptionBase {

public:

  /** The traits class for the template argument class. */
  typedef ClassTraits<T> Traits;

public:

  /**
   * Default constructor. If \a abst is true then the corresponding
   * class is abstract.
   */
  ClassDescriptionTBase(bool abst)
    : ClassDescriptionBase(Traits::className(), typeid(T), Traits::version(),
			   Traits::library(), abst) 
  {
    DescriptionList::Register(*this);
    T::Init();
  }


  /**
   * The descructor.
   */
  virtual ~ClassDescriptionTBase() {}

  /**
   * Set up the base class information for this object.
   */
  virtual void setup() {
    DescriptionVector bases;
    ClassDescriptionHelper<T,1>::addBases(bases);
    baseClasses(bases.begin(), bases.end());
  }

};

/**
 * A concreate implementation of ClassDescriptionBase describing an
 * abstract class with persistent data.
 */
template <typename T>
class AbstractClassDescription: public ClassDescriptionTBase<T> {

public:

  /** The traits class for the template argument class. */
  typedef ClassTraits<T> Traits;

public:

  /**
   * Default constructor.
   */
  AbstractClassDescription() : ClassDescriptionTBase<T>(true) {}

  /**
   * Do not create an object of the described class (which is
   * abstract). Just return the null pointer.
   */
  virtual BPtr create() const {
    throw std::logic_error("Tried to instantiate virtual class " + Named::name());
  }

  /**
   * Output the members of an object of the described class to a
   * persistent stream.
   * @param b the object to be written.
   * @param os the persistent stream.
   */
  virtual void output(tcBPtr b, PersistentOStream & os) const {
    Traits::output(Traits::cast(b), os);
  }

  /**
   * Read the members of an object of the described class from a
   * persistent stream.
   * @param b the object to be read.
   * @param is the persistent stream.
   * @param oldVersion the version number of the object when it was written.
   */
  virtual void input(tBPtr b, PersistentIStream & is,
		     int oldVersion) const {
    Traits::input(Traits::cast(b), is, oldVersion);
  }

};

/**
 * A concreate implementation of ClassDescriptionBase describing a
 * concrete class with persistent data.
 */
template <typename T>
class ClassDescription: public ClassDescriptionTBase<T> {

public:

  /** The traits class for the template argument class. */
  typedef ClassTraits<T> Traits;

public:

  /**
   * Default constructor.
   */
  ClassDescription() : ClassDescriptionTBase<T>(false) {}

  /**
   * Create an object of the described class.
   */
  virtual BPtr create() const { return Traits::create(); }

  /**
   * Output the members of an object of the described class to a
   * persistent stream.
   * @param b the object to be written.
   * @param os the persistent stream.
   */
  virtual void output(tcBPtr b, PersistentOStream & os) const {
    Traits::output(Traits::cast(b), os);
  }

  /**
   * Read the members of an object of the described class from a
   * persistent stream.
   * @param b the object to be read.
   * @param is the persistent stream.
   * @param oldVersion the version number of the object when it was written.
   */
  virtual void input(tBPtr b, PersistentIStream & is,
		     int oldVersion) const {
    Traits::input(Traits::cast(b), is, oldVersion);
  }

};

/**
 * A concreate implementation of ClassDescriptionBase describing a
 * concrete class without persistent data.
 */
template <typename T>
class NoPIOClassDescription: public ClassDescriptionTBase<T> {

public:

  /** The traits class for the template argument class. */
  typedef ClassTraits<T> Traits;

public:

  /**
   * Default constructor.
   */
  NoPIOClassDescription() : ClassDescriptionTBase<T>(false) {}

  /**
   * Create an object of the described class.
   */
  virtual BPtr create() const { return Traits::create(); }

  /**
   * Do nothing since the described class has no persistent data.
   */
  virtual void output(tcBPtr, PersistentOStream &) const {}

  /**
   * Do nothing since the described class has no persistent data.
   */
  virtual void input(tBPtr, PersistentIStream &, int) const {}

};

/**
 * A concreate implementation of ClassDescriptionBase describing an
 * abstract class without persistent data.
 */
template <typename T>
class AbstractNoPIOClassDescription: public ClassDescriptionTBase<T> {

public:

  /** The traits class for the template argument class. */
  typedef ClassTraits<T> Traits;

public:

  /**
   * Default constructor.
   */
  AbstractNoPIOClassDescription() : ClassDescriptionTBase<T>(true) {}

  /**
   * Do not create an object of the described class (which is
   * abstract). Just return the null pointer.
   */
  virtual BPtr create() const {
    throw std::logic_error("Tried to instantiate virtual class " + Named::name());
  }

  /**
   * Do nothing since the described class has no persistent data.
   */
  virtual void output(tcBPtr, PersistentOStream & ) const {}

  /**
   * Do nothing since the described class has no persistent data.
   */
  virtual void input(tBPtr, PersistentIStream &, int) const {}

};

}

#define ThePEG_DECLARE_CLASS_DESCRIPTION(Class)                    \
/** Describe a concrete class with persistent data. */             \
static ClassDescription<Class> init ## Class                       \

#define ThePEG_DECLARE_ABSTRACT_CLASS_DESCRIPTION(Class)           \
/** Describe an abstract class with persistent data. */            \
static AbstractClassDescription<Class> init ## Class               \

#define ThePEG_DECLARE_NOPIO_CLASS_DESCRIPTION(Class)              \
/** Describe a concrete class without persistent data. */          \
static NoPIOClassDescription<Class> init ## Class                  \

#define ThePEG_DECLARE_ABSTRACT_NOPIO_CLASS_DESCRIPTION(Class)     \
/** Describe an abstract class without persistent data. */         \
static AbstractNoPIOClassDescription<Class> init ## Class          \

#define ThePEG_IMPLEMENT_CLASS_DESCRIPTION(Class)                  \
ClassDescription<Class> Class::init ## Class                       \

#define ThePEG_IMPLEMENT_ABSTRACT_CLASS_DESCRIPTION(Class)         \
AbstractClassDescription<Class> Class::init ## Class               \

#define ThePEG_IMPLEMENT_NOPIO_CLASS_DESCRIPTION(Class)            \
NoPIOClassDescription<Class> Class::init ## Class                  \

#define ThePEG_IMPLEMENT_ABSTRACT_NOPIO_CLASS_DESCRIPTION(Class)   \
AbstractNoPIOClassDescription<Class> Class::init ## Class          \

#endif /* ThePEG_ClassDescription_H */
