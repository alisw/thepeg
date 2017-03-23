// -*- C++ -*-
//
// DescribeClass.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_DescribeClass_H
#define ThePEG_DescribeClass_H

#include "ThePEG/Utilities/ClassDescription.h"

namespace ThePEG {

/**
 * Helper class for multiple base classes. If a class to be described
 * by DescribeClass has multiple base classes, B1 and B2, they should
 * be encoded by BaseClasses<B1,B2>. This can be done for up to four
 * base classes.
 */
template <typename BaseT1, typename BaseT2 = int,
          typename BaseT3 = int, typename BaseT4 = int>
struct BaseClasses {};

/**
 * Traits class used by DescribeCLassT for transparent handling of one
 * base class or a several base classes encoded by BaseClasses.
 */
template <typename T>
struct BaseClassesTraits {

  /** The first base class */
  typedef T Base1;

  /** The second base class */
  typedef int Base2;

  /** The third base class */
  typedef int Base3;

  /** The fourth base class */
  typedef int Base4;

};

/**
 * Traits class used by DescribeCLassT for transparent handling of one
 * base class or a several base classes encoded by BaseClasses.
 */
template <typename BaseT1, typename BaseT2, typename BaseT3, typename BaseT4>
struct BaseClassesTraits< BaseClasses<BaseT1,BaseT2,BaseT3,BaseT4> > {

  /** The first base class */
  typedef BaseT1 Base1;

  /** The second base class */
  typedef BaseT2 Base2;

  /** The third base class */
  typedef BaseT3 Base3;

  /** The fourth base class */
  typedef BaseT4 Base4;

};

/**
 * Helper class used by DescribeClassT for transparent handling of
 * classes with and without persistent I/O functions.
 */
template <typename T, bool NoPIO>
struct DescribeClassPIOHelper {

  /**
   * Call standard output function.
   */
  static void output(const T & t, PersistentOStream & os) {
    t.persistentOutput(os);
  }

  /**
   * Call standard input function.
   */
  static void input(T & t, PersistentIStream & is, int oldVersion) {
    t.persistentInput(is, oldVersion);
  }

};

/**
 * Helper class used by DescribeClassT for transparent handling of
 * classes with and without persistent I/O functions.
 */
template <typename T>
struct DescribeClassPIOHelper<T,true> {

  /**
   * Do nothing as T has no persistent I/O functions.
   */
  static void output(const T &, PersistentOStream &) {}

  /**
   * Do nothing as T has no persistent I/O functions.
   */
  static void input(T &, PersistentIStream &, int) {}

};

/** 
 * Helper class used by DescribeClassT for transparent handling of
 * abstract and concrete classes.
 */
template <typename T, bool abstract>
struct DescribeClassAbstractHelper {

  /**
   * Default-creat an object.
   */
  static typename ThePEG::Ptr<T>::pointer create() {
    return new_ptr(T());
  }
};

/** 
 * Helper class used y DescribeClassT for transparent handling of
 * abstract and concrete classes.
 */
template <typename T>
struct DescribeClassAbstractHelper<T,true> {

  /**
   * Throw an exception as the class is bastract.
   */
  static typename ThePEG::Ptr<T>::pointer create() {
    throw std::logic_error("Tried to instantiate abstract class " +
			   ClassTraits<T>::className());
   }

};




/**
 * DescribeClassT and its derived companion classes DescribeClass
 * DescribeAbstractClass, DescribeNoPIOClass and
 * DescribeAbstractNoPIOClass, is a simplified interface to the type
 * information system in ThePEG. For simple classes there is no need
 * to specialize the ClassTraits and BaseClassTrait classes and to
 * have a static member variable of ClassDescription in the class (as
 * in the full ThePEG type info system). Instead it is enough to have
 * one statically initialized variable of one of the DescraibeClass
 * classes for each class. The Abstract and NoPIO versions of this
 * class should be used for abstract classes and classes without
 * persistent I/O functions respectively.
 */
template <typename T, typename BaseT, bool Abstract = false, bool NoPIO = false>
class DescribeClassT: public ClassDescriptionBase {

public:

  ThePEG_DECLARE_TEMPLATE_POINTERS(T,TPtr);
  ThePEG_DECLARE_POINTERS(Base,BPtr);

  /**
   * Constructor taking the name of the class, the dynamic library
   * where it is located and an optional version number as argument.
   */
  DescribeClassT(string cname, string lib, int vers = 0)
    : ClassDescriptionBase(cname, typeid(T), vers, lib, Abstract) {
    DescriptionList::Register(*this);
    T::Init();
  }

  /**
   * The descructor.
   */
  virtual ~DescribeClassT() {}

  /**
   * Set up the base class information for this object.
   */
  virtual void setup() {
    DescriptionVector bases;
    const ClassDescriptionBase * b =
      DescriptionList::find(typeid(typename BaseClassesTraits<BaseT>::Base1));
    if ( b ) bases.push_back(b);
    b = DescriptionList::find(typeid(typename BaseClassesTraits<BaseT>::Base2));
    if ( b ) bases.push_back(b);
    b = DescriptionList::find(typeid(typename BaseClassesTraits<BaseT>::Base3));
    if ( b ) bases.push_back(b);
    b = DescriptionList::find(typeid(typename BaseClassesTraits<BaseT>::Base4));
    if ( b ) bases.push_back(b);
    baseClasses(bases.begin(), bases.end());
  }

  /**
   * Default-create an object.
   */
  virtual BPtr create() const {
    return DescribeClassAbstractHelper<T,Abstract>::create();
  }

  /**
   * Call standard output function.
   */
  virtual void output(tcBPtr o, PersistentOStream & os) const {
    tcTPtr t = dynamic_ptr_cast<tcTPtr>(o);
    DescribeClassPIOHelper<T,NoPIO>::output(*t, os);
  }

  /**
   * Call standard input function.
   */
  virtual void input(tBPtr o, PersistentIStream & is, int oldVersion) const {
    tTPtr t = dynamic_ptr_cast<tTPtr>(o);
    DescribeClassPIOHelper<T,NoPIO>::input(*t, is, oldVersion);
  }

};
   
/**
 * DescribeClass and its companion classes DescribeAbstractClass,
 * DescribeNoPIOClass and DescribeAbstractNoPIOClass, is a simplified
 * interface to type information system in ThePEG. For simple classes
 * there is no need to specialize the ClassTraits and BaseClassTrait
 * classes and to have a static member variable of ClassDescription in
 * the class (as in the full ThePEG type info system). Instead it is
 * enough to have one statically initialized variable of one of the
 * DescraibeClass classes for each class. The Abstract and NoPIO
 * versions of this class may be used for abstract classes and
 * classes without persistent I/O functions respectively.
 *
 * The template arguments are as follows:
 * T       : the class being described. 
 * BaseT   : The base class of T. If several base classes these should be
 *           encoded s template arguments to the BaseClasses class
 *           (at most four base classes can be encoded).
 * Abstract: shoule be set to true if class T is abstract.
 * NoPIO   : shoule be set to true if class T does not persistent I/O functions.
 */
template <typename T, typename BaseT = int,
  bool Abstract = false, bool NoPIO = false>
class DescribeClass: public DescribeClassT<T,BaseT,Abstract,NoPIO> {

public:

  /**
   * Constructor taking the name of the class, the dynamic library
   * where it is located and an optional version number as argument.
   */
  DescribeClass(string cname, string lib, int vers = 0)
    : DescribeClassT<T,BaseT,Abstract,NoPIO>(cname, lib, vers) {}

};
   
/**
 * DescribeClass and its companion classes DescribeAbstractClass,
 * DescribeNoPIOClass and DescribeAbstractNoPIOClass, is a simplified
 * interface to type information system in ThePEG. For simple classes
 * there is no need to specialize the ClassTraits and BaseClassTrait
 * classes and to have a static member variable of ClassDescription in
 * the class (as in the full ThePEG type info system). Instead it is
 * enough to have one statically initialized variable of one of the
 * DescraibeClass classes for each class. The Abstract and NoPIO
 * versions of this class may be used for abstract classes and
 * classes without persistent I/O functions respectively.
 *
 * The template arguments are as follows:
 * T       : the class being described. 
 * BaseT   : The base class of T. If several base classes these should be
 *           encoded s template arguments to the BaseClasses class
 *           (at most four base classes can be encoded).
 */
template <typename T, typename BaseT = int>
class DescribeNoPIOClass: public DescribeClassT<T,BaseT,false,true> {

public:

  /**
   * Constructor taking the name of the class, the dynamic library
   * where it is located and an optional version number as argument.
   */
  DescribeNoPIOClass(string cname, string lib, int vers = 0)
    : DescribeClassT<T,BaseT,false,true>(cname, lib, vers) {}

};
   
/**
 * DescribeClass and its companion classes DescribeAbstractClass,
 * DescribeNoPIOClass and DescribeAbstractNoPIOClass, is a simplified
 * interface to type information system in ThePEG. For simple classes
 * there is no need to specialize the ClassTraits and BaseClassTrait
 * classes and to have a static member variable of ClassDescription in
 * the class (as in the full ThePEG type info system). Instead it is
 * enough to have one statically initialized variable of one of the
 * DescraibeClass classes for each class. The Abstract and NoPIO
 * versions of this class may be used for abstract classes and
 * classes without persistent I/O functions respectively.
 *
 * The template arguments are as follows:
 * T       : the class being described. 
 * BaseT   : The base class of T. If several base classes these should be
 *           encoded s template arguments to the BaseClasses class
 *           (at most four base classes can be encoded).
 */
template <typename T, typename BaseT = int>
class DescribeAbstractClass: public DescribeClassT<T,BaseT,true,false> {

public:

  /**
   * Constructor taking the name of the class, the dynamic library
   * where it is located and an optional version number as argument.
   */
  DescribeAbstractClass(string cname, string lib, int vers = 0)
    :DescribeClassT<T,BaseT,true,false>(cname, lib, vers) {}

};
   
/**
 * DescribeClass and its companion classes DescribeAbstractClass,
 * DescribeNoPIOClass and DescribeAbstractNoPIOClass, is a simplified
 * interface to type information system in ThePEG. For simple classes
 * there is no need to specialize the ClassTraits and BaseClassTrait
 * classes and to have a static member variable of ClassDescription in
 * the class (as in the full ThePEG type info system). Instead it is
 * enough to have one statically initialized variable of one of the
 * DescraibeClass classes for each class. The Abstract and NoPIO
 * versions of this class may be used for abstract classes and
 * classes without persistent I/O functions respectively.
 *
 * The template arguments are as follows:
 * T       : the class being described. 
 * BaseT   : The base class of T. If several base classes these should be
 *           encoded s template arguments to the BaseClasses class
 *           (at most four base classes can be encoded).
 */
template <typename T, typename BaseT = int>
class DescribeAbstractNoPIOClass: public DescribeClassT<T,BaseT,true,true> {

public:

  /**
   * Constructor taking the name of the class, the dynamic library
   * where it is located and an optional version number as argument.
   */
  DescribeAbstractNoPIOClass(string cname, string lib, int vers = 0)
    :DescribeClassT<T,BaseT,true,true>(cname, lib, vers) {}

};
   
}


#endif /* ThePEG_DescribeClass_H */
