// -*- C++ -*-
#ifndef ARIADNE_HardRemnant_H
#define ARIADNE_HardRemnant_H
//
// This is the declaration of the HardRemnant class.
//

#include "RemnantParton.h"
#include "HardRemnant.fh"

namespace Ariadne {

using namespace ThePEG;

/**
 * The HardRemnant class represent the struct quark in DIS.
 */
class HardRemnant: public RemnantParton {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline HardRemnant();

  /**
   * The copy constructor.
   */
  inline HardRemnant(const HardRemnant &);

  /**
   * The destructor.
   */
  virtual ~HardRemnant();
  //@}

public:

  /**
   * Return Q2 associated with the parton.
   */
  inline Energy2 Q2() const;

  /**
   * Set Q2 associated with the parton.
   */
  inline void setQ2(Energy2 Q2);

private:

  /**
   * Q2 associated with the parton.
   */
  Energy2 theQ2;

protected:

  /** @name Functions relating to the DipoleState to which this belongs. */
  //@{
  /**
   * Return a simple clone of this object. Should be implemented as
   * <code>return new_ptr(*this);</code> by a derived class.
   */
  virtual ClonePtr clone() const;

  /**
   * Fill the provided set with all pointers to CloneBase objects used
   * in this object.
   */
  virtual void fillReferences(CloneSet &) const;

  /**
   * Rebind pointers to other CloneBase objects. Called after a number
   * of interconnected CloneBase objects have been cloned, so that
   * the cloned objects will refer to the cloned copies afterwards.
   *
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   */
  virtual void rebind(const TranslationMap & trans);
  //@}

public:

  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void persistentInput(PersistentIStream & is, int version);
  //@}

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

public:

  /**
   * Print out debugging information on std::cerr.
   */
  virtual void debugme() const;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<HardRemnant> initHardRemnant;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  HardRemnant & operator=(const HardRemnant &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of HardRemnant. */
template <>
struct BaseClassTrait<Ariadne::HardRemnant,1> {
  /** Typedef of the first base class of HardRemnant. */
  typedef Ariadne::RemnantParton NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the HardRemnant class and the shared object where it is defined. */
template <>
struct ClassTraits<Ariadne::HardRemnant>
  : public ClassTraitsBase<Ariadne::HardRemnant> {
  /** Return a platform-independent class name */
  static string className() { return "Ariadne::HardRemnant"; }
  /**
   * The name of a file containing the dynamic library where the class
   * HardRemnant is implemented. It may also include several, space-separated,
   * libraries if the class HardRemnant depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "libArCascade.so"; }
};

/** @endcond */

}

#include "HardRemnant.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "HardRemnant.tcc"
#endif

#endif /* ARIADNE_HardRemnant_H */
