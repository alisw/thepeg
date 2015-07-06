// -*- C++ -*-
#ifndef ARIADNE_RemnantParton_H
#define ARIADNE_RemnantParton_H
//
// This is the declaration of the RemnantParton class.
//

#include "Parton.h"
#include "RemnantParton.fh"

namespace Ariadne {

using namespace ThePEG;

/**
 * The RemnantParton class represents the remnant left after a parton
 * has been extracted from an incoming particle.
 */
class RemnantParton: public Parton {

public:

  /**
   * The DipoleState is a friend.
   */
  friend class DipoleState;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline RemnantParton();

  /**
   * The copy constructor.
   */
  inline RemnantParton(const RemnantParton &);

  /**
   * The destructor.
   */
  inline virtual ~RemnantParton();
  //@}

public:

  /** @name Simple access functions. */
  //@{
  /**
   * The ParticleData object corresponding to the particle from which
   * this parton was extracted.
   */
  inline tcPDPtr parentDataPtr() const;

  /**
   * The ParticleData object corresponding to the particle from which
   * this parton was extracted.
   */
  inline const ParticleData & parentData() const;

  /**
   * The momentum of the particle from which this parton was
   * extracted.
   */
  inline const Lorentz5Momentum & parentMomentum() const;

  /**
   * Return the inverse extension of this remant.
   */
  inline Energy mu() const;

  /**
   * Return the dimensionality of the extension of this remnant.
   */
  inline double alpha() const;

  /**
   * Set the inverse extension of this remant.
   */
  inline void mu(Energy);

  /**
   * Set the dimensionality of the extension of this remnant.
   */
  inline void alpha(double);
  //@}

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

private:

  /**
   * The inverse extension of this remant.
   */
  Energy theMu;

  /**
   * The dimensionality of the extension of this remnant.
   */
  double theAlpha;

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
  static ClassDescription<RemnantParton> initRemnantParton;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  RemnantParton & operator=(const RemnantParton &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of RemnantParton. */
template <>
struct BaseClassTrait<Ariadne::RemnantParton,1> {
  /** Typedef of the first base class of RemnantParton. */
  typedef Ariadne::Parton NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the RemnantParton class and the shared object where it is defined. */
template <>
struct ClassTraits<Ariadne::RemnantParton>
  : public ClassTraitsBase<Ariadne::RemnantParton> {
  /** Return a platform-independent class name */
  static string className() { return "Ariadne::RemnantParton"; }
  /** Return the name of the shared library be loaded to get
   *  access to the RemnantParton class and every other class it uses
   *  (except the base class). */
  static string library() { return "libArCascade.so"; }
};

/** @endcond */

}

#include "RemnantParton.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "RemnantParton.tcc"
#endif

#endif /* ARIADNE_RemnantParton_H */
