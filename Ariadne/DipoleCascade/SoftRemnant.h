// -*- C++ -*-
#ifndef ARIADNE_SoftRemnant_H
#define ARIADNE_SoftRemnant_H
//
// This is the declaration of the SoftRemnant class.
//

#include "RemnantParton.h"
#include "SoftRemnant.fh"

namespace Ariadne {

using namespace ThePEG;

/**
 * The SoftRemnant class represents the remnant left after a parton
 * has been extracted from an incoming particle.
 */
class SoftRemnant: public RemnantParton {

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
  inline SoftRemnant();

  /**
   * The copy constructor.
   */
  inline SoftRemnant(const SoftRemnant &);

  /**
   * The destructor.
   */
  inline virtual ~SoftRemnant();
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
   * The momentum of the incoming parton
   */
  inline LorentzMomentum incomingMomentum() const;

  /**
   * The momentum of the particle from which this parton was
   * extracted.
   */
  inline void parentMomentum(const Lorentz5Momentum & p);

  /**
   * Set the data pointer for the particle from which this parton was
   * extracted.
   */
  void parentData(tcPDPtr);

  /**
   * Get the pdf associated with the parent.
   */
  inline const PDF & pdf() const;

  /**
   * Set the pdf associated with the parent.
   */
  inline void pdf(PDF);
  //@}

  /**
   * Calculate x of the extracted parton.
   */
  double partonx() const;

  /**
   * Produce a ThePEG::Particle corresponding to this parton. The
   * momentum of the produced particle is rotated with \a r w.r.t. the
   * parton.
   */
  virtual tPPtr produceParticle(const LorentzRotation & r);

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
   * The corresponding ParticleData object.
   */
  tcPDPtr theParentDataPtr;

  /**
   * The momentum of this Parton.
   */
  Lorentz5Momentum theParentMomentum;

  /**
   * The PDF associated with the parent.
   */
  PDF thePdf;

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
  static ClassDescription<SoftRemnant> initSoftRemnant;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SoftRemnant & operator=(const SoftRemnant &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of SoftRemnant. */
template <>
struct BaseClassTrait<Ariadne::SoftRemnant,1> {
  /** Typedef of the first base class of SoftRemnant. */
  typedef Ariadne::RemnantParton NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the SoftRemnant class and the shared object where it is defined. */
template <>
struct ClassTraits<Ariadne::SoftRemnant>
  : public ClassTraitsBase<Ariadne::SoftRemnant> {
  /** Return a platform-independent class name */
  static string className() { return "Ariadne::SoftRemnant"; }
  /** Return the name of the shared library be loaded to get
   *  access to the SoftRemnant class and every other class it uses
   *  (except the base class). */
  static string library() { return "libArCascade.so"; }
};

/** @endcond */

}

#include "SoftRemnant.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "SoftRemnant.tcc"
#endif

#endif /* ARIADNE_SoftRemnant_H */
