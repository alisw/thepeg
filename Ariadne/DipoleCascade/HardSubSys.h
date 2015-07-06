// -*- C++ -*-
#ifndef ARIADNE_HardSubSys_H
#define ARIADNE_HardSubSys_H
//
// This is the declaration of the HardSubSys class.
//

#include "ThePEG/Config/ThePEG.h"
#include "CascadeBase.h"
#include "DipoleState.fh"
#include "Emitter.fh"
#include "Parton.fh"
#include "String.fh"
#include "HardRemnant.fh"

namespace Ariadne {

using namespace ThePEG;

/**
 * The HardSubSys class is used by DipoleState to store information
 * about the hard sub-system (i.e. outgoing particles from the hard
 * sub-process). The class also keeps track of added particles
 * (coloured or not) and possible Lorentz rotations.
 */
class HardSubSys: public CascadeBase {

public:

  /**
   * The DipoleState is a friend.
   */
  friend class DipoleState;

  /**
   * A set of Parton pointers.
   */
  typedef set<tParPtr> PartonSet;

  /**
   * A pair of hard remnants.
   */
  typedef pair<tHardRemPtr,tHardRemPtr> tHardRemPair;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline HardSubSys();

  /**
   * The copy constructor.
   */
  inline HardSubSys(const HardSubSys &);

  /**
   * The destructor.
   */
  virtual ~HardSubSys();
  //@}

public:

  /**
   * Add a new particle to the set of active partons or the set of
   * produced colour-singlet particles.
   */
  inline void add(tParPtr p);

  /**
   * Remove a particle from the set of active partons or the set of
   * produced colour-singlet particles.
   */
  inline bool remove(tParPtr p);

  /**
   * Add a final state (colour-singlet) particle from the initial
   * sub-process.
   */
  void add(tPPtr p);

  /**
   * Add an intermediate particle from the initial sub-process.
   */
  inline void addIntermediate(tPPtr p);

  /**
   * Return the total momentum of this hard sub-system.
   */
  inline const Lorentz5Momentum & momentum() const;

  /**
   * Set the total momentum of this hard sub-system. If posZDir is true
   * the value of the incoming particle with positive momentum along the
   * z axis is modified.
   */
  void setMomentum(const LorentzMomentum & q, bool posZDir, double x = 1.0);

  /**
   * Perform a Lorentz transformation on the sub-system.
   */
  void transform(const LorentzRotation & r);

  /**
   * Access the set of partons which are active in the cascade.
   */
  inline const PartonSet & active() const;

  /**
   * Access the set of created (colour-singlet) particles from the cascade.
   */
  inline const PartonSet & produced() const;

  /**
   * Access the vector of final state (colour-singlet) particles in
   * the initial sub-process.
   */
  inline const tPVector & initial() const;

  /**
   * Access the vector of intermediate particles in the initial
   * sub-process.
   */
  inline const tPVector & intermediates() const;

  /**
   * Return the aggregated LorentzRotation to be used on the
   * initial() and intermediates() to obtain their final momentum.
   */
  inline const LorentzRotation & totalRotation() const;

  /**
   * Indicate whether the initial() particles have been rotated.
   */
  inline bool rotated();

  /**
   * Get the pair of hard remnant partons.
   */
  inline const tHardRemPair & hardRemnants() const;

  /**
   * Set the a \a hard remnant partons. If it is not the \a first hard
   * remnant, it is the second.
   */
  inline void hardRemnant(const tHardRemPtr & hard, bool first);

  /**
   * Return the sum of Q2 for the hard remnants.
   */
  Energy2 Q2();

protected:

  /**
   * Recalculate the total momentum. Is called from momentum() when
   * necessary.
   */
  void sumMomentum() const;

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
  virtual void rebind(const CloneBase::TranslationMap & trans);
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
   * The set of partons which are active in the cascade.
   */
  PartonSet theActivePartons;

  /**
   * The set of created (colour-singlet) particles from the cascade.
   */
  PartonSet theProducedParticles;

  /**
   * A vector of final state (colour-singlet) particles in the initial
   * sub-process. Coloured particles are stored in theActivePartons.
   */
  tPVector theInitialParticles;

  /**
   * A vector of intermediate particles in the initial sub-process.
   */
  tPVector theIntermediateParticles;

  /**
   * A pair of hard remnants.
   */
  tHardRemPair theHardRemnants;

  /**
   * The total momentum of the hard sub-system.
   */
  mutable Lorentz5Momentum theMomentum;

  /**
   * Indicate if the total momentum needs to be recalculated.
   */
  mutable bool isModified;

  /**
   * The sum of the momenta of theInitialParticles. If the sub-system
   * is transformed only this is changed, while theInitialParicles are
   * untouched. Their momenta can instead be modified afterwards using
   * theRotation.
   */
  Lorentz5Momentum theInitMomentum;

  /**
   * The accumulated LorentzRotation to be applied to
   * theInitialParticles and theIntermediateParticles.
   */
  LorentzRotation theRotation;

  /**
   * Indicate if rotations have been done.
   */
  bool isRotated;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<HardSubSys> initHardSubSys;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  HardSubSys & operator=(const HardSubSys &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of HardSubSys. */
template <>
struct BaseClassTrait<Ariadne::HardSubSys,1> {
  /** Typedef of the first base class of HardSubSys. */
  typedef PersistentBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the HardSubSys class and the shared object where it is defined. */
template <>
struct ClassTraits<Ariadne::HardSubSys>
  : public ClassTraitsBase<Ariadne::HardSubSys> {
  /** Return a platform-independent class name */
  static string className() { return "Ariadne::HardSubSys"; }
  /**
   * The name of a file containing the dynamic library where the class
   * HardSubSys is implemented. It may also include several, space-separated,
   * libraries if the class HardSubSys depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "libArCascade.so"; }
};

/** @endcond */

}

#include "HardSubSys.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "HardSubSys.tcc"
#endif

#endif /* ARIADNE_HardSubSys_H */
