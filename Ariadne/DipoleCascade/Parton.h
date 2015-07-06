// -*- C++ -*-
#ifndef ARIADNE_Parton_H
#define ARIADNE_Parton_H
//
// This is the declaration of the Parton class.
//

#include "CascadeBase.h"
#include "Parton.fh"
#include "Dipole.fh"
#include "String.fh"
#include "ThePEG/Vectors/Lorentz5Vector.h"

namespace Ariadne {

using namespace ThePEG;

/**
 * The Parton class represents partons which are able to cascade
 * according to the Dipole Cascade Model.
 */
class Parton: public CascadeBase {

public:

  /**
   * The DipoleState is a friend.
   */
  friend class DipoleState;

  /**
   * A pair of partons.
   */
  typedef pair<tParPtr,tParPtr> tParPair;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  Parton();

  /**
   * The copy constructor.
   */
  Parton(const Parton &);

  /**
   * The destructor.
   */
  virtual ~Parton();
  //@}

public:

  /** @name Simple access functions. */
  //@{
  /**
   * If this was an original parton, this points to the corresponding
   * Particle object.
   */
  inline tcPPtr orig() const;

  /**
   * The final ThePEG::Particle produced for this parton with the
   * produceParticle() function.
   */
  inline tPPtr particle() const;

  /**
   * Fill the given iterator with the original Particle objects
   * corresponding to the ultimate parents.
   */
  template <typename OIterator>
  void getOriginalParents(OIterator it) const;

  /**
   * If this parton was created during the cascade, return the
   * pointers to the parent partons.
   */
  inline tParPair parents() const;

  /**
   * The corresponding ParticleData object.
   */
  inline tcPDPtr dataPtr() const;

  /**
   * The corresponding ParticleData object.
   */
  inline const ParticleData & data() const;

  /**
   * Return true if this parton is a gluon.
   */
  inline bool isG() const;

  /**
   * Return true if this parton is coloured.
   */
  inline bool coloured() const;

  /**
   * The momentum of this Parton.
   */
  inline const Lorentz5Momentum & momentum() const;

  /**
   * The momentum of this Parton.
   */
  inline Lorentz5Momentum & momentum();

  /**
   * The string to which this parton belongs.
   */
  inline tStrPtr string() const;

  /**
   * The dipole connecting to the anti-colour line.
   */
  inline tDipPtr iDip() const;

  /**
   * The previous parton in the string. Returns iDip()->iPart().
   */
  tParPtr prev() const;

  /**
   * The dipole connecting to the colour line.
   */
  inline tDipPtr oDip() const;

  /**
   * The next parton in the string. Returns oDip()->oPart().
   */
  tParPtr next() const;

  /**
   * If this was an original parton, set the pointer to the
   * corresponding Particle object. Also set the momentum and the
   * ParticleData pointer.
   */
  void orig(tcPPtr);

  /**
   * If this parton was created during the cascade, set the pointers
   * to the parent partons.
   */
  inline void parents(tParPair);

  /**
   * Set the data pointer for this parton.
   */
  void data(tcPDPtr);

  /**
   * Set the string to which this parton belongs.
   */
  inline void string(tStrPtr);

  /**
   * The dipole connecting to the anti-colour line.
   */
  inline void iDip(tDipPtr);

  /**
   * The dipole connecting to the colour line.
   */
  inline void oDip(tDipPtr);
  //@}

  /**
   * Produce a ThePEG::Particle corresponding to this parton. The
   * momentum of the produced particle is rotated with \a r w.r.t. the
   * parton.
   */
  virtual tPPtr produceParticle(const LorentzRotation & r);

  /**
   * Calculate the invatiant pt2 of the parton.
   */
  Energy2 invPT2() const;

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
  virtual void fillReferences(CloneSet &) const ;

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
   * If this was an original parton, this points to the corresponding
   * Particle object.
   */
  tcPPtr theOrig;

  /**
   * The final ThePEG::Particle produced for this parton with the
   * produceParticle() function.
   */
  PPtr theParticle;

  /**
   * If this parton was created during the cascade, these points to
   * the parent partons.
   */
  tParPair theParents;

  /**
   * The corresponding ParticleData object.
   */
  tcPDPtr theDataPtr;

  /**
   * True if this parton is a gluon.
   */
  bool isGluon;

  /**
   * The momentum of this Parton.
   */
  Lorentz5Momentum theMomentum;

  /**
   * The string to which this parton belongs.
   */
  tStrPtr theString;

  /**
   * The dipole connecting to the incoming (anti-)colour line.
   */
  tDipPtr theIDip;

  /**
   * The dipole connecting to the outgoing colour line.
   */
  tDipPtr theODip;

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
  static ClassDescription<Parton> initParton;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  Parton & operator=(const Parton &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of Parton. */
template <>
struct BaseClassTrait<Ariadne::Parton,1> {
  /** Typedef of the first base class of Parton. */
  typedef Ariadne::CascadeBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the Parton class and the shared object where it is defined. */
template <>
struct ClassTraits<Ariadne::Parton>
  : public ClassTraitsBase<Ariadne::Parton> {
  /** Return a platform-independent class name */
  static string className() { return "Ariadne::Parton"; }
  /** Return the name of the shared library be loaded to get
   *  access to the Parton class and every other class it uses
   *  (except the base class). */
  static string library() { return "libArCascade.so"; }
};

/** @endcond */

}

#include "Parton.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "Parton.tcc"
#endif

#endif /* ARIADNE_Parton_H */
