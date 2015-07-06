// -*- C++ -*-
#ifndef Ariadne5_Parton_H
#define Ariadne5_Parton_H
//
// This is the declaration of the Parton class.
//

#include "CascadeBase.h"
#include "Parton.fh"
#include "QCDDipole.fh"
#include "Emission.h"
#include "ThePEG/Vectors/Lorentz5Vector.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/ParticleTraits.h"

namespace Ariadne5 {

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

  /**
   * Typedef for position in transverse coordinate space.
   */
  typedef Transverse<InvEnergy> BPos;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor. Derived classes which are complex
   * partons must set the \a spec argument to true.
   */
  Parton(bool spec = false);

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
  inline tPPtr orig() const {
    return theOrig;
  }

  /**
   * The Emission that produced this parton, or null if original.
   */
  inline tcEmPtr emission() const {
    return theEmission;
  }

  /**
   * The final ThePEG::Particle produced for this parton with the
   * produceParticle() function.
   */
  inline tPPtr particle() const {
    return theParticle;
  }

  /**
   * Fill the given iterator with the original Particle objects
   * corresponding to the ultimate parents.
   */
  template <typename OIterator>
  void getOriginalParents(OIterator it) const {
    for ( int i = 0, N = adoptedOriginals.size(); i < N; ++i )
      *it++ = adoptedOriginals[i];
    if ( orig() ) *it++ = orig();
    else {
      if ( emission()->colourParent )
	emission()->colourParent->getOriginalParents(it);
      if ( emission()->antiColourParent )
	emission()->antiColourParent->getOriginalParents(it);
    }
  }

  /**
   * If this parton was created during the cascade, return the
   * pointers to the parent partons. The first parent is the coloured
   * parent and the second is the anti-colour one.
   */
  tParPair parents() const;

  /**
   * The corresponding ParticleData object.
   */
  inline tcPDPtr dataPtr() const {
    return theDataPtr;
  }

  /**
   * The corresponding ParticleData object.
   */
  inline const ParticleData & data() const {
    return *theDataPtr;
  }

  /**
   * Return the colour charge of this parton.
   */
  inline bool hasColour() const {
    return dataPtr() && data().iColour() != PDT::Colour0;
  }

  /**
   * Return the colour charge of this parton.
   */
  inline PDT::Colour iColour() const {
    return dataPtr()? data().iColour(): PDT::Colour0;
  }

  /**
   * The original incoming anti-colour line corresponding to this
   * parton.
   */
  tColinePtr origICol() const {
    return theICol;
  }

  /**
   * The original outgoing colour line corresponding to this
   * parton.
   */
  tColinePtr origOCol() const {
    return theOCol;
  }

  /**
   * The original incoming anti-colour line corresponding to this
   * parton.
   */
  void origICol(tColinePtr cl) {
    theICol = cl;
  }

  /**
   * The original outgoing colour line corresponding to this
   * parton.
   */
  void origOCol(tColinePtr cl) {
    theOCol = cl;
  }

  /**
   * Return the charge of this parton in umits of e.
   */
  inline double charge() const {
    return dataPtr()? double(data().charge()/Units::eplus): 0.0;
  }

  /**
   * Return the charge of this parton in umits of e.
   */
  inline bool hasCharge() const {
    return dataPtr() && data().iCharge() != PDT::Charge0;
  }

  /**
   * Return true if this parton is a gluon.
   */
  inline bool isG() const {
    return iColour() == PDT::Colour8;
  }

  /**
   * Return true if this parton is coloured.
   */
  inline bool coloured() const {
    return dataPtr() && data().coloured();
  }

  /**
   * Return true if this is not a normal parton in a dipole together
   * with an \a other parton.
   */
  virtual bool special(tcParPtr other = tcParPtr()) const {
    return isSpecial;
  }

  /**
   * Return true if this parton was emitted in a failsafe Emission.
   */
  bool failsafe() const {
    return emission() && emission()->failsafe;
  }

  /**
   * Return true if this parton has information about neighboring
   * dipoles or partons and if any of these hwve been touched.
   */
  virtual bool touchedNeighbours(tcQCDPtr d) const {
    return false;
  }

  /**
   * The momentum of this Parton.
   */
  inline const Lorentz5Momentum & momentum() const {
    return theMomentum;
  }

  /**
   * The momentum of this Parton.
   */
  inline Lorentz5Momentum & momentum() {
    return theMomentum;
  }

  /**
   * Set the momentum of this parton.
   */
  virtual void setMomentum(const Lorentz5Momentum & p) {
    theMomentum = p;
  }

  /**
   * The production vertex position of this Parton.
   */
  inline const LorentzPoint & vertex() const {
    return theVertex;
  }

  /**
   * The production vertex position of this Parton.
   */
  inline LorentzPoint & vertex() {
    return theVertex;
  }

  /**
   * Set the production vertex position of this Parton.
   */
  virtual void setVertex(const LorentzPoint p) {
    theVertex = p;
  }

  /**
   * Get the position in impact parameter space in inverse energy units
   */
  BPos bPos() const {
    return BPos(vertex().x()/Constants::hbarc, vertex().y()/Constants::hbarc);
  }

  /**
   * Set the position in impact parameter space in inverse energy units
   */
  const BPos & bPos(const BPos & p) {
    theVertex.setX(p.x()*hbarc);
    theVertex.setY(p.y()*hbarc);
    return p;
  }

  /**
   * If this was an original parton, set the pointer to the
   * corresponding Particle object. Also set the momentum and the
   * ParticleData pointer.
   */
  void orig(tPPtr);

  /**
   * If this parton has absorbed a gluon, it may also adopt that gluons original parton.
   */
  void adoptOriginals(tParPtr op);

  /**
   * Set the data pointer for this parton.
   */
  void data(tcPDPtr);

  /**
   * Set the Emission which produced this parton, or null if original.
   */
  void emission(tcEmPtr e);
  //@}

  /**
   * Produce a ThePEG::Particle corresponding to this parton. The
   * momentum of the produced particle is rotated with \a r w.r.t. the
   * parton.
   */
  virtual tPPtr produceParticle(const LorentzRotation & r = LorentzRotation());

  /**
   * Calculate the invatiant pt2 of the parton.
   */
  Energy2 invPT2() const;

  /**
   * Some sub classes may need information about emissions which may
   * influence their state.
   */
  virtual void notify(const Emission &);

  /**
   * Return true if this parton is present in the final state.
   */
  bool finalState() const {
    return inFinalState;
  }

  /**
   * Set true if this parton is present in the final state.
   */
  void setFinalState(bool yes = true) {
    inFinalState = yes;
  }

  /**
   * Return the dipole system to which this parton currently belongs.
   */
  int system() const {
    return theSystem;
  }

  /**
   * Set the dipole system to which this parton currently belongs.
   */
  void system(int s) {
    theSystem = s;
  }

  /**
   * Return the dipole system to which this parton originally belonged.
   */
  int origSystem() const {
    return theOrigSystem;
  }

  /**
   * Return the dipole system to which this parton originally belonged.
   */
  void origSystem(int s) {
    theOrigSystem = theSystem = s;
  }

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
   * If this was an original parton, this points to the corresponding
   * Particle object.
   */
  tPPtr theOrig;

  /**
   * The Emission which produced this parton, or null if original.
   */
  cEmPtr theEmission;

  /**
   * The final ThePEG::Particle produced for this parton with the
   * produceParticle() function.
   */
  PPtr theParticle;

  /**
   * The corresponding ParticleData object.
   */
  tcPDPtr theDataPtr;

  /**
   * The momentum of this Parton.
   */
  Lorentz5Momentum theMomentum;

  /**
   * The production vertex position of this parton.
   */
  LorentzPoint theVertex;

  /**
   * The original incoming anti-colour line corresponding to this
   * parton.
   */
  tColinePtr theICol;

  /**
   * The original outgoing colour line corresponding to this
   * parton.
   */
  tColinePtr theOCol;

  /**
   * Set true if this parton need to be notified about other emission
   * which may influence its state.
   */
  bool doNeedNotification;

  /**
   * Set true if this parton is present in the final state.
   */
  bool inFinalState;

  /**
   * Is true if this is a derived complex class.
   */
  bool isSpecial;

  /**
   * The dipole system to which this parton currently belongs.
   */
  int theSystem;

  /**
   * The dipole system to which this parton originally belonged.
   */
  int theOrigSystem;

  /**
   * If this parton has absorbed a gluon, it should also adopt that
   * gluons original parton.
   */
  vector<tPPtr> adoptedOriginals;

public:

  /**
   * Print out debugging information on std::cerr.
   */
  virtual void debugme() const;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  Parton & operator=(const Parton &);

};

}

#endif /* Ariadne5_Parton_H */
