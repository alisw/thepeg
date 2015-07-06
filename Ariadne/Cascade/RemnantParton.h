// -*- C++ -*-
#ifndef Ariadne5_RemnantParton_H
#define Ariadne5_RemnantParton_H
//
// This is the declaration of the RemnantParton class.
//

#include "Parton.h"
#include "RemnantParton.fh"
#include "Models/RemnantModel.fh"
#include "ThePEG/PDF/PartonBinInstance.h"
#include "ThePEG/PDF/PDF.h"

namespace Ariadne5 {

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
  RemnantParton();

  /**
   * The destructor.
   */
  inline virtual ~RemnantParton() {}
  //@}

public:

  /** @name Simple access functions. */
  //@{
  /**
   * The object responsible for modeling some aspects of the remnant.
   */
  inline const RemnantModel & model() const {
    return *theModel;
  }

  /**
   * Return the inverse extension of this remant.
   */
  inline Energy mu() const {
    return theMu;
  }

  /**
   * Return the dimensionality of the extension of this remnant.
   */
  inline double alpha() const {
    return theAlpha;
  }

  /**
   * Return the supression power of emission scales above the maximum.
   */
  inline double beta() const {
    return theBeta;
  }

  /**
   * The original incoming particle. May be null if no information was
   * available in the SubProcess.
   */
  tPPtr parent() const {
    return theParent;
  }

  /**
   * The original incoming particle type.
   */
  const ParticleData & parentData() const {
    return *theParentData;
  }

  /**
   * The incoming particle momentum.
   */
  const Lorentz5Momentum & parentMomentum() const {
    return theParentMomentum;
  }

  /**
   * The original extracted parton. May be null if no information was
   * available in the SubProcess.
   */
  tPPtr originalExtracted() const {
    return theOriginalExtracted;
  }

  /**
   * The extracted parton type.
   */
  const ParticleData & extractedData() const {
    return *theExtractedData;
  }

  /**
   * The momentum of the extracted parton.
   */
  const Lorentz5Momentum & extractedMomentum() const {
    return theExtractedMomentum;
  }

  /**
   * The current x-value.
   */
  double x() const {
    return theX;
  }

  /**
   * Get the current factorization scale.
   */
  Energy2 muF2() const {
    return theMuF2;
  }

  /**
   * Set a new x-value. If negative, use the information in
   * extractedMomentum() and the reference particle to calculate.
   */
  void x(double xin);

  /**
   * Set the current factorization scale.
   */
  void muF2(Energy2 Q2) {
    theMuF2 = Q2;
  }

  /**
   * For hard remnants, the energy fraction of the virtual boson in DIS.
   */
  double y() const {
    return theY;
  }

  /**
   * The final extracted parton as produced by produceParticle()
   * incase this is a soft remnant.
   */
  tPPtr extracted() const {
    return theExtracted;
  }

  /**
   * The PDF associated with the parent. May be empty for hard remnants.
   */
  const PDF & pdf() const {
    return thePDF;
  }

  /**
   * Return the PDF value (momentum density) for the currently
   * extracted parton, given a scale). Return zero if not a soft
   * remnant.
   */
  double xfx(Energy2 scale) const {
    return pdf().pdf()? pdf().xfx(theExtractedData, scale, x()): 0.0;
  }

  /**
   * Return the PDF value (momentum density) for the currently
   * extracted parton, given a a scale and momentum fraction, \a xx. Return
   * zero if not a soft remnant.
   */
  double xfx(Energy2 scale, double xx) const {
    return pdf().pdf()? pdf().xfx(theExtractedData, scale, xx): 0.0;
  }

  /**
   * Return the PDF value (momentum density) for the given \a
   * extracted parton, given a a scale and momentum fraction, \a xx. Return
   * zero if not a soft remnant.
   */
  double xfx(tcPDPtr extracted, Energy2 scale, double xx) const {
    return pdf().pdf()? pdf().xfx(extracted, scale, xx): 0.0;
  }

  /**
   * Return the ratio of pdfs when evolving backwards from the current
   * parton at the current x, to a \a newparton at x/\a z at the given
   * \a scale.
   *
   * @return 0 if no PDF object has been assigned, or a negative value
   * if the density for the current parton is zero.
   */
  double xfratio(tcPDPtr newparton, Energy2 scale, double z) const;

  /**
   * Return the effective mass for this parton. Optionally supply
   * anothe parton to be extracted.
   */
  Energy effectiveMass(tcPDPtr extracted = tcPDPtr()) const;

  /**
   * For an emission where this remnant would acquire a transverse
   * momentum recoil, use the RemnantModel in the controlling
   * AriadneHandler to calculate a weight related to the subsequent
   * emission of a recoil gluon. \a ph is the momentum of the hard
   * subsystem and \a pr is the momentum of the remnant after the
   * emission under consideration. If the emission is an initial state
   * g->q or q->g splitting, the new \a extracted parton should be given.
   */
  double recoilWeight(const LorentzMomentum & ph, const LorentzMomentum & pr,
		      tcPDPtr extracted = tcPDPtr()) const;

  /**
   * If this remnant has acquired a transverse momentum recoil in a
   * previous emission, calculate the weight related to the subsequent
   * emission of a recoil gluon.
   */
  double recoilWeight() const;

  /**
   * Return the soft suppression weight for an emission at scale \a
   * rho, taking a fraction x of the positive momentum of the incoming
   * particle.
   */
  double softSuppression(Energy rho, double xplus) const;


  /**
   * Calculate a weight of an emission generated from a dipole with
   * this remnant giving a parton with momentum \a pem. \a prem is the
   * momentum of the remnant after the emission. Veto if it is not
   * considered a final-state splitting, otherwise reweight if too
   * much momentum is taken from the remnant. Optionally an \a
   * extracted flavour may be supplied if the emission was an initial
   * state g->q or q->g splitting.
   */
  double reweightFS(Energy rho, LorentzMomentum pem,
		    LorentzMomentum prem, tcPDPtr extracted = tcPDPtr()) const;


  /**
   * Get the rotation to transform the momentum \a phold to \a
   * phnew. It is done in a way such that the momenta which are
   * opposite in direction from the remnant are left untouched.
   */
  LorentzRotation
  getHardTransform(const LorentzMomentum & phold,
		   const LorentzMomentum & phnew) const;
  

  /**
   * Return true if this is a hard remnant.
   */
  bool hard() const {
    return !pdf().pdf();
  }

  /**
   * Return the direction of this remnant. Negative means incoming
   * from the seconf beam, positive the first.
   */
  int dir() const {
    return theDirection;
  }

  /**
   * Set the direction of this remnant. Negative means incoming
   * from the seconf beam, positive the first.
   */
  void dir(int d) {
    theDirection = d;
  }

  /**
   * Set the transformation needed to go to the rest system of this
   * remnants incoming particle and the other incoming particle, with
   * this incoming particle along the positive z-axis.
   */
  const LorentzRotation & setBoost();

  /**
   * Return the transformation needed to go to the rest system of this
   * remnants incoming particle and the other incoming particle, with
   * this incoming particle along the positive z-axis.
   */
  const LorentzRotation & getBoost() const {
    return theBoost;
  }

  /**
   * Return the inverse of getBoost().
   */
  const LorentzRotation & invBoost() const {
    return theInvBoost;
  }

  /**
   * Return the squared transverse momentum of the given momentum in
   * the frame specified by getBoost().
   */
  Energy2 getPT2Kick(const LorentzMomentum & p) const {
    return (getBoost()*p).perp2();
  }

  /**
   * Return the squared transverse momentum of this remnant in the
   * frame specified by getBoost().
   */
  Energy2 getPT2Kick() const {
    return getPT2Kick(momentum());
  }

  /**
   * Produce a ThePEG::Particle corresponding to this parton. The
   * momentum of the produced particle is rotated with \a r w.r.t. the
   * parton. If this is a hard remnant, the corresponding particle is
   * produced, while if this is a soft remnant, the corresponding
   * extracted particle is produced.
   */
  virtual tPPtr produceParticle(const LorentzRotation & r = LorentzRotation());
  //@}

public:

  /** @name Setup functions. */
  //@{
  /**
   * Setup remnant from a PartonBinInstance.
   */
  void setup(const PartonBinInstance & pb, int indir);

  /**
   * Setup parent.
   */
  void setupParent(tPPtr p) {
    theParent = p;
    setupParent(p->dataPtr(), p->momentum());
  }

  /**
   * Setup parent.
   */
  void setupParent(tcPDPtr pd, const Lorentz5Momentum & p) {
    theParentData = pd;
    theParentMomentum = p;
  }

  /**
   * Set the particle type for this remnant. Is always taken to be the
   * anti-partner of the extracted parton.
   */
  void setData() {
    data(extractedData().CC()? tcPDPtr(extractedData().CC()): theExtractedData);
  }

  /**
   * Setup extracted parton.
   */
  void setupExtracted(tPPtr p) {
    theOriginalExtracted = p;
    setExtracted(p->dataPtr(), p->momentum());
  }

  /**
   * Set extracted parton.
   */
  void setExtracted(tcPDPtr pd, const Lorentz5Momentum & p) {
    theExtractedData = pd;
    setData();
    theExtractedMomentum = p;
    x(-1.0);
  }

  void setExtracted(tcPDPtr pd = tcPDPtr()) {
    if ( !pd ) pd = theExtractedData;
    setExtracted(pd, parentMomentum() - momentum());
  }

  virtual void setMomentum(const Lorentz5Momentum & p);

  void setMomentum(const Lorentz5Momentum & p, tcPDPtr pd) {
    momentum() = p;
    setExtracted(pd);
  }

  bool untouched() const {
    return originalExtracted() &&
      theExtractedData == originalExtracted()->dataPtr() &&
      theExtractedMomentum == originalExtracted()->momentum();
  }

  /**
   * Set parameters for hard remnant.
   */
  void setupHard(int indir);

  /**
   * Setup pdf.
   */
  void setupPDF(const PDF & f) {
    thePDF = f;
  }

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

protected:

  /**
   * The object modeling some aspects of remnants.
   */
  tcRemnantModelPtr theModel;

  /**
   * The inverse extension of this remant. To be set by sub-classes.
   */
  Energy theMu;

  /**
   * The dimensionality of the extension of this remnant. To be set by
   * sub-classes.
   */
  double theAlpha;

  /**
   * The supression power of emission scales above the maximum.
   */
  double theBeta;

private:

  /**
   * The original incoming particle. May be null if no information was
   * available in the SubProcess.
   */
  tPPtr theParent;

  /**
   * The original incoming particle type.
   */
  tcPDPtr theParentData;

  /**
   * The incoming particle momentum.
   */
  Lorentz5Momentum theParentMomentum;

  /**
   * The original extracted parton. May be null if no information was
   * available in the SubProcess.
   */
  tPPtr theOriginalExtracted;

  /**
   * The extracted parton type.
   */
  tcPDPtr theExtractedData;

  /**
   * The momentum of the extracted parton.
   */
  Lorentz5Momentum theExtractedMomentum;

  /**
   * The final extracted parton produced in produceParticle() if this
   * was a soft remnant.
   */
  PPtr theExtracted;

  /**
   * The PDF associated with the parent. May be empty for hard remnants.
   */
  PDF thePDF;
  
  /**
   * The current x-value.
   */
  double theX;

  /**
   * The current factorization scale.
   */
  Energy2 theMuF2;

  /**
   * For hard remnants, the energy fraction of the virtual boson in DIS.
   */
  double theY;

  /**
   * The direction of this remnant. Negative means incoming from the
   * seconf beam, positive the first.
   */
  int theDirection;

  /**
   * The transformation needed to go to the rest system of this
   * remnants incoming particle and the other incoming particle, with
   * this incoming particle along the positive z-axis.
   */
  LorentzRotation theBoost;

  /**
   * The inverse of theBoost.
   */
  LorentzRotation theInvBoost;

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
  RemnantParton & operator=(const RemnantParton &);

};

}

#endif /* Ariadne5_RemnantParton_H */
