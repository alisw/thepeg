// -*- C++ -*-
//
// ParticleData.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_ParticleData_H
#define ThePEG_ParticleData_H
// This is the declaration of the ParticleData class.

#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/PDT/PDT.h"
#include "ThePEG/PDT/PID.h"
#include "ThePEG/Vectors/LorentzVector.h"
#include "ThePEG/Vectors/ThreeVector.h"
#include "ThePEG/Interface/Interfaced.h"
#include "ThePEG/Utilities/Selector.h"
#include "ThePEG/PDT/WidthGenerator.h"
#include "ThePEG/PDT/MassGenerator.h"
#include "ThePEG/PDT/DecayMode.fh"
#include "ThePEG/Utilities/ClassTraits.h"
#include "ThePEG/Utilities/ClassDescription.h"

namespace ThePEG {

/**
 * ParticleData inherits from InterfacedBase and represents the
 * properties of a particle type. It is also able to produce instances
 * of this Particle type and, among other things, to decay them.
 *
 * @see \ref ParticleDataInterfaces "The interfaces"
 * defined for ParticleData.
 */
class ParticleData: public Interfaced {

public:

  /** The Repository is a friend. */
  friend class Repository;

  /** The EventGenerator is a friend. */
  friend class EventGenerator;

  /** DecayMode is a friend. */
  friend class DecayMode;

  /** A selector of DecayMode objects. */
  typedef Selector<tDMPtr> DecaySelector;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  ParticleData();

  /**
   * Destructor.
   */
  virtual ~ParticleData();
  //@}

  /** @name The Create methods are special interfaces for ParticleData
      classes. */
  //@{
  /**
   * Create a Particle which is its own anti-particle.
   */
  static PDPtr Create(PID newId, const string & newPDGName);

  /**
   * Create a particle - anti particle pair.
   */
  static PDPair Create(PID newId, const string & newPDGName, const string & newAntiPDGName);
  //@}

public:

  /** @name Acces to number and name. */
  //@{
  /**
   * Return the PDG id number.
   */
  long id() const { return theId; }

  /**
   * Return the generic PDG name. Note that this is not really
   * standardised.
   */
  const string & PDGName() const { return thePDGName; }

  /**
   * Return the generic PDG name. Note that this is not really
   * standardised.
   */
  const string & genericName() const { return thePDGName; }
  //@}

  /** @name Functions used for producing Particle instances. */
  //@{
  /**
   * Produce a particle specifying momentum.
   */
  PPtr produceParticle(const Lorentz5Momentum &) const;

  /**
   * Produce a particle specifying momentum.
   */
  PPtr produceParticle(const LorentzMomentum &) const;

  /**
   * Produce a particle specifying 4-momentum and a mass.
   */
  PPtr produceParticle(const LorentzMomentum &, Energy m) const;

  /**
   * Produce a particle specifying 3-momentum.
   */
  PPtr produceParticle(const Momentum3 & pp = Momentum3()) const;

  /**
   * Produce a particle specifying mass and 3-momentum.
   */
  PPtr produceParticle(Energy m, const Momentum3 & pp = Momentum3()) const;

  /**
   * Produce a particle specifying light-cone momentum components and
   * transverse momentum components.
   */
  PPtr produceParticle(Energy plus, Energy minus, Energy px, Energy py) const;

  /**
   * Generate a mass for an instance of this particle type.
   */
  Energy generateMass() const;

  /**
   * Generate a width for an instance of this particle type. Given a
   * \a mass of an instance of this particle type, calculate its width.
   */
  Energy generateWidth(Energy mass) const;

  /**
   * Generate a mass for an instance of this particle type. Given a \a
   * mass and a \a width of an instance of this particle type,
   * generate a life time.
   */
  Length generateLifeTime(Energy mass, Energy width) const;
  // Given a mass and a width of an instance of this particle type,
  // generate a life time.
  //@}

  /** @name Access the decay modes. */
  //@{
  /**
   * Return the nominal decay selector for this particle. Ie. the
   * decay modes weighted by their nominal branching ratios.
   */
  const DecaySelector & decaySelector() const { return theDecaySelector; }

  /**
   * Selects a decay mode randomly according to the branching
   * ratios. The nominal branching ratios may be changed for the
   * particular Particle instance \a p, iether by an assigned
   * WidthGenerator or the respective Decayers.
   */
  tDMPtr selectMode(Particle & p) const;

  /**
   * Access all the decay modes, including those which are
   * switched off, or have zero branching ratio
   */
  const DecaySet & decayModes() const { return theDecayModes; }
  //@}

  /**
   * Set the nominal mass
   */
  Energy mass(Energy);

  /**
   * Return the nominal mass.
   */
  Energy mass() const { return theMass; }

  /**
   * Return the mass to be used when evaluating hard process cross sections.
   */
  Energy hardProcessMass() const { 
    return 
      hardProcessMassSet ? theHardProcessMass : mass();
  }

  /**
   * Return the maximum possible mass of this particle type.
   */
  Energy massMax() const { return mass() + widthUpCut(); }

  /**
   * Return the minimum possible mass of this particle type.
   */
  Energy massMin() const { return max(mass() - widthLoCut(), ZERO); }

  /**
   * Return the constituent mass of this particle if relevant. This
   * version simply returns the nominal mass.
   */
  virtual Energy constituentMass() const { return mass(); }

  /**
   * Set the width.
   */
  Energy width(Energy);

  /**
   * Get the width. If no width is specified, it is calculated from
   * the lifetime.
   */
  Energy width() const {
    return theWidth >= ZERO ? theWidth :
      ( theCTau > Length() ? hbarc/theCTau :
	( theCTau == Length() ? Constants::MaxEnergy : ZERO ) );
  }

  /**
   * Set the width cut. Both upper and lower cut is set.
   */
  Energy widthCut(Energy wci) {
    widthUpCut(wci);
    return widthLoCut(wci);
  }

  /**
   * Get the width cut.
   */
  Energy widthCut() const { return max(widthUpCut(), widthLoCut()); }

  /**
   * Set the upper width cut.
   */
  Energy widthUpCut(Energy);

  /**
   * Get the upper width cut.
   */
  Energy widthUpCut() const {
    return theWidthUpCut >= ZERO? theWidthUpCut: Constants::MaxEnergy;
  }

  /**
   * Set the lower width cut.
   */
  Energy widthLoCut(Energy);

  /**
   * Get the lower width cut.
   */
  Energy widthLoCut() const {
    return theWidthLoCut >= ZERO? theWidthLoCut: Constants::MaxEnergy;
  }

  /**
   * Set the life time cTau.
   */
  Length cTau(Length);

  /**
   * Get the life time cTau cTau. If no life time is specified, it is
   * calculated from the width. If the width is also not specified,
   * the lifetime is assumed to be zero for ustable particles and
   * infinite for stable ones.
   */
  Length cTau() const {
    return theCTau > Length() ? theCTau :
      ( theWidth > ZERO ? hbarc/theWidth :
	( stable() ? Constants::MaxLength : Length() ) );
  }

  /**
   * Return the width to be used when evaluating hard process cross sections.
   */
  Energy hardProcessWidth() const { 
    return 
      hardProcessWidthSet ? theHardProcessWidth : width();
  }

  /**
   * Set the charge. The charge should be given
   * in units of e/3 using the PDT::Charge enum.
   */
  PDT::Charge iCharge(PDT::Charge);

  /**
   * Get the charge. The charge is returned in standard units and in
   * iCharge the charge is returned in units of e/3.
   */
  Charge charge() const { return eplus*double(theCharge)/3.0; }

  /**
   * Get the charge. The charge is returned in units of e/3.
   */
  PDT::Charge iCharge() const { return theCharge; }

  /**
   * Return true if charged.
   */
  bool charged() const { return PDT::charged(theCharge); }

  /**
   * Return true if positively charged.
   */
  bool positive() const { return PDT::positive(theCharge); }

  /**
   * Return true if negatively charged.
   */
  bool negative() const { return PDT::negative(theCharge); }

  /**
   * Set the spin. The spin should be given as 2J+1 (in units of
   * hbar/2) using the PDT::Spin enum.
   */
  PDT::Spin iSpin(PDT::Spin);

  /**
   * Get the spin.The spin is returned in standard units.
   */
  AngularMomentum spin() const { return hbar_Planck*double(theSpin-1)*0.5; }

  /**
   * Get the spin. The spin is returned as 2J+1 in units of hbar/2.
   */
  PDT::Spin iSpin() const { return theSpin; }

  /**
   * Set the colour of the particle in units of PDT::Colour.
   */
  PDT::Colour iColour(PDT::Colour);

  /**
   * Get the colour of the particle in units of PDT::Colour.
   */
  PDT::Colour iColour() const { return theColour; }

  /**
   * Return true if coloured.
   */
  bool coloured() const { return PDT::coloured(iColour()); }

  /**
   * Return true if (\a anti) coloured or colour-octet.
   */
  bool hasColour(bool anti = false) const {
    return anti? hasAntiColour():
      ( iColour() == PDT::Colour3 || iColour() == PDT::Colour6 ||
	iColour() == PDT::Colour8 );
  }

  /**
   * Return true if anti coloured or colour-octet.
   */
  bool hasAntiColour() const {
    return iColour() == PDT::Colour3bar || iColour() == PDT::Colour6bar || 
      iColour() == PDT::Colour8;
  }

  /**
   * Specify if particle is to be considered stable according to \a
   * stab.
   */
  void stable(bool stab);

  /**
   * Return true if particle is to be considered stable. If the decay
   * table is empty the function always returns true, even if the
   * member variable is false.
   */
  bool stable() const { return isStable; }

  /**
   * Get the pointer to the corresponding anti partner.
   */
  tPDPtr CC() const { return theAntiPartner; }

  /**
   * Specify if the anti partner chould be changed automatically when
   * this object is changed according to \a sync.
   */
  void synchronized(bool sync);

  /**
   * Return true if the anti partner chould be changed automatically
   * when this object is changed.
   */
  bool synchronized() const { return syncAnti; }

  /**
   * If there is an anti-partner, update this object to have correct
   * anti-properties.
   */
  void synchronize();

  /**
   * Set the mass generator object.
   */
  void massGenerator(tMassGenPtr);

  /**
   * Get the mass generator object.
   */
  tMassGenPtr massGenerator() const { return theMassGenerator; }

  /**
   * Set the width generator object.
   */
  void widthGenerator(tWidthGeneratorPtr);

  /**
   * Get the width generator object.
   */
  tWidthGeneratorPtr widthGenerator() const { return theWidthGenerator; }

  /**
   * Specify if the branching ratio of the Particle instances should vary with their
   * masses. 
   */
  void variableRatio(bool varRatio);

  /**
   * Return true if the branching ratio should vary with the mass of the Particle
   * instance.
   */
  bool variableRatio() const { return theVariableRatio; }

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

  static void Init();

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}

  /**
   * Special clone function used by the Repository. Also makes copies
   * the decay modes and the anti-partner if it exists and if
   * synchronized() is true.
   */
  virtual PDPtr pdclone() const;

  /**
   * Protected constructor only to be used by subclasses or by the
   * Create method.
   */
  ParticleData(PID newId, const string & newPDGName);

  /**
   * Read setup info from a standard stream. The following information
   * must be supplied in a white-space separated list: PDG number,
   * generic name, default mass (GeV), default width (GeV), width cut
   * (GeV), the lifetime ctau (mm), the charge, the colour, the spin,
   * stable (true) or not (false). Note that if a minus sign is given
   * instead of a generic name, the name of the object will be used
   * instead.
   */
  virtual void readSetup(istream & is);

  /**
   * Used by subclasses or by the Create method to setup
   * anti-relationship.
   */
  static void antiSetup(const PDPair & pap);


protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Check sanity of the object during the setup phase.
   */
  virtual void doupdate();

  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

  /**
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given
   * pointer.
   */
  virtual void rebind(const TranslationMap & trans)
   ;

  /**
   * Return a vector of all pointers to Interfaced objects used in this
   * object.
   * @return a vector of pointers.
   */
  virtual IVector getReferences();

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();
  //@}

protected:

  /**
   * Add a decay mode for this particle.
   */
  void addDecayMode(tDMPtr);

  /**
   * Remove a decay mode for this particle.
   */
  void removeDecayMode(tDMPtr);

private:

  /**
   * Id number according to the STDHEP/PDG standard.
   */
  PID theId;

  /**
   * Name and Id number according to the STDHEP/PDG standard.
   */
  string thePDGName;

  /**
   * Nominal mass.
   */
  Energy theMass;

  /**
   * Width.
   */
  Energy theWidth;

  /**
   * The mass to be used when evaluating hard process cross sections.
   */
  Energy theHardProcessMass;

  /**
   * True, if a hard process mass has been set.
   */
  bool hardProcessMassSet;

  /**
   * The width to be used when evaluating hard process cross sections.
   */
  Energy theHardProcessWidth;

  /**
   * True, if a hard process width has been set.
   */
  bool hardProcessWidthSet;

  /**
   * Upper width cut.
   */
  Energy theWidthUpCut;

  /**
   * Lower width cut.
   */
  Energy theWidthLoCut;

  /**
   * Lifetime.
   */
  Length theCTau;

  /**
   * Three times the charge.
   */
  PDT::Charge theCharge;

  /**
   * 2 times the spin plus one.
   */
  PDT::Spin theSpin;

  /**
   * The colour for this particle.
   */
  PDT::Colour theColour;

  /**
   * A pointer to an object capable to generate a mass for a particle
   * of this type.
   */
  MassGenPtr theMassGenerator;

  /**
   * True if the particle is considered stable.
   */
  bool isStable;

  /**
   * A selector of decay modes weighted by the nominal branching
   * ratios.
   */
  DecaySelector theDecaySelector;

  /**
   * The set of all decay modes.
   */
  DecaySet theDecayModes;

  /**
   * A pointer to an object capable to generate the branching
   * fractions for different decay modes for this particle type. The
   * object will be asked to generate branching fractions every time
   * the ParticleData object it updated and will modify the branching
   * fractions for every particle instance if variableRatio is true.
   */
  WidthGeneratorPtr theWidthGenerator;

  /**
   * Determine whether the branching fractions are allowed to change
   * on a particle-by-particle basis.
   */
  bool theVariableRatio;

  /**
   * Pointer to the object corresponding to the antiparticle. Set to
   * null if it is its own antiparticle.
   */
  tPDPtr theAntiPartner;

  /**
   * If syncAnti is true all changes to this object will be transfered
   * to the antiParticle.
   */
  bool syncAnti;

  /**
   * Helper variable to keep track of the default mass.
   */
  Energy theDefMass;

  /**
   * Helper variable to keep track of the default width.
   */
  Energy theDefWidth;

  /**
   * Helper variable to keep track of the default width cut.
   */
  Energy theDefCut;

  /**
   * Helper variable to keep track of the default lifetime.
   */
  Length theDefCTau;

  /**
   * Helper variable to keep track of the default charge.
   */
  PDT::Charge theDefCharge;

  /**
   * Helper variable to keep track of the default spin.
   */
  PDT::Spin theDefSpin;

  /**
   * Helper variable to keep track of the default colour.
   */
  PDT::Colour theDefColour;

  /**
   * Utility function for the interface.
   */
  void setMass(Energy);

  /**
   * Utility function for the interface.
   */
  void setHardProcessMass(Energy);

  /**
   * Reset the hard process mass
   */
  string doUnsetHardProcessMass(string);

  /**
   * Adjust the nominal mass to the hard process mass if a reshuffling
   * is not desirable.
   */
  string doAdjustNominalMass(string);

  /**
   * Utility function for the interface.
   */
  Energy defMass() const;

  /**
   * Utility function for the interface.
   */
  void setWidth(Energy);

  /**
   * Utility function for the interface.
   */
  void setHardProcessWidth(Energy);

  /**
   * Reset the hard process mass
   */
  string doUnsetHardProcessWidth(string);

  /**
   * Utility function for the interface.
   */
  Energy getWidth() const;

  /**
   * Utility function for the interface.
   */
  Energy defWidth() const;

  /**
   * Utility function for the interface.
   */
  void setCut(Energy);

  /**
   * Utility function for the interface.
   */
  Energy getCut() const;

  /**
   * Utility function for the interface.
   */
  Energy defCut() const;

  /**
   * Utility function for the interface.
   */
  void setUpCut(Energy);

  /**
   * Utility function for the interface.
   */
  Energy getUpCut() const;

  /**
   * Utility function for the interface.
   */
  void setLoCut(Energy);

  /**
   * Utility function for the interface.
   */
  Energy getLoCut() const;

  /**
   * Utility function for the interface.
   */
  void setCTau(Length);

  /**
   * Utility function for the interface.
   */
  Length getCTau() const;

  /**
   * Utility function for the interface.
   */
  Length defCTau() const;

  /**
   * Utility function for the interface.
   */
  void setStable(long);

  /**
   * Utility function for the interface.
   */
  long getStable() const;

  /**
   * Utility function for the interface.
   */
  void setSync(long);

  /**
   * Utility function for the interface.
   */
  long getSync() const;

  /**
   * Utility function for the interface.
   */
  void setVariableRatio(long);

  /**
   * Utility function for the interface.
   */
  long getVariableRatio() const;

  /**
   * Utility function for the interface.
   */
  string doSync(string);

  /**
   * Utility function for the interface.
   */
  void setMassGenerator(MassGenPtr);

  /**
   * Utility function for the interface.
   */
  void setWidthGenerator(WidthGeneratorPtr);

  /**
   * Utility function for the interface.
   */
  void setCharge(int);

  /**
   * Utility function for the interface.
   */
  string ssetCharge(string);

  /**
   * Utility function for the interface.
   */
  int getCharge() const;

  /**
   * Utility function for the interface.
   */
  int defCharge() const;

  /**
   * Utility function for the interface.
   */
  void setSpin(int);

  /**
   * Utility function for the interface.
   */
  int getSpin() const;

  /**
   * Utility function for the interface.
   */
  int defSpin() const;

  /**
   * Utility function for the interface.
   */
  void setColour(long);

  /**
   * Utility function for the interface.
   */
  long getColour() const;

  /**
   * Utility function for the interface.
   */
  long defColour() const;

  /**
   * Utility function for the interface.
   */
  void insDecayModes(DMPtr dm, int);

  /**
   * Utility function for the interface.
   */
  void delDecayModes(int i);

  /**
   * Utility function for the interface.
   */
  vector<DMPtr> getDecayModes() const;

  /**
   * Utility function for the interface.
   */
  string doSelectDecayModes(string);

  /**
   * Utility function for the interface.
   */
  string doPrintDecayModes(string);

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<ParticleData> initParticleData;

};

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the base classes
 *  of ParticleData. */
template <>
struct BaseClassTrait<ParticleData,1>: public ClassTraitsType {
  /** Typedef of the first base class of ParticleData. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of the
 *  ParticleData class. */
template <>
struct ClassTraits<ParticleData>: public ClassTraitsBase<ParticleData> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::ParticleData"; }
};

/** @endcond */

}

#endif /* ThePEG_ParticleData_H */
