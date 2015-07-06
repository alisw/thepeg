// -*- C++ -*-
//
// DecayMode.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_DecayMode_H
#define ThePEG_DecayMode_H
// This is the declaration of the DecayMode class.

#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/Interface/Interfaced.h"
#include "DecayMode.fh"
#include "MatcherBase.h"
#include "Decayer.h"

namespace ThePEG {

ThePEG_DECLARE_MULTISET(tPDPtr,ParticleMSet);
ThePEG_DECLARE_MULTISET(tPMPtr,MatcherMSet);
ThePEG_DECLARE_MULTISET(tDMPtr,ModeMSet);

/**
 * The DecayMode class describes a decay channel of a particle. In its
 * simplest form it contains simply a parent ParticleData object and a
 * list of decay products, but it can also specify a set of
 * MatcherBase objects each representing one of a set of possible
 * decay products. A matcher can also be specified to represents an
 * unlimited set of decay products. Decay chains can be represented by
 * specifying other decay channels where the parents are taken to be
 * intermediate resonances. It is also possible to specify the absence
 * of intermediate resonances.
 *
 * Each decay mode can be uniquely described by a character string on
 * the form
 * <code>{decaying-particle-name}->{decay-product-specifier}[,{decay-product-specifier},...];</code>
 * where no spaces are allowed anywhere. The decaying-particle-name
 * should be a path to a particle in the Repository or a
 * PDG-standardized particle name (see ParticleData::PDGName()) and
 * the decay-product-specifier can be one of the following:
 *
 * <ul>
 *
 * <li> a path to a particle in the Repository or a PDG-standardized
 * particle name if it is a specific decay product,
 *
 * <li> a question mark followed by the name of a particle matcher
 * object in the Repository if representing one of several different
 * alternative decay products,
 *
 * <li> a star followed by the name of a particle matcher object in
 * the Repository if representing any number of several different
 * alternative decay products (note that only one of these
 * <i>wild-card</i> matchers may be specified, if several are given
 * only the last one will be taken into account),
 *
 * <li> a whole decay mode string enclosed in square brackets if
 * representing a resonance decay product with a specified decay mode,
 *
 * <li> an exclamation mark followed by a path to a particle in the
 * Repository or a PDG-standardized particle name representing the
 * exclusion of an intermediate resonance decaying into the other
 * specified decay products.
 *
 * <li> two paths to particles in the Repository or PDG-standardized
 * particle names with an equal sign between, indicating two coloured
 * particles which are in a colour-singlet state.
 *
 * </ul>
 *
 * Note that the order of the specified decay products will be
 * preserved when the corresponding particles are produced in a decay.
 *
 * The possibility of specifying matchers as decay products means that
 * one decay mode may overlap with another one. When an EventGenerator
 * is initialized all decay modes are checked so that a given decay
 * mode afterwards will report if there are other modes which are
 * matched by its matchers through the overlap() function.
 *
 * @see \ref DecayModeInterfaces "The interfaces"
 * defined for DecayMode.
 * @see ParticleData
 * @see MatcherBase
 */
class DecayMode: public Interfaced {

public:

  /** ParticleData is a friend. */
  friend class ParticleData;

  /** RemnantData is a friend. */
  friend class RemnantData;

  /** The EventGenerator is a friend. */
  friend class EventGenerator;

public:

  /** A vector of DecayMode pointers. */
  typedef vector<tcDMPtr> ModeVector;
  /** A vector of pairs of ParticleData pointers. */
  typedef vector<tPDPair> LinkVector;

public:

  /**
   * Create a decay mode from a given tag. This function is used
   * directly by the Repository. If name of the decaying particle is a
   * valid path to a particle object, the decaymode will be added to
   * that particle, otherwise it will be added to the default particle
   * of that name.
   */
  static DMPtr constructDecayMode(string & tag, vector<DMPtr> * save = 0);

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  DecayMode();

  /**
   * Copy-constructor.
   */
  DecayMode(const DecayMode &);

  /**
   * Destructor.
   */
  ~DecayMode();
  //@}

  /**
   * Return a clone of this decay mode with \a pd as the decaying
   * particle.
   */
  virtual DMPtr clone(tPDPtr pd) const;

public:

  /**
   * Return the tag for this decay mode. This string is a unique
   * identifier for this decay mode.
   */
  const string & tag() const {
    return theTag.size() ? theTag : ( theTag = makeTag() );
  }

  /**
   * Get a pointer to the particle data object corresponding to
   * the decaying particle.
   */
  tcPDPtr parent() const { return theParent; }

  /**
   * The set of identified decay products.
   */
  const ParticleMSet & products() const { return theProducts; }

  /**
   * The set of identified decay products in the order they were specified.
   */
  const tPDVector & orderedProducts() const { return theOrderedProducts; }

  /**
   * Produce particles corresponding to the identified decay
   * products. They will be orderd in the same order they were
   * sspecified.
   */
  PVector produceProducts() const;

  /**
   * The set of identified resonance products with specified decay
   * modes
   */
  const ModeMSet & cascadeProducts() const { return theCascadeProducts; }

  /**
   * The set of matchers each corresponding to one decay product.
   */
  const MatcherMSet & productMatchers() const { return theMatchers; }

  /**
   * The pointer to a matcher corresponding to any number of decay
   * products
   */
  tPMPtr wildProductMatcher() const { return theWildMatcher; }

  /**
   * The set particles corresponding to excluded intermediate
   * resonances.
   */
  const ParticleMSet & excluded() const { return theExcluded; }

  /**
   * Return the branching ratio to be used.
   */
  double brat() const;

  /**
   * Calculate the branching ratio for a particular particle instance.
   */
  double brat(const Particle &) const;

  /**
   * Get the decayer assigned to this mode.
   */
  tDecayerPtr decayer() const { return theDecayer; }

  /**
   * Check if another decay mode is included in this one.
   */
  bool includes(const DecayMode &) const;

  /**
   * Return a pointer to the corresponding decaymode for the
   * antiparticle decay.
   */
  tDMPtr CC() const { return theAntiPartner; }

  /**
   * Check if another decay mode has the same final state as this
   * one.
   */
  bool operator == (const DecayMode & d) const { 
    return tag() == d.tag() ;
  }

  /**
   * Return a vector of pairs of decay products which are linked
   * together (e.g. colourless q-qbar pairs).
   */
  const LinkVector & links() const { return theLinks; }

  /**
   * Return the list of overlapping decay modes.
   */
  const ModeVector & overlap() const { return theOverlap; }

  /**
   * Modify this mode to have properties corresponding to its anti-partner.
   */
  void synchronize();

  /**
   *  Check whether this decay mode is switched on
   */
  bool on() const { return isOn; }

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
   * Standard Init function used to initialize the interface.
   */
  static void Init();

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Check sanity of the object during the setup phase.
   */
  virtual void doupdate();

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
  //@}

protected:

  /**
   * Set a pointer to the particle data object corresponding to
   * the decaying particle.
   */
  void parent(tPDPtr pd) { theParent = pd; }

  /**
   * Set the branching ratio to be used.
   */
  void brat(double);

  /**
   * Switch on this decay mode.
   */
  void switchOn();

  /**
   * Switch off this decay mode.
   */
  void switchOff();

  /**
   * Set the decayer. The set method returns false if the decayer
   * does not claim to be able to handle the decay.
   */
  void decayer(tDecayerPtr);

  /**
   * Add identified decay products.
   */
  void addProduct(tPDPtr);

  /**
   * Add a pair of decay products which are linked together
   * (e.g. colourless q-qbar pairs).
   */
  void addLink(tPDPtr a, tPDPtr b);

  /**
   * Add identified resonant product with specified decay mode.
   */
  void addCascadeProduct(tDMPtr);

  /**
   * Add a mathcer corresponding to one decay product.
   */
  void addProductMatcher(tPMPtr);

  /**
   * Add a matcher corresponding to any number of decay products.
   */
  void setWildMatcher(tPMPtr);

  /**
   * Add a particle corresponding to an excluded intermediate
   * resonance.
   */
  void addExcluded(tPDPtr);

  /**
   * Protected creation and clone methods.
   */
  static DMPtr Create(tPDPtr newParent, double newBrat = 0.0,
		      bool newOn = false);
  /**
   * Protected constructor.
   */
  DecayMode(tPDPtr newParticle, double newBrat, bool newOn);

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
   * Protected special clone function.
   */
  DMPtr dmclone() const;

  /**
   * Read setup info from a standard stream used by the
   * Repository. The following information must be supplied in a
   * white-space separated list: the branching ratio, on or off (true
   * or false), and the name of a Decayer.
   */
  virtual void readSetup(istream & is);

  /**
   * The set of identified decay products.
   */
  ParticleMSet & products() { return theProducts; }

  /**
   * The set of identified resonant products with specified decay
   * modes
   */
  ModeMSet & cascadeProducts() { return theCascadeProducts; }

  /**
   * The set of matchers each corresponding to one decay product.
   */
  MatcherMSet & productMatchers() { return theMatchers; }

  /**
   * The pointer to a matcher corresponding to any number of decay
   * products
   */
  tPMPtr & wildProductMatcher() { return theWildMatcher; }

  /**
   * The set particles corresponding to excluded intermediate
   * resonances.
   */
  ParticleMSet & excluded() { return theExcluded; }

  /**
   *  Set the pointer to the corresponding decaymode for the
   * antiparticle decay.
   */
  void CC(tDMPtr cc) {theAntiPartner = cc;}

private:

  /**
   * Add a decay mode to the list of overlapping modes if included.
   */
  bool addOverlap(tcDMPtr);

  /**
   * Remove all decay modes from the list of overlapping modes.
   */
  void resetOverlap();

  /**
   * Check if two sets of particles have equivalent types.
   */
  bool compareId(const ParticleMSet &, const ParticleMSet &) const;

  /**
   * Check if a particle set contains a given particle ID.
   */
  ParticleMSet::const_iterator findId(const ParticleMSet &,
				      const ParticleData &) const;
  
  /**
   * Use the members in this decay channel and generate the
   * corresponding tag.
   */
  string makeTag() const;

  /**
   * Delete the tag (it will be regenerated later if asked for).
   */
  void resetTag() {
    theTag = "";
    if ( CC() ) CC()->theTag = "";
  }


private:

  /**
   * Utility function for the interface.
   */
  void setOn(long);

  /**
   * Utility function for the interface.
   */
  long getOn() const;

  /**
   * Utility function for the interface.
   */
  void setDecayer(DecayerPtr);

private:

  /**
   * The tag.
   */
  mutable string theTag;

  /**
   * The branching ratio.
   */
  double theBrat;

  /**
   * True if this mode is switched on.
   */
  bool isOn;

  /**
   * Pointer to a particle data object corresponding to the decaying
   * particle.
   */
  tPDPtr theParent;

  /**
   * The set of specified decay particles.
   */
  ParticleMSet theProducts;

  /**
   * The set of specified decay particles in the order they was specified.
   */
  tPDVector theOrderedProducts;

  /**
   * The set of matching decay channels corresponding to a specified
   * with a specified subsequent decay mode.
   */
  ModeMSet theCascadeProducts;

  /**
   * The set of matching decay products. Each of the matchers
   * correspond to one particle.
   */
  MatcherMSet theMatchers;

  /**
   * A particle matcher which corresponds to zero or more particles.
   */
  tPMPtr theWildMatcher;

  /**
   * A set of particles which are not allowed as intermediate
   * resonances.
   */
  ParticleMSet theExcluded;

  /**
   * A list of decay modes which are included in this one.
   */
  ModeVector theOverlap;

  /**
   * The decayer object responsible for performing the decay.
   */
  DecayerPtr theDecayer;

  /**
   * The corresponding decay mode of the anti particle.
   */
  tDMPtr theAntiPartner;

  /**
   * The vector of pairs of decay products which are linked together
   * (e.g. colourless q-qbar pairs).
   */
  LinkVector theLinks;

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<DecayMode> initDecayMode;

  /**
   *  Private and non-existent assignment operator.
   */
  DecayMode & operator=(const DecayMode &);

};

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of DecayMode. */
template <>
struct BaseClassTrait<DecayMode,1>: public ClassTraitsType {
  /** Typedef of the first base class of DecayMode. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of the
 *  DecayMode class. */
template <>
struct ClassTraits<DecayMode>:
    public ClassTraitsBase<DecayMode> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::DecayMode"; }
};

/** @endcond */

}

#endif /* ThePEG_DecayMode_H */
