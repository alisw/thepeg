// -*- C++ -*-
//
// MatcherBase.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_MatcherBase_H
#define ThePEG_MatcherBase_H
// This is the declaration of the MatcherBase class.


#include "ParticleData.h"
#include "ThePEG/EventRecord/Particle.h"

namespace ThePEG {

/**
 * MatcherBase is an abstract base class to be used for objects
 * representing groups of ParticleData objects. Concrete
 * implementations will typically use the templated Matcher class for
 * easy building of a full sub-class.
 *
 * @see ParticleData
 * @see Matcher
 * 
 */
class MatcherBase: public Interfaced {

public:

  /** Repository needs to be a friend. */
  friend class Repository;

  /**
   * Convenient typedef.
   */
  typedef set<tPDPtr> tPDSet;

  /**
   * Convenient typedef.
   */
  typedef set<tPMPtr> tPMSet;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  MatcherBase();

  /**
   * Copy-constructor.
   */
  MatcherBase(const MatcherBase &);

  /**
   * Destructor.
   */
  virtual ~MatcherBase();
  //@}

public:

  /** @name Virtual functions to be overridden by sub-classes. */
  //@{
  /**
   * Check if a particle type meets the criteria.
   */
  virtual bool check(const ParticleData &) const = 0;

  /**
   * Specialized clone method for MatcherBase used by the
   * Repository. A sub class must make sure that also the MatcherBase
   * object corresponding to the complex conjugate of this is cloned.
   */
  virtual PMPtr pmclone() const = 0;
  //@}

  /** @name Check if something is matched. */
  //@{
  /**
   * Check if a Particle meets the criteria.
   */
  bool checkp(const Particle & p) const { return check(p.data()); }

  /**
   * Check if a given particle type belongs to the set of
   * matches. This function looks for the same ParticleData object in
   * the set of all particles matched by this matcher. May be quicker
   * than to go through the check proceedure.
   */
  bool matches(const ParticleData & pd) const {
    return member(matchingParticles, PDPtr(const_cast<ParticleData *>(&pd)));
  }


  /**
   * Check if a given particle belongs to the set of matches. This
   * function looks for the corresponding ParticleData object in the
   * set of all particles matched by this matcher. May be quicker than
   * to go through the check proceedure.
   */
  bool matches(const Particle & p) const { return matches(p.data()); }
  
  /**
   * Check if a given particle matcher belongs to the set of
   * matches. This function looks for the same MatcherBase object in
   * the set of all matchers matched by this matcher.
   */
  bool matches(const MatcherBase & pm) const {
    return member(matchingMatchers, PMPtr(const_cast<MatcherBase *>(&pm)));
  }
  //@}

  /** @name Access the sets of matching particles and matchers. */
  //@{
  /**
   * Access to the set of matching particles.
   */
  const tPDSet & particles() const { return matchingParticles; }
  /**
   * Access to the set of matching matchers.
   */
  const tPMSet & matchers() const { return matchingMatchers; }
  //@}

  /** @name Access common properties of all matched particles. */
  //@{
  /**
   * Returns the minimum mass of the matching particles.
   */
  Energy minMass() const { return theMinMass; }

  /**
   * Returns the maximum mass of the matching particles.
   */
  Energy maxMass() const { return theMaxMass; }

  /**
   * Returns the common mass of the matching particles. If all matching
   * particles do not have exactly the same mass, -1.0 GeV is returned.
   */
  Energy mass() const { return commonMass; }

  /**
   * Returns the common width of the matching particles. If all matching
   * particles do not have exactly the same width, -1.0 GeV is returned.
   */
  Energy width() const { return commonWidth; }

  /**
   * Returns the common decay length of the matching particles. If all
   * matching particles do not have exactly the same decay length -1.0
   * mm is returned.
   */
  Length cTau() const { return commonCTau; }

  /**
   * Return common charge. If all matching particles have the same
   * charge the common charge is returned. Otherwise if all are
   * positive (negative), PDT::Positive (PDT::Negative) is
   * returned. Otherwise if all are charged, PDT::Charged is
   * returned. Otherwise PDT::ChargeUndefined is returned.
   */
  PDT::Charge iCharge() const { return commonCharge; }

  /**
   * Are the particles charged? If all matching particles are charged, return
   * true, otherwise false.
   */
  bool charged() const { return PDT::charged(commonCharge); }

  /**
   * Are the particles positively charged? If all matching particles
   * are positively charged, return true, otherwise false.
   */
  bool positive() const { return PDT::positive(commonCharge); }

  /**
   * Are the particles negatively charged? If all matching particles
   * are negatively charged, return true, otherwise false.
   */
  bool negative() const { return PDT::negative(commonCharge); }

  /**
   * Return common spin. If all matching particles have the same spin,
   * the common spin is returned. Otherwise PDT::SpinUndefined is
   * returned.
   */
  PDT::Spin iSpin() const { return commonSpin; }

  /**
   * If all matching particles have the same colour, the common colour
   * is returned. Otherwise if all are coloured, PDT::Coloured is
   * returned. Otherwise PDT::ColourUndefined is returned.
   */
  PDT::Colour iColour() const { return commonColour; }

  /**
   * Are the particles coloured? If all matching particles are
   * coloured, return true, otherwise false.
   */
  bool coloured() const { return PDT::coloured(commonColour); }

  /**
   * Are the particles stable? Returns (0)1 if all matching particles
   * are (un)stable. Otherwise -1 is returned.
   */
  int stable() const { return commonStable; }
  //@}

  /**
   * Get the matcher object matching the antiparticles of this. If
   * no-one exists null is returned.
   */
  tPMPtr CC() const { return theAntiPartner; }

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
  //@}

protected:

  /**
   * Add a particle to the set of matching particles if it meets the
   * criteria.
   */
  void addPIfMatch(tPDPtr);

  /**
   * Add a particle matcher to the set of matching matchers if it
   * meets the criteria.
   */
  void addMIfMatch(tPMPtr);

  /**
   * Add a number of particles to the set of matching particles if
   * they meets the criteria.
   */
  template <typename Iterator>
  void addPIfMatch(Iterator first, Iterator last) {
    for ( ; first != last; ++first ) addPIfMatch(*first);
  }

  /**
   * Add a number of particles to the set of matching particles if
   * they meets the criteria.
   */
  template <typename Cont>
  void addPIfMatchFrom(const Cont & c) {
    addPIfMatch(c.begin(), c.end());
  }

  /**
   * Add a number of particle matchers to the set of matching
   * matchers if they meets the criteria.
   */
  template <typename Iterator>
  void addMIfMatch(Iterator first, Iterator last) {
    for ( ; first != last; ++first ) addMIfMatch(*first);
  }

  /**
   * Add a number of particle matchers to the set of matching
   * matchers if they meets the criteria.
   */
  template <typename Cont>
  void addMIfMatchFrom(const Cont & c) {
    addMIfMatch(c.begin(), c.end());
  }

  /**
   * Clear information about matching particles and matchers.
   */
  void clear();

  /**
   * Set antipartner.
   */
  static void setCC(tPMPtr pm, tPMPtr apm) {
    pm->theAntiPartner = apm;
    apm->theAntiPartner = pm;
  }

private:

  /**
   * The set of particle data objects matched by this matcher.
   */
  tPDSet matchingParticles;

  /**
   * A set of matchers which matches a subset of this matcher.
   */
  tPMSet matchingMatchers;

  /**
   * The maximum mass of all matching particles.
   */
  Energy theMaxMass;

  /**
   * The minimum mass of all matching particles.
   */
  Energy theMinMass;

  /**
   * The common mass of all matching particles.
   */
  Energy commonMass;

  /**
   * The common width of all matching particles.
   */
  Energy commonWidth;

  /**
   * The common decay length of all matching particles.
   */
  Length commonCTau;

  /**
   * The common charge of all matching particles.
   */
  PDT::Charge commonCharge;

  /**
   * The common spin of all matching particles.
   */
  PDT::Spin commonSpin;

  /**
   * The common colour of all matching particles.
   */
  PDT::Colour commonColour;

  /**
   * The common stability of all matching particles.
   */
  int commonStable;

  /**
   * Pointer to a matcher object which matches all anti particles
   * which are matched by this matcher.
   */
  tPMPtr theAntiPartner;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class with persistent data.
   */
  static AbstractClassDescription<MatcherBase> initMatcherBase;

  /**
   *  Private and non-existent assignment operator.
   */
  MatcherBase & operator=(const MatcherBase &);

};

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the base classes
 *  of MatcherBase. */
template <>
struct BaseClassTrait<MatcherBase,1>: public ClassTraitsType {
  /** Typedef of the first base class of MatcherBase. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of the
 *  MatcherBase class. */
template <>
struct ClassTraits<MatcherBase>:
    public ClassTraitsBase<MatcherBase> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::MatcherBase"; }
};

/** @endcond */

}

#endif /* ThePEG_MatcherBase_H */
