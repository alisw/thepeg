// -*- C++ -*-
//
// SimpleFlavour.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_SimpleFlavour_H
#define THEPEG_SimpleFlavour_H
// This is the declaration of the SimpleFlavour class.

#include "ThePEG/Handlers/FlavourGenerator.h"
#include "ThePEG/Utilities/VSelector.h"
// #include "SimpleFlavour.fh"
// #include "SimpleFlavour.xh"

namespace ThePEG {

/**
 * SimpleFlavour is a simple class to generate hadrons given the quark
 * flavours. It implements a simplified version of the model of the
 * old fortran version of Pythia.
 *
 * @see \ref SimpleFlavourInterfaces "The interfaces"
 * defined for SimpleFlavour.
 */
class SimpleFlavour: public FlavourGenerator {

public:

  /** A map of <code>Selector</code>s. */
  typedef map<long, VSelector< pair<long,long> > > ProbabilityMap;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  SimpleFlavour();

  /**
   * Destructor.
   */
  virtual ~SimpleFlavour();
  //@}

public:

  /** @name Virtual functions mandated by the FlavourGenerator base class. */
  //@{
  /**
   * Generate a hadron from a quark. Given a quark(antiquark, diquark
   * or antidiquark), choose a quark-antiquark (or
   * antidiquark-diquark) pair. Return (first) a hadron formed by the
   * original quark and the antiquark together with (second) the
   * generated quark. Returns null pointers if the generation failed.
   * @param quark a quark, antiquark, diquark or antidiquark.
   * @return a pair of ParticleData pointers. The \a first is the
   * hadron produced and the \a second is the anti-partner of the
   * (anti-)(di-)quark generated to form the hadron.
   */
  virtual tcPDPair generateHadron(tcPDPtr quark) const;

  /**
   * Get hadron from flavours. Return a hadron with the flavour
   * content given by the (anti-)(di-)quarks in the argument. The
   * arguments are given as PDG codes.
   * @param iq1 the PDG code of the first flavour.
   * @param iq2 the PDG code of the second flavour.
   * @return the corresponding hadron type or null if none could be
   * generated.
   */
  virtual tcPDPtr getHadron(long iq1, long iq2) const;
  using FlavourGenerator::getHadron;

  /**
   * Return a baryon with the flavour content given by the
   * (anti)quarks in the argument.  The arguments are given as
   * particle data pointers.
   * @param q1 the PDG code of the first flavour.
   * @param q2 the PDG code of the second flavour.
   * @param q3 the PDG code of the third flavour.
   * @return the corresponding baryon type or null if none could be
   * generated.
   */
  virtual tcPDPtr getBaryon(long q1, long q2, long q3) const;
  using FlavourGenerator::getBaryon;

  /**
   * Generate a random quark flavour.
   */
  virtual long selectQuark() const;

  /**
   * Generate a random (di)quark flavour.
   */
  virtual long selectFlavour() const;
  //@}

public:

  /** @name Access the parameters controlling the generation. */
  //@{
  /**
   * Return the suppression factor of strange quarks w.r.t. u and d.
   */
  double sSup() const { return theSSup; }

  /**
   * Return the suppression factor for di-quarks w.r.t. quarks
   */
  double diSup() const { return theDiSup; }

  /**
   * Return the suppression of spin-1 di-quarks w.r.t. spin-0 ones;
   */
  double di1Sup() const { return theDi1Sup; }

  /**
   * Return the suppression of strange di-quarks w.r.t. u and d ones
   * in addition to the standard strangness suppression of quarks.
   */
  double diSSup() const  { return theDiSSup; }

  /**
   * Return the extra suppression of eta's
   */
  double etaSup() const  { return theEtaSup; }

  /**
   * Return the extra suppression of ets-prime's
   */
  double etaPSup() const  { return theEtaPSup; }

  /**
   * Return the extra suppression for baryons of the spin 3/2
   * decuplet.
   */
  double baryon10Sup() const { return theBaryon10Sup; }

  /**
   * Return the probability that light (u/d) mesons has spin 1;
   */
  double pSpin1() const { return thePSpin1; }

  /**
   * Return the probability that strange mesons has spin 1;
   */
  double pSpinS1() const { return thePSpinS1; }

  /**
   * Return the probability that charmed and heavier mesons has spin
   * 1;
   */
  double pSpinC1() const { return thePSpinC1; }
  //@}

protected:

  /**
   * Calculate the probabilities for generateHadron for the given
   * flavour.
   */
  virtual void setProbabilities(long iq) const;

  /**
   * Return the probability that the given quark flavours end up in a
   * vector meson rather than in a pseudo scalar meson.
   */
  virtual double vectorMesonProbability(long iq1, long iq2) const;

  /**
   * Return the probability that the given quark and diquark flavours
   * end up in a spin 3/2 decuplet baryon rather than in a spin 1/2
   * octet baryon.
   */
  virtual double baryonDecupletProbability(long iq1, long iq2) const;

  /**
   * Return a pseudo scalar meson formed by the two quark flavours.
   */
  virtual tcPDPtr pseudoScalarMeson(long iq, long iqbar) const;

  /**
   * Return a vector meson formed by the two quark flavours.
   */
  virtual tcPDPtr vectorMeson(long iq, long iqbar) const;

  /**
   * Return a spin 1/2 octet baryon formed by the given quark and
   * diquark flavours.
   */
  virtual tcPDPtr baryonOctet(long iq, long idq) const;

  /**
   * Return a spin 3/2 decuplet baryon formed by the given quark and
   * diquark flavours.
   */
  virtual tcPDPtr baryonDecuplet(long iq, long idq) const;

  /**
   * Return the PDG code of a pseudo scalar meson formed by the two
   * quark flavours for \a iqh >= \a iql > 0.
   */
  virtual long pseudoScalarId(long iqh, long iql) const;

  /**
   * Return the PDG code of a vector meson formed by the two quark
   * flavours for \a iqh >= \a iql > 0.
   */
  virtual long vectorId(long iqh, long iql) const;

  /**
   * Return the PDG code for a spin 1/2 octet baryon formed by the
   * given quark flavours (\a iqa >= \a iqb >= \a iqc > 0). iq is one
   * of the flavours and the other two are assumed to be in a diquark
   * (in a spin-1 state if \a dqs1).
   */
  virtual long baryonOctetId(long iqa, long iqb, long iqc,
			     long iq, bool dqs1) const;

  /**
   * Return the PDG code for a spin 3/2 decuplet baryon formed by the
   * given quark flavours (\a iqa >= \a iqb >= \a iqc > 0).
   */
  virtual long baryonDecupletId(long iqa, long iqb, long iqc) const;

  /**
   * Return the PDG code of pseudo scalar mesons formed by the two
   * quark flavours (for \a iqh >= \a iql > 0), together with suitable
   * weights.
   */
  virtual vector< pair<long,double> >
  pseudoScalarIds(long iqh, long iql) const;

  /**
   * Return the PDG codes of vector mesons formed by the two quark
   * flavours (for \a iqh >= \a iql > 0), together with
   * suitable weights.
   */
  virtual vector< pair<long,double> > vectorIds(long iqh, long iql) const;

  /**
   * Return the PDG codes for spin 1/2 octet baryons formed by the
   * given quark flavours (\a iqa >= \a iqb >= \a iqc > 0) together
   * with suitable weights. iq is one of the flavours and the other
   * two are assumed to be in a diquark (in a spin-1 state if \a dqs1).
   */
  virtual vector< pair<long,double> >
  baryonOctetIds(long iqa, long iqb, long iqc,
		 long iq, bool dqs1) const;

  /**
   * Return the PDG codes for spin 3/2 decuplet baryons formed by
   * the given quark flavours (\a iqa >= \a iqb >= \a iqc > 0) together with
   * suitable weights.
   */
  virtual vector< pair<long,double> >
  baryonDecupletIds(long iqa, long iqb, long iqc) const;

  /**
   * Clear all cashed weights.
   */
  void clear();

  /**
   * Return the SU(6) weight for the given quark and di-quark flavours
   * to end up with in a baryon with the given spin (2S+1).
   */
  static double weightSU6QDiQSpin(long iq, long idq, int spin);

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
   * Standard Init function used to initialize the interfaces.
   */
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

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();
  //@}

private:

  /**
   * Suppression factor of strange quarks w.r.t. u and d.
   */
  double theSSup;

  /**
   * Suppression factor for di-quarks w.r.t. quarks.
   */
  double theDiSup;

  /**
   * Suppression of spin-1 di-quarks w.r.t. spin-0 ones.
   */
  double theDi1Sup;

  /**
   * Suppression of strange di-quarks w.r.t. u and d ones in addition
   * to the standard strangness suppression of quarks.
   */
  double theDiSSup;

  /**
   * Extra suppression of eta's.
   */
  double theEtaSup;

  /**
   * Extra suppression of ets-prime's.
   */
  double theEtaPSup;

  /**
   * Extra suppression for baryons of the spin 3/2 decuplet.
   */
  double theBaryon10Sup;

  /**
   * Probability that light (u/d) mesons has spin 1.
   */
  double thePSpin1;

  /**
   * Probability that strange mesons has spin 1.
   */
  double thePSpinS1;

  /**
   * Probability that charmed and heavier mesons has spin 1.
   */
  double thePSpinC1;

  /**
   * A selector used to weight the creation of (di)quark-anti(di)quark
   * pairs.
   */
  mutable VSelector<long> theFlavourSelector;

  /**
   * A map of selectors to cash probabilities for generateHadron.
   */
  mutable ProbabilityMap theProbabilities;


private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<SimpleFlavour> initSimpleFlavour;

  /**
   * Private and non-existent assignment operator.
   */
  SimpleFlavour & operator=(const SimpleFlavour &);

};

}


namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * This template specialization informs ThePEG about the base class of
 * SimpleFlavour.
 */
template <>
struct BaseClassTrait<SimpleFlavour,1>: public ClassTraitsType {
  /** Typedef of the base class of SimpleFlavour. */
  typedef FlavourGenerator NthBase;
};

template <>
/**
 * This template specialization informs ThePEG about the name of the
 * SimpleFlavour class and the shared object where it is defined.
 */
struct ClassTraits<SimpleFlavour>
  : public ClassTraitsBase<SimpleFlavour> {
  /** Return the class name.  */
  static string className() { return "ThePEG::SimpleFlavour"; }
  /**
   * Return the name of the shared library to be loaded to get access
   * to the SimpleFlavour class and every other class it uses (except
   * the base class).
   */
  static string library() { return "SimpleFlavour.so"; }

};

/** @endcond */

}

#endif /* THEPEG_SimpleFlavour_H */
