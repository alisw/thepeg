// -*- C++ -*-
//
// FlavourGenerator.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_FlavourGenerator_H
#define ThePEG_FlavourGenerator_H
// This is the declaration of the FlavourGenerator class.

#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/Handlers/HandlerBase.h"

namespace ThePEG {

/**
 * FlavourGenerator is an abstract base class to be used to implement
 * models describing the quark content of hadrons. FlavourGenerator
 * inherits from the HandlerBase class.
 *
 * The interface is based on the flavour generation implementation in
 * Pythia but is general enough to be used in other situations. The
 * main virtual functions to be overridden in subclasses are
 * generateHadron(tcPDPtr), getHadron(tcPDPtr, tcPDPtr),
 * getHadron(long, long) getBaryon(tcPDPtr, tcPDPtr, tcPDPtr),
 * getBaryon(long, long, long), selectQuark() and selectFlavour(). In
 * this base class the getHadron(tcPDPtr, tcPDPtr) and getHadron(long,
 * long) are implemented to call eachother, so a subclass must
 * implement at least one of them. The same thing is true for
 * getBaryon(tcPDPtr, tcPDPtr, tcPDPtr) and getBaryon(long, long,
 * long)
 *
 * @see \ref FlavourGeneratorInterfaces "The interfaces"
 * defined for FlavourGenerator.
 * @see HandlerBase
 */
class FlavourGenerator: public HandlerBase {

public:

  /** @name Virtual functions to be overridden by subclasses. */
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
  virtual tcPDPair generateHadron(tcPDPtr quark) const = 0;

  /**
   * Same as generateHadron(tcPDPtr), but throws an exception if no
   * hadron could be produced.
   */
  tcPDPair alwaysGenerateHadron(tcPDPtr quark) const;

  /**
   * Get hadron from flavours. Return a hadron with the flavour
   * content given by the (anti-)(di-)quarks in the argument. The
   * arguments are given as ParticleData pointers. The default
   * versions will call the getHadron(long, long).
   * @param q1 the first flavour.
   * @param q2 the second flavour.
   * @return the corresponding hadron type or null if none could be
   * generated.
   */
  virtual tcPDPtr getHadron(tcPDPtr q1, tcPDPtr q2) const;

  /**
   * Get hadron from flavours. Return a hadron with the flavour
   * content given by the (anti-)(di-)quarks in the argument. The
   * arguments are given as PDG codes. The default
   * versions will call the getHadron(tcPDPtr, tcPDPtr).
   * @param iq1 the PDG code of the first flavour.
   * @param iq2 the PDG code of the second flavour.
   * @return the corresponding hadron type or null if none could be
   * generated.
   */
  virtual tcPDPtr getHadron(long iq1, long iq2) const;

  /**
   * Same as getHadron(tcPDPtr, tcPDPtr) but thows an exception if no
   * hadron could be produced.
   */
  tcPDPtr alwaysGetHadron(tcPDPtr q1, tcPDPtr q2) const;

  /**
   * Same as getHadron(long, long) but thows an exception if no hadron
   * could be produced.
   */
  tcPDPtr alwaysGetHadron(long iq1, long iq2) const;

  /**
   * Return a baryon with the flavour content given by the
   * (anti)quarks in the argument.  The arguments are given as
   * particle data pointers. The default versions will call
   * getBaryon(long, long, long). If no corresponding hadron was
   * formed it should return the null pointer.
   * @param q1 the first flavour.
   * @param q2 the second flavour.
   * @param q3 the third flavour.
   * @return the corresponding baryon type or null if none could be
   * generated.
   */
  virtual tcPDPtr getBaryon(tcPDPtr q1, tcPDPtr q2, tcPDPtr q3) const;

  /**
   * Return a baryon with the flavour content given by the
   * (anti)quarks in the argument.  The arguments are given as
   * particle data pointers. The default versions will call
   * getBaryon(tcPDPtr, tcPDPtr, tcPDPtr). If no corresponding hadron was
   * formed it should return the null pointer.
   * @param q1 the PDG code of the first flavour.
   * @param q2 the PDG code of the second flavour.
   * @param q3 the PDG code of the third flavour.
   * @return the corresponding baryon type or null if none could be
   * generated.
   */
  virtual tcPDPtr getBaryon(long q1, long q2, long q3) const;

  /**
   * Same as getBaryon(tcPDPtr, tcPDPtr, tcPDPtr), but throws an
   * exception if no baryon could be produced.
   */
  tcPDPtr alwaysGetBaryon(tcPDPtr q1, tcPDPtr q2, tcPDPtr q3) const;
  /**
   * Same as getBaryon(long, long, long), but throws an exception if
   * no baryon could be produced.
   */
  tcPDPtr alwaysGetBaryon(long q1, long q2, long q3) const;

  /**
   * Generate a random quark flavour.
   */
  virtual long selectQuark() const = 0;

  /**
   * Generate a random (di)quark flavour.
   */
  virtual long selectFlavour() const = 0;
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
   * Standard Init function used to initialize the interface.
   */
  static void Init();


private:

  /**
   * Describe aa abstract class without persistent data.
   */
  static AbstractNoPIOClassDescription<FlavourGenerator> initFlavourGenerator;

  /**
   *  Private and non-existent assignment operator.
   */
  FlavourGenerator & operator=(const FlavourGenerator &);

};

/** @cond EXCEPTIONCLASSES */
/** An Exception class used by FlavourGenerator classes if no hadrons
    could be generated. */
class FlavourGeneratorException: public Exception {};
/** @endcond */

/** @cond TRAITSPECIALIZATIONS */

/**
 * This template specialization informs ThePEG about the
 * base class of FlavourGenerator.
 */
template <>
struct BaseClassTrait<FlavourGenerator,1>: public ClassTraitsType {
  /** Typedef of the base class of FlavourGenerator. */
  typedef HandlerBase NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * DecayHandler class.
 */
template <>
struct ClassTraits<FlavourGenerator>:
    public ClassTraitsBase<FlavourGenerator> {
  /** Return the class name. */
  static string className() {
    return "ThePEG::FlavourGenerator";
  }
};

/** @endcond */

}

#endif /* ThePEG_FlavourGenerator_H */
