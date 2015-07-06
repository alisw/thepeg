// -*- C++ -*-
#ifndef PYTHIA7_LundZGenerator_H
#define PYTHIA7_LundZGenerator_H
// This is the declaration of the LundZGenerator class.

#include "FragConfig.h"
// #include "LundZGenerator.fh"
// #include "LundZGenerator.xh"
#include "ThePEG/Handlers/ZGenerator.h"


namespace Pythia7 {

/**
 * The LundZGenerator generates longitudinal scaling variable z of
 * hadron produced in the Lund string fragmentation scheme, according
 * to the Lund Symmetric Fragmentation Function,
 * \f$f(z)=(1/z^{c})(1-z)^{a}exp(-bm_\perp^2/z)\f$, where the \f$a\f$,
 * \f$b\f$ and \f$c\f$ factors are determined by the setShape()
 * function.
 *
 * For diuark production there is an enhancement of the effective
 * \f$a\f$ parameter (see deltaDQ).
 *
 * A heavy endpoint quark (above c quark) is treated according to the
 * Bowler modification of the Lund Symmetric fragmentation function.
 * (see the rQ() method).
 *
 * The LundZGenerator inherits all the necessary interfaces from the 
 * ZGenerator class.
 *
 * @see \ref LundZGeneratorInterfaces "The interfaces"
 * defined for LundZGenerator.
 */
class LundZGenerator: public ThePEG::ZGenerator {

public:


  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  inline LundZGenerator();

  /**
   * Copy-constructor.
   */
  inline LundZGenerator(const LundZGenerator &);

  /**
   * Destructor.
   */
  virtual ~LundZGenerator();
  //@}

public:

  /** @name Virtual functions mandated by the ZGenerator base class. */
  //@{
  /**
   * generate the scaling variable z of hadrons created at each step
   * of the fragmentation procedure, given the two hadron constituents
   * \a lastPD leftover at the previous step, and the \a newPD newly
   * created in the current step and the hadron transverse mass
   * squared \a mT2.
   */
  virtual double generate(cPDPtr lastPD, cPDPtr newPD, Energy2 mT2) const;

  /**
   * Return the default value of the \f$a\f$ parameter of the Lund Symmetric 
   * fragmentation function
   */
  inline double aSym() const;

  /**
   * Return the default value of the \f$b\f$ parameter of the Lund
   * Symmetric fragmentation function
   */
  inline InvEnergy2 bSym() const;
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


protected:

  /**
   * Determine the \f$a\f$, \f$b\f$ and \f$c\f$ parameters of the
   * fragmentation function parametrization, given the flavour content
   * of the two end-points \a lastPD and \a newPD and the hadron \a
   * mT2.
   */
  virtual void setShape(cPDPtr lastPD, cPDPtr newPD, Energy2 mT2) const;

  /**
   * Return the value of the \f$r_Q\f$ parameter of the Bowler
   * parametrization given the the heavy quark id number, \a hf. Note
   * that the same \f$r_Q\f$ value is returned for all \a hf above the
   * c-quark.
   */
  inline virtual double rQ(long hf) const;


protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
  //@}


protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Check sanity of the object during the setup phase.
   */
  inline virtual void doupdate() throw(UpdateException);

  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  inline virtual void doinit() throw(InitException);

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  inline virtual void dofinish();
  //@}

private:

  /**
   * Return the \f$a\f$ factor of the fragmentation function
   * parametrization.
   */
  inline double af() const;

  /**
   * Return the \f$b m_\perp^2\f$ factor of the fragmentation function
   * parametrization.
   */
  inline double bf() const;

  /**
   * Return the \f$c\f$ factor of the fragmentation function
   * parametrization.
   */
  inline double cf() const;

  /**
   * Return true when \a zm is below the minimal
   * value for Zmax.
   */
  inline bool smallZmax(double zm) const; 

  /**
   * Return true when \a zm is above the maximal value for Zmax.
   */
  inline bool largeZmax(double zm) const;   


  /**
   * Return true when both \a zm and \a bmt2 are large
   */
  inline bool largeZmax(double zm, double bmt2) const;   


  /**
   * Return true when \f$c\f$ is close to 1.
   */
  inline bool cCloseToOne() const;

  /**
   * Return true when \f$c\f$ is close to \f$a\f$.
   */
  inline bool cCloseToA() const;

private:

  /**
   * The default \f$a\f$ parameter of the Lund Symmetric fragmentation
   * function
   */
  double asym;         // PARJ(41) 

  /**
   * The default \f$b\f$ parameter of the Lund Symmetric fragmentation
   * function
   */
  InvEnergy2 bsym;     // PARJ(42)


  /**
   * The amount by which the effective a\f$a\f$ parameter in the Lund
   * flavour dependent symmetric fragmentation function is assumed to
   * be larger than the default \f$a\f$ when diquarks are produced.
   */
  double deltaDQ;      // PARJ(45) 

  /**
   * The value of the \f$r_Q\f$ parameter of the Bowler
   * parametrization for all quarks heavier than the c-quark.
   */
  double rQc;          // PARJ(46)


  /**
   * Current value of \f$bm_\perp^2\f$.
   */
  mutable double bmT2;

  /**
   * Current value of \f$a\f$.
   */
  mutable double theAf;

  /**
   * Current value of \f$c\f$.
   */
  mutable double theCf;

  /**
   * Init Interface Description/Map 
   */
  static ClassDescription<LundZGenerator> initLundZGenerator;


  /**
   *  Private and non-existent assignment operator.
   */
  LundZGenerator & operator=(const LundZGenerator &);

};

}

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * This template specialization informs ThePEG about the base class of
 * Pythia7::LundZGenerator.
 */
template <>
struct BaseClassTrait<Pythia7::LundZGenerator,1>: public ClassTraitsType {
  /** Typedef of the base class of Pythia7::LundZGenerator. */
  typedef ZGenerator NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * Pythia7::LundZGenerator class and the shared object where it is
 * defined.
 */
template <>
struct ClassTraits<Pythia7::LundZGenerator>
  : public ClassTraitsBase<Pythia7::LundZGenerator> {
  /** Return the class name.  */
  static string className() { return "Pythia7::LundZGenerator"; }
  /**
   * Return the name of the shared library to be loaded to get access
   * to the Pythia7::LundZGenerator class and every other class it
   * uses (except the base class).
   */
  static string library() { return "libP7String.so"; }
};

/** @endcond */

}

#include "LundZGenerator.icc"
#ifndef PYTHIA7_TEMPLATES_IN_CC_FILE
// #include "LundZGenerator.tcc"
#endif

#endif /* PYTHIA7_LundZGenerator_H */











