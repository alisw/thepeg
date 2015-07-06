// -*- C++ -*-
#ifndef PYTHIA7_LundPtGenerator_H
#define PYTHIA7_LundPtGenerator_H
// This is the declaration of the LundPtGenerator class.

#include "FragConfig.h"
// #include "LundPtGenerator.fh"
// #include "LundPtGenerator.xh"
#include "ThePEG/Handlers/PtGenerator.h"


namespace Pythia7 {

/**
 * LundPtGenerator is the transverse momentum generator used in the
 * Lund string fragmentation scheme.  It generates the \f$(p_x,p_y)\f$
 * components of the transverse momentum of a quark (diquark) in a
 * \f$q-\bar{q}\f$ (or a \f$(qq)-(\overline{qq})\f$) pair created in
 * the string colour field, according to the flavour independent
 * Gaussian distribution in \f$p_x\f$ and \f$p_y\f$ including the
 * possibility of non-Gaussian tails.
 * 
 * LundPtGenerator inherits the from the PtGenerator class and
 * overrides the generate() function.
 *
 * Note that the 4-vector version of the transverse momentum is
 * obtained given the transverse vectors of the StringRegion where the
 * pair is produced.
 *
 * @see \ref LundPtGeneratorInterfaces "The interfaces" defined
 * for LundPtGenerator.
 *
 */
class LundPtGenerator: public ThePEG::PtGenerator {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  inline LundPtGenerator();

  /**
   * Copy-constructor.
   */
  inline LundPtGenerator(const LundPtGenerator &);

  /**
   * Destructor.
   */
  virtual ~LundPtGenerator();
  //@}

public:

  /** @name Virtual functions required by the PtGenerator class. */
  //@{
  /**
   * Generate (\f$k_x, k_y\f$) components of the transverse
   * momentum. They will be distributed as
   * \f$\exp(-k_\perp^2/\sigma^2)k_\perp dk_\perp\f$ with
   * \f$k_\perp^2=k_x^2+k_y^2\f$ and \f$\sigma=\f$ Sigma().
   */
  virtual TransverseMomentum generate() const ;
  //@}

  /** @name Access to parameters. */
  //@{
  /**
   * Get the gaussian width.
   */
  inline Energy Sigma() const;

  /**
   * Get the non-Gaussian fraction of the Gaussian transverse momentum
   * distribution to be enhanced by the factor nGaussFactor().
   */
  inline double NGaussFraction() const;

  /**
   * Get the non-Gaussian tails enhancement factor.
   */
  inline double NGaussFactor() const;
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
   * The gaussian width.
   */
  Energy sigma;             //PARJ(21)

  /**
   * The non-Gaussian fraction of the Gaussian transverse momentum
   * distribution to be enhanced by the factor nGaussFactor().
   */
  double nGaussfraction;    //PARJ(23) 

  /**
   * The non-Gaussian tails enhancement factor.
   */
  double nGaussfactor;      //PARJ(24)

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<LundPtGenerator> initLundPtGenerator;

  /**
   * Private and non-existent assignment operator.
   */
  LundPtGenerator & operator=(const LundPtGenerator &);
};

}

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * This template specialization informs ThePEG about the base class of
 * Pythia7::LundPtGenerator.
 */
template <>
struct BaseClassTrait<Pythia7::LundPtGenerator,1>: public ClassTraitsType {
  /** Typedef of the base class of Pythia7::LundPtGenerator. */
  typedef PtGenerator NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * Pythia7::LundPtGenerator class and the shared object where it is
 * defined.
 */
template <>
struct ClassTraits<Pythia7::LundPtGenerator>
  : public ClassTraitsBase<Pythia7::LundPtGenerator> {
  /** Return the class name. */
  static string className() { return "Pythia7::LundPtGenerator"; }
  /** Return the name of the shared library to be loaded to get access
   * to the Pythia7::LundPtGenerator class and every other class it
   * uses (except the base class). */
  static string library() { return "libP7String.so"; }
};

/** @endcond */

}

#include "LundPtGenerator.icc"
#ifndef PYTHIA7_TEMPLATES_IN_CC_FILE
// #include "LundPtGenerator.tcc"
#endif

#endif /* PYTHIA7_LundPtGenerator_H */
