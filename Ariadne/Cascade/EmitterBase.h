// -*- C++ -*-
#ifndef Ariadne5_EmitterBase_H
#define Ariadne5_EmitterBase_H
//
// This is the declaration of the EmitterBase class.
//

#include "ThePEG/Handlers/HandlerBase.h"
#include "DipoleBase.h"
#include "QCDDipole.fh"
#include "ReweightBase.h"
#include "EmitterBase.fh"
#include "Emission.fh"
#include "ThePEG/Utilities/Triplet.h"
#include "ThePEG/Utilities/Current.h"
#include "AriadneHandler.h"

namespace Ariadne5 {

using namespace ThePEG;

/**
 * EmitterBase is the base class of all Ariadne classes implementing
 * a specific model for emission from different kinds of dipoles. A
 * sub-class must implement a sub-class of Emission to store the
 * result of an emission.
 *
 * @see \ref EmitterBaseInterfaces "The interfaces"
 * defined for EmitterBase.
 */
class EmitterBase: public HandlerBase {

public:

  /** Convenient typedef. */
  typedef Triplet<Lorentz5Momentum,Lorentz5Momentum,Lorentz5Momentum> Trip;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  EmitterBase(bool reco = false): canReconstruct(reco), willReconstruct(reco) {}

  /**
   * The destructor.
   */
  virtual ~EmitterBase() {}
  //@}

public:

  /** @name Virtual functions to be overridden by sub-classes. */
  //@{
  /**
   * Return true if and only if this emitter can handle the given \a
   * dipole.
   */
  virtual bool canHandle(const DipoleBase & dipole) const = 0;

  /**
   * If the given \a emitter overlaps with this model for the given \a
   * dipole, return true if this object should take precedence. Must
   * only be called for a \a dipole for which canHandle() is true.
   */
  virtual bool overrides(const EmitterBase & emitter,
			 DipoleBase & dipole) const = 0;

  /**
   * Check if objects related to the given \a dipole have been touched
   * in a way such that emissions must be regenerated.
   */
  virtual bool touched(const DipoleBase & dipole) const;

  /**
   * Generate the a phase space point for an emission corresponding to
   * this model. Must only be called for a \a dipole for which
   * canHandle() is true.
   */
  virtual EmPtr generate(const DipoleBase & dipole,
			 Energy rhomin, Energy rhomax) const = 0;

  /**
   * Perform the \a emission previously generated.
   * @return true if the emission was successful
   */
  virtual bool perform(const Emission & emission) const = 0;

  /**
   * Reverse a previously performed emission. Sub-classes which has
   * signalled that they can revert an emission but fails to do so,
   * must throw a Exception::runerror.
   */
  virtual void revert(const Emission & emission) const;
  
  /**
   * Return true if this emitter should be used to reconstruct
   * emissions in the CKKW-L algorithm.
   */
  bool reconstructor() const {
    return willReconstruct;
  }

  /**
   * Return a list of inverse emissions which this model could have
   * performed to arrive at the given \a state.
   */
  virtual vector<EmPtr> inverseEmissions(const DipoleState & state) const;

  /**
   * Check if an inverse \a emission suggested by the given \a emitter
   * should rather be handled by this model.
   */
  virtual bool overrideInverse(const Emission & emission) const;

  /**
   * Perform the inverse \a emission, previously reported by
   * inverseEmissions().
   * @return false of the reconstruction failed.
   */
  virtual bool performInverse(const Emission &, DipoleState &) const;

  /**
   * Possibility to override the general non-perturbative cutoff for
   * this process. The default is AriadneHandler::pTCut() which is
   * typically used for QCD radiation.
   */
  virtual Energy rhoCut() const;

  //@}

  /** @name These functions are related to external reweighting of
   *  emissions and should normally not be overridden. */
  //@{
  /**
   * Return the product of preweights from all reweighters registered
   * for this model.
   */
  virtual double preweight(const Emission & emission) const;

  /**
   * Return the product of reweights from all reweighters registered
   * for this emission model.
   */
  virtual double reweight(const Emission & emission) const;

  /**
   * Return true if any of the reweighters registered for this
   * emission model says so.
   */
  virtual bool finalVeto(const Emission & emission) const;

  /**
   * Return true if any of the reweighters registered for this
   * emission model may veto in finalVeto().
   */
  virtual bool hasFinalVeto() const;
  //@}

  /**
   * Access the AriadneHandler object currently in use.
   */
  static AriadneHandler & handler() {
    return Current<AriadneHandler>::current();
  }

  /**
   * Generate an squared invariant transverse momentum. The
   * distribution used is \f$C\alpha_S\frac{dp_\perp^2}{p_\perp^2}\f$
   * multiplied with the corresponding Sudakov formfactor. The maximum
   * \f$p_\perp^2\f$ is \a pt2max and \a C is the constant. If we are
   * using a running \f$\alpha_S\f$ runrdnsud() will be called
   * automatically.
   */
  static Energy2 rndsud(double C, Energy2 pt2max, Energy2 pt2min);

  /**
   * Generate an squared invariant transverse momentum. The
   * distribution used is \f$C\alpha_S\frac{dp_\perp^2}{p_\perp^2}\f$
   * multiplied with the corresponding Sudakov formfactor. The maximum
   * \f$p_\perp^2\f$ is \a pt2max and \a C is the constant. This
   * function takes care of the running of \f$\alpha_S\f$ or whatever.
   */
  static Energy2 runrndsud(double C, Energy2 pt2max, Energy2 pt2min);

  /**
   * Check consistency of generated phase-space point for standard
   * dipole variables.
   */
  static bool check(double x1, double x3, double y1, double y2, double y3);

  /**
   * Return the invariant transverse momentum of parton \a p2
   * w.r.t. partons \a p1 and \a p3.
   */
  static Energy2 invPT2(tcParPtr p1, tcParPtr p2, tcParPtr p3);

  /**
   * Return the invariant rapidity of parton \a p2 w.r.t. partons \a
   * p1 and \a p3.
   */
  static double invY(tcParPtr p1, tcParPtr p2, tcParPtr p3);

  /**
   * Return the momenta of three particles with masses \a m1, \a m2
   * and \a m3 respectively and where the first particle takes an
   * energy fraction \a x1 and the third an energy fraction \a x3
   * (x1+x2+x3=2). If \a nr1 is true the first particle will b along
   * the positive z-axis. Conversely if \a nr3 is true the third
   * particle will always be along the negative z-axis. If \a nr1 and
   * \a nr3 is true, the summed squared transverse energy of the
   * first and third momenta will be minimized, and if niether \a nr1
   * or \a nr3 is true, \a nr1 (\a nr3) will be set to true with
   * probability \f$x_{1(3)}^2/(x_1^2+x_3^2)\f$. If \a userho is true,
   * the recoils will be determined by the momentum fractions rather
   * than the energy fractions.
   */
  static Trip getMomenta(Energy2 s, double x1, double x3,
			 Energy m1, Energy m2, Energy m3,
			 bool nr1, bool nr3, bool userho = false);

  /**
   * Same as
   * getMomenta(Energy2,double,double,Energy,Energy,Energy,bool,bool,bool),
   * but also rotate an angle phi around z-axis
   */
  static Trip getMomenta(Energy2 s, double x1, double x3,
			 Energy m1, Energy m2, Energy m3,
			 bool nr1, bool nr3, double phi, bool userho = false);

  /**
   * Calls
   * getMomenta(Energy2,double,double,Energy,Energy,Energy,bool,bool,bool)
   * But also rotates the resulting momenta an angle \a phi around the
   * z-axis and transforms them from the rest frame of \a p1 and \a
   * p3). \a s must be the squared invariant mass of \a p1 and \a p3.
   */
  static Trip getMomenta(Energy2 s, double x1, double x3,
			 Energy m1, Energy m2, Energy m3,
			 bool nr1, bool nr3, double phi,
			 const Lorentz5Momentum & p1,
			 const Lorentz5Momentum & p3, bool userho = false);

protected:

  /**
   * Exception class used if weights larger that one is encountered.
   */
  struct WeightException: Exception {};

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

public:

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).

private:

  /**
   * A list of reweighting objects which may be applied to
   * emissions generated by this emission model.
   */
  vector<DipoleRWPtr> theReweighters;

  /**
   * If true, this emitter can be used to perform inverse emissions in
   * the CKKW-L algorithm.
   */
  bool canReconstruct;

  /**
   * If true, this emitter will be used to perform inverse emissions in
   * the CKKW-L algorithm.
   */
  bool willReconstruct;

private:

  /**
   * Utility function for the interface.
   */
  void setReconstruct(bool t);

  /**
   * Utility function for the interface.
   */
  bool defReconstruct() const;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  EmitterBase & operator=(const EmitterBase &);

};

}

#endif /* Ariadne5_EmitterBase_H */
