// -*- C++ -*-
#ifndef DIPSY_SimpleNucleus_H
#define DIPSY_SimpleNucleus_H
//
// This is the declaration of the SimpleNucleus class.
//

#include "Ariadne/DIPSY/WaveFunction.h"
#include "ThePEG/Analysis/FactoryBase.h"

namespace DIPSY {

using namespace ThePEG;

/**
 * The SimpleNucleus class represents the unevolved proton wave
 * function described in terms of an equilatteral triangle of dipoles
 * with a size distributed as a Gaussian.
 *
 * @see \ref SimpleNucleusInterfaces "The interfaces"
 * defined for SimpleNucleus.
 */
class SimpleNucleus: public WaveFunction {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  SimpleNucleus()
  : A(208), Z(82), R(6.68*femtometer), a(0.546*femtometer), w(0.0),
    wf(WaveFunctionPtr()), wfp(WaveFunctionPtr()), wfn(WaveFunctionPtr()),
    Rn(1.3*femtometer), useInterNucleonSwing(true), nTry(1), doRecenter(false),
    rdist(0) {}

  /**
   * The destructor.
   */
  virtual ~SimpleNucleus();
  //@}

public:

  /** @name Main virtual functions to be overridden in sub-classes. */
  //@{
  /**
   * Initialize the wavefunction for the given DipoleEventHandler.
   */
  virtual void initialize(const DipoleEventHandler &);

  /**
   * Generate a dipole state according to this wavefunction, given an
   * event handler and a positive light-cone momentum component.
   */
  virtual DipoleStatePtr generate(const DipoleEventHandler & eh, Energy plus);

  /**
   * Return the invariant mass squared of the particle.
   */
  virtual Energy2 m2() const;
  //@}

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  virtual void dofinish();

  //@}

public:

  /** @Name Simple access functions. */
  //@{
  /**
   * Return the nucleon wave function.
   */
  tWaveFunctionPtr nucleon() const {
    return wf;
  }

  /**
   * Return the proton wave function.
   */
  tWaveFunctionPtr proton() const {
    return wfp;
  }

  /**
   * Return the number of protons.
   */
  int nProtons() const {
    return abs(Z);
  }

  /**
   * Return the neutron wave function.
   */
  tWaveFunctionPtr neutron() const {
    return wfn;
  }

  /**
   * Return the number of protons.
   */
  int nNeutrons() const {
    return abs(A) - nProtons();
  }

  /**
   * Return the mass number.
   */
  int massNumber() const {
    return A;
  }

  /**
   * Return true if we want inter-nucleon swings.
   */
  bool interNucleonSwing() const {
    return useInterNucleonSwing;
  }

  /**
   * Return true if the nucleons should be recentered to put the
   * center of mass at zero impact parameter after distribution of
   * nucleons.
   */
  inline bool recenter() const {
    return doRecenter;
  }

  /**
   * Return the Wood-Saxon potential.
   */
  double ws(Length r) const {
    Length Rc = min(R - Rn/2.0, R);
    return (1.0 + w*sqr(r)/sqr(Rc))/(1.0 + exp((r - Rc)/a));
  }

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

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const {
    return new_ptr(*this);
  }

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const {
    return new_ptr(*this);
  }
  //@}


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).

private:

  /**
   * The mass number of the nucleaus.
   */
  int A;

  /**
   * The charge of the nucleaus.
   */
  int Z;

  /**
   * The size of the nucleus.
   */
  Length R;

  /**
   * The suppression at the edge of the nucleus.
   */
  Length a;

  /**
   * Optional quadratic term in Wood-Saxon potential
   */
  double w;

  /**
   * Wave function of the individual nucleons.
   */
  WaveFunctionPtr wf;

  /**
   * Wave function of the individual protons.
   */
  WaveFunctionPtr wfp;

  /**
   * Wave function of the individual neutrons.
   */
  WaveFunctionPtr wfn;

  /**
   * Size of an individual nucleon for the exclusion mechanism.
   */
  Length Rn;

  /**
   * Allow swings between nucleons in the nuclei.
   */
  bool useInterNucleonSwing;

  /**
   * Number of times to try new angles for the same r-value if
   * overlapping nucleons are found. If zero, the whole set of
   * generated positions is discarded if anoverlap is found.
   */
  int nTry;

  /**
   * Return true if the nucleons should be recentered to put the
   * center of mass at zero impact parameter after distribution of
   * nucleons.
   */
  bool doRecenter;

  /** The distribution of nucleon centers. */
  FactoryBase::tH1DPtr rdist;


  /**
   * Helper function to set measured parametrs of a nucleus given by
   * the name.
   */
  string setParameters(string);

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SimpleNucleus & operator=(const SimpleNucleus &);

};

}

#endif /* DIPSY_SimpleNucleus_H */
