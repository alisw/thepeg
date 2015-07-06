// -*- C++ -*-
#ifndef DIPSY_SimpleProton_H
#define DIPSY_SimpleProton_H
//
// This is the declaration of the SimpleProton class.
//

#include "Ariadne/DIPSY/WaveFunction.h"
#include "ThePEG/PDT/SimpleBaryonRemnantDecayer.h"

namespace DIPSY {

using namespace ThePEG;

/**
 * The SimpleProton class represents the unevolved proton wave
 * function described in terms of an equilatteral triangle of dipoles
 * with a size distributed as a Gaussian.
 *
 * @see \ref SimpleProtonInterfaces "The interfaces"
 * defined for SimpleProton.
 */
class SimpleProton: public WaveFunction {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  SimpleProton();

  /**
   * The destructor.
   */
  virtual ~SimpleProton();
  //@}

public:

  /**
   * Get the width of the Gaussian distribution.
   */
  InvEnergy R() const;

  /**
   * Get the shift in the average of the Gaussian distribution.
   */
  InvEnergy r0() const;

  /**
   * Get the width of the distribution of the angle of the proton.
   */
  inline double angleWidth() const {
    return theAngleWidth;
  }

  /**
   * Get the width of the distribution of the angle of the proton.
   */
  inline double rapidityWidth() const {
    return theRapidityWidth;
  }

  /**
   * Get the maximum mass difference between the valence system and
   * the hadron allowed for collapsing into the original hadron after
   * the evolution.
   */
  inline Energy collapseTolerance() const {
    return theCollapseTolerance;
  }

  /**
   * Controls if the ends of the dipoles are connected or not.
   */
  inline int connectOption() const {
    return theConnected;
  }

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

  /**
   * Generate the shifted gaussian distribution.
   */
  virtual InvEnergy rndShiftGauss() const;

  /**
   * Set the corresponding particle. Check that it is actually a
   * proton or at least a baryon.
   */
  virtual void setParticle(PDPtr);

  /**
   * Fix up valens partons if they were not of the correct flavour or
   * if they should collapse into a hadron.
   */
  virtual void fixValence(Step & step, tPPtr particle, const vector<PPtr> & valence) const;

  /**
   * try to collapse the valence gluons into the original
   * particle. @return true if succeeded.
   */
  virtual bool collapseToProton(Step &, tPPtr particle, const vector<PPtr> &) const;

  /**
   * Pick a valence gluon and split it into a quark-diquark pair.
   */
  virtual bool splitValence(Step & step, tPPtr particle, const vector<PPtr> & valence) const;
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
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).

private:

  /**
   * The width of the Gaussian distribution.
   */
  InvEnergy theR;

  /**
   * The shift in the average of the Gaussian distribution.
   */
  InvEnergy theR0;

  /**
   * The width of the gaussian in rapidity for the individual partons.
   */
  double theRapidityWidth;

  /**
   * The width of the gaussian in shape for the individual partons.
   */
  double theAngleWidth;

  /**
   * Controls if the ends of the dipoles are connected or not.
   */
  int theConnected;

  /**
   * The maximum mass difference between the valence system and the
   * hadron allowed for collapsing into the original hadron after the
   * evolution.
   */
  Energy theCollapseTolerance;

  /**
   * A pointer to a RemnantDecayer object which is able to handle the
   * remnants of the particle.
   */
  Ptr<SimpleBaryonRemnantDecayer>::pointer remdec;

private:

  /**
   * Excleption class
   */
  struct RemnantException: public Exception {};
  
  /**
   * Utility function for the interface.
   */
  void setDecayer(RemDecPtr rd);

  /**
   * Simple helper function.
   */
  static void swapColourLines(tPPtr gluon);

  /**
   * Simple helper function.
   */
  static void insertIntermediate(Step & step, tPPtr parent, tPPtr children);

  /**
   * Simple helper function.
   */
  static void eraseIntermediate(Step & step, tPPtr p);

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SimpleProton & operator=(const SimpleProton &);

};

}

#endif /* DIPSY_SimpleProton_H */
