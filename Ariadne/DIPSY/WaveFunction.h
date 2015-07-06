// -*- C++ -*-
#ifndef DIPSY_WaveFunction_H
#define DIPSY_WaveFunction_H
//
// This is the declaration of the WaveFunction class.
//

#include "ThePEG/Handlers/HandlerBase.h"
#include "ThePEG/PDT/ParticleData.h"
#include "WaveFunction.fh"
#include "DipoleState.h"
#include "DipoleEventHandler.fh"

namespace DIPSY {

using namespace ThePEG;

/**
 * WaveFunction is the base class for wavefunction objects capable of
 * generating initial DipoleState objects.
 *
 * @see \ref WaveFunctionInterfaces "The interfaces"
 * defined for WaveFunction.
 */
class WaveFunction: public HandlerBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline WaveFunction();

  /**
   * The copy constructor.
   */
  inline WaveFunction(const WaveFunction &);

  /**
   * The destructor.
   */
  virtual ~WaveFunction();
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
  virtual DipoleStatePtr
  generate(const DipoleEventHandler & eh, Energy plus) = 0;

  /**
   * Fix up valens partons if they were not of the correct flavour or
   * if they should collapse into a hadron.
   */
  virtual void fixValence(Step & step, tPPtr particle, const vector<PPtr> & valence) const;

  /**
   * Return the invariant mass squared of the particle.
   */
  virtual Energy2 m2() const = 0;
  //@}

  /**
   * The DipoleEventHandler in charge of the generation.
   */
  inline tcDipoleEventHandlerPtr eventHandler() const;

  /**
   * Get the corresponding particle.
   */
  inline tcPDPtr particle() const;

  /**
   * Set the corresponding particle. May be overridden by sub classes
   * to check the particle is relevant for the particular wave
   * function.
   */
  virtual void setParticle(PDPtr);

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


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


private:

  /**
   * The DipoleEventHandler in charge of the generation.
   */
  tcDipoleEventHandlerPtr theEventHandler;

  /**
   * The corresponding particle.
   */
  PDPtr theParticle;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  WaveFunction & operator=(const WaveFunction &);

};

}

#include "WaveFunction.icc"

#endif /* DIPSY_WaveFunction_H */
