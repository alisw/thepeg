// -*- C++ -*-
#ifndef DIPSY_WFInfo_H
#define DIPSY_WFInfo_H
//
// This is the declaration of the WFInfo class.
//

#include "ThePEG/Config/ThePEG.h"
#include "WFInfo.fh"
#include "WaveFunction.fh"
#include "ThePEG/EventRecord/EventInfoBase.h"
#include "ThePEG/EventRecord/Particle.h"

namespace DIPSY {

using namespace ThePEG;

/**
 * WFInfo is the base class for additional information about the
 * WaveFunction used for constructing an intitial Dipole state. This
 * base class keeps information about the size, r, of the initial
 * DipoleSystem. The class does not have any virtual
 * functions. Instead dynaic_cast should be used to access information
 * from sub-classes.
 */
class WFInfo: public ThePEG::EventInfoBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor giving the wave function that produces
   * this object and the size as argument.
   */
  inline WFInfo(tcWaveFunctionPtr wfin = tcWaveFunctionPtr(),
		InvEnergy rini = 0.0*InvGeV)
    : theWF(wfin), theR(rini) {}

  /**
   * The destructor.
   */
  virtual ~WFInfo();
  //@}

public:

  /**
   * Get the wave function that produced this info object.
   */
  inline const WaveFunction & wf() const {
    return *theWF;
  }

  /**
   * The initial size of the DipoleSystem.
   */
  inline InvEnergy r() const {
    return theR;
  }

  /**
   * Return the WFInfo object stored with the given particle.
   */
  static tcWFInfoPtr getWFInfo(const Particle & particle);

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

private:

  /**
   * The wave function that produced this info object.
   */
  tcWaveFunctionPtr theWF;

  /**
   * The initial size of the DipoleSystem.
   */
  InvEnergy theR;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  WFInfo & operator=(const WFInfo &);

};

}

#endif /* DIPSY_WFInfo_H */
