// -*- C++ -*-
#ifndef ARIADNE5_RRGluonEmission_H
#define ARIADNE5_RRGluonEmission_H
//
// This is the declaration of the RRGluonEmission class.
//

#include "FSGluonEmission.h"
#include "Ariadne/Cascade/EmitterBase.h"
#include "ThePEG/Utilities/UtilityBase.h"

namespace Ariadne5 {

using namespace ThePEG;

/**
 * The RRGluonEmission class contains all information about a
 * generated and performed final state gluon emission from a remnant
 * parton.
 */
class RRGluonEmission: public FSGluonEmission {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The constructor relevant for real emissions.
   */
  RRGluonEmission(const EmitterBase & inmodel, const DipoleBase & indipole,
		 tRemParPtr inrem1, tRemParPtr inrem3,
		  const Lorentz5Momentum & inph)
    : FSGluonEmission(inmodel, indipole, 0.0, 0.0),
      rem1(inrem1), rem3(inrem3), mh2(inph.mass2()), oph(inph),
      ophr1(inrem1->getBoost()*inph), ophr3(inrem3->getBoost()*inph), ph(inph) {
    pold = make_pair(rem1->momentum(), rem3->momentum());
    Rrcm = Utilities::getBoostToCM(pold);
    Rircm = Rrcm.inverse();
    ophr = Rrcm*oph;
    opr1 = rem1->getBoost()*pold.first;
    opr3 = rem3->getBoost()*pold.second;
  }

  /**
   * The deault constructor should not normally be used.
   */
  RRGluonEmission() {}

  /**
   * The destructor.
   */
  virtual ~RRGluonEmission();
  //@}

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual ClonePtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual ClonePtr fullclone() const;
  //@}


public:

  /**
   * Remnant 1
   */
  tRemParPtr rem1;

  /**
   * Remnant 3
   */
  tRemParPtr rem3;

  /**
   * The squared mass of the hard sub-system.
   */
  Energy2 mh2;

  /**
   * The momentum of the hard sub-system before the emission.
   */
  Lorentz5Momentum oph;

  /**
   * The momentum of the hard sub-system before the emission in the
   * remnants rest system.
   */
  Lorentz5Momentum ophr;

  /**
   * The momentum of the hard sub-system before the emission in the
   * system of remnant 1.
   */
  Lorentz5Momentum ophr1;

  /**
   * The momentum of the hard sub-system before the emission in the
   * system of remnant 3.
   */
  Lorentz5Momentum ophr3;

  /**
   * The momentum of the remnant 1 in its own system
   */
  Lorentz5Momentum opr1;

  /**
   * The momentum of the remnant 3 in its own system
   */
  Lorentz5Momentum opr3;

  /**
   * The momentum of the hard sub-system after the emission.
   */
  Lorentz5Momentum ph;

  /**
   * The rotation of the hard sub-system giving oph -> ph
   */
  mutable LorentzRotation Rh;

  /**
   * The rotation from the lab system to the remnants rest frame.
   */
  LorentzRotation Rrcm;

  /**
   * The inverse rotation from the lab system to the remnants rest frame.
   */
  LorentzRotation Rircm;

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
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  RRGluonEmission & operator=(const RRGluonEmission &);

};

}

#endif /* ARIADNE_RRGluonEmission_H */
