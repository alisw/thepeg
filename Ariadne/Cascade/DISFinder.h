// -*- C++ -*-
#ifndef Ariadne5_DISFinder_H
#define Ariadne5_DISFinder_H
//
// This is the declaration of the DISFinder class.
//

#include "ThePEG/Handlers/HandlerBase.h"
#include "DipoleState.fh"
#include "RemnantParton.fh"
#include "DISFinder.fh"

namespace Ariadne5 {

using namespace ThePEG;

/**
 * The DISFinder class is responsible for finding DIS-type scattered
 * leptons in a SubProcess. This base class only looks at the
 * mother-daughter relationships in the SubProcess. The case where
 * these relationships are not available is treated very crudely and
 * more sophisticated methods may implemented in sub-classes.
 *
 * @see \ref DISFinderInterfaces "The interfaces"
 * defined for DISFinder.
 */
class DISFinder: public HandlerBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  DISFinder();

  /**
   * The destructor.
   */
  virtual ~DISFinder();
  //@}

public:

  /** @name Main virtual function. */
  //@{
  /**
   * If DIS-like scattered leptons are found in the given SubProcess,
   * return properly setup corresponding hard remnants. The
   * DipoleState is supplied to allow for proper creation of
   * RemnantParton objects.
   */
  virtual pair<tRemParPtr,tRemParPtr>
  findDISLeptons(SubProcess &, DipoleState &) const;

  /**
   * If findDISLeptons() has found scattered leptons, find also the
   * scattered quarks in the SubProcess and return them as hard
   * remnants. The DipoleState is supplied to allow for proper
   * creation of RemnantParton objects.
   */
  virtual pair<tRemParPtr,tRemParPtr>
  findDISQuarks(pair<tRemParPtr,tRemParPtr> leptons,
		SubProcess &, DipoleState &) const;
  //@}

protected:

  /** @name Additional helper function. */
  //@{
  /**
   * Check if the given particle is a t-channel electro-weak boson.
   */
  bool virtualBoson(tPPtr) const;

  /**
   * Return the scattered lepton and virtual boson in a SubProcess if
   * the given \a incoming particle is a lepton.
   */
  virtual pair<tPPtr,tPPtr> getScattered(tPPtr incoming, SubProcess & sub) const;

  /**
   * Create and setup a remnant given the incoming and scattered
   * lepton, and the boson (if found). The DipoleState is supplied to
   * allow for proper creation of RemnantParton objects.
   */
  virtual tRemParPtr
  createRemnant(tPPtr inc, tPPtr lep, tPPtr bos, DipoleState & state) const;

  /**
   * Find the quark in a SubProcess corresponding to the scattered \a
   * lepton.
   */
  virtual tPPtr getQuark(tRemParPtr lepton, SubProcess & sub) const;

  /**
   * Create a hard remnant corresponding to the scattered \a quark,
   * given a scattered \a lepton. The DipoleState is supplied to
   * allow for proper creation of RemnantParton objects.
   */
  tRemParPtr createQuarkRemnant(tRemParPtr lepton, tPPtr quark,
				DipoleState & state) const;
  //@}

public:

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
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DISFinder & operator=(const DISFinder &);

};

}

#endif /* Ariadne5_DISFinder_H */
