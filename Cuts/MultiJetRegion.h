// -*- C++ -*-
//
// MultiJetRegion.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
// Copyright (C) 2009-2019 Simon Platzer
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_MultiJetRegion_H
#define ThePEG_MultiJetRegion_H
//
// This is the declaration of the MultiJetRegion class.
//

#include "ThePEG/Cuts/JetRegion.h"

namespace ThePEG {

/**
 * MultiJetRegion implements pairwise constraints on jets matching several jet regions.
 *
 * @see JetRegion
 * @see JetCuts
 *
 * @see \ref MultiJetRegionInterfaces "The interfaces"
 * defined for MultiJetRegion.
 */
class MultiJetRegion: public HandlerBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  MultiJetRegion();

  /**
   * The destructor.
   */
  virtual ~MultiJetRegion();
  //@}

public:

  /**
   * Return the jet regions to act on.
   */
  const vector<Ptr<JetRegion>::ptr>& regions() const { return theRegions; }

  /**
   * Return the minimum jet-jet invariant mass.
   */
  Energy massMin() const { return theMassMin; }

  /**
   * Return the maximum jet-jet invariant mass.
   */
  Energy massMax() const { return theMassMax; }

  /**
   * Return the minimum jet-jet lego-plot separation.
   */
  double deltaRMin() const { return theDeltaRMin; }

  /**
   * Return the maximum jet-jet lego-plot separation.
   */
  double deltaRMax() const { return theDeltaRMax; }

  /**
   * Return the minimum jet-jet rapidity separation.
   */
  double deltaYMin() const { return theDeltaYMin; }

  /**
   * Return the maximum jet-jet rapidity separation.
   */
  double deltaYMax() const { return theDeltaYMax; }

  /**
   * Return the cut weight encountered from the last call to matches()
   */
  double cutWeight() const { return theCutWeight; }

public:

  /**
   * Describe the currently active cuts in the log file.
   */
  virtual void describe() const;

  /**
   * Return true, if the requirements on the jet regions are fullfilled.
   */
  virtual bool matches();

  /**
   * Return true, if the requirements on two jet regions are fullfilled.
   */
  virtual bool matches(int i, int j);

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
   * The jet regions to act on.
   */
  vector<Ptr<JetRegion>::ptr> theRegions;

  /**
   * The minimum jet-jet invariant mass.
   */
  Energy theMassMin;

  /**
   * The maximum jet-jet invariant mass.
   */
  Energy theMassMax;

  /**
   * The minimum jet-jet lego-plot separation.
   */
  double theDeltaRMin;

  /**
   * The maximum jet-jet lego-plot separation.
   */
  double theDeltaRMax;

  /**
   * The minimum jet-jet rapidity separation.
   */
  double theDeltaYMin;

  /**
   * The maximum jet-jet rapidity separation.
   */
  double theDeltaYMax;

  /**
   * The cut weight encountered from the last call to matches()
   */
  double theCutWeight;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MultiJetRegion & operator=(const MultiJetRegion &) = delete;

};

}

#endif /* ThePEG_MultiJetRegion_H */
