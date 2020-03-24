// -*- C++ -*-
//
// JetRegion.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
// Copyright (C) 2009-2019 Simon Platzer
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_JetRegion_H
#define ThePEG_JetRegion_H
//
// This is the declaration of the JetRegion class.
//

#include "ThePEG/Handlers/HandlerBase.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Cuts/Cuts.h"

namespace ThePEG {

/**
 * JetRegion implements the requirement of finding a jet inside a
 * given range of transverse momenta, and (pseudo-)rapidity.
 *
 * @see JetPairRegion
 * @see JetCuts
 *
 * @see \ref JetRegionInterfaces "The interfaces"
 * defined for JetRegion.
 */
class JetRegion: public HandlerBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  JetRegion();

  /**
   * The destructor.
   */
  virtual ~JetRegion();
  //@}

public:

  /**
   * Return the minimum pt.
   */
  Energy ptMin() const { return thePtMin; }

  /**
   * Return the maximum pt.
   */
  Energy ptMax() const { return thePtMax; }

  /**
   * Return the rapidity ranges.
   */
  const vector<pair<double,double> >& yRanges() const { return theYRanges; }

  /**
   * Return the jets accepted by this region (with respect to the
   * ordering imposed by the JetCuts object). If empty, any jet will
   * be accepted.
   */
  const vector<int>& accepts() const { return theAccepts; }

  /**
   * Return true, if this jet region is fuzzy
   */
  bool fuzzy() const { return theFuzzy; }

  /**
   * Return the cut weight encountered from the last call to matches()
   */
  double cutWeight() const { return theCutWeight; }

  /**
   * Perform a (potentially) fuzzy check on energy-type quantities
   */
  bool lessThanEnergy(Energy a, Energy b, double& weight) const;

  /**
   * Perform a (potentially) fuzzy check on angular-type quantities
   */
  bool lessThanRapidity(double a, double b, double& weight) const;

public:

  /**
   * Describe the currently active cuts in the log file.
   */
  virtual void describe() const;

  /**
   * Return true, if the given jet matches this region.
   */
  virtual bool matches(tcCutsPtr parent, int n, const LorentzMomentum& p,
		       double yHat = 0.0);

  /**
   * Return true, if this region matched a jet in the last call to matches().
   */
  bool didMatch() { return theDidMatch; }

  /**
   * Reset this region to act  on a new event.
   */
  virtual void reset() { theDidMatch = false; }

  /**
   * Return the number of the last jet matching this region.
   */
  int lastNumber() const { return theLastNumber; }

  /**
   * Return the momentum of the last jet matching this region.
   */
  const LorentzMomentum& lastMomentum() const { return theLastMomentum; }

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
   * Command to insert a rapidity range
   */
  string doYRange(string);

  /**
   * The minimum pt.
   */
  Energy thePtMin;

  /**
   * The maximum pt.
   */
  Energy thePtMax;

  /**
   * The rapidity ranges.
   */
  vector<pair<double,double> > theYRanges;

  /**
   * The jets accepted by this region (with respect to the ordering
   * imposed by the JetCuts object). If empty, any jet will be
   * accepted.
   */
  vector<int> theAccepts;

  /**
   * True, if this region matched a jet in the last call to matches().
   */
  bool theDidMatch;

  /**
   * The number of the last jet matching this region.
   */
  int theLastNumber;

  /**
   * Return the momentum of the last jet matching this region.
   */
  LorentzMomentum theLastMomentum;

  /**
   * True if this region is fuzzy
   */
  bool theFuzzy;

  /**
   * The cut weight encountered from the last call to matches()
   */
  double theCutWeight;

  /**
   * The smearing width for the pt or mass cuts, if fuzzy
   */
  Energy theEnergyCutWidth;

  /**
   * The smearing width for the rapidity cut, if fuzzy
   */
  double theRapidityCutWidth;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  JetRegion & operator=(const JetRegion &) = delete;

};

}

#endif /* ThePEG_JetRegion_H */
