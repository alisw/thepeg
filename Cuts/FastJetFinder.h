// -*- C++ -*-
//
// FastJetFinder.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
// Copyright (C) 2009-2017 Simon Platzer
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_FastJetFinder_H
#define THEPEG_FastJetFinder_H
//
// This is the declaration of the FastJetFinder class.
//

#include "ThePEG/Cuts/JetFinder.h"

namespace ThePEG {

/**
 * FastJetFinder implements the class of longitudinally invariant kt
 * jet clustering algorithms.
 *
 * @see \ref FastJetFinderInterfaces "The interfaces"
 * defined for FastJetFinder.
 */
class FastJetFinder: public JetFinder {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  FastJetFinder();

  /**
   * The destructor.
   */
  virtual ~FastJetFinder();
  //@}

public:

  /**
   * Perform jet clustering on the given outgoing particles.
   * Optionally, information on the incoming particles is provided.
   * Return true, if a clustering has been performed.
   */
  virtual bool cluster(tcPDVector & ptype, vector<LorentzMomentum> & p,
		       tcCutsPtr parent, tcPDPtr t1 = tcPDPtr(),
		       tcPDPtr t2 = tcPDPtr()) const;

  /**
   * Describe this jet fined.
   */
  virtual void describe() const;

public:
  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object. Called in the read phase.
   */
  virtual void doinit();

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();
  //@}

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
   * The resolution cut.
   */
  Energy2 theDCut;

  /**
   * The `cone radius' R.
   */
  double theConeRadius;

  /**
   * The possible variants.
   */
  enum variants {
    kt = 1,
    CA = 2,
    antiKt = 3,
    sphericalKt = 4,
    sphericalCA = 5,
    sphericalAntiKt = 6
  };

  /**
   * The variant.
   */
  int theVariant;

  /**
   * The possible modes.
   */
  enum modes {
    inclusive = 1,
    exclusive = 2
  };

  /**
   * The mode.
   */
  int theMode;

  /**
   * The possible recombination schemes.
   */
  enum recombinations {
    recoPt = 1,
    recoE = 2
  };

  /**
   * The recombination scheme
   */
  int theRecombination;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  FastJetFinder & operator=(const FastJetFinder &);

};

}

#endif /* THEPEG_FastJetFinder_H */
