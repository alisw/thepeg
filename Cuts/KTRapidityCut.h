// -*- C++ -*-
#ifndef THEPEG_KTRapidityCut_H
#define THEPEG_KTRapidityCut_H
//
// This is the declaration of the KTRapidityCut class.
//

#include "ThePEG/Cuts/OneCutBase.h"

namespace ThePEG {

/**
 * The KTRapidityCut class is a simple concrete sub-class of OneCutbase simply
 * requiring a minimum transverse momentum of any outgoing
 * particle. It is also possible to require a minimum and maximum
 * rapidity. Optionally the restrictions only apply to particles
 * matching a specific matcher object.
 *
 * @see \ref KTRapidityCutInterfaces "The interfaces"
 * defined for KTRapidityCut.
 * @see SimpleKtCut
 */
class KTRapidityCut: public OneCutBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  KTRapidityCut(Energy minKT=10*GeV) 
    : theMinKT(minKT), theMaxKT(Constants::MaxEnergy),
      theMinRapidity(-Constants::MaxRapidity),
      theMaxRapidity(Constants::MaxRapidity) {}

  /**
   * The destructor.
   */
  virtual ~KTRapidityCut();
  //@}

public:

  /** @name Overwritten virtual functions defined in the base class. */
  //@{
  /**
   * Return the minimum allowed value of the transverse momentum of an
   * outgoing parton.
   */
  virtual Energy minKT(tcPDPtr p) const;

  /**
   * Return the minimum allowed pseudo-rapidity of an outgoing parton
   * of the given type. The pseudo-rapidity is measured in the lab
   * system.
   */
  virtual double minEta(tcPDPtr p) const;

  /**
   * Return the maximum allowed pseudo-rapidity of an outgoing parton
   * of the given type. The pseudo-rapidity is measured in the lab
   * system.
   */
  virtual double maxEta(tcPDPtr p) const;

  /**
   * Return true if a particle with type \a ptype and momentum \a p
   * passes the cuts. The \a parent contains information about the
   * kinematics of the hard sub-process.
   */
  virtual bool passCuts(tcCutsPtr parent,
			tcPDPtr ptype, LorentzMomentum p) const;
  //@}

  /**
   * Describe the currently active cuts in the log file.
   */
  virtual void describe() const;

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

private:

  /**
   * Helper function used by the interface.
   */
  Energy maxKTMin() const;

  /**
   * Helper function used by the interface.
   */
  Energy minKTMax() const;

  /**
   * Helper function used by the interface.
   */
  double maxRapidityMin() const;

  /**
   * Helper function used by the interface.
   */
  double minRapidityMax() const;

  /**
   * Return the minimum allowed rapidity of an outgoing parton
   * of the given type. The rapidity is measured in the lab
   * system.
   */
  virtual double minRapidityMax(tcPDPtr p) const;

  /**
   * Return the maximum allowed rapidity of an outgoing parton
   * of the given type. The rapidity is measured in the lab
   * system.
   */
  virtual double maxRapidityMin(tcPDPtr p) const;

private:

  /**
   * The minimum allowed value of the transverse momentum of an
   * outgoing parton.
   */
  Energy theMinKT;

  /**
   * The maximum allowed value of the transverse momentum of an
   * outgoing parton.
   */
  Energy theMaxKT;

  /**
   * The minimum allowed rapidity of an outgoing parton. The
   * rapidity is measured in the lab system.
   */
  double theMinRapidity;

  /**
   * The maximum allowed rapidity of an outgoing parton. The
   * rapidity is measured in the lab system.
   */
  double theMaxRapidity;

  /**
   * If non-null only particles matching this object will be affected
   * by this cut.
   */
  PMPtr theMatcher;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<KTRapidityCut> initKTRapidityCut;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  KTRapidityCut & operator=(const KTRapidityCut &) = delete;

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of KTRapidityCut. */
template <>
struct BaseClassTrait<KTRapidityCut,1> {
  /** Typedef of the first base class of KTRapidityCut. */
  typedef OneCutBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the KTRapidityCut class and the shared object where it is defined. */
template <>
struct ClassTraits<KTRapidityCut>
  : public ClassTraitsBase<KTRapidityCut> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::KTRapidityCut"; }
  /**
   * The name of a file containing the dynamic library where the class
   * KTRapidityCut is implemented. It may also include several, space-separated,
   * libraries if the class KTRapidityCut depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "KTRapidityCut.so"; }
};

/** @endcond */

}

#endif /* THEPEG_KTRapidityCut_H */
