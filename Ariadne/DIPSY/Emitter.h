// -*- C++ -*-
#ifndef DIPSY_Emitter_H
#define DIPSY_Emitter_H
//
// This is the declaration of the Emitter class.
//

#define THEPEG_NEW_CLASS_DESCRIPTION

#include "ThePEG/Handlers/HandlerBase.h"
#include "Emitter.fh"
#include "Dipole.fh"
#include "Parton.h"

namespace DIPSY {

using namespace ThePEG;

/**
 * The Emitter class is responsible for generating and performing
 * emissions from dipoles. This base class does the standard thing. Any
 * non-standard thing must be implemented in a sub-class overriding
 * the generate() and/or emit() functions.
 *
 * @see \ref EmitterInterfaces "The interfaces"
 * defined for Emitter.
 */
class Emitter: public HandlerBase {

public:

  /** Convenietn typedef. */
  typedef Parton::Point Point;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline Emitter()
    : testingPS(false), fixY(0.0), thePSInflation(1.0), thePMinusOrdering(1.0),
      thePTScale(1.0), theRScale(1.0*InvGeV), theRMax(0.0*InvGeV),
      theBothOrderedEvo(true), theBothOrderedInt(true), theBothOrderedFS(true),
      theRangeMode(0),theMinusOrderingMode(0), thestMode(false) {}

  /**
   * The destructor.
   */
  virtual ~Emitter();
  //@}

public:

  /** @name Virtual functions which may be overridden by sub-classes. */
  //@{

  /**
   * Return the running coupling for the given size (scale). /CF
   */
  double alphaS(InvEnergy r) const;

  /**
   * if pT is 1/r or 2/r.
   */
  inline double pTScale() const {
    return thePTScale;
  }

  /**
   * Alpha bar = alphas*Nc/pi /CF
   */
  inline double alphaBar(InvEnergy r) const {
    return alphaS(r)*3.0/M_PI;
  }

  /**
   * Sets the shape of the overestimate in the genreation.
   * The cross section does not (should not) depend on rScale. /CF
   */
  inline void rScale(InvEnergy rScale) {
    theRScale = rScale;
  }

  /**
   * Get bothordered.
   */
  inline bool bothOrderedEvo() const {
    return theBothOrderedEvo;
  };
  inline bool bothOrderedInt() const {
    return theBothOrderedInt;
  };
  inline bool bothOrderedFS() const {
    return theBothOrderedFS;
  };

  /**
   * Set PSInfaltion
   */
  inline void PSInflation(double PSInflation) {
    thePSInflation = PSInflation;
  }

  /**
   * Get PSInfaltion
   */
  inline double PSInflation() const {
    return thePSInflation;
  }

  /**
   * Set PMinusOrdering
   */
  inline void PMinusOrdering(double x) {
    thePMinusOrdering = x;
  }

  /**
   * Get PMinusOrdering
   */
  inline double PMinusOrdering() const {
    return thePMinusOrdering;
  }

  /**
   * get the rangeMode
   **/
  inline int rangeMode() const {
    return theRangeMode;
  }

  /**
   * get the minusOrderingMode
   **/
  inline int minusOrderingMode() const {
    return theMinusOrderingMode;
  }

  /**
   * The confinement scale.
   */
  InvEnergy rMax() const;

  /**
   * Generate a possible emission or a swing from a given dipole in the
   * given rapidity interval [\a miny,\a maxy].
   */
  virtual void generate(Dipole & dipole, double miny, double maxy) const;

  /**
   * Generate a possible emission or a swing from a given dipole using
   * shadows in the given rapidity interval [\a miny,\a maxy].
   */
  virtual void generateWithShadows(Dipole & dipole, double miny, double maxy) const;

  /**
   * Perform the emission previously generated for the given \a
   * dipole. If no emission has been generated a runtime_error is
   * thrown.
   */
   virtual void emit(Dipole & dipole) const;
  /**
   * Perform the emission previously generated for the given \a
   * dipole using shadows. If no emission has been generated a runtime_error is
   * thrown.
   */
   virtual void emitWithShadows(Dipole & dipole) const;
  //@}

protected:

  /**
   * Internal function.
   */
  virtual bool OEVeto(DipolePtr dip, double y0, Point p) const;

  /**
   * Internal function.
   */
  virtual bool OEVeto(DipolePtr dip, tSPartonPtr sp1, tSPartonPtr sp2, double y0, Point p) const;

  /**
   * Internal function.
   */
  virtual bool YVeto(double y, DipolePtr dip, double Coe, double rateOE) const;

  /**
   * Internal function.
   */
  virtual bool YVeto(double y, DipolePtr dip, tSPartonPtr sp1, tSPartonPtr sp2,
		     double Coe, double rateOE) const;

  /**
   * Internal function.
   */
  virtual double generateY(DipolePtr dip, double ymin, double ymax) const;

  /**
   * Internal function.
   */
  virtual double generateY(DipolePtr dip, tSPartonPtr sp1, tSPartonPtr sp2,
			   double ymin, double ymax) const;

  /**
   * Internal function.
   */
  virtual Point generateXT(DipolePtr dip, double y0) const;

  /**
   * Internal function.
   */
  virtual Point generateXT(DipolePtr dip, tSPartonPtr sp1, tSPartonPtr sp2, double y0) const;

public:
  /**
   * Debugging and testing
   */
  mutable bool testingPS;
  mutable double fixY;

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

protected:  //@{

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr clone() const {
    return new_ptr(*this);
  }

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const {
    return new_ptr(*this);
  }
  //@}


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).

protected:

  /**
   * How much P+- ordering should be required. 1 is nomral ordering
   * high values are no ordering, low number are strong ordering.
   */
  double thePSInflation;

  /**
   * Controls the strength of the p- ordering in the evolution.
   */
  double thePMinusOrdering;

  /**
   * if pt is 1/r or 2/r.
   */
  double thePTScale;

  /**
   * the scale in the distribution of transverse r. /CF
   */
  InvEnergy theRScale;

  /**
   * The confinement scale.
   */
  InvEnergy theRMax;

  /**
   * If the plus/minus ordering should be full strenght
   * for both parents in teh cascade, interaction and final state.
   */
  bool theBothOrderedEvo;
  bool theBothOrderedInt;
  bool theBothOrderedFS;

  /**
   * How the range is determined.
   **/
  int theRangeMode;

  /**
   * How the minus ordering is made in the cascade.
   **/
  int theMinusOrderingMode;

  /**
   * If test-modification to evolution should be on. /CF
   */
  bool thestMode;

public:

  /**
   * Exception class.
   */
  struct EmitterException: public Exception {};

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  Emitter & operator=(const Emitter &);

};

}

#ifndef THEPEG_NEW_CLASS_DESCRIPTION

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of Emitter. */
template <>
struct BaseClassTrait<DIPSY::Emitter,1> {
  /** Typedef of the first base class of Emitter. */
  typedef HandlerBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the Emitter class and the shared object where it is defined. */
template <>
struct ClassTraits<DIPSY::Emitter>
  : public ClassTraitsBase<DIPSY::Emitter> {
  /** Return a platform-independent class name */
  static string className() { return "DIPSY::Emitter"; }
  /**
   * The name of a file containing the dynamic library where the class
   * Emitter is implemented. It may also include several, space-separated,
   * libraries if the class Emitter depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "libAriadne5.so libDIPSY.so"; }
};

/** @endcond */

}

#endif

#endif /* DIPSY_Emitter_H */
