// -*- C++ -*-
#ifndef Ariadne5_QCDDipole_H
#define Ariadne5_QCDDipole_H
//
// This is the declaration of the QCDDipole class.
//

#include "DipoleBase.h"
#include "QCDDipole.fh"
#include "Parton.fh"
#include "DipoleState.fh"
#include "ColourIndex.h"

namespace Ariadne5 {

using namespace ThePEG;

/**
 * The QCDDipole class represents QCD dipoles between
 * <code>Parton</code>s. The QCDDipole also defines the colour flow
 * from the incoming to the outgoing Parton, where the incoming Parton
 * carries anti-colour and the outgoing Parton carries colour.
 *
 * NOTE! The dipole rest frame is defined such that the
 * anti-colour-carrying parton is along the positive z-axis.
 */
class QCDDipole: public DipoleBase {

public:

  /**
   * The QCDDipoleState is a friend.
   */
  friend class DipoleState;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  QCDDipole();

  /**
   * The destructor.
   */
  inline virtual ~QCDDipole() {};
  //@}

public:

  /**
   * Get (or create and return) the colour line corresponding
   * corresponding to this dipole, if present. Is only meaningful if
   * Parton::produceParticle() has been called for the connected
   * partons. The colour line should be properly connected to the
   * produced particle.
   */
  virtual ColinePtr colourLine() const;

  /** @name Simple access functions. */
  //@{
  /**
   * The anti-colour-carrying parton in this dipole.
   */
  inline tParPtr iPart() const {
    return theIPart;
  }

  /**
   * The previous dipole in the string. Returns the dipole whose
   * oPart() is this iPart() if any.
   */
  tQCDPtr prev() const {
    return thePrev;
  }

   /**
   * Set the previous dipole in the string.
   */
  void prev(tQCDPtr d) {
    thePrev = d;
  }

 /**
   * The colour-carrying parton in this dipole.
   */
  inline tParPtr oPart() const {
    return theOPart;
  }

  /**
   * The next dipole in the string. Returns the dipole whose
   * iPart() is this oPart() if any.
   */
  tQCDPtr next() const {
    return theNext;
  }

  /**
   * Set the next dipole in the string.
   */
  void next(tQCDPtr d) {
    theNext = d;
  }

  /**
   * Return all partons in this string (or string piece if it is in a junction).
   */
  vector<tParPtr> string() const;

  /**
   * The colour index of this dipole.
   */
  inline const ColourIndex & colourIndex() const {
    return theColourIndex;
  }

  /**
   * Set the anti-colour-carrying parton in this dipole.
   */
  inline void iPart(tParPtr x) {
    theIPart = x;
  }

  /**
   * Set the colour-carrying parton in this dipole.
   */
  inline void oPart(tParPtr x) {
    theOPart = x;
  }

  /**
   * Set the colour index of this dipole.
   */
  inline void colourIndex(const ColourIndex & x) {
    theColourIndex = x;
  }

  /**
   * Generte a colour index of this dipole.
   */
  void generateColourIndex();

  /**
   * Return the squared invariant mass of this dipole.
   */
  Energy2 sdip() const;

  // *** ATTENTION *** Do we need this?
  /**
   * The particle type if this dipole comes directly from a decay of a
   * resonance. If not, null is returned.
   */
  inline tcPDPtr resonance() const {
    return theResonance;
  }

  /**
   * The particle type if this dipole comes directly from a decay of a
   * resonance.
   */
  inline void resonance(tcPDPtr x) {
    theResonance = x;
  }
  //@}

protected:

  /** @name Functions relating to the DipoleState to which this belongs. */
  //@{
  /**
   * Return a simple clone of this object. Should be implemented as
   * <code>return new_ptr(*this);</code> by a derived class.
   */
  virtual ClonePtr clone() const;

  /**
   * Fill the provided set with all pointers to CloneBase objects used
   * in this object.
   */
  virtual void fillReferences(CloneSet &) const;

  /**
   * Rebind pointers to other CloneBase objects. Called after a number
   * of interconnected CloneBase objects have been cloned, so that
   * the cloned objects will refer to the cloned copies afterwards.
   *
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   */
  virtual void rebind(const TranslationMap & trans);
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

private:

  /**
   * The anti-colour-carrying parton in this dipole.
   */
  tParPtr theIPart;

  /**
   * The colour-carrying parton in this dipole.
   */
  tParPtr theOPart;

  /**
   * The next dipole in the string.
   */
  tQCDPtr theNext;

  /**
   * The previous dipole in the string.
   */
  tQCDPtr thePrev;

  /**
   * The colour index of this dipole. 0 means no index has been assigned.
   */
  ColourIndex theColourIndex;

  /**
   * The particle type if this dipole comes directly from a decay of a
   * resonance.
   */
  tcPDPtr theResonance;

public:

  /**
   * Print out debugging information on std::cerr.
   */
  virtual void debugme() const;

  /**
   * Check integrity of the emitter. Return false if error is found.
   */
  virtual bool checkIntegrity();

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  QCDDipole & operator=(const QCDDipole &);

};

}

#endif /* Ariadne5_QCDDipole_H */
