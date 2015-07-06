// -*- C++ -*-
#ifndef Ariadne5_ReweightBase_H
#define Ariadne5_ReweightBase_H
//
// This is the declaration of the ReweightBase class.
//

#include "ThePEG/Handlers/HandlerBase.h"
#include "ReweightBase.fh"
#include "DipoleBase.fh"
#include "Parton.fh"
#include "Emission.fh"

namespace Ariadne5 {

using namespace ThePEG;

/**
 * ReweightBase is the base class for implementing different ways of
 * reweighting the dipole emissions in Ariadne. There are three
 * virtual functions which may be overridden in concrete sub-classes:
 * preweight() should return a factor multiplying the standard dipole
 * emission before reweighting and reweight() should return the actual
 * reweight. Alternatively the finalVeto() function can be used to
 * accept/reject an already performed dipole emission. Object of a
 * concrete sub-class can inserted to a list in
 * Ariadne5::CascadeHandler and will then be applied to all dipole
 * emissions. The object could possibly be used for several different
 * dipoles, and should therefore not store any information in the
 * preweigth() function to be used in the subsequent reweigth() and
 * finalVeto() calls.
 *
 * @see \ref ReweightBaseInterfaces "The interfaces"
 * defined for ReweightBase.
 */
class ReweightBase: public HandlerBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor. \a mayveto must be set to true by
   * subclasses which implements the finalVeto function.
   */
  ReweightBase(bool mayveto = false) : mayVeto(mayveto) {};

  /**
   * The destructor.
   */
  virtual ~ReweightBase();
  //@}

public:

  /** @name Virtual functions to be overridden in sub-classes. */
  //@{
  /**
   * Return a factor to multiply the basic emission emission
   * probability for a given \a emission to ensure that a subsequent
   * call to reweight() will give a number less than unity.
   */
  virtual double preweight(const Emission & emission) const;
 
  /**
   * Return the weight associated with an \a emission. Must return a
   * number between zero and one.
   */
  virtual double reweight(const Emission & emission) const;

  /**
   * In addition to the reweight function a final hit/miss veto may be
   * given after the \a emission has been performed. Will only be
   * called if hasFinalVeto() returns true.
   */
  virtual bool finalVeto(const Emission & emission) const;

  /**
   * Indicate that this object may veto a performed emission.
   */
  bool hasFinalVeto() const {
    return mayVeto;
  }
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


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).

private:

  /**
   * Indicate that this object may veto a performed emission.
   */
  bool mayVeto;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ReweightBase & operator=(const ReweightBase &);

};

}

#endif /* Ariadne5_ReweightBase_H */
