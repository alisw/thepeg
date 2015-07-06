// -*- C++ -*-
#ifndef Ariadne5_DipoleBase_H
#define Ariadne5_DipoleBase_H
//
// This is the declaration of the DipoleBase class.
//

#include "CascadeBase.h"
#include "DipoleBase.fh"

namespace Ariadne5 {

using namespace ThePEG;

/**
 * DipoleBase is the base class of all Dipole classes used in the
 * Dipole Cascade Model and its extensions. Currently it does not
 * introduce any functionality.
 */
class DipoleBase: public CascadeBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  DipoleBase() {}

  /**
   * The destructor.
   */
  virtual ~DipoleBase() {}
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

  /**
   * The cutoff scale for this dipole.
   */
  Energy rhoCut() const {
    return theRhoCut;
  }

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
   * Check the integrity of this dipole. Return false if error is
   * found.
   */
  virtual bool checkIntegrity();

public:

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

protected:

  /**
   * The cutoff scale for this dipole. To be set by sub-classes on
   * construction.
   */
  Energy theRhoCut;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DipoleBase & operator=(const DipoleBase &);

};

}


#endif /* ARIADNE_DipoleBase_H */
