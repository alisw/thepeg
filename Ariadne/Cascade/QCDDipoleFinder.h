// -*- C++ -*-
#ifndef Ariadne5_QCDDipoleFinder_H
#define Ariadne5_QCDDipoleFinder_H
//
// This is the declaration of the QCDDipoleFinder class.
//

#include "ThePEG/Handlers/HandlerBase.h"
#include "QCDDipoleFinder.fh"
#include "QCDDipole.fh"
#include "DipoleState.fh"
#include "Parton.fh"

namespace Ariadne5 {

using namespace ThePEG;

/**
 * The QCDDipoleFinder class and its sub-classes are responsible for
 * identifying and introducing of QCD dipoles in the setup phase of a
 * given DipoleState.
 *
 * @see \ref QCDDipoleFinderInterfaces "The interfaces"
 * defined for QCDDipoleFinder.
 */
class QCDDipoleFinder: public HandlerBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  QCDDipoleFinder();

  /**
   * The destructor.
   */
  virtual ~QCDDipoleFinder();
  //@}

public:

  /**
   * Find, create and return the QCD dipoles in the given
   * DipoleState.
   */
  virtual vector<tQCDPtr> findDipoles(DipoleState &) const;

  /**
   * Create a QCDDipole in the given DipoleState between the Parton
   * objects \a pi and \a po carrying the anti-colour and colour
   * respectively.
   */
  tQCDPtr createDipole(tParPtr pi, tParPtr po, DipoleState & state) const;
  
public:

  /**
   * Exception class used when no consistent set of dipoles can be found.
   */
  struct QCDFinderException: public Exception {};

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
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  QCDDipoleFinder & operator=(const QCDDipoleFinder &);

};

}

#endif /* Ariadne5_QCDDipoleFinder_H */
