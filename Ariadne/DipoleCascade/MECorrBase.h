// -*- C++ -*-
#ifndef ARIADNE_MECorrBase_H
#define ARIADNE_MECorrBase_H
//
// This is the declaration of the MECorrBase class.
//

#include "ThePEG/Handlers/HandlerBase.h"
#include "MECorrBase.fh"
#include "Dipole.h"
#include "DipoleState.fh"

namespace Ariadne {

using namespace ThePEG;

/**
 * MECorrBase is the base class for implementing leading order
 * matrix-element corrections to the dipole emissions in
 * Ariadne. There are three virtual functions which may be overridden
 * in concrete sub-classes: canHandle() should check if the
 * matrix-element correction is applicable to a given Dipole in a
 * given DipoleState; preweight() should return a factor multiplying
 * the standard dipole emission before reweighting; and reweight()
 * should return the actual matrix-element weight. Alternatively the
 * finalVeto() function can be used to accept/reject an already
 * performed dipole emission. An object of a concrete sub-class can
 * inserted to a list in Ariadne::CascadeHandler. The object could
 * possibly be used for several different dipoles, and should
 * therefore not store any information in the preweigth() function to
 * be used in the subsequent reweigth() and finalVeto() calls.
 *
 * @see \ref MECorrBaseInterfaces "The interfaces"
 * defined for MECorrBase.
 */
class MECorrBase: public HandlerBase {

public:

  /** Typedef for convenience. */
  typedef Emitter::EmissionType EmissionType;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline MECorrBase() {};

  /**
   * The destructor.
   */
  virtual ~MECorrBase();
  //@}

public:

  /** @name Virtual functions to be overridden in sub-classes. */
  //@{
  /**
   * Check if the matrix-element correction can be applied to the
   * given \a dipole in the given dipole \a state.
   * @return true if this correction is applicable, false otherwise.
   */
  virtual bool
  canHandle(tcEmiPtr dipole, tcDipoleStatePtr state) const = 0;

  /**
   * Return a factor to multiply the basic dipole emission probability
   * for an emission of a parton with PDG number \a id from the given
   * \a dipole to ensure that a subsequent call to reweight() will
   * give a number less than unity. Besides the \a id, an emission \a
   * type need to be specified.
   */
  virtual double
  preweight(tcEmiPtr dipole, tcDipoleStatePtr state, long id,
	    const EmissionType & type) const;
 
  /**
   * Return the matrix-element correction weight associated with an
   * emission of a parton with PDG number \a id and mass emass from
   * the given dipole, \a dip. The transverse momentum squared, \a
   * pt2, and the vector of doubles that will be stored in genVar must
   * be given. Besides the \a id, an emission \a type need to be
   * specified.
   * 
   */
  virtual double
  reweight(tcEmiPtr dip, tcDipoleStatePtr state, long id,
      Energy2 pt2, vector<double> & genVar,
      const EmissionType & type) const;

  /**
   * In addition to the reweight function a final hit/miss veto may be
   * given after the emission has been performed. The arguments are
   * the original \a dipole, the final dipole \a state and the emitted
   * \a parton. Also the scale (\a pt2) and \a type of the emission
   * must be supplied.
   */
  virtual bool
  finalVeto(tcEmiPtr dipole, tcDipoleStatePtr state, tcParPtr parton,
      Energy2 pt2, const EmissionType & type) const;
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
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class with persistent data.
   */
  static AbstractClassDescription<MECorrBase> initMECorrBase;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MECorrBase & operator=(const MECorrBase &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MECorrBase. */
template <>
struct BaseClassTrait<Ariadne::MECorrBase,1> {
  /** Typedef of the first base class of MECorrBase. */
  typedef HandlerBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MECorrBase class and the shared object where it is defined. */
template <>
struct ClassTraits<Ariadne::MECorrBase>
  : public ClassTraitsBase<Ariadne::MECorrBase> {
  /** Return a platform-independent class name */
  static string className() { return "Ariadne::MECorrBase"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MECorrBase is implemented. It may also include several, space-separated,
   * libraries if the class MECorrBase depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "libArCascade.so"; }
};

/** @endcond */

}

#endif /* ARIADNE_MECorrBase_H */
