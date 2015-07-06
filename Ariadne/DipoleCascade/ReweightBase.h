// -*- C++ -*-
#ifndef ARIADNE_ReweightBase_H
#define ARIADNE_ReweightBase_H
//
// This is the declaration of the ReweightBase class.
//

#include "ThePEG/Handlers/HandlerBase.h"
#include "ReweightBase.fh"
#include "Dipole.h"
#include "DipoleState.fh"

namespace Ariadne {

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
 * Ariadne::CascadeHandler and will then be applied to all dipole
 * emissions, even those which are already correced by an MECorrBase
 * object. The object could possibly be used for several different
 * dipoles, and should therefore not store any information in the
 * preweigth() function to be used in the subsequent reweigth() and
 * finalVeto() calls.
 *
 * @see \ref ReweightBaseInterfaces "The interfaces"
 * defined for ReweightBase.
 */
class ReweightBase: public HandlerBase {

public:

  /** Typedef for convenience. */
  typedef Emitter::EmissionType EmissionType;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  ReweightBase() {};

  /**
   * The destructor.
   */
  virtual ~ReweightBase();
  //@}

public:

  /** @name Virtual functions to be overridden in sub-classes. */
  //@{
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
   * Return the weight associated with an emission of a parton with
   * PDG number \a id and mass emass from the given dipole, \a
   * dip. The transverse momentum squared, \a pt2, and the vector of
   * doubles that will be stored in genVar must be given. Besides the
   * \a id, an emission \a type need to be specified. Must return a
   * number between zero and one.
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
  static AbstractClassDescription<ReweightBase> initReweightBase;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ReweightBase & operator=(const ReweightBase &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of ReweightBase. */
template <>
struct BaseClassTrait<Ariadne::ReweightBase,1> {
  /** Typedef of the first base class of ReweightBase. */
  typedef HandlerBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the ReweightBase class and the shared object where it is defined. */
template <>
struct ClassTraits<Ariadne::ReweightBase>
  : public ClassTraitsBase<Ariadne::ReweightBase> {
  /** Return a platform-independent class name */
  static string className() { return "Ariadne::ReweightBase"; }
  /**
   * The name of a file containing the dynamic library where the class
   * ReweightBase is implemented. It may also include several, space-separated,
   * libraries if the class ReweightBase depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "libArCascade.so"; }
};

/** @endcond */

}

#endif /* ARIADNE_ReweightBase_H */
