// -*- C++ -*-
#ifndef ARIADNE_DYCorr_H
#define ARIADNE_DYCorr_H
//
// This is the declaration of the DYCorr class.
//

#include "Ariadne/DipoleCascade/MECorrBase.h"
#include "DYCorr.fh"

namespace Ariadne {

using namespace ThePEG;

/**
 * DYCorr is a class implementing the first order matrix element
 * corrections to the W/Z/gamma production in hadron collisions.
 *
 * @see \ref DYCorrInterfaces "The interfaces"
 * defined for DYCorr.
 */
class DYCorr: public MECorrBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline DYCorr();

  /**
   * The copy constructor.
   */
  inline DYCorr(const DYCorr &);

  /**
   * The destructor.
   */
  virtual ~DYCorr();
  //@}

public:

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

public:

  /**
   * Check if the matrix-element correction can be applied to the
   * given \a dipole in the given dipole \a state.
   * @return true if this correction is applicable, false otherwise.
   */
  virtual bool
  canHandle(tcEmiPtr dipole, tcDipoleStatePtr state) const;

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

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
  //@}


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).


private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an concrete class without persistent data.
   */
  static NoPIOClassDescription<DYCorr> initDYCorr;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  DYCorr & operator=(const DYCorr &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of DYCorr. */
template <>
struct BaseClassTrait<Ariadne::DYCorr,1> {
  /** Typedef of the first base class of DYCorr. */
  typedef Ariadne::MECorrBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the DYCorr class and the shared object where it is defined. */
template <>
struct ClassTraits<Ariadne::DYCorr>
  : public ClassTraitsBase<Ariadne::DYCorr> {
  /** Return a platform-independent class name */
  static string className() { return "Ariadne::DYCorr"; }
  /**
   * The name of a file containing the dynamic library where the class
   * DYCorr is implemented. It may also include several, space-separated,
   * libraries if the class DYCorr depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "libArCascade.so"; }
};

/** @endcond */

}

#include "DYCorr.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "DYCorr.tcc"
#endif

#endif /* ARIADNE_DYCorr_H */
