// -*- C++ -*-
#ifndef ARIADNE5_SpecialQQEmission_H
#define ARIADNE5_SpecialQQEmission_H
//
// This is the declaration of the SpecialQQEmission class.
//

#include "FSQQEmission.h"
#include "PseudoParton.h"
#include "Ariadne/Cascade/EmitterBase.h"

namespace Ariadne5 {

using namespace ThePEG;

/**
 *The SpecialQQEmission class contains all information about a
 * generated and performed final state gluon emission from special
 * partons.
 */
class SpecialQQEmission: public FSQQEmission {

public:
  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The only relevant constructor.
   */
  SpecialQQEmission(const EmitterBase & inmodel,
		       const DipoleBase & indipole,
		       int infl, Energy inmq, Energy2 Sin,
		       const PseudoParton & inip, const PseudoParton & inop)
    : FSQQEmission(inmodel, indipole, infl, inmq), S(Sin),
      pip(inip), pop(inop), wrem1(1.0), wrem3(1.0) {}

  /**
   * The deault constructor should not normally be used.
   */
  SpecialQQEmission() {}

  /**
   * The destructor.
   */
  virtual ~SpecialQQEmission();
  //@}

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual ClonePtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual ClonePtr fullclone() const;
  //@}


public:

  /**
   * The total squared invariant mass
   */
  Energy2 S;

  /**
   * Pseudo parton 1
   */
  PseudoParton pip;

  /**
   * Pseudo parton 3
   */
  PseudoParton pop;

  /**
   * The remnant weight for parton 1
   */
  double wrem1;
  /**
   * The remnant weight for parton 3
   */
  double wrem3;

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
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  SpecialQQEmission & operator=(const SpecialQQEmission &);

};

}

#endif /* ARIADNE_SpecialQQEmission_H */
