// -*- C++ -*-
#ifndef DIPSY_ImpactParameterGenerator_H
#define DIPSY_ImpactParameterGenerator_H
//
// This is the declaration of the ImpactParameterGenerator class.
//

#include "ThePEG/Handlers/HandlerBase.h"
#include "ImpactParameterGenerator.fh"
#include "ImpactParameters.h"
#include "ThePEG/Repository/UseRandom.h"

namespace DIPSY {

using namespace ThePEG;

/**
 * The objective of the ImpactParameterGenerator class is to generate
 * ImpactParameters objects to be used when two
 * <code>DipoleState</code>s have been generated w.r.t.  origo in the
 * transverse to displace and rotate one of them before
 * collision. This base class will generate impact parameters
 * according to a Gaussian distribution, and the weigth in the
 * produced ImpactParameters objects is set accordingly. Sub-classes
 * may override the generate() function to use any distribution as
 * long as the weight is set accordingly.
 *
 * @see \ref ImpactParameterGeneratorInterfaces "The interfaces"
 * defined for ImpactParameterGenerator.
 */
class ImpactParameterGenerator: public HandlerBase {

public:

  /**
   * Use the same Point class as the Parton.
   */
  typedef Parton::Point Point;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  ImpactParameterGenerator();

  /**
   * The copy constructor.
   */
  ImpactParameterGenerator(const ImpactParameterGenerator &);

  /**
   * The destructor.
   */
  virtual ~ImpactParameterGenerator();
  //@}

public:

  /**
   * Generate an ImpactParameters object.
   */
  virtual ImpactParameters generate(double seed = UseRandom::rnd()) const;

  /**
   * Generate an ImpactParameters object, with
   * the two points at the centre of the distribution.
   */
  virtual ImpactParameters generateDynamic(vector<pair<Parton::Point, InvEnergy> > points1,
					   vector<pair<Parton::Point, InvEnergy> > points2) const;

  /**
   * The width of the generated distribution.
   */
  InvEnergy width() const {
    return theWidth;
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
   * The width of the generated distribution.
   */
  InvEnergy theWidth;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ImpactParameterGenerator & operator=(const ImpactParameterGenerator &);

};

}

#endif /* DIPSY_ImpactParameterGenerator_H */
