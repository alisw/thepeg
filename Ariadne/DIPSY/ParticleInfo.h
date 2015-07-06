// -*- C++ -*-
#ifndef DIPSY_ParticleInfo_H
#define DIPSY_ParticleInfo_H
//
// This is the declaration of the ParticleInfo class.
//

#include "ThePEG/EventRecord/EventInfoBase.h"
#include "ThePEG/EventRecord/Particle.h"
#include "Parton.h"

namespace DIPSY {

using namespace ThePEG;

/**
 * Here is the documentation of the ParticleInfo class.
 */
class ParticleInfo: public ThePEG::EventInfoBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  ParticleInfo(tcPartonPtr p = tcPartonPtr());

  /**
   * The destructor.
   */
  virtual ~ParticleInfo();
  //@}

public:

  /**
   * Return the parton from which the particle was created.
   */
  tcPartonPtr parton() const {
    return theParton;
  }

  /**
   * Return the parton from which the given particle was created.
   */
  static tcPartonPtr getParton(const Particle & particle);

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
   * The DIPSY::Parton from which the particle was created.
   */
  cPartonPtr theParton;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  ParticleInfo & operator=(const ParticleInfo &);

};

}

#endif /* DIPSY_ParticleInfo_H */
