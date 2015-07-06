// -*- C++ -*-
#ifndef DIPSY_NucleusData_H
#define DIPSY_NucleusData_H
//
// This is the declaration of the NucleusData class.
//

#include "ThePEG/PDT/ParticleData.h"

namespace DIPSY {

using namespace ThePEG;

/**
 * Here is the documentation of the NucleusData class.
 *
 * @see \ref NucleusDataInterfaces "The interfaces"
 * defined for NucleusData.
 */
class NucleusData: public ThePEG::ParticleData {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  NucleusData();

  /**
   * The destructor.
   */
  virtual ~NucleusData();
  //@}

  /** @name The Create methods are special interfaces for ParticleData
      classes. */
  //@{
  /**
   * Create a Particle which is its own anti-particle.
   */
  static PDPtr Create(long newId, string newPDGName);

  /**
   * Create a particle - anti particle pair.
   */
  static PDPair Create(long newId, string newPDGName, string newAntiPDGName);
  //@}

public:

  /**
   * Get the atomic number for this nucleus.
   */
  inline unsigned int A() const {
    return theA;
  }

  /**
   * Get the number of protons in this nucleus.
   */
  inline int Z() const {
    return theZ;
  }

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

  /**
   * Protected constructor only to be used by subclasses or by the
   * Create method.
   */
  NucleusData(long newId, string newPDGName);

  /**
   * Read setup info from a standard stream. The information that must
   * be supplied is the same as for ParticleData::readSetup with an
   * additional constituent mass (in GeV) added in the end.
   */
  virtual void readSetup(istream & is);

  /**
   * ParticleData clone method
   */
  virtual PDPtr pdclone() const;

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
   * The atomic number for this nucleus.
   */
  unsigned int theA;

  /**
   * The number of protons in this nucleus.
   */
  int theZ;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  NucleusData & operator=(const NucleusData &);

};

}

#endif /* DIPSY_NucleusData_H */
