// -*- C++ -*-
#ifndef PYTHIA7_Hadron_H
#define PYTHIA7_Hadron_H
// This is the declaration of the Hadron class.

#include "FragConfig.h"
// #include "Hadron.fh"
// #include "Hadron.xh"
#include "Oriented.h"
#include "EndPoint.h"

namespace Pythia7 {

/**
 * Hadron describe the minimal information of a hadron created during
 * the String Fragmentation process.  Each hadron will produce the
 * corresponding full Particle objects at the end of a succefull
 * generation.
 *
 */
class Hadron {

public:

  /** @name Standard constructors, assignment and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline Hadron();

  /**
   * The copy constructor.
   */
  inline Hadron(const Hadron &);

  /**
   * The destructor.
   */
  inline ~Hadron();

  /**
   * Assignment operator.
   */
  inline Hadron & operator=(const Hadron &);
  //@}

  /**
   * Return a pointer the corresponding Particle object .
   */
  inline PPtr createParticle();

  /**
   * Set the Hadron transverse mass from two endpoints.
   */
  inline void mT2(const EndPoint& , const EndPoint&);

  /**
   * Get the Hadron transverse mass squared. 
   */
  inline Energy2 mT2() const; 

  /**
   * Get the Hadron transverse mass. 
   */
  inline Energy mT() const;

  /**
   * get Particle mass squared.
   */
  inline Energy2 m2() const;

  /**
   * get Particle mass.
   */
  inline Energy mass() const;

  /**
   * get Particle energy.
   */
  inline Energy e() const;

  /**
   * Set the final momentum of the Hadron after the stepping process
   * is done.
   */
  inline void storeMomentum();

  /**
   * Set the Particle type.
   */
  inline void PData(tcPDPtr);

  /** Pointer to the actual Particle */
  PPtr theParticle;

  /** Pointer to the corresponding ParticleData object. */
  cPDPtr  Data;

  /** The transverse mass squared. */
  Energy2 mt2;

  /** The forward light-cone fractions. */
  double Xf;

  /** The backward light-cone fractions. */
  double Xb;

  /** The momentum */
  LorentzMomentum P;

  /** The mass. */
  Energy theParticleMass;

  /** From which side was this particle produced? */
  int ProductionSide;
  
};

}


#include "Hadron.icc"
#ifndef PYTHIA7_TEMPLATES_IN_CC_FILE
// #include "Hadron.tcc"
#endif

#endif /* PYTHIA7_Hadron_H */
