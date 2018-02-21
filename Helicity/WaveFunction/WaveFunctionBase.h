// -*- C++ -*-
//
// WaveFunctionBase.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2017 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_WaveFunctionBase_H
#define ThePEG_WaveFunctionBase_H
//
// This is the declaration of the WaveFunctionBase class.

#include <ThePEG/Vectors/Lorentz5Vector.h>
#include <ThePEG/Vectors/LorentzVector.h>
#include <ThePEG/PDT/ParticleData.h>
#include <ThePEG/Helicity/HelicityDefinitions.h>

namespace ThePEG {

namespace Helicity {



/** \ingroup Helicity
 *  Definition of the enumerated values used for the direction of the 
 *  particles in the calculation of the wavefunction.
 */
enum Direction 
{
  incoming, /**< An incoming particle. */
  outgoing, /**< An outgoing particle. */
  intermediate /**< An intermediate particle. */
};

/** \ingroup Helicity
 *  \author Peter Richardson
 *
 * This class is the base class for all wavefunctions for use in helicity amplitude
 * calculations. The general approach is to use a similar philosophy 
 * to the FORTRAN HELAS code but with additional structure.
 *
 * This class contains the storage of the particle type and 5-momentum 
 * and methods to set/access this information.
 *
 * The methods for the wavefunction itself will be implemented in the classes
 * derived from this one for the specific spin type, for example scalar, spinor,
 * vector and tensor. 
 *
 *  @see ScalarWaveFunction
 *  @see SpinorWaveFunction
 *  @see SpinorBarWaveFunction
 *  @see VectorWaveFunction
 *  @see RSSpinorWaveFunction
 *  @see RSSpinorBarWaveFunction
 *  @see TensorWaveFunction
 */
class WaveFunctionBase{

public:

  /// Constructors
  //@{
  /**
   * Default constructor
   */
  WaveFunctionBase() 
    : _particle(), _momentum(), _dir(intermediate) 
  {}

  /**
   * 
   */
  WaveFunctionBase(const Lorentz5Momentum & p,
		   tcPDPtr pd, Direction dir = intermediate) 
    : _particle(pd), _momentum(p), _dir(dir) 
  {
    if(_dir==outgoing) _momentum *= -1.0; 
    if ( dir != outgoing ) {
      tcPDPtr anti = pd->CC();
      if ( anti ) _particle = anti;
    }
  }
  //@}


  /**
   * Access to the momentum components and mass
   */
  //@{
  /**
   * Get the x component of the momentum.
   */
  Energy px() const {return _momentum.x();}

  /**
   * Get the y component of the momentum.
   */
  Energy py() const {return _momentum.y();}

  /**
   * Get the z component of the momentum.
   */
  Energy pz() const {return _momentum.z();}

  /**
   * Get the energy.
   */
  Energy e()  const {return _momentum.e();}

  /**
   * Get the mass.
   */
  Energy mass() const {return _momentum.mass();}

  /**
   * Get off-shell mass squared.
   */
  Energy2 m2() const {return _momentum.m2();}

  /**
   *  Access to the 5-momentum
   */
  const Lorentz5Momentum & momentum() const {return _momentum;}
  //@}

  /**
   *  Access to the particle properties
   */
  //@{
  /** 
   * Get the particle id.
   */
  long id() const {return _particle->id();}

  /** 
   * Get 2s+1 for the particle.
   */
  PDT::Spin iSpin() const {return _particle->iSpin();}

  /**
   * Get the particle pointer.
   */
  tcPDPtr particle() const {return _particle;}

  /** 
   * Get the direction of particle.
   */
  ThePEG::Helicity::Direction direction() const {return _dir;}

  /** 
   * Set the direction of the particle 
   */ 
  void direction(ThePEG::Helicity::Direction in) {_dir=in;} 
  //@}

protected:
  
  /**
   *  Perform the Lorentz transformation of the wave function
   */
  void transformMomentum(const LorentzRotation & r) {
    _momentum.transform(r);
  }

private:

  /**
   * Constant pointer to the particle info.
   */
  tcPDPtr _particle;

  /**
   * Lorentz 5 momentum.
   */
  Lorentz5Momentum _momentum;

  /**
   * Incoming or outgoing.
   */
  Direction _dir;
};
}
}

#endif
