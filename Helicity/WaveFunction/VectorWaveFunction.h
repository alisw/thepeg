// -*- C++ -*-
//
// VectorWaveFunction.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2019 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_VectorWaveFunction_H
#define ThePEG_VectorWaveFunction_H
//
// This is the declaration of the VectorWaveFunction class.
//
#include "WaveFunctionBase.h"
#include <ThePEG/Helicity/LorentzPolarizationVector.h>
#include <ThePEG/Helicity/VectorSpinInfo.h>
#include <ThePEG/EventRecord/RhoDMatrix.h>
#include <ThePEG/EventRecord/Particle.h>

namespace ThePEG {
namespace Helicity {

/** \ingroup Helicity
 *
 *  \author Peter Richardson
 *
 *  The VectorWaveFunction class is designed to store the wavefunction
 *  of a vector in a form suitable for use in helicity amplitude calculations 
 *  of the matrix element using a similar philosophy to the FORTRAN HELAS code.
 *
 *  In addition to storing the vector using the LorentzPolarizationVector class
 *  it inherits from the WaveFunctionBase class to provide storage of the 
 *  momentum and ParticleData for the vector boson.
 *
 *  This class also contains the code which does the actually calculation of the
 *  vector wavefunction.
 *
 *  There are two choices available for the calculation of the wavefunction.
 *  These are set using the VectorPhase enumeration which specifies a default choice.
 *  The first choice, vector_phase, includes a phase factor \f$\exp(\pm i \phi)\f$
 *  for the \f$\pm\f$ helicity states while the second, vector_nophase, does not.
 *
 *  N.B. In our convention 0 is the \f$-1\f$ helicity state and 
 *        1 is the \f$0\f$ helicity state
 *        2 is the \f$+1\f$ helicity state
 *
 *  @see WaveFunctionBase
 *  @see LorentzPolarizationVector
 */
class VectorWaveFunction : public WaveFunctionBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Constructor, set the momentum and Wavefunction, the direction can also
   * be specified. 
   * @param p The momentum.
   * @param part The ParticleData pointer
   * @param wave The wavefunction, \e i.e. the polarization vector.
   * @param dir The direction of the particle.
   */
  VectorWaveFunction(const Lorentz5Momentum & p,tcPDPtr part,
		     const LorentzPolarizationVector & wave,
		     Direction  dir=intermediate) 
    : WaveFunctionBase(p,part,dir), _wf(wave)
  {
    assert(iSpin()==3);
  }

  /**
   * Constructor, set the momentum and components of the wavefunction.
   * @param p The momentum.
   * @param part The ParticleData pointer
   * @param x The x component of the polarization vector
   * @param y The y component of the polarization vector
   * @param z The z component of the polarization vector
   * @param t The t component of the polarization vector
   */
  VectorWaveFunction(const Lorentz5Momentum & p,tcPDPtr part,const Complex & x,
		     const Complex & y,const Complex & z, const Complex & t) 
    : WaveFunctionBase(p,part), _wf(x,y,z,t)
  {
    assert(iSpin()==3);
  }
  
  /**
   * Constructor, set the momentum, helicity and direction, optionally the choice
   * of the phase.
   * @param p The momentum.
   * @param part The ParticleData pointer.
   * @param ihel The helicity (0,1,2 as described above.)
   * @param dir The direction.
   * @param phase The phase choice.
   */
  VectorWaveFunction(const Lorentz5Momentum & p,tcPDPtr part,
		     unsigned int ihel,Direction dir,
		     VectorPhase phase=default_vector_phase) 
    : WaveFunctionBase(p,part,dir)
  {
    assert(iSpin()==3);
    calculateWaveFunction(ihel,phase);
  }
  
  /**
   * Constructor, set the 5-momentum and direction, zero the wavefunction.
   * @param p The 5-momentum.
   * @param part The ParticleData pointer.
   * @param dir The direction.
   */
  VectorWaveFunction(const Lorentz5Momentum &p,
		     tcPDPtr part,Direction dir)  
    : WaveFunctionBase(p,part,dir), _wf()
  {
    assert(iSpin()==3);
  }
  
  /**
   * Default constructor.
   */
  VectorWaveFunction() {}

  /**
   *  Special for spin correlations \todo make static?
   */
  VectorWaveFunction(vector<VectorWaveFunction> & wave,
		     tPPtr part,Direction dir,bool time,bool massless,
		     bool=true,
		     VectorPhase phase=default_vector_phase) {
    calculateWaveFunctions(wave,part,dir,massless,phase);
    constructSpinInfo(wave,part,dir,time,massless);
  }
  //@}

  /**
   *  Access to the wavefunction and its components.
   */
  //@{
  /**
   * Return wavefunction as polarization vector. 
   */
  const LorentzPolarizationVector & wave() const { return _wf;}
  
  /**
   * Get x component.
   */
  Complex x() const {return _wf.x();}
  
  /**
   * Get y component.
   */
  Complex y() const {return _wf.y();}
  
  /**
   * Get z component.
   */
  Complex z() const {return _wf.z();}
  
  /**
   * Get t component.
   */
  Complex t() const {return _wf.t();}

  /**
   * Reset functions.
   */
  //@{
  /**
   * Reset the helicity (recalculation the polarization vector).
   * @param ihel The new helicity (0,1,2 as described above.)
   * @param phase The phase choice.
   */
  void reset(unsigned int ihel,VectorPhase phase=default_vector_phase) {
    calculateWaveFunction(ihel,phase);
  }
  //@}

public:

  /**
   *  Perform the Lorentz transformation of the wave function
   */
  void transform(const LorentzRotation & r) {
    _wf.transform(r);
    transformMomentum(r);
  }

public:

  /**
   *  Calculate the wavefunctions
   */
  static void calculateWaveFunctions(vector<LorentzPolarizationVector> & waves,
				     tPPtr particle,Direction,bool massless,
				     VectorPhase phase=default_vector_phase);

  /**
   *  Calculate the wavefunctions
   */
  static void calculateWaveFunctions(vector<VectorWaveFunction> & waves,
				     tPPtr particle,Direction,bool massless,
				     VectorPhase phase=default_vector_phase);

  /**
   *  Calculate the wavefunctions
   */
  static void calculateWaveFunctions(vector<VectorWaveFunction> & waves,
				     const Lorentz5Momentum & momentum,
				     tcPDPtr parton, Direction,bool massless,
				     VectorPhase phase=default_vector_phase);

  /**
   *  Calculate the wavefunctions
   */
  static void calculateWaveFunctions(vector<LorentzPolarizationVector> & waves,
				     RhoDMatrix & rho,
				     tPPtr particle,Direction,bool massless,
				     VectorPhase phase=default_vector_phase);

  /**
   *  Calculate the wavefunctions
   */
  static void calculateWaveFunctions(vector<VectorWaveFunction> & waves,
				     RhoDMatrix & rho,
				     tPPtr particle,Direction,bool massless,
				     VectorPhase phase=default_vector_phase);

  /**
   *  Construct the SpinInfo object
   */
  static void constructSpinInfo(const vector<LorentzPolarizationVector> & waves,
				tPPtr part,Direction dir, bool time,bool massless);

  /**
   *  Construct the SpinInfo object
   */
  static void constructSpinInfo(const vector<VectorWaveFunction> & waves,
				tPPtr part,Direction dir, bool time,bool massless);

private:
  
  /**
   * Calculate the wavefunction
   * @param ihel The helicity  (0,1,2 as described above.)
   * @param phase The phase choice.
   */
  void calculateWaveFunction(unsigned int ihel,
			     VectorPhase phase=default_vector_phase);

private:
  
  /**
   * Storage of the wavefunction as a Lorentz Vector.
   */
  LorentzPolarizationVector _wf;
  
};

}
}

#endif
