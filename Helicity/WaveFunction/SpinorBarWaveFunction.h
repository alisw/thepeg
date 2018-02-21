// -*- C++ -*-
//
// SpinorBarWaveFunction.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2017 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_SpinorBarWaveFunction_H
#define ThePEG_SpinorBarWaveFunction_H
//
// This is the declaration of the SpinorBarWaveFunction class.

#include "WaveFunctionBase.h"
#include <ThePEG/Helicity/LorentzSpinorBar.h>
#include <ThePEG/Helicity/FermionSpinInfo.h>
#include <ThePEG/EventRecord/Particle.h>
#include <ThePEG/EventRecord/RhoDMatrix.h>

namespace ThePEG {

namespace Helicity {

/**
 *  Forward declaration of the SpinorWaveFunction class
 */
class SpinorWaveFunction;

/** \ingroup Helicity
 *  \author Peter Richardson
 *
 *  The SpinorBarWaveFunction class is designed to store the wavefunction
 *  of a barred spinor in a form suitable for use in helicity amplitude 
 *  calculations of the matrix element using a similar philosophy to the 
 *  FORTRAN HELAS code.
 *
 *  In addition to storing the spinor using the LorentzSpinorBar class
 *  it inherits from the WaveFunctionBase class to provide storage of
 *  the momentum and ParticleData for the fermion.
 *
 *  This class also contains the code which does the actually calculation 
 *  of the barred spinor for an external particle.
 *
 *  When calculating the wavefunction the direction of the particle is used,
 *
 *  \e i.e. 
 *  - incoming calculates a \f$\bar{v}\f$ spinor.
 *  - outgoing calculates a \f$\bar{u}\f$ spinor.
 *
 *  N.B. In our convention 0 is the \f$-\frac12\f$ helicity state and 
 *        1 is the \f$+\frac12\f$ helicity state
 *
 *  @see WaveFunctionBase
 *  @see LorentzSpinorBar
 *  @see HelicityDefinitions
 */
class SpinorBarWaveFunction : public WaveFunctionBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Constructor, set the momentum and the components of the spinor.
   * @param p The momentum.
   * @param part The ParticleData pointer.
   * @param s1 The first component
   * @param s2 The second component
   * @param s3 The third component
   * @param s4 The fourth component
   */
  SpinorBarWaveFunction(const Lorentz5Momentum & p,tcPDPtr part,
			complex<double> s1,complex<double> s2,
			complex<double> s3,complex<double> s4)
    : WaveFunctionBase(p,part), _wf(s1,s2,s3,s4)
  {
    assert(iSpin()==2);
  }
  
  
  /**
   * Constructor, set the momentum and the wavefunction.
   * @param p The momentum.
   * @param part The ParticleData pointer.
   * @param wave The wavefunction.
   * @param dir The direction of the particle
   */
  SpinorBarWaveFunction(const Lorentz5Momentum & p,tcPDPtr part,
			const LorentzSpinorBar<double> & wave,
			Direction dir=intermediate) 
    : WaveFunctionBase(p,part,dir), _wf(wave)
  {
    assert(iSpin()==2);
  }
  
  SpinorBarWaveFunction(const tPPtr & p, 
			const LorentzSpinorBar<SqrtEnergy> & wave,
			Direction dir=intermediate) 
    : WaveFunctionBase(p->momentum(),p->dataPtr(),dir), 
      _wf(wave.Type())
  {
    assert(iSpin()==2);
    for (unsigned int i=0; i<4; ++i)
      _wf[i]=wave[i]*UnitRemoval::InvSqrtE;
  }
 
  /**
   * Constructor, set the momentum, helicity, direction.
   * @param p The momentum.
   * @param part The ParticleData pointer.
   * @param ihel The helicity (0,1 as described above.)
   * @param dir The direction.
   */
  SpinorBarWaveFunction(const Lorentz5Momentum & p,tcPDPtr part,
			unsigned int ihel,Direction dir)
    : WaveFunctionBase(p,part,dir)
  {
    assert(iSpin()==2);
    calculateWaveFunction(ihel);
  }

  /**
   * Constructor, set the momentum, direction, zero the 
   * wavefunction.
   * @param p The momentum.
   * @param part The ParticleData pointer.
   * @param dir The direction.
   */
  SpinorBarWaveFunction(const Lorentz5Momentum & p,tcPDPtr part,
			Direction dir) 
    : WaveFunctionBase(p,part,dir), _wf()
  {
    assert(iSpin()==2);
  }

  /**
   * Default constructor.
   */
  SpinorBarWaveFunction() 
    : WaveFunctionBase(), _wf()
  {}

  /**
   *  Special for spin correlations
   */
  SpinorBarWaveFunction(vector<SpinorBarWaveFunction> & wave,
			tPPtr part,Direction dir,bool time,bool=true) {
    calculateWaveFunctions(wave,part,dir);
    constructSpinInfo(wave,part,dir,time);
  }
  //@}

  /**
   *  Access to the wavefunction and its components.
   */
  //@{
  /**
   * Subscript operator for the wavefunction.
   */
  complex<double> operator ()(int i) const {
    assert(i>=0 &&i<=3);
    return _wf(i);
  }
    
  /**
   * Return wavefunction as LorentzSpinorBar<double>.
   */
  const LorentzSpinorBar<double> & wave() const {return _wf;}
  
  /// Return wavefunction as LorentzSpinorBar<SqrtEnergy>
  LorentzSpinorBar<SqrtEnergy> dimensionedWave() const {
    return dimensionedWf();
  }
  
  /**
   * Get the first spin component component.
   */
  complex<double> s1() const {return _wf.s1();}

  /**
   * Get the second spin component component.
   */
  complex<double> s2() const {return _wf.s2();}

  /**
   * Get the third spin component component.
   */
  complex<double> s3() const {return _wf.s3();}

  /**
   * Get the fourth spin component component.
   */
  complex<double> s4() const {return _wf.s4();}

  /**
   * Take the conjugate of the spinor \f$u_c=C\bar{u}^T\f$. This operation
   * transforms u-spinors to v-spinors and vice-versa and is required when
   * dealing with majorana particles.
   */
  void conjugate();

  /**
   * Return the barred spinor
   */
  SpinorWaveFunction bar();

  /**
   * Reset functions.
   */
  //@{
  /**
   * Reset the helicity (calculates the new spinor).
   * @param ihel The helicity (0,1 as described above.)
   */
  void reset(unsigned int ihel) {
    calculateWaveFunction(ihel);
  }
  //@}

private:

  /**
   * Calcuate the wavefunction.
   * @param ihel The helicity (0,1 as described above.)
   */
  void calculateWaveFunction(unsigned int ihel);


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
  static void calculateWaveFunctions(vector<LorentzSpinorBar<SqrtEnergy> > & waves,
				     tPPtr particle,Direction);

  /**
   *  Calculate the wavefunctions
   */
  static void calculateWaveFunctions(vector<SpinorBarWaveFunction> & waves,
				     tPPtr particle,Direction);

  /**
   *  Calculate the wavefunctions
   */
  static void calculateWaveFunctions(vector<LorentzSpinorBar<SqrtEnergy> > & waves,
				     RhoDMatrix & rho,
				     tPPtr particle,Direction);

  /**
   *  Calculate the wavefunctions
   */
  static void calculateWaveFunctions(vector<SpinorBarWaveFunction> & waves,
				     RhoDMatrix & rho,
				     tPPtr particle,Direction);

  /**
   *  Construct the SpinInfo object
   */
  static void constructSpinInfo(const vector<LorentzSpinorBar<SqrtEnergy> > & waves,
				tPPtr part,Direction dir, bool time);

  /**
   *  Construct the SpinInfo object
   */
  static void constructSpinInfo(const vector<SpinorBarWaveFunction> & waves,
				tPPtr part,Direction dir, bool time);

private:
  
  /**
   * Storage of the Lorentz SpinorBar wavefunction.
   */
  LorentzSpinorBar<double> _wf;

  /// Return wavefunction as LorentzSpinorBar<SqrtEnergy>
  LorentzSpinorBar<SqrtEnergy> dimensionedWf() const {
    LorentzSpinorBar<SqrtEnergy> temp(_wf.Type());
    for (unsigned int i=0; i<4; ++i)
	temp(i) = _wf(i)*UnitRemoval::SqrtE;
    return temp;
  }
};
}
}

#endif




