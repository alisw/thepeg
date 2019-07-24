// -*- C++ -*-
//
// RSSpinorWaveFunction.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2017 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_RSSpinorWaveFunction_H
#define ThePEG_RSSpinorWaveFunction_H
// This is the declaration of the RSSpinorWaveFunction class.

#include "WaveFunctionBase.h"
#include <ThePEG/Helicity/LorentzRSSpinor.h>
#include <ThePEG/Helicity/RSFermionSpinInfo.h>
#include <ThePEG/EventRecord/Particle.h>
#include <ThePEG/EventRecord/RhoDMatrix.h>

namespace ThePEG {

namespace Helicity {

/** \ingroup Helicity
 *
 *  The RSSpinorWaveFunction class is designed to store the wavefunction
 *  of a spin-3/2 particle in a form suitable for use in helicity amplitude
 *  calculations of the matrix element using a similar philosophy to the 
 *  FORTRAN HELAS code.
 *
 *  In addition to storing the spinor using the LorentzRSSpinor class
 *  it inherits from the <code>WaveFunctionBase</code> class to provide storage of
 *  the momentum and ParticleData for the fermion.
 *
 *  This class also contains the code which does the actually calculation of the
 *  spinor for an external particle.
 *
 *  When calculating the wavefunction the direction of the particle is used,
 *
 *  \e i.e. 
 *  - incoming calculates a \f$u\f$ spinor.
 *  - outgoing calculates a \f$v\f$ spinor.
 *
 *  The spinors are calculated using a Clebsch-Gordon decomposition in the rest-frame
 *  for a massive particle and boosted to the lab-frame. For massless particles the
 *  calculation is performed in the lab-frame (N.B. there are only two helicities
 *  \f$\pm\frac32\f$ in this case.)
 *
 *  N.B. In our convention 0 is the \f$-\frac32\f$ helicity state,
 *        1 is the \f$-\frac12\f$ helicity state,
 *        2 is the \f$+\frac12\f$ helicity state
 *        3 is the \f$+\frac32\f$ helicity state and 
 *
 * @see WaveFunctionBase
 * @see LorentzRSSpinor
 * @see HelicityDefinitions
 * 
 */
class RSSpinorWaveFunction: public WaveFunctionBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Constructor, set the momentum and the components of the spinor.
   * @param p The momentum.
   * @param part The ParticleData pointer.
   * @param xs1 The first  spinor component of the \f$x\f$ vector.
   * @param xs2 The second spinor component of the \f$x\f$ vector.
   * @param xs3 The third  spinor component of the \f$x\f$ vector.
   * @param xs4 The fourth spinor component of the \f$x\f$ vector.
   * @param ys1 The first  spinor component of the \f$y\f$ vector.
   * @param ys2 The second spinor component of the \f$y\f$ vector.
   * @param ys3 The third  spinor component of the \f$y\f$ vector.
   * @param ys4 The fourth spinor component of the \f$y\f$ vector.
   * @param zs1 The first  spinor component of the \f$z\f$ vector.
   * @param zs2 The second spinor component of the \f$z\f$ vector.
   * @param zs3 The third  spinor component of the \f$z\f$ vector.
   * @param zs4 The fourth spinor component of the \f$z\f$ vector.
   * @param ts1 The first  spinor component of the \f$t\f$ vector.
   * @param ts2 The second spinor component of the \f$t\f$ vector.
   * @param ts3 The third  spinor component of the \f$t\f$ vector.
   * @param ts4 The fourth spinor component of the \f$t\f$ vector.
   */
  RSSpinorWaveFunction(const Lorentz5Momentum & p,tcPDPtr part,
		       complex<double> xs1, complex<double> xs2,
		       complex<double> xs3, complex<double> xs4,
		       complex<double> ys1, complex<double> ys2,
		       complex<double> ys3, complex<double> ys4,
		       complex<double> zs1, complex<double> zs2,
		       complex<double> zs3, complex<double> zs4,
		       complex<double> ts1, complex<double> ts2,
		       complex<double> ts3, complex<double> ts4)
    : WaveFunctionBase(p,part), _wf(xs1,xs2,xs3,xs4,
				    ys1,ys2,ys3,ys4,
				    zs1,zs2,zs3,zs4,
				    ts1,ts2,ts3,ts4)
  {
    assert(iSpin()==4);
  }

  /**
   * Constructor, set the momentum and the wavefunction.
   * @param p The momentum.
   * @param part The ParticleData pointer.
   * @param wave The wavefunction.
   * @param dir The direction of the particle
   */
  RSSpinorWaveFunction(const Lorentz5Momentum & p, tcPDPtr part,
		       const LorentzRSSpinor<double> & wave,
		       Direction dir=intermediate) 
    : WaveFunctionBase(p,part,dir), _wf(wave)
  {
    assert(iSpin()==4);
  }
  
  /**
   * Constructor, set the momentum and the wavefunction.
   * @param p The Particle pointer.
   * @param wave The wavefunction.
   * @param dir The direction of the particle
   */
  RSSpinorWaveFunction(const tPPtr & p,
		       const LorentzRSSpinor<SqrtEnergy> & wave,
		       Direction dir=intermediate) 
    : WaveFunctionBase(p->momentum(),p->dataPtr(),dir), _wf(wave.Type())
  {
    assert(iSpin()==4);
    for (unsigned int i=0; i<4; ++i)
      for(unsigned int j=0; j<4; ++j)
	_wf(i,j)=Complex(wave(i,j)*UnitRemoval::InvSqrtE);
  }

  /**
   * Constructor, set the momentum, helicity, direction.
   * @param p The momentum.
   * @param part The ParticleData pointer.
   * @param ihel The helicity (0,1,2,3 as described above.)
   * @param dir The direction.
   */
  RSSpinorWaveFunction(const Lorentz5Momentum & p,tcPDPtr part,
		       unsigned int ihel, Direction dir) 
    : WaveFunctionBase(p,part,dir)
  {
    assert(iSpin()==4);
    calculateWaveFunction(ihel);
  }

  /**
   * Constructor, set the momentum, direction, zero the 
   * wavefunction.
   * @param p The momentum.
   * @param part The ParticleData pointer.
   * @param dir The direction.
   */
  RSSpinorWaveFunction(const Lorentz5Momentum & p,tcPDPtr part,Direction dir)
    : WaveFunctionBase(p,part,dir), _wf()
  {
    assert(iSpin()==4);
  }

  /**
   * Default constructor
   */
  RSSpinorWaveFunction() 
    : WaveFunctionBase(), _wf()
  {}
  //@}

  /**
   *  Access to the wavefunction and its components.
   */
  //@{
  /**
   * subscript operator for the wavefunction
   * Set components by index.
   */
  complex<double> operator ()(int i, int j) const {
    assert(  i>=0 && i<=3 && j>=0 && j<=3);
    return _wf(i,j);
  }

  /**
   * return wavefunction as LorentzRSSpinor
   */
  const LorentzRSSpinor<double> & wave() const {return _wf;}

  /**
   * Get first spinor component for the x vector
   */
  complex<double> xs1() const {return _wf.xs1();}

  /**
   * Get second spinor component for the x vector
   */
  complex<double> xs2() const {return _wf.xs2();}

  /**
   * Get third  spinor component for the x vector
   */
  complex<double> xs3() const {return _wf.xs3();}

  /**
   * Get fourth  spinor component for the x vector
   */
  complex<double> xs4() const {return _wf.xs4();}

  /**
   * Get first spinor component for the y vector
   */
  complex<double> ys1() const {return _wf.ys1();}

  /**
   * Get second spinor component for the y vector
   */
  complex<double> ys2() const {return _wf.ys2();}
  
  /**
   * Get third spinor component for the y vector
   */
  complex<double> ys3() const {return _wf.ys3();}
  
  /**
   * Get fourth spinor component for the y vector
   */
  complex<double> ys4() const {return _wf.ys4();}
  
  /**
   * Get first spinor component for the z vector
   */
  complex<double> zs1() const {return _wf.zs1();}
  
  /**
   * Get second spinor component for the z vector
   */
  complex<double> zs2() const {return _wf.zs2();}
  
  /**
   * Get third spinor component for the z vector
   */
  complex<double> zs3() const {return _wf.zs3();}
  
  /**
   * Get fourth spinor component for the z vector
   */
  complex<double> zs4() const {return _wf.zs4();}
  
  /**
   * Get first spinor component for the t vector
   */
  complex<double> ts1() const {return _wf.ts1();}
  
  /**
   * Get second spinor component for the t vector
   */
  complex<double> ts2() const {return _wf.ts2();}
  
  /**
   * Get third spinor component for the t vector
   */
  complex<double> ts3() const {return _wf.ts3();}
  
  /**
   * Get fourth spinor component for the t vector
   */
  complex<double> ts4() const {return _wf.ts4();}
  //@}

  /**
   * reset functions
   */
  //@{
  /**
   * Reset the helicity (calculates the new spinor).
   * @param ihel The helicity (0,1,2,3 as described above.)
   */
  void reset(unsigned int ihel) {
    calculateWaveFunction(ihel);
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
  static void calculateWaveFunctions(vector<LorentzRSSpinor<SqrtEnergy> > & waves,
				     tPPtr particle,Direction);

  /**
   *  Calculate the wavefunctions
   */
  static void calculateWaveFunctions(vector<RSSpinorWaveFunction> & waves,
				     tPPtr particle,Direction);

  /**
   *  Calculate the wavefunctions
   */
  static void calculateWaveFunctions(vector<LorentzRSSpinor<SqrtEnergy> > & waves,
				     RhoDMatrix & rho,
				     tPPtr particle,Direction);

  /**
   *  Calculate the wavefunctions
   */
  static void calculateWaveFunctions(vector<RSSpinorWaveFunction> & waves,
				     RhoDMatrix & rho,
				     tPPtr particle,Direction);

  /**
   *  Construct the SpinInfo object
   */
  static void constructSpinInfo(const vector<LorentzRSSpinor<SqrtEnergy> > & waves,
				tPPtr part,Direction dir, bool time);

  /**
   *  Construct the SpinInfo object
   */
  static void constructSpinInfo(const vector<RSSpinorWaveFunction> & waves,
				tPPtr part,Direction dir, bool time);

private:

  /**
   * Calcuate the wavefunction.
   * @param ihel The helicity (0,1,2,3 as described above.)
   */
  void calculateWaveFunction(unsigned int ihel);

private:

  /**
   * storage of the Lorentz RSSpinor
   */
  LorentzRSSpinor<double> _wf;

  /// Return wavefunction as LorentzRSSpinor<SqrtEnergy>
  LorentzRSSpinor<SqrtEnergy> dimensionedWf() const {
    LorentzRSSpinor<SqrtEnergy> temp(_wf.Type());
    for (unsigned int i=0; i<4; ++i)
      for (unsigned int j=0; j<4; ++j)
	temp(i,j) = _wf(i,j)*UnitRemoval::SqrtE;
    return temp;
  }
};

}
}

#endif /* ThePEG_RSSpinorWaveFunction_H */
