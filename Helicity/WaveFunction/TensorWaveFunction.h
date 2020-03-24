// -*- C++ -*-
//
// TensorWaveFunction.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2019 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_TensorWaveFunction_H
#define ThePEG_TensorWaveFunction_H
//
// This is the declaration of the TensorWaveFunction class.
//
#include "WaveFunctionBase.h"
#include "VectorWaveFunction.h"
#include <ThePEG/Helicity/LorentzTensor.h>
#include <ThePEG/Helicity/TensorSpinInfo.h>
#include <ThePEG/EventRecord/Particle.h>
#include <ThePEG/EventRecord/RhoDMatrix.h>

namespace ThePEG {
namespace Helicity {

/**\ingroup Helicity
 * Definition of the enumerated values of the phase to include in the 
 * calculation of the polarization tensor.
 */
enum TensorPhase {
  tensor_phase, /**< Include the phase factor.*/
  tensor_nophase, /**< No phase-factor. */
  default_tensor_phase=tensor_nophase /**< Default option.*/
};

/** \ingroup Helicity
 *  \author Peter Richardson
 *
 *  The TensorWaveFunction class is designed to store the wavefunction
 *  of a tensor in a form suitable for use in helicity amplitude 
 *  calculations of the matrix element using a similar philosophy to the 
 *  FORTRAN HELAS code.
 * 
 *  In addition to storing the tensor using the LorentzTensor class
 *  it inherits from the WaveFunctionBase class to provide storage of
 *  the momentum and ParticleData for the tensor particle.
 *
 *  This class also contains the code which does the actually 
 *  calculation of the tensor wavefunction.
 *
 *  There are two choices available for the calculation of the 
 *  wavefunction. These are set using the TensorPhase enumeration 
 *  which specifies a default choice.
 *  The first choice, tensor_phase, includes a phase factor 
 *  \f$\exp(\pm i \phi)\f$ for the \f$\pm\f$ helicity states while the second, 
 *  tensor_nophase, does not.
 *
 *  N.B. In our convention 
 *        0 is the \f$-2\f$ helicity state,
 *        1 is the \f$-1\f$ helicity state,
 *        2 is the \f$ 0\f$ helicity state,
 *        3 is the \f$+1\f$ helicity state and
 *        4 is the \f$+2\f$ helicity state.
 *
 *  @see WaveFunctionBase
 *  @see LorentzTensor
 *  @see VectorWaveFunction
 */
class TensorWaveFunction : public WaveFunctionBase {

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
  TensorWaveFunction(const Lorentz5Momentum & p,tcPDPtr part,
		     const LorentzTensor<double> & wave,
		     Direction  dir=intermediate) 
    : WaveFunctionBase(p,part,dir), _wf(wave)
  {
    assert(iSpin()==PDT::Spin2);
  }

  /**
   * Constructor, set the momentum and the components of the tensor.
   * @param p The momentum.
   * @param part The ParticleData pointer
   * @param xx The \f$xx\f$ component.
   * @param xy The \f$xy\f$ component.
   * @param xz The \f$xz\f$ component.
   * @param xt The \f$xt\f$ component.
   * @param yx The \f$yx\f$ component.
   * @param yy The \f$yy\f$ component.
   * @param yz The \f$yz\f$ component.
   * @param yt The \f$yt\f$ component.
   * @param zx The \f$zx\f$ component.
   * @param zy The \f$zy\f$ component.
   * @param zz The \f$zz\f$ component.
   * @param zt The \f$zt\f$ component.
   * @param tx The \f$tx\f$ component.
   * @param ty The \f$ty\f$ component.
   * @param tz The \f$tz\f$ component.
   * @param tt The \f$tt\f$ component.
   */
  TensorWaveFunction(const Lorentz5Momentum & p,tcPDPtr part,
		     Complex xx,Complex xy,Complex xz,Complex xt,Complex yx,
		     Complex yy,Complex yz,Complex yt,Complex zx,Complex zy,
		     Complex zz,Complex zt,Complex tx,Complex ty,Complex tz,
		     Complex tt) 
    : WaveFunctionBase(p,part), _wf(xx,xy,xz,xt,
				    yx,yy,yz,yt,
				    zx,zy,zz,zt,
				    tx,ty,tz,tt)
  {
    assert(iSpin()==PDT::Spin2);
  }

  /**
   * Constructor, set the momentum, helicity, direction and optionally the phase
   * @param p The momentum.
   * @param part The ParticleData pointer
   * @param ihel The helicity (0,1,2,3,4 as described above.)
   * @param dir The direction.
   * @param phase The phase choice.
   */
  TensorWaveFunction(const Lorentz5Momentum & p,tcPDPtr part,
		     unsigned int ihel,Direction dir,
		     TensorPhase phase=default_tensor_phase) 
    : WaveFunctionBase(p,part,dir)
  {
    assert(iSpin()==PDT::Spin2);
    calculateWaveFunction(ihel,phase);
  }

  /**
   * Constructor, set the 5-momentum and direction, zero the wavefunction.
   * @param p The momentum.
   * @param part The ParticleData pointer.
   * @param dir The direction.
   */
  TensorWaveFunction(const Lorentz5Momentum & p,
		     tcPDPtr part,Direction dir) 
    : WaveFunctionBase(p,part,dir), _wf()
  {
    assert(iSpin()==PDT::Spin2);
  }

  /** 
   * Default constructor.
   */
  TensorWaveFunction() {}

  /**
   *  Special for spin correlations \todo make static?
   */
  TensorWaveFunction(vector<TensorWaveFunction> & wave,
		     tPPtr part,Direction dir,bool time,bool massless,
		     bool=true,
		     TensorPhase phase=default_tensor_phase) {
    calculateWaveFunctions(wave,part,dir,massless,phase);
    constructSpinInfo(wave,part,dir,time,massless);
  }
  //@}

  /**
   *  Access to the wavefunction and its components.
   */
  //@{
  /**
   * Subscript operator for the wavefunction.
   */
  Complex operator ()(int i, int j) const {
    return _wf(i,j);
  }

  /**
   * Set components by index.
   */
  Complex & operator () (int i, int j) {
    return _wf(i,j);
  }

  /**
   * Return wavefunction as polarization vector.
   */
  const LorentzTensor<double> & wave() const {return _wf;}

  /**
   * Get the \f$xx\f$ component.
   */
  Complex xx() const {return _wf.xx();}

  /**
   * Get the \f$yx\f$ component.
   */
  Complex yx() const {return _wf.yx();}

  /**
   * Get the \f$zx\f$ component.
   */
  Complex zx() const {return _wf.zx();}

  /**
   * Get the \f$tx\f$ component.
   */
  Complex tx() const {return _wf.tx();}

  /**
   * Get the \f$xy\f$ component.
   */
  Complex xy() const {return _wf.xy();}

  /**
   * Get the \f$yy\f$ component.
   */
  Complex yy() const {return _wf.yy();}

  /**
   * Get the \f$zy\f$ component.
   */
  Complex zy() const {return _wf.zy();}

  /**
   * Get the \f$ty\f$ component.
   */
  Complex ty() const {return _wf.ty();}

  /**
   * Get the \f$xz\f$ component.
   */
  Complex xz() const {return _wf.xz();}

  /**
   * Get the \f$yz\f$ component.
   */
  Complex yz() const {return _wf.yz();}

  /**
   * Get the \f$zz\f$ component.
   */
  Complex zz() const {return _wf.zz();}

  /**
   * Get the \f$tz\f$ component.
   */
  Complex tz() const {return _wf.tz();}

  /**
   * Get the \f$xt\f$ component.
   */
  Complex xt() const {return _wf.xt();}

  /**
   * Get the \f$yt\f$ component.
   */
  Complex yt() const {return _wf.yt();}

  /**
   * Get the \f$zt\f$ component.
   */
  Complex zt() const {return _wf.zt();}

  /**
   * Get the \f$tt\f$ component.
   */
  Complex tt() const {return _wf.tt();}
  //@}

  /**
   * Reset functions.
   */
  //@{

  /**
   * Reset helicity (recalculate the tensor ).
   * @param ihel The new helicity (0,1,2,3,4 as described above.)
   * @param phase The phase choice.
   */
  void reset(unsigned int ihel,TensorPhase phase=default_tensor_phase) {
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
  static void calculateWaveFunctions(vector<LorentzTensor<double> > & waves,
				     tPPtr particle,Direction,bool massless,
				     TensorPhase phase=default_tensor_phase);

  /**
   *  Calculate the wavefunctions
   */
  static void calculateWaveFunctions(vector<TensorWaveFunction> & waves,
				     tPPtr particle,Direction,bool massless,
				     TensorPhase phase=default_tensor_phase);

  /**
   *  Calculate the wavefunctions
   */
  static void calculateWaveFunctions(vector<LorentzTensor<double> > & waves,
				     RhoDMatrix & rho,
				     tPPtr particle,Direction,bool massless,
				     TensorPhase phase=default_tensor_phase);

  /**
   *  Calculate the wavefunctions
   */
  static void calculateWaveFunctions(vector<TensorWaveFunction> & waves,
				     RhoDMatrix & rho,
				     tPPtr particle,Direction,bool massless,
				     TensorPhase phase=default_tensor_phase);

  /**
   *  Construct the SpinInfo object
   */
  static void constructSpinInfo(const vector<LorentzTensor<double> > & waves,
				tPPtr part,Direction dir, bool time,bool massless);

  /**
   *  Construct the SpinInfo object
   */
  static void constructSpinInfo(const vector<TensorWaveFunction> & waves,
				tPPtr part,Direction dir, bool time,bool massless);

private:

  /**
   * Calculate the wavefunction.
   * @param ihel The helicity (0,1,2,3,4 as described above.)
   * @param phase The phase choice.
   */
  void calculateWaveFunction(unsigned int ihel,
			     TensorPhase phase=default_tensor_phase);

private:

  /**
   * Storage of the wavefunction as a Lorentz Tensor.
   */
  LorentzTensor<double> _wf;

};
}
}

#endif
