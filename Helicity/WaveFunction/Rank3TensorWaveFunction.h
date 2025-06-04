// -*- C++ -*-
//
// Rank3TensorWaveFunction.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2019 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_Rank3TensorWaveFunction_H
#define ThePEG_Rank3TensorWaveFunction_H
//
// This is the declaration of the Rank3TensorWaveFunction class.
//
#include "WaveFunctionBase.h"
#include "VectorWaveFunction.h"
#include <ThePEG/Helicity/LorentzRank3Tensor.h>
#include <ThePEG/Helicity/Rank3TensorSpinInfo.h>
#include <ThePEG/EventRecord/Particle.h>
#include <ThePEG/EventRecord/RhoDMatrix.h>

namespace ThePEG {
namespace Helicity {

/** \ingroup Helicity
 *  \author Peter Richardson
 *
 *  The Rank3TensorWaveFunction class is designed to store the wavefunction
 *  of a tensor in a form suitable for use in helicity amplitude 
 *  calculations of the matrix element using a similar philosophy to the 
 *  FORTRAN HELAS code.
 * 
 *  In addition to storing the tensor using the LorentzRank3Tensor class
 *  it inherits from the WaveFunctionBase class to provide storage of
 *  the momentum and ParticleData for the tensor particle.
 *
 *  This class also contains the code which does the actually 
 *  calculation of the tensor wavefunction.
 *
 *  N.B. In our convention 
 *        0 is the \f$-3\f$ helicity state,
 *        1 is the \f$-2\f$ helicity state,
 *        2 is the \f$-1\f$ helicity state,
 *        3 is the \f$ 0\f$ helicity state,
 *        4 is the \f$+1\f$ helicity state and
 *        5 is the \f$+2\f$ helicity state.
 *        6 is the \f$+3\f$ helicity state.
 *
 *  @see WaveFunctionBase
 *  @see LorentzRank3Tensor
 *  @see VectorWaveFunction
 */
class Rank3TensorWaveFunction : public WaveFunctionBase {

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
  Rank3TensorWaveFunction(const Lorentz5Momentum & p,tcPDPtr part,
			  const LorentzRank3Tensor<double> & wave,
			  Direction  dir=intermediate) 
    : WaveFunctionBase(p,part,dir), _wf(wave) {
    assert(iSpin()==PDT::Spin3);
  }

  /**
   * Constructor, set the momentum, helicity, direction and optionally the phase
   * @param p The momentum.
   * @param part The ParticleData pointer
   * @param ihel The helicity (0,1,2,3,4 as described above.)
   * @param dir The direction.
   * @param phase The phase choice.
   */
  Rank3TensorWaveFunction(const Lorentz5Momentum & p,tcPDPtr part,
			  unsigned int ihel,Direction dir) 
    : WaveFunctionBase(p,part,dir) {
    assert(iSpin()==PDT::Spin3);
    calculateWaveFunction(ihel);
  }

  /**
   * Constructor, set the 5-momentum and direction, zero the wavefunction.
   * @param p The momentum.
   * @param part The ParticleData pointer.
   * @param dir The direction.
   */
  Rank3TensorWaveFunction(const Lorentz5Momentum & p,
			  tcPDPtr part,Direction dir) 
    : WaveFunctionBase(p,part,dir), _wf() {
    assert(iSpin()==PDT::Spin3);
  }

  /** 
   * Default constructor.
   */
  Rank3TensorWaveFunction() {}

  /**
   *  Special for spin correlations \todo make static?
   */
  Rank3TensorWaveFunction(vector<Rank3TensorWaveFunction> & wave,
			  tPPtr part,Direction dir,bool time,bool massless,
			  bool=true) {
    calculateWaveFunctions(wave,part,dir,massless);
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
  Complex operator ()(int i, int j, int k) const {
    return _wf(i,j,k);
  }

  /**
   * Set components by index.
   */
  Complex & operator () (int i, int j, int k) {
    return _wf(i,j,k);
  }

  /**
   * Return wavefunction as polarization rank-3 tensor .
   */
  const LorentzRank3Tensor<double> & wave() const {return _wf;}
  //@}

  /**
   * Reset functions.
   */
  //@{

  /**
   * Reset helicity (recalculate the tensor ).
   * @param ihel The new helicity (0,1,2,3,4 as described above.)
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
  static void calculateWaveFunctions(vector<LorentzRank3Tensor<double> > & waves,
				     tPPtr particle,Direction,bool massless);

  /**
   *  Calculate the wavefunctions
   */
  static void calculateWaveFunctions(vector<Rank3TensorWaveFunction> & waves,
				     tPPtr particle,Direction,bool massless);

  /**
   *  Calculate the wavefunctions
   */
  static void calculateWaveFunctions(vector<LorentzRank3Tensor<double> > & waves,
				     RhoDMatrix & rho,
				     tPPtr particle,Direction,bool massless);

  /**
   *  Calculate the wavefunctions
   */
  static void calculateWaveFunctions(vector<Rank3TensorWaveFunction> & waves,
				     RhoDMatrix & rho,
				     tPPtr particle,Direction,bool massless);

  /**
   *  Construct the SpinInfo object
   */
  static void constructSpinInfo(const vector<LorentzRank3Tensor<double> > & waves,
				tPPtr part,Direction dir, bool time,bool massless);

  /**
   *  Construct the SpinInfo object
   */
  static void constructSpinInfo(const vector<Rank3TensorWaveFunction> & waves,
				tPPtr part,Direction dir, bool time,bool massless);

private:

  /**
   * Calculate the wavefunction.
   * @param ihel The helicity (0,1,2,3,4 as described above.)
   */
  void calculateWaveFunction(unsigned int ihel);

private:

  /**
   * Storage of the wavefunction as a Lorentz Rank3Tensor.
   */
  LorentzRank3Tensor<double> _wf;

};
}
}

#endif
