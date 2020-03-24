// -*- C++ -*-
//
// FFVTVertex.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2019 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_FFVTVertex_H
#define ThePEG_FFVTVertex_H
//
// This is the declaration of the FFVTVertex class.
//
#include "ThePEG/Helicity/Vertex/AbstractFFVTVertex.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "FFVTVertex.fh"

namespace ThePEG {
namespace Helicity {

/** \ingroup Helicity
 *  
 *  The FFVTVertex class is the implementation of the 
 *  fermion-fermion--vector-tensor vertex. 
 *  It inherits from the AbstractFFVTVertex class for the storage of the particles
 *  interacting at the vertex and implements the helicity amplitude calculations.
 *
 *  All implementations of this vertex should inherit from it and implement the
 *  virtual setCoupling member.
 *
 *  The form of the vertex is
 * \f[\frac{ig\kappa}4t^a_{nm}\bar{f_2}(C_{\mu\nu,\rho\sigma}-g_{\mu\nu}g_{\rho\sigma})
 * \gamma^\sigma f_1\epsilon_{3\rho}\epsilon_4^{\mu\nu}\f]
 *  where
 *  -\f$C_{\mu\nu,\rho\sigma}=g_{\mu\rho}g_{\nu\sigma}+g_{\mu\sigma}g_{\nu\rho}
 *         -g_{\mu\nu}g_{\rho\sigma}\f$.
 *
 *  @see AbstractFFVTVertex
 */
class FFVTVertex: public AbstractFFVTVertex {
      
public:

  /**
   *  Default constructor
   */
  FFVTVertex() : left_(1.), right_(1.) {}

  /**
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();
  
public:
  
  /**
   * Members to calculate the helicity amplitude expressions for vertices
   * and off-shell particles.
   */
  //@{
  /**
   * Evalulate the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param sp1   The wavefunction for the ferimon.
   * @param sbar2 The wavefunction for the antifermion.
   * @param vec3  The wavefunction for the vector.
   * @param ten4  The wavefunction for the tensor.
   */
  Complex evaluate(Energy2 q2,const SpinorWaveFunction & sp1,
		   const SpinorBarWaveFunction & sbar2,
		   const VectorWaveFunction & vec3, const TensorWaveFunction & ten4);
  //@}

  /**
   *   Set coupling methods
   */
  //@{
  /**
   * Dummy for a three point interaction.
   */
  virtual void setCoupling(Energy2,tcPDPtr,tcPDPtr,tcPDPtr) {
    assert(false);
  }

  /**
   * Calculate the couplings for a four point interaction.
   * This method is virtual and must be implemented in 
   * classes inheriting from this.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param part1 The ParticleData pointer for the first  particle.
   * @param part2 The ParticleData pointer for the second particle.
   * @param part3 The ParticleData pointer for the third  particle.
   * @param part4 The ParticleData pointer for the fourth particle.
   */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1,tcPDPtr part2,tcPDPtr part3,
			   tcPDPtr part4)=0;
  //@}

  /**
   *  Left and right couplings
   */
  //@{
  /**
   *  Get left
   */
  Complex left() const {return left_;}

  /**
   *  Set left
   */
  void left(Complex in) {left_ = in;}

  /**
   *  Get right
   */
  Complex right() const {return right_;}

  /**
   *  Set right
   */
  void right(Complex in) {right_ = in;}
  //@}
    
private:
  
  /**
   * Private and non-existent assignment operator.
   */
  FFVTVertex & operator=(const FFVTVertex &) = delete;
  
private:

  /**
   *  Left coupling
   */
  Complex left_;

  /**
   *  Right coupling
   */
  Complex right_;

};
}

}
#endif /* ThePEG_FFVTVertex_H */
