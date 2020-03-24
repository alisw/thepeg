// -*- C++ -*-
//
// VVVTVertex.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2019 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_VVVTVertex_H
#define ThePEG_VVVTVertex_H
//
// This is the declaration of the VVVTVertex class.
//
#include "ThePEG/Helicity/Vertex/AbstractVVVTVertex.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "VVVTVertex.fh"

namespace ThePEG {
namespace Helicity {
  
/** \ingroup Helicity
 *
 *  The VVTVertex class is the implementation of the 
 *  vector-vector-vector-tensor vertex. 
 *  It inherits from the AbstractVVVTVertex class for the storage of the particles
 *  interacting at the vertex and implements the helicity amplitude calculations.
 *
 *  All implementations of this vertex should inherit from it and implement the
 *  virtual setCoupling member.
 *
 *  The vertex has the form
 *  \f[
 *  g\frac\kappa2f^{abc}\left[
 *  C_{\mu\nu,\rho\sigma}(k_1-k_2)_\lambda+C_{\mu\nu,\rho\lambda}(k_3-k_1)_\sigma
 * +C_{\mu\nu,\sigma\lambda}(k_2-k_3)_\rho+F_{\mu\nu,\rho\sigma\lambda}
 *   \right]\epsilon_1^\rho\epsilon^\sigma_2\epsilon^\lambda_3
 *   \epsilon^{\mu\nu}_4
 *  \f]
 *  where
 *  -\f$C_{\mu\nu,\rho\sigma}=g_{\mu\rho}g_{\nu\sigma}+g_{\mu\sigma}g_{\nu\rho}
 *         -g_{\mu\nu}g_{\rho\sigma}\f$
 *  -\f$F_{\mu\nu,\rho\sigma\lambda} = 
 *      g_{\mu\rho}g_{\sigma\lambda}(k_2-k_3)_\nu
 *     +g_{\mu\sigma}g_{\rho\lambda}(k_3-k_1)_\nu
 *     +g_{\mu\lambda}g_{\rho\sigma}(k_1-k_2)_\nu+(\mu\leftrightarrow\nu)
 *   \f$
 *
 *  @see AbstractVVVTVertex
 */
class VVVTVertex: public AbstractVVVTVertex {
    
public:

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
   * @param vec1  The wavefunction for the first  vector.
   * @param vec2  The wavefunction for the second vector.
   * @param vec3  The wavefunction for the third  vector.
   * @param ten4  The wavefunction for the tensor.
   */
  Complex evaluate(Energy2 q2,const VectorWaveFunction & vec1,
		   const VectorWaveFunction & vec2,
		   const VectorWaveFunction & vec3, const TensorWaveFunction & ten4);

  /**
   * Evaluate the off-shell tensor coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell tensor.
   * @param out The ParticleData pointer for the off-shell tensor.
   * @param vec1  The wavefunction for the first  vector.
   * @param vec2  The wavefunction for the second vector.
   * @param vec3  The wavefunction for the third  vector.
   * @param mass The mass of the off-shell particle if not taken from the ParticleData
   * object
   * @param width The width of the off-shell particle if not taken from the ParticleData
   * object
   */
  TensorWaveFunction evaluate(Energy2 q2,int iopt, tcPDPtr out,
			      const VectorWaveFunction & vec1,
			      const VectorWaveFunction & vec2,
			      const VectorWaveFunction & vec3,
			      complex<Energy> mass=-GeV, complex<Energy> width=-GeV);

  /**
   * Evaluate the off-shell vector coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell vector.
   * @param out The ParticleData pointer for the off-shell vector.
   * @param vec1  The wavefunction for the first  vector.
   * @param vec2  The wavefunction for the second vector.
   * @param ten4  The wavefunction for the tensor.
   * @param mass The mass of the off-shell particle if not taken from the ParticleData
   * object
   * @param width The width of the off-shell particle if not taken from the ParticleData
   * object
   */
  VectorWaveFunction evaluate(Energy2 q2,int iopt, tcPDPtr out,
			      const VectorWaveFunction & vec1,
			      const VectorWaveFunction & vec2,
			      const TensorWaveFunction & ten4,
			      complex<Energy> mass=-GeV, complex<Energy> width=-GeV);
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
  
private:
  
  /**
   * Private and non-existent assignment operator.
   */
  VVVTVertex & operator=(const VVVTVertex &) = delete;
  
};

}

}
#endif /* ThePEG_VVVTVertex_H */
