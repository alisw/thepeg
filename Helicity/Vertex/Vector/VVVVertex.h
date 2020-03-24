// -*- C++ -*-
//
// VVVVertex.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2019 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_VVVVertex_H
#define ThePEG_VVVVertex_H
//
// This is the declaration of the VVVVertex class.

#include "ThePEG/Helicity/Vertex/AbstractVVVVertex.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "VVVVertex.fh"

namespace ThePEG {
namespace Helicity{
  
/** \ingroup Helicity
 *
 *  The VVVVertex class is the base class for triple vector vertices
 *  using the perturbative form. 
 *  It inherits from the AbstractVVVVertex class for the storage of the 
 *  particles allowed at the vertex.
 *
 *  Classes which implement a specific vertex should inherit from this and
 *  implement the virtual setCoupling member.
 *
 *  The form of the vertex is
 *  \f[ig\left[  (p_1-p_2)^\gamma g^{\alpha\beta }
 *              +(p_2-p_3)^\alpha g^{\beta \gamma}
 *              +(p_3-p_1)^\beta  g^{\alpha\gamma}
 *   \right]\epsilon_{1\alpha}\epsilon_{2\beta}\epsilon_{3\gamma}\f]
 *
 *  @see AbstractVVVVertex
 */
class VVVVertex: public AbstractVVVVertex {
    
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
   * Evaluate the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param vec1 The wavefunction for the first  vector.
   * @param vec2 The wavefunction for the second vector.
   * @param vec3 The wavefunction for the third  vector.
   */
  Complex evaluate(Energy2 q2, const VectorWaveFunction & vec1,
		   const VectorWaveFunction & vec2, const VectorWaveFunction & vec3);

  /**
   * Evaluate the off-shell vector coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell vector.
   * @param out The ParticleData pointer for the off-shell vector.
   * @param vec2 The wavefunction for the second vector.
   * @param vec3 The wavefunction for the third  vector.
   * @param mass The mass of the off-shell particle if not taken from the ParticleData
   * object
   * @param width The width of the off-shell particle if not taken from the ParticleData
   * object
   */
  VectorWaveFunction evaluate(Energy2 q2,int iopt, tcPDPtr out,
			      const VectorWaveFunction & vec2,
			      const VectorWaveFunction & vec3,
			      complex<Energy> mass=-GeV, complex<Energy> width=-GeV);
  //@}

  /**
   *   Set coupling methods
   */
  //@{
  /**
   * Calculate the couplings for a three point interaction.
   * This method is virtual and must be implemented in 
   * classes inheriting from this.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param part1 The ParticleData pointer for the first  particle.
   * @param part2 The ParticleData pointer for the second particle.
   * @param part3 The ParticleData pointer for the third  particle.
   */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1,
			   tcPDPtr part2,tcPDPtr part3)=0;

  /**
   * Dummy setCouplings for a four point interaction 
   * This method is virtual and must be implemented in 
   * classes inheriting from this.
   */
  virtual void setCoupling(Energy2,tcPDPtr,tcPDPtr,tcPDPtr,tcPDPtr) {
    assert(false);
  }
  //@}
  
private:
  
  /**
   * Private and non-existent assignment operator.
   */
  VVVVertex & operator=(const VVVVertex &) = delete;
  
};

}

}
#endif /* ThePEG_VVVVertex_H */
