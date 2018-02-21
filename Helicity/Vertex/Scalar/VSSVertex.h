// -*- C++ -*-
//
// VSSVertex.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2017 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_VSSVertex_H
#define ThePEG_VSSVertex_H
//
// This is the declaration of the VSSVertex class.
//
#include "ThePEG/Helicity/Vertex/AbstractVSSVertex.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "VSSVertex.fh"

namespace ThePEG {
namespace Helicity {

/** \ingroup Helicity
 *
 *  The VSSVertex class is the implementation of the vector-scalar-scalar
 *  vertex. It inherits from the AbstractVSSVertex class for storage of the particles 
 *  and implements the helicity calculations.
 *
 *  All such vertices should inherit from this class and implement the virtual
 *  setCoupling member
 *
 *  The form of the vertex is
 * \f[-ic\left(p_2-p_3\right)\cdot\epsilon_1\phi_2\phi_3\f]
 *
 *  @see AbstractVSSVertex
 */
class VSSVertex: public AbstractVSSVertex {
    
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
   * @param vec1 The wavefunction for the vector.
   * @param sca2 The wavefunction for the first  scalar.
   * @param sca3 The wavefunction for the second scalar.
   */
  Complex evaluate(Energy2 q2,const VectorWaveFunction & vec1,
		   const ScalarWaveFunction & sca2, const ScalarWaveFunction & sca3);

  /**
   * Evaluate the off-shell vector coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell vector.
   * @param out The ParticleData pointer for the off-shell vector.
   * @param sca2 The wavefunction for the first  scalar.
   * @param sca3 The wavefunction for the second scalar.
   * @param mass The mass of the off-shell particle if not taken from the ParticleData
   * object
   * @param width The width of the off-shell particle if not taken from the ParticleData
   * object
   */
  VectorWaveFunction evaluate(Energy2 q2,int iopt,tcPDPtr out,
			      const ScalarWaveFunction & sca2,
			      const ScalarWaveFunction & sca3,
			      complex<Energy> mass=-GeV, complex<Energy> width=-GeV);

  /**
   * Evaluate the off-shell scalar coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell scalar.
   * @param out The ParticleData pointer for the off-shell scalar.
   * @param vec1 The wavefunction for the vector.
   * @param sca2 The wavefunction for the scalar.
   * @param mass The mass of the off-shell particle if not taken from the ParticleData
   * object
   * @param width The width of the off-shell particle if not taken from the ParticleData
   * object
   */
  ScalarWaveFunction evaluate(Energy2 q2,int iopt, tcPDPtr out,
			      const VectorWaveFunction & vec1,
			      const ScalarWaveFunction & sca2,
			      complex<Energy> mass=-GeV, complex<Energy> width=-GeV);

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
   * Describe an abstract class with persistent data.
   */
  static AbstractNoPIOClassDescription<VSSVertex> initVSSVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  VSSVertex & operator=(const VSSVertex &);
  
};

}

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of VSSVertex.
 */
template <>
struct BaseClassTrait<ThePEG::Helicity::VSSVertex,1> {
  /** Typedef of the base class of VSSVertex. */
  typedef ThePEG::Helicity::AbstractVSSVertex NthBase;
};
  
/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<ThePEG::Helicity::VSSVertex>
  : public ClassTraitsBase<ThePEG::Helicity::VSSVertex> {
  
  /**
   * Return the class name.
   */
  static string className() { return "ThePEG::VSSVertex"; }
};

/** @endcond */
  
}

#endif /* ThePEG_VSSVertex_H */
