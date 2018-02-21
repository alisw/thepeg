// -*- C++ -*-
//
// VVSVertex.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2017 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_VVSVertex_H
#define ThePEG_VVSVertex_H
//
// This is the declaration of the VVSVertex class.

#include "ThePEG/Helicity/Vertex/AbstractVVSVertex.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "VVSVertex.fh"

namespace ThePEG {
namespace Helicity {

/** \ingroup Helicity
 *
 * The VVSVertex class is the implementation of the vector-vector-scalar.
 * It inherits from the AbstractVVSVertex class for the storage of the particles and
 * implements the helicity calculations.
 *
 * All interactions of this type should inherit from it and implement the virtual
 * setCoupling member.
 *
 *  The form of the vertex is
 *  \f[icg^{\mu\nu}\epsilon_{1\mu}\epsilon_{2\nu}\f]
 *
 * @see AbstractVVSVertex
 */
class VVSVertex: public AbstractVVSVertex {
  
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
   * @param sca3 The wavefunction for the scalar.
   */
  Complex evaluate(Energy2 q2,const VectorWaveFunction & vec1,
		   const VectorWaveFunction & vec2, const ScalarWaveFunction & sca3);

  /**
   * Evaluate the off-shell vector coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell vector.
   * @param out The ParticleData pointer for the off-shell vector.
   * @param vec2 The wavefunction for the vector.
   * @param sca3 The wavefunction for the scalar.
   * @param mass The mass of the off-shell particle if not taken from the ParticleData
   * object
   * @param width The width of the off-shell particle if not taken from the ParticleData
   * object
   */
  VectorWaveFunction evaluate(Energy2 q2,int iopt,tcPDPtr out,
			      const VectorWaveFunction & vec2,
			      const ScalarWaveFunction & sca3,
			      complex<Energy> mass=-GeV, complex<Energy> width=-GeV);

  /**
   * Evaluate the off-shell scalar coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell scalar.
   * @param out The ParticleData pointer for the off-shell scalar.
   * @param vec1 The wavefunction for the first  vector.
   * @param vec2 The wavefunction for the second vector.
   * @param mass The mass of the off-shell particle if not taken from the ParticleData
   * object
   * @param width The width of the off-shell particle if not taken from the ParticleData
   * object
   */
  ScalarWaveFunction evaluate(Energy2 q2,int iopt, tcPDPtr out,
			      const VectorWaveFunction & vec1,
			      const VectorWaveFunction & vec2,
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
   * Describe an abstract base class with persistent data.
   */
  static AbstractNoPIOClassDescription<VVSVertex> initVVSVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  VVSVertex & operator=(const VVSVertex &);
  
};

}

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of VVSVertex.
 */
template <>
struct BaseClassTrait<ThePEG::Helicity::VVSVertex,1> {
  /** Typedef of the base class of VVSVertex. */
  typedef ThePEG::Helicity::AbstractVVSVertex NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<ThePEG::Helicity::VVSVertex>
  : public ClassTraitsBase<ThePEG::Helicity::VVSVertex> {

  /**
   * Return the class name.
   */
  static string className() { return "ThePEG::VVSVertex"; }
};

/** @endcond */

}


#endif /* ThePEG_VVSVertex_H */
