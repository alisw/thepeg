// -*- C++ -*-
//
// SSSVertex.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2011 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_SSSVertex_H
#define ThePEG_SSSVertex_H
//
// This is the declaration of the SSSVertex class.
//
#include "ThePEG/Helicity/Vertex/AbstractSSSVertex.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "SSSVertex.fh"

namespace ThePEG {
namespace Helicity {

/** \ingroup Helicity
 *  
 *  The SSSVertex class is the implementation of the interaction of
 *  three scalars. It inherits from the AbstractSSSVertex class for the storage of the
 *  particles interacting at the vertex and implements the helicity calculations.
 *
 *  Any classes implementating the vertex should inherit from it and implement
 *  the virtual set Coupling member.
 *
 *  The form of the vertex is
 * \f[ic\phi_1\phi_2\phi_3\f]
 *
 *  @see AbstractSSSVertex
 */
class SSSVertex: public AbstractSSSVertex {

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
   * @param sca1 The wavefunction for the first  scalar.
   * @param sca2 The wavefunction for the second scalar.
   * @param sca3 The wavefunction for the third  scalar.
   */
  Complex evaluate(Energy2 q2,const ScalarWaveFunction & sca1,
		   const ScalarWaveFunction & sca2,const ScalarWaveFunction & sca3);

  /**
   * Evaluate the off-shell scalar coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell scalar.
   * @param out The ParticleData pointer for the off-shell scalar.
   * @param sca1 The wavefunction for the first  scalar.
   * @param sca2 The wavefunction for the second scalar.
   * @param mass The mass of the off-shell particle if not taken from the ParticleData
   * object
   * @param width The width of the off-shell particle if not taken from the ParticleData
   * object
   */
  ScalarWaveFunction evaluate(Energy2 q2,int iopt, tcPDPtr out, 
			      const ScalarWaveFunction & sca1,
			      const ScalarWaveFunction & sca2,
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
  static AbstractNoPIOClassDescription<SSSVertex> initSSSVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  SSSVertex & operator=(const SSSVertex &);
  
};

}

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of SSSVertex.
 */
template <>
struct BaseClassTrait<ThePEG::Helicity::SSSVertex,1> {
  /** Typedef of the base class of SSSVertex. */
  typedef ThePEG::Helicity::AbstractSSSVertex NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<ThePEG::Helicity::SSSVertex>
  : public ClassTraitsBase<ThePEG::Helicity::SSSVertex> {

  /**
   * Return the class name.
   */
  static string className() { return "ThePEG::SSSVertex"; }
};

/** @endcond */

}


#endif /* ThePEG_SSSVertex_H */
