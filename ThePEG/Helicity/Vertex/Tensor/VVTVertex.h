// -*- C++ -*-
//
// VVTVertex.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2011 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_VVTVertex_H
#define ThePEG_VVTVertex_H
//
// This is the declaration of the VVTVertex class.

#include "ThePEG/Helicity/Vertex/AbstractVVTVertex.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "VVTVertex.fh"

namespace ThePEG {
namespace Helicity {

/** \ingroup Helicity
 *
 *  The VVTVertex class is the implementation of the 
 *  vector-vector-tensor vertex. 
 *  It inherits from the AbstractVVTVertex class for the storage of the particles
 *  interacting at the vertex and implements the helicity amplitude calculations.
 *
 *  All implementations of this vertex should inherit from it and implement the
 *  virtual setCoupling member.
 *
 *  The vertex has the form
 *  \f[
 *    \left[m^2_v+k_1\cdot k_2)C_{\mu\nu,\rho\sigma}+D_{\mu\nu,\rho\sigma}\right]
 *    \epsilon_1^\rho\epsilon_2^\sigma \epsilon_3^{\mu\nu}
 *  \f]
 *  where
 * - \f$C_{\mu\nu,\rho\sigma}=g_{\mu\rho}g_{\nu\sigma}+g_{\mu\sigma}g_{\nu\rho}
 *         -g_{\mu\nu}g_{\rho\sigma}\f$
 * - \f$D_{\mu\nu,\rho\sigma}=
 *    g_{\mu\nu}k_{1\sigma}k_{2\rho}-\left[g_{\mu\sigma}k_{1\nu}k_{2\rho}+g_{\mu\rho}k_{1\sigma}k_{2\nu}-g_{\rho\sigma}k_{1\mu}k_{2\nu}+(\mu\leftrightarrow\nu)\right]
 *   \f$
 *
 *  @see AbstractVVTVertex
 */
class VVTVertex: public AbstractVVTVertex {
  
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
   * @param ten3  The wavefunction for the tensor.
   * @param vmass The mass of the vector boson.
   */
  Complex evaluate(Energy2 q2,const VectorWaveFunction & vec1,
		   const VectorWaveFunction & vec2,
		   const TensorWaveFunction & ten3,
		   Energy vmass=-GeV);

  /**
   * Evaluate the off-shell tensor coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell tensor.
   * @param out The ParticleData pointer for the off-shell tensor.
   * @param vec1  The wavefunction for the first  vector.
   * @param vec2  The wavefunction for the second vector.
   * @param vmass The mass of the vector boson.
   * @param mass The mass of the off-shell particle if not taken from the ParticleData
   * object
   * @param width The width of the off-shell particle if not taken from the ParticleData
   * object
   */
  TensorWaveFunction evaluate(Energy2 q2,int iopt, tcPDPtr out,
			      const VectorWaveFunction & vec1,
			      const VectorWaveFunction & vec2,
			      Energy vmass=-GeV,
			      complex<Energy> mass=-GeV, complex<Energy> width=-GeV);

  /**
   * Evaluate the off-shell vector coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell vector.
   * @param out The ParticleData pointer for the off-shell vector.
   * @param vec1  The wavefunction for the first vector.
   * @param ten3  The wavefunction for the tensor.
   * @param mass The mass of the off-shell particle if not taken from the ParticleData
   * object
   * @param width The width of the off-shell particle if not taken from the ParticleData
   * object
   */
  VectorWaveFunction evaluate(Energy2 q2,int iopt, tcPDPtr out,
			      const VectorWaveFunction & vec1,
			      const TensorWaveFunction & ten3,
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
  
private:
  
  /**
   * Describe a concrete class with persistent data.
   */
  static AbstractNoPIOClassDescription<VVTVertex> initVVTVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  VVTVertex & operator=(const VVTVertex &);
  
};

}

/** @cond TRAITSPECIALIZATIONS */
  
/**
 * The following template specialization informs ThePEG about the
 * base class of VVTVertex.
 */
template <>
struct BaseClassTrait<ThePEG::Helicity::VVTVertex,1> {
  /** Typedef of the base class of VVTVertex. */
  typedef ThePEG::Helicity::AbstractVVTVertex NthBase;
};
  
/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<ThePEG::Helicity::VVTVertex>
  : public ClassTraitsBase<ThePEG::Helicity::VVTVertex> {
  
  /**
   * Return the class name.
   */
  static string className() { return "ThePEG::VVTVertex"; }
};

/** @endcond */
  
}


#endif /* ThePEG_VVTVertex_H */
