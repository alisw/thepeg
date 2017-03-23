// -*- C++ -*-
//
// FFTVertex.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2011 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_FFTVertex_H
#define ThePEG_FFTVertex_H
//
// This is the declaration of the FFTVertex class.
//
#include "ThePEG/Helicity/Vertex/AbstractFFTVertex.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "FFTVertex.fh"

namespace ThePEG {
namespace Helicity {

/** \ingroup Helicity
 *
 *  The FFTVertex class is the implementation of the fermion-fermion-tensor
 *  vertex. It inherits from the AbstractFFTVertex class for the storage of the particles
 *  interacting at the vertex and implements the helicity amplitude calculations.
 *
 *  All implementations of this vertex should inherit from it and implement the
 *  virtual setCoupling member.
 *
 *  The vertex has the form
 * \f[-\frac{i\kappa}8\bar{f_2}\left[
 *  \gamma_\mu(k_1-k_2)_\nu+\gamma_\nu(k_1-k_2)_\mu
 * -2g_{\mu\nu}(k\!\!\!\!\!\not\,\,\,_1-k\!\!\!\!\!\not\,\,\,_2)+4g_{\mu\nu}m_{f}
 * \right]f_1\epsilon^{\mu\nu}_3\f] 
 *
 *  @see AbstractFFTVertex
 */
class FFTVertex: public AbstractFFTVertex {
      
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
   * @param sp1   The wavefunction for the ferimon.
   * @param sbar2 The wavefunction for the antifermion.
   * @param ten3  The wavefunction for the tensor.
   */
  Complex evaluate(Energy2 q2,const SpinorWaveFunction & sp1,
		   const SpinorBarWaveFunction & sbar2, 
		   const TensorWaveFunction & ten3);

  /**
   * Evaluate the off-shell tensor coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell tensor.
   * @param out The ParticleData pointer for the off-shell tensor.
   * @param sp1   The wavefunction for the ferimon.
   * @param sbar2 The wavefunction for the antifermion.
   * @param mass The mass of the off-shell particle if not taken from the ParticleData
   * object
   * @param width The width of the off-shell particle if not taken from the ParticleData
   * object
   */
  TensorWaveFunction evaluate(Energy2 q2, int iopt,tcPDPtr out,
			      const SpinorWaveFunction & sp1,
			      const SpinorBarWaveFunction & sbar2,
			      complex<Energy> mass=-GeV, complex<Energy> width=-GeV);

  /**
   * Evaluate the off-shell spinor coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell spinor.
   * @param out The ParticleData pointer for the off-shell spinor.
   * @param sp1   The wavefunction for the ferimon.
   * @param ten3  The wavefunction for the tensor.
   * @param mass The mass of the off-shell particle if not taken from the ParticleData
   * object
   * @param width The width of the off-shell particle if not taken from the ParticleData
   * object
   */
  SpinorWaveFunction evaluate(Energy2 q2,int iopt,tcPDPtr out,
			      const SpinorWaveFunction & sp1,
			      const TensorWaveFunction & ten3,
			      complex<Energy> mass=-GeV, complex<Energy> width=-GeV);

  /**
   * Evaluate the off-shell barred spinor coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell barred spinor.
   * @param out The ParticleData pointer for the off-shell barred spinor.
   * @param sbar2 The wavefunction for the antifermion.
   * @param ten3  The wavefunction for the tensor.
   * @param mass The mass of the off-shell particle if not taken from the ParticleData
   * object
   * @param width The width of the off-shell particle if not taken from the ParticleData
   * object
   */
  SpinorBarWaveFunction evaluate(Energy2 q2,int iopt, tcPDPtr out,
				 const SpinorBarWaveFunction & sbar2,
				 const TensorWaveFunction& ten3,
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
   * Describe an abstract class with persistent data.
   */
  static AbstractNoPIOClassDescription<FFTVertex> initFFTVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  FFTVertex & operator=(const FFTVertex &);
  
};

}

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of FFTVertex.
 */
template <>
struct BaseClassTrait<ThePEG::Helicity::FFTVertex,1> {
  /** Typedef of the base class of FFTVertex. */
  typedef ThePEG::Helicity::AbstractFFTVertex NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<ThePEG::Helicity::FFTVertex>
  : public ClassTraitsBase<ThePEG::Helicity::FFTVertex> {
  
  /**
   * Return the class name.
   */
  static string className() { return "ThePEG::FFTVertex"; }
 
};
  
/** @endcond */
  
}

#endif /* ThePEG_FFTVertex_H */
