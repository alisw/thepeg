// -*- C++ -*-
//
// RFVVertex.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2017 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_RFVVertex_H
#define ThePEG_RFVVertex_H
//
// This is the declaration of the RFVVertex class.

#include <ThePEG/Helicity/Vertex/AbstractRFVVertex.h>
#include <ThePEG/Helicity/WaveFunction/VectorWaveFunction.h>
#include <ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h>
#include <ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h>
#include "RFVVertex.fh"

namespace ThePEG {

namespace Helicity{

/** \ingroup Helicity
 *
 *  The RFVVertex class is the base class for all helicity amplitude
 *  vertices which use the renormalisable form for the 
 *  spin-3/2 fermion-fermion-vector vertex. 
 *
 *  Any such vertices should inherit from this class and implement the virtual
 *  setcoupling member function. The base AbstractRFVVertex class is used to store the
 *  particles allowed to interact at the vertex.
 *
 *  The form of the vertex is
 *  \f[ic\bar{f_2}_\mu \left[ g^{\mu\nu} a^{1\lambda} P_\lambda
 *                           +\gamma^\nu p_{1}^\mu a^{2\lambda} P_\lambda
 *                           +p_{1}^\mu  p_{2}^\nu a^{3\lambda} P_\lambda 
 * \right]f_1\epsilon_{3\nu}\f]
 *
 *  @see AbstractRFVVertex
 */
class RFVVertex: public AbstractRFVVertex {
  
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
   * @param vec3  The wavefunction for the vector.
   */
  Complex evaluate(Energy2 q2,const RSSpinorWaveFunction & sp1,
		   const SpinorBarWaveFunction & sbar2,
		   const VectorWaveFunction & vec3);

  /**
   * Evalulate the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param sp1   The wavefunction for the ferimon.
   * @param sbar2 The wavefunction for the antifermion.
   * @param vec3  The wavefunction for the vector.
   */
  Complex evaluate(Energy2 q2,const SpinorWaveFunction & sp1,
		   const RSSpinorBarWaveFunction & sbar2,
		   const VectorWaveFunction & vec3);

  /**
   * Evaluate the off-shell barred spinor coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell barred spinor.
   * @param out The ParticleData pointer for the off-shell barred spinor.
   * @param sbar2 The wavefunction for the antifermion.
   * @param vec3  The wavefunction for the vector.
   * @param mass The mass of the off-shell particle if not taken from the ParticleData
   * object
   * @param width The width of the off-shell particle if not taken from the ParticleData
   * object
   */
  SpinorBarWaveFunction evaluate(Energy2 q2,int iopt,tcPDPtr out,
				 const RSSpinorBarWaveFunction & sbar2,
				 const VectorWaveFunction & vec3,
				 complex<Energy> mass=-GeV, complex<Energy> width=-GeV);

  /**
   * Evaluate the off-shell barred spinor coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell barred spinor.
   * @param out The ParticleData pointer for the off-shell barred spinor.
   * @param sbar2 The wavefunction for the antifermion.
   * @param vec3  The wavefunction for the vector.
   * @param mass The mass of the off-shell particle if not taken from the ParticleData
   * object
   * @param width The width of the off-shell particle if not taken from the ParticleData
   * object
   */
  RSSpinorBarWaveFunction evaluate(Energy2 q2,int iopt,tcPDPtr out,
				   const SpinorBarWaveFunction & sbar2,
				   const VectorWaveFunction & vec3,
				   complex<Energy> mass=-GeV, complex<Energy> width=-GeV);

  /**
   * Evaluate the off-shell vector coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell vector.
   * @param out The ParticleData pointer for the off-shell vector.
   * @param sp1   The wavefunction for the ferimon.
   * @param sbar2 The wavefunction for the antifermion.
   * @param mass The mass of the off-shell particle if not taken from the ParticleData
   * object
   * @param width The width of the off-shell particle if not taken from the ParticleData
   * object
   */
  VectorWaveFunction evaluate(Energy2 q2,int iopt,tcPDPtr out,
			      const RSSpinorWaveFunction & sp1,
			      const SpinorBarWaveFunction & sbar2,
			      complex<Energy> mass=-GeV, complex<Energy> width=-GeV);

  /**
   * Evaluate the off-shell vector coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell vector.
   * @param out The ParticleData pointer for the off-shell vector.
   * @param sp1   The wavefunction for the ferimon.
   * @param sbar2 The wavefunction for the antifermion.
   * @param mass The mass of the off-shell particle if not taken from the ParticleData
   * object
   * @param width The width of the off-shell particle if not taken from the ParticleData
   * object
   */
  VectorWaveFunction evaluate(Energy2 q2,int iopt,tcPDPtr out,
			      const SpinorWaveFunction & sp1,
			      const RSSpinorBarWaveFunction & sbar2,
			      complex<Energy> mass=-GeV, complex<Energy> width=-GeV);

  /**
   * Evaluate the off-shell spinor coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell spinor.
   * @param out The ParticleData pointer for the off-shell spinor.
   * @param sp1   The wavefunction for the ferimon.
   * @param vec3  The wavefunction for the vector.
   * @param mass The mass of the off-shell particle if not taken from the ParticleData
   * object
   * @param width The width of the off-shell particle if not taken from the ParticleData
   * object
   */
  RSSpinorWaveFunction evaluate(Energy2 q2,int iopt,tcPDPtr out,
				const SpinorWaveFunction & sp1,
				const VectorWaveFunction & vec3,
				complex<Energy> mass=-GeV, complex<Energy> width=-GeV);

  /**
   * Evaluate the off-shell spinor coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell spinor.
   * @param out The ParticleData pointer for the off-shell spinor.
   * @param sp1   The wavefunction for the ferimon.
   * @param vec3  The wavefunction for the vector.
   * @param mass The mass of the off-shell particle if not taken from the ParticleData
   * object
   * @param width The width of the off-shell particle if not taken from the ParticleData
   * object
   */
  SpinorWaveFunction evaluate(Energy2 q2,int iopt,tcPDPtr out,
			      const RSSpinorWaveFunction & sp1,
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

  /**
   * Get the Couplings 
   */
  //@{
  /**
   * Get the left coupling.
   */
  const vector<Complex> &  left() const { return _left; }

  /**
   * Get the right coupling.
   */
  const vector<Complex> & right() const { return _right; }
  //@}

protected:

  /**
   *  Set the couplings
   */
  //@{
  /**
   * Set the left coupling.
   */
  void left(const vector<Complex> & in) { _left = in; }

  /**
   * Set the right coupling.
   */
  void right(const vector<Complex> & in) { _right = in; }
  //@}

private:
  
  /**
   * Describe an abstract class with persistent data.
   */
  static AbstractNoPIOClassDescription<RFVVertex> initRFVVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  RFVVertex & operator=(const RFVVertex &);
  
private:

  /**
   * Left coupling.
   */
  vector<Complex> _left;

  /**
   * Right coupling.
   */
  vector<Complex> _right;
  
};
}

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of RFVVertex.
 */
template <>
struct BaseClassTrait<ThePEG::Helicity::RFVVertex,1> {
  /** Typedef of the base class of RFVVertex. */
  typedef ThePEG::Helicity::AbstractRFVVertex NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<ThePEG::Helicity::RFVVertex>
  : public ClassTraitsBase<ThePEG::Helicity::RFVVertex> {

  /**
   * Return the class name.
   */
  static string className() { return "ThePEG::RFVVertex"; }
};

/** @endcond */

}

#endif /* ThePEG_RFVVertex_H */
