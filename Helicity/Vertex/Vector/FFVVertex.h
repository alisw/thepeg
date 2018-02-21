// -*- C++ -*-
//
// FFVVertex.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2017 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_FFVVertex_H
#define ThePEG_FFVVertex_H
//
// This is the declaration of the FFVVertex class.

#include <ThePEG/Helicity/Vertex/AbstractFFVVertex.h>
#include <ThePEG/Helicity/WaveFunction/VectorWaveFunction.h>
#include <ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h>
#include <ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h>
#include "FFVVertex.fh"

namespace ThePEG {

namespace Helicity{

/** \ingroup Helicity
 *
 *  The FFVVertex class is the base class for all helicity amplitude
 *  vertices which use the renormalisable form for the 
 *  fermion-fermion-vector vertex. 
 *
 *  Any such vertices should inherit from this class and implement the virtual
 *  setcoupling member function. The base AbstractFFVVertex class is used to store the
 *  particles allowed to interact at the vertex.
 *
 *  The form of the vertex is
 *  \f[ic\bar{f_2}\gamma^\mu a^\lambda P_\lambda f_1\epsilon_{3\mu}\f]
 *
 *  @see AbstractFFVVertex
 */
class FFVVertex: public AbstractFFVVertex {
  
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
  virtual Complex evaluate(Energy2 q2,const SpinorWaveFunction & sp1,
			   const SpinorBarWaveFunction & sbar2,
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
  virtual SpinorBarWaveFunction evaluate(Energy2 q2,int iopt,tcPDPtr out,
					 const SpinorBarWaveFunction & sbar2,
					 const VectorWaveFunction & vec3,
					 complex<Energy> mass=-GeV,
					 complex<Energy> width=-GeV);

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
  virtual VectorWaveFunction evaluate(Energy2 q2,int iopt,tcPDPtr out,
				      const SpinorWaveFunction & sp1,
				      const SpinorBarWaveFunction & sbar2,
				      complex<Energy> mass=-GeV,
				      complex<Energy> width=-GeV);

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
  virtual SpinorWaveFunction evaluate(Energy2 q2,int iopt,tcPDPtr out,
				      const SpinorWaveFunction & sp1,
				      const VectorWaveFunction & vec3,
				      complex<Energy> mass=-GeV,
				      complex<Energy> width=-GeV);
  //@}

  /**
   *  Special members for off-shell fermion wavefunctions with massless
   *  gauge bosons at small angles in the small angle limit for
   *  numerical accuracy. In order to get sufficient accuracy it is
   *  assumed that the fermion lies along either the positive or negative z
   *  axis.
   */
  //@{
  /** Small angle approx for an off-shell spinor
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell spinor.
   * @param out The ParticleData pointer for the off-shell spinor.
   * @param sp1   The wavefunction for the ferimon.
   * @param vec3  The wavefunction for the vector.
   * @param fhel Helicity of the fermion
   * @param vhel Helicity of the vector
   * @param ctheta   The cosine of the 
   *                 polar angle of the photon with respect to the fermion
   * @param phi      The azimuthal angle of the photon with respect to the fermion
   * @param stheta   The sine of the
   *                 polar angle of the photon with respect to the fermion
   * @param includeEikonal Whether or not to include the eikonal piece
   * @param direction Whether fermion along + or - z direction
   * @param mass The mass of the off-shell particle if not taken from the ParticleData
   * object
   * @param width The width of the off-shell particle if not taken from the ParticleData
   * object
   */
  virtual SpinorWaveFunction evaluateSmall(Energy2 q2,int iopt, tcPDPtr out,
					   const SpinorWaveFunction & sp1,
					   const VectorWaveFunction & vec3,
					   unsigned int fhel, unsigned int vhel,
					   double ctheta, double phi, double stheta,
					   bool includeEikonal = true,
					   SmallAngleDirection direction = PostiveZDirection,
					   Energy mass=-GeV, Energy width=-GeV);

  /** Small angle approx for an off-shell spinor
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell spinor.
   * @param out The ParticleData pointer for the off-shell spinor.
   * @param sbar2 The wavefunction for the antifermion.
   * @param vec3  The wavefunction for the vector.
   * @param fhel Helicity of the fermion
   * @param vhel Helicity of the vector
   * @param ctheta   The cosine of the 
   *                 polar angle of the photon with respect to the fermion
   * @param phi      The azimuthal angle of the photon with respect to the fermion
   * @param stheta   The sine of the
   *                 polar angle of the photon with respect to the fermion
   * @param includeEikonal Whether or not to include the eikonal piece
   * @param direction Whether fermion along + or - z direction
   * @param mass The mass of the off-shell particle if not taken from the ParticleData
   * object
   * @param width The width of the off-shell particle if not taken from the ParticleData
   * object
   */
  virtual SpinorBarWaveFunction evaluateSmall(Energy2 q2,int iopt, tcPDPtr out,
					      const SpinorBarWaveFunction & sbar2,
					      const VectorWaveFunction & vec3,
					      unsigned int fhel, unsigned int vhel,
					      double ctheta, double phi, double stheta,
					      bool includeEikonal = true,
					      SmallAngleDirection direction = PostiveZDirection,
					      Energy mass=-GeV, Energy width=-GeV);
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
  const Complex & left() const { return _left; }

  /**
   * Get the right coupling.
   */
  const Complex & right() const { return _right; }
  //@}

protected:

  /**
   *  Set the couplings
   */
  //@{
  /**
   * Set the left coupling.
   */
  void left(const Complex & in) { _left = in; }

  /**
   * Set the right coupling.
   */
  void right(const Complex & in) { _right = in; }
  //@}

private:
  
  /**
   * Describe an abstract class with persistent data.
   */
  static AbstractNoPIOClassDescription<FFVVertex> initFFVVertex;
  
  /**
   * Private and non-existent assignment operator.
   */
  FFVVertex & operator=(const FFVVertex &);
  
private:

  /**
   * Left coupling.
   */
  Complex _left;

  /**
   * Right coupling.
   */
  Complex _right;
  
};
}

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of FFVVertex.
 */
template <>
struct BaseClassTrait<ThePEG::Helicity::FFVVertex,1> {
  /** Typedef of the base class of FFVVertex. */
  typedef ThePEG::Helicity::AbstractFFVVertex NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<ThePEG::Helicity::FFVVertex>
  : public ClassTraitsBase<ThePEG::Helicity::FFVVertex> {

  /**
   * Return the class name.
   */
  static string className() { return "ThePEG::FFVVertex"; }
};

/** @endcond */

}


#endif /* ThePEG_FFVVertex_H */
