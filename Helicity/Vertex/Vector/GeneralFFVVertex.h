// -*- C++ -*-
//
// GeneralFFVVertex.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2019 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_GeneralFFVVertex_H
#define ThePEG_GeneralFFVVertex_H
//
// This is the declaration of the GeneralFFVVertex class.

#include <ThePEG/Helicity/Vertex/AbstractFFVVertex.h>
#include <ThePEG/Helicity/WaveFunction/VectorWaveFunction.h>
#include <ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h>
#include <ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h>
#include "GeneralFFVVertex.fh"

namespace ThePEG {

namespace Helicity{

/** \ingroup Helicity
 *
 *  The GeneralFFVVertex class is the base class for all helicity amplitude
 *  vertices which use a general form of the fermion-fermion-vector vertex.
 *
 *  Any such vertices should inherit from this class and implement the virtual
 *  setcoupling member function. The base AbstractFFVVertex class is
 *  used to store the particles allowed to interact at the vertex.
 *
 *  The form of the vertex is
 *  \f[ic\bar{f_2}\gamma^\mu a^\lambda P_\lambda f_1\epsilon_{3\mu}\f]
 *
 *  @see AbstractFFVVertex
 */
class GeneralFFVVertex: public AbstractFFVVertex {
  
public:
  
  /**
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();

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
  Complex evaluate(Energy2 q2,const SpinorWaveFunction & sp1,
		   const SpinorBarWaveFunction & sbar2,const VectorWaveFunction & vec3);

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
			      const SpinorWaveFunction & sp1,
			      const SpinorBarWaveFunction & sbar2,
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
			      const SpinorWaveFunction & sp1,
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
  const Complex & getLeft() { return _left; }

  /**
   * Get the right coupling.
   */
  const Complex & getRight() { return _right; }

  /**
   *  Get the left coupling for the \f$\sigma^{\mu\nu}\f$ term
   */
  const complex<InvEnergy> getLeftSigma() { return _leftSigma; }

  /**
   *  Get the right coupling for the \f$\sigma^{\mu\nu}\f$ term
   */
  const complex<InvEnergy> getRightSigma() { return _rightSigma; }
  //@}

protected:

  /**
   *  Set the couplings
   */
  //@{
  /**
   * Set the left coupling.
   */
  void setLeft(const Complex & in) { _left = in; }

  /**
   * Set the right coupling.
   */
  void setRight(const Complex & in) { _right = in; }

  /**
   *  Set the left coupling for the \f$\sigma^{\mu\nu}\f$ term
   */
  void  setLeftSigma(const complex<InvEnergy> & in) { _leftSigma  = in; }

  /**
   *  Set the right coupling for the \f$\sigma^{\mu\nu}\f$ term
   */
  void setRightSigma(const complex<InvEnergy> & in) { _rightSigma = in; }
  //@}

private:
  
  /**
   * Private and non-existent assignment operator.
   */
  GeneralFFVVertex & operator=(const GeneralFFVVertex &) = delete;
  
private:

  /**
   *  Left \f$\gamma^\mu\f$ coupling.
   */
  Complex _left;

  /**
   * Right \f$\gamma^\mu\f$ coupling.
   */
  Complex _right;

  /**
   * Left \f$\sigma^{\mu\nu}\f$ coupling
   */
  complex<InvEnergy> _leftSigma;

  /**
   * Right \f$\sigma^{\mu\nu}\f$ coupling
   */
  complex<InvEnergy> _rightSigma;
  
};
}

}
#endif /* ThePEG_GeneralFFVVertex_H */
