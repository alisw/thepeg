// -*- C++ -*-
//
// VVVSVertex.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2017 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_VVVSVertex_H
#define ThePEG_VVVSVertex_H
//
// This is the declaration of the VVVSVertex class.

#include "ThePEG/Helicity/Vertex/AbstractVVVSVertex.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "VVVSVertex.fh"

namespace ThePEG {
namespace Helicity{
  
/** \ingroup Helicity
 *
 *  The VVVSVertex class is the base class for triple vector scalar vertices.
 *  It inherits from the AbstractVVVSVertex class for the storage of the 
 *  particles allowed at the vertex. Given this vertice has dimenison
 *  5 the two forms possible for either a scalar or pseudoscalar particle.
 *
 *  Classes which implement a specific vertex should inherit from this and
 *  implement the virtual setCoupling member.
 *
 *  The form of the vertex is
 *  \f[ig\left[  (p_1-p_2)^\gamma g^{\alpha\beta }
 *              +(p_2-p_3)^\alpha g^{\beta \gamma}
 *              +(p_3-p_1)^\beta  g^{\alpha\gamma}
 *   \right]\epsilon_{1\alpha}\epsilon_{2\beta}\epsilon_{3\gamma}\f]
 *  for a scalar particle and
 *  \f[ig\epsilon^{\delta\alpha\beta\gamma}\epsilon_{1\alpha}\epsilon_{2\beta}\epsilon_{3\gamma}
 *     (p_1+p_2+p_3)_\delta\f]
 *  @see AbstractVVVSVertex
 */
class VVVSVertex: public AbstractVVVSVertex {

public :
  
  /**
   * Default constructor
   */
  VVVSVertex() : scalar_(true) {}
    
public:
 
  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void persistentInput(PersistentIStream & is, int version);
  //@}

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
   * @param sca  The wavefunction for the scalar particle
   */
  Complex evaluate(Energy2 q2, const VectorWaveFunction & vec1,
		   const VectorWaveFunction & vec2,
		   const VectorWaveFunction & vec3,
		   const ScalarWaveFunction & sca);

  /**
   * Evaluate the off-shell vector coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell vector.
   * @param out The ParticleData pointer for the off-shell vector.
   * @param vec2 The wavefunction for the second vector.
   * @param vec3 The wavefunction for the third  vector.
   * @param sca  The wavefunction for the scalar particle
   * @param mass The mass of the off-shell particle if not taken from the ParticleData
   * object
   * @param width The width of the off-shell particle if not taken from the ParticleData
   * object
   */
  VectorWaveFunction evaluate(Energy2 q2,int iopt, tcPDPtr out,
			      const VectorWaveFunction & vec2,
			      const VectorWaveFunction & vec3,
			      const ScalarWaveFunction & sca,
			      complex<Energy> mass=-GeV, complex<Energy> width=-GeV);

  /**
   * Evaluate the off-shell scalar coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell scalar.
   * @param out The ParticleData pointer for the off-shell scalar.
   * @param vec1 The wavefunction for the first vector.
   * @param vec2 The wavefunction for the second vector.
   * @param vec3 The wavefunction for the third  vector.
   * @param mass The mass of the off-shell particle if not taken from the ParticleData
   * object
   * @param width The width of the off-shell particle if not taken from the ParticleData
   * object
   */
  ScalarWaveFunction evaluate(Energy2 q2,int iopt, tcPDPtr out,
			      const VectorWaveFunction & vec1,
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
   * @param part1 The ParticleData pointer for the first  vector.
   * @param part2 The ParticleData pointer for the second vector.
   * @param part3 The ParticleData pointer for the third  vector.
   * @param part4 The ParticleData pointer for the scalar
   */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1,
			   tcPDPtr part2,tcPDPtr part3,
			   tcPDPtr part4)=0;

  /**
   * Dummy setCouplings for a three point interaction 
   * This method is virtual and must be implemented in 
   * classes inheriting from this.
   */
  virtual void setCoupling(Energy2,tcPDPtr,tcPDPtr,tcPDPtr) {
    assert(false);
  }
  //@}

protected:

  /**
   *  Set the type of the vertex
   */
  void scalar(bool in) {scalar_=in;}

private:
  
  /**
   * Private and non-existent assignment operator.
   */
  VVVSVertex & operator=(const VVVSVertex &) = delete;

  /**
   *   Whether or ont the vertex has a scalar or pseudoscalr particle
   */
  bool scalar_;
  
};

}
}

#endif /* ThePEG_VVVSVertex_H */
