// -*- C++ -*-
#ifndef Helicity_GeneralVSSVertex_H
#define Helicity_GeneralVSSVertex_H
//
// This is the declaration of the GeneralVSSVertex class.
//

#include "ThePEG/Helicity/Vertex/AbstractVSSVertex.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "GeneralVSSVertex.fh"

namespace ThePEG {
namespace Helicity {

using namespace ThePEG;

/**
 * Here is the documentation of the GeneralVSSVertex class.
 *
 * @see \ref GeneralVSSVertexInterfaces "The interfaces"
 * defined for GeneralVSSVertex.
 */
class GeneralVSSVertex: public AbstractVSSVertex {

public:

  /**
   * The default constructor.
   */
  GeneralVSSVertex() : a_(1.), b_(-1.)
  {}

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

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

protected:

  /**
   *  Member functions top set/get the coefficents
   */
  //@{
  /**
   * Access coefficient of \f$p_2\f$
   */
  Complex a() const {return a_;}
  
  /**
   * Access coefficient of \f$p_3\f$
   */
  Complex b() const {return b_;}
  
  /**
   * Set coefficient of \f$p_2\f$
   */
  void a(Complex in) {a_=in;}
  
  /**
   * Set coefficient of \f$p_3\f$
   */
  void b(Complex in) {b_=in;}
  //@}
  
private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  GeneralVSSVertex & operator=(const GeneralVSSVertex &) = delete;

private:

  /**
   *  The coefficent of \f$p_2\f$
   */
  Complex a_;

  /**
   *  The coefficent of \f$p_3\f$
   */
  Complex b_;

};

}
}

#endif /* Helicity_GeneralVSSVertex_H */
