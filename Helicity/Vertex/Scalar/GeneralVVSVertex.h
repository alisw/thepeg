// -*- C++ -*-
#ifndef HELICITY_GeneralVVSVertex_H
#define HELICITY_GeneralVVSVertex_H
//
// This is the declaration of the GeneralVVSVertex class.
//

#include "ThePEG/Helicity/Vertex/AbstractVVSVertex.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "GeneralVVSVertex.fh"

namespace ThePEG {
namespace Helicity {

using namespace ThePEG;

/**
 * The GeneralVVSVertex class implements a
 * general Vector-Vector-Scalar vertex allowing for decay modes 
 * that only enter at the one-loop level
 * 
 * The loop integral is calculated by Passarino-Veltman reduction 
 * and the coefficients are stored here. They must be calculated
 * in the inheriting class along with implementation of the
 * setCoupling member.
 *
 *  The vertex takes the form
 *  \f[  a_{00}g^{\mu\nu}p_1\cdot p_2 + a_{11}p_1^\mu p_1^\nu
 *     + a_{12}p_1^\mu p_2^\nu + a_{21}p_2^\mu p_1^\nu
 *     + a_{22}p_2^\mu p_2^\nu +a_e\epsilon^{p_1 \mu \nu p_2}\f]
 *
 */
class GeneralVVSVertex: public AbstractVVSVertex {

public:

  /**
   * The default constructor.
   */
  GeneralVVSVertex() : _a00(1), 
		       _a11(0), _a12(0), 
		       _a21(0), _a22(0),
		       _aEp(0)
  {}

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
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
  virtual Complex evaluate(Energy2 q2,const VectorWaveFunction & vec1,
			   const VectorWaveFunction & vec2,
			   const ScalarWaveFunction & sca3);

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
  virtual VectorWaveFunction evaluate(Energy2 q2,int iopt,tcPDPtr out,
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
  virtual ScalarWaveFunction evaluate(Energy2 q2,int iopt, tcPDPtr out,
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

public:

  /**@name Set and Get tensor coefficients.*/
  //@{
  /**
   * Access coefficient of \f$g^{\mu\nu}\f$
   */
  Complex a00() const {return _a00;}
      
  /**
   * Access coefficient of \f$p_1^\mu p_1^\nu\f$
   */
  Complex a11() const {return _a11;}
      
  /**
   * Access coefficient of \f$p_1^\mu p_2^\nu\f$
   */
  Complex a12() const {return _a12;}
      
  /**
   * Access coefficient of \f$p_2^\mu p_1^\nu\f$
   */
  Complex a21() const {return _a21;}
      
  /**
   * Access coefficient of \f$p_2^\mu p_2^\nu\f$
   */
  Complex a22() const {return _a22;}
      
  /**
   * Access coefficient of \f$\epsilon^{\mu\nu\alpha\beta}p_1\alpha p_2\beta\f$
   */
  Complex aEp() const {return _aEp;}
      
  /**
   * Set tensor coefficient of \f$g^{\mu\nu}\f$
   */
  void a00(const Complex & val) {_a00 = val;}
      
  /**
   * Set tensor coefficient of \f$p_1^\mu p_1^\nu\f$
   */
  void a11(const Complex & val) {_a11 = val;}
      
  /**
   * Set tensor coefficient of \f$p_1^\mu p_2^\nu\f$
   */
  void a12(const Complex & val) {_a12 = val;}
      
  /**
   * Set tensor coefficient of \f$p_2^\mu p_1^\nu\f$
   */
  void a21(const Complex & val) {_a21 = val;}
      
  /**
   * Set tensor coefficient of \f$p_2^\mu p_2^\nu\f$
   */ 
  void a22(const Complex & val) {_a22 = val;}
      
  /**
   * Set tensor coefficient of \f$\epsilon^{\mu\nu\alpha\beta}p_1\alpha p_2\beta\f$
   */
  void aEp(const Complex & val) {_aEp = val;}
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  GeneralVVSVertex & operator=(const GeneralVVSVertex &) = delete;

      
private:
  
  /** @name Store tensor coefficients.*/
  //@{
  /**
   * Coefficient of \f$g^{\mu\nu}\f$
   */
  Complex _a00;

  /**
   * Coefficient of \f$p_1^\mu p_1^\nu\f$
   */
  Complex _a11;

  /**
   * Coefficient of \f$p_1^\mu p_2^\nu\f$
   */
  Complex _a12;

  /**
   * Coefficient of \f$p_2^\mu p_1^\nu\f$
   */
  Complex _a21;

  /**
   * Coefficient of \f$p_2^\mu p_2^\nu\f$
   */
  Complex _a22;

  /**
   * Coefficient of \f$\epsilon^{\mu\nu\alpha\beta}p_1\alpha p_2\beta\f$
   */
  Complex _aEp;
  //@}
};

}
}
#endif /* HELICITY_GeneralVVSVertex_H */
