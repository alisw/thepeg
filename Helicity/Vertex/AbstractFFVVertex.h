// -*- C++ -*-
#ifndef HELICITY_AbstractFFVVertex_H
#define HELICITY_AbstractFFVVertex_H
//
// This is the declaration of the AbstractFFVVertex class.
//

#include "VertexBase.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "AbstractFFVVertex.fh"

namespace ThePEG {
namespace Helicity {

/**
 * The AbstractFFVVertex class provides a base class for all
 * fermion-fermion-vector vertices in ThePEG.
 */
class AbstractFFVVertex: public VertexBase {

public:

  /**
   *  Enum for the direction in the small angle limit
   */
  enum SmallAngleDirection {
    NegativeZDirection = -1, ///< Along -z 
    PostiveZDirection  =  1  ///< Along +z 
  };
  
public:


  /**
   * Default constructor
   */
  AbstractFFVVertex() : VertexBase(VertexType::FFV) {}

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
			   const VectorWaveFunction & vec3) = 0;

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
					 complex<Energy> mass=-GeV, complex<Energy> width=-GeV) = 0;

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
				      complex<Energy> mass=-GeV, complex<Energy> width=-GeV) = 0;

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
				      complex<Energy> mass=-GeV, complex<Energy> width=-GeV) = 0;
  //@}

  /**
   *  Special members for off-shell fermion wavefunctions with massless
   *  gauge bosons at small angles in the small angle limit for
   *  numerical accuracy. In order to get sufficient accuracy it is
   *  assumed that the fermion lies along either the positive or negative z
   *  axis
   *
   *  Unlike the other members this is not required to be implemented in
   *  all inheriting classes and a default implementation which returns
   *  the general result is included.
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

public:

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  AbstractFFVVertex & operator=(const AbstractFFVVertex &) = delete;

};

}
}


namespace ThePEG {

}
#endif /* HELICITY_AbstractFFVVertex_H */
