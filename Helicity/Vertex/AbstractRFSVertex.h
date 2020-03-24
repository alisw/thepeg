// -*- C++ -*-
#ifndef HELICITY_AbstractRFSVertex_H
#define HELICITY_AbstractRFSVertex_H
//
// This is the declaration of the AbstractRFSVertex class.
//

#include "VertexBase.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/RSSpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/RSSpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "AbstractRFSVertex.fh"

namespace ThePEG {
namespace Helicity {

/**
 * The AbstractRFSVertex class provides a base class for all
 * spin-3/2 fermion-fermion-scalar vertices in ThePEG.
 */
class AbstractRFSVertex: public VertexBase {

public:      

  /**
   * Default constructor
   */
  AbstractRFSVertex() : VertexBase(VertexType::RFS) {}

  /**
   * Members to calculate the helicity amplitude expressions for vertices
   * and off-shell particles.
   */
  //@{
  /**
   * Evalulate the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param sp1   The wavefunction for the RS ferimon.
   * @param sbar2 The wavefunction for the antifermion.
   * @param sca3  The wavefunction for the scalar.
   */
  virtual Complex evaluate(Energy2 q2,const RSSpinorWaveFunction & sp1,
			   const SpinorBarWaveFunction & sbar2,
			   const ScalarWaveFunction & sca3) = 0;
  /**
   * Evalulate the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param sp1   The wavefunction for the ferimon.
   * @param sbar2 The wavefunction for the RS antifermion.
   * @param sca3  The wavefunction for the scalar.
   */
  virtual Complex evaluate(Energy2 q2,const SpinorWaveFunction & sp1,
			   const RSSpinorBarWaveFunction & sbar2,
			   const ScalarWaveFunction & sca3) = 0;

  /**
   * Evaluate the off-shell spinor coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell spinor.
   * @param out The ParticleData pointer for the off-shell spinor.
   * @param sp1   The wavefunction for the ferimon.
   * @param sca3  The wavefunction for the scalar.
   * @param mass The mass of the off-shell particle if not taken from the ParticleData
   * object
   * @param width The width of the off-shell particle if not taken from the ParticleData
   * object
   */
  virtual SpinorWaveFunction evaluate(Energy2 q2,int iopt,tcPDPtr out,
				      const RSSpinorWaveFunction & sp1, 
				      const ScalarWaveFunction & sca3,
				      complex<Energy> mass=-GeV, complex<Energy> width=-GeV) = 0;

  /**
   * Evaluate the off-shell RS spinor coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell spinor.
   * @param out The ParticleData pointer for the off-shell spinor.
   * @param sp1   The wavefunction for the ferimon.
   * @param sca3  The wavefunction for the scalar.
   * @param mass The mass of the off-shell particle if not taken from the ParticleData
   * object
   * @param width The width of the off-shell particle if not taken from the ParticleData
   * object
   */
  virtual RSSpinorWaveFunction evaluate(Energy2 q2,int iopt,tcPDPtr out,
					const SpinorWaveFunction & sp1, 
					const ScalarWaveFunction & sca3,
					complex<Energy> mass=-GeV, complex<Energy> width=-GeV) = 0;

  /**
   * Evaluate the off-shell barred spinor coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell barred spinor.
   * @param out The ParticleData pointer for the off-shell barred spinor.
   * @param sbar2 The wavefunction for the antifermion.
   * @param sca3  The wavefunction for the scalar.
   * @param mass The mass of the off-shell particle if not taken from the ParticleData
   * object
   * @param width The width of the off-shell particle if not taken from the ParticleData
   * object
   */
  virtual SpinorBarWaveFunction evaluate(Energy2 q2,int iopt,tcPDPtr out,
					 const RSSpinorBarWaveFunction & sbar2,
					 const ScalarWaveFunction & sca3,
					 complex<Energy> mass=-GeV, complex<Energy> width=-GeV) = 0;

  /**
   * Evaluate the off-shell barred RS spinor coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell barred spinor.
   * @param out The ParticleData pointer for the off-shell barred spinor.
   * @param sbar2 The wavefunction for the antifermion.
   * @param sca3  The wavefunction for the scalar.
   * @param mass The mass of the off-shell particle if not taken from the ParticleData
   * object
   * @param width The width of the off-shell particle if not taken from the ParticleData
   * object
   */
  virtual RSSpinorBarWaveFunction evaluate(Energy2 q2,int iopt,tcPDPtr out,
					   const SpinorBarWaveFunction & sbar2,
					 const ScalarWaveFunction & sca3,
					   complex<Energy> mass=-GeV, complex<Energy> width=-GeV) = 0;

  /**
   * Evaluate the off-shell scalar coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell scalar.
   * @param out The ParticleData pointer for the off-shell scalar.
   * @param sp1   The wavefunction for the RS ferimon.
   * @param sbar2 The wavefunction for the antifermion.
   * @param mass The mass of the off-shell particle if not taken from the ParticleData
   * object
   * @param width The width of the off-shell particle if not taken from the ParticleData
   * object
   */
  virtual ScalarWaveFunction evaluate(Energy2 q2,int iopt,tcPDPtr out,
				      const RSSpinorWaveFunction & sp1, 
				      const SpinorBarWaveFunction & sbar2,
				      complex<Energy> mass=-GeV, complex<Energy> width=-GeV) = 0;

  /**
   * Evaluate the off-shell scalar coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell scalar.
   * @param out The ParticleData pointer for the off-shell scalar.
   * @param sp1   The wavefunction for the ferimon.
   * @param sbar2 The wavefunction for the RS antifermion.
   * @param mass The mass of the off-shell particle if not taken from the ParticleData
   * object
   * @param width The width of the off-shell particle if not taken from the ParticleData
   * object
   */
  virtual ScalarWaveFunction evaluate(Energy2 q2,int iopt,tcPDPtr out,
				      const SpinorWaveFunction & sp1, 
				      const RSSpinorBarWaveFunction & sbar2,
				      complex<Energy> mass=-GeV, complex<Energy> width=-GeV) = 0;
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
  AbstractRFSVertex & operator=(const AbstractRFSVertex &) = delete;

};

}
}

#endif /* HELICITY_AbstractRFSVertex_H */
