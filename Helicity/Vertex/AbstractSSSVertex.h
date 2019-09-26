// -*- C++ -*-
#ifndef HELICITY_AbstractSSSVertex_H
#define HELICITY_AbstractSSSVertex_H
//
// This is the declaration of the AbstractSSSVertex class.
//

#include "VertexBase.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "AbstractSSSVertex.fh"

namespace ThePEG {
namespace Helicity {

/**
 * The AbstractSSSVertex class is the base class for all scalar-scalar-scalar
 * interactions in ThePEG
 */
class AbstractSSSVertex: public VertexBase {

public:

  /**
   * Default constructor
   */
  AbstractSSSVertex() : VertexBase(VertexType::SSS) {}

  /**
   * Members to calculate the helicity amplitude expressions for vertices
   * and off-shell particles.
   */
  //@{  
  /**
   * Evaluate the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param sca1 The wavefunction for the first  scalar.
   * @param sca2 The wavefunction for the second scalar.
   * @param sca3 The wavefunction for the third  scalar.
   */
  virtual Complex evaluate(Energy2 q2,const ScalarWaveFunction & sca1,
			   const ScalarWaveFunction & sca2,
			   const ScalarWaveFunction & sca3) = 0;

  /**
   * Evaluate the off-shell scalar coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell scalar.
   * @param out The ParticleData pointer for the off-shell scalar.
   * @param sca1 The wavefunction for the first  scalar.
   * @param sca2 The wavefunction for the second scalar.
   * @param mass The mass of the off-shell particle if not taken from the ParticleData
   * object
   * @param width The width of the off-shell particle if not taken from the ParticleData
   * object
   */
  virtual ScalarWaveFunction evaluate(Energy2 q2,int iopt, tcPDPtr out, 
				      const ScalarWaveFunction & sca1,
				      const ScalarWaveFunction & sca2,
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
  AbstractSSSVertex & operator=(const AbstractSSSVertex &) = delete;

};

}
}


namespace ThePEG {

}
#endif /* HELICITY_AbstractSSSVertex_H */
