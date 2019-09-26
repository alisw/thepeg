// -*- C++ -*-
#ifndef HELICITY_AbstractSSSSVertex_H
#define HELICITY_AbstractSSSSVertex_H
//
// This is the declaration of the AbstractSSSSVertex class.
//

#include "VertexBase.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "AbstractSSSSVertex.fh"

namespace ThePEG {
namespace Helicity {

/**
 * The AbstractSSSSVertex class is the base class for all scalar-scalar-scalar
 * interactions in ThePEG
 */
class AbstractSSSSVertex: public VertexBase {

public:

  /**
   * Default constructor
   */
  AbstractSSSSVertex() : VertexBase(VertexType::SSSS) {}

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
   * @param sca4 The wavefunction for the fourth scalar.
   */
  virtual Complex evaluate(Energy2 q2,const ScalarWaveFunction & sca1,
			   const ScalarWaveFunction & sca2,
			   const ScalarWaveFunction & sca3,
			   const ScalarWaveFunction & sca4) = 0;

  /**
   * Evaluate the off-shell scalar coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell scalar.
   * @param out The ParticleData pointer for the off-shell scalar.
   * @param sca1 The wavefunction for the first  scalar.
   * @param sca2 The wavefunction for the second scalar.
   * @param sca3 The wavefunction for the third  scalar.
   */
  virtual ScalarWaveFunction evaluate(Energy2 q2,int iopt, tcPDPtr out, 
				      const ScalarWaveFunction & sca1,
				      const ScalarWaveFunction & sca2,
				      const ScalarWaveFunction & sca3) = 0;
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
  AbstractSSSSVertex & operator=(const AbstractSSSSVertex &) = delete;

};

}
}


namespace ThePEG {

}
#endif /* HELICITY_AbstractSSSSVertex_H */
