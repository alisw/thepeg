// -*- C++ -*-
#ifndef HELICITY_AbstractSSSTVertex_H
#define HELICITY_AbstractSSSTVertex_H
//
// This is the declaration of the AbstractSSSTVertex class.
//

#include "VertexBase.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "AbstractSSSTVertex.fh"

namespace ThePEG {
namespace Helicity {

/**
 * The AbstractSSSTVertex class is the base class for all scalar-scalar-scalar-tensor
 * interactions in ThePEG.
 */
class AbstractSSSTVertex: public VertexBase {

public:

  /**
   * Default constructor
   */
  AbstractSSSTVertex() : VertexBase(VertexType::SSST) {}

  /**
   * Members to calculate the helicity amplitude expressions for vertices
   * and off-shell particles.
   */
  //@{
  /**
   * Evalulate the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param sca1  The wavefunction for the first  scalar.
   * @param sca2  The wavefunction for the second scalar.
   * @param sca3  The wavefunction for the third  scalar.
   * @param ten4  The wavefunction for the tensor.
   */
  virtual Complex evaluate(Energy2 q2,const ScalarWaveFunction & sca1,
			   const ScalarWaveFunction & sca2,
			   const ScalarWaveFunction & sca3, 
			   const TensorWaveFunction & ten4) = 0;
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
  AbstractSSSTVertex & operator=(const AbstractSSSTVertex &) = delete;

};

}
}


namespace ThePEG {

}
#endif /* HELICITY_AbstractSSSTVertex_H */
