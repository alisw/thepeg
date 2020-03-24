// -*- C++ -*-
#ifndef HELICITY_AbstractVVSTVertex_H
#define HELICITY_AbstractVVSTVertex_H
//
// This is the declaration of the AbstractVVSTVertex class.
//

#include "VertexBase.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "AbstractVVSTVertex.fh"

namespace ThePEG {
namespace Helicity {

/**
 * The AbstractVVSTVertex class is the base class for all vector-vector-scalar-tensor
 * interactions in ThePEG.
 */
class AbstractVVSTVertex: public VertexBase {

public:

  /**
   * Default constructor
   */
  AbstractVVSTVertex() : VertexBase(VertexType::VVST) {}

  /**
   * Members to calculate the helicity amplitude expressions for vertices
   * and off-shell particles.
   */
  //@{
  /**
   * Evalulate the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param vec1  The wavefunction for the first  vector.
   * @param vec2  The wavefunction for the second vector.
   * @param sca3  The wavefunction for the third  vector.
   * @param ten4  The wavefunction for the tensor.
   */
  virtual Complex evaluate(Energy2 q2,const VectorWaveFunction & vec1,
			   const VectorWaveFunction & vec2,
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
  AbstractVVSTVertex & operator=(const AbstractVVSTVertex &) = delete;

};

}
}


namespace ThePEG {

}
#endif /* HELICITY_AbstractVVSTVertex_H */
