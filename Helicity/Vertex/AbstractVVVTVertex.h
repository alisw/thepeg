// -*- C++ -*-
#ifndef HELICITY_AbstractVVVTVertex_H
#define HELICITY_AbstractVVVTVertex_H
//
// This is the declaration of the AbstractVVVTVertex class.
//

#include "VertexBase.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "AbstractVVVTVertex.fh"

namespace ThePEG {
namespace Helicity {

/**
 * The AbstractVVVTVertex class is the base class for all vector-vector-vector-tensor
 * interactions in ThePEG.
 */
class AbstractVVVTVertex: public VertexBase {

public:

  /**
   * Default constructor
   */
  AbstractVVVTVertex() : VertexBase(VertexType::VVVT) {}

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
   * @param vec3  The wavefunction for the third  vector.
   * @param ten4  The wavefunction for the tensor.
   */
  virtual Complex evaluate(Energy2 q2,const VectorWaveFunction & vec1,
			   const VectorWaveFunction & vec2,
			   const VectorWaveFunction & vec3, 
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
  AbstractVVVTVertex & operator=(const AbstractVVVTVertex &) = delete;

};

}
}


namespace ThePEG {

}
#endif /* HELICITY_AbstractVVVTVertex_H */
