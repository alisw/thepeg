// -*- C++ -*-
#ifndef HELICITY_AbstractVVVVVertex_H
#define HELICITY_AbstractVVVVVertex_H
//
// This is the declaration of the AbstractVVVVVertex class.
//

#include "VertexBase.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "AbstractVVVVVertex.fh"

namespace ThePEG {
namespace Helicity {

/**
 * The AbstractVVVVVertex class is the base class for
 * vector-vector-vector-vector interactions in ThePEG
 */
class AbstractVVVVVertex: public VertexBase {

public:

  /**
   * Default constructor
   */
  AbstractVVVVVertex() : VertexBase(VertexType::VVVV) {}
  
  /**
   * Members to calculate the helicity amplitude expressions for vertices
   * and off-shell particles.
   */
  //@{
  /**
   * Evaluate the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Evaluation option, 0 just evaluate the four point vertex, 1
   * include all the three point diagrams as well.
   * @param vec1 The wavefunction for the first  vector.
   * @param vec2 The wavefunction for the second vector.
   * @param vec3 The wavefunction for the third  vector.
   * @param vec4 The wavefunction for the fourth vector.
   */
  virtual Complex evaluate(Energy2 q2, int iopt,
			   const VectorWaveFunction & vec1, 
			   const VectorWaveFunction & vec2,
			   const VectorWaveFunction & vec3, 
			   const VectorWaveFunction & vec4) = 0;
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
  AbstractVVVVVertex & operator=(const AbstractVVVVVertex &) = delete;

};

}
}


namespace ThePEG {

}
#endif /* HELICITY_AbstractVVVVVertex_H */
