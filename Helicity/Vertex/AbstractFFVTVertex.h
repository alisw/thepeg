// -*- C++ -*-
#ifndef HELICITY_AbstractFFVTVertex_H
#define HELICITY_AbstractFFVTVertex_H
//
// This is the declaration of the AbstractFFVTVertex class.
//

#include "VertexBase.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "AbstractFFVTVertex.fh"

namespace ThePEG {
namespace Helicity {

/**
 * The AbstractFFVTVertex class is the base class for all 
 * fermion-fermion-vector-tensor interactions in ThePEG.
 */
class AbstractFFVTVertex: public VertexBase {

public:
  

  /**
   * Default constructor
   */
  AbstractFFVTVertex() : VertexBase(VertexType::FFVT) {}

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
   * @param ten4  The wavefunction for the tensor.
   */
  virtual Complex evaluate(Energy2 q2,const SpinorWaveFunction & sp1,
			   const SpinorBarWaveFunction & sbar2,
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
  AbstractFFVTVertex & operator=(const AbstractFFVTVertex &) = delete;

};

}
}


namespace ThePEG {

}
#endif /* HELICITY_AbstractFFVTVertex_H */
