// -*- C++ -*-
#ifndef HELICITY_AbstractFFSTVertex_H
#define HELICITY_AbstractFFSTVertex_H
//
// This is the declaration of the AbstractFFSTVertex class.
//

#include "VertexBase.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "AbstractFFSTVertex.fh"

namespace ThePEG {
namespace Helicity {

/**
 * The AbstractFFSTVertex class is the base class for all 
 * fermion-fermion-scalar-tensor interactions in ThePEG.
 */
class AbstractFFSTVertex: public VertexBase {

public:
  

  /**
   * Default constructor
   */
  AbstractFFSTVertex() : VertexBase(VertexType::FFST) {}

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
   * @param sca3  The wavefunction for the vector.
   * @param ten4  The wavefunction for the tensor.
   */
  virtual Complex evaluate(Energy2 q2,const SpinorWaveFunction & sp1,
			   const SpinorBarWaveFunction & sbar2,
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
  AbstractFFSTVertex & operator=(const AbstractFFSTVertex &) = delete;

};

}
}


namespace ThePEG {

}
#endif /* HELICITY_AbstractFFSTVertex_H */
