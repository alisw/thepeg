// -*- C++ -*-
#ifndef HELICITY_AbstractFFVVVertex_H
#define HELICITY_AbstractFFVVVertex_H
//
// This is the declaration of the AbstractFFVVVertex class.
//

#include "VertexBase.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "AbstractFFVVVertex.fh"

namespace ThePEG {
namespace Helicity {

/**
 * The AbstractFFVVVertex class provides a base class for all
 * fermion-fermion-vector-vector vertices in ThePEG.
 */
class AbstractFFVVVertex: public VertexBase {

public:

  /**
   * Default constructor
   */
  AbstractFFVVVertex() : VertexBase(VertexType::FFVV) {}

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
   * @param vec4  The wavefunction for the vector.
   */
  virtual Complex evaluate(Energy2 q2,const SpinorWaveFunction & sp1,
			   const SpinorBarWaveFunction & sbar2,
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
  AbstractFFVVVertex & operator=(const AbstractFFVVVertex &) = delete;

};

}
}

#endif /* HELICITY_AbstractFFVVVertex_H */
