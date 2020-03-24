// -*- C++ -*-
#ifndef Helicity_AbstractFFFFVertex_H
#define Helicity_AbstractFFFFVertex_H
//
// This is the declaration of the AbstractFFFFVertex class.
//

#include "VertexBase.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "AbstractFFFFVertex.fh"

namespace ThePEG {
namespace Helicity {

/**
 * The AbstractFFFFVertex class provides a base class for all 4-fermion vertices
 * in ThePEG
 */
class AbstractFFFFVertex: public VertexBase {

public:

  /**
   * The default constructor.
   */
  AbstractFFFFVertex() : VertexBase(VertexType::FFFF) {}

  /**
   * Members to calculate the helicity amplitude expressions for vertices
   * and off-shell particles.
   */
  //@{
  /**
   * Evalulate the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param sp1   The wavefunction for the first ferimon.
   * @param sp2   The wavefunction for the second ferimon.
   * @param sbar1 The wavefunction for the first antifermion.
   * @param sbar2 The wavefunction for the second antifermion.
   */
  virtual Complex evaluate(Energy2 q2,
			   const SpinorWaveFunction & sp1,
			   const SpinorBarWaveFunction & sbar1,
			   const SpinorWaveFunction & sp2,
			   const SpinorBarWaveFunction & sbar2) = 0;
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
  AbstractFFFFVertex & operator=(const AbstractFFFFVertex &) = delete;

};

}
}


#endif /* Helicity_AbstractFFFFVertex_H */
