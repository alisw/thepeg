// -*- C++ -*-
#ifndef ThePEG_AbstractFFVSVertex_H
#define ThePEG_AbstractFFVSVertex_H
//
// This is the declaration of the AbstractFFVSVertex class.
//

#include "VertexBase.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "AbstractFFVSVertex.fh"

namespace ThePEG {
namespace Helicity {

/**
 * Here is the documentation of the AbstractFFVSVertex class.
 */
class AbstractFFVSVertex: public VertexBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  AbstractFFVSVertex() : VertexBase(VertexType::FFVS) {}
  
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
   * @param sca4  The wavefunction for the scalar.
   */
  virtual Complex evaluate(Energy2 q2,const SpinorWaveFunction & sp1,
			   const SpinorBarWaveFunction & sbar2,
			   const VectorWaveFunction & vec3,
			   const ScalarWaveFunction & sca4) = 0;
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
  AbstractFFVSVertex & operator=(const AbstractFFVSVertex &) = delete;

};

}
}

#endif /* ThePEG_AbstractFFVSVertex_H */
