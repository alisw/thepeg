// -*- C++ -*-
#ifndef ThePEG_AbstractRFVSVertex_H
#define ThePEG_AbstractRFVSVertex_H
//
// This is the declaration of the AbstractRFVSVertex class.
//

#include "VertexBase.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/RSSpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/RSSpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "AbstractRFVSVertex.fh"

namespace ThePEG {
namespace Helicity {
  

/**
 * The AbstractRFSVertex class provides a base class for all
 * spin-3/2 fermion-fermion-vector-scalar vertices in ThePEG.
 */
class AbstractRFVSVertex: public VertexBase {

public:

  /**
   * The default constructor.
   */
  AbstractRFVSVertex() : VertexBase(VertexType::RFVS) {}

  /**
   * Members to calculate the helicity amplitude expressions for vertices
   * and off-shell particles.
   */
  //@{
  /**
   * Evalulate the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param sp1   The wavefunction for the RS ferimon.
   * @param sbar2 The wavefunction for the antifermion.
   * @param vec3  The wavefunction for the vector.
   * @param sca4  The wavefunction for the scalar.
   */
  virtual Complex evaluate(Energy2 q2,const RSSpinorWaveFunction & sp1,
			   const SpinorBarWaveFunction & sbar2,
			   const VectorWaveFunction & vec3,
			   const ScalarWaveFunction & sca4) = 0;
  /**
   * Evalulate the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param sp1   The wavefunction for the ferimon.
   * @param sbar2 The wavefunction for the RS antifermion.
   * @param vec3  The wavefunction for the vector.
   * @param sca4  The wavefunction for the scalar.
   */
  virtual Complex evaluate(Energy2 q2,const SpinorWaveFunction & sp1,
			   const RSSpinorBarWaveFunction & sbar2,
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
  AbstractRFVSVertex & operator=(const AbstractRFVSVertex &) = delete;

};

}
}

#endif /* ThePEG_AbstractRFVSVertex_H */
