// -*- C++ -*-
#ifndef ThePEG_AbstractRFSSVertex_H
#define ThePEG_AbstractRFSSVertex_H
//
// This is the declaration of the AbstractRFSSVertex class.
//

#include "ThePEG/Helicity/Vertex/VertexBase.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/RSSpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/RSSpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "AbstractRFSSVertex.fh"

namespace ThePEG {
namespace Helicity {
/**
 * The AbstractRFSVertex class provides a base class for all
 * spin-3/2 fermion-fermion-scalar-scalar vertices in ThePEG.
 */
class AbstractRFSSVertex: public VertexBase {

public:
  
  /**
   * The default constructor.
   */
  AbstractRFSSVertex() : VertexBase(VertexType::RFSS) {}

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
   * @param sca3  The wavefunction for the 1st scalar.
   * @param sca4  The wavefunction for the 2nd scalar.
   */
  virtual Complex evaluate(Energy2 q2,const RSSpinorWaveFunction & sp1,
			   const SpinorBarWaveFunction & sbar2,
			   const ScalarWaveFunction & sca3,
			   const ScalarWaveFunction & sca4) = 0;
  /**
   * Evalulate the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param sp1   The wavefunction for the ferimon.
   * @param sbar2 The wavefunction for the RS antifermion.
   * @param sca3  The wavefunction for the 1st scalar.
   * @param sca4  The wavefunction for the 2nd scalar.
   */
  virtual Complex evaluate(Energy2 q2,const SpinorWaveFunction & sp1,
			   const RSSpinorBarWaveFunction & sbar2,
			   const ScalarWaveFunction & sca3,
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
  AbstractRFSSVertex & operator=(const AbstractRFSSVertex &) = delete;

};

}
}

#endif /* ThePEG_AbstractRFSSVertex_H */
