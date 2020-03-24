// -*- C++ -*-
#ifndef ThePEG_AbstractRFVVVertex_H
#define ThePEG_AbstractRFVVVertex_H
//
// This is the declaration of the AbstractRFVVVertex class.
//

#include "VertexBase.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/RSSpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/RSSpinorBarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "AbstractRFVVVertex.fh"

namespace ThePEG {
namespace Helicity {
  

/**
 * The AbstractRFSVertex class provides a base class for all
 * spin-3/2 fermion-fermion-vector-vector vertices in ThePEG.
 */
class AbstractRFVVVertex: public VertexBase {

public:

  /**
   * The default constructor.
   */
  AbstractRFVVVertex() : VertexBase(VertexType::RFVV) {}

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
   * @param vec3  The wavefunction for the 1st vector.
   * @param vec4  The wavefunction for the 2nd vector.
   */
  virtual Complex evaluate(Energy2 q2,const RSSpinorWaveFunction & sp1,
			   const SpinorBarWaveFunction & sbar2,
			   const VectorWaveFunction & vec3,
			   const VectorWaveFunction & vec4) = 0;
  /**
   * Evalulate the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param sp1   The wavefunction for the ferimon.
   * @param sbar2 The wavefunction for the RS antifermion.
   * @param vec3  The wavefunction for the 1st vector.
   * @param vec4  The wavefunction for the 2nd vector.
   */
  virtual Complex evaluate(Energy2 q2,const SpinorWaveFunction & sp1,
			   const RSSpinorBarWaveFunction & sbar2,
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
  AbstractRFVVVertex & operator=(const AbstractRFVVVertex &) = delete;

};

}
}

#endif /* ThePEG_AbstractRFVVVertex_H */
