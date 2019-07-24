// -*- C++ -*-
#ifndef HELICITY_AbstractVVTVertex_H
#define HELICITY_AbstractVVTVertex_H
//
// This is the declaration of the AbstractVVTVertex class.
//

#include "VertexBase.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "AbstractVVTVertex.fh"

namespace ThePEG {
namespace Helicity {

/**
 * Here is the documentation of the AbstractVVTVertex class.
 */
class AbstractVVTVertex: public VertexBase {

public:


  /**
   * Default constructor
   */
  AbstractVVTVertex() : VertexBase(VertexType::VVT) {}

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
   * @param ten3  The wavefunction for the tensor.
   * @param vmass The mass of the vector boson
   */
  virtual Complex evaluate(Energy2 q2,const VectorWaveFunction & vec1,
			   const VectorWaveFunction & vec2,
			   const TensorWaveFunction & ten3,
			   Energy vmass=-GeV) = 0;

  /**
   * Evaluate the off-shell tensor coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell tensor.
   * @param out The ParticleData pointer for the off-shell tensor.
   * @param vec1  The wavefunction for the first  vector.
   * @param vec2  The wavefunction for the second vector.
   * @param mass The mass of the off-shell particle if not taken from the ParticleData
   * object
   * @param width The width of the off-shell particle if not taken from the ParticleData
   * object
   * @param vmass The mass of the vector boson
   */
  virtual TensorWaveFunction evaluate(Energy2 q2,int iopt, tcPDPtr out,
				      const VectorWaveFunction & vec1,
				      const VectorWaveFunction & vec2,
				      Energy vmass=-GeV,
				      complex<Energy> mass=-GeV, complex<Energy> width=-GeV) = 0;

  /**
   * Evaluate the off-shell vector coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell vector.
   * @param out The ParticleData pointer for the off-shell vector.
   * @param vec1  The wavefunction for the first vector.
   * @param ten3  The wavefunction for the tensor.
   * @param mass The mass of the off-shell particle if not taken from the ParticleData
   * object
   * @param width The width of the off-shell particle if not taken from the ParticleData
   * object
   */
  virtual VectorWaveFunction evaluate(Energy2 q2,int iopt, tcPDPtr out,
				      const VectorWaveFunction & vec1,
				      const TensorWaveFunction & ten3,
				      complex<Energy> mass=-GeV, complex<Energy> width=-GeV) = 0;
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
  AbstractVVTVertex & operator=(const AbstractVVTVertex &) = delete;

};

}
}


namespace ThePEG {

}
#endif /* HELICITY_AbstractVVTVertex_H */
