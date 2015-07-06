// -*- C++ -*-
#ifndef HELICITY_AbstractVVVTVertex_H
#define HELICITY_AbstractVVVTVertex_H
//
// This is the declaration of the AbstractVVVTVertex class.
//

#include "VertexBase.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "AbstractVVVTVertex.fh"

namespace ThePEG {
namespace Helicity {

/**
 * The AbstractVVVTVertex class is the base class for all vector-vector-vector-tensor
 * interactions in ThePEG.
 */
class AbstractVVVTVertex: public VertexBase {

public:

  /**
   * Default constructor
   */
  AbstractVVVTVertex() : VertexBase(VertexType::VVVT) {}

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
   * @param vec3  The wavefunction for the third  vector.
   * @param ten4  The wavefunction for the tensor.
   */
  virtual Complex evaluate(Energy2 q2,const VectorWaveFunction & vec1,
			   const VectorWaveFunction & vec2,
			   const VectorWaveFunction & vec3, 
			   const TensorWaveFunction & ten4) = 0;
  
  /**
   * Evaluate the off-shell tensor coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell tensor.
   * @param out The ParticleData pointer for the off-shell tensor.
   * @param vec1  The wavefunction for the first  vector.
   * @param vec2  The wavefunction for the second vector.
   * @param vec3  The wavefunction for the third  vector.
   * @param mass The mass of the off-shell particle if not taken from the ParticleData
   * object
   * @param width The width of the off-shell particle if not taken from the ParticleData
   * object
   */
  virtual TensorWaveFunction evaluate(Energy2 q2,int iopt, tcPDPtr out,
				      const VectorWaveFunction & vec1,
				      const VectorWaveFunction & vec2,
				      const VectorWaveFunction & vec3,
				      complex<Energy> mass=-GeV, complex<Energy> width=-GeV) = 0;

  /**
   * Evaluate the off-shell vector coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell vector.
   * @param out The ParticleData pointer for the off-shell vector.
   * @param vec1  The wavefunction for the first  vector.
   * @param vec2  The wavefunction for the second vector.
   * @param ten4  The wavefunction for the tensor.
   * @param mass The mass of the off-shell particle if not taken from the ParticleData
   * object
   * @param width The width of the off-shell particle if not taken from the ParticleData
   * object
   */
  virtual VectorWaveFunction evaluate(Energy2 q2,int iopt, tcPDPtr out,
				      const VectorWaveFunction & vec1,
				      const VectorWaveFunction & vec2,
				      const TensorWaveFunction & ten4,
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
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class without persistent data.
   */
  static AbstractNoPIOClassDescription<AbstractVVVTVertex> initAbstractVVVTVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  AbstractVVVTVertex & operator=(const AbstractVVVTVertex &);

};

}
}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of AbstractVVVTVertex. */
template <>
struct BaseClassTrait<Helicity::AbstractVVVTVertex,1> {
  /** Typedef of the first base class of AbstractVVVTVertex. */
  typedef Helicity::VertexBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the AbstractVVVTVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Helicity::AbstractVVVTVertex>
  : public ClassTraitsBase<Helicity::AbstractVVVTVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Helicity::AbstractVVVTVertex"; }
};

/** @endcond */

}

#endif /* HELICITY_AbstractVVVTVertex_H */
