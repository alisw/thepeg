// -*- C++ -*-
#ifndef HELICITY_AbstractVVSSVertex_H
#define HELICITY_AbstractVVSSVertex_H
//
// This is the declaration of the AbstractVVSSVertex class.
//

#include "VertexBase.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "AbstractVVSSVertex.fh"

namespace ThePEG {
namespace Helicity {

/**
 * The AbstractVVSSVertex class is the base class for vector-vector-scalar-scalar
 * interactions in ThePEG
 */
class AbstractVVSSVertex: public VertexBase {

public:

  /**
   * Default constructor
   */
  AbstractVVSSVertex() : VertexBase(VertexType::VVSS) {}

  /**
   * Members to calculate the helicity amplitude expressions for vertices
   * and off-shell particles.
   */
  //@{
  /**
   * Evaluate the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param vec1 The wavefunction for the first  vector.
   * @param vec2 The wavefunction for the second vector.
   * @param sca3 The wavefunction for the first  scalar.
   * @param sca4 The wavefunction for the second scalar.
   */
  virtual Complex evaluate(Energy2 q2, const VectorWaveFunction & vec1,
			   const VectorWaveFunction & vec2,
			   const ScalarWaveFunction & sca3,
			   const ScalarWaveFunction & sca4) = 0;

  /**
   * Evaluate the off-shell vector coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell vector.
   * @param out The ParticleData pointer for the off-shell vector.
   * @param vec2 The wavefunction for the second vector.
   * @param sca3 The wavefunction for the first  scalar.
   * @param sca4 The wavefunction for the second scalar.
   * @param mass The mass of the off-shell particle if not taken from the ParticleData
   * object
   * @param width The width of the off-shell particle if not taken from the ParticleData
   * object
   */
  virtual VectorWaveFunction evaluate(Energy2 q2, int iopt,tcPDPtr out,
				      const VectorWaveFunction & vec2,
				      const ScalarWaveFunction & sca3,
				      const ScalarWaveFunction & sca4,
				      complex<Energy> mass=-GeV, complex<Energy> width=-GeV) = 0;

  /**
   * Evaluate the off-shell scalar coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell scalar.
   * @param out The ParticleData pointer for the off-shell scalar.
   * @param vec1 The wavefunction for the first  vector.
   * @param vec2 The wavefunction for the second vector.
   * @param sca3 The wavefunction for the second scalar.
   * @param mass The mass of the off-shell particle if not taken from the ParticleData
   * object
   * @param width The width of the off-shell particle if not taken from the ParticleData
   * object
   */
  virtual ScalarWaveFunction evaluate(Energy2 q2, int iopt,tcPDPtr out,
				      const VectorWaveFunction & vec1,
				      const VectorWaveFunction & vec2,
				      const ScalarWaveFunction & sca3,
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
  static AbstractNoPIOClassDescription<AbstractVVSSVertex> initAbstractVVSSVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  AbstractVVSSVertex & operator=(const AbstractVVSSVertex &);

};

}
}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of AbstractVVSSVertex. */
template <>
struct BaseClassTrait<Helicity::AbstractVVSSVertex,1> {
  /** Typedef of the first base class of AbstractVVSSVertex. */
  typedef Helicity::VertexBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the AbstractVVSSVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Helicity::AbstractVVSSVertex>
  : public ClassTraitsBase<Helicity::AbstractVVSSVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Helicity::AbstractVVSSVertex"; }
};

/** @endcond */

}

#endif /* HELICITY_AbstractVVSSVertex_H */
