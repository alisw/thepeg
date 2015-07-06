// -*- C++ -*-
#ifndef HELICITY_AbstractVVVVertex_H
#define HELICITY_AbstractVVVVertex_H
//
// This is the declaration of the AbstractVVVVertex class.
//

#include "VertexBase.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "AbstractVVVVertex.fh"

namespace ThePEG {
namespace Helicity {

/**
 * The AbstractVVVVertex class provides the base class for all
 * vector-vector-vector interactions in ThePEG
 */
class AbstractVVVVertex: public VertexBase {

public:
  
  /**
   * Default constructor
   */
  AbstractVVVVertex() : VertexBase(VertexType::VVV) {}

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
   * @param vec3 The wavefunction for the third  vector.
   */
  virtual Complex evaluate(Energy2 q2, const VectorWaveFunction & vec1,
			   const VectorWaveFunction & vec2,
			   const VectorWaveFunction & vec3) = 0;

  /**
   * Evaluate the off-shell vector coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell vector.
   * @param out The ParticleData pointer for the off-shell vector.
   * @param vec2 The wavefunction for the second vector.
   * @param vec3 The wavefunction for the third  vector.
   * @param mass The mass of the off-shell particle if not taken from the ParticleData
   * object
   * @param width The width of the off-shell particle if not taken from the ParticleData
   * object
   */
  virtual VectorWaveFunction evaluate(Energy2 q2,int iopt, tcPDPtr out,
				      const VectorWaveFunction & vec2,
				      const VectorWaveFunction & vec3,
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
  static AbstractNoPIOClassDescription<AbstractVVVVertex> initAbstractVVVVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  AbstractVVVVertex & operator=(const AbstractVVVVertex &);

};

}
}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of AbstractVVVVertex. */
template <>
struct BaseClassTrait<Helicity::AbstractVVVVertex,1> {
  /** Typedef of the first base class of AbstractVVVVertex. */
  typedef Helicity::VertexBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the AbstractVVVVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Helicity::AbstractVVVVertex>
  : public ClassTraitsBase<Helicity::AbstractVVVVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Helicity::AbstractVVVVertex"; }
};

/** @endcond */

}

#endif /* HELICITY_AbstractVVVVertex_H */
