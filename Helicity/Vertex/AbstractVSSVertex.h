// -*- C++ -*-
#ifndef HELICITY_AbstractVSSVertex_H
#define HELICITY_AbstractVSSVertex_H
//
// This is the declaration of the AbstractVSSVertex class.
//

#include "VertexBase.h"
#include "ThePEG/Helicity/WaveFunction/ScalarWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "AbstractVSSVertex.fh"

namespace ThePEG {
namespace Helicity {

/**
 * The AbstractVSSVertex class is the base class for vector-scalar-scalar
 * interactions in ThePEG
 */
class AbstractVSSVertex: public VertexBase {

public:

  /**
   * Default constructor
   */
  AbstractVSSVertex() : VertexBase(VertexType::VSS) {}

  /**
   * Members to calculate the helicity amplitude expressions for vertices
   * and off-shell particles.
   */
  //@{  
  /**
   * Evaluate the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param vec1 The wavefunction for the vector.
   * @param sca2 The wavefunction for the first  scalar.
   * @param sca3 The wavefunction for the second scalar.
   */
  virtual Complex evaluate(Energy2 q2,const VectorWaveFunction & vec1,
			   const ScalarWaveFunction & sca2, 
			   const ScalarWaveFunction & sca3) = 0;

  /**
   * Evaluate the off-shell vector coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell vector.
   * @param out The ParticleData pointer for the off-shell vector.
   * @param sca2 The wavefunction for the first  scalar.
   * @param sca3 The wavefunction for the second scalar.
   * @param mass The mass of the off-shell particle if not taken from the ParticleData
   * object
   * @param width The width of the off-shell particle if not taken from the ParticleData
   * object
   */
  virtual VectorWaveFunction evaluate(Energy2 q2,int iopt,tcPDPtr out,
				      const ScalarWaveFunction & sca2,
				      const ScalarWaveFunction & sca3,
				      complex<Energy> mass=-GeV, complex<Energy> width=-GeV) = 0;

  /**
   * Evaluate the off-shell scalar coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell scalar.
   * @param out The ParticleData pointer for the off-shell scalar.
   * @param vec1 The wavefunction for the vector.
   * @param sca2 The wavefunction for the scalar.
   * @param mass The mass of the off-shell particle if not taken from the ParticleData
   * object
   * @param width The width of the off-shell particle if not taken from the ParticleData
   * object
   */
  virtual ScalarWaveFunction evaluate(Energy2 q2,int iopt, tcPDPtr out,
				      const VectorWaveFunction & vec1,
				      const ScalarWaveFunction & sca2,
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
  static AbstractNoPIOClassDescription<AbstractVSSVertex> initAbstractVSSVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  AbstractVSSVertex & operator=(const AbstractVSSVertex &);

};

}
}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of AbstractVSSVertex. */
template <>
struct BaseClassTrait<Helicity::AbstractVSSVertex,1> {
  /** Typedef of the first base class of AbstractVSSVertex. */
  typedef Helicity::VertexBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the AbstractVSSVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Helicity::AbstractVSSVertex>
  : public ClassTraitsBase<Helicity::AbstractVSSVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Helicity::AbstractVSSVertex"; }
};

/** @endcond */

}

#endif /* HELICITY_AbstractVSSVertex_H */
