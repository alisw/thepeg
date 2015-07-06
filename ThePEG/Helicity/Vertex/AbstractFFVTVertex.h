// -*- C++ -*-
#ifndef HELICITY_AbstractFFVTVertex_H
#define HELICITY_AbstractFFVTVertex_H
//
// This is the declaration of the AbstractFFVTVertex class.
//

#include "VertexBase.h"
#include "ThePEG/Helicity/WaveFunction/VectorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/TensorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorWaveFunction.h"
#include "ThePEG/Helicity/WaveFunction/SpinorBarWaveFunction.h"
#include "AbstractFFVTVertex.fh"

namespace ThePEG {
namespace Helicity {

/**
 * The AbstractFFVTVertex class is the base class for all 
 * fermion-fermion-vector-tensor interactions in ThePEG.
 */
class AbstractFFVTVertex: public VertexBase {

public:
  

  /**
   * Default constructor
   */
  AbstractFFVTVertex() : VertexBase(VertexType::FFVT) {}

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
   * @param ten4  The wavefunction for the tensor.
   */
  virtual Complex evaluate(Energy2 q2,const SpinorWaveFunction & sp1,
			   const SpinorBarWaveFunction & sbar2,
			   const VectorWaveFunction & vec3,
			   const TensorWaveFunction & ten4) = 0;

  /**
   * Evaluate the off-shell tensor coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell tensor.
   * @param out The ParticleData pointer for the off-shell tensor.
   * @param sp1   The wavefunction for the ferimon.
   * @param sbar2 The wavefunction for the antifermion.
   * @param vec3  The wavefunction for the vector.
   * @param mass The mass of the off-shell particle if not taken from the ParticleData
   * object
   * @param width The width of the off-shell particle if not taken from the ParticleData
   * object
   */
  virtual TensorWaveFunction evaluate(Energy2 q2,int iopt, tcPDPtr out,
				      const SpinorWaveFunction & sp1,
				      const SpinorBarWaveFunction & sbar2,
				      const VectorWaveFunction & vec3,
				      complex<Energy> mass=-GeV, complex<Energy> width=-GeV) = 0;

  /**
   * Evaluate the off-shell vector coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell vector.
   * @param out The ParticleData pointer for the off-shell vector.
   * @param sp1   The wavefunction for the ferimon.
   * @param sbar2 The wavefunction for the antifermion.
   * @param ten4  The wavefunction for the tensor.
   * @param mass The mass of the off-shell particle if not taken from the ParticleData
   * object
   * @param width The width of the off-shell particle if not taken from the ParticleData
   * object
   */
  virtual VectorWaveFunction evaluate(Energy2 q2,int iopt, tcPDPtr out,
				      const SpinorWaveFunction & sp1,
				      const SpinorBarWaveFunction & sbar2, 
				      const TensorWaveFunction & ten4,
				      complex<Energy> mass=-GeV, complex<Energy> width=-GeV) = 0;

  /**
   * Evaluate the off-shell spinor coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell spinor.
   * @param out The ParticleData pointer for the off-shell spinor.
   * @param sp1   The wavefunction for the ferimon.
   * @param vec3  The wavefunction for the vector.
   * @param ten4  The wavefunction for the tensor.
   * @param mass The mass of the off-shell particle if not taken from the ParticleData
   * object
   * @param width The width of the off-shell particle if not taken from the ParticleData
   * object
   */
  virtual SpinorWaveFunction evaluate(Energy2 q2,int iopt, tcPDPtr out,
				      const SpinorWaveFunction & sp1,
				      const VectorWaveFunction & vec3,
				      const TensorWaveFunction & ten4,
				      complex<Energy> mass=-GeV, complex<Energy> width=-GeV) = 0;

  /**
   * Evaluate the off-shell barred spinor coming from the vertex.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param iopt Option of the shape of the Breit-Wigner for the off-shell barred spinor.
   * @param out The ParticleData pointer for the off-shell barred spinor.
   * @param sbar2 The wavefunction for the antifermion.
   * @param vec3  The wavefunction for the vector.
   * @param ten4  The wavefunction for the tensor.
   * @param mass The mass of the off-shell particle if not taken from the ParticleData
   * object
   * @param width The width of the off-shell particle if not taken from the ParticleData
   * object
   */
  virtual SpinorBarWaveFunction evaluate(Energy2 q2,int iopt, tcPDPtr out,
					 const SpinorBarWaveFunction & sbar2,
					 const VectorWaveFunction & vec3,
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
  static AbstractNoPIOClassDescription<AbstractFFVTVertex> initAbstractFFVTVertex;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  AbstractFFVTVertex & operator=(const AbstractFFVTVertex &);

};

}
}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of AbstractFFVTVertex. */
template <>
struct BaseClassTrait<Helicity::AbstractFFVTVertex,1> {
  /** Typedef of the first base class of AbstractFFVTVertex. */
  typedef Helicity::VertexBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the AbstractFFVTVertex class and the shared object where it is defined. */
template <>
struct ClassTraits<Helicity::AbstractFFVTVertex>
  : public ClassTraitsBase<Helicity::AbstractFFVTVertex> {
  /** Return a platform-independent class name */
  static string className() { return "Helicity::AbstractFFVTVertex"; }
};

/** @endcond */

}

#endif /* HELICITY_AbstractFFVTVertex_H */
