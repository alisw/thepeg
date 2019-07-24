// -*- C++ -*-
#ifndef THEPEG_MixedParticleData_H
#define THEPEG_MixedParticleData_H
//
// This is the declaration of the MixedParticleData class.
//

#include "ParticleData.h"
#include "MixedParticleData.fh"

namespace ThePEG {

/**
 * The MixedParticleData class is designed to store the 
 * particle data for particles which undergo mixing with their
 * anti-particle, e.g. \f$B^0-\bar{B}^0\f$, \f$B_s^0-\bar{B}_s^0\f$, 
 * \f$D^0-\bar{D}^0\f$.
 *
 *  The key parameters are the mass difference between the heavier 
 *  and lighter states, \f$\Delta m = m_H-m_L\f$, the width difference
 *  \f$\Delta\Gamma = \Gamma_H-\Gamma_L\f$. While the mass difference
 *  is positive by definition the sign of the width difference must be
 *  determined experimentally. It should be noted that because the width
 *  difference is expected to be negative in the Standard Model many
 *  experiment results use the convention 
 *  \f$\Delta\Gamma = \Gamma_L - \Gamma_H\f$ so care should be taken.
 *
 *  In order to be as general as possible we have included the option
 *  of CPT violation, in addition to the violation of CP, in the mixing
 *  in this case the light and heavy eigenstates are given by
 * \f[|M_L\rangle \propto p \sqrt{1-z}|     M^0 \rangle 
 *                    + q \sqrt{1+z}|\bar{M}^0\rangle,\f]
 * \f[|M_H\rangle \propto p \sqrt{1+z}|     M^0 \rangle 
 *                    - q \sqrt{1-z}|\bar{M}^0\rangle,\f]
 *  where \f$p\f$, \f$q\f$ and \f$z\f$ are complex parameters with \f$z=0\f$ if CPT
 *  is conserved.
 *
 *  This gives the time evolution of the states
 * \f[|M^0_{\rm phys}(t)\rangle = (g_+(t)+zg_-(t))|M^0\rangle
 *                                -sqrt{1-z^2}\frac{q}{p}g_-(t)|\bar{M}^0\rangle,\f]
 * \f[\bar{M}^-_{\rm phys}\rangle = (g_+(t)+zg_-(t))|\bar{M}^0\rangle
 *                                -sqrt{1-z^2}\frac{p}{q}g_-(t)|M^0\rangle,\f]
 * for particles initially in a pure particle or antiparticle state
 * where
 * \f[g_\pm(t)=\frac12\left(e^{-im_Ht-\frac12\Gamma_Ht} \pm
 *                       e^{-im_Lt-\frac12\Gamma_Lt}\right).\f]
 * @see \ref MixedParticleDataInterfaces "The interfaces"
 * defined for MixedParticleData.
 */
class MixedParticleData: public ParticleData {

public:

  /**
   * The default constructor.
   */
  MixedParticleData() : ParticleData(), _deltam(0.*GeV), 
			_deltagamma(0.*GeV), _pqmag(1.), _pqphase(0.),
			_pq(1.,0.), _zmag(0.), _zphase(0.), _z(0.), _x(0.), _y(0.),
			_prob(make_pair(1.,0.))
  {}

  /** @name The Create methods are special interfaces for ParticleData
      classes. */
  //@{
  /**
   * Create a Particle which is its own anti-particle.
   */
  static PDPtr Create(long newId, string newPDGName);

  /**
   * Create a particle - anti particle pair.
   */
  static PDPair Create(long newId, string newPDGName, string newAntiPDGName);
  //@}

public:

  /**
   *  Mixing parameters
   */
  //@{
  /**
   *  The mass difference
   */
  Energy deltaM() const {return _deltam;}

  /**
   *  The width difference
   */
  Energy deltaGamma() const {return _deltagamma;}

  /**
   *  \f$p/q\f$
   */
  Complex pq() const {return _pq;}

  /**
   *  \f$z\f$
   */
  Complex z() const {return _z;}

  /**
   *  The \f$x\f$ mixing variable
   */
  double x() const {return _x;}

  /**
   *  the \f$y\f$ mixing variable
   */
  double y() const {return _y;}

  /**
   *  The time-integrated mixing probabilities
   */
  pair<double,double> prob() const {return _prob;}

  /**
   *  For a given paricle decide is it undergoes mixing
   *  and generate the lifetime
   */
  pair<bool,Length> generateLifeTime() const;

  /**
   *  The amplitudes for the different states
   */
  pair<Complex,Complex> mixingAmplitudes(Length,bool) const;

public:

  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void persistentInput(PersistentIStream & is, int version);
  //@}

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

protected:

  /**
   * Protected constructor only to be used by subclasses or by the
   * Create method.
   */
  MixedParticleData(long newId, string newPDGName);

  /**
   * ParticleData clone method
   */
  virtual PDPtr pdclone() const;

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}

private:

  /**
   *  Function for the interface to set the mass difference
   */
  void setDeltaM(Energy);

  /**
   *  Function for the interface to set the width difference
   */
  void setDeltaGamma(Energy);

  /**
   *  Function for the interface to set the magnitude of p/q
   */
  void setPQMagnitude(double);

  /**
   *  Function for the interface to set the phase of p/q
   */
  void setPQPhase(double);

  /**
   *  Function for the interface to set the magnitude of z
   */
  void setZMagnitude(double);

  /**
   *  Function for the interface to set the phase of z
   */
  void setZPhase(double);

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<MixedParticleData> initMixedParticleData;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MixedParticleData & operator=(const MixedParticleData &) = delete;

private:

  /**
   *  Mixing parameters
   */
  //@{
  /**
   *  The mass difference
   */
  Energy _deltam;

  /**
   *  The width difference
   */
  Energy _deltagamma;

  /**
   *  The magnitude of \f$p/q\f$
   */
  double _pqmag;

  /**
   *  The phase of \f$p/q\f$
   */
  double _pqphase;

  /**
   *  \f$p/q\f$
   */
  Complex _pq;

  /**
   *  The magnitude of \f$z\f$
   */
  double _zmag;

  /**
   *  The phase of \f$z\f$
   */
  double _zphase;

  /**
   * The \f$z\f$ parameter
   */
  Complex _z;

  /**
   *  The \f$x\f$ mixing variable
   */
  double _x;

  /**
   *  the \f$y\f$ mixing variable
   */
  double _y;

  /**
   *  The time-integrated mixing probabilities
   */
  pair<double,double> _prob;
  //@}

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MixedParticleData. */
template <>
struct BaseClassTrait<MixedParticleData,1> {
  /** Typedef of the first base class of MixedParticleData. */
  typedef ParticleData NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MixedParticleData class and the shared object where it is defined. */
template <>
struct ClassTraits<MixedParticleData>
  : public ClassTraitsBase<MixedParticleData> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::MixedParticleData"; }
};

/** @endcond */

}

#endif /* THEPEG_MixedParticleData_H */
