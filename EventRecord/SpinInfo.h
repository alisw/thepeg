// -*- C++ -*-
//
// SpinInfo.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2017 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_SpinInfo_H
#define ThePEG_SpinInfo_H
// This is the declaration of the SpinInfo class.

#include "ThePEG/EventRecord/EventInfoBase.h"
#include "ThePEG/PDT/PDT.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "HelicityVertex.h"

namespace ThePEG {

/**
 *   The SpinInfo is the base class for the spin information for the
 *   spin correlation algorithm. The implementations for different spin
 *   states inherit from this.
 *
 *   The class contains pointers to the vertex where the particle is
 *   produced and where it decays, together with methods to set/get
 *   these.
 *
 *   There are two flags decayed which store information on the state
 *   of the particle.
 *
 *   The decayed() members provides access to the _decay data member
 *   which is true if the spin density matrix required to perform the
 *   decay of a timelike particle has been calculated (this would be a
 *   decay matrix for a spacelike particle.) This is set by the
 *   decay() method which calls a method from the production vertex to
 *   calculate this matrix. The decay() method should be called by a
 *   decayer which uses spin correlation method before it uses the
 *   spin density matrix to calculate the matrix element for the
 *   decay.
 *
 *   The developed() member provides access to the _developed data
 *   member which is true if the decay matrix required to perform the
 *   decays of the siblings of a particle has been calculated (this
 *   would a spin density matrix for a spacelike particle.) This is
 *   set by the developed() method which calls a method from the decay
 *   vertex to calculate the matrix. The developed() method is called
 *   by a DecayHandler which is capable of performing spin
 *   correlations after all the unstable particles produced by a
 *   decaying particle are decayed.
 *
 *   Methods are also provided to access the spin density and decay
 *   matrices for a particle.
 *
 * @author Peter Richardson
 *
 */
class SpinInfo: public EventInfoBase {

public:

  /**
   *  Status for the implementation of spin correlations
   */
  enum DevelopedStatus {
    Undeveloped=0, /**< Not developed. */
    Developed=1,   /**< Developed. */
    NeedsUpdate=2  /**< Developed but needs recalculating due to some change. */
  };

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  SpinInfo() 
    : _timelike(false), _prodloc(-1), _decayloc(-1), 
      _decayed(false), _developed(Undeveloped) {}

  /**
   * Standard Constructor.
   * @param s the spin.
   * @param p the production momentum.
   * @param time true if the particle is time-like.
   */
  SpinInfo(PDT::Spin s, 
	   const Lorentz5Momentum & p = Lorentz5Momentum(), 
	   bool time = false)
    : _timelike(time), _prodloc(-1), _decayloc(-1),
      _decayed(false), _developed(Undeveloped),
      _rhomatrix(s), _Dmatrix(s), _spin(s),
      _productionmomentum(p), _currentmomentum(p) {}

  /**
   * Copy-constructor.
   */
  SpinInfo(const SpinInfo &);
  //@}

public:

  /**
   * Returns true if the polarization() has been implemented in a
   * subclass. This default version returns false.
   */
  virtual bool hasPolarization() const { return false; }

  /**
   * Return the angles of the polarization vector as a pair of
   * doubles. first is the polar angle and second is the azimuth
   * wrt. the particles direction. This default version of the
   * function returns 0,0, and if a subclass implements a proper
   * function it should also implement 'hasPolarization()' to return
   * true.
   */
  virtual DPair polarization() const { return DPair(); }

public:

  /**
   * Standard Init function.
   */
  static void Init();

  /**
   * Rebind to cloned objects. If a FermionSpinInfo is cloned together
   * with a whole Event and this has pointers to other event record
   * objects, these should be rebound to their clones in this
   * function.
   */
  virtual void rebind(const EventTranslationMap & trans);

  /**
   * Standard clone method.
   */
  virtual EIPtr clone() const;

  /**
   * Method to handle the delelation
   */
  void update() const;

  /**
   * Perform a lorentz rotation of the spin information
   */
  virtual void transform(const LorentzMomentum & m, const LorentzRotation & r) {
    _currentmomentum = m;
    _currentmomentum.transform(r);
  }

public:


  /** @name Access the vertices. */
  //@{
  /**
   * Set the vertex at which the particle was produced.
   */
  void productionVertex(VertexPtr in) const {
    _production=in;
    // add to the list of outgoing if timelike
    int temp(-1);
    if(_timelike) in->addOutgoing(this,temp); 
    // or incoming if spacelike
    else          in->addIncoming(this,temp);
    _prodloc=temp;
  }

  /**
   * Get the vertex at which the particle was produced.
   */
  tcVertexPtr productionVertex() const { return _production; }

  /**
   * Set the vertex at which the particle decayed or branched.
   */
  void decayVertex(VertexPtr in) const {
    if(in) {
      _decay=in;
      if(_timelike) {
	int temp(-1);
	in->addIncoming(this,temp);
	_decayloc=temp;
	assert(temp==0);
      }
      else {
	int temp(-1);
	in->addOutgoing(this,temp);
	_decayloc=temp;
      }
    }
    else {
      _decay=VertexPtr();
      _decayloc=-1;
    }
  }

  /**
   * Get the vertex at which the particle decayed or branched.
   */
  tcVertexPtr decayVertex() const { return _decay; }
  //@}

  /** @name Access information about the associated particle. */
  //@{
  /**
   * Has the particle decayed?
   */
  bool decayed() const { return _decayed; }

  /**
   * Set if the particle has decayed.
   */
  void decayed(bool b) const { _decayed = b; }

  /**
   * Return true if the decay matrix required to perform the decays of
   * the siblings of a particle has been calculated.
   */
  DevelopedStatus developed() const { return _developed; }

  /**
   * Calculate the rho matrix for the decay if not already done.
   */
  void decay(bool recursive=false) const ;

  /**
   * Set the developed flag and calculate the D matrix for the decay.
   */
  void develop() const ;

  /**
   *  Needs update
   */
  void needsUpdate() const {_developed=NeedsUpdate;}

  /**
   * Return 2s+1 for the particle
   */
  PDT::Spin iSpin() const { return _spin; }

  /**
   * Return the momentum of the particle when it was produced.
   */
  const Lorentz5Momentum & productionMomentum() const {
    return _productionmomentum;
  }

  /**
   *  The current momentum of the particle
   */
  const Lorentz5Momentum & currentMomentum() const {
    return _currentmomentum;
  }

  /**
   * Return true if particle is timelike (rather than spacelike).
   */
  bool timelike() const { return _timelike; }
  //@}

  /**
   *  Access to the locations
   */
  //@{
  /**
   *  Production Location
   */
  int productionLocation() const {return _prodloc;}

  /**
   *  Decay Location
   */
  int decayLocation() const {return _decayloc;}
  //@}

public:

  /** @name Access the rho and D matrices. */
  //@{
  /**
   * Access the rho matrix.
   */
  RhoDMatrix rhoMatrix() const { return _rhomatrix; }

  /**
   * Access the rho matrix.
   */
  RhoDMatrix & rhoMatrix() { return _rhomatrix; }

  /**
   * Access the D matrix.
   */
  RhoDMatrix DMatrix() const { return _Dmatrix; }

  /**
   * Access the D matrix.
   */
  RhoDMatrix & DMatrix() { return _Dmatrix; }
  //@}

protected:

  /**
   *  Check if momentum is near to the current momentum
   */
  bool isNear(const Lorentz5Momentum & p) {
    return currentMomentum().isNear(p,_eps);
  }

private:

  /**
   * Describe a concrete class without persistent data.
   */
  static NoPIOClassDescription<SpinInfo> initSpinInfo;

  /**
   * Private and non-existent assignment operator.
   */
  SpinInfo & operator=(const SpinInfo &);

private:

  /**
   * Set the developed flag and calculate the D matrix for the decay,
   * and all decays further up the chain.
   */
  void redevelop() const ;

  /**
   * Recursively recalulate all the rho matrices from the top of the chain
   */
  void redecay() const ;

private:

  /**
   * Pointer to the production vertex for the particle
   */
  mutable VertexPtr _production;

  /**
   * Pointers to the decay vertex for the particle
   */
  mutable VertexPtr _decay;

  /**
   * Is this is timelike (true) or spacelike (false ) particle?  This
   * is used to decide if the particle is incoming or outgoing at the
   * production vertex
   */
  bool _timelike;

  /**
   * Location in the hard vertex array at production.
   */
  mutable int _prodloc;

  /**
   * Location in the hard vertex array at decay.
   */
  mutable int _decayloc;

  /**
   * Has the particle been decayed?  (I.e. has the rho matrix for the
   * decay been calculated.)
   */
  mutable bool _decayed;

  /**
   * Has the particle been developed?  (I.e. has the D matrix encoding
   * the info about the decay been calculated)
   */
  mutable DevelopedStatus _developed;

  /**
   * Storage of the rho matrix.
   */
  mutable RhoDMatrix _rhomatrix;

  /**
   * Storage of the decay matrix
   */
  mutable RhoDMatrix _Dmatrix;

  /**
   * The spin of the particle
   */
  PDT::Spin _spin;

  /**
   * Momentum of the particle when it was produced
   */
  Lorentz5Momentum _productionmomentum;

  /**
   * Momentum of the particle when it decayed
   */
  mutable Lorentz5Momentum _decaymomentum;
  
  /**
   * Current momentum of the particle
   */
  Lorentz5Momentum _currentmomentum;

  /**
   *  A small energy for comparing momenta to check if Lorentz Transformations
   *  should be performed
   */
  static const double _eps;
};

}


namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * This template specialization informs ThePEG about the base class of
 * SpinInfo.
 */
template <>
struct BaseClassTrait<ThePEG::SpinInfo,1>: public ClassTraitsType {
  /** Typedef of the base class of SpinInfo. */
  typedef EventInfoBase NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * SpinInfo class and the shared object where it is defined.
 */
template <>
struct ClassTraits<ThePEG::SpinInfo>
  : public ClassTraitsBase<ThePEG::SpinInfo> {
  /**
   * Return the class name.
   */
  static string className() { return "ThePEG::SpinInfo"; }
};

/** @endcond */

}

#endif /* ThePEG_SpinInfo_H */
