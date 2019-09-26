// -*- C++ -*-
//
// HelicityVertex.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2017 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_HelicityVertex_H
#define ThePEG_HelicityVertex_H
// This is the declaration of the HelicityVertex class.

#include "HelicityVertex.fh"
#include "ThePEG/EventRecord/EventConfig.h"
#include "ThePEG/Utilities/ClassDescription.h"
#include "RhoDMatrix.h"
#include "ThePEG/Utilities/Rebinder.h"
#include "ThePEG/Helicity/HelicityDefinitions.h"

namespace ThePEG {

/**
 *  The HelicityVertex class is designed to store the helicity
 *  amplitude expression for the matrix element for use by the spin
 *  correlation algorithm. It implements the storage of the pointers
 *  to the incoming and outgoing particles at the vertex and virtual
 *  methods for calculating the rho and D matrices. The concrete
 *  implementations of the vertices for specific processes, eg
 *  production or decay, inherit from this and implement the storage
 *  of the matrix element together with the set and get methods.
 *
 *  These methods are then called by the SpinInfo class to perform the
 *  calculations.
 *
 *
 * @see SpinInfo
 *
 * @author Peter Richardson
 *
 */
class HelicityVertex: public EventInfoBase {

public:

  /**
   * Output the spin density matrix for debugging purposes.
   */
  friend ostream & operator<<(ostream & os, const HelicityVertex & vert);

public:

  /** A vector of SpinInfo objects. */
  typedef vector<tcSpinPtr> SpinVector;

public:

  /**
   * Standard Init function.
   */
  static void Init();

  /**
   * Rebind to cloned objects. If a HelicityVertex is cloned together
   * with a whole Event and this has pointers to other event record
   * objects, these should be rebound to their clones in this
   * function.
   */
  virtual void rebind(const EventTranslationMap & trans);

public:

  /** @name Access the incoming and outgoing particles. */
  //@{
  /**
   * Access the spin of the incoming particles.
   */
  const SpinVector & incoming() const {return _incoming;}

  /**
   * Access the spin of the outgoing particles.
   */
  const SpinVector & outgoing() const {return _outgoing;}

  /**
   * Add the spin of an incoming particle.
   * @param spin the spin of the particle.
   * @param loc is set to the position in the list of incoming spins.
   */
  void addIncoming(tcSpinPtr spin, int & loc) {
    if(loc<0) {
      _incoming.push_back(spin);
      loc=_incoming.size()-1;
    }
    else {
      _incoming[loc] = spin;
    }
  }

  /**
   * Add the spin of an outgoing particle.
   * @param spin the spin of the particle.
   * @param loc is set to the position in the list of outgoing spins.
   */
  void addOutgoing(tcSpinPtr spin, int & loc) {
    if(loc<0) {
      _outgoing.push_back(spin);
      loc=_outgoing.size()-1;
    }
    else {
      _outgoing[loc]= spin;
    }
  }

  /**
   * Reset the \a spin of the incoming particle at position \a loc.
   */
  void resetIncoming(tcSpinPtr spin, int loc) {
    assert( loc < int(_incoming.size()) && loc >= 0 );
    _incoming[loc]=spin;
  }

  /**
   * Reset the \a spin of the outgoing particle at position \a loc.
   */
  void resetOutgoing(tcSpinPtr spin, int loc) {
    assert( loc < int(_outgoing.size()) && loc >= 0 );
    _outgoing[loc]=spin;
  }
  //@}

public:

  /** @name Mthods to calculate rho and D matrices. */
  //@{
  /**
   * Get the rho matrix for the outgoing particle at position \a loc.
   */
  virtual RhoDMatrix getRhoMatrix(int loc,bool recursive) const = 0;

  /**
   * Get the D matrix for the incoming particle at position \a loc.
   */
  virtual RhoDMatrix getDMatrix(int loc) const = 0;
  //@}

private:

  /**
   * Describe an abstract base class without persistent data.
   */
  static AbstractNoPIOClassDescription<HelicityVertex> initHelicityVertex;

  /**
   * Private and non-existent assignment operator.
   */
  HelicityVertex & operator=(const HelicityVertex &) = delete;

private:

  /**
   * Pointers to the incoming particle spins at the vertex.
   */
  SpinVector _incoming;

  /**
   * Pointers to the outgoing particle spins at the vertex.
   */
  SpinVector _outgoing;

};

/**
 *  Output operator
 */
inline ostream & operator<<(ostream & os, const HelicityVertex & vert) {
  os << "the incoming particles at the vertex are" << endl;
  for(unsigned int ix=0;ix<vert._incoming.size();++ix) {
    os << "the " << ix << " th incoming particle " << vert._incoming[ix] << "\n";
  }
  os << "the outgoing particles at the vertex are" << endl;
  for(unsigned int ix=0;ix<vert._outgoing.size();++ix) {
    os << "the " << ix << " th outgoing particle " << vert._outgoing[ix] << "\n";
  }
  return os;
}
  
}


namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

  /**
   * This template specialization informs ThePEG about the
   * base class of HelicityVertex.
   */
  template <>
  struct BaseClassTrait<ThePEG::HelicityVertex,1>
    : public ClassTraitsType {
  /** Typedef of the base class of HelicityVertex. */
    typedef EventInfoBase NthBase;
  };

  /**
   * This template specialization informs ThePEG about the name of
   * the HelicityVertex class and the shared object where it is defined.
   */
template <>
struct ClassTraits<ThePEG::HelicityVertex>
  : public ClassTraitsBase<ThePEG::HelicityVertex> {
  /**
   * Return the class name.
   */
  static string className() { return "ThePEG::HelicityVertex"; }
};

/** @endcond */

}

#endif /* ThePEG_HelicityVertex_H */
