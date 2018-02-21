// -*- C++ -*-
//
// VertexBase.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2017 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_VertexBase_H
#define ThePEG_VertexBase_H
//
// This is the declaration of the VertexBase class.

#include <ThePEG/Interface/Interfaced.h>
#include <ThePEG/PDT/ParticleData.h>
#include <ThePEG/Helicity/HelicityDefinitions.h>
#include <ThePEG/Repository/EventGenerator.h>
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "VertexBase.fh"
#include <array>


namespace ThePEG {
namespace Helicity {

/**
 * Namespace for naming of vertices. Each child class should extend this
 * with its own spin configuration.
 */
namespace VertexType {
  typedef unsigned T;
  /**
   *  Undefined Enum for the Lorentz structures
   */
  const T UNDEFINED = 0;
}

/** \ingroup Helicity
 * 
 *  The VertexBase class is the base class for all helicity amplitude
 *  vertices. In implements the storage of the particles 
 *  which are allowed to interact at the vertex and some simple functions 
 *  which are often needed by the classes which implement the specific 
 *  vertices.
 *
 *  In practice little use is made of this information and it is mainly
 *  included for future extensions. It can also be used at the development
 *  and debugging stage.
 *
 */
class VertexBase  : public Interfaced {
/**
 *  The output operator is a friend to avoid the data being public.
 */
friend ostream & operator<<(ostream &, const VertexBase &);

public:
  
  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Constructor for \f$n\f$-point vertices.
   * @param name The type of vertex
   * @param kine Whether the kinematic invariants should be calculated.
   */
  VertexBase(VertexType::T name, bool kine=false);
  //@}

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
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();
  
public:

  /**
   *  Access to the particle information
   */
  //@{
  /**
   * Number of different particle combinations allowed.
   */
  unsigned int size() const { return _particles.size(); }

public:
  /**
   * Is a particle allowed as an incoming particle?
   * @param p The ParticleData pointer
   */
  bool isIncoming(tPDPtr p) const {
    return _inpart.find(p) != _inpart.end();
  }

  /**
   * Is a particle allowed as an outgoing particle?
   * @param p The ParticleData pointer
   */
  bool isOutgoing(tPDPtr p) const {
    return _outpart.find(p) != _outpart.end();
  }

  /**
   * Get the list of incoming particles.
   */
  const set<tPDPtr> & incoming() const { return _inpart; }

  /**
   * Get the list of outgoing particles.
   */
  const set<tPDPtr> & outgoing() const { return _outpart; }

  /**
   * Get the coupling.
   */
  Complex norm() const { return _norm; }

  /**
   * Function to search the list.
   * @param ilist Which list to search
   * @param id The PDG code to look for.
   */
  vector<long> search(unsigned int ilist,long id) const;

  /**
   * Function to search the list.
   * @param ilist Which list to search
   * @param id The particle to look for.
   */
  vector<tPDPtr> search(unsigned int ilist,tcPDPtr id) const;

  /**
   * Is a given combination allowed.
   * @param id1 PDG code of the first particle.
   * @param id2 PDG code of the second particle.
   * @param id3 PDG code of the third particle.
   * @param id4 PDG code of the fourth particle.
   */
  bool allowed(long id1, long id2, long id3, long id4 = 0) const;

  /**
   * Get name of Vertex
   */
  VertexType::T getName() const { return _theName; }

  /**
   * Get number of lines on Vertex
   */
  unsigned int getNpoint() const { return _npoint; }

  /**
   * Get the order in \f$g_EM\f$
   */
  unsigned int orderInGem() const { return _ordergEM; }

  /**
   * Get the order in \f$g_s\f$
   */
  unsigned int orderInGs() const { return _ordergS; }
  //@}

public:

  /**
   * @name Calculation of the strong, electromagnetic and weak couplings
   */
  //@{
  /**
   *  Strong coupling
   */
  double strongCoupling(Energy2 q2) const {
    if(_coupopt==0) {
      double val = 4.0*Constants::pi*generator()->standardModel()->alphaS(q2);
      assert(val>=0.);
      return sqrt(val);
    }
    else if(_coupopt==1)
      return sqrt(4.0*Constants::pi*generator()->standardModel()->alphaS());
    else
      return _gs;
  }

  /**
   *  Electromagnetic coupling
   */
  double electroMagneticCoupling(Energy2 q2) const {
    if(_coupopt==0)
      return sqrt(4.0*Constants::pi*generator()->standardModel()->alphaEMME(q2));
    else if(_coupopt==1)
      return sqrt(4.0*Constants::pi*generator()->standardModel()->alphaEMMZ());
    else
      return _ee;
  }

  /**
   *  Weak coupling
   */
  double weakCoupling(Energy2 q2) const {
    if( _coupopt == 0 )
      return sqrt(4.0*Constants::pi*generator()->standardModel()->alphaEMME(q2)/
		  generator()->standardModel()->sin2ThetaW());
    else if( _coupopt == 1 )
      return sqrt(4.0*Constants::pi*generator()->standardModel()->alphaEMMZ()/
		  generator()->standardModel()->sin2ThetaW());
    else
      return _ee/_sw;
  }

  double sin2ThetaW() const {
    if( _coupopt == 0 || _coupopt  == 1)
      return generator()->standardModel()->sin2ThetaW();
    else
      return sqr(_sw);
  }
  //@}

public:

  /**
   *   Set coupling methods
   */
  //@{
  /**
   * Calculate the couplings for a three point interaction.
   * This method is virtual and must be implemented in 
   * classes inheriting from this.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param part1 The ParticleData pointer for the first  particle.
   * @param part2 The ParticleData pointer for the second particle.
   * @param part3 The ParticleData pointer for the third  particle.
   */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1,
			   tcPDPtr part2,tcPDPtr part3)=0;

  /**
   * Calculate the couplings for a four point interaction.
   * This method is virtual and must be implemented in 
   * classes inheriting from this.
   * @param q2 The scale \f$q^2\f$ for the coupling at the vertex.
   * @param part1 The ParticleData pointer for the first  particle.
   * @param part2 The ParticleData pointer for the second particle.
   * @param part3 The ParticleData pointer for the third  particle.
   * @param part4 The ParticleData pointer for the fourth particle.
   */
  virtual void setCoupling(Energy2 q2,tcPDPtr part1,tcPDPtr part2,tcPDPtr part3,
			   tcPDPtr part4)=0;
  //@}

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

  /**
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given
   * pointer.
   */
  virtual void rebind(const TranslationMap & trans);

  /**
   * Return a vector of all pointers to Interfaced objects used in this
   * object.
   * @return a vector of pointers.
   */
  virtual IVector getReferences();
  //@}

protected:
  /**
   *  Members to set-up the particles
   */
  //@{
  /**
   * Set up the lists of outer particles for the vertex.
   * @param ids A vector of PDG codes for the particles.
   */
  void addToList(const vector<long> & ids);

  /**
   * Set up the lists of outer particles for the three-/four-point vertex.
   * For small vertices, this form is much easier to use.
   * @param ida The PDG codes for the first  set of particles.
   * @param idb The PDG codes for the second set of particles.
   * @param idc The PDG codes for the third  set of particles.
   * @param idd The PDG codes for the fourth set of particles.
   */
  void addToList(long ida, long idb, long idc, long idd = 0);
  //@}

protected:
  /**
   *  Members for the amplitude calculations
   */
  //@{
  /**
   * Set the coupling.
   * @param coup The coupling.
   */
  void norm(const Complex & coup) { _norm = coup; }

  /**
   * Calculate the propagator for a diagram.
   * @param iopt The option for the Breit-Wigner shape
   * @param q2 The scale
   * @param part The ParticleData pointer for the off-shell particle.
   * @param mass The mass if not to be taken from the ParticleData object
   * @param width The width if not to be taken from the ParticleData object
   */
  virtual Complex propagator(int iopt, Energy2 q2,tcPDPtr part,
			     complex<Energy> mass=-GeV,
			     complex<Energy> width=-GeV);

  /**
   * Calculate propagator multiplied by coupling.
   * @param iopt The option for the Breit-Wigner shape
   * @param q2 The scale
   * @param part The ParticleData pointer for the off-shell particle.
   * @param mass The mass if not to be taken from the ParticleData object
   * @param width The width if not to be taken from the ParticleData object
   */
  Complex normPropagator(int iopt, Energy2 q2,tcPDPtr part,
			 complex<Energy> mass=-GeV, 
			 complex<Energy> width=-GeV) {
    return _norm*propagator(iopt,q2,part,mass,width);
  }
  //@}    

public:
  /** @name Kinematic invariants for loop diagrams */
  //@{

  /**
   * Whether or not to calculate the kinematics invariants
   */
  bool kinematics() const { return _calckinematics; }

  /**
   * Set whether or not to calculate the kinematics invariants
   */
  void kinematics(bool kine ) { _calckinematics=kine; }

  /**
   *  Calculate the kinematics for a 3-point vertex
   */
  void calculateKinematics(const Lorentz5Momentum & p0,
			   const Lorentz5Momentum & p1,
			   const Lorentz5Momentum & p2) {
    _kine[0][0]=p0*p0;
    _kine[1][1]=p1*p1;
    _kine[2][2]=p2*p2;
    _kine[0][1]=p0*p1;_kine[1][0]=_kine[0][1];
    _kine[0][2]=p0*p2;_kine[2][0]=_kine[0][2];
    _kine[1][2]=p1*p2;_kine[2][1]=_kine[1][2];
  }
  
  /**
   *  Calculate the kinematics for a 4-point vertex
   */
  void calculateKinematics(const Lorentz5Momentum & p0,
			   const Lorentz5Momentum & p1,
			   const Lorentz5Momentum & p2,
			   const Lorentz5Momentum & p3) {
    _kine[0][0]=p0*p0;
    _kine[1][1]=p1*p1;
    _kine[2][2]=p2*p2;
    _kine[3][3]=p3*p3;
    _kine[0][1]=p0*p1;_kine[1][0]=_kine[0][1];
    _kine[0][2]=p0*p2;_kine[2][0]=_kine[0][2];
    _kine[0][3]=p0*p3;_kine[3][0]=_kine[0][3];
    _kine[1][2]=p1*p2;_kine[2][1]=_kine[1][2];
    _kine[1][3]=p1*p3;_kine[3][1]=_kine[1][3];
    _kine[2][3]=p2*p3;_kine[3][2]=_kine[2][3];
  }
  
  /**
   *  Calculate the kinematics for a n-point vertex
   */
  void calculateKinematics(const vector<Lorentz5Momentum> & p) {
    for(size_t ix=0;ix<p.size();++ix) {
      for(size_t iy=0;iy<=ix;++ix) {
	_kine[ix][iy]=p[ix]*p[iy];
	_kine[iy][ix]=_kine[ix][iy];
      }
    }
  }

  /**
   * Get one of the kinematic invariants
   */
  Energy2 invariant(unsigned int ix ,unsigned int iy) const {
    assert ( ix < _npoint && iy < _npoint );
    return _kine[ix][iy];
  }
  //@}
  
protected:

  /**
   * Set the order in \f$g_EM\f$
   * @param order The order of the vertex in \f$g_EM\f$
   */
  void orderInGem(unsigned int order) { _ordergEM = order; }

  /**
   * Set the order in \f$g_s\f$
   * @param order The order of the vertex in \f$g_s\f$
   */
  void orderInGs(unsigned int order) { _ordergS = order; }
  
private:
  
  /**
   * Describe a concrete class with persistent data.
   */
  static AbstractClassDescription<ThePEG::Helicity::VertexBase> initVertexBase;
  
  /**
   * Private and non-existent assignment operator.
   */
  VertexBase & operator=(const VertexBase &);
  
private:

  /**
   * Storage of the particles.   
   */
  //@{
  /**
   *  Particles interacting at the vertex
   */
  vector<vector<PDPtr> > _particles;

  /**
   *  Number of particles at the vertex
   */
  unsigned int _npoint;

  /**
   * ParticleData pointers for the allowed incoming particles.
   */
  set<tPDPtr> _inpart;

  /**
   * ParticleData pointers for the allowed outgoing particles.
   */
  set<tPDPtr> _outpart;
  //@}

  /**
   * The overall coupling.
   */
  Complex _norm;

  /**
   * Whether or not to calculate the kinematic invariants for the vertex
   */
  bool _calckinematics;

  /**
   * Kinematica quantities needed for loop vertices
   */
  std::array<std::array<Energy2,5>,5> _kine;

  /**
   * Name of vertex
   */
  VertexType::T _theName;

  /**
   * Order of vertex in \f$g_EM\f$
   */
  unsigned int _ordergEM;

  /**
   * Order of vertex in \f$g_s\f$
   */
  unsigned int _ordergS;

  /**
   *  option for the coupling
   */
  unsigned int _coupopt;

  /**
   *  Fixed value of strong coupling to use
   */
  double _gs;

  /**
   *  Fixed value of the electromagentic coupling to use
   */
  double _ee;

  /**
   *  Fixed value of \f$\sin\theta_W\f$ to use
   */
  double _sw;
};
  
/**
 * Output the information on the vertex.
 */
ostream & operator<<(ostream &, const VertexBase &);
  
}
}

#endif /* ThePEG_VertexBase_H */
