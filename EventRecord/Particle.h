// -*- C++ -*-
//
// Particle.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_Particle_H
#define ThePEG_Particle_H
// This is the decalaration of the Particle class.

#include "EventConfig.h"
#include "ThePEG/Vectors/Lorentz5Vector.h"
#include "ThePEG/Vectors/LorentzRotation.h"
#include "ThePEG/Utilities/ClassDescription.h"
#include "ThePEG/EventRecord/MultiColour.h"
#include "ThePEG/EventRecord/SpinInfo.h"
#include "ThePEG/PDT/ParticleData.h"

namespace ThePEG {

/**
 * The Particle class is used to describe an instance of a
 * particle. Properties of the corresponding particle type can be
 * accessed through a pointer to a ParticleData object.
 *
 * A Particle object contains pointers to other particles, such as a
 * list of parents and a list of children. It may also contain a
 * pointer to the previous or next instance of the same physical
 * particle if the properties of a given particle has been changed
 * during the generation. Coloured particles contains pointers to
 * ColourLine defining the colour connections to other particles.
 *
 * The Particle also has a pointer to the Step object where it was
 * first introduced in the Event.
 *
 * When printing a particle the format of the output is governed by
 * the static <code>outputFormat</code> string. When a particle is sent
 * to an <code>ostream</code>, the format string is written but with
 * keys prefixed by the <code>\%</code> character replaced with
 * infromation about the particle as follows:<BR> <code>\%\%</code> is
 * replaced by a singel <code>\%</code><BR> <code>\%C</code> sets a flag so
 * that subsequent output of children and parents etc. will contain
 * colour information.<BR> <code>\%n</code> is replaced by the particles
 * number in a fied ow width<BR> <code>\%s</code> is replaced by the name
 * of the particle<BR> <code>\%i</code> is replaced by the id number of
 * the particle type<BR><code>\%x, \%y, \%z, \%e, \%m</code> is replaced by
 * the x-, y-, z-, energy- and mass-component of the particles
 * momentum respectively<BR><code>\%dx, \%dy, \%dz, \%dt, \%dT</code> is
 * replaced by the x-, y-, z-, time- and invariant time-component of
 * the particles lifeLength respectively<BR><code>\%Vx, \%Vy, \%Vz,
 * \%Vt</code> is replaced by the x-, y-, z- and time-component of the
 * creation point relative to the vertex of the
 * collision.<BR><code>\%Lx, \%Ly, \%Lz, \%Lt</code> is replaced by the x-,
 * y-, z- and time-component of the creation point in the current
 * lab-system.<BR> <code>\%p[,]</code> is replaced by a list (of numbers)
 * of the particles parents enclosed by <code>[</code> and <code>]</code>
 * and separated by <code>,</code>, The parent from which the particle
 * inherits its (anti-) colour is marked by a (-)+<BR>
 * <code>\%c(,)</code> is replaced by a list (of numbers) of the
 * particles children enclosed by <code>(</code> and <code>)</code> and
 * separated by <code>,</code>, The child which inherits the particles
 * (anti-) colour is marked by a (-)+<BR> <code>\%&gt;</code> is replaced
 * by the number of the colour neighbor<BR> <code>\%&lt;</code> is
 * replaced by the number of the anti-colour neighbor<BR>
 * <code>\%^</code> is replaced by the number of the previous instance of
 * the same physical particle<BR> <code>\%v</code> is replaced by the
 * number of the next instance of the same physical particle.<BR>
 * <code>\%l{,}</code> is replaced by the indices of the colour lines to
 * which this particle is connected enclosed by <code>{</code> and
 * <code>}</code> and separated by <code>,</code>, The line corresponding
 * to the (anti-) colour of the particle is prefixed by a (-)+
 *
 * @see Event
 * @see Collision
 * @see Step
 * @see SubProcess
 * @see Lorentz5Vector
 * @see ColourLine
 * @see ColourBase 
 */
class Particle: public EventRecordBase {

public:

  /** Most of the Event classes are friends with each other. */
  friend class Event;
  /** Most of the Event classes are friends with each other. */
  friend class Collision;
  /** Most of the Event classes are friends with each other. */
  friend class Step;
  /** Most of the Event classes are friends with each other. */
  friend class SubProcess;
  /** ParticleData needs to be a friend. */
  friend class ParticleData;

  struct ParticleRep;

public:

  /**
   * Standard Constructor. Note that the default constructor is
   * private - there is no particle without a pointer to a
   * ParticleData object.
   */
  Particle(tcEventPDPtr newData) : theData(newData), theRep(0) {}

  /**
   * Copy constructor.
   */
  Particle(const Particle &);

  /**
   * Destructor.
   */
  virtual ~Particle();
  //@}

  /** @name Functions relating to ancestry of particles. */
  //@{
  /**
   * Returns true if and only if this particle has decayed.
   */
  bool decayed() const {
    return hasRep() && !rep().theChildren.empty();
  }

  /**
   * The list of decay products.
   */
  const ParticleVector & children() const {
    static const ParticleVector null;
    return hasRep() ? rep().theChildren : null;
  }

  /**
   * Add a child (the childs parent pointer will be set accordingly).
   */
  void addChild(tPPtr c) {
      rep().theChildren.push_back(c);
      (c->rep()).theParents.push_back(this);
  }

  /**
   * Remove the given child from the list of children of this particle
   * (the corresponding parent pointer of the child will also be
   * removed).
   */
  void abandonChild(tPPtr child) {
    removeChild(child);
    child->removeParent(this);
  }

  /**
   * The list of parent particles.
   */
  const tParticleVector & parents() const {
    static const tParticleVector null;
    return hasRep() ? rep().theParents : null;
  }

  /**
   * Return a set of neighboring particles coming from the same decay
   * as this one. The return value is a newly recalculated set
   * every time. It must be stored to be used further, do not directly call
   * e.g. siblings().begin() or siblings().end()!
   */
  tParticleSet siblings() const;

  /**
   * Undo the decay of this particle, removing all children (and
   * grand children ...) from the event record
   */
  void undecay() {
    if ( hasRep() ) {
      rep().theChildren.clear();
      rep().theNext = tPPtr();
    }
  }

  /**
   * If this particle has decayed set the corresponding decay mode.
   */
  void decayMode(tDMPtr dm) { rep().theDecayMode = dm; }

  /**
   * If this particle has decayed get the corresponding decay mode.
   */
  tDMPtr decayMode() const {
    return hasRep() ? rep().theDecayMode : tDMPtr();
  }

  /**
   * Next instance. Pointer to another instance of the same
   * physical particle in later steps.
   */
  tPPtr next() const {
    return hasRep() ? rep().theNext : PPtr();
  }

  /**
   * Previous instance. Pointer to another instance of the same
   * physical particle in earlier steps.
   */
  tPPtr previous() const {
    return hasRep() ? rep().thePrevious : tPPtr();
  }

  /**
   * Original instance. If there exists another previous instance of
   * this particle return this instance (recursively).
   */
  tcPPtr original() const {
    return previous() ? tcPPtr(previous()->original()) : tcPPtr(this);
  }

  /**
   * Original instance. If there exists another previous instance of
   * this particle return this instance (recursively).
   */
  tPPtr original() {
    return previous() ? previous()->original() : tPPtr(this);
  }


  /**
   * Final instance. If there exists another subsequent instance of
   * this particle return this instance (recursively).
   */
  tcPPtr final() const {
    return next() ? tcPPtr(next()->final()) : tcPPtr(this);
  }

  /**
   * Final instance. If there exists another subsequent instance of
   * this particle return this instance (recursively).
   */
  tPPtr final() {
    return next() ? next()->final() : tPPtr(this);
  }

  //@}

  /** @name Relations to the Event and Step. */
  //@{
  /**
   * Get the first Step object where this particle occurred.
   */
  tStepPtr birthStep() const { 
    return hasRep() ? rep().theBirthStep : tStepPtr();
  }

  /**
   * Get the order-number for this particle in the current event.
   */
  int number() const { 
    return hasRep() ? rep().theNumber : 0; 
  }
  //@}

  /** @name Access the underlying ParticleData object. */
  //@{
  /**
   * Access the ParticleData object of this particle type
   */
  const ParticleDataClass & data() const { return *theData; }

  /**
   * Access the ParticleData object of this particle type
   */
  tcEventPDPtr dataPtr() const { return theData; }

  /**
   * Return the PDG name of this particle.
   */
  const string & PDGName() const { return data().PDGName(); }

  /**
   * Return the PDG id number of this particle.
   */
  long id() const { return data().id(); }
  //@}

  /** @name Functions to access the momentum. */
  //@{
  /**
   * Return the momentum of this particle.
   */
  const Lorentz5Momentum & momentum() const { return theMomentum; }

  /**
   * Set the 3-momentum of this particle. The energy is set to be
   * consistent with the mass.
   */
  void set3Momentum(const Momentum3 & p) {
    theMomentum.setVect(p);
    theMomentum.rescaleEnergy();
  }

  /**
   * Set the momentum of this particle. Afterwards, the underlying
   * Lorentz5Momentum may have inconsistent mass.
   */
  void setMomentum(const LorentzMomentum & p) {
    theMomentum = p;
  }

  /**
   * Set the momentum and mass.
   */
  void set5Momentum(const Lorentz5Momentum & p) {
    theMomentum = p;
  }

  /**
   * Acces the mass of this particle.
   */
  Energy mass() const { return momentum().mass(); }

  /**
   * Acces the mass of this particle type.
   */
  Energy nominalMass() const { return data().mass(); }

  /**
   * Get the scale at which this particle is considered resolved.
   */
  Energy2 scale() const { 
    return hasRep() ? rep().theScale : -1.0*GeV2;
  }

  /**
   * Set the scale at which this particle is considered resolved.
   */
  void scale(Energy2 q2) { rep().theScale = q2; }

  /**
   * Get the scale above which this particle should
   * not radiate.
   */
  Energy2 vetoScale() const { 
    return hasRep() ? rep().theVetoScale : -1.0*GeV2;
  }

  /**
   * Set the scale above which this particle should
   * not radiate.
   */
  void vetoScale(Energy2 q2) { rep().theVetoScale = q2; }

  /**
   * Return the transverse mass (squared), calculated from the energy
   * and the longitudinal momentum.
   */
  Energy2 mt2() const { return sqr(momentum().t()) - sqr(momentum().z()); }

  /**
   * Return the transverse mass (squared), calculated from the energy
   * and the longitudinal momentum.
   */
  Energy mt() const { return sqrt(mt2()); }

  /**
   * Return the transverse mass (squared), calculated from the mass
   * and the transverse momentum.
   */
  Energy2 perpmass2() const { return momentum().perp2() + momentum().mass2(); }

  /**
   * Return the transverse mass (squared), calculated from the mass
   * and the transverse momentum.
   */
  Energy perpmass() const { return sqrt(perpmass2()); }

  /**
   * Return the (pseudo) rapidity.
   */
  double rapidity() const {
    return ( Pplus() > ZERO && Pminus() > ZERO )?
      0.5*log(Pplus()/Pminus()) : Constants::MaxFloat;
  }

  /**
   * Return the (pseudo) rapidity.
   */
  double eta() const {
    Energy rho = momentum().rho();
    return rho > abs(momentum().z())?
      0.5*log((rho+momentum().z())/(rho-momentum().z())) : Constants::MaxFloat;
  }

  /**
   * Return the positive and negative light-cone momenta.
   */
  Energy Pplus() const { return momentum().plus(); }
  /**
   * Return the positive and negative light-cone momenta.
   */
  Energy Pminus() const { return momentum().minus(); }
  //@}

  /** @name Functions to access the position. */
  //@{
  /**
   * The creation vertex of this particle. The point is given
   * relative to the collision vertex.
   */
  const LorentzPoint & vertex() const {
    static const LorentzPoint null;
    return hasRep() ? rep().theVertex : null;
  }

  /**
   * The creation vertex of this particle. The absolute
   * position in the lab is given.
   */
  LorentzPoint labVertex() const;

  /**
   * The decay vertex of this particle. The point is given
   * relative to the collision vertex.
   */
  LorentzPoint decayVertex() const { 
    return vertex() + lifeLength();
  }

  /**
   * The decay vertex of this particle. The absolute
   * position in the lab is given.
   */
  LorentzPoint labDecayVertex() const {
    return labVertex() + lifeLength();
  }

  /**
   * The life time/length. Return the Lorentz vector connecting the
   * creation to the decay vertes.
   */
  const Lorentz5Distance & lifeLength() const {
    static const Lorentz5Distance null;
    return hasRep() ? rep().theLifeLength : null;
  }

  /**
   * Set the creation vertex relative to the collision vertex.
   */
  void setVertex(const LorentzPoint & p) {
    rep().theVertex = p;
  }

  /**
   * Set the creation vertex in the lab frame of this particle.
   */
  void setLabVertex(const LorentzPoint &);

  /**
   * Set the life length of this particle. The life time will be
   * automatically rescaled to be consistent with the invariant
   * distance.
   */
  void setLifeLength(const Distance & d) {
    rep().theLifeLength.setVect(d);
    rep().theLifeLength.rescaleEnergy();
  }

  /**
   * Set the life time/length of a particle. The invariant distance
   * may become inconsistent.
   */
  void setLifeLength(const LorentzDistance & d) {
    rep().theLifeLength = d;
  }

  /**
   * Set the life time/length of a particle.
   */
  void setLifeLength(const Lorentz5Distance & d) {
    rep().theLifeLength = d;
  }

  /**
   * The invariant life time of this particle.
   */
  Time lifeTime() const { return lifeLength().m(); }

  //@}

  /** @name Functions for (Lorentz) transformations. */
  //@{
  /**
   * Do Lorentz transformations on this particle.
   */
  void transform(const LorentzRotation & r);

  /**
   * Do Lorentz transformations on this particle. \a bx, \a by and \a
   * bz are the boost vector components.
   */
  void boost(double bx, double by, double bz) {
    transform(LorentzRotation(Boost(bx, by, bz)));
  }

  /**
   * Do Lorentz transformations on this particle. \a b is the boost
   * vector.
   */
  void boost(const Boost & b) { transform(LorentzRotation(b)); }

  /**
   * Rotate around the x-axis.
   */
  void rotateX(double a);

  /**
   * Rotate around the y-axis.
   */
  void rotateY(double a);

  /**
   * Rotate around the z-axis.
   */
  void rotateZ(double a);

  /**
   * Rotate around the given \a axis.
   */
  void rotate(double a, const Axis & axis);

  /**
   * Mirror in the xy-plane.
   */
  void mirror() { theMomentum.setZ(-theMomentum.z()); }

  /**
   * Do Lorentz transformations on this particle and its decendants.
   */
  void deepTransform(const LorentzRotation & r);

  /**
   * Do Lorentz transformations on this particle and its
   * decendants. \a bx, \a by and \a bz are the boost vector
   * components.
   */
  void deepBoost(double bx, double by, double bz) {
    deepTransform(LorentzRotation(Boost(bx, by, bz)));
  }

  /**
   * Do Lorentz transformations on this particle and its
   * decendants. \a b is the boost vector.
   */
  void deepBoost(const Boost & b) { deepTransform(LorentzRotation(b)); }

  /**
   * Rotate this particle and its decendants around the x-axis.
   */
  void deepRotateX(double a);

  /**
   * Rotate this particle and its decendants around the y-axis.
   */
  void deepRotateY(double a);

  /**
   * Rotate this particle and its decendants around the z-axis.
   */
  void deepRotateZ(double a);

  /**
   * Rotate this particle and its decendants around the given \a axis.
   */
  void deepRotate(double a, const Axis & axis);

  //@}

  /** @name Functions controlling possible mass/momentum inconsistencies. */
  //@{
  /**
   * Return the relative inconsistency in the mass component.
   */
  double massError() const { return theMomentum.massError(); }

  /**
   * Return the relative inconsistency in the energy component.
   */
  double energyError() const { return theMomentum.energyError(); }

  /**
   * Return the relative inconsistency in the spatial components.
   */
  double rhoError() const { return theMomentum.rhoError(); }

  /**
   * Rescale energy, so that the invariant length/mass of the
   * LorentzVector agrees with the current one.
   */
  void rescaleEnergy() { theMomentum.rescaleEnergy(); }

  /**
   * Rescale spatial component, so that the invariant length/mass of
   * the LorentzVector agrees with the current one.
   */
  void rescaleRho() { theMomentum.rescaleRho(); }

  /**
   * Set the invariant length/mass member, so that it agrees with the
   * invariant length/mass of the LorentzVector.
   */
  void rescaleMass() { theMomentum.rescaleMass(); }
  //@}

  /** @name Acces incormation about colour connections */
  //@{
  /**
   * True if this particle has colour information. To determine if
   * this particle is actually coloured, the coloured(), hasColour() or
   * hasAntiColour() methods should be used instead.
   */
  bool hasColourInfo() const {
    return hasRep() && rep().theColourInfo;
  }

  /**
   * Return the colour lines to which this particles anti-colour is
   * connected.
   */
  tColinePtr antiColourLine() const {
    return hasColourInfo() ? colourInfo()->antiColourLine() : tColinePtr();
  }

  /**
   * Return the colour lines to which this particles (\a anti-)colour
   * is connected.
   */
  tColinePtr colourLine(bool anti = false) const {
    if ( anti ) return antiColourLine();
    return hasColourInfo() ? colourInfo()->colourLine() : tColinePtr();
  }

  /**
   * Return true if the particle is connected to the given (\a anti-)
   * colour \a line.
   */
  bool hasColourLine(tcColinePtr line, bool anti = false) const {
    return hasColourInfo() ? colourInfo()->hasColourLine(line, anti) : false;
  }

  /**
   * Return true if the particle is connected to the given anti-colour
   * \a line.
   */
  bool hasAntiColourLine(tcColinePtr line) const {
    return hasColourLine(line, true);
  }

  /**
   * True if this particle type is not a colour singlet.
   */
  bool coloured() const { return data().coloured(); }

  /**
   * True if this particle type carries (\a anti-)colour.
   */
  bool hasColour(bool anti = false) const { return data().hasColour(anti); }

  /**
   * True if this particle type carries anti-colour.
   */
  bool hasAntiColour() const { return data().hasAntiColour(); }

  /**
   * Get the ColourBase object.
   */
  tcCBPtr colourInfo() const {
    return hasRep() ? rep().theColourInfo : CBPtr();
  }

  /**
   * Get the ColourBase object.
   */
  tCBPtr colourInfo() {
    if ( !rep().theColourInfo ) {
      switch(theData->iColour()) {
      case PDT::Colour6: 
      case PDT::Colour6bar:
	rep().theColourInfo = new_ptr(MultiColour());
	break;
      default:
	rep().theColourInfo = new_ptr(ColourBase());
      }
    }
    return rep().theColourInfo;
  }

  /**
   * Set the ColourBase object.
   */
  void colourInfo(tCBPtr c) {
    rep().theColourInfo = c;
  }

  /**
   * Get a pointer to the colour neighbor. Returns a particle in the
   * range \a first to \a last which colour is connected to the same
   * line as this particles anti-colour. If \a anti is true return
   * antiColourNeighbour().
   */
  template <typename Iterator>
  typename std::iterator_traits<Iterator>::value_type
  colourNeighbour(Iterator first, Iterator last, bool anti = false) const;

  /**
   * Get a pointer to the anti-colour neighbor. Returns a particle in
   * the range \a first to \a last which anti-colour is
   * connected to the same line as this particles colour.
   */
  template <typename Iterator>
  typename std::iterator_traits<Iterator>::value_type
  antiColourNeighbour(Iterator first, Iterator last) const {
    return colourNeighbour(first, last, true);
  }

  /**
   * Set the colour neighbor. Connects the given particles colour to
   * the same colour line as this particles anti-colour. If \a anti is
   * true call antiColourNeighbour(tPPtr).
   */
  void colourNeighbour(tPPtr, bool anti = false);

  /**
   * Set the anti-colour neighbor. Connects the given particles
   * anti-colour to the same colour line as this particles colour.
   */
  void antiColourNeighbour(tPPtr p) { colourNeighbour(p, true); }

  /**
   * Connect colour. Create a colour line connecting to it this
   * particles colour and the given particles anti-colour.
   */
  void antiColourConnect(tPPtr neighbour) {
    colourConnect(neighbour, true);
  }

  /**
   * Connect colour. Create a colour line connecting to it this
   * particles anti-colour and the given particles colour. If \a anti
   * is true call antiColourConnect(tPPtr).
   */
  void colourConnect(tPPtr neighbour, bool anti = false) {
    colourNeighbour(neighbour, anti);
  }

  /**
   * Incoming colour. Return the parent particle which colour is
   * connected to the same colour line as this particle. If \a anti is
   * true return incomingAntiColour().
   */
  tPPtr incomingColour(bool anti = false) const;

  /**
   * Incoming anti-colour. Return the parent particle which
   * anti-colour is connected to the same colour line as this
   * particle.
   */
  tPPtr incomingAntiColour() const { return incomingColour(true); }

  /**
   * Set incoming colour. Connect this particles colour to the same
   * colour line as the given particle. If \a anti
   * is true call incomingAntiColour(tPPtr).
   */
  void incomingColour(tPPtr p, bool anti = false) { p->outgoingColour(this, anti); }

  /**
   * Set incoming anti-colour. Connect this particles anti colour to
   * the same colour line as the given particle.
   */
  void incomingAntiColour(tPPtr p) { p->outgoingColour(this, true); }

  /**
   * Outgoing colour. Return the daughter particle which colour is
   * connected to the same colour line as this particle. If \a anti is
   * true return outgoingAntiColour().
   */
  tPPtr outgoingColour(bool anti = false) const;
  /**
   * Outgoing anti-colour. Return the daughter particle which
   * anti-colour is connected to the same colour line as this
   * particle.
   */
  tPPtr outgoingAntiColour() const { return outgoingColour(true); }

  /**
   * Set outgoing colour. Connect this particles colour to the same
   * colour line as the given particle. If \a anti
   * is true call outgoingAntiColour(tPPtr).
   */
  void outgoingColour(tPPtr, bool anti = false);

  /**
   * Set outgoing anti-colour. Connect this particles anti-colour to
   * the same colour line as the given particle.
   */
  void outgoingAntiColour(tPPtr p) { outgoingColour(p, true); }

  /**
   * Specify colour flow. Calls outgoingColour(tPPtr,bool).
   */
  void colourFlow(tPPtr child, bool anti = false) {
    outgoingColour(child, anti);
  }

  /**
   * Specify anticolour flow. Calls outgoingAntiColour(tPPtr,bool).
   */
  void antiColourFlow(tPPtr child) { colourFlow(child, true); }

  /**
   * Remove all colour information;
   */
  void resetColour() {
    if ( hasColourInfo() ) rep().theColourInfo = CBPtr();
  }

  //@}

  /** @name Functions to access spin. */
  //@{
  /**
   * Return the Spin object.
   */
  tcSpinPtr spinInfo() const { 
    return hasRep() ? rep().theSpinInfo : SpinPtr(); 
  }

  /**
   * Return the Spin object.
   */
  tSpinPtr spinInfo() {
    return hasRep() ? rep().theSpinInfo : SpinPtr(); 
  }

  /**
   * Set the Spin object.
   */
  void spinInfo(tSpinPtr s) { rep().theSpinInfo = s; }
  //@}

  /** @name Accessing user-defined information. */
  //@{
  /**
   * Access user-defined information as a vector of EventInfoBase pointers.
   */
  const EIVector & getInfo() const {
    static const EIVector null;
    return hasRep() ? rep().theExtraInfo : null;
  }

  /**
   * Access user-defined information as a vector of EventInfoBase pointers.
   */
  EIVector & getInfo() { return rep().theExtraInfo; }
  //@}

public:

  /** @name Accessing user-defined information. */
  //@{
  /**
   * True if this particle has instantiated the object with
   * information other than type and momentum.
   */
  bool hasRep() const { return theRep; }

  /**
   * If this particle has only a type and momentum, instantiate the
   * rest of the information.
   */
  void initFull();

  //@}

public:

  /** @name Input and output functions. */
  //@{
  /**
   * Standard function for writing to a persistent stream.
   */
  void persistentOutput(PersistentOStream &) const;

  /**
   * Standard function for reading from a persistent stream.
   */
  void persistentInput(PersistentIStream &, int);

  //@}

  /**
   * Print particle info to a stream \a os. The \a step is used to
   * access information about colour neighbors and other struff.
   */
  ostream & print(ostream & os, tcStepPtr step = tcStepPtr()) const;

  /**
   * Print a range of particles.
   */
  template <typename Iterator>
  static void PrintParticles(ostream & os, Iterator first, Iterator last,
			     tcStepPtr step = tcStepPtr());

  /**
   * Print a container of particles.
   */
  template <typename Cont>
  static inline void PrintParticles(ostream & os, const Cont & c,
			     tcStepPtr step = tcStepPtr()) {
    PrintParticles(os, c.begin(), c.end(), step);
  }
   
  /**
   * Standard Init function. @see Base::Init().
   */
  static void Init();

  /**
   * Specify how to print particles. The format string is analogous to
   * the one used by eg. the unix 'date' command as described above.
   */
  static string outputFormat;

private:

  /**
   * Standard clone function.
   */
  virtual PPtr clone() const;

  /**
   * Rebind to cloned objects. When an Event is cloned, a shallow
   * copy is done first, then all <code>Particle</code>s etc, are
   * cloned, and finally this method is used to see to that the
   * pointers in the cloned Particle points to the cloned objects.
   */
  virtual void rebind(const EventTranslationMap &);

  /**
   * Set the order-number for this particle in the current event.
   */
  void number(int n) { rep().theNumber = n; }

  /**
   * Remove the given particle from the list of children.
   */
  void removeChild(tPPtr c) {
    if ( hasRep() )
      rep().theChildren.erase(remove(rep().theChildren.begin(),
				     rep().theChildren.end(), c),
			      rep().theChildren.end());
  }

  /**
   * Remove the given particle from the list of parents.
   */
  void removeParent(tPPtr p) {
    if ( hasRep() )
      rep().theParents.erase(remove(rep().theParents.begin(),
				    rep().theParents.end(), p),
			     rep().theParents.end());
  }

  /**
   * Set the mass of this particle.
   */
  void mass(Energy m) { theMomentum.setMass(m); }

  /**
   * Set the invaiant life time of this particle.
   */
  void lifeTime(Length t) { rep().theLifeLength.setTau(t); }

  /**
   * Return a reference to the bulk information of this particle. if
   * no ParticleRep object exists, one is created.
   */
  ParticleRep & rep() {
    if ( !hasRep() ) initFull();
    return *theRep;
  }

  /**
   * Return a reference to the bulk information of this particle. if
   * no ParticleRep object exists, we return the default values.
   */
  const ParticleRep & rep() const {
    static const ParticleRep null;
    return hasRep() ? *theRep : null;
  }

  /**
   * The pointer to the ParticleData object
   */
  cEventPDPtr theData;

  /**
   * The momentum.
   */
  Lorentz5Momentum theMomentum;

  /**
   * The rest of the information in this particle is only instantiated
   * if needed.
   */
  ParticleRep * theRep;

public:

  /**
   * This class is used internally in the Particle class to represent
   * information besides momentum and type. A corresponding object
   * will only be instantiated if needed to save memory and time when
   * temporarily creating particles.
   */
  struct ParticleRep {

    /**
     * Default constructor.
     */
    ParticleRep() : theScale(-1.0*GeV2), theVetoScale(-1.0*GeV2), theNumber(0) {}

    /**
     * Copy constructor.
     */
    ParticleRep(const ParticleRep &);

    /**
     * The pointers to the parents.
     */
    tParticleVector theParents;

    /**
     * The pointers to the children.
     */
    ParticleVector theChildren;

    /**
     * The pointer to the previous instance.
     */
    tPPtr thePrevious;

    /**
     * The pointer to the next instance.
     */
    PPtr theNext;

    /**
     * If this particle has decayed this is the pointer to the
     * corresponding decay mode.
     */
    tDMPtr theDecayMode;

    /**
     * The pointer to the first step where this particle occurred.
     */
    tStepPtr theBirthStep;

    /**
     * The creation point.
     */
    LorentzPoint theVertex;

    /**
     * The life time/length.
     */
    Lorentz5Distance theLifeLength;

    /**
     * the resolution scale.
     */
    Energy2 theScale;

    /**
     * the veto scale.
     */
    Energy2 theVetoScale;

    /**
     * The order-number for this particle in the current event.
     */
    int theNumber;

    /**
     * A pointer to the colour information object.
     */
    CBPtr theColourInfo;

    /**
     * Spin information
     */
    SpinPtr theSpinInfo;

    /**
     * Additional used-defined information.
     */
    EIVector theExtraInfo;

  };

public:

  /**
   * Print out debugging information for this object on std::cerr. To
   * be called from within a debugger via the debug() function.
   */
  virtual void debugme() const;

protected:

  /**
   * Private default constructor must only be used by the
   * PersistentIStream class via the ClassTraits<Particle> class.
   */
  Particle() : theRep(0) {}

  /**
   * The ClassTraits<Particle> class must be a friend to be able to
   * use the private default constructor.
   */
  friend struct ClassTraits<Particle>;

private:

  /**
   * Private and non-existent assignment.
   */
  Particle & operator=(const Particle &);

  /**
   * Describe concrete class with persistent data.
   */
  static ClassDescription<Particle> initParticle;

};

/**
 * Write a Particle object to a stream.
 */
ostream & operator<<(ostream &, const Particle &);


/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base class of Particle. */
template <>
struct BaseClassTrait<Particle,1>: public ClassTraitsType {
  /** Typedef of the first base class of Collision. */
  typedef EventRecordBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the Particle class and how to create it. */
template <>
struct ClassTraits<Particle>: public ClassTraitsBase<Particle> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::Particle"; }
  /** Create a Particle object. */
  static TPtr create() { return TPtr::Create(Particle()); }
};

/** @endcond */

}

#ifndef ThePEG_TEMPLATES_IN_CC_FILE
#include "Particle.tcc"
#endif

#endif /* ThePEG_Particle_H */
