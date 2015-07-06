// -*- C++ -*-
#ifndef Ariadne5_Junction_H
#define Ariadne5_Junction_H
//
// This is the declaration of the Junction class.
//

#include "Parton.h"
#include "Junction.fh"
namespace Ariadne5 {

using namespace ThePEG;

/**
 * The Junction class inherits from Parton, but is only a place holder
 * so that Another Parton can be colour-connected to it. The class
 * keeps track of the other partons which are connected to the
 * Junction.
 */
class Junction: public Parton {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  Junction();

  /**
   * The destructor.
   */
  virtual ~Junction();
  //@}


public:

  /**
   * Return true if this is a source.
   */
  bool source() const;

  /**
   * Return true if this is a sink.
   */
  bool sink() const;

  /**
   * Get the Parton which is assumed to be connected to the given
   * Parton.
   */
  tParPtr getNeighbour(tcParPtr p) const;

  /**
   * Return true if any of the neighboring dipoles or partons have
   * been touched.
   */
  virtual bool touchedNeighbours(tcQCDPtr d) const;

  /**
   * Get the neighboring dipoles.
   */
  pair<tcQCDPtr,tcQCDPtr> getNeighbours(tcQCDPtr d) const;

  /**
   * Randomly choose a neighboring parton in the junction, recursing
   * if this itself was a junction.
   */
  tParPtr getRandomRecoiler(tcQCDPtr d) const;

  /**
   * Replace the \a oldd QCDDipole with the \a newd QCDDipole as connected
   * to this Junction.
   * @return false if \a oldd was not a member.
   */
  bool replace(tQCDPtr oldd, tQCDPtr newd);

  /**
   * Set the three partons connected to this Junction.
   */
  void set(tParPtr pin1, tParPtr pin2, tParPtr pin3) {
    p1 = pin1;
    p2 = pin2;
    p3 = pin3;
  }

  /**
   * Do not produce a ThePEG::Particle corresponding to this parton,
   * as it is not really a parton.
   */
  virtual tPPtr produceParticle(const LorentzRotation & r = LorentzRotation());

protected:

  /**
   * Helper function for recursing touchedNeighbours().
   */
  bool touchedBranch(tcQCDPtr d) const;

  /** @name Functions relating to the DipoleState to which this belongs. */
  //@{
  /**
   * Return a simple clone of this object. Should be implemented as
   * <code>return new_ptr(*this);</code> by a derived class.
   */
  virtual ClonePtr clone() const;

  /**
   * Fill the provided set with all pointers to CloneBase objects used
   * in this object.
   */
  virtual void fillReferences(CloneSet &) const;

  /**
   * Rebind pointers to other CloneBase objects. Called after a number
   * of interconnected CloneBase objects have been cloned, so that
   * the cloned objects will refer to the cloned copies afterwards.
   *
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   */
  virtual void rebind(const TranslationMap & trans);
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
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

private:

  /**
   * The first dipole.
   */
  tQCDPtr d1;

  /**
   * The second dipole.
   */
  tQCDPtr d2;

  /**
   * The third dipole.
   */
  tQCDPtr d3;

  /**
   * The first parton.
   */
  tParPtr p1;

  /**
   * The second parton.
   */
  tParPtr p2;

  /**
   * The third parton.
   */
  tParPtr p3;


public:

  /**
   * Print out debugging information on std::cerr.
   */
  virtual void debugme() const;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  Junction & operator=(const Junction &);

};

}

#endif /* Ariadne5_Junction_H */
