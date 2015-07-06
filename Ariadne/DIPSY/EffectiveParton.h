// -*- C++ -*-
#ifndef DIPSY_EffectiveParton_H
#define DIPSY_EffectiveParton_H
//
// This is the declaration of the EffectiveParton class.
//

#include "Ariadne/DIPSY/Parton.h"
#include "Ariadne/DIPSY/EffectiveParton.fh"

namespace DIPSY {

using namespace ThePEG;

/**
 * An effective parton is a colour chain of partons that are close
 * enough to each other so that the neighbouring dipole can't really
 * resolve them.
 */
class EffectiveParton: public Parton {

public:

  /**
   * Caching info for Partons.
   */
  struct RangeInfo {

    /**
     * The range at which the number of partons included changes.
     */
    InvEnergy2 range;

    /**
     * The effective transverse momentum at this range;
     */
    TransverseMomentum pT;

    /**
     * The effective positive light-cone momentum at this range.
     */
    Energy plus;

    /**
     * The effective negative light-cone momentum at this range.
     */
    Energy minus;

    /**
     * The effective rapidity for this range.
     */
    double y;

    /**
     * Included partons.
     */
    vector<tPartonPtr> partons;

    /**
     * Included dipoles.
     */
    vector<tDipolePtr> dipoles;

    /**
     * The external connecting dipoles.
     */
    pair<tDipolePtr,tDipolePtr> dipair;

  };



public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  EffectiveParton();

  /**
   * Constructs an effective parton for \a p.
   */
  EffectiveParton(Parton & p);

  /**
   * Static creation method that also performs
   * initialization. Constructs an effective parton around p,
   * including other partons within range.
   */
  static EffectivePartonPtr create(Parton & p, InvEnergy range);

  /**
   * The destructor.
   */
  virtual ~EffectiveParton();
  //@}

protected:

  /**
   * Setup the and access cache of range information.
   */
  void setupCache(InvEnergy maxrange);
  void uncache();
  void recache();
  /** @name The virtual functions to be overridden in sub-classes. */
  //@{
  /**
   * Return a simple clone of this object. Should be implemented as
   * <code>return new_ptr(*this);</code> by a derived class.
   */
  virtual Ariadne5::ClonePtr clone() const;

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

  /**
   * Returns the internal dipoles in the nonresolved chain.
   */
  inline const vector<tDipolePtr> & internalDipoles() {
    return cache.empty()? theInternalDipoles: cache[currange].dipoles;
  }

  /**
   * Returns the internal dipoles in the nonresolved chain.
   */
  inline tPartonPtr originalParton() {
    return theOriginalParton;
  }

  /**
   * recalculates the partons and dipoles within the range, and updates
   * theInternalDipoles and momenta accordingly.
   */
  void printData();

  /**
   * recalculates the partons and dipoles within the range, and updates
   * theInternalDipoles and momenta accordingly.
   */
    void setRange( InvEnergy range );

  /**
   * returns true if any of the internal dipoles, or their neighbours,
   * has emitted.
   */
  bool changed();

  /**
   * Finds and returns the partons belonging to the internal dipoles.
   */
  const set<tPartonPtr> & internalPartons() const {
    return theInternalPartons;
  };

  /**
   * Finds and returns the partons belonging to the internal dipoles.
   */
  const vector<tPartonPtr> & cachedPartons() const {
    return cache.back().partons;
  };

  /**
   * Check if a recoil of pT and plus will push any of the internal partons above
   * rapidity ymax.
   */
  bool recoilsOverYMax(TransverseMomentum pT, Energy plus, double ymax) const;

  /**
   * recoils the effective parton from providing an emission with
   * pT transverse momentum, and plus p+. Goes with the pT definition
   * where every recoil in every emission is remembered.
   */
  void recoil( TransverseMomentum pT, Energy plus );

  /**
   * recoils the effective parton from providing an emission with
   * plus p+. pT is decided from each partons two dipole neighbours.
   * Goes with the pT determined only from the two current
   * colour neighbours, without history dependence.
   */
  void newRecoil( Energy plus );

protected:

  /**
   * recursively add partons within range.
   */
  void addPartons( InvEnergy range );

  /**
   * recursively add partons within range. USes the "Relatives" mode of adding partons.
   */
  void addRelatives( InvEnergy range, tPartonPtr p );

  /**
   * Adds p to the internal partons, and updates the momentum of the effective parton.
   */
  void merge(tPartonPtr p);

  /**
   * checks that the sum ofh the partons momenta is the momenta of the effective parton.
   */
  bool checkSums() const;

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
   * Exception class used when effective parton gets (almost) zero pt)
   */
  struct EffectivePartonPTException: public Exception {};

private:

  /**
   * The original parton.
   */
  PartonPtr theOriginalParton;

  /**
   * The internal dipoles in the nonresolved chain.
   */
  vector<tDipolePtr> theInternalDipoles;

  /**
   * The internal partons in the nonresolved effective parton.
   */
  set<tPartonPtr> theInternalPartons;

  /**
   * Cache informations about different ranges.
   */
  vector<RangeInfo> cache;

  /**
   * The current range 
   */
  int currange;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  EffectiveParton & operator=(const EffectiveParton &);

};

}

#endif /* DIPSY_EffectiveParton_H */
