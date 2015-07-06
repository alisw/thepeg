// -*- C++ -*-
#ifndef DIPSY_Dipole_H
#define DIPSY_Dipole_H
//
// This is the declaration of the Dipole class.
//

#include "ThePEG/Config/ThePEG.h"
#include "Ariadne/Config/CloneBase.h"
#include "Dipole.fh"
#include "Parton.h"
#include "EffectiveParton.h"
#include "DipoleState.fh"

namespace DIPSY {

using namespace ThePEG;

/**
 * Here is the documentation of the Dipole class.
 */
class Dipole: public Ariadne5::CloneBase {

public:

  /**
   * A pair of dipoles
   */
  typedef pair<tDipolePtr,tDipolePtr> tDipolePair;

  /**
   * A pair of partons.
   */
  typedef pair<PartonPtr,PartonPtr> PartonPair;

  /**
   * The DipoleState is a friend.
   */
  friend class DipoleState;

    /**
     * A map for stoing effective partons.
     */
    typedef vector< pair<InvEnergy,EffectivePartonPtr> > EffMap;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline Dipole();

  /**
   * The copy constructor.
   */
  inline Dipole(const Dipole &);

  /**
   * The destructor.
   */
  virtual ~Dipole();
  //@}

protected:

  /** @name The virtual functions to be overridden in sub-classes. */
  //@{
  /**
   * Return a simple clone of this object. Should be implemented as
   * <code>return new_ptr(*this);</code> by a derived class.
   */
  virtual Ariadne5::ClonePtr clone() const;
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

  /** @name Simple access functions. */
  //@{
  /**
   * rMax for the dipoles handler.
   */
  InvEnergy rMax() const;

  /**
   * The DipoleState to which this Dipole belongs.
   */
  inline DipoleState & dipoleState() const;

  /**
   * Set the DipoleState to which this Dipole belongs.
   * added by CF to access from emitter.
   */
  inline void dipoleState(tDipoleStatePtr);

  /**
   */
  inline const PartonPair & partons() const;

  /**
   * Get the effective parton.
   */
  inline const pair< EffectivePartonPtr, EffectivePartonPtr > &
  effectivePartons() const;

  /**
   * Get an effective parton corresponding to one of the partons in
   * the original dipole given a range.
   */
  tEffectivePartonPtr getEff(tPartonPtr p, InvEnergy range) const;


  /**
   * Set the effective parton.
   */
  inline void effectivePartons( tEffectivePartonPtr, tEffectivePartonPtr );

  /**
   * controls if any of the effective partons have been changed by emissions.
   */
  inline const bool effectivePartonsChanged() const;

  /**
   * Set the neighboring dipoles.
   */
  inline const tDipolePair & neighbors() const;

  /**
   * Get the child dipoles.
   */
  inline const tDipolePair & children() const;

  /**
   * Get the child dipoles.
   */
  inline tDipolePair & children();

  /**
   * Get the generated gluon to be emitted. Returns null if none has
   * been generated.
   */
  inline tPartonPtr generatedGluon() const;

  /**
   * Set the generated gluon to be emitted.
   * Added by CF to get access from emitter.
   */
  inline void generatedGluon(tPartonPtr);

  /**
   * Get the dipole to swing with. Return null if no swing has been
   * generated.
   */
  inline tDipolePtr swingDipole() const;

  /**
   * Set the dipole to swing with.
   * Aded by CF
   */
  inline  void swingDipole(tDipolePtr);

  /**
   * Return true if this dipole has generated an emission or a swing.
   */
  inline bool hasGen();

  /**
   * The rapidity of the generated emission or swing.
   */
  inline double generatedY() const;

  /**
   * Set the rapidity of the generated emission or swing.
   * Added by CF to access from emitter.
   */
  inline void generatedY(double);

  /**
   * Get a pointer to a Dipole in another DipoleState with which this
   * Dipole has interacted. @return null if Dipole has not interacted.
   */
  inline tDipolePtr interacted() const;

  /**
   * Returns the lengths that the new dipoles will have when it has
   * swinged with the interacted() dipole. returns 0 if not interacting.
   */
  pair <InvEnergy, InvEnergy> interactionLengths() const;

  /**
   * Set the partons.
   */
  inline void partons(PartonPair);

  /**
   * Set the neighboring dipoles.
   */
  inline void neighbors(tDipolePair);

  /**
   * Set the neighboring dipole.
   */
  inline void firstNeighbor(tDipolePtr);

  /**
   * Set the neighboring dipole.
   */
  inline void secondNeighbor(tDipolePtr);

  /**
   * Get the colour.
   */
  inline int colour() const;

  /**
   * Set the colour.
   */
  inline void colour(int c, int sys = 0);

  /**
   * Get the colour system.
   */
  int colourSystem() const;

  /**
   * Set the colour system.
   */
  void colourSystem(int sys);

  /**
   * says if the dipole is part of a DGLAP chain, and thus should not be absorbed.
   */
  inline bool DGLAPsafe();

  /**
   * Sets if the dipole is part of a DGLAP chain or not.
   */
  inline void DGLAPsafe(bool);

  /**
   * says if the (initial) dipole is participating in the interaction.
   */
  inline bool participating();

  /**
   * Sets if the (initial) dipole is participating in the interaction.
   */
  inline void participating(bool);

  /**
   * gets if the generated swing is a gluon exchange.
   */
  inline bool recoilSwing() const {
    return isRecoilSwing;
  }

  /**
   * sets if the generated swing is a gluon exchange.
   */
  inline void recoilSwing(bool recoil) {
    isRecoilSwing = recoil;
  }

  /**
   * Indicate that the Dipole has interacted with the given Dipole in
   * another DipoleState.
   */
  void interact(Dipole &);

  /**
   * Get the squared transverse size of this Dipole.
   */
  inline InvEnergy2 size2() const;

  /**
   * Get the transverse size of this Dipole.
   */
  inline InvEnergy size() const;
  //@}

  /** @name Functions for extracting final dipoles. */
  //@{
  /**
   * Given an output iterator fill the corresponding container with
   * final, undecayed dipoles. This is done recursively, ie. if this
   * dipole has decayed, the extract methods of the two children are
   * called instead.
   */
  template <typename OutputIterator>
  void extract(OutputIterator it);

  /**
   * Extract all undecayed dipoles into the vector according to colour
   * index.
   */
  void colExtract(vector< vector<tDipolePtr> > &);

  //@}

  /**
   * Recursively find the dipole which is next to swing in final state.
   */
  tDipolePtr getFSSwinger(double miny, double maxy);

  /**
   * Generate a possible emission or a swing from this Dipole in the
   * given rapidity interval [\a miny,\a maxy].
   */
  void generateFSRec(double miny, double maxy);

  /** @name Functions for shower generation. */
  //@{
  /**
   * Recursively find the dipole which is next to emit or swing in the
   * given rapidity interval [\a miny,\a maxy].
   */
  tDipolePtr getEmitter(double miny, double maxy);

  /**
   * Generate a possible emission or a swing from this Dipole in the
   * given rapidity interval [\a miny,\a maxy].
   */
  void generate(double miny, double maxy);

  /**
   * Generate a swing for this Dipole in the given rapidity interval
   * [\a miny,\a maxy]. If \a force is true, always generate a swing,
   * otherwise only check if a swing is possible with dipoles which
   * are new or has changed.
   */
  void generateRec(double miny, double maxy, bool force);

  /**
   * Try really hard to generates a swing, but not with higher rapdity step
   * than /a ymax. /a ymax = 0 means no limit.
   * If no swing found, return false.
   */
  bool forceGenerateRec(double ymax);

  /**
   * Perform a recombination previously generated by generateRec().
   */
  void recombine();

  /**
   * absorb the dipole, replacing it with a single parton.
   */
  void absorb();
  
  /**
   * Perform the emission previously generated for this \a dipole. If
   * no emission has been generated a runtime_error is thrown.
   */
   void emit();
  //@}


public:

  /**
   * Split this dipole by emitting the generated gluon. \a colsel
   * gives the probability that the first child dipole inherits the
   * colour index of the mother.
   */
  void splitDipole(double colsel);

  /**
   * Reset any emission or swong generated for this dipole.
   */
  void reset();

  /**
   * Check if conditions for this dipole has changed since last
   * generated emission or swing.
   */
  bool touched() const {
    return isTouched;
  }

  /**
   * Flag that conditions for this dipole has changed since last
   * generated emission or swing.
   */
  void touch() {
    isTouched = true;
  }

  /**
   * Flag that this dipole no has generated emission or swing.
   */
  void untouch() {
    isTouched = false;
  }

  /**
   * Flag that this dipole should no longer emit.
   */
  void turnOff() {
    isOn = false;
  }


protected:

  struct NothingGenerated: public Exception {};

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
   * The DipoleState to which this Dipole belongs.
   */
  tDipoleStatePtr theDipoleState;

  /**
   * The partons.
   */
  PartonPair thePartons;

  /**
   * The neighboring dipoles.
   */
  tDipolePair theNeighbors;

  /**
   * The child dipoles.
   */
  tDipolePair theChildren;

  /**
   * The generated gluon to be emitted. Null if none has been generated.
   */
  PartonPtr theGeneratedGluon;

  /**
   * the effective partons.
   */
  pair< EffectivePartonPtr, EffectivePartonPtr > theEffectivePartons;

  /**
   * The dipole to swing with. Null if no swing has been generated.
   */
  tDipolePtr theSwingDipole;

  /**
   * The rapidity of the generated emission or swing.
   */
  double theGeneratedY;

  /**
   * The colour index of this Dipole.
   */
  int theColour;

  /**
   * If non-null, pointer to a Dipole in another DipoleState with which
   * this Dipole has interacted.
   */
  tDipolePtr theInteracted;

  /**
   * If the dipole is part of a DGLAP chain or not in the final state.
   */
  bool theGLAPsafe;

  /**
   * if the (initial) dipole is part of a participating nucleon in a
   * heavy ion collision.
   */
  bool isParticipating;

  /**
   * If the generated swing is gluon exchange or quadrupole swing.
   */
  bool isRecoilSwing;

  /**
   * Indicate if conditions for this dipole has changed since last
   * generated emission or swing.
   */
  bool isTouched;

  /**
   * This can be set to false to manually turn a dipole of, ie disallow
   * emissions or swings.
   **/
  bool isOn;

public:

  /**
   * Cache values needed for (FS) swing.
   */
  mutable InvEnergy2 swingCache;

  /**
   * A pair of maps where the effective partons are stored.
   */
  mutable pair<EffMap,EffMap> * effmap;

  /**
   * Build the map of effective partons.
   */
  void buildEffMap() const;
  
private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  Dipole & operator=(const Dipole &);

};

}

#include "Dipole.icc"

#endif /* DIPSY_Dipole_H */
