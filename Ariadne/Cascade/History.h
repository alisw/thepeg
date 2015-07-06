// -*- C++ -*-
#ifndef Ariadne5_History_H
#define Ariadne5_History_H
//
// This is the declaration of the History class.
//

#include "DipoleState.h"
#include "Emission.h"
#include "EmitterBase.h"
#include "AriadneHandler.fh"
#include "History.fh"
#include "ThePEG/Utilities/Selector.h"

namespace Ariadne5 {

using namespace ThePEG;

/**
 * Here is the documentation of the History class.
 */
class History: public Base {

public:

  /**
   * Enumerate different ways of reordering scales in unordered paths.
   */
  enum Reordering {
    UseMax = 1, /**< Always take the maximum. */
    UseMin = -1 /**< Always take the minimum. */
  };

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor, si never useful.
   */
  History();

  /**
   * The main constructor, which recursively constructs a tree
   * structure of History objects. \a depth is the maximum number of
   * clusterings requested. \a statein is the DipoleState to be
   * clustered.
   */
  History(int depth, DipoleStatePtr statein);

protected:

  /**
   * The constructor used when creating recursive nodes. The \a depth
   * is the maximum number of clusterings requested. \a scalein is the
   * scale at which the \a statein was clustered. If \a isOrdered is
   * true, the previous clusterings has been ordered. \a probin is the
   * accumulated probabilities for the previous clusterings, and \a
   * mothin is the parent history node. \a emin is the inverse
   * emission which produced the \a statein.
   */
  History(int depth, Energy scalein, DipoleStatePtr statein,
	  bool isOrdered, double probin, tHistoryPtr mothin, tcEmPtr emin);

public:

  /**
   * The destructor.
   */
  virtual ~History();
  //@}

public:

  /**
   * In the initial history node, select one of the paths according to
   * the probabilities. This function must only be called for the
   * initial history node. All reconstructed paths, except the chosen
   * one, are discarded.
   */
  tHistoryPtr select(Reordering ord = UseMax);

  /**
   * Find the weight calculated from the ratio of couplings, the
   * no-emission probabilities, and possible PDF ratios. \a as0 is the
   * value of \f$\alpha_s\f$ and \a muF is the factorization scale
   * used in the matrix element generator. This function must only be
   * called from the original node.
   */
  double weight(double as0, Energy2 muF2) const;

  /**
   * The scale of this step, corresponding to clustering which
   * constructed the corresponding state. For the initial node, this
   * is the total energy before selection of the path and the
   * factorization scale of the fully resconstructed state after the
   * selection.
   */
  Energy scale() const {
    return theScale;
  }

  /**
   * The state of the event correponding to this step in the
   * reconstruction.
   */
  tDipoleStatePtr state() const {
    return theState;
  }

private:

  /**
   * Do the actual recursive clustering, with a maximum number of
   * steps given by \a depth. If \a isOrdered is false, this belongs
   * to a path which is already unordered.
   */
  void recurse(int depth, bool isOrdered);

  /**
   * Check if an ordered (and complete) path has been found in the
   * initial node, in which case we will no longer be interested in
   * any unordered paths.
   */
  bool onlyOrderedPaths() const;

  /**
   * When a full path has been found, register it with the initial
   * history node.
   */
  bool registerPath(tHistoryPtr h, Energy rho0,
		    bool isOrdered, bool isComplete);

  /**
   * Helper class to order potential reconstructions in scale.
   */
  struct ScaleCmp {
    /** Compare two reconstructions. */
  bool operator()(tcEmPtr e1, tcEmPtr e2) { return e1->rho < e2->rho; }
  };

  /**
   * Return the set of all possible single clusterings of the current
   * state.
   */
  set<EmPtr,ScaleCmp> getClusterings();

  /**
   * Perform the clustering of the current state and return the
   * clustered state.
   */
  DipoleStatePtr cluster(tEmPtr);

  /**
   * Perform a trial shower between \a maxscale down to this scale and
   * return the corresponding Sudakov form factor.
   */
  double trialShower(Energy maxscale) const;

  /**
   * Remove all children from this node, except the given one.
   */
  void clean(HistoryPtr h = HistoryPtr());

  /**
   * Find the weight calculated from the ratio of couplings and
   * possible PDF ratios. \a as0 is the value of \f$\alpha_s\f$ used
   * in the matrix element generator. \a hdl is the resonsible
   * handler.
   */
  double weight(double as0, const AriadneHandler & hdl) const;

  /**
   * Reorder scales for an unordered path.
   */
  void reorder(Reordering);

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
   * The state of the event correponding to this step in the
   * reconstruction.
   */
  mutable DipoleStatePtr theState;

  /**
   * The previous step from which this step has been clustered. If
   * null, this is the initial step with the n-jet state generated by
   * the matrix element.
   */
  tHistoryPtr mother;

  /**
   * The inverse emission which reconstructed this state.
   */
  tcEmPtr emission;

  /**
   * The different steps corresponding to possible clusterings of this
   * state.
   */
  vector<HistoryPtr> children;

  /**
   * Selector keeping track of all reconstructed paths according to
   * their probabilities. This is empty unless this is the initial
   * step (mother == null).
   */
  Selector<tHistoryPtr, double> paths;

  /**
   * A map keeping track of the maximum scales for each path.
   */
  map<tHistoryPtr,Energy> scales;
  
  /**
   * This is set true if an ordered (and complete) path has been found
   * and inserted in paths.
   */
  mutable bool foundOrderedPath;

  /**
   * This is set true if a complete (with the required number of
   * clusterings) path has been found and inserted in paths.
   */
  bool foundCompletePath;

  /**
   * The scale of this step, corresponding to clustering which
   * constructed the corresponding state is zero for the initial node.
   */
  Energy theScale;

  /**
   * The probability associated with this step and the previous steps.
   */
  double prob;

  /**
   * The last node in the selected path.
   */
  HistoryPtr sel;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  History & operator=(const History &);

};

}

#endif /* Ariadne5_History_H */
