// -*- C++ -*-
//
// Graphviz.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_GraphvizPlot_H
#define THEPEG_GraphvizPlot_H
//
// This is the declaration of the GraphvizPlot class.
//

#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Handlers/AnalysisHandler.h"

namespace ThePEG {

/** \ingroup Analysis
 * The GraphvizPlot class generates
 * an output of the tree structure of the event which can be viewed using dot.
 *
 * @see \ref GraphvizPlotInterfaces "The interfaces"
 * defined for GraphvizPlot.
 */
class GraphvizPlot: public AnalysisHandler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  GraphvizPlot() : _eventNumber(1), _quiet(false) {}
  //@}

public:

  /** @name Virtual functions required by the AnalysisHandler class. */
  //@{
  /**
   * Analyze a given Event. Note that a fully generated event
   * may be presented several times, if it has been manipulated in
   * between. The default version of this function will call transform
   * to make a lorentz transformation of the whole event, then extract
   * all final state particles and call analyze(tPVector) of this
   * analysis object and those of all associated analysis objects. The
   * default version will not, however, do anything on events which
   * have not been fully generated, or have been manipulated in any
   * way.
   * @param event pointer to the Event to be analyzed.
   * @param ieve the event number.
   * @param loop the number of times this event has been presented.
   * If negative the event is now fully generated.
   * @param state a number different from zero if the event has been
   * manipulated in some way since it was last presented.
   */
  virtual void analyze(tEventPtr event, long ieve, int loop, int state);
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

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const {return new_ptr(*this);}

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const {return new_ptr(*this);}
  //@}

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  virtual void dofinish();
  //@}

  /**
   * Helper function to obtain the name of a particle.
   */
  string particleName(tcPPtr) const;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<GraphvizPlot> initGraphvizPlot;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  GraphvizPlot & operator=(const GraphvizPlot &) = delete;

private:

  /**
   * Event number that should be drawn 
   */
  long _eventNumber;

  /**
   * Tell the object not to write out messages to the standard output.
   */
  bool _quiet;


};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of GraphvizPlot. */
template <>
struct BaseClassTrait<GraphvizPlot,1> {
  /** Typedef of the first base class of GraphvizPlot. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the GraphvizPlot class and the shared object where it is defined. */
template <>
struct ClassTraits<GraphvizPlot>
  : public ClassTraitsBase<GraphvizPlot> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::GraphvizPlot"; }
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the GraphvizPlot class and any other class on which it depends
   *  (except the base class). */
  static string library() { return "GraphvizPlot.so"; }
};

/** @endcond */

}

#endif /* THEPEG_GraphvizPlot_H */
