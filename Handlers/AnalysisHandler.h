// -*- C++ -*-
//
// AnalysisHandler.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_AnalysisHandler_H
#define ThePEG_AnalysisHandler_H
// This is the declaration of the AnalysisHandler class.

#include "HandlerBase.h"
#include "AnalysisHandler.fh"
#include "ThePEG/Vectors/LorentzRotation.h"
#include "ThePEG/Analysis/FactoryBase.h"
#include "ThePEG/EventRecord/Event.h"
#include <stdexcept>

namespace ThePEG {

/**
 * The AnalysisHandler is the base class of all analysis objects which
 * may be handled by the FullEventGenerator. The main function is the
 * virtual <code>analyze()</code> method which which is called for
 * each analysis handler after each event. The method may be called
 * several times for each event - this may be checked by the analysis
 * handler by looking at the <code>ieve</code>, <code>loop</code> and
 * <code>state</code> arguments to the <code>analyze</code> method.
 *
 * Initialization of histograms etc. should be made in the
 * <code>doinitrun()</code> function, while writing out of histograms
 * and analysis results should be performed in the
 * <code>dofinish()</code> function.
 *
 * @see \ref AnalysisHandlerInterfaces "The interfaces"
 * defined for AnalysisHandler.
 * @see FullEventGenerator
 * @see Event
 */
class AnalysisHandler: public HandlerBase {

public:

  /**
   * Convenient typedef for pointer to AIDA::IHistogram1D.
   */
  typedef FactoryBase::tH1DPtr tH1DPtr;

  /**
   * Convenient typedef for pointer to AIDA::IHistogram1D.
   */
  typedef FactoryBase::tcH1DPtr tcH1DPtr;

  /**
   * Convenient typedef for pointer to AIDA::IHistogram1D.
   */
  typedef FactoryBase::tH2DPtr tH2DPtr;

  /**
   * Convenient typedef for pointer to AIDA::IHistogram1D.
   */
  typedef FactoryBase::tcH2DPtr tcH2DPtr;

public:

  /** @name Virtual functions required by the AnalysisHandler class. */
  //@{
  /**
   * Analyze a given Event. Note that a fully generated event may be
   * presented several times, if it has been manipulated in
   * between. The default version of this function will extract all
   * final state particles, temporarily boost them according to the
   * transform(tEventPtr) function and call analyze(tPVector) of this
   * analysis object and those of all associated analysis objects. The
   * default version will not, however, do anything on events which
   * have not been fully generated, or have been manipulated in any
   * way.
   * @param event pointer to the Event to be analyzed.
   * @param ieve  the event number.
   * @param loop  the number of times this event has been presented.
   *              If negative the event is now fully generated.
   * @param state a number different from zero if the
   *              event has been manipulated in some way since it was last
   *              presented.
   */
  virtual void analyze(tEventPtr event, long ieve, int loop, int state);

  /**
   * Transform the event to the desired Lorentz frame and return the
   * corresponding LorentzRotation.
   * @param event a pointer to the Event to be transformed.
   * @return the LorentzRotation used in the transformation.
   * @deprecated Use transform(tcEventPtr) instead. This method is no
   *    longer used automatically.
   */
  virtual LorentzRotation transform(tEventPtr event) const;

  /**
   * Return a LorentzTransform which would put the event in the
   * desired Lorentz frame.
   * @param event a pointer to the Event to be considered.
   * @return the LorentzRotation used in the transformation.
   */
  virtual LorentzRotation transform(tcEventPtr event) const;

  /**
   * Analyze the given vector of particles. The default version calls
   * analyze(tPPtr) for each of the particles.
   * @param particles the vector of pointers to particles to be analyzed
   * @deprecated Use analyze(const tPVector &, double) instead.
   */
  virtual void analyze(const tPVector & particles);

  /**
   * Analyze the given vector of particles. The default version calls
   * analyze(tPPtr) for each of the particles.
   * @param particles the vector of pointers to particles to be analyzed
   * @param weight the weight of the current event.
   */
  virtual void analyze(const tPVector & particles, double weight);

  /**
   * Analyze the given particle.
   * @param particle pointer to the particle to be analyzed.
   * @deprecated us analyze(tPPtr, double) instead.
   */
  virtual void analyze(tPPtr particle);

  /**
   * Analyze the given particle.
   * @param particle pointer to the particle to be analyzed.
   * @param weight the weight of the current event.
   */
  virtual void analyze(tPPtr particle, double weight);

  //@}

  /** @name Functions to access histograms. */
  //@{
  /**
   * Check if the associated EventGenerator has been assigned a
   * histogram factory. If \a warn is true also emit a warning saying
   * that no histograms will be generated.
   */
  bool checkHistogramFactory(bool warn = false) const;

  /**
   * Access the HistogramFactory from the EventGenerator.
   */
  FactoryBase & histogramFactory();

  /**
   * Access the HistogramFactory from the EventGenerator.
   */
  const FactoryBase & histogramFactory() const;

  /**
   * Access the underlying AIDA::IHistogramFactory in the
   * HistogramFactory from the EventGenerator.
   */
  AIDA::IHistogramFactory & iHistogramFactory() const {
    return histogramFactory().histogramFactory();
  }

  /**
   * Normalize the histogran \a h using the collected statistics from
   * the EventGenerator. If the histogram has been filled with the
   * Event::weight() as weight a plotted function will correspond to a
   * proper cross section distribution in units of \a unit. This only
   * works for evenly binned histograms. If not evenly binned, nothing
   * will be done.
   */
  void normalize(tH1DPtr h, CrossSection unit = picobarn) const;

  /**
   * Normalize the histogran \a h to unit integral.
   */
  void unitNormalize(tH1DPtr h) const;
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
   * Standard Init function used to initialize the interface.
   */
  static void Init();

protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;

  //@}

private:

  /**
   * A list of slave analysis objects which are called for the same
   * extracted particles and in the same Lorentz frame as this one.
   */
  AnalysisVector theSlaves;
  //@}

public:

  /** Exception class used if no histogram factory was found. */
  class NoHistFactory: public InitException {};

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<AnalysisHandler> initAnalysisHandler;

};

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of AnalysisHandler. */
template <>
struct BaseClassTrait<AnalysisHandler,1>: public ClassTraitsType {
  /** Typedef of the first base class of AnalysisHandler. */
  typedef HandlerBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the AnalysisHandler class and the shared object where it is defined. */
template <>
struct ClassTraits<AnalysisHandler>: public ClassTraitsBase<AnalysisHandler> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::AnalysisHandler"; }
};

/** @endcond */

}

#endif /* ThePEG_AnalysisHandler_H */
