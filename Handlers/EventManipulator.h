// -*- C++ -*-
//
// EventManipulator.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_EventManipulator_H
#define ThePEG_EventManipulator_H
// This is the declaration of the EventManipulator class.

#include "ThePEG/Interface/Interfaced.h"
#include "ThePEG/Handlers/EventHandler.fh"
#include <stdexcept>

namespace ThePEG {

/**
 * An object of the EventManipulator class may be assigned to a
 * FullEventGenerator object. The manipulate() method is called for
 * each event generated, after the AnalysisHandlers have been called,
 * and may manipulate the event in any way needed. The manipulator may
 * alseo add StepHandlers to the EventHandler which produced the
 * event. The manipulate() method returns an integer which should be
 * zero if nothing was done to the event. If the EventHandler has
 * steps left to do, these are performed, after which the
 * <code>AnalysisHandler</code>s are called with the return value from
 * the previous manipulate() call. Then manipulate is called again and
 * the procedure is repeated until the EventHandler has no more steps
 * to do.
 *
 * @see \ref EventManipulatorInterfaces "The interfaces"
 * defined for EventManipulator.
 * @see FullEventGenerator
 * @see AnalysisHandler
 * @see EventHandler
 * @see StepHandler
 * 
 */
class EventManipulator: public Interfaced {

public:

  /**
   * Manipulate an event and the event handler.
   * @param eh the EventHandler in charge of the generation.
   * @param event the Event to be manipulated.
   * @return zero if the event was not manipulated. Otherwise return
   * an integer which will be given to the
   * <code>AnalysisHandler</code>s of the current FullEventGenerator.
   */
  virtual int manipulate(tEHPtr eh, tEventPtr event) = 0;

public:

  /**
   * Standard Init function used to initialize the interface.
   */
  static void Init();

protected:

private:

  /**
   * Describe an abstract class without persistent data.
   */
  static AbstractNoPIOClassDescription<EventManipulator> initEventManipulator;

  /**
   *  Private and non-existent assignment operator.
   */
  EventManipulator & operator=(const EventManipulator &) = delete;

};

/** @cond TRAITSPECIALIZATIONS */

/**
 * This template specialization informs ThePEG about the
 * base class of EventManipulator.
 */
template <>
struct BaseClassTrait<EventManipulator,1>: public ClassTraitsType {
  /** Typedef of the base class of EventManipulator. */
  typedef Interfaced NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * DecayHandler class and the shared object where it is defined.
 */
template <>
struct ClassTraits<EventManipulator>:
    public ClassTraitsBase<EventManipulator> {
  /** Return the class name. */
  static string className() {  return "ThePEG::EventManipulator"; }
};

/** @endcond */

}

#endif /* ThePEG_EventManipulator_H */
