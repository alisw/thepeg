// -*- C++ -*-
//
// EventInfoBase.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_EventInfoBase_H
#define ThePEG_EventInfoBase_H
// This is the declaration of the EventInfoBase class.

#include "ThePEG/EventRecord/EventConfig.h"
#include "ThePEG/Utilities/ClassDescription.h"

namespace ThePEG {

/**
 * EventInfoBase is a base class for information objects. It is used
 * as a base class for classes representing user-defined information
 * which may be associated with a Particle. The class itself is
 * practically empty. Information added in sub-classes can be accessed
 * from a Particle by the Particle::getInfo() function and the
 * resulting pointers need to be dynamically cast to check if they are
 * of a desired class.
 */
class EventInfoBase: public EventRecordBase {

public:

  /**
   * Rebind to cloned objects. If an EventInfoBase is cloned together
   * with a whole Event and this has pointers to other event record
   * objects, these should be rebound to their clones in this
   * function.
   */
  virtual void rebind(const EventTranslationMap & ) {}

  /**
   * Standard Init function. @see Base::Init().
   */
  static void Init() {}

  /**
   * Standard clone method.
   */
  virtual EIPtr clone() const { return new_ptr(*this); }

private:

  /**
   * Describe concrete class without persistent data.
   */
  static NoPIOClassDescription<EventInfoBase> initEventInfoBase;

  /**
   *  Private and non-existent assignment operator.
   */
  EventInfoBase & operator=(const EventInfoBase &) = delete;

};


/** @cond TRAITSPECIALIZATIONS */
ThePEG_DECLARE_CLASS_TRAITS(EventInfoBase,EventRecordBase);
/** @endcond */

}

#endif /* ThePEG_EventInfoBase_H */
