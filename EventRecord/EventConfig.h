// -*- C++ -*-
//
// EventConfig.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_EventConfig_H
#define ThePEG_EventConfig_H

/** \file EventConfig.h

 * This is the main config header file for the Event classes. Do not
 * make changes in this file. If these classes are used outside of
 * ThePEG you probably need to need to modify base classes persistency
 * scheme etc. If so, edit a copy of the file which can be included
 * instead of this file using the macro
 * <code>ThePEG_ALTERNATIVE_EVENT_RECORD</code>.
 */

#ifndef ThePEG_NOT_ThePEG

#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/Utilities/Rebinder.fh"
#include "ThePEG/Persistency/PersistentOStream.fh"
#include "ThePEG/Persistency/PersistentIStream.fh"

#ifndef ThePEG_ALTERNATIVE_EVENT_RECORD

#include "ThePEG/Vectors/Lorentz5Vector.h"
#include "ThePEG/Vectors/LorentzRotation.h"

namespace ThePEG {

/** EventRecordBase is the base class of all event record classes. It
 *  must be typedefed to a class which allows for garbage
 *  collection, since none of the event record classes are deleted
 *  explicitly.
 */
typedef Base EventRecordBase;
/** Alias for a reference counted pointer to EventRecordBase. */
typedef Ptr<EventRecordBase>::pointer EventBasePtr;
/** Alias for a reference counted pointer to const EventRecordBase. */
typedef Ptr<EventRecordBase>::const_pointer cEventBasePtr;
/** Alias for a transient pointer to EventRecordBase. */
typedef Ptr<EventRecordBase>::transient_pointer tEventBasePtr;
/** Alias for a transient pointer to const EventRecordBase. */
typedef Ptr<EventRecordBase>::transient_const_pointer tcEventBasePtr;

/** Alias for a rebinder object able to relate pointers to original
 *  objects to pointers to their clones. */
typedef Rebinder<EventRecordBase> EventTranslationMap;


/** ParticleClass is the name used for Particle in the event record classes. */
typedef Particle ParticleClass;
/** ParticleDataClass is the name used for ParticleData in the event
 *  record classes. */
typedef ParticleData ParticleDataClass;

/** Alias for a reference counted pointer to ParticleDataClass. */
typedef Ptr<ParticleDataClass>::pointer EventPDPtr;
/** Alias for a reference counted pointer to const ParticleDataClass. */
typedef Ptr<ParticleDataClass>::const_pointer cEventPDPtr;
/** Alias for a transient pointer to ParticleDataClass. */
typedef Ptr<ParticleDataClass>::transient_pointer tEventPDPtr;
/** Alias for a transient pointer to const ParticleDataClass. */
typedef Ptr<ParticleDataClass>::transient_const_pointer tcEventPDPtr;

/** A vector of transient pointers to Particle. */
typedef vector<tPPtr> tParticleVector;
/** A vector of pointers to Particle. */
typedef vector<PPtr> ParticleVector;
/** A set of pointers to Particle. */
typedef set<PPtr, less<PPtr> > ParticleSet;
/** A set of transient pointers to Particle. */
typedef set<tPPtr, less<tPPtr> > tParticleSet;
/** A set of transient pointers to const Particle. */
typedef set<tcPPtr, less<tcPPtr> > tcParticleSet;
/** A vector of pointers to Step. */
typedef vector<StepPtr> StepVector;
/** A vector of pointers to SubProcess. */
typedef vector<SubProPtr> SubProcessVector;
/** A vector of transient pointers to SubProcess. */
typedef vector<tSubProPtr> tSubProcessVector;
/** A vector of pointers to Collision. */
typedef vector<CollPtr> CollisionVector;
/** A set of pointers to Step. */
typedef set<StepPtr, less<StepPtr> > StepSet;
/** A set of pointers to SubProcess. */
typedef set<SubProPtr, less<SubProPtr> > SubProcessSet;

/** A helper class to facilitate persistent input and output. */
struct EventConfig {

  /**
   * Optional pointer to the current EventGenerator.
   * If \a currentGenerator is set during persistent output, only the
   * PDG number of a particle type is written rather than the full
   * ParticleData object. Also only the name of handlers is written
   * rather than the full objects. When this is read back in again,
   * the \a currentGenerator must be set so that conversion from
   * name/number back to objects can be done.
   */
  static tcEventBasePtr currentGenerator;

  /** Write a handler object to a persistent stream. */
  static void putHandler(PersistentOStream & os, tcEventBasePtr h);
  /** Read a handler object from a persistent stream. */
  static void getHandler(PersistentIStream & is, tcEventBasePtr & h);
  /** Write a ParticleData object to a persistent stream. */
  static void putParticleData(PersistentOStream & os, tcEventPDPtr pd);
  /** Read a ParticleData object from a persistent stream. */
  static void getParticleData(PersistentIStream & is, cEventPDPtr & pd);
  /** Return the name of a handler object. */
  static string nameHandler(tcEventBasePtr h);

};

}

#else

#include ThePEG_ALTERNATIVE_EVENT_RECORD

#endif

#endif /* ThePEG_NOT_ThePEG */


#endif /* ThePEG_EventConfig_H */

