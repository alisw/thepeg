// -*- C++ -*-
//
// HandlerGroup.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_HandlerGroup_H
#define ThePEG_HandlerGroup_H
// This is the declaration of the HandlerGroup class.

#include "ThePEG/Config/ThePEG.h"
// #include "HandlerGroup.fh"
// #include "HandlerGroup.xh"

namespace ThePEG {

/**
 * HandlerGroupBase is the base class for the templated HandlerGroup
 * utility class to manage a group of <code>StepHandler</code>s.
 *
 * The derived StepHandler has a main StepHandler (CascadeHandler,
 * MultipleInteractionHandler, HadronizationHandler or DecayHandler)
 * while this bease class has a list of pre-hadlers and a list of
 * post-handlers.
 *
 * The <code>HandlerGroup</code> class is used in the
 * EventHandler and SubProcessHandler to manage the
 * post-sub-process handler, the cascade, multiple interaction,
 * hadronization and decay handler groups. When an event is generated,
 * after the main sub-process is performed, all handler groups are
 * processed in turn. In each group the pre-hadnlers are run first,
 * followed by the main handler (which may be run several times is
 * more than one Hint has been specified) and finally the
 * post-handlers are run.
 *
 * When a group is initialised before each run, an auxilliary
 * HandlerGroupBase object may be specified to override the default
 * handlers in this group.
 *
 * @see HandlerGroup
 */
class HandlerGroupBase {

public:

  /** Associate a StepHandler with a Hint object. */
  typedef pair<StepHdlPtr, HintPtr> StepWithHint;

  /** A vector of StepHandler objects. */
  typedef vector<StepHdlPtr> StepVector;

  /** A vector of StepHandler objects associated with Hint objects. */
  typedef vector<StepWithHint> StepHintVector;

  /** A vector of Hint objects. */
  typedef deque<HintPtr> HintVector;

public:

  /**
   * Default constructor.
   */
  HandlerGroupBase();

  /**
   * Destructor.
   */
  virtual ~HandlerGroupBase();

  /**
   * Returns true if current selections in this group is empty.
   */
  bool empty() const { return isEmpty; }

  /**
   * Initialize, taking the default StepHandlers as the current ones,
   * possibly overridden by the default ones in the auxilliary group
   * supplied in the argument.
   */
  void init(const HandlerGroupBase & ext) {
    clear();
    refillDefaults(ext);
  }

  /**
   * Return the next step;
   */
  StepWithHint next();

  /**
   * Add a step handler, \a sh to the current list of
   * pre-handlers. Optionally a \a hint may be specified. If the main
   * handler has already been executed, the object is reinitialized
   * using \a ext to override defaults.
   */
  void addPreHandler(tStepHdlPtr sh, tHintPtr hint,
		     const HandlerGroupBase & ext);

  /**
   * Add a step handler, \a sh, to the current list of
   * post-handlers. Optionally a \a hint may be specified. If the main
   * handler has already been executed, the object is reinitialized
   * using \a ext to override defaults.
   */
  void addPostHandler(tStepHdlPtr sh, tHintPtr hint,
		      const HandlerGroupBase &);

  /**
   * Add a \a hint to the currently selected main handler. If the main
   * handler has already been executed, the object is reinitialized
   * using \a ext to override defaults.
   */
  void addHint(tHintPtr hint, const HandlerGroupBase & ext);

  /**
   * Return a reference to the list of default pre-handlers.
   */
  StepVector & preHandlers() { return theDefaultPreHandlers; }

  /**
   * Return a reference to the list of default pre-handlers.
   */
  const StepVector & preHandlers() const { return theDefaultPreHandlers; }

  /**
   * Return a pointer to the default main handler.
   */
  virtual tStepHdlPtr defaultHandler() const = 0;

  /**
   * Return a reference to the list of default post-handlers.
   */
  StepVector & postHandlers() { return theDefaultPostHandlers; }

  /**
   * Return a reference to the list of default post-handlers.
   */
  const StepVector & postHandlers() const { return theDefaultPostHandlers; }

  /**
   * Return a pointer to the current main handler.
   */
  virtual tStepHdlPtr handler() const = 0;

  /**
   * Unset the current main handler.
   */
  virtual void setHandler() = 0;

  /**
   * Set the current main handler, but also refill the current pre-
   * and post- handlers with the defaults from \a ext.
   */
  virtual bool setHandler(tStepHdlPtr, const HandlerGroupBase & ext) = 0;

  /**
   * Set the current main handler. If the null pointer use the default
   * main handler.
   */
  virtual void refillDefaultHandler(tStepHdlPtr) = 0;

  /**
   * Fill main, pre- and post- handlers with the default ones. The
   * default handlers in the argument takes precedence to this.
   */
  void refillDefaults(const HandlerGroupBase &);

  /**
   * Clear all current handlers, but don't touch the default ones.
   */
  virtual void clear();

  /**
   * Return the base class name of the main handler type.
   */
  virtual string handlerClass() const = 0;

  /**
   * Utility function used for the interface.
   */
  void interfaceSetPrehandler(StepHdlPtr p, int i);

  /**
   * Utility function used for the interface.
   */
  void interfaceInsertPrehandler(StepHdlPtr p, int i);

  /**
   * Utility function used for the interface.
   */
  void interfaceErasePrehandler(int i);

  /**
   * Utility function used for the interface.
   */
  vector<StepHdlPtr> interfaceGetPrehandlers() const;

  /**
   * Utility function used for the interface.
   */
  void interfaceSetPosthandler(StepHdlPtr p, int i);

  /**
   * Utility function used for the interface.
   */
  void interfaceInsertPosthandler(StepHdlPtr p, int i);

  /**
   * Utility function used for the interface.
   */
  void interfaceErasePosthandler(int i);

  /**
   * Utility function used for the interface.
   */
  vector<StepHdlPtr> interfaceGetPosthandlers() const;

  /**
   * Write to persistent streams.
   */
  virtual void write(PersistentOStream &) const;

  /**
   * Read from persistent streams.
   */
  virtual void read(PersistentIStream &);

protected:

  /**
   * The copy constructor is only used via subclasses.
   */
  HandlerGroupBase(const HandlerGroupBase &);

  /**
   * True if the current handlers are empty.
   */
  bool isEmpty;

private:

  /**
   * Add handlers from the def vector to the current, supplying them
   * with default hints.
   */
  void checkInsert(StepHintVector & current, const StepVector & def);

protected:

  /**
   * The default pre-handlers with hints.
   */
  StepVector theDefaultPreHandlers;

  /**
   * The default post-handlers with hints.
   */
  StepVector theDefaultPostHandlers;

  /**
   * The current pre-handlers with hints.
   */
  StepHintVector thePreHandlers;

  /**
   * The current hints for the main handler.
   */
  HintVector theHints;

  /**
   * The current post-handlers with hints.
   */
  StepHintVector thePostHandlers;

private:

  /**
   * Assignment is private.
   */
  HandlerGroupBase & operator=(const HandlerGroupBase &);
  
};

/**
 * HandlerGroup is a templated utility class to manage a
 * group of <code>StepHandler</code>s. All HandlerGroup
 * classes are derived from the <code>HandlerGroupBase</code> class. As
 * an example the specialization
 * <code>HandlerGroup<CascadeHandler></code> keeps a
 * CascadeHandler object and associated pre- and
 * post- StepHandlers, defining shich steps should be
 * performed before the perturbative cascade, which object should be
 * used for the cascade and which steps should be performed after.
 *
 * The <code>HandlerGroup</code> keesp both a default main handler and
 * the corresponding default pre- and post- handlers as well as the
 * main handler and pre/post hadlers chosen for the current event. The
 * current handlers are accompanied by Hints. Handlers which are
 * copied from the default ones are accompanied by the default Hint,
 * while handlers supplied from the outside may be accompanied by any
 * kind of hint. The main handler can be supplied with several hints,
 * the pre- and post- handlers may only have one hint each.
 *
 * The <code>HandlerGroup</code> class is used in the
 * EventHandler and SubProcessHandler to manage the
 * post-sub-process handler, the cascade, multiple interaction,
 * hadronization and decay handler groups.
 * 
 * @see EventHandler
 * @see SubProcessHandler
 * @see StepHandler
 * @see CascadeHandler
 * @see MultipleInteractionHandler
 * @see HadronizationHandler
 * @see DecayHandler
 * 
 */
template <typename HDLR>
class HandlerGroup: public HandlerGroupBase {

public:

  /** A pointer to the template argument class. */
  typedef typename Ptr<HDLR>::pointer HdlPtr;

  /** A transient pointer to the template argument class. */
  typedef typename Ptr<HDLR>::transient_pointer tHdlPtr;

public:

  /**
   * Destructor.
   */
  virtual ~HandlerGroup();

  /**
   * Set the current main handler. Also refill the current pre- and
   * post- handlers with the defaults from \a ext.
   */
  virtual bool setHandler(tStepHdlPtr, const HandlerGroupBase & ext);

  /**
   * Unset the current main handler.
   */
  virtual void setHandler() { theHandler = HdlPtr(); }

  /**
   * Return a pointer to the current main handler.
   */
  virtual tStepHdlPtr handler() const {
    return dynamic_ptr_cast<tStepHdlPtr>(theHandler);
  }

  /**
   * Return a pointer to the default main handler.
   */
  virtual tStepHdlPtr defaultHandler() const {
    return dynamic_ptr_cast<tStepHdlPtr>(theDefaultHandler);
  }

  /**
   * Set the current main handler. If the null pointer use the default
   * main handler.
   */
  virtual void refillDefaultHandler(tStepHdlPtr);

  /**
   * Clear all current handlers, but don't touch the default ones.
   */
  virtual void clear();

  /**
   * Return the base class name of the main handler type.
   */
  virtual string handlerClass() const;

  /**
   * Utility function used for the interface.
   */
  void interfaceSetHandler(HdlPtr);

  /**
   * Utility function used for the interface.
   */
  HdlPtr interfaceGetHandler() const;

  /**
   * Write to persistent streams.
   */
  virtual void write(PersistentOStream & os) const {
    os << theDefaultHandler << theHandler;
    HandlerGroupBase::write(os);
  }

  /**
   * Read from persistent streams.
   */
  virtual void read(PersistentIStream & is) {
    is >> theDefaultHandler >> theHandler;
    HandlerGroupBase::read(is);
  }

private:


  /**
   * The default main handler.
   */
  HdlPtr theDefaultHandler;

  /**
   * The current main handler.
   */
  HdlPtr theHandler;

private:

  /**
   * Assignment is private.
   */
  HandlerGroup<HDLR> & operator=(const HandlerGroup<HDLR> &);
  
};

/** Namespace to encapsulate enums related to <code>HandlerGroup</code>s. */
namespace Group {

/**
 * Enumeration for the type of <code>HandlerGroup</code>s.
 */
enum Handler {
  subproc, /**< The sub-process group. */
  cascade, /**< The CascadeHandler group. */
  multi,   /**< The MultipleInteractionHandler group. */
  hadron,  /**< The HadronizationHandler group. */
  decay    /**< The DecayHandler group. */
};

/** Enumeration for the type of step handler */
enum Level {
  before, /**< A pre-handler. */
  main,   /**< The mainhandler. */
  after   /**< A post-handler. */
};
}

/** Output a HandlerGroup to a PersistentOStream. */
template <typename HDLR>
inline PersistentOStream & operator<<(PersistentOStream & os,
				      const HandlerGroup<HDLR> & hg) {
  hg.write(os);
  return os;
}

/** Input a HandlerGroup from a PersistentIStream. */
template <typename HDLR>
inline PersistentIStream & operator>>(PersistentIStream & is,
				      HandlerGroup<HDLR> & hg) {
  hg.read(is);
  return is;
}

}

/** Macro for declaring a prepost group */
#define ThePEG_DECLARE_PREPOST_GROUP(HandlerClass,prepost)                    \
/** Utility function for the interface. */                                    \
void interfaceSet##prepost##HandlerClass(StepHdlPtr, int);                    \
/** Utility function for the interface. */                                    \
void interfaceInsert##prepost##HandlerClass(StepHdlPtr, int);                 \
/** Utility function for the interface. */                                    \
void interfaceErase##prepost##HandlerClass(int);                              \
/** Utility function for the interface. */                                    \
vector<StepHdlPtr> interfaceGet##prepost##HandlerClass() const

/** Macro for declaring a group interface */
#define ThePEG_DECLARE_GROUPINTERFACE(HandlerClass,ptr)                       \
ThePEG_DECLARE_PREPOST_GROUP(HandlerClass,Pre);                               \
/** Utility function for the interface. */                                    \
void interfaceSet##HandlerClass(ptr);                                         \
/** Utility function for the interface. */                                    \
ptr interfaceGet##HandlerClass() const;                                       \
ThePEG_DECLARE_PREPOST_GROUP(HandlerClass,Post)

/** Macro for implementing a prepost group. */
#define ThePEG_IMPLEMENT_PREPOST_GROUP(ThisClass,HandlerClass,member,pp)      \
void ThisClass::interfaceSet##pp##HandlerClass(StepHdlPtr p , int i) {     \
  member.interfaceSet##pp##handler(p,i);                                     \
}                                                                              \
void ThisClass::interfaceInsert##pp##HandlerClass(StepHdlPtr p, int i) {   \
  member.interfaceInsert##pp##handler(p,i);                                  \
}                                                                              \
void ThisClass::interfaceErase##pp##HandlerClass(int i) {                  \
  member.interfaceErase##pp##handler(i);                                     \
}                                                                              \
vector<StepHdlPtr> ThisClass::interfaceGet##pp##HandlerClass() const {     \
  return member.interfaceGet##pp##handlers();                                  \
}

/** Macro for implementing a group interface. */
#define ThePEG_IMPLEMENT_GROUPINTERFACE(ThisClass,HandlerClass,member,ptr)     \
ThePEG_IMPLEMENT_PREPOST_GROUP(ThisClass,HandlerClass,member,Pre)              \
void ThisClass::interfaceSet##HandlerClass(ptr p) {                            \
  member.interfaceSetHandler(p);                                               \
}                                                                              \
ptr ThisClass::interfaceGet##HandlerClass() const {                            \
  return member.interfaceGetHandler();                                         \
}                                                                              \
ThePEG_IMPLEMENT_PREPOST_GROUP(ThisClass,HandlerClass,member,Post)             \

/** Macro for declaring prepost objects. */
#define ThePEG_DECLARE_PREPOST_OBJECTS(ThisClass,HandlerClass,pp,ba)           \
static RefVector<ThisClass,StepHandler> interface##pp##HandlerClass            \
(#pp #HandlerClass "s",                                                        \
 "A list of handlers to be called " #ba " the " #HandlerClass ". "             \
 "If handler objects are specified in a EventHandler and "                     \
 "the SubProcessHandler chosen in a given collision also specifies some, "     \
 "the latter will caled first.",                                               \
 0, 0, false, false, true, false,                                              \
 &ThisClass::interfaceSet##pp##HandlerClass,                                   \
 &ThisClass::interfaceInsert##pp##HandlerClass,                                \
 &ThisClass::interfaceErase##pp##HandlerClass,                                 \
 &ThisClass::interfaceGet##pp##HandlerClass)

/** Macro for declaring group interface objects. */
#define ThePEG_DECLARE_GROUPINTERFACE_OBJECTS(ThisClass,HandlerClass)          \
ThePEG_DECLARE_PREPOST_OBJECTS(ThisClass,HandlerClass,Pre, before);            \
static Reference<ThisClass,HandlerClass> interface ## HandlerClass             \
(#HandlerClass,                                                                \
 "The " #HandlerClass " object used in this " #ThisClass ". "                  \
 "If a " #HandlerClass " object is specified in a EventHandler and "           \
 "the SubProcessHandler chosen in a given collision also specifies one,"       \
 "the latter will be used.",                                                   \
 0, false, false, true, true,                                                  \
 &ThisClass::interfaceSet ## HandlerClass,                                     \
 &ThisClass::interfaceGet ## HandlerClass);                                    \
ThePEG_DECLARE_PREPOST_OBJECTS(ThisClass,HandlerClass,Post, after)

#ifndef ThePEG_TEMPLATES_IN_CC_FILE
#include "HandlerGroup.tcc"
#endif

#endif /* ThePEG_HandlerGroup_H */
