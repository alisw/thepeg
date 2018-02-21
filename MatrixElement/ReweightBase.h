// -*- C++ -*-
//
// ReweightBase.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_ReweightBase_H
#define ThePEG_ReweightBase_H
// This is the declaration of the ReweightBase class.

#include "ThePEG/Handlers/HandlerBase.h"
#include "ThePEG/EventRecord/SubProcess.h"
#include "ThePEG/Handlers/LastXCombInfo.h"
#include "ThePEG/Handlers/StandardXComb.fh"

namespace ThePEG {

/**
 * The ReweightBase class is the base class of all objects
 * representing external biases to matrix elements. These can be used
 * to enhance certain matrix elements or certain phase space
 * regions. They can be used in two ways, either to completely change
 * the matrix element (re-weight), in which case the total cross
 * section will be affected or, when using weighted events in an
 * EventHandler, to pre-weight certain events but leaving the cross
 * section unchanged
 *
 * There is only one virtual function which must be overridden in
 * derived classes: weight().
 *
 * @see \ref ReweightBaseInterfaces "The interfaces"
 * defined for ReweightBase.
 * @see MEBase
 * @see EventHandler
 * @see SubProcessHandler
 */
class ReweightBase: public HandlerBase, public LastXCombInfo<> {

public:
 
  /** @name Standard constructors and destructors. */
  //@{
  
  /**
   * Destructor.
   */
  virtual ~ReweightBase();
  //@}
  
public:

  /**
   * Return the wieght for the kinematical configuation provided by
   * the assigned XComb object (in the LastXCombInfo base class).
   */
  virtual double weight() const = 0;

  /**
   * Assigne an XComb object with information about the sub-process to
   * be used in the reweighting.
   */
  void setXComb(tXCombPtr xc);

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
   * Standard Init function used to initialize the interfaces.
   */
  static void Init();

private:

  /**
   * Describe an abstract base class with persistent data.
   */
  static AbstractClassDescription<ReweightBase> initReweightBase;

  /**
   *  Private and non-existent assignment operator.
   */
  ReweightBase & operator=(const ReweightBase &);

};

}


namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * This template specialization informs ThePEG about the base class of
 * ReweightBase.
 */
template <>
struct BaseClassTrait<ReweightBase,1>: public ClassTraitsType {
  /** Typedef of the base class of ReweightBase. */
  typedef HandlerBase NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * ReweightBase class.
 */
template <>
struct ClassTraits<ReweightBase>: public ClassTraitsBase<ReweightBase> {
  /** Return the class name. */
  static string className() { return "ThePEG::ReweightBase"; }
};

/** @endcond */

}

#endif /* ThePEG_ReweightBase_H */
