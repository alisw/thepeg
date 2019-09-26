// -*- C++ -*-
//
// HIHepMCFile.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_HIHepMCFile_H
#define THEPEG_HIHepMCFile_H
//
// This is the declaration of the HIHepMCFile class.
//
#include <iostream>
#include <fstream>
#include "ThePEG/Handlers/AnalysisHandler.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "HepMC/IO_BaseClass.h"

namespace ThePEG {

/** \ingroup Analysis
 * The HIHepMCFile class outputs ThePEG events in HepMC format with additional Heavy Ion information.
 *
 * @see \ref HIHepMCFileInterfaces "The interfaces"
 * defined for HIHepMCFile.
 */
class HIHepMCFile: public AnalysisHandler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  HIHepMCFile();

  /**
   * The copy constructor.
   */
  HIHepMCFile(const HIHepMCFile &);
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
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
   virtual void doinitrun();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  virtual void dofinish();
  //@}

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<HIHepMCFile> initHIHepMCFile;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  HIHepMCFile & operator=(const HIHepMCFile &) = delete;

private:

  /**
   * Last event that should be written out as HepMC format
   */
  long _eventNumber;

  /**
   *  The HepMC format
   */
  int _format;

  /**
   *  The HepMC filename
   */
  string _filename;

  /**
   *  The HepMC I/O handler
   */
  HepMC::IO_BaseClass *_hepmcio;

  /**
   *  The HepMC dump file
   */
  ofstream _hepmcdump;

  /**
   * Selector for the choice of units
   */
  int _unitchoice;

  /**
   * Choice of output precision in GenEvent format
   */
  unsigned int _geneventPrecision;
};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of HIHepMCFile. */
template <>
struct BaseClassTrait<HIHepMCFile,1> {
  /** Typedef of the first base class of HIHepMCFile. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the HIHepMCFile class and the shared object where it is defined. */
template <>
struct ClassTraits<HIHepMCFile>
  : public ClassTraitsBase<HIHepMCFile> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::HIHepMCFile"; }
  /** Return the name(s) of the shared library (or libraries) be loaded to get
   *  access to the HIHepMCFile class and any other class on which it depends
   *  (except the base class). */
  static string library() { return "HepMCAnalysis.so"; }
};

/** @endcond */

}

#endif /* THEPEG_HIHepMCFile_H */
