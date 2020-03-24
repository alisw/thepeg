// -*- C++ -*-
#ifndef ThePEG_NLOHepMCFile_H
#define ThePEG_NLOHepMCFile_H
//
// This is the declaration of the NLOHepMCFile class.
//
#include <config.h>
#include "ThePEG/Handlers/AnalysisHandler.h"
#include "ThePEG/Config/HepMCHelper.h"
namespace ThePEG {

/**
 * \author Simon Platzer
 *
 * \brief Write hard sub processes or sub process groups to HepMC.
 *
 */
class NLOHepMCFile: public AnalysisHandler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  NLOHepMCFile();

  /**
   * The copy constructor.
   */
  NLOHepMCFile(const NLOHepMCFile &);

  /**
   * The destructor.
   */
  virtual ~NLOHepMCFile() {}
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

  /**
   * Produca a HepMC event for the given subprocess
   */
  HepMC::GenEvent * makeEvent(tEventPtr event, tSubProPtr sub, long no,
			      Energy eUnit, Length lUnit,
			      CrossSection xsec, CrossSection xsecErr) const;

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


// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).

private:

  /**
   * The PDG ID to be used for remnants
   */
  long _remnantId;

  /**
   *  The HepMC format
   */
  int _format;

  /**
   *  The HepMC filename
   */
  string _filename;

#ifdef HAVE_HEPMC_ROOTIO
  /**
   *  The name of TTRee in ROOT file
   */
  string _ttreename;

  /**
   *  The name of branch in ROOT file
   */
  string _tbranchname;

#endif
  /**
   *  The HepMC I/O handler
   */
#ifdef HAVE_HEPMC3
  HepMC::Writer *_hepmcio;
#else
  HepMC::IO_BaseClass *_hepmcio;
#endif

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

  /**
   * Count events
   */
  unsigned long _eventNumber;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  NLOHepMCFile & operator=(const NLOHepMCFile &) = delete;

};

}

#endif /* ThePEG_NLOHepMCFile_H */
