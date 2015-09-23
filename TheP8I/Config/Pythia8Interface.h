// -*- C++ -*-


#ifndef THEP8I_Pythia8Interface_H
#define THEP8I_Pythia8Interface_H
//
// This is the declaration of the Pythia8Interface class.
//

#include "ThePEG/Utilities/Throw.h"
#include "TheP8I/Config/TheP8I.h"
#include "TheP8I/Config/RndmEngine.h"
#include "ThePEG/Utilities/ObjectIndexer.h"
#include "RopeUserHooks.h"
#include "Pythia.h"
namespace TheP8I {

/**
 * The Pythia8Interface class is a wrapper around a static Pythia8
 * object. All communication between ThePEG and Pythia8 is handled by
 * static functions of the Pythia8Interface class.
 */
class Pythia8Interface: public Base {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  Pythia8Interface();

  /**
   * The destructor.
   */
  ~Pythia8Interface();
  //@}

public:

  /**
   * Initialize the main Pythia8 object setting all parameters which
   * are not default by copying them from the corresponding
   * interfaces..
   */
  void init(const Interfaced &, const vector<string> &);

  /**
   * Set parameters which are not default in the main Pythia8 object
   * by copying them from the corresponding interfaces..
   */
  void setParameters(const Interfaced &, const vector<string> &);

  /**
   * Access the local Pythia8 object.
   */
  Pythia8::Pythia & operator()() {
    return *pythia;
  }

  /**
   * Access the event object of the local Pythia8 object.
   */
  Pythia8::Event & event() {
    return pythia->event;
  }

  /**
   * Check if the local Pythia8 object has been created.
   */
  bool created() const {
    return pythia != 0;
  }

  /**
   * Prepare for introducing a new event.
   */
  void clearEvent() {
    event().reset();
    colourIndex.clear();
    colourIndex(0, tColinePtr());
    particleIndex.clear();
  }

  /**
   * Add a ThePEG::Particle to the Pythia8 event (if it hasn't been
   * added already) and return its index.
   */
  int addParticle(tPPtr p, int status, int mother1, int mother2);

  /**
   * Add a ThePEG::Particle to the Pythia8 event (if it hasn't been
   * added already) and return its index.
   */
  int addParticle(tcPPtr p, int status, int mother1, int mother2) {
    return addParticle(const_ptr_cast<tPPtr>(p), status, mother1, mother2);
  }

  /**
   * Get the colour Pythia8 colour index for a given colour
   * line. Generate a new index if the line has not been seen before.
   */
  int addColourLine(tColinePtr c);

  /**
   * Return index of the given particle \a p.
   * @return -1 if not found.
   */
  int indexOf(tPPtr p) {
    return particleIndex.included(p)? particleIndex(p): -1;
  }

  /**
   * Return a particle corresponding to the given entry in the Pythia8 event.
   */
  PPtr getParticle(int idx);

  /**
   * After the event has been filled, let Pythia8 do its thing. Return
   * true if everything was OK.
   */
  bool go();

  /*
  * Get a pointer to the UserHooks object (for ropes), pointer is NULL if not initialized
  */
  RopeUserHooks * getRopeUserHooksPtr(){
    return hooks;
  }

  void enableHooks(){
    doHooks = true;
  }
  
  void errorlist();
  /**
   * Retrieve the version number of Pythia.
   */
  double version() const {
    return pythia? pythia->settings.parm("Pythia:versionNumber"): -1.0;
  }

  /**
   * Print out debugging information.
   */
  void debug();

  /**
   * Perform all static initializations.
   */
  static void Init();

private:

  /**
   * The main Pythia8 object.
   */
  Pythia8::Pythia * pythia;

  /**
   * Pythia 8 user hooks object for fooling around with ropes.
   */  
  RopeUserHooks * hooks;

  bool doHooks;
  /**
   * Association between ColourLines and colour indices.
   */
  ObjectIndexer<int,ColourLine> colourIndex;

  /**
   * Association between Particles and indices.
   */
  ObjectIndexer<int,Particle> particleIndex;

  /**
   * True if the static members have been initialized.
   */
  static bool initialized;

  /**
   * The place where Pythia 8 is installed.
   */
  static string xmlDir;

  /**
   * The interface to the random generatorof ThePEG.
   */
  static RndmEngine rnd;

  /**
   * Internal Exception class.
   */
  struct Pythia8InitException: public Exception {};

  /**
   * Internal Exception class.
   */
  struct Pythia8ExecException: public Exception {};

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  Pythia8Interface & operator=(const Pythia8Interface &);

};

}


#endif /* THEP8I_Pythia8Interface_H */
