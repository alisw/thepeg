// -*- C++ -*-
#ifndef PYTHIA7_LundFlavourHandler_H
#define PYTHIA7_LundFlavourHandler_H
// This is the declaration of the LundFlavourHandler class.

#include "FragConfig.h"
#include "Oriented.h"
#include "LundFlavourGenerator.h"
// #include "LundFlavourHandler.fh"
// #include "LundFlavourHandler.xh"

namespace Pythia7 {


/**
 * LundFlavourHandler is an helper class to administrate flavour
 * generation procedure that requires extra handling that simple calls
 * to the LundFlavourGenerator.
 *   
 * Particularly, the LundFlavourHandler is responsible for handling
 * the popcorn procedure for Baryon production in the Lund
 * fragmentation scheme. It uses the LundFlavourGenerator of the
 * LundFragHandler for flavour generation.
 * 
 * The LundFragHandler invokes the flavour generation through the
 * generateHadron() and delegates to the LundFlavourHandler the
 * responsibility for handling the <i>popcorn</i> meson production.
 *
 * The LundFlavourHandler is not intended to be interfaced. It should
 * be created given the pointer to the interfaced LundFlavourGenerator
 * used in the LundFragHandler.
 *
 * The generateHadron() is made virtual to provide users with the
 * possibility to implement separate strategies for popcorn meson
 * production.
 *
 * @see LundFragHandler
 * @see LundFlavourGenerator
 */
class LundFlavourHandler{

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Standard constructor given a reference to a LunfFlavourGenerator
   * object.
   */
  inline LundFlavourHandler(const LFlGenPtr& );

  /**
   * Destructor.
   */
  virtual ~LundFlavourHandler();
  //@}

  /**
   * Handle the generation of new flavour according  to the simple popcorn 
   * algorithm using the Flavour Generator of the LundFragHandler
   */
  virtual PDPtr generateHadron(tcPDPtr inPD, cPDPtr& newPD);

  /**
   * Set the flavour generator used by the fragmentation handler.
   */
  inline void FlGen(const LFlGenPtr& );

  /**
   * Get the flavour generator used by the fragmentation handler.
   */
  inline LFlGenPtr FlGen() const;

  /**
   * Initialise the number of popcorn meson and the curtain quark Id
   */
  inline void initialize();

  /**
   * Return the number of remaining popcorn mesons to be produced  
   * in the current popcorn generation started in the direction 
   * of fragmentation fixed by the LundFragHandler
   */
  inline int& PopN();

  /**
   * Return the curtain quark Id, given by the flavour generator 
   * at the beginning of the popcorn generation started in the direction 
   * of fragmentation fixed by the LundFragHandler
   */
  inline long& curtainQId();

private:

  /**
   * Pointer to the flavour generator 
   */
  LFlGenPtr theFlGen;

  /**
   * number of popcorn meson in between the Baryon and the antiBaryon.  
   * for the generations started from one or another side of the string. 
   */
  pair<int, int> thePopN;

  /**
   * the curtain quark Id in the popcorn generation.
   */
  pair<long, long> theCurtainQId;

};
}

#include "LundFlavourHandler.icc"
#ifndef PYTHIA7_TEMPLATES_IN_CC_FILE
// #include "LundFlavourHandler.tcc"
#endif

#endif /* PYTHIA7_LundFlavourHandler_H */


