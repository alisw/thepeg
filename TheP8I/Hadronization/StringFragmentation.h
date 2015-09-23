// -*- C++ -*-

#ifndef THEP8I_StringFragmentation_H
#define THEP8I_StringFragmentation_H
//
// This is the declaration of the StringFragmentation class.
//

#include "ThePEG/Handlers/HadronizationHandler.h"
#include "ThePEG/Handlers/ClusterCollapser.h"
#include "TheP8I/Config/Pythia8Interface.h"
#include "OverlapPythiaHandler.h"
#include "TheP8EventShapes.h"
#include <sstream>
//#include "YODA/Histo1D.h"
//#include "YODA/Histo2D.h"
//#include "YODA/Writer.h"
//#include "YODA/WriterYODA.h"

namespace TheP8I {

      bool sorting_principle(pair<ThePEG::ColourSinglet*, vector<TheP8I::Ropewalk::Dipole*> > c1, pair<ThePEG::ColourSinglet*, vector<TheP8I::Ropewalk::Dipole*> > c2){
	return (c1.first->momentum().perp() < c2.first->momentum().perp());
  }


using namespace ThePEG;

/**
 * Here is the documentation of the StringFragmentation class.
 *
 * @see \ref StringFragmentationInterfaces "The interfaces"
 * defined for StringFragmentation.
 */
class StringFragmentation: public HadronizationHandler {

typedef map<string, double> PytPars;

friend class Ropewalk::Dipole;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  StringFragmentation();

  /**
   * The destructor.
   */
  virtual ~StringFragmentation();
  //@}

public:

  /**
    * The main function called by the EventHandler class to
    * perform the Hadronization step.
    * @param eh the EventHandler in charge of the Event generation.
    * @param tagged if not empty these are the only particles which should
    * be considered by the StepHandler.
    * @param hint a Hint object with possible information from previously
    * performed steps.
    * @throws Veto if the StepHandler requires the current step to be
    * discarded.
    * @throws Stop if the generation of the current Event should be stopped
    * after this call.
    * @throws Exception if something goes wrong.
    */


  virtual void handle(EventHandler & eh, const tPVector & tagged,
		      const Hint & hint);

  /** @name Functions used by the persistent I/O system. */
  //@{

   /**
   * Let the given Pythia8Interface hadronize the ColourSinglet
   * systems provided.
   */
  bool hadronizeSystems(Pythia8Interface & pyt, const vector<ColourSinglet> & singlets,
                        const tPVector & all);
  

  void hadronizeTweak(Pythia8Interface & pyt, const ColourSinglet & singlet,
                        const tPVector & all);
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;

  /**
   * Helper function to get the maximum transverse momenta of any
   * parton in a string. Optionally only consider partons within a
   * rapidity interval |eta| < deltaY.
   */
  static Energy maxPT(const ColourSinglet & cs, double deltaY = 0.0);

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
  //@}

private:

  /**
   * The interface to the Pythia 8 object,
   */
  Pythia8Interface pythia;

   /**
   * Internal switch for the fragmentation scheme 
   */
  int fScheme;

  /**
  * Object that handles calculation of new Pythia parameters
  */
  ParameterHandler phandler;

  /**
   * The intrinsic string radius
   */
  Length stringR0;

   /**
   * The assumed mass for calculating rapidities in the reeperbahn scheme
   */
  Energy stringm0;

 /**
   * Supression of diquarks emerging from breaking junctions
   */
  double junctionDiquark;

  double alpha;

  bool average;
  
  
  /**
   * If **window** is set, this determines the abs(y) window for which no enhancement will be made
   */
  double stringyCut;
 
  /**
   * If **window** is set, this determines the pT window for which no enhancement will be made
   */
  Energy stringpTCut;

  /**
   * Assumed m0 value in calculation of normalization of f(z) the Lund splitting function.
   */
  double fragmentationMass;

  /**
   * Parameter surpressing baryonic content further in calculation of effective parameters
   */
  double baryonSuppression;

  // Force hadronize systems to run to the end, regardless of windows
  bool forcerun;

  /**
   * Can provide a window in y and pT where no enhancement is made, in order to not include
   * very hard gluons in central rapidity region in a rope
   */
  bool window;

  /**
   * Throw away some string entirely while making ropes 
   */
  bool throwaway;

  /**
   * Do rap. overlap
   **/
  bool rapidityOverlap;
  
  TheP8EventShapes * _eventShapes;
  /**
   * Use thrust axis to calculate rapidities; useful for LEP validation
   */
  int useThrustAxis;
 
  /**
   * Additional interfaces to Pythia8 objects in case the string
   * overlap model is to be used.
   */
  OverlapPythiaHandler * opHandler;

  /**
   * Sometimes Pythia gives up on an event too easily. We allow it to
   * re-try a couple of times.
   */
  int maxTries;

  // Keep track of how many events we already hadronized using this object, in order to clean up
  int nev;

  // Can do a small analysis
  int doAnalysis;
  //vector<YODA::Histo1D*> _histograms;
  //vector<YODA::Histo2D*> _histograms2D;
  string analysisPath;

  void dofinish();

  // Can print the strings of the event in a matlab friendly format
  string PrintStringsToMatlab(vector<ColourSinglet>& singlets);

  string convert(double d);
#include "StringFragmentation-var.h"

  /**
   * The object used to avoid too small strings in the hadronization.
   */
  ClusterCollapserPtr theCollapser;

public:

  /**
   * Exception class declaration.
   */
  struct StringFragError: public Exception {};

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<StringFragmentation> initStringFragmentation;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  StringFragmentation & operator=(const StringFragmentation &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of StringFragmentation. */
template <>
struct BaseClassTrait<TheP8I::StringFragmentation,1> {
  /** Typedef of the first base class of StringFragmentation. */
  typedef HadronizationHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the StringFragmentation class and the shared object where it is defined. */
template <>
struct ClassTraits<TheP8I::StringFragmentation>
  : public ClassTraitsBase<TheP8I::StringFragmentation> {
  /** Return a platform-independent class name */
  static string className() { return "TheP8I::StringFragmentation"; }
  /**
   * The name of a file containing the dynamic library where the class
   * StringFragmentation is implemented. It may also include several, space-separated,
   * libraries if the class StringFragmentation depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "libTheP8I.so"; }
};

/** @endcond */

}

#endif /* THEP8I_StringFragmentation_H */
