// -*- C++ -*-
#ifndef THEPEG_HZTool_H
#define THEPEG_HZTool_H
//
// This is the declaration of the HZTool class.
//

#include "ThePEG/Handlers/AnalysisHandler.h"

namespace ThePEG {

/**
 * This class wraps the HZTool fortran library. The specified
 * functions of in the HZTool library will be called for each event
 * (after the ThePEG::Event is first converted to a HepMC::GenEvent
 * and then translated to the HEPEVT fortran common block). Note that
 * only one HZTool AnalysisHandler can be used for a given
 * EventGenerator, otherwise the result is undefined.
 *
 * @see \ref HZToolInterfaces "The interfaces"
 * defined for HZTool.
 */
class HZTool: public AnalysisHandler {

public:

  /**
   * Typedef of pointer to standard HZTool analysis routine.
   */
  typedef void (*HZAnaFn)(const int &);

  /**
   * A vector of function pointers.
   */
  typedef vector<HZAnaFn> HZAnaFnVector;

  /**
   * A vector of strings.
   */
  typedef vector<string> StringVector;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline HZTool();

  /**
   * The copy constructor.
   */
  inline HZTool(const HZTool &);

  /**
   * The destructor.
   */
  virtual ~HZTool();
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
   * Transform the event to the desired Lorentz frame and return the
   * corresponding LorentzRotation.
   * @param event a pointer to the Event to be transformed.
   * @return the LorentzRotation used in the transformation.
   */
  virtual LorentzRotation transform(tEventPtr event) const;

  /**
   * Analyze the given vector of particles. The default version calls
   * analyze(tPPtr) for each of the particles.
   * @param particles the vector of pointers to particles to be analyzed
   */
  virtual void analyze(const tPVector & particles);

  /**
   * Analyze the given particle.
   * @param particle pointer to the particle to be analyzed.
   */
  virtual void analyze(tPPtr particle);
  //@}

public:

  /**
   * The filename to which the HZTool histograms are written. If empty
   * the name of the controlling EventGenerator is used instead. The
   * standard '.rz' suffix is added to the name.
   */
  inline string filename() const;

  /**
   * Add the exchanged boson in a DIS event to the HEPEVT event
   * record. For some reason HZTool needs this, rather than figuring
   * it out from the scattered lepton. Also changes the direction of
   * the beams to fit the assumptions of HZTool. If this is not a DIS
   * event nothing is modified.
   *
   * @param sub the primary sub-process.
   */
  static void addDISLepton(const SubProcess & sub);

protected:

  /**
   * Return a function pointer to the standard HZTool analysis routine
   * called \a name. Returns null if no routine is found.
   */
  static HZAnaFn getFunctionPointer(string name);

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
  inline virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  inline virtual IBPtr fullclone() const;
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
   * The filename to which the HZTool histograms are written. If empty
   * the name of the controlling EventGenerator is used instead. The
   * standard '.rz' suffix is added to the name.
   */
  string theFilename;

  /**
   * The names of the HZTool analysis routines to be called in this
   * analysis handler.
   */
  StringVector functionNames;

  /**
   * The function pointers corresponding to the functionNames.
   */
  HZAnaFnVector functions;

  /**
   * The sum of all event weights in a run.
   */
  double sumweight;

public:

  /**
   * Utility function for the interface.
   */
  void insertFunction(string name, int pos);

  /**
   * Utility function for the interface.
   */
  void setFunction(string name, int pos);

  /**
   * Utility function for the interface.
   */
  void delFunction(int pos);

  /**
   * Exception class used if non existing HZTool function is requested.
   */
  struct MissingFunction: public InterfaceException {};

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<HZTool> initHZTool;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  HZTool & operator=(const HZTool &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of HZTool. */
template <>
struct BaseClassTrait<HZTool,1> {
  /** Typedef of the first base class of HZTool. */
  typedef AnalysisHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the HZTool class and the shared object where it is defined. */
template <>
struct ClassTraits<HZTool>
  : public ClassTraitsBase<HZTool> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::HZTool"; }
  /**
   * The name of a file containing the dynamic library where the class
   * HZTool is implemented. It may also include several, space-separated,
   * libraries if the class HZTool depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "HZTool.so"; }
};

/** @endcond */

}

#include "HZTool.icc"
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
// #include "HZTool.tcc"
#endif

#endif /* THEPEG_HZTool_H */
