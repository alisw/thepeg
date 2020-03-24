// -*- C++ -*-
//
// CascadeHandler.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_CascadeHandler_H
#define ThePEG_CascadeHandler_H
// This is the declaration of the CascadeHandler class.

#include "StepHandler.h"
#include "ThePEG/Handlers/LastXCombInfo.h"
#include "ThePEG/PDF/PDF.h"

namespace ThePEG {


/**
 * The CascadeHandler is the base class of all handlers implementing
 * perturbative partonic cascade models. It is derived from the more
 * general StepHandler class, and implements the handle() function to
 * do some standard initialization before calling the main cascade()
 * function.
 *
 * @see \ref CascadeHandlerInterfaces "The interfaces"
 * defined for CascadeHandler.
 * @see StepHandler
 * @see EventHandler
 * @see SubProcessHandler
 */
class CascadeHandler: public StepHandler, public LastXCombInfo<> {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The destructor.
   */
  virtual ~CascadeHandler();
  //@}

public:

  /** @name Virtual functions required by the StepHandler class. */
  //@{
  /**
    * The main function called by the EventHandler class to
    * perform a step.
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
  //@}

  /**
   * The main function to be overwritten by sub-classes. It is called
   * by handle() after storing some information which is then
   * available through simple access functions.
   */
  virtual void cascade() = 0;

  /**
   * The CascadeHandler can be used inside the process generation to
   * do so-called CKKW reweighting of the hard sub-process. In this
   * case this function is called after information about the
   * sub-process is made available through the LastXCombInfo base
   * class. Only the function belonging to the primary CascadeHandler
   * for the event to be generated is called. Sub-classes may
   * implement it to give a suitable weight in return. The
   * CascadeHandler may store information about the generated
   * sub-process to be used in the subsequent cascade. It is however
   * not guaranteed that the reweightCKKW() will have been called for
   * the subprocess handed to the handle() function. This default
   * implementation of the function simply return one. The current
   * sub-process is mixed together with other processes with a
   * multiplicity of outgoing particles between \a minMult and \a
   * maxMult.
   */
  virtual double reweightCKKW(int minMult, int maxMult);

public:

  /** @name Access information stored by the handle() function. */
  //@{
  /**
   * Return the vector of tagged particles which should be
   * showered. It the vector is empty, the patons from the current
   * sub-process is supposed to be showered.
   */
  const tPVector & tagged() const { return *theTagged; }

  /**
   * Return the int provided in the current call to handle().
   */
  const Hint & hint() const { return *theHint; }

  /**
   * Return references to the PDF used by the first incoming particle.
   */
  const PDF & firstPDF() const { return pdfs().first; }

  /**
   * Return references to the PDF used by the first incoming particle.
   */
  const PDF & secondPDF() const { return pdfs().second; }

  /**
   * Return references to the currently used PDF's.
   */
  const pair<PDF,PDF> & pdfs() const { return thePDFs; }

  /**
   * Set alternative PDFBase objects to be used for cascade.
   */
  void resetPDFs(const pair<tcPDFPtr,tcPDFPtr> & pdfpair);

  /**
   * Set alternative PDFBase objects to be used for cascade.
   */
  void resetPDFs(const pair<tcPDFPtr,tcPDFPtr> & pdfpair, PBPair ppair);

  /**
   * Set the XComb object with information about the sub-process
   * generation.
   */
  void setXComb(tXCombPtr xc);

  /**
   * Return true, if this cascade handler will perform reshuffling from hard
   * process masses.
   */
  virtual bool isReshuffling() const { return false; }
  
  /**
   * For multiple cascade calls, this flag tells
   * if cascade was called before.
   */
  bool didRunCascade() const {return theDidRunCascade;}
  
  /**
   * Set the flag to inform if prior cascades had been called.
   */
  static void setDidRunCascade(bool c){theDidRunCascade=c;}

  //@}

public:

  /**
   * Standard Init function used to initialize the interface.
   */
  static void Init();

private:

  /**
   * Store the tagged argument given to handle().
   */
  const tPVector * theTagged;

  /**
   * Store the Hint arguments given to handle().
   */
  const Hint * theHint;

  /**
   * The pdfs used to extract the incoming partons.
   */
  pair<PDF,PDF> thePDFs;
  
  /**
   * If there are multiple cascade calls, this flag tells
   * if cascade was called before.
   */
  static bool theDidRunCascade;
  

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class without persistent data.
   */
  static AbstractNoPIOClassDescription<CascadeHandler> initCascadeHandler;

  /**
   *  Private and non-existent assignment operator.
   */
  CascadeHandler & operator=(const CascadeHandler &) = delete;

};

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of CascadeHandler. */
template <>
struct BaseClassTrait<CascadeHandler,1>: public ClassTraitsType {
  /** Typedef of the first base class of CascadeHandler. */
  typedef StepHandler NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the CascadeHandler class and the shared object where it is defined. */
template <>
struct ClassTraits<CascadeHandler>: public ClassTraitsBase<CascadeHandler> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::CascadeHandler"; }
};

/** @endcond */

}

#endif /* ThePEG_CascadeHandler_H */
