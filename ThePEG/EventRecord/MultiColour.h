// -*- C++ -*-
//
// MultiColour.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_MultiColour_H
#define THEPEG_MultiColour_H
//
// This is the declaration of the MultiColour class.
//

#include "ThePEG/EventRecord/ColourBase.h"

namespace ThePEG {

/**
 * This class is used to store colour information of RemnantParticle
 * objects and other particle classes with complicated colour
 * structures. Rather than just having a 
 */
class MultiColour: public ColourBase {

public:

  using ColourBase::colourLine;
  using ColourBase::antiColourLine;

  /**
   * Return the anti-colour lines to which this particle is
   * connected.
   */
  virtual vector<tcColinePtr> antiColourLines() const;

  /**
   * Return the colour lines to which this particle is
   * connected.
   */
  virtual vector<tcColinePtr> colourLines() const;

  /**
   * Add the given (\a anti-) colour \a line to the particle. If the base
   * class has no (anti-) colour line, it will also be set.
   */
  virtual void colourLine(tColinePtr line, bool anti = false);

  /**
   * Add the given (\a anti-) colour \a line to the particle. If the base
   * class has no (anti-) colour line, it will also be set.
   */
  virtual void colourLine(tColinePtr line, int index, bool anti = false);

  /**
   * Add the given anti-colour \a line to the particle. If the base
   * class has no anti-colour line, it will also be set.
   */
  virtual void antiColourLine(tColinePtr line);

  /**
   * Add the given anti-colour \a line to the particle. If the base
   * class has no anti-colour line, it will also be set.
   */
  virtual void antiColourLine(tColinePtr line, int index);

  /**
   * Remove the given (\a anti-) colour \a line from the particle. If
   * the line is the colourLine() of the base class, it will be
   * removed there as well.
   */
  virtual void removeColourLine(tcColinePtr line, bool anti = false);

  /**
   * Remove the given anti-colour \a line from the particle. If the
   * line is the antiColourLine() of the base class, it will be
   * removed there as well.
   */
  virtual void removeAntiColourLine(tcColinePtr line);

  /**
   * Return true if the particle is connected to the given (\a anti-)
   * colour \a line.
   */
  virtual bool hasColourLine(tcColinePtr line, bool anti = false) const;

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

  /**
   * Standard clone method.
   */
  virtual EIPtr clone() const { return new_ptr(*this); }

private:

  /**
   * The set of colour lines to which a particle is attached.
   */
  list<cColinePtr> theColourLines;

  /**
   * The set of anti-colour lines to which a particle is attached.
   */
  list<cColinePtr> theAntiColourLines;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is a concrete class with persistent data.
   */
  static ClassDescription<MultiColour> initMultiColour;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  MultiColour & operator=(const MultiColour &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of MultiColour. */
template <>
struct BaseClassTrait<MultiColour,1> {
  /** Typedef of the first base class of MultiColour. */
  typedef ColourBase NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the MultiColour class and the shared object where it is defined. */
template <>
struct ClassTraits<MultiColour>
  : public ClassTraitsBase<MultiColour> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::MultiColour"; }
  /**
   * The name of a file containing the dynamic library where the class
   * MultiColour is implemented. It may also include several, space-separated,
   * libraries if the class MultiColour depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return "MultiColour.so"; }
};

/** @endcond */

}

#endif /* THEPEG_MultiColour_H */
