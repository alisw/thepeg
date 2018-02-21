// -*- C++ -*-
//
// BlobDiagram.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_BlobDiagram_H
#define ThePEG_BlobDiagram_H
// This is the declaration of the BlobDiagram class.

#include "ThePEG/MatrixElement/DiagramBase.h"
#include "ThePEG/MatrixElement/ColourLines.h"
#include "ThePEG/Handlers/StandardXComb.fh"

namespace ThePEG {

/**
 * The BlobDiagram class inherits from DiagramBase and represents a general
 * Feynman diagram of which no further substructure is assumed.
 *
 * @see DiagramBase
 * @see ColourLines
 * 
 */
class BlobDiagram: public DiagramBase {

public:

  /** The integer type reresenting vector sizes. */
  typedef cPDVector::size_type size_type;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor
   */
  BlobDiagram() 
    : DiagramBase() {}

  /**
   * Constructor specifiying incoming partons
   */
  BlobDiagram(int id, tcPDPtr first, tcPDPtr second) 
    : DiagramBase() {
    addParton(first);
    addParton(second);
    diagramInfo(2,id);
  }

  /**
   * Destructor.
   */
  ~BlobDiagram();
  //@}

public:

  /**
   * Add a space- or time-like parton.
   */
  BlobDiagram& operator,(PDPtr pd) { addParton(pd); return *this; }

  /**
   * Add a space- or time-like parton.
   */
  BlobDiagram& operator,(cPDPtr pd) { addParton(pd); return *this; }

  /**
   * Add a space- or time-like parton.
   */
  BlobDiagram& operator,(tPDPtr pd) { addParton(pd); return *this; }

  /**
   * Add a space- or time-like parton.
   */
  BlobDiagram& operator,(tcPDPtr pd) { addParton(pd); return *this; }

  /**
   * Construct a sub process corresponding to this diagram. The
   * incoming partons, and the momenta of the outgoing ones, are given
   * by the XComb object. All parent/children pointers should be set
   * correspondingly and the partons should be colour connected as
   * specified by the ColourLines object.
   */
  virtual tPVector construct(SubProPtr sb, const StandardXComb&,
			     const ColourLines&) const;

  /**
   * Return the types of the incoming partons.
   */
  tcPDPair incoming() const {
    return tcPDPair(partons()[0],partons()[1]); 
  }

  /**
   * Return the outgoing parton types of this tree diagram.
   */
  tcPDVector outgoing() const {
    return tcPDVector(partons().begin()+2,partons().end());
  }

  /**
   * Return the incoming followed by the outgoing parton types of this
   * tree diagram.
   */
  tcPDVector external() const {
    return tcPDVector(partons().begin(),partons().end());
  }

  /**
   * Return the number of outgoing partons.
   */
  size_type nOutgoing() const { return partons().size() - 2; }

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

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<BlobDiagram> initBlobDiagram;

  /**
   *  Private and non-existent assignment operator.
   */
  BlobDiagram & operator=(const BlobDiagram &);

};

}

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * This template specialization informs ThePEG about the
 * base class of BlobDiagram.
 */
template <>
struct BaseClassTrait<BlobDiagram,1>: public ClassTraitsType {
  /** Typedef of the base class of BlobDiagram. */
  typedef DiagramBase NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * BlobDiagram class.
 */
template <>
struct ClassTraits<BlobDiagram>: public ClassTraitsBase<BlobDiagram> {
  /** Return the class name. */
  static string className() { return "ThePEG::BlobDiagram"; }
};

/** @endcond */

}

#endif /* ThePEG_BlobDiagram_H */
