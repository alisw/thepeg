// -*- C++ -*-
//
// Tree2toNDiagram.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_Tree2toNDiagram_H
#define ThePEG_Tree2toNDiagram_H
// This is the declaration of the Tree2toNDiagram class.

#include "ThePEG/MatrixElement/DiagramBase.h"
#include "ThePEG/MatrixElement/ColourLines.h"
#include "ThePEG/Handlers/StandardXComb.fh"
#include "Tree2toNDiagram.xh"

namespace ThePEG {

/**
 * The Tree2toNDiagram class inherits from DiagramBase and represents
 * a Feynman tree diagram. It is represented by a chain of \f$n\f$
 * space-like propagators, where one incoming particle has index 1 and
 * other incoming one index \f$n\f$. For adiagram with in total
 * \f$m\f$ propagators the timelike propagators are then numbered
 * \f$n+1\f$ through \f$m\f$. The vector of type of the propagators
 * are accessible from the partons() method, and the parents of
 * propagator \f$i\f$ form the parents(int) method. The parent of a
 * space-like propagator is simply the previous space-like one. The
 * parent of a time-like propagator is either a previous time-like
 * propagator or the first of the connecting space-like ones.
 *
 * A Tree2toNDiagram is created by first constructing it with an
 * integer corresponding to the number of space-like propagators. Then
 * the comma operator is used to add first the particle data objects
 * corresponding to the space-like propagators, then the time-like
 * ones preceeded with the index of their parents. To complete a
 * Tree2toNDiagram, a negative integer is added with the comma
 * operator. This number is then used as an identifier. Note that the
 * parent must have been added before a child is. As an example, the
 * s-channel diagram \f$e \nu_e \rightarrow u \bar{d}\f$ is created
 * thus:<br>
 * <code>Tree2toNDiagram(2),eplus,nue,1,Wplus,3,u,3,dbar</code>.<br>
 * Similarly the t-channel diagram \f$e d \rightarrow \nu_e u\f$ is
 * created thus:<br>
 * <code>Tree2toNDiagram(3),eplus,Wplus,d,1,nu,2,u</code>.  Note that
 * only two chidren are allowed per propagator. This means that
 * four-propagator vertices are not allowed, but must be represented
 * by two three-propagator ones.
 *
 * Please note that for technical reasons, when specifying the
 * diagrams with the comma operator the numbering of the particles is
 * \f$1\ldots m\f$, while the internal representation (in the
 * parent(int) and children(int) function) is using \f$0\ldots m-1\f$
 *
 * @see DiagramBase
 * @see ColourLines
 * 
 */
class Tree2toNDiagram: public DiagramBase {

public:

  /** The integer type reresenting vector sizes. */
  typedef cPDVector::size_type size_type;
  /** A multi-set of particle data objects. */
  typedef multiset<tcPDPtr> PDMSet;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  Tree2toNDiagram()
    : theNSpace(0), theNOutgoing(0), nextOrig(0) {}

  /**
   * Destructor.
   */
  ~Tree2toNDiagram();

  /**
   * The standard constructor giving the number of \a space-like
   * propagators.
   */
  explicit Tree2toNDiagram(int space) 
    : theNSpace(space), theNOutgoing(0), nextOrig(-1) {}
  //@}

public:

  /**
   * If less than zero indicate that this tree is competed. Otherwise
   * signal the parent of the next added parton.
   */
  Tree2toNDiagram & operator,(int o) {
    nextOrig = o - 1;
    if ( o < 0 ) check();
    return *this;
  }

  /**
   * Add a space- or time-like parton.
   */
  Tree2toNDiagram & operator,(PDPtr pd) { return add(pd); }

  /**
   * Add a space- or time-like parton.
   */
  Tree2toNDiagram & operator,(cPDPtr pd) { return add(pd); }

  /**
   * Add a space- or time-like parton.
   */
  Tree2toNDiagram & operator,(tPDPtr pd) { return add(pd); }

  /**
   * Add a space- or time-like parton.
   */
  Tree2toNDiagram & operator,(tcPDPtr pd) { return add(pd); }

  /**
   * Construct a sub process corresponding to this diagram. The
   * incoming partons, and the momenta of the outgoing ones, are given
   * by the XComb object. All parent/children pointers should be set
   * correspondingly and the partons should be colour connected as
   * specified by the ColourLines object.
   */
  virtual tPVector construct(SubProPtr sb, const StandardXComb &,
			     const ColourLines &) const;

  /**
   * Return the types of the incoming partons.
   */
  tcPDPair incoming() const;

  /**
   * Return the complete vector of partons in this tree diagram.
   */
  const cPDVector & allPartons() const { return thePartons; }

  /**
   * Return the outgoing parton types of this tree diagram.
   */
  tcPDVector outgoing() const;

  /**
   * Return the incoming followed by the outgoing parton types of this
   * tree diagram.
   */
  tcPDVector external() const;

  /**
   * Return the index of the parent of the given parton.
   */
  int parent(int i) const { return theParents[i]; }

  /**
   * Return the indices of the children of the given parton.
   */
  pair<int,int> children(int) const;

  /**
   * Return the number of space-like partons
   */
  int nSpace() const { return theNSpace; }

  /**
   * Extend this diagram to accomodate the given number of space-like lines
   */
  void resize(size_type nSpace) {
    theNSpace = max(nSpace,theNSpace);
  }

  /**
   * Return the number of outgoing partons.
   */
  int nOutgoing() const { return theNOutgoing; }

private:

  /**
   * Check the consistency of this tree diagram.
   */
  void check();

  /**
   * Add a space-like parton to this diagram.
   */
  void addSpacelike(tcPDPtr pd) {
    if ( thePartons.size() >= theNSpace ) throw Tree2toNDiagramError();
    theParents.push_back(thePartons.size() - 1);
    thePartons.push_back(pd);
  }
  /**
   * Add a time-like parton to this diagram.
   */
  void addTimelike(tcPDPtr);

  /**
   * Add a time-like parton to this diagram indicating its \a origin.
   */
  void addTimelike(tcPDPtr, size_type origin);

  /**
   * Add a parton to this diagram.
   */
  Tree2toNDiagram & add(tcPDPtr);

public:

  /**
   * Compare this diagram's topology to another one.
   */
  virtual bool isSame(tcDiagPtr) const;

  /**
   * Compare this diagram's topology to another one modulo
   * permutations of external legs; provide a map of this diagram's
   * external legs to the other diagram's external legs.
   */
  virtual bool isSame(tcDiagPtr, map<int,int>&) const;

  /**
   * Check for equality.
   */
  bool equals(Ptr<Tree2toNDiagram>::tcptr, int start=0, int startCmp=0) const;

  /**
   * Check for equality modulo permutations of external legs.
   */
  bool equals(Ptr<Tree2toNDiagram>::tcptr, 
	      map<int,int>&,
	      int start=0, int startCmp=0) const;

  /**
   * Merge the two external partons referred to by indices as in the
   * partons() vector returned by DiagramBase. If both are timelike,
   * the parent will become the new outgoing parton, if one is space-
   * and the other timelike, the spacelike child will become the new
   * incoming parton.  Return the position of the merged parton in the
   * resulting diagram or -1 if the merging is not possible. In
   * addition, return a mapping of a certain (non-merged) external leg
   * id to the id in the merged diagram.
   */
  int mergeEmission(int emitter, int id, map<int,int>& remap);

  /**
   * Translate a parton's id in the diagram to a parton's id in a
   * vector of incoming followed by outgoing partons.
   */
  int externalId(int id) const;

  /**
   * Translate a parton's id in a vector of incoming followed by
   * outgoing partons to a parton's id in the diagram.
   */
  int diagramId(int id) const;

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
   * The number of space-like partons
   */
  size_type theNSpace;

  /**
   * The number of outgoing partons.
   */
  int theNOutgoing;

  /**
   * The parent of the next added parton.
   */
  int nextOrig;

  /**
   * The complete vector of partons in this tree diagram.
   */
  cPDVector thePartons;

  /**
   * The index of the parents.
   */
  vector<int> theParents;

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<Tree2toNDiagram> initTree2toNDiagram;

  /**
   *  Private and non-existent assignment operator.
   */
  Tree2toNDiagram & operator=(const Tree2toNDiagram &) = delete;

};

}

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * This template specialization informs ThePEG about the
 * base class of Tree2toNDiagram.
 */
template <>
struct BaseClassTrait<Tree2toNDiagram,1>: public ClassTraitsType {
  /** Typedef of the base class of Tree2toNDiagram. */
  typedef DiagramBase NthBase;
};

/**
 * This template specialization informs ThePEG about the name of the
 * Tree2toNDiagram class.
 */
template <>
struct ClassTraits<Tree2toNDiagram>: public ClassTraitsBase<Tree2toNDiagram> {
  /** Return the class name. */
  static string className() { return "ThePEG::Tree2toNDiagram"; }
};

/** @endcond */

}

#endif /* ThePEG_Tree2toNDiagram_H */
