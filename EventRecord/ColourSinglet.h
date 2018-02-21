// -*- C++ -*-
//
// ColourSinglet.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_ColourSinglet_H
#define ThePEG_ColourSinglet_H

#include "ThePEG/EventRecord/EventConfig.h"

namespace ThePEG {

/**
 * A <code>ColourSinglet</code> object represents a set of
 * colour-connected partons in a total colour-singlet state. The
 * simplest form is a single string represented by a vector of partons
 * starting from a colour-triplet parton followed by a number of
 * colour-octet partons and ending in a anti-triplet. A closed string
 * is represented by a vector of colour-octet partons. If junctions
 * are present, there will be a vector of string pieces containing the
 * associated partons, with associated information about the string
 * pieces connected via sinks or sources.
 *
 * @see ColourLine
 * @see Particle
 * 
 */
class ColourSinglet {

public:

  /** A piece of string. */
  typedef deque<tcPPtr> StringPiece;
  /** An integer type used as index in a string pieece. */
  typedef StringPiece::size_type Index;
  /** Representaion of a junction. */
  typedef pair<Index,Index> Junction;

public:

  /**
   * Default constructor.
   */
  ColourSinglet() {}

  /**
   * Constructor taking an initial colour line and a set of partons to
   * select from.
   */
  ColourSinglet(tcColinePtr cl, tcParticleSet & left);

protected:

  /**
   * Internal constructor. When splitting an internal line, create a
   * new colour singlet starting from the string piece \a si in the
   * singlet \a cs.
   */
  ColourSinglet(const ColourSinglet & cs, Index si);

public:

  /**
   * Access the vector of partons.
   */
  tcPVector & partons() { return thePartons; }

  /**
   * Access the vector of partons.
   */
  const tcPVector & partons() const { return thePartons; }

  /**
   * Access a parton giving an index.
   */
  tcPPtr parton(tcPVector::size_type i) const { return thePartons[i]; }

  /**
   * Return the total momentum for the partons in the colour singlet.
   */
  LorentzMomentum momentum() const;

  /**
   * Return the number of string pieces.
   */
  Index nPieces() const { return thePieces.size(); }

  /**
   * Return the partons belonging to string piece with index i (note
   * that the indices starts with 1).
   */
  StringPiece & piece(Index i) { return thePieces[i - 1]; }

  /**
   * Return the partons belonging to string piece with index i (note
   * that the indices starts with 1).
   */
  const StringPiece & piece(Index i) const { return thePieces[i - 1]; }

  /**
   * Return the sink (or source if \a forward is false) neighbors of
   * a string piece with index i (note that the indices starts with 1).
   */
  Junction & junction(Index i, bool forward) {
    return forward ? sink(i) : source(i);
  }

  /**
   * Return the sink (or source if \a forward is false) neighbors of
   * a string piece with index i (note that the indices starts with 1).
   */
  const Junction & junction(Index i, bool forward) const {
    return forward ? sink(i) : source(i);
  }

  /**
   * Return the sink or source neighbors of a string piece with index
   * i (note that the indices starts with 1).
   */
  Junction & sink(Index i) { return theSinks[i - 1]; }

  /**
   * Return the sink or source neighbors of a string piece with index
   * i (note that the indices starts with 1).
   */
  const Junction & sink(Index i) const { return theSinks[i - 1]; }

  /**
   * Return the sink or source neighbors of a string piece with index
   * i (note that the indices starts with 1).
   */
  Junction & source(Index i) { return theSources[i - 1]; }
  /**
   * Return the sink or source neighbors of a string piece with index
   * i (note that the indices starts with 1).
   */
  const Junction & source(Index i) const { return theSources[i - 1]; }

  /**
   * Extract colour-singlet strings/clusters of partons from the given
   * set. The set will be empty afterwards if all went well - even
   * colour-singlet particles will be removed.
   */
  static vector<ColourSinglet> getSinglets(tcParticleSet & left);

  /**
   * Extract colour-singlet strings/clusters of partons from the given
   * range of particles.
   */
  template <typename Iterator>
  static vector<ColourSinglet>
  getSinglets(Iterator first, Iterator last) {
    tcParticleSet pset(first, last);
    return getSinglets(pset);
  }

  /**
   * Return a vector with data objects corresponding to all triplet
   * partons in this singlet.
   */
  tcPDVector getTripletData() const;

  /**
   * Split a junction string. Split an internal string piece \a sp
   * (one that starts from a source and ends in a sink). Randomly pair
   * together the connectin string pieces and randomly assign possible
   * colour octet partons to the two joined string pieces. The singlet
   * will then split into two - one is returned and this will contain
   * the other. If no argument is given, a random internal line will
   * be picked.
   */
  ColourSinglet splitInternal(Index sp = 0);

  /**
   * Split a junction string. Split an internal string piece \a sp
   * (one that starts from a source and ends in a sink). The source
   * neighbour \a sa will be joined with the sink neighbour \a
   * sc. Possible colour octet partons in the split string piece will
   * be assigned to the joined \a sa-\a sc string piece if the
   * corresponing \a assign element is true. The singlet will then
   * split into two - one is returned and this will contain the other.
   */
  ColourSinglet splitInternal(Index sp, Index sa, Index sc,
			      const vector<bool> & assing);

  /**
   * Split a junction string. If a sting piece \a sp is connected to a
   * junction and ends in a diquark diq, split the diquark into two
   * quarks, remove the junction, adding the other partons to the
   * first neighbouring string piece if the corresponding \a assign
   * element is true. The singlet will then split into two - one is
   * returned and this will contain the other. Note that the
   * parent-children relationship and colour flows of the diquark and
   * the two quarks is not handled by this function.
   */
  ColourSinglet
  splitDiQuarkJunction(Index sp, tcPPtr diq, tcPPair qq,
		       const vector<bool> & assign = vector<bool>());

  /**
   * Split a di-quark string. If this singlet is a simple string with
   * (anti-)diquarks in both ends, split the diquarks into the given
   * quark pairs and create two new string pieces the betwen the quark
   * and anit-quarks, adding the remaining partons to the first piece
   * if the corresponding \a assign element is true. The singlet will
   * then split into two - one is returned and this will contain the
   * other. Note that the parent-children relationship and colour
   * flows of the diquark and the two quarks is not handled by this
   * function.
   */
  ColourSinglet splitDiDiQuark(tcPPair qq1, tcPPair qq2,
			       const vector<bool> & assign = vector<bool>());

  /**
   * Swap this colour singlet for the argument.
   */
  void swap(ColourSinglet & x) {
    thePartons.swap(x.thePartons);
    thePieces.swap(x.thePieces);
    theSources.swap(x.theSources);
    theSinks.swap(x.theSinks);
  }

private:

  /**
   * Fill a string piece. Follow a colour line \a forward in colour
   * removing partons from the \a left set and adding them to the
   * string piece \a s0 as we go along. If a line ends in a colour
   * sink/source, follow the other two colour lines in turn with the
   * value of \a forward reversed.
   */
  bool fill(Index s0, bool forward, tcColinePtr first, tcParticleSet & left);

  /**
   * Fill a string piece. When creating a new singlet from an old one
   * which has been split, add the string piece \a i1. If it ends in a
   * sink or a source add the neighbors recursively.
   */
  void fill(Index i0, bool forward, const ColourSinglet & cs, Index i1);

  /**
   * Add a new string piece and return its index.
   */
  Index addPiece() {
    thePieces.push_back(StringPiece());
    theSinks.push_back(Junction());
    theSources.push_back(Junction());
    return nPieces();
  }

  /**
   * Add a junction which is a source (or sink if \a forward) to the
   * string piece \a s0. Return the indices of the two connected
   * string pieces.
   */
  Junction addJunction(Index s0, bool forward);

public:

  /** @cond EXCEPTIONCLASSES */
  /** An exception class to indicate that something went wrong with a
   *  ColourSinglet. */
  class ColourSingletException: public Exception {};
  /** @endcond */

private:

  /**
   * The vector of all partons.
   */
  tcPVector thePartons;

  /**
   * The different string pieces.
   */
  vector<StringPiece> thePieces;

  /**
   * The source neighbours of the string pieces.
   */
  vector<Junction> theSources;

  /**
   * The sink neighbours of the string pieces.
   */
  vector<Junction> theSinks;

};

}

#endif /* ThePEG_ColourSinglet_H */
