// -*- C++ -*-
//
// ColourLine.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_ColourLine_H
#define ThePEG_ColourLine_H
// This is the declaration of the ColourLine class.

#include "EventConfig.h"
#include "ThePEG/Utilities/ClassDescription.h"
#include "ThePEG/EventRecord/ColourSinglet.h"

namespace ThePEG {

/**
 * The ColourLine class represents colour lines connecting
 * <code>Particle</code>s. A <code>ColourLine</code> keeps track on
 * the particles connected to it. To connect a particle to a colour
 * line the <code>addColoured()</code> and
 * <code>addAntiColoured()</code> functions should be used - these
 * will automatically set up the Particle correctly. There is no
 * method in a Particle to directly set its colour lines.
 *
 * If a colour line stems from a colour source or ends in a colour
 * sink, it is possible to obtain the neighbouring colour lines. This
 * is also the way junction strings and sinks and sources are
 * implemented.
 *
 * @see Particle
 * @see ColourBase
 */
class ColourLine: public EventRecordBase {

public:

  /** @name Creation functions. */
  //@{
  /**
   * Create a colour line. Set a pair of colour - anticolour particles
   * in a newly created colour line.
   */
  static tColinePtr create(tPPtr col, tPPtr anti);

  /**
   * Create a colour line. Set a particle for which the created object
   * is a (anti-)colour line .
   * @param p the particle to be connected.
   * @param anti if true, the created object is the anti-colour line
   * of \a p.
   */
  static tColinePtr create(tPPtr p, bool anti = false);

  /**
   * Create a colour line. Set a particle for which the created object
   * is a anti-colour line .
   * @param p the particle to be connected.
   */
  static tColinePtr createAnti(tPPtr p) { return create(p, true); }

  /**
   * Create a coloue line which is a connector between two junctions,
   * a source junction with neigboring colour lines \a son1 and \a
   * son2 and a sink junction with neigboring colour lines \a sin1 and
   * \a sin2.
   */
  static tColinePtr create(tColinePtr son1, tColinePtr son2,
			   tColinePtr sin1, tColinePtr sin2);
  //@}

  /**
   * Destructor.
   */
  virtual ~ColourLine();

public:

  /** @name Access particles connected to the colour line. */
  //@{
  /**
   * Return the vectors of particles connected to this line with their
   * colours.
   */
  const tPVector & coloured() const { return theColoured; }

  /**
   * Return the vectors of particles connected to this line with their
   * anti-colours.
   */
  const tPVector & antiColoured() const { return theAntiColoured; }

  /**
   * Return the first particle on this colour line. Returns null if
   * this line stems from a colour source. If the particle is
   * outgoing, its anti colour is connected, otherwise its colour is
   * connected.
   */
  tPPtr startParticle() const;

  /**
   * Return the last particle on this colour line. Returns null if
   * this line ends in a colour sink. If the particle is outgoing, its
   * colour is connected, otherwise its anti colour is connected.
   */
  tPPtr endParticle() const;

  //@}

  /** @name Add and remove particles in a colour line. */
  //@{
  /**
   * Add a particle having this as a anti-colour line.
   */
  void addAntiColoured(tPPtr);

  /**
   * Add a particle having this as a (anti-)colour line.
   * @param p the particle to be connected.
   * @param anti if true, this is the anti-colour line of \a p.
   */
  void addColoured(tPPtr p, bool anti = false);

  /**
   * Add a particle having this as a anti-colour line at a given index.
   */
  void addAntiColouredIndexed(tPPtr p, int index);

  /**
   * Add a particle having this as a (anti-)colour line at a given index.
   * @param p the particle to be connected.
   * @param anti if true, this is the anti-colour line of \a p.
   */
  void addColouredIndexed(tPPtr p, int index, bool anti=false);

  /**
   * Remove a particle having this as an anti-colour line.
   */
  void removeAntiColoured(tPPtr);

  /**
   * Remove a particle having this as a (anti-)colour line.
   * @param p the particle to be removed.
   * @param anti if true, this is the anti-colour line of \a p.
   */
  void removeColoured(tPPtr p, bool anti = false);

  //@}

  /** @name Functions for junction strings. */
  //@{
  /**
   * If this colour line ends in a colour sink, these two colour lines
   * ends in the same.
   */
  tColinePair sinkNeighbours() const { return theSinkNeighbours; }

  /**
   * If this colour line stems from a colour source (sink), these two colour
   * lines stems from (ends in) the same.
   * @param anti if true return sinkNeighbours().
   */
  tColinePair sourceNeighbours(bool anti = false) const {
    return anti? theSinkNeighbours: theSourceNeighbours;
  }

  /**
   * Add two colour lines as neighbours to this line. Afterwards all
   * three will end in the same sink. Also the neighbors are set up
   * correspondingly.
   */
  void setSinkNeighbours(tColinePtr l1, tColinePtr l2) {
    theSinkNeighbours.second = l1->theSinkNeighbours.second = l2;
    l2->theSinkNeighbours.second = theSinkNeighbours.first = l1;
    l1->theSinkNeighbours.first = l2->theSinkNeighbours.first = this;
  }

  /**
   * Add two colour lines as neighbours to this line. Afterwards all
   * three will stem from the same source. Also the neighbors are set
   * up correspondingly.
   */
  void setSourceNeighbours(tColinePtr l1, tColinePtr l2) {
    theSourceNeighbours.second = l1->theSourceNeighbours.second = l2;
    l2->theSourceNeighbours.second = theSourceNeighbours.first = l1;
    l1->theSourceNeighbours.first = l2->theSourceNeighbours.first = this;
  }

  //@}

  /**
   * Join with the given ColourLine. The colour of the given \a line
   * is joined so that it will flow into this line, ie. the
   * anti-coloured particle in the end of the \a line will become
   * connected to the coloured particle in the of this line.  After
   * the joining the given \a line will not be connected to
   * anything.
   */
  bool join(ColinePtr line);

  /**
   * Return the first (anti-)coloured parton among the given range of
   * particles which is on this colour line.
   */
  template <typename Iterator>
  typename std::iterator_traits<Iterator>::value_type
  getColouredParticle(Iterator first, Iterator last, bool anti = false) const {
    typedef typename std::iterator_traits<Iterator>::value_type ParticlePointer;
    for ( ; first != last; ++first )
      if ( (**first).coloured() && (**first).hasColourLine(this, anti) )
	return *first;
    return ParticlePointer();
  }

  /**
   * Write out information about this colour line to the stream.
   */
  void write(ostream & os, tcEventPtr event, bool anti) const;

public:

  /**
   * Standard function for writing to a persistent stream.
   */
  void persistentOutput(PersistentOStream &) const;

  /**
   * Standard function for reading from a persistent stream.
   */
  void persistentInput(PersistentIStream &, int);

private:

  /**
   * The particles connecting to this colour line, following the
   * incoming colour flow.
   */
  tPVector theColoured;

  /**
   * The particles connecting to this colour line, following the
   * outgoing colour flow.
   */
  tPVector theAntiColoured;

  /**
   * If this colour line stems from a colour source, these two colour
   * lines stems from the same.
   */
  tColinePair theSourceNeighbours;

  /**
   * If this colour line ends in a colour sink, these two colour lines
   * ends in the same.
   */
  tColinePair theSinkNeighbours;

  /**
   * Colour lines which are connectors between two junctions do not
   * have a particle which owns it, instead it is owned by one of the
   * source neighbours.
   */
  vector<ColinePtr> orphanedConnectors;

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<ColourLine> initColourLine;

  /**
   *  Private and non-existent assignment operator.
   */
  ColourLine & operator=(const ColourLine &) = delete;

};

}


namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/**
 * The following template specialization informs ThePEG about the
 * base class of ColourLine.
 */
template <>
struct BaseClassTrait<ColourLine,1>: public ClassTraitsType {
  /** Typedef of the first base class of ColourLine. */
  typedef EventRecordBase NthBase;
};

/**
 * The following template specialization informs ThePEG about the
 * name of this class and the shared object where it is defined.
 */
template <>
struct ClassTraits<ColourLine>: public ClassTraitsBase<ColourLine> {
  /** Return the class name. */
  static string className() { return "ThePEG::ColourLine"; }
  /** Return the name of the shared library to be loaded to get
   * access to this class. */
  static string library() { return "ColourLine.so"; }
};

/** @endcond */

}

#endif /* ThePEG_ColourLine_H */
