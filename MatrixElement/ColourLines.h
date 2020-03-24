// -*- C++ -*-
//
// ColourLines.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_ColourLines_H
#define ThePEG_ColourLines_H
// This is the declaration of the ColourLines class.

#include "ThePEG/Config/ThePEG.h"

namespace ThePEG {

/**
 * The ColourLines class defines the colour flow in a SubProcess. It
 * defines a number of colour lines and specifies which particles are
 * connected to them.
 * 
 */
class ColourLines: public Base {

public:

  /** A single colour line */
  typedef vector<pair<int,int> > Line;
  /** A vector of colour lines. */
  typedef vector<Line> LineVector;
  /** A vector of <code>ColourLine</code>. */
  typedef vector<ColinePtr> Vertex;
  /** A vector of vertices. */
  typedef vector<Vertex> VertexVector;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  ColourLines() {}

  /**
   * The standard constructor. The string \a s should contain a
   * comma-separated sequence of integers. Each sequence of numbers
   * indicates a colour line and the integer represents a parton
   * connected to it. If the integer is negative, it means that the
   * line is the corresponding partons anti-colour. Note that the
   * partons are numbered from 1: The first spacelike particle is 1, the second
   * is 2 and the internal time-like and outgoing are numbered after all
   * the spacelike particles.
   */
  ColourLines(string s);
  //@}

  /** 
   * Reset this ColourLines object. The string \a s should contain a 
   * comma-separated sequence of integers. Each sequence of numbers 
   * indicates a colour line and the integer represents a parton 
   * connected to it. If the integer is negative, it means that the 
   * line is the corresponding partons anti-colour. Note that the 
   * partons are numbered from 1: The first incoming is 1, the second 
   * is 2 and the internal and outgoing are numbered 3 and upwards. 
   */ 
  void reset(string s); 

public:

  /**
   * Create the corresponding <code>ColourLine</code>s and connect the
   * given \a partons. The partons are assumed to be in the same order
   * as the numbers specified in the constructor.
   */
 void connect(const tPVector & partons) const;

private:

  /**
   * The vector of colour lines.
   */
  LineVector theLines;

};

}

#endif /* ThePEG_ColourLines_H */
