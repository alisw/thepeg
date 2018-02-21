// -*- C++ -*-
//
// Particle.tcc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the templated member functions of
// the Particle class.
//


namespace ThePEG {

struct ParticleNumLess {
  bool operator()(tPPtr p1, tPPtr p2)
    { return p1->number() < p2->number(); }
};

template <typename Iterator>
void Particle::
PrintParticles(ostream & os, Iterator first, Iterator last, tcStepPtr step) {
  typedef multiset<PPtr,ParticleNumLess> NumSet;
  NumSet pset(first, last);
  for ( NumSet::iterator i = pset.begin(); i != pset.end(); ++i )
    (**i).print(os, step);
}

template <typename Iterator>
typename std::iterator_traits<Iterator>::value_type Particle::
colourNeighbour(Iterator first, Iterator last, bool anti) const {
  if ( !hasColourInfo() || !hasColour(!anti) ) return tPPtr();
  tColinePtr line = colourLine(!anti);
  for ( ; first != last; ++first )
    if ( (**first).hasColourInfo() && (**first).hasColourLine(line, anti) )
      return *first;
  return tPPtr();
}

}

