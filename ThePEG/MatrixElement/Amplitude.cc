// -*- C++ -*-
//
// Amplitude.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Amplitude class.
//

#include "Amplitude.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;

AbstractNoPIOClassDescription<Amplitude> Amplitude::initAmplitude;
// Definition of the static class description member.

void Amplitude::Init() {
  static ClassDocumentation<Amplitude> documentation
    ("This the abstract class from which any amplitude class, associated with ",
     "a vertex inherits from.");
}



Complex Amplitude::overestimateValue( const tcPDVector & particles,
				      const vector<Lorentz5Momentum> & momenta, 
				      const vector<int> & helicities ) {
  return value(particles,momenta,helicities);
}


Complex Amplitude::value( const PVector & particles,
			  const vector<int> & helicities ) {
  tcPDVector dataParticles;
  vector<Lorentz5Momentum> momenta;
  for ( PVector::const_iterator cit = particles.begin();
	cit != particles.end(); ++cit ) {
    dataParticles.push_back( (*cit)->dataPtr() );
    momenta.push_back( (*cit)->momentum() );
  }
  return value(dataParticles,momenta,helicities);
}


Complex Amplitude::overestimateValue( const PVector & particles,
				      const vector<int> & helicities ) {
  tcPDVector dataParticles;
  vector<Lorentz5Momentum> momenta;
  for ( PVector::const_iterator cit = particles.begin();
	cit != particles.end(); ++cit ) {
    dataParticles.push_back( (*cit)->dataPtr() );
    momenta.push_back( (*cit)->momentum() );
  }
  return overestimateValue(dataParticles,momenta,helicities);
}


