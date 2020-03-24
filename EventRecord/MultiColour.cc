// -*- C++ -*-
//
// MultiColour.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MultiColour class.
//

#include "MultiColour.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include <algorithm>
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;

vector<tcColinePtr> MultiColour::antiColourLines() const {
  return vector<tcColinePtr>(theAntiColourLines.begin(),
			     theAntiColourLines.end());
}

vector<tcColinePtr> MultiColour::colourLines() const {
  return vector<tcColinePtr>(theColourLines.begin(),
			     theColourLines.end());
}

bool MultiColour::hasColourLine(tcColinePtr line, bool anti) const {
  return ( anti? ( find(theAntiColourLines.begin(),theAntiColourLines.end(),line) != 
		   theAntiColourLines.end() ):
	   ( find(theColourLines.begin(),theColourLines.end(),line) !=
	     theColourLines.end() ) );
}

void MultiColour::colourLine(tColinePtr line, bool anti) {
  if ( anti ) antiColourLine(line);
  else {
    if ( !colourLine() ) ColourBase::colourLine(line);
    if(find(theColourLines.begin(),theColourLines.end(),line)==
       theColourLines.end())
      theColourLines.push_back(line);
  }
}

void MultiColour::colourLine(tColinePtr line, int index, bool anti) {
  if ( anti ) {
    antiColourLine(line,index);
    return;
  }
  if ( !colourLine() ) ColourBase::colourLine(line);
  if(find(theColourLines.begin(),theColourLines.end(),line)!=theColourLines.end())
    return;
  int ix=0;
  for(list<cColinePtr>::iterator it=theColourLines.begin();
      it!=theColourLines.end();++it) {
    ++ix;
    if(ix==index) {
      it = theColourLines.insert(it,line);
      ++it;
      removeColourLine(*it);
      if ( !colourLine() ) ColourBase::colourLine(line);
      return;
    }
  }
  for(;ix<index-1;++ix) {
    theColourLines.push_back(ColinePtr());
  }
  theColourLines.push_back(line);
  if ( !colourLine() ) ColourBase::colourLine(line);
}

void MultiColour::antiColourLine(tColinePtr line, int index) {
  if ( !antiColourLine() ) ColourBase::antiColourLine(line);
  if(find(theAntiColourLines.begin(),theAntiColourLines.end(),line)!=theAntiColourLines.end())
    return;
  int ix=0;
  for(list<cColinePtr>::iterator it=theAntiColourLines.begin();
      it!=theAntiColourLines.end();++it) {
    ++ix;
    if(ix==index) {
      it = theAntiColourLines.insert(it,line);
      ++it;
      removeAntiColourLine(*it);
      if ( !antiColourLine() ) ColourBase::antiColourLine(line);
      return;
    }
  } 
  for(;ix<index-1;++ix) {
    theAntiColourLines.push_back(ColinePtr());
  }
  theAntiColourLines.push_back(line);
  if ( !antiColourLine() ) ColourBase::antiColourLine(line);
}

void MultiColour::antiColourLine(tColinePtr line) {
  if ( !antiColourLine() ) ColourBase::antiColourLine(line);
  if(find(theAntiColourLines.begin(),theAntiColourLines.end(),line)==
     theAntiColourLines.end())
    theAntiColourLines.push_back(line);
}

void MultiColour::removeColourLine(tcColinePtr line, bool anti) {
  if ( anti ) removeAntiColourLine(line);
  else {
    ColourBase::removeColourLine(line);
    theColourLines.remove(line);
  }
}

void MultiColour::removeAntiColourLine(tcColinePtr line) {
  ColourBase::removeAntiColourLine(line);
  theAntiColourLines.remove(line);
}

void MultiColour::persistentOutput(PersistentOStream & os) const {
  os << theColourLines << theAntiColourLines;
}

void MultiColour::persistentInput(PersistentIStream & is, int) {
  is >> theColourLines >> theAntiColourLines;
}

ClassDescription<MultiColour> MultiColour::initMultiColour;
// Definition of the static class description member.

void MultiColour::Init() {

  static ClassDocumentation<MultiColour> documentation
    ("There is no documentation for the MultiColour class");

}

