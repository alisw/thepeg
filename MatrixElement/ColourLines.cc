// -*- C++ -*-
//
// ColourLines.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ColourLines class.
//

#include "ColourLines.h"
#include "ColourLines.xh"
#include "ThePEG/EventRecord/ColourLine.h"
#include "ThePEG/EventRecord/MultiColour.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Utilities/StringUtils.h"

using namespace ThePEG;

ColourLines::ColourLines(string s) {
  reset(s); 
} 

void ColourLines::reset(string s) { 
  theLines.clear(); 
  while ( true ) {
    string line = StringUtils::car(s, ",");
    line = StringUtils::stripws(line);
    Line l;
    while (line!="") {
      string loc = StringUtils::car(line);
      string first  = StringUtils::car(loc,":");
      string second = StringUtils::cdr(loc,":");
      if(second!="") {
	int i;
	istringstream is(first);
	is >> i;
	int j;
	istringstream is2(second);
	is2 >> j;
	l.push_back(make_pair(i,j));
      }
      else {
	int i;
	istringstream is(first);
	is >> i;
	l.push_back(make_pair(i,0));
      }
      line = StringUtils::cdr(line);
    };
    if ( l.empty() ) return;
    theLines.push_back(l);
    s = StringUtils::cdr(s, ",");
  }
}

void ColourLines::connect(const tPVector & partons) const {
  VertexVector sinks;
  VertexVector sources;
  long np = partons.size();

  // Create each line and connect the specified partons to them. Save
  // all lines coming from a source or ending in a sink.
  for ( LineVector::size_type il = 0; il < theLines.size(); ++il ) {
    const Line & line = theLines[il];
    ColinePtr cline = new_ptr(ColourLine());
    for ( Line::size_type i = 0; i < line.size(); ++i ) {
      if ( line[i].first > np ) {
	// this is a colour source.
	int is = line[i].first - np;
	sources.resize(is);
	sources[is - 1].push_back(cline);
      } else if ( -line[i].first > np ) {
	// this is a colour sink.
	int is = -line[i].first - np;
	sources.resize(is);
	sources[is - 1].push_back(cline);
      } else if ( line[i].first > 0 ) {
	// This is a coloured particle.
	if ( !partons[line[i].first - 1]->hasColour() )
	  throw ColourGeometryException(partons, line);
	if(line[i].second==0) {
	  cline->addColoured(partons[line[i].first - 1]);
	}
	else {
	  Ptr<MultiColour>::pointer colour = 
	    dynamic_ptr_cast<Ptr<MultiColour>::pointer>(partons[line[i].first - 1]->colourInfo());
	  assert(colour);
	  colour->colourLine(cline,line[i].second);
	}
      } else {
	if ( !partons[-line[i].first - 1]->hasAntiColour() )
	  throw ColourGeometryException(partons, line);
	if(line[i].second==0) {
	  cline->addAntiColoured(partons[-line[i].first - 1]);
	}
	else {
	  Ptr<MultiColour>::pointer colour = 
	    dynamic_ptr_cast<Ptr<MultiColour>::pointer>(partons[-line[i].first - 1]->colourInfo());
	  assert(colour);
	  colour->antiColourLine(cline,line[i].second);
	}
      }
    }
  }

  // Now connect up all lines steming from sources.
  for ( VertexVector::size_type i = 0; i < sources.size(); ++i ) {
    if ( sources[i].empty() ) continue;
    if ( sources[i].size() != 3 ) throw ColourGeometryException(partons,
                                                                vector<pair<int,int> >() );
    sources[i][0]->setSourceNeighbours(sources[i][1], sources[i][2]);
  }

  // Now connect up all lines ending in sinks.
  for ( VertexVector::size_type i = 0; i < sinks.size(); ++i ) {
    if ( sinks[i].empty() ) continue;
    if ( sinks[i].size() != 3 ) throw ColourGeometryException(partons,
                                                              vector<pair<int,int> >());
    sinks[i][0]->setSinkNeighbours(sinks[i][1], sinks[i][2]);
  }
}

ColourGeometryException::
ColourGeometryException(const tPVector & p, const vector<pair<int,int> > & c) {
  if ( c.empty() )
    theMessage << "The number of colour lines steming from one colour source "
	       << "or ending in one colour sink was not equal to three.\n";
  else {
    theMessage << "Cannot connect the following partons:\n";
    for ( unsigned i = 0; i < p.size(); ++i )
      theMessage << " " << p[i]->PDGName();
    theMessage << "\n to the following colour line:\n";
    for ( unsigned i = 0; i < c.size(); ++i ) theMessage << " (" << c[i].first << ","
							 << c[i].second << ") ";
    theMessage << endl;
  }
  severity(maybeabort);
}

