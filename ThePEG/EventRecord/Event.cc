// -*- C++ -*-
//
// Event.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined member functions of
// the Event class.
//
#include "Event.h"
#include "ThePEG/EventRecord/SubProcess.h"
#include "ThePEG/Utilities/Rebinder.h"
#include "ThePEG/EventRecord/ParticleTraits.h"
#include "ThePEG/Config/algorithm.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/DIterator.h"
#include <iostream>

using namespace ThePEG;

Event::Event(const PPair & newIncoming, tcEventBasePtr newHandler,
	     string newName, long newNumber, double newWeight)
  : Named(newName), theIncoming(newIncoming), theHandler(newHandler),
    theNumber(newNumber), theWeight(newWeight), theParticleNumber(0) {
  addParticle(incoming().first);
  addParticle(incoming().second);
}

Event::Event(const Event & e)
  : Base(e), Named(e), 
    theIncoming(e.theIncoming), theCollisions(e.theCollisions),
    allSteps(e.allSteps), allSubProcesses(e.allSubProcesses),
    allParticles(e.allParticles), theHandler(e.theHandler),
    theNumber(e.theNumber), theWeight(e.theWeight),
    theParticleNumber(e.theParticleNumber) {}

Event::~Event() {
  for ( int i = 0, N = theCollisions.size(); i < N; ++i )
    if ( theCollisions[i]->event() == this )
      theCollisions[i]->theEvent = tEventPtr();
  theIncoming = PPair();
  theCollisions.clear();
  allSteps.clear();
  allSubProcesses.clear();
  allParticles.clear();
  theHandler = tcEventBasePtr();
  theColourLines.clear();
  theNumber = -1;
  theWeight = 0.0;
}

void Event::setInfo(tcEventBasePtr newHandler, string newName,
		    long newNumber, double weight) {
  theHandler = newHandler;
  name(newName);
  theNumber = newNumber;
  theWeight = weight;
}

double Event::optionalWeight(const string& name) const {
  map<string,double>::const_iterator w =
    theOptionalWeights.find(name);
  if ( w == theOptionalWeights.end() )
    return 0.;
  return w->second;
}

void Event::optionalWeight(const string& name, double value) {
  theOptionalWeights[name] = value;
}

tCollPtr Event::newCollision() {
  theCollisions.push_back(new_ptr(Collision(incoming(), this)));
  return theCollisions.back();
}

tStepPtr Event::newStep() {
  if ( theCollisions.empty() ) newCollision();
  return theCollisions.back()->newStep();
}


void Event::addCollision(tCollPtr c) {
  if ( !c ) return;
  theCollisions.push_back(c);
  addParticles(c->all().begin(), c->all().end());
  allSubProcesses.insert(c->subProcesses().begin(), c->subProcesses().end()); 
}

void Event::addParticle(tPPtr p) {
  if ( !p ) return;
  if ( member(allParticles, p) ) return;
  allParticles.insert(p);
  p->number(++theParticleNumber);
}

void Event::transform(const LorentzRotation & r) {
  for_each(allParticles, Transformer(r));
}

int Event::colourLineIndex(tcColinePtr line) const {
  ColourLineMap::const_iterator found = theColourLines.find(line);
  if ( found != theColourLines.end() ) return found->second;
  int index = theColourLines.size() + 1;
  theColourLines[line] = index;
  return index;
}

void Event::primaryCollision(tCollPtr c) {
  if ( !c ) return;
  if ( theCollisions.empty() ) theCollisions.push_back(c);
  else {
    if ( theCollisions[0] )
      for ( ParticleSet::const_iterator it = theCollisions[0]->all().begin();
	    it != theCollisions[0]->all().end(); ++it ) 
	allParticles.erase(*it);
    theCollisions[0] = c;
  }
  addParticles(c->all().begin(), c->all().end());
}

void Event::removeDecay(tPPtr p) {
  while ( !p->children().empty() ) removeParticle(p->children().back());
}

void Event::removeEntry(tPPtr p) {
  ParticleSet::iterator it = allParticles.find(p);
  if ( it == allParticles.end() ) return;
  for ( DIterator<CollisionVector::iterator> cit = theCollisions.begin();
	cit != theCollisions.end(); ++cit ) cit->removeEntry(p);
  allParticles.erase(it);
}

void Event::removeParticle(tPPtr p) {
  if ( p->next() ) removeParticle(p->next());
  while ( !p->children().empty() ) removeParticle(p->children().back());
  removeEntry(p);
}

void Event::cleanSteps() {
  for ( DIterator<CollisionVector::iterator> cit = theCollisions.begin();
	cit != theCollisions.end(); ++cit ) cit->cleanSteps();
}

EventPtr Event::clone() const {
  EventPtr newEvent = ptr_new<EventPtr>(*this);
  EventTranslationMap trans;
  trans[this] = newEvent;
  for ( CollisionVector::const_iterator cit = theCollisions.begin();
	cit != theCollisions.end(); ++cit ) trans[*cit] = (**cit).clone();
  for ( SubProcessSet::const_iterator spit = allSubProcesses.begin();
	spit != allSubProcesses.end(); ++spit ) trans[*spit] = (**spit).clone();
  for ( StepSet::const_iterator sit = allSteps.begin();
        sit != allSteps.end(); ++sit ) trans[*sit] = (**sit).clone();
  for ( ParticleSet::const_iterator pit = allParticles.begin();
	pit != allParticles.end(); ++pit )
    trans[*pit] = (**pit).clone();
  newEvent->rebind(trans);
  return newEvent;
}

void Event::rebind(const EventTranslationMap & trans) {
  theIncoming.first = trans.translate(theIncoming.first);
  theIncoming.second = trans.translate(theIncoming.second);
  for ( CollisionVector::iterator cit = theCollisions.begin();
	cit != theCollisions.end(); ++cit )
    (*cit = trans.translate(*cit))->rebind(trans);
  SubProcessSet newSubProcesses;
  for ( SubProcessSet::const_iterator spit = allSubProcesses.begin();
	spit != allSubProcesses.end(); ++spit )
    newSubProcesses.insert(trans.translate(*spit));
  allSubProcesses.swap(newSubProcesses);
  StepSet newSteps;
  for ( StepSet::const_iterator sit = allSteps.begin();
	sit != allSteps.end(); ++sit )
    newSteps.insert(trans.translate(*sit));
  allSteps.swap(newSteps);
  ParticleSet newParticles;
  for ( ParticleSet::const_iterator pit = allParticles.begin();
	pit != allParticles.end(); ++pit )
    newParticles.insert(trans.translate(*pit));
  allParticles.swap(newParticles);
}
  
ostream & ThePEG::operator<<(ostream & os, const Event & e) {
  os << string(78, '*') << endl
     << "Event number " << e.number() << " (id: " << e.name() << ") ";
  if ( e.handler() ) os << "performed by "
			<<  EventConfig::nameHandler(e.handler());
  os << endl;
  for ( unsigned int i = 0; i < e.collisions().size(); ++i ) {
    os << string(78, '=') << endl;
    if ( e.collisions().size() != 1 ) {
      if ( i ) {
	os << "Secondary Collision " << i;
	if ( e.collisions()[i]->handler() )
	  os << " performed by "
	     << EventConfig::nameHandler(e.collisions()[i]->handler());
      } else
	os << "Primary Collision";
    }
    os << endl << *e.collisions()[i];
  }
  return os;
}

// Helpers for the Graphviz output
namespace {
  static const string header = "digraph test {\nrankdir=LR;\nranksep=1.5;\n";

  inline long startnode(tcPPtr p) {
    return p->parents().empty() ? -p->uniqueId : p->parents()[0]->uniqueId;
  }

  inline long endnode(tcPPtr p) {
    return p->children().empty() ? p->uniqueId : startnode( p->children()[0] );
  }

  static const char * colours[] = {
    "red",
    "green",
    "blue",
    "orange",
    "aquamarine",
    "deeppink",
    "darkviolet",
    "darkolivegreen",
    "cyan"
  };

  struct Vertex {
    tcParticleSet in;
    tcParticleSet out;
  };

  typedef map<long, Vertex> VertexMap;

  template<typename Iter>
  Lorentz5Momentum sumP(Iter a, Iter b) {
    Lorentz5Momentum sum;
    for ( Iter it = a; it != b; ++it ) {
      sum += (*it)->momentum();
    }
    return sum;
  }
}

void ThePEG::Event::printGraphviz() const {
  ThePEG::printGraphviz(cout, this);
}


void ThePEG::printGraphviz(ostream & os, tcEventPtr ev) {
  os << header
     << "node [width=0.03,height=0.03,shape=point,label=\"\"];\n";

  tcParticleSet all;
  ev->select(inserter(all), SelectAll());

  VertexMap vertices;
  for (tcParticleSet::const_iterator it = all.begin();
       it != all.end(); ++it) {
    tcPPtr p = (*it);

    long start = startnode(p);
    long end = endnode(p);

    vertices[start].out.insert(p);
    vertices[end  ].in .insert(p);

    os << start << " -> " << end
       << " [label=\"" << p->number() << ' '
       << p->PDGName() << "\\n"
       << p->momentum().e()/GeV << "\\n"
       << p->momentum().mass()/GeV
       << "\"";

    if ( p->hasColourInfo() &&
	 ( p->colourLine() || p->antiColourLine() )) {
      os << ",penwidth=2,color=\"";

      const vector<tcColinePtr> & clines = p->colourInfo()->colourLines();
      for ( int i = 0, N = clines.size(); i < N; ++i ) {
	int colindex = ev->colourLineIndex(clines[i]) % 9;
	if ( i > 0 )  os << ':';
	os << colours[colindex];
      }
      const vector<tcColinePtr> & aclines = p->colourInfo()->antiColourLines();
      for ( int i = 0, N = aclines.size(); i < N; ++i ) {
	int colindex = ev->colourLineIndex(aclines[i]) % 9;
	if ( i > 0 || !clines.empty() )  os << ':';
	os << colours[colindex];
      }
      os << '"';
    }
    os << "];\n";
  }

  int label = 0;
  for ( VertexMap::const_iterator v = vertices.begin();
	v != vertices.end(); ++v ) {

    const long vertexId = v->first;
    const tcParticleSet & in  = v->second.in;
    const tcParticleSet & out = v->second.out;

    if ( in.empty() || out.empty() ) 
      continue;

    Lorentz5Momentum diff 
      = sumP(out.begin(), out.end()) - sumP(in.begin(), in.end());

    if ( abs(diff.e()) > 1.0*GeV ) {
      ++label;
      std::stringstream tail;
      tail << " [label=\"" << std::setprecision(4) 
	   << abs(diff.e()/GeV)
	   << " GeV\",arrowsize=3,penwidth=5,"
	   << "color=\"#ff000010\",fontcolor=\"#ff000010\"];\n";
      if ( diff.e() > ZERO )
	os << "mom" << label << " -> " << vertexId << tail.str();
      else
	os << vertexId << " -> " << "mom" << label << tail.str();
    }
  }

  os << '}' << endl;
}

void Event::debugme() const {
  cerr << *this;
  EventRecordBase::debugme();
}

void Event::persistentOutput(PersistentOStream & os) const {
  os << theIncoming << theCollisions << allSteps << allSubProcesses
     << allParticles << theNumber << theWeight << theOptionalWeights 
     << theParticleNumber;
  EventConfig::putHandler(os, theHandler);
}

void Event::persistentInput(PersistentIStream & is, int) {
  is >> theIncoming >> theCollisions >> allSteps >> allSubProcesses
     >> allParticles >> theNumber >> theWeight >> theOptionalWeights
     >> theParticleNumber;
  EventConfig::getHandler(is, theHandler);
}

ClassDescription<Event> Event::initEvent;

void Event::Init() {}

