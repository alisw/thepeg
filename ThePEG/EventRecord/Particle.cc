// -*- C++ -*-
//
// Particle.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined member functions of
// the Particle class.
//
#include "Particle.h"
#include "ThePEG/EventRecord/Step.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/EventRecord/ColourLine.h"
#include "ThePEG/Utilities/Rebinder.h"
#include "ThePEG/Config/algorithm.h"
#include "ThePEG/EventRecord/ParticleTraits.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/DecayMode.h"
#include <iostream>
#include <iomanip>
#include <cctype>

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
#include "Particle.tcc"
#endif

using namespace ThePEG;

Particle::ParticleRep::ParticleRep(const ParticleRep & p)
  : theParents(p.theParents), theChildren(p.theChildren),
    thePrevious(p.thePrevious), theNext(p.theNext),
    theBirthStep(p.theBirthStep), theVertex(p.theVertex),
    theLifeLength(p.theLifeLength), theScale(p.theScale),
    theVetoScale(p.theVetoScale), theNumber(p.theNumber), 
    theExtraInfo(p.theExtraInfo.size()) {
  if ( p.theColourInfo )
    theColourInfo = dynamic_ptr_cast<CBPtr>(p.theColourInfo->clone());
  if ( p.theSpinInfo )
    theSpinInfo = dynamic_ptr_cast<SpinPtr>(p.theSpinInfo->clone());
  for ( int i = 0, N = p.theExtraInfo.size(); i < N; ++i )
    theExtraInfo[i] = p.theExtraInfo[i]->clone();
}

Particle::Particle(const Particle & p)
  : Base(p), theData(p.theData), theMomentum(p.theMomentum), theRep(p.theRep) {
  if ( p.theRep ) {
    theRep = new ParticleRep(*p.theRep);
    theRep->theParents.clear();
  }
}

Particle::~Particle() {
  if ( theRep ) {
    if ( colourLine() ) colourLine()->removeColoured(this);
    if ( antiColourLine() ) antiColourLine()->removeAntiColoured(this);
    delete theRep;
  }
  theRep = 0;
  theData = cEventPDPtr();
}

void Particle::initFull() {
  if ( theRep ) return;
  theRep = new ParticleRep;

  Energy width = data().generateWidth(mass());

  if ( width > ZERO ) {
    Time lifetime = data().generateLifeTime(mass(), width);
    theRep->theLifeLength.setTau(lifetime);
    theRep->theLifeLength.
      setVect((momentum().vect()*(lifetime /
				max(mass(), Constants::epsilon*GeV))));
    theRep->theLifeLength.rescaleEnergy();
  }
}

PPtr Particle::clone() const {
  return ptr_new<PPtr>(*this);
}

void Particle::rebind(const EventTranslationMap & trans) {
  for ( ParticleVector::iterator pit = rep().theChildren.begin();
	pit != rep().theChildren.end(); ++pit ) *pit = trans.translate(*pit);
  for ( tParticleVector::iterator pit = rep().theParents.begin();
	pit != rep().theParents.end(); ++pit ) *pit = trans.translate(*pit);
  rep().thePrevious = trans.translate(rep().thePrevious);
  rep().theNext = trans.translate(rep().theNext);
  if ( hasColourInfo() ) colourInfo()->rebind(trans);
  if ( spinInfo() ) spinInfo()->rebind(trans);
  rep().theBirthStep = trans.translate(rep().theBirthStep);
  for ( EIVector::const_iterator ie = rep().theExtraInfo.begin();
	ie != rep().theExtraInfo.end(); ++ie ) (**ie).rebind(trans);
}

tParticleSet Particle::siblings() const {
  tParticleSet theSiblings;
  for ( tParticleVector::const_iterator pit = parents().begin();
	pit != parents().end(); ++pit )
    theSiblings.insert((*pit)->children().begin(), (*pit)->children().end());
  theSiblings.erase(const_cast<Particle *>(this));
  return theSiblings;
}

void Particle::colourNeighbour(tPPtr p, bool anti) {
  tColinePtr line = colourLine(!anti);
  if ( !line ) line = ColourLine::create(this, !anti);
  line->addColoured(p, anti);  
}

void Particle::outgoingColour(tPPtr p, bool anti) {
  tColinePtr line = colourLine(anti);
  if ( !line ) line = ColourLine::create(this, anti);
  line->addColoured(p, anti);  
}

tPPtr Particle::incomingColour(bool anti) const {
  if ( !hasColourInfo() ) return tPPtr();
  tColinePtr line = colourLine(anti);
  if ( !line ) return tPPtr();
  for ( int i = 0, N = parents().size(); i < N; ++i )
    if ( parents()[i]->hasColourLine(line, anti) ) return parents()[i];
  return tPPtr();
}

tPPtr Particle::outgoingColour(bool anti) const {
  if ( !hasColourInfo() ) return tPPtr();
  tColinePtr line = colourLine(anti);
  if ( !line ) return tPPtr();
  for ( int i = 0, N = children().size(); i < N; ++i )
    if ( children()[i]->hasColourLine(line, anti) ) return children()[i];
  return tPPtr();
}

LorentzPoint Particle::labVertex() const {
  LorentzPoint r(rep().theBirthStep && rep().theBirthStep->collision()?
		 vertex() + rep().theBirthStep->collision()->vertex():
		 vertex());
  return r;
}

void Particle::setLabVertex(const LorentzPoint & p) {
  rep().theVertex = ( rep().theBirthStep && rep().theBirthStep->collision()?
		      p - rep().theBirthStep->collision()->vertex() : p );
}

void Particle::transform(const LorentzRotation & r) {
  if ( hasRep() && spinInfo() ) spinInfo()->transform(momentum(), r);
  theMomentum.transform(r);
  if ( !hasRep() ) return;
  rep().theVertex.transform(r);
  rep().theLifeLength.transform(r);
}

void Particle::deepTransform(const LorentzRotation & r) {
  transform(r);
  if ( !theRep ) return;
  for ( int i = 0, N = children().size(); i < N; ++i )
    rep().theChildren[i]->deepTransform(r);
  if ( rep().theNext ) rep().theNext->deepTransform(r);
}

void Particle::rotateX(double a) {
  LorentzRotation r;
  r.rotateX(a);
  transform(r);
}

void Particle::deepRotateX(double a) {
  LorentzRotation r;
  r.rotateX(a);
  deepTransform(r);
}

void Particle::rotateY(double a) {
  LorentzRotation r;
  r.rotateY(a);
  transform(r);
}

void Particle::deepRotateY(double a) {
  LorentzRotation r;
  r.rotateY(a);
  deepTransform(r);
}

void Particle::rotateZ(double a) {
  LorentzRotation r;
  r.rotateZ(a);
  transform(r);
}

void Particle::deepRotateZ(double a) {
  LorentzRotation r;
  r.rotateZ(a);
  deepTransform(r);
}

void Particle::rotate(double a, const Axis & axis) {
  LorentzRotation r;
  r.rotate(a, axis);
  transform(r);
}

void Particle::deepRotate(double a, const Axis & axis) {
  LorentzRotation r;
  r.rotate(a, axis);
  deepTransform(r);
}

string Particle::outputFormat =
"%n3%s10 %i7 %p[,]0 %c(,) %^^0%vv0 %>>0%<>0 %l{,}0\n"
"                            %x10.3%y10.3%z10.3%e10.3%m10.3\n";

int getNumber(string::const_iterator & pos, int def) {
  if ( !isdigit(*pos) ) return def;
  def = *pos++ - '0';
  while ( isdigit(*pos) ) def = 10*def + *pos++ - '0';
  return def;
}

void writePrecision(ostream & os, string::const_iterator & pos,
			      int defw, int defp, double x) {
  defw = getNumber(pos, defw);
  if ( *pos == '.' ) defp = getNumber(++pos, defp);
  int oldp = os.precision();
  os << setprecision(defp) << setw(defw) << x << setprecision(oldp);
}

void writeStringAdjusted(ostream & os, bool left, int w, string str) {
  while ( !left && w-- > int(str.size()) ) os << ' ';
  os << str;
  while ( left && w-- > int(str.size()) ) os << ' ';
}

template <typename Container>
void writeParticleRanges(ostream & os, const Container & co, char sep, int w) {
  set<int> cnum;
  for ( typename Container::const_iterator it = co.begin();
	it != co.end(); ++it) cnum.insert((**it).number());

  bool elipsis = false;
  int last = -10;
  for ( set<int>::iterator it = cnum.begin(); it != cnum.end(); ++it) {
    int n = *it;
    int next = 0;
    set<int>::iterator itn = it;
    if ( ++itn != cnum.end() ) next = *itn;
    bool writeit = true;
    bool writesep = false;
    if ( elipsis && ( n != last + 1 || n != next - 1 ) )
      elipsis = false;
    else if ( !elipsis && n == last + 1 && n == next -1 ) {
      os << "..";
      elipsis = true;
      writeit = false;
    }
    else if ( elipsis && n == last + 1 && n == next -1 )
      writeit = false;
    else if ( it != cnum.begin() )
      writesep = true;
    if ( writeit ) {
      if ( writesep ) os << sep;
      os << setw(w) << n;
    }
    last = n;
  }
}

ostream & ThePEG::operator<<(ostream & os, const Particle & p) {
  return p.print(os, p.birthStep());
}

ostream & Particle::print(ostream & os, tcStepPtr step) const {
  if ( !step ) step = birthStep();
  tCollPtr coll = step? step->collision(): tCollPtr();
  tEventPtr event = coll? coll->event(): tEventPtr();
  string::const_iterator pos = Particle::outputFormat.begin();
  ios::fmtflags saveflags = os.setf(ios::fixed, ios::floatfield);
  while ( pos != Particle::outputFormat.end() ) {
    if ( *pos == '%' && ++pos != Particle::outputFormat.end() ) {
      bool left = false;
      if ( *pos == '-' ) {
	left = true;
	os.setf(ios::left, ios::adjustfield);
	++pos;
      } else {
	os.setf(ios::right, ios::adjustfield);
      }
      char mark;
      char open;
      char close;
      char sep;
      int w;
      string str;
      string fill;
      if ( pos == Particle::outputFormat.end() ) break;
      bool fullColour = false;
      switch ( *pos ) {
      case 'n':
	os << setw(getNumber(++pos, 3)) << number();
	break;
      case 'i':
	os << setw(getNumber(++pos, 8)) << id();
	break;
      case 's':
	writeStringAdjusted(os, left, getNumber(++pos, 8), PDGName());
	break;
      case 'x':
	writePrecision(os, ++pos, 10, 3, momentum().x()/GeV);
	break;
      case 'y':
	writePrecision(os, ++pos, 10, 3, momentum().y()/GeV);
	break;
      case 'z':
	writePrecision(os, ++pos, 10, 3, momentum().z()/GeV);
	break;
      case 'e':
	writePrecision(os, ++pos, 10, 3, momentum().e()/GeV);
	break;
      case 'm':
	writePrecision(os, ++pos, 10, 3, momentum().mass()/GeV);
	break;
      case 'P':
	fullColour = true;
      case 'p':
	open = *++pos;
	sep = *++pos;
	close = *++pos;
	w = getNumber(++pos, 0);
	if ( parents().empty() ) break;
	if ( open ) os << open;
	writeParticleRanges(os, parents(), sep, w);
	if ( fullColour && hasColourInfo() &&
	     ( incomingColour() || incomingAntiColour() ) ) {
	  if ( close ) os << open;
	  if ( incomingColour() )
	    os << "+" << incomingColour()->number();
	  if ( incomingAntiColour() )
	    os << "-" << incomingAntiColour()->number();
	  if ( close ) os << close;
	}
	if ( close ) os << close;
	break;
      case 'l':
	open = *++pos;
	sep = *++pos;
	close = *++pos;
	w = getNumber(++pos, 0);
	if ( hasColourInfo() &&
	     ( colourLine() || antiColourLine() ) && event) {
	  if ( open ) os << open;
	  vector<tcColinePtr> clines = colourInfo()->colourLines();
	  for ( int i = 0, N = clines.size(); i < N; ++i ) {
	    if ( i > 0 && sep )  os << sep;
	    clines[i]->write(os, event, false);
	  }
	  vector<tcColinePtr> aclines = colourInfo()->antiColourLines();
	  for ( int i = 0, N = aclines.size(); i < N; ++i ) {
	    if ( ( i > 0 || clines.size() ) && sep )  os << sep;
	    aclines[i]->write(os, event, true);
	  }
	  if ( close ) os << close;
	}
	break;
      case 'C':
	fullColour = true;
      case 'c':
	open = *++pos;
	sep = *++pos;
	close = *++pos;
	w = getNumber(++pos, 0);
	if ( children().empty() ) break;
	if ( open ) os << open;
	writeParticleRanges(os, children(), sep, w);
	if ( fullColour && hasColourInfo() &&
	     ( outgoingColour() || outgoingAntiColour() ) ) {
	  if ( close ) os << open;
	  if ( outgoingColour() )
	    os << "+" << outgoingColour()->number();
	  if ( outgoingAntiColour() )
	    os << "-" << outgoingAntiColour()->number();
	  if ( close ) os << close;
	}
	if ( close ) os << close;
	break;
      case '>':
	mark = *++pos;
	w = getNumber(++pos, 0);
	if ( hasColourInfo() && step && step->colourNeighbour(this) ) {
	  os << setw(w-1) << step->colourNeighbour(this)->number() << mark;
	}
	break;
      case '<':
	mark = *++pos;
	w = getNumber(++pos, 0);
	if ( hasColourInfo() && step && step->antiColourNeighbour(this) ) {
	  int n = step->antiColourNeighbour(this)->number();
	  ostringstream oss;
	  oss << mark << n;
	  writeStringAdjusted(os, left, w, oss.str());
	}
	break;
      case 'v':
	mark = *++pos;
	w = getNumber(++pos, 0);
	if ( next() ) {
	  if ( left && mark ) os << mark;
	  os << setw(w) << next()->number();
	  if ( !left && mark ) os << mark;
	}
	break;
      case '^':
	mark = *++pos;
	w = getNumber(++pos, 0);
	if ( previous() ) {
	  if ( left && mark ) os << mark;
	  os << setw(w) << previous()->number();
	  if ( !left && mark ) os << mark;
	}
	break;
      case 'd':
	switch ( *++pos ) {
	case 'x':
	  writePrecision(os, ++pos, 10, 3, lifeLength().x()/mm);
	  break;
	case 'y':
	  writePrecision(os, ++pos, 10, 3, lifeLength().y()/mm);
	  break;
	case 'z':
	  writePrecision(os, ++pos, 10, 3, lifeLength().z()/mm);
	  break;
	case 't':
	  writePrecision(os, ++pos, 10, 3, lifeLength().e()/mm);
	  break;
	case 'T':
	  writePrecision(os, ++pos, 10, 3, lifeLength().tau()/mm);
	  break;
	}
	break;
      case 'V':
	switch ( *++pos ) {
	case 'x':
	  writePrecision(os, ++pos, 10, 3, vertex().x()/mm);
	  break;
	case 'y':
	  writePrecision(os, ++pos, 10, 3, vertex().y()/mm);
	  break;
	case 'z':
	  writePrecision(os, ++pos, 10, 3, vertex().z()/mm);
	  break;
	case 't':
	  writePrecision(os, ++pos, 10, 3, vertex().e()/mm);
	  break;
	}
      case 'L':
	switch ( *++pos ) {
	case 'x':
	  writePrecision(os, ++pos, 10, 3, labVertex().x()/mm);
	  break;
	case 'y':
	  writePrecision(os, ++pos, 10, 3, labVertex().y()/mm);
	  break;
	case 'z':
	  writePrecision(os, ++pos, 10, 3, labVertex().z()/mm);
	  break;
	case 't':
	  writePrecision(os, ++pos, 10, 3, labVertex().e()/mm);
	  break;
	}
	break;
      default:
	os << *pos++;
      }
    } else {
      if ( pos != Particle::outputFormat.end() ) os << *pos++;
    }
  }
  os.flags(saveflags);
  return os;
}

void Particle::debugme() const {
  cerr << *this;
  EventRecordBase::debugme();
}

void Particle::persistentOutput(PersistentOStream & os) const {
  EventConfig::putParticleData(os, theData);
  os << ounit(theMomentum, GeV) << bool( theRep != 0 );
  if ( !theRep ) return;
  os << rep().theParents << rep().theChildren
     << rep().thePrevious << rep().theNext << rep().theBirthStep
     << ounit(rep().theVertex, mm) << ounit(rep().theLifeLength, mm)
     << ounit(rep().theScale, GeV2) << ounit(rep().theVetoScale, GeV2) 
     << rep().theNumber << rep().theDecayMode
     << rep().theColourInfo << rep().theSpinInfo << rep().theExtraInfo;
}

void Particle::persistentInput(PersistentIStream & is, int) {
  bool readRep;
  EventConfig::getParticleData(is, theData);
  is >> iunit(theMomentum, GeV) >> readRep;
  if ( !readRep ) return;
  if ( !hasRep() ) theRep = new ParticleRep;

  is >> rep().theParents >> rep().theChildren
     >> rep().thePrevious >> rep().theNext >> rep().theBirthStep
     >> iunit(rep().theVertex, mm) >> iunit(rep().theLifeLength, mm)
     >> iunit(rep().theScale, GeV2) >> iunit(rep().theVetoScale, GeV2) 
     >> rep().theNumber >> rep().theDecayMode
     >> rep().theColourInfo >> rep().theSpinInfo >> rep().theExtraInfo;
}

ClassDescription<Particle> Particle::initParticle;

void Particle::Init() {}

ThePEG_IMPLEMENT_SET(PPtr,ParticleSet)

