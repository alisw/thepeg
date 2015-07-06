// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QCDDipole class.
//

#include "QCDDipole.h"
#include "Junction.h"
#include "RemnantParton.h"
#include "ResonanceParton.h"
#include "DipoleState.h"
#include "AriadneHandler.h"
#include "Emission.h"

#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Utilities/Current.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/EventRecord/Particle.h"

using namespace Ariadne5;

QCDDipole::QCDDipole() {
  theRhoCut = Current<AriadneHandler>()->pTCut();
}

ClonePtr QCDDipole::clone() const {
  return new_ptr(*this);
}

ColinePtr QCDDipole::colourLine() const {
  tJunctionPtr ji = dynamic_ptr_cast<tJunctionPtr>(iPart());
  tJunctionPtr jo = dynamic_ptr_cast<tJunctionPtr>(oPart());
  // If both partons are junctions then we have no colour line to report.
  if ( ji && jo ) return ColinePtr();

  // First check if we are attached to an original colour line.
  ColinePtr cl;
  // This must always be the case if we have a juction.
  if ( ji || jo ) {
    if ( ji ) cl = oPart()->origICol();
    if ( jo ) cl = iPart()->origOCol();
    if ( !cl )
      Throw<Exception>()
	<< "Could not find colour line of Junction in Ariadne. This is a "
	<< "serious error, please contact the author." << Exception::runerror;
    return cl;
  } else {
  // Now simply check the partons and see if there is already a colour line
    cl = iPart()->origOCol();
    if ( oPart()->origICol() ) {
      if ( cl && oPart()->origICol() != cl )
	Throw<Exception>()
	  << "Found Conflicting colour lines in Ariadne. This is a "
	  << "serious error, please contact the author." << Exception::runerror;
      cl = oPart()->origICol();
    }
    // If no colour line was found we may create one.
    if ( !cl ) cl = new_ptr(ColourLine());

    // Now connect the produced particles.
    tRemParPtr ri = dynamic_ptr_cast<tRemParPtr>(iPart());
    if ( ri && ri->hard() ) ri = tRemParPtr();
    tRemParPtr ro = dynamic_ptr_cast<tRemParPtr>(oPart());
    if ( ro && ro->hard() ) ro = tRemParPtr();

    if ( ri ) cl->addAntiColoured(ri->extracted());
    else if ( !ji ) cl->addColoured(iPart()->particle());
    if ( ro ) cl->addColoured(ro->extracted());
    else if ( !jo ) cl->addAntiColoured(oPart()->particle());

  }
  return cl;
}

Energy2 QCDDipole::sdip() const {
  return (iPart()->momentum() + oPart()->momentum()).m2();
}

void QCDDipole::generateColourIndex() {
  theColourIndex.generate(next()? next()->colourIndex(): ColourIndex(),
			  prev()? prev()->colourIndex(): ColourIndex());
}

vector<tParPtr> QCDDipole::string() const {
  vector<tParPtr> ret;
  pair<tcQCDPtr,tcQCDPtr> range = DipoleState::StringEnds(this);
  tParPtr p1 = range.first->iPart();
  tParPtr p2 = range.second->oPart();
  while ( range.first ) {
    ret.push_back(range.first->iPart());
    if ( range.first == range.second )
      range.first = tcQCDPtr();
    else
      range.first = range.first->next();
  }
  if ( p1 != p2 ) ret.push_back(p2);
  return ret;
}


void QCDDipole::fillReferences(CloneSet & cset) const {
  DipoleBase::fillReferences(cset);
  cset.insert(oPart());
  cset.insert(iPart());
}

void QCDDipole::rebind(const TranslationMap & trans) {
  DipoleBase::rebind(trans);
  theOPart = trans.translate(oPart());
  theIPart = trans.translate(iPart());
  theNext =  trans.translate(theNext);
  thePrev =  trans.translate(thePrev);
}

void QCDDipole::persistentOutput(PersistentOStream & os) const {
  os << theIPart << theOPart << theNext << thePrev 
     << theColourIndex << theResonance;
}

void QCDDipole::persistentInput(PersistentIStream & is, int) {
  is >> theIPart >> theOPart >> theNext >> thePrev
     >> theColourIndex >> theResonance;
}

DescribeClass<QCDDipole,DipoleBase>
describeAriadne5QCDDipole("Ariadne5::QCDDipole", "libAriadne5.so");

void QCDDipole::Init() {}

void QCDDipole::debugme() const {
  DipoleBase::debugme();
  cerr << "D"
       << setw(3) << state()->index(this) << " "
       << setw(3) << state()->index(iPart()) << " "
       << setw(3) << state()->index(oPart())
       << setw(10) << setprecision(3)
       << (touched()? " *": "  ")
       << "[" << colourIndex() << "]";
}

bool QCDDipole::checkIntegrity() {
  if ( oPart()->origICol() && iPart()->origOCol() &&
       oPart()->origICol() != iPart()->origOCol() )
    return false;
  return DipoleBase::checkIntegrity() &&
    iPart() && (!next() || next()->prev() == this) &&
    oPart() && (!prev() || prev()->next() == this);
}

