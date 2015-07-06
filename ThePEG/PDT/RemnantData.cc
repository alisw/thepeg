// -*- C++ -*-
//
// RemnantData.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RemnantData class.
//

#include "RemnantData.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/RemnantDecayer.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;

RemnantData::
RemnantData(tcPDPtr particle, RemDecPtr dec)
  //  : ParticleData(particle->id(), "Rem:" + particle->PDGName()),
  : ParticleData(82, "Rem:" + particle->PDGName()),
    parentPD(particle), theDecayer(dec) {
  synchronized(false);
  width(ZERO);
  widthCut(ZERO);
  iCharge(particle->iCharge());
  iSpin(particle->iSpin());
  iColour(particle->iColour());
  stable(true);
  fixColour();
}

IBPtr RemnantData::clone() const {
  return new_ptr(*this);
}

IBPtr RemnantData::fullclone() const {
  return new_ptr(*this);
}

bool RemnantData::extract(tcPDPtr parton) {
  if ( !decayMode ) {
    stable(false);
    decayMode = new_ptr(DecayMode());
    decayMode->parent(this);
    decayMode->brat(1.0);
    decayMode->decayer(theDecayer);
    decayMode->switchOn();
    addDecayMode(decayMode);
  }
  if ( !parton ) return false;
  //  if ( !extracted.empty() && !decayer->multiCapable() ) return false;
  if ( !theDecayer->canHandle(parentPD, parton) ) return false;
  extracted.insert(parton);
  iCharge(PDT::Charge(iCharge() - parton->iCharge()));
  return fixColour();
}

bool RemnantData::
reextract(tcPDPtr oldp, tcPDPtr newp) {
  multiset<tcPDPtr>::iterator it = extracted.find(oldp);
  if ( it == extracted.end() ) return false;
  if ( !theDecayer->canHandle(parentPD, newp) ) return false;
  extracted.erase(it);
  extracted.insert(newp);
  iCharge(PDT::Charge(iCharge() + oldp->iCharge() - newp->iCharge()));
  return fixColour();
}

bool RemnantData::
remove(tcPDPtr oldp) {
  multiset<tcPDPtr>::iterator it = extracted.find(oldp);
  if ( it == extracted.end() ) return false;
  extracted.erase(it);
  iCharge(PDT::Charge(iCharge() + oldp->iCharge()));
  return fixColour();
}

bool RemnantData::fixColour() {
  bool col = parentPD->hasAntiColour();
  bool acol = parentPD->hasColour();
  for ( multiset<tcPDPtr>::const_iterator it = extracted.begin();
	it != extracted.end(); ++it ) {
    if ( (**it).hasColour() ) col = true;
    if ( (**it).hasAntiColour() ) acol = true;
  }
  if ( col )
    iColour(acol? PDT::Colour8: PDT::Colour3bar);
  else
    iColour(acol? PDT::Colour3: PDT::Colour0);

  return true;
}

void RemnantData::persistentOutput(PersistentOStream & os) const {
  os << parentPD << theDecayer << decayMode << extracted;
}

void RemnantData::persistentInput(PersistentIStream & is, int) {
  is >> parentPD >> theDecayer >> decayMode >> extracted;
}

ClassDescription<RemnantData> RemnantData::initRemnantData;
// Definition of the static class description member.

void RemnantData::Init() {

  static ClassDocumentation<RemnantData> documentation
    ("There is no documentation for the RemnantData class");

}

