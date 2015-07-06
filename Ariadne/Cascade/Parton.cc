// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Parton class.
//

#include "Parton.h"
#include "Junction.h"
#include "QCDDipole.h"
#include "Emission.h"
#include "DipoleState.h"
#include "AriadneHandler.h"
#include "Ariadne/Config/UnitFO.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/DescribeClass.h"

using namespace Ariadne5;

Parton::Parton(bool spec)
  : inFinalState(false), isSpecial(spec), theSystem(0), theOrigSystem(0) {}

Parton::~Parton() {}

ClonePtr Parton::clone() const {
  return new_ptr(*this);
}

void Parton::orig(tPPtr x) {
  theOrig = x;
  data(x->dataPtr());
  theMomentum = x->momentum();
  theVertex = x->vertex();
  theICol = x->antiColourLine();
  theOCol = x->colourLine();
}

void Parton::adoptOriginals(tParPtr orphan) {
  if ( orphan->orig() ) adoptedOriginals.push_back(orphan->orig());
  adoptedOriginals.insert(adoptedOriginals.end(),
			  orphan->adoptedOriginals.begin(), orphan->adoptedOriginals.end());
}
  

void Parton::emission(tcEmPtr e) {
  theEmission = e;
}

tPPtr Parton::produceParticle(const LorentzRotation & r) {
  theParticle = data().produceParticle(r*momentum());
  theParticle->setVertex(vertex());
  particle()->scale(sqr(coloured()? handler()->pTCut(): handler()->pTCutEM()));
  return particle();
}

void Parton::data(tcPDPtr x) {
  theDataPtr = x;
}

Energy2 Parton::invPT2() const {
  // *** ATTENTION***
  return ZERO;
}

void Parton::notify(const Emission &) {}

void Parton::fillReferences(CloneSet & cset) const {
  CascadeBase::fillReferences(cset);
  cset.insert(theEmission);
}

void Parton::rebind(const TranslationMap & trans) {
  CascadeBase::rebind(trans);
  theEmission = trans.translate(theEmission);
}

Parton::tParPair Parton::parents() const {
  tParPair ret;
  if ( emission() )
    ret = make_pair(emission()->colourParent, emission()->antiColourParent);
  return ret;
}

void Parton::persistentOutput(PersistentOStream & os) const {
  os << theOrig << theEmission << theParticle << theDataPtr
     << ounit(theMomentum, GeV) << ounit(theVertex, femtometer)
     << theICol << theOCol << inFinalState << isSpecial
     << theSystem << theOrigSystem << adoptedOriginals;
}

void Parton::persistentInput(PersistentIStream & is, int) {
  is >> theOrig >> theEmission >> theParticle >> theDataPtr
     >> iunit(theMomentum, GeV) >> iunit(theVertex, femtometer)
     >> theICol >> theOCol >> inFinalState >> isSpecial
     >> theSystem >> theOrigSystem >> adoptedOriginals;
}

DescribeClass<Parton,CascadeBase>
describeAriadne5Parton("Ariadne5::Parton", "libAriadne5.so");

void Parton::Init() {}

void Parton::debugme() const {
  CascadeBase::debugme();
  cerr << "P"
       << setw(3) << state()->index(this);
  if ( dataPtr() )
    cerr << setw(4) << data().id();
  else
    cerr << "*** NO DATA ***";
  cerr << (touched()? " *": "  ") << (emission()? "n": "o")
       << founit(momentum(), GeV, 10, 4);
}

