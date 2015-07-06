// -*- C++ -*-
//
// ParticleData.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ParticleData class.
//

#include "ParticleData.h"
#include "ParticleData.xh"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/Utilities/HoldFlag.h"
#include "ThePEG/Utilities/Rebinder.h"
#include "ThePEG/Utilities/StringUtils.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Repository/Repository.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Config/algorithm.h"
#include "ThePEG/Utilities/Exception.h"
#include "ThePEG/Utilities/EnumIO.h"
#include "ThePEG/Repository/UseRandom.h"

namespace ThePEG {

ParticleData::ParticleData()
  : theId(0), thePDGName(""), theMass(-1.0*GeV), theWidth(-1.0*GeV),
    theWidthUpCut(-1.0*GeV), theWidthLoCut(-1.0*GeV), theCTau(-1.0*mm),
    theCharge(PDT::ChargeUnknown),
    theSpin(PDT::SpinUnknown), theColour(PDT::ColourUnknown), isStable(true),
    theVariableRatio(false), syncAnti(false), theDefMass(-1.0*GeV),
    theDefWidth(-1.0*GeV), theDefCut(-1.0*GeV), theDefCTau(-1.0*mm),
    theDefCharge(PDT::ChargeUnknown), theDefSpin(PDT::SpinUnknown),
    theDefColour(PDT::ColourUnknown) {}

ParticleData::
ParticleData(PID newId, const string & newPDGName)
  : theId(newId), thePDGName(newPDGName), theMass(-1.0*GeV), theWidth(-1.0*GeV),
    theWidthUpCut(-1.0*GeV), theWidthLoCut(-1.0*GeV), theCTau(-1.0*mm),
    theCharge(PDT::ChargeUnknown),
    theSpin(PDT::SpinUnknown), theColour(PDT::ColourUnknown), isStable(true),
    theVariableRatio(false), syncAnti(false), theDefMass(-1.0*GeV),
    theDefWidth(-1.0*GeV), theDefCut(-1.0*GeV), theDefCTau(-1.0*mm),
    theDefCharge(PDT::ChargeUnknown), theDefSpin(PDT::SpinUnknown),
    theDefColour(PDT::ColourUnknown) {}

ParticleData::~ParticleData() {}

PDPtr ParticleData::Create(PID newId, const string & newPDGName) {
  return new_ptr(ParticleData(newId, newPDGName));
}

PDPair ParticleData::
Create(PID newId, const string & newPDGName, const string & newAntiPDGName) {
  PDPair pap;
  pap.first = new_ptr(ParticleData(newId, newPDGName));
  pap.second = new_ptr(ParticleData(-newId, newAntiPDGName));
  antiSetup(pap);
  return pap;
}

void ParticleData::readSetup(istream & is) {
  long id;
  is >> id >> thePDGName >> iunit(theDefMass, GeV) >> iunit(theDefWidth, GeV)
     >> iunit(theDefCut, GeV) >> iunit(theDefCTau, mm) >> ienum(theDefCharge)
     >> ienum(theDefColour) >> ienum(theDefSpin) >> ienum(isStable);
  theId = id;
  theMass = theDefMass;
  theWidth = theDefWidth;
  theWidthUpCut = theDefCut;
  theWidthLoCut = theDefCut;
  theCTau = theDefCTau;
  theCharge = theDefCharge;
  theColour = theDefColour;
  theSpin = theDefSpin;
  if ( PDGName() == "-" ) thePDGName = name();
  return;
}
  
void ParticleData::antiSetup(const PDPair & pap) {
  pap.first->theAntiPartner = pap.second;
  pap.second->theAntiPartner = pap.first;
  pap.first->syncAnti = pap.second->syncAnti = true;
}

PDPtr ParticleData::pdclone() const {
  return new_ptr(*this);
}

IBPtr ParticleData::clone() const {
  return pdclone();
}
  
IBPtr ParticleData::fullclone() const {
  PDPtr pd = pdclone();
  Repository::Register(pd);
  pd->theDecaySelector.clear();
  pd->theDecayModes.clear();
  pd->isStable = true;
  PDPtr apd;
  if ( CC() ) {
    apd = CC()->pdclone();
    Repository::Register(apd);
    apd->theDecaySelector.clear();
    apd->theDecayModes.clear();
    apd->isStable = true;
    pd->theAntiPartner = apd;
    apd->theAntiPartner = pd;
    pd->syncAnti = syncAnti;
    apd->syncAnti = CC()->syncAnti;
  }
  HoldFlag<> dosync(pd->syncAnti, true);
  for ( DecaySet::const_iterator it = theDecayModes.begin();
	it != theDecayModes.end(); ++it )
    pd->addDecayMode(*it);
  return pd;
}

Energy ParticleData::width(Energy wi) {
  theWidth = wi;
  if ( synchronized() && CC() ) CC()->theWidth = theWidth;
  return theWidth;
}

Energy ParticleData::widthUpCut(Energy wci) {
  theWidthUpCut = wci;
  if ( synchronized() && CC() ) CC()->theWidthUpCut = theWidthUpCut;
  return theWidthUpCut;
}

Energy ParticleData::widthLoCut(Energy wci) {
  theWidthLoCut = wci;
  if ( synchronized() && CC() ) CC()->theWidthLoCut = theWidthLoCut;
  return theWidthLoCut;
}

Length ParticleData::cTau(Length ti) {
  theCTau = ti;
  if ( synchronized() && CC() ) CC()->theCTau = theCTau;
  return theCTau;
}

PDT::Charge ParticleData::iCharge(PDT::Charge ci) {
  theCharge = ci;
  if ( synchronized() && CC() ) CC()->theCharge = PDT::Charge(-ci);
  return theCharge;
}

PDT::Spin ParticleData::iSpin(PDT::Spin si) {
  theSpin = si;
  if ( synchronized() && CC() ) CC()->theSpin = si;
  return si;
}

PDT::Colour ParticleData::iColour(PDT::Colour ci) {
  theColour = ci;
  if ( synchronized() && CC() ) CC()->theColour = PDT::Colour(-ci);
  return theColour;
}

void ParticleData::stable(bool s) {
  isStable = s;
  if ( synchronized() && CC() ) CC()->isStable = s;
}

void ParticleData::synchronized(bool h) {
  syncAnti = h;
  if ( CC() ) CC()->syncAnti = h;
}

void ParticleData::variableRatio(bool varRatio) {
  theVariableRatio=varRatio;
}

void ParticleData::addDecayMode(tDMPtr dm) {
  if ( member(theDecayModes, dm) ) return;
  cPDPtr parent = dm->parent();
  if ( !parent ) parent = this;
  if ( parent != this ) {
    dm = dm->clone(this);
  }
  theDecayModes.insert(dm);
  theDecaySelector.insert(dm->brat(), dm);
  if ( CC() ) {
    if ( !synchronized() ) dm->CC()->switchOff();
    CC()->theDecayModes.insert(dm->CC());
    CC()->theDecaySelector.insert(dm->CC()->brat(), dm->CC());
  }
}

void ParticleData::removeDecayMode(tDMPtr dm) {
  theDecayModes.erase(theDecayModes.find(dm));
  if(theDecayModes.empty()) isStable = true;
  theDecaySelector.erase(dm);
  if ( !CC() ) return;
  CC()->theDecayModes.erase(dm->CC());
  if(CC()->theDecayModes.empty()) CC()->isStable = true;
  CC()->theDecaySelector.erase(dm->CC());
}

void ParticleData::synchronize() {
  if ( !CC() ) return;
  isStable = CC()->isStable;
  theMass = CC()->theMass;
  theWidth = CC()->theWidth;
  theWidthUpCut = CC()->theWidthUpCut;
  theWidthLoCut = CC()->theWidthLoCut;
  theCTau = CC()->theCTau;
  theCharge = PDT::Charge(-CC()->theCharge);
  theSpin = CC()->theSpin;
  theColour = PDT::antiColour(CC()->theColour);
  theMassGenerator = CC()->theMassGenerator;
  theWidthGenerator = CC()->theWidthGenerator;
  syncAnti = CC()->syncAnti;
  theDecaySelector.clear();
  for ( DecaySet::iterator it = theDecayModes.begin();
	it != theDecayModes.end(); ++it ) {
    (*it)->synchronize();
    theDecaySelector.insert((*it)->brat(), *it);
  }
}

void ParticleData::doupdate() {
  Interfaced::doupdate();
  bool redo = touched();
  for_each(theDecayModes, UpdateChecker(redo));
  UpdateChecker::check(theMassGenerator, redo);
  UpdateChecker::check(theWidthGenerator, redo);
  if ( !redo ) return;

  theDecaySelector.clear();
  for ( DecaySet::const_iterator dit = theDecayModes.begin();
	dit != theDecayModes.end(); ++dit ) {
    tDMPtr dm = *dit;
    dm->resetOverlap();
    for ( DecaySet::const_iterator dit2 = theDecayModes.begin();
	  dit2 != theDecayModes.end(); ++dit2 )
      if ( dit2 != dit ) dm->addOverlap(dm);
    if ( dm->brat() > 0.0 ) theDecaySelector.insert(dm->brat(), dm);
  }
  if ( theMassGenerator && !theMassGenerator->accept(*this) )
    throw UpdateException();
  if ( theWidthGenerator &&
       !theWidthGenerator->accept(*this) )
    throw UpdateException();
  if ( theWidthGenerator ) theDecaySelector = theWidthGenerator->rate(*this);
  touch();
}

tDMPtr ParticleData::selectMode(Particle & p) const {
  if ( &(p.data()) != this ) return tDMPtr();
  try {
    if ( !theWidthGenerator || !theVariableRatio )
      return theDecaySelector.select(UseRandom::current());
    DecaySelector local;
    if ( theWidthGenerator )
      local = theWidthGenerator->rate(p);
    else
      for ( DecaySet::const_iterator mit = theDecayModes.begin();
	    mit != theDecayModes.end(); ++mit  )
	local.insert((*mit)->brat(p), *mit);
    return local.select(UseRandom::current());
  }
  catch (range_error) {
    return tDMPtr();
  }
}

void ParticleData::rebind(const TranslationMap & trans) {
  if ( CC() ) theAntiPartner = trans.translate(theAntiPartner);
  DecaySet newModes;
  DecaySelector newSelector;
  for ( DecaySet::iterator it = theDecayModes.begin();
	it != theDecayModes.end(); ++it ) {
    DMPtr dm;
    dm = trans.translate(*it);
    if ( !dm ) throw RebindException();
    newModes.insert(dm);
    newSelector.insert(dm->brat(), dm);
  }
  theDecayModes.swap(newModes);
  theDecaySelector.swap(newSelector);
}

IVector ParticleData::getReferences() {
  IVector refs = Interfaced::getReferences();
  if ( CC() ) refs.push_back(CC());
  refs.insert(refs.end(), theDecayModes.begin(), theDecayModes.end());
  return refs;
}



void ParticleData::massGenerator(tMassGenPtr mg) {
  if ( mg && !mg->accept(*this) ) return;
  if ( mg && synchronized() && CC() && !mg->accept(*CC()) ) return;
  theMassGenerator = mg;
  if ( synchronized() && CC() ) CC()->theMassGenerator = mg;
}

void ParticleData::widthGenerator(tWidthGeneratorPtr newGen) {
  if ( newGen && !newGen->accept(*this) ) return;
  if ( newGen && synchronized() && CC() && !newGen->accept(*CC()) ) return;
  theWidthGenerator = newGen;
  if ( synchronized() && CC() ) CC()->theWidthGenerator = newGen;
}

Energy ParticleData::generateMass() const {
  return massGenerator()? massGenerator()->mass(*this): mass();
}

Energy ParticleData::generateWidth(Energy m) const {
  return widthGenerator()? widthGenerator()->width(*this, m): width();
}

Length ParticleData::generateLifeTime(Energy m, Energy w) const {
  return widthGenerator() ? 
    widthGenerator()->lifeTime(*this, m, w) :
    UseRandom::rndExp(cTau());
}

PPtr ParticleData::produceParticle(const Lorentz5Momentum & pp) const {
  PPtr p = new_ptr(Particle(this));
  p->set5Momentum(pp);
  return p;
}

PPtr ParticleData::produceParticle(const LorentzMomentum & pp) const {
  PPtr p(produceParticle(Lorentz5Momentum(pp)));
  return p;
}

PPtr ParticleData::produceParticle(const LorentzMomentum & pp, Energy m) const {
  PPtr p(produceParticle(Lorentz5Momentum(pp, m)));
  return p;
}

PPtr ParticleData::produceParticle(Energy m, const Momentum3 & pp) const {
  PPtr p(produceParticle(Lorentz5Momentum(m, pp)));
  return p;
}

PPtr ParticleData::produceParticle(const Momentum3 & pp) const {
  PPtr p(produceParticle(Lorentz5Momentum(generateMass(), pp)));
  return p;
}

PPtr ParticleData::
produceParticle(Energy plus, Energy minus, Energy px, Energy py) const {
  PPtr p(produceParticle(LorentzMomentum(px, py, 0.5*(plus-minus),
					 0.5*(plus+minus))));
  return p;
}

void ParticleData::setMass(Energy mi) {
  theMass = mi;
  ParticleData * apd = CC().operator->();
  if ( synchronized() && apd ) apd->theMass = theMass;
}

Energy ParticleData::defMass() const {
  return theDefMass;
}

void ParticleData::setWidth(Energy wi) {
  width(wi);
}

Energy ParticleData::getWidth() const {
  return width();
}

Energy ParticleData::defWidth() const {
  return theDefWidth;
}

void ParticleData::setCut(Energy ci) {
  widthCut(ci);
}

Energy ParticleData::getCut() const {
  return (theWidthUpCut >= ZERO && theWidthLoCut >= ZERO)?
    max(theWidthUpCut, theWidthLoCut): min(theWidthUpCut, theWidthLoCut);
}

Energy ParticleData::defCut() const {
  return theDefCut;
}

void ParticleData::setUpCut(Energy ci) {
  widthUpCut(ci);
}

Energy ParticleData::getUpCut() const {
  return theWidthUpCut;
}

void ParticleData::setLoCut(Energy ci) {
  widthLoCut(ci);
}

Energy ParticleData::getLoCut() const {
  return theWidthLoCut;
}

void ParticleData::setCTau(Length ti) {
  cTau(ti);
}

Length ParticleData::getCTau() const {
  return cTau();
}

Length ParticleData::defCTau() const {
  return theDefCTau;
}

void ParticleData::setStable(long is) {
  stable(is);
}

long ParticleData::getStable() const {
  return stable();
}

void ParticleData::setSync(long is) {
  synchronized(is);
}

long ParticleData::getSync() const {
  return synchronized();
}

void ParticleData::setVariableRatio(long is) {
  variableRatio(is);
}

long ParticleData::getVariableRatio() const {
  return variableRatio();
}

string ParticleData::doSync(string) {
  synchronize();
  return "";
}

void ParticleData::setMassGenerator(MassGenPtr gi) {
  massGenerator(gi);
}

void ParticleData::setWidthGenerator(WidthGeneratorPtr wg) {
  widthGenerator(wg);
}

void ParticleData::setColour(long c) {
  theColour = PDT::Colour(c);
}

long ParticleData::getColour() const {
  return theColour;
}

long ParticleData::defColour() const {
  return theDefColour;
}

void ParticleData::setCharge(int c) {
  theCharge = PDT::Charge(c);
}

string ParticleData::ssetCharge(string arg) {
  istringstream is(arg);
  long i;
  if ( is >> i ) {
    theCharge = PDT::Charge(i);
    return "New charge is " + arg;
  }
  if ( arg == "unknown" )
    theCharge = PDT::ChargeUnknown;
  else if ( arg == "charged" )
    theCharge = PDT::Charged;
  else if ( arg == "positive" )
    theCharge = PDT::Positive;
  else if ( arg == "negative" )
    theCharge = PDT::Negative;
  else throw ParticleChargeCommand(*this, arg);
  return  "New charge is " + arg;
}

int ParticleData::getCharge() const {
  return theCharge;
}

int ParticleData::defCharge() const {
  return theDefCharge;
}

void ParticleData::setSpin(int s) {
  theSpin = PDT::Spin(s);
}

int ParticleData::getSpin() const {
  return theSpin;
}

int ParticleData::defSpin() const {
  return theDefSpin;
}

ClassDescription<ParticleData> ParticleData::initParticleData;

struct ParticleOrdering {
  bool operator()(tcPDPtr p1, tcPDPtr p2) {
    return abs(p1->id()) > abs(p2->id()) ||
      ( abs(p1->id()) == abs(p2->id()) && p1->id() > p2->id() ) ||
      ( p1->id() == p2->id() && p1->fullName() > p2->fullName() );
  }
};

struct ModeOrdering {
  bool operator()(const tcDMPtr & d1, const tcDMPtr & d2) {
    ParticleOrdering ord;
    return ord(d1->parent(), d2->parent()) ||
      ( !ord(d2->parent(), d1->parent()) &&
	( d1->tag() < d2->tag() ||
	  ( d1->tag() == d2->tag() && d1->fullName() < d2->fullName() ) ) );
  }
};

void ParticleData::persistentOutput(PersistentOStream & os) const {
  multiset<tcDMPtr,ModeOrdering>
    modes(theDecayModes.begin(), theDecayModes.end());
  os << long(theId) << thePDGName << ounit(theMass, GeV) << ounit(theWidth, GeV)
     << ounit(theWidthUpCut, GeV) << ounit(theWidthLoCut, GeV)
     << ounit(theCTau, mm) << oenum(theCharge) << oenum(theSpin)
     << oenum(theColour);
  os << theMassGenerator << isStable << modes << theDecaySelector
     << theWidthGenerator << theVariableRatio << theAntiPartner << syncAnti
     << ounit(theDefMass, GeV) << ounit(theDefWidth, GeV)
     << ounit(theDefCut, GeV) << ounit(theDefCTau, mm) << oenum(theDefColour)
     << oenum(theDefCharge) << oenum(theDefSpin);
}

void ParticleData::persistentInput(PersistentIStream & is, int) {
  long id;
  is >> id >> thePDGName >> iunit(theMass, GeV) >> iunit(theWidth, GeV)
     >> iunit(theWidthUpCut, GeV) >> iunit(theWidthLoCut, GeV)
     >> iunit(theCTau, mm) >> ienum(theCharge) >> ienum(theSpin)
     >> ienum(theColour) >> theMassGenerator >> isStable
     >> theDecayModes >> theDecaySelector >> theWidthGenerator >> theVariableRatio
     >> theAntiPartner >> syncAnti >> iunit(theDefMass, GeV)
     >> iunit(theDefWidth, GeV) >> iunit(theDefCut, GeV)
     >> iunit(theDefCTau, mm) >> ienum(theDefColour) >> ienum(theDefCharge)
     >> ienum(theDefSpin);
  theId = id;
}

void ParticleData::Init() {

  static ClassDocumentation<ParticleData> documentation
    ("There is no documentation for the ThePEG::ParticleData class");

  static Parameter<ParticleData,Energy> interfaceMass
    ("NominalMass",
     "The nominal mass in GeV of the particle. The actual mass "
     "of a particle instance is generated depending on the "
     "nominal mass and the width and is generated by the "
     "<interface>Mass_generator</interface> object associated with the "
     "particle.",
     &ParticleData::theMass, GeV, ZERO, ZERO, Constants::MaxEnergy,
     false, false, Interface::lowerlim,
     &ParticleData::setMass, 0, 0, 0, &ParticleData::defMass);

  static Parameter<ParticleData,Energy> interfaceDefMass
    ("DefaultMass",
     "The default nominal mass in GeV of the particle. The actual mass "
     "of a particle instance is generated depending on the "
     "nominal mass and the width and is generated by the "
     "<interface>Mass_generator</interface> object associated with the "
     "particle.",
     &ParticleData::theDefMass, GeV, ZERO, ZERO, Constants::MaxEnergy,
     false, true, Interface::lowerlim);
  interfaceDefMass.setHasDefault(false);

  static Parameter<ParticleData,Energy> interfaceWidth
    ("Width",
     "The width of the particle in GeV.",
     0, GeV, ZERO, ZERO, Constants::MaxEnergy,
     false, false, Interface::lowerlim,
     &ParticleData::setWidth, &ParticleData::getWidth,
     0, 0, &ParticleData::defWidth);

  static Parameter<ParticleData,Energy> interfaceDefWidth
    ("DefaultWidth",
     "The default width of the particle in GeV.",
     &ParticleData::theDefWidth, GeV, ZERO, ZERO, Constants::MaxEnergy,
     false, true, Interface::lowerlim);
  interfaceDefWidth.setHasDefault(false);

  static Parameter<ParticleData,Energy> interfaceWidthUpCut
    ("WidthUpCut",
     "The upper hard cutoff in GeV in generated mass, which is the maximum "
     "allowed upwards deviation from the nominal mass. A negative value "
     "corresponds to no limit.",
     0, GeV, ZERO, -1.0*GeV, Constants::MaxEnergy,
     false, false, Interface::lowerlim,
     &ParticleData::setUpCut, &ParticleData::getUpCut,
     0, 0, &ParticleData::defCut);
  
  static Parameter<ParticleData,Energy> interfaceWidthLoCut
    ("WidthLoCut",
     "The lower hard cutoff in GeV in generated mass, which is the maximum "
     "allowed downwards deviation from the nominal mass. A negative value "
     "corresponds to no limit.",
     0, GeV, ZERO, -1.0*GeV, Constants::MaxEnergy,
     false, false, Interface::lowerlim,
     &ParticleData::setLoCut, &ParticleData::getLoCut,
     0, 0, &ParticleData::defCut);
  
  static Parameter<ParticleData,Energy> interfaceWidthCut
    ("WidthCut",
     "The hard cutoff in GeV in generated mass, which is the maximum "
     "allowed deviation from the nominal mass. Sets both the upper and lower "
     "cut. (The displayed value is the maximium of the upper and lower cut.) "
     "A negative value corresponds to no limit.",
     0, GeV, ZERO, -1.0*GeV, Constants::MaxEnergy,
     false, false, Interface::lowerlim,
     &ParticleData::setCut, &ParticleData::getCut,
     0, 0, &ParticleData::defCut);
  interfaceWidthCut.setHasDefault(false);
  
  static Parameter<ParticleData,Energy> interfaceDefWidthCut
    ("DefaultWidthCut",
     "The default hard cutoff in GeV in generated mass, which is the maximum "
     "allowed deviation from the nominal mass. For the actual cutoff, the "
     "upper and lower cut can be set separately.",
     &ParticleData::theDefCut, GeV, ZERO, ZERO, Constants::MaxEnergy,
     false, true, Interface::lowerlim);
  interfaceDefWidthCut.setHasDefault(false);
  
  static Parameter<ParticleData,Length> interfaceCTau
    ("LifeTime",
     "c times the average lifetime of the particle measuerd in mm."
     "The actual lifetime of a particle instance is generated "
     "from this number by the <interface>Mass_generator</interface> "
     "object associated with the particle.",
     0, mm, ZERO, ZERO, Constants::MaxLength,
     false, false, Interface::lowerlim,
     &ParticleData::setCTau, &ParticleData::getCTau,
     0, 0, &ParticleData::defCTau);
  interfaceCTau.setHasDefault(false);

  static Parameter<ParticleData,Length> interfaceDefCTau
    ("DefaultLifeTime",
     "c times the default average lifetime of the particle measuerd in mm."
     "The actual lifetime of a particle instance is generated "
     "from this number by the <interface>Mass_generator</interface> "
     "object associated with the particle.",
     &ParticleData::theDefCTau, mm, ZERO, ZERO, Constants::MaxLength,
     false, true, Interface::lowerlim);
  interfaceDefCTau.setHasDefault(false);

  static Switch<ParticleData> interfaceColour
    ("Colour",
     "The colour quantum number of this particle type.",
     0, -1, false, false, &ParticleData::setColour, &ParticleData::getColour,
     &ParticleData::defColour);
  static SwitchOption interfaceColourUndefined
    (interfaceColour, "Undefined", "The coulur is undefined.", -1);
  static SwitchOption interfaceColourNeutral
    (interfaceColour, "Neutral", "This particle is colour neutral.", 0);
  static SwitchOption interfaceColour3
    (interfaceColour, "Triplet", "This particle is a colour triplet.", 3);
  static SwitchOption interfaceColour3bar
    (interfaceColour, "AntiTriplet",
     "This particle is a colour anti-triplet.", -3);
  static SwitchOption interfaceColour6
    (interfaceColour, "Sextet", "This particle is a colour sextet.", 6);
  static SwitchOption interfaceColour6bar
    (interfaceColour, "AntiSextet",
     "This particle is a colour anti-sextet.", -6);
  static SwitchOption interfaceColour8
    (interfaceColour, "Octet", "This particle is a colour octet.", 8);
  
  static Switch<ParticleData,PDT::Colour> interfaceDefColour
    ("DefaultColour",
     "The default colour quantum number of this particle type.",
     &ParticleData::theDefColour, PDT::Colour(-1), false, true);
  static SwitchOption interfaceDefColourUndefined
    (interfaceDefColour, "Undefined", "The coulur is undefined.", -1);
  static SwitchOption interfaceDefColourNeutral
    (interfaceDefColour, "Neutral", "This particle is colour neutral.", 0);
  static SwitchOption interfaceDefColour3
    (interfaceDefColour, "Triplet", "This particle is a colour triplet.", 3);
  static SwitchOption interfaceDefColour3bar
    (interfaceDefColour, "AntiTriplet",
     "This particle is a colour anti-triplet.", -3);
  static SwitchOption interfaceDefColour6
    (interfaceDefColour, "Sextet", "This particle is a colour sextet.", 6);
  static SwitchOption interfaceDefColour6bar
    (interfaceDefColour, "AntiSextet",
     "This particle is a colour anti-sextet.", -6);
  static SwitchOption interfaceDefColour8
    (interfaceDefColour, "Octet", "This particle is a colour octet.", 8);
  interfaceDefColour.setHasDefault(false);  

  static Parameter<ParticleData, int> interfaceCharge
    ("Charge",
     "The charge of this particle in units of e/3. "
     "See also the command interface <interface>SetCharge</interface>.",
     0, 0, -24, 24, false, false, true,
     &ParticleData::setCharge, &ParticleData::getCharge, 0, 0,
     &ParticleData::defCharge);

  static Parameter<ParticleData, PDT::Charge> interfaceDefCharge
    ("DefaultCharge",
     "The default charge of this particle in units of e/3. "
     "See also the command interface <interface>SetCharge</interface>.",
     &ParticleData::theDefCharge, PDT::Charge(0), PDT::Charge(-24),
     PDT::Charge(24), false, true, true);
  interfaceDefCharge.setHasDefault(false);

  static Command<ParticleData> interfaceSetCharge
    ("SetCharge",
     "Set the charge of this particle. The argument should be given as an "
     "interger giving three times the unit charge, or 'unknown', "
     "'charged', 'positive' or 'negative'", &ParticleData::ssetCharge);

  static Parameter<ParticleData, int> interfaceSpin
    ("Spin",
     "The spin quantim number of this particle on the form 2j+1.",
     0, 0, 0, 9, false, false, true,
     &ParticleData::setSpin, &ParticleData::getSpin, 0, 0,
     &ParticleData::defSpin);

  static Parameter<ParticleData, PDT::Spin> interfaceDefSpin
    ("DefaultSpin",
     "The default spin quantim number of this particle on the form 2j+1.",
     &ParticleData::theDefSpin, PDT::Spin(0), PDT::Spin(0), PDT::Spin(9),
     false, true, true);
  interfaceDefSpin.setHasDefault(false);

  static Switch<ParticleData> interfaceStable
    ("Stable",
     "Indicates if the particle is stable or not.",
     0, 0, false, false,
     &ParticleData::setStable, &ParticleData::getStable, 0);
  static SwitchOption interfaceStableYes
    (interfaceStable,
     "Stable",
     "This particle is stable",
     1);
  static SwitchOption interfaceStableNo
    (interfaceStable,
     "Unstable",
     "This particle is not stable",
     0);
  interfaceStable.setHasDefault(false);

  static Switch<ParticleData> interfaceVariableRatio
    ("VariableRatio",
     "Indicates if the branching ratios of the particle are allowed"
     " to vary for given Particle instances depending on the mass of the instance.",
     0, 0, false, false,
     &ParticleData::setVariableRatio, &ParticleData::getVariableRatio, 0);
  static SwitchOption interfaceVariableRatioYes
    (interfaceVariableRatio,
     "Yes",
     "The branching ratio varies.",
     1);
  static SwitchOption interfaceVariableRatioNo
    (interfaceVariableRatio,
     "No",
     "The branching ratio does not vary.",
     0);

  static Switch<ParticleData> interfaceSync
    ("Synchronized",
     "Indicates if the changes to this particle is propagated to "
     "its anti-partner or not. Note that setting this switch does not "
     "actually synchronize the properties with the anti-partner, "
     "it only assures that following changes are propagated. "
     "To sync the particle with its anti-particle, use the "
     "<interface>Synchronize</interface> command.",
     0, 1, false, false,
     &ParticleData::setSync, &ParticleData::getSync, 0);
  static SwitchOption interfaceSyncYes
    (interfaceSync,
     "Synchronized",
     "Changes to this particle will propagate to its "
     "anti-partner",
     1);
  static SwitchOption interfaceSyncNo
    (interfaceSync,
     "Not_synchronized",
     "Changes to this particle will propagate to its "
     "anti-partner",
     0);
  interfaceSync.setHasDefault(false);

  static Command<ParticleData> interfaceSynchronize
    ("Synchronize",
     "Synchronizes this particle so that all its properties "
     "correspond to those of its anti-partner",
     &ParticleData::doSync, false);

  static Reference<ParticleData,MassGenerator> interfaceMassGenerator
    ("Mass_generator",
     "An object derived from the ThePEG::MassGenerator"
     "class, which is able to generate a mass for a given "
     "particle instance",
     &ParticleData::theMassGenerator, false, false, true, true,
     &ParticleData::setMassGenerator, 0, 0);

  static Reference<ParticleData,WidthGenerator> interfaceWidthGenerator
    ("Width_generator",
     "An object derived from the ThePEG::WidthGenerator class, "
     "which is able to calculate the full and partial widths for"
     "this particle type and for a given instance of this "
     "particle type.",
     &ParticleData::theWidthGenerator, false, false, true, true,
     &ParticleData::setWidthGenerator, 0, 0);

  static RefVector<ParticleData,DecayMode> interfaceDecayModes
    ("DecayModes",
     "The list of decay modes defined for this particle type.",
     0, -1, false, false, false, false, 0,
     &ParticleData::insDecayModes, &ParticleData::delDecayModes,
     &ParticleData::getDecayModes);


  static Command<ParticleData> interfaceSelectDecayModes
    ("SelectDecayModes",
     "Only the decay modes which are given as (white-space separated) "
     "decay tags will be switched on, all others will be switched off. "
     "If no argument or 'none' is given, all decay modes are switched off. "
     "If the argument is 'all', all decay modes are switched on.",
     &ParticleData::doSelectDecayModes, false);


  static Command<ParticleData> interfacePrintDecayModes
    ("PrintDecayModes",
     "Print all decay modes of this particle.",
     &ParticleData::doPrintDecayModes, true);


  interfaceStable.rank(14);
  interfaceDecayModes.rank(13);
  interfaceMass.rank(12);
  interfaceWidth.rank(11);
  interfaceWidthCut.rank(10);
  interfaceCTau.rank(9);
  interfaceMassGenerator.rank(8);
  interfaceWidthGenerator.rank(7);
  interfaceWidthUpCut.rank(-0.1);
  interfaceWidthLoCut.rank(-0.1);


}

string ParticleData::doPrintDecayModes(string) {
  multimap<double,tDMPtr, std::greater<double> > sorted;
  for ( DecaySet::iterator it = decayModes().begin();
	it != decayModes().end(); ++it )
    sorted.insert(make_pair((**it).brat(), *it));
  ostringstream os;
  for ( multimap<double,tDMPtr, 
	  std::greater<double> >::iterator it = sorted.begin();
	it != sorted.end(); ++it )
    os << it->second->tag()
       << (it->second->on()? " ": " (off) ")
       << it->first << endl;
  return os.str();
}


string ParticleData::doSelectDecayModes(string args) {
  DecaySet on;
  while ( !args.empty() ) {
    string arg = StringUtils::car(args);
    if ( arg == "all" ) {
      on = decayModes();
      break;
    }
    if ( arg == "none" ) {
      on.clear();
      break;
    }
    string name = arg;
    args = StringUtils::cdr(args);
    if ( arg.empty() ) continue;
    if ( arg[0] != '/' ) arg = fullName() + "/" + arg;
    DMPtr dm = Repository::GetPtr<DMPtr>(arg);
    if ( !dm ) return "Error: No decay mode with tag '" + name + "' exists.";
    on.insert(dm);
  }
  for ( DecaySet::iterator it = decayModes().begin();
	it != decayModes().end(); ++it ) {
    if ( on.find(*it) != on.end() ) {
      (**it).switchOn();
      on.erase(*it);
    } else {
      (**it).switchOff();
    }
  }
  if ( !on.empty() )
    return "Error: decay mode '" + (**on.begin()).tag() + "'was not available.";
  return "";
}

void ParticleData::insDecayModes(DMPtr dm, int) {
  addDecayMode(dm);
}

void ParticleData::delDecayModes(int i) {
  vector<DMPtr> mv = getDecayModes();
  if ( i >= 0 && static_cast<unsigned int>(i) < mv.size() ) removeDecayMode(mv[i]);
}

vector<DMPtr> ParticleData::getDecayModes() const {
  return vector<DMPtr>(theDecayModes.begin(), theDecayModes.end());
}

ParticleChargeCommand::
ParticleChargeCommand(const ParticleData & pd, string arg) {
  theMessage << "Cannot set the charge of particle '" << pd.name()
	     << "' to '" << arg << "'.";
  severity(warning);
}

void ParticleData::doinit() {
  Interfaced::doinit();
  if( theMassGenerator )  theMassGenerator->init();
  if( theWidthGenerator ) theWidthGenerator->init();
}

void ParticleData::doinitrun() {
  Interfaced::doinitrun();
  if( theMassGenerator )  theMassGenerator->initrun();
  if( theWidthGenerator ) theWidthGenerator->initrun();
}

}
