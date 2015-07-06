// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SpecialGluonEmitter class.
//

#include "SpecialGluonEmission.h"
#include "SpecialGluonEmitter.h"
#include "FSGluonEmitter.h"
#include "RemnantModel.h"
#include "ColourResonanceModel.h"
#include "Ariadne/Cascade/QCDDipole.h"
#include "Ariadne/Cascade/DipoleState.h"
#include "Ariadne/Cascade/AriadneHandler.h"
#include "Ariadne/Cascade/Junction.h"
#include "Ariadne/Cascade/RemnantParton.h"
#include "Ariadne/Cascade/ResonanceParton.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Utilities/UtilityBase.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "Ariadne/Config/UnitFO.h"


#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Ariadne5;

SpecialGluonEmitter::SpecialGluonEmitter() {}

SpecialGluonEmitter::~SpecialGluonEmitter() {}

IBPtr SpecialGluonEmitter::clone() const {
  return new_ptr(*this);
}

IBPtr SpecialGluonEmitter::fullclone() const {
  return new_ptr(*this);
}


bool SpecialGluonEmitter::canHandle(const DipoleBase & e) const {
  tcQCDPtr d = dynamic_ptr_cast<tcQCDPtr>(&e);
  if ( !d ) return false;
  if ( d->iPart()->special(d->oPart()) || d->oPart()->special(d->iPart()) )
    return true;
  return false;
}

bool SpecialGluonEmitter::touched(const DipoleBase & dipole) const {
  if ( dipole.touched() ) return true;
  tcQCDPtr d = dynamic_ptr_cast<tcQCDPtr>(&dipole);
  return d->iPart()->touched() || d->oPart()->touched() ||
    d->iPart()->touchedNeighbours(d) || d->oPart()->touchedNeighbours(d);
}  

bool SpecialGluonEmitter::
overrides(const EmitterBase & em, DipoleBase &) const {
  if ( typeid(em) == typeid(FSGluonEmitter) ) return true;
  return false;
}

EmPtr SpecialGluonEmitter::
generate(const DipoleBase & dipole, Energy rhomin, Energy rhomax) const {
  const QCDDipole & d = dynamic_cast<const QCDDipole &>(dipole);

  tParPtr ip = d.iPart();
  tParPtr op = d.oPart();
  bool ijunc = false;
  bool ojunc = false;
  if ( tJunctionPtr j = dynamic_ptr_cast<tJunctionPtr>(ip) ) {
    ijunc = true;
    ip = j->getRandomRecoiler(&d);
  }
  if ( tJunctionPtr j = dynamic_ptr_cast<tJunctionPtr>(op) ) {
    if ( ijunc ) return EmPtr(); // ignore dipoles between two junctions.
    op = j->getRandomRecoiler(&d);
    ojunc = true;
  }

  // Set pseudo-partons
  PseudoParton pip(ip->isG(), !ip->isG(), PseudoParton::normal,
		   ip->momentum(), ip);
  PseudoParton pop(op->isG(), !op->isG(), PseudoParton::normal,
		   op->momentum(), op);

  // Set soft remnant pointers
  // *** ATTENTION *** What about recoils for hard remnants?
  tRemParPtr remip = dynamic_ptr_cast<tRemParPtr>(ip);
  tRemParPtr remop = dynamic_ptr_cast<tRemParPtr>(op);
  if ( remip && remip->hard() ) remip = tRemParPtr();
  if ( remop && remop->hard() ) remop = tRemParPtr();

  // Set pointers to possible resonance products.
  tResParPtr resip = dynamic_ptr_cast<tResParPtr>(ip);
  tResParPtr resop = dynamic_ptr_cast<tResParPtr>(op);

  const RemnantModel & remmod = Current<AriadneHandler>()->remnantModel();

  tColourResonanceModelPtr resmod =
    Current<AriadneHandler>()->colourResonanceModel();
  if ( ( resip || resop ) && !resmod )
    Throw<MissingModel>()
      << "SpecialGluonEmitter '" << name() << "' could not find a "
      << "ColourResonanceModel object assigned to the AriadneHandler class."
      << Exception::runerror;

  // Now setup pseudo partons for the main gluon emission.
  if ( remip )
    pip = remmod.getPseudoParton(remip);
  else if ( resip && resmod->stillSpecial(resip) )
    pip = resmod->getPseudoParton(resip);

  if ( remop )
    pop = remmod.getPseudoParton(remop);
  else if ( resop && resmod->stillSpecial(resop) )
    pop = resmod->getPseudoParton(resop);

  if ( remip && remop ) return EmPtr();
  // Now setup the main gluon emission.
  if ( remip ) pop.noRecoil = false;
  if ( remop ) pip.noRecoil = false;

  double C = (pip.isGluon || pop.isGluon? 3.0/4.0: 2.0/3.0)/
    Constants::pi;
  int n1 =  pip.isGluon? 3: 2;
  int n3 =  pop.isGluon? 3: 2;
  Energy2 S = (pip.p + pop.p).m2();
  if ( S <= sqr(rhomin*2.0) ) return EmPtr();
  Energy W = sqrt(S);
  double y1 = max(pip.p.mass2()/S, 0.0);
  double y3 = max(pop.p.mass2()/S, 0.0);
  SpecialGluonEmission e(*this, dipole, y1, y3, S, pip, pop);

  C *= preweight(e);

  rhomax = min(rhomax, W/2.0);
  if ( rhomax <= rhomin ) return EmPtr();
  double yint = 2.0*acosh(0.5*W/rhomin);

  while (true) {
    if ( rhomax <= rhomin ) return EmPtr();
    double weight =
      FSGluonEmitter::gluonEngine(e, rhomin, rhomax, W, C, yint, n1, n3);
    if ( weight < 0.0 ) return EmPtr();
    rhomax = e.rho;
    if ( weight == 0.0 ) continue;
    weight *= reweight(e);

    // Only do specieal stuff if we have remnants.

    try {

      if ( remip || remop )
	e.genmom = getMomenta(S, e.x1, e.x3, pip.p.mass(), ZERO, pop.p.mass(),
			      pip.noRecoil, pop.noRecoil,
			      UseRandom::rnd(2.0*Constants::pi), pip.p, pop.p);
      if ( remip )
	weight *=
	  (e.wrem1 = remmod.reweightFS(remip, e.rho,
				       e.genmom.second, e.genmom.first));

      if ( remop )
	weight *=
	  (e.wrem3 = remmod.reweightFS(remop, e.rho,
				       e.genmom.second, e.genmom.third));

      e.colourParent = ip;
      e.mainParent = e.antiColourParent = op;
      // The main parent is the one which looses most energy.
      if ( sqr(e.x3) > UseRandom::rnd()*(sqr(e.x1) + sqr(e.x3)) )
	e.mainParent = ip;
      // Do not allow a junction side to be the the main parent.
      if ( ijunc && e.mainParent == ip ) continue;
      if ( ojunc && e.mainParent == op ) continue;

    } catch ( ImpossibleKinematics ) {
      continue;
    }

    if ( weight > UseRandom::rnd() ) break;

  }

  // Save information in the Emission object.
  e.ymax = log(W/Current<AriadneHandler>()->pTCut());
  e.radiators.push_back(op);
  e.radiators.push_back(ip);
  e.pold = make_pair(pip.p, pop.p);
  if ( e.mainParent == ip ) swap(e.radiators[0], e.radiators[1]);

  return new_ptr(e);

}


bool SpecialGluonEmitter::
perform(const Emission & emission) const {
  tColourResonanceModelPtr resmod =
    Current<AriadneHandler>()->colourResonanceModel();
  QCDDipole & d = dynamic_cast<QCDDipole &>(*emission.dipole);
  const SpecialGluonEmission & e =
    dynamic_cast<const SpecialGluonEmission &>(emission);
  tParPtr ip = e.colourParent;
  tParPtr op = e.antiColourParent;
  tParPtr g = FSGluonEmitter::insertGluon(d, e.mainParent).first;

  // Set soft remnant pointers;
  tRemParPtr remip = dynamic_ptr_cast<tRemParPtr>(ip);
  tRemParPtr remop = dynamic_ptr_cast<tRemParPtr>(op);
  if ( remip && remip->hard() ) remip = tRemParPtr();
  if ( remop && remop->hard() ) remop = tRemParPtr();

  tResParPtr resip = dynamic_ptr_cast<tResParPtr>(ip);
  tResParPtr resop = dynamic_ptr_cast<tResParPtr>(op);
  PseudoParton pip = e.pip;
  PseudoParton pop = e.pop;
  // If none of the pseudo-partons were remnants we cannot be sure the
  // momenta of the pseudo-partons hasn't changed. Also the new
  // momenta were not pre-generated.
  if ( !remip && !remop ) {
    pip = PseudoParton(ip->isG(), !ip->isG(), PseudoParton::normal,
		       ip->momentum(), ip);
    pop = PseudoParton(op->isG(), !op->isG(), PseudoParton::normal,
		       op->momentum(), op);
    if ( resip ) pip = resmod->getPseudoParton(resip);
    if ( resop ) pop = resmod->getPseudoParton(resop);
    e.pold = make_pair(pip.p, pop.p);
    e.genmom = getMomenta(e.S, e.x1, e.x3, pip.p.mass(), ZERO, pop.p.mass(),
			  pip.noRecoil, pop.noRecoil,
			  UseRandom::rnd(2.0*Constants::pi), pip.p, pop.p);
  }

  g->momentum() = e.genmom.second;
  g->setVertex(e.mainParent->vertex());

  g->emission(&e);
  d.state()->addHadronicFS(g);

  e.partons.push_back(g);

  if ( remip ) remip->setMomentum(e.genmom.first);
  else if ( resip ) resip->setMomentum(e.genmom.first);
  else ip->momentum() = e.genmom.first;

  if ( remop ) remop->setMomentum(e.genmom.third);
  else if ( resop ) resop->setMomentum(e.genmom.third);
  else op->momentum() = e.genmom.third;

  return true;
}

void SpecialGluonEmitter::revert(const Emission & emission) const {
  tColourResonanceModelPtr resmod =
    Current<AriadneHandler>()->colourResonanceModel();
  QCDDipole & d = dynamic_cast<QCDDipole &>(*emission.dipole);
  const SpecialGluonEmission & e =
    dynamic_cast<const SpecialGluonEmission &>(emission);

  FSGluonEmitter::removeGluon(d, d.oPart(), e.mainParent);

  tParPtr ip = e.colourParent;
  tParPtr op = e.antiColourParent;

  e.pip.realParton->setMomentum(e.pold.first);
  e.pop.realParton->setMomentum(e.pold.second);
  d.state()->untouchHadronicState();

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void SpecialGluonEmitter::persistentOutput(PersistentOStream &) const {}

void SpecialGluonEmitter::persistentInput(PersistentIStream &, int) {}

DescribeClass<SpecialGluonEmitter,FSGluonEmitter>
describeAriadne5SpecialGluonEmitter("Ariadne5::SpecialGluonEmitter",
				    "libAriadne5.so");

void SpecialGluonEmitter::Init() {

  static ClassDocumentation<SpecialGluonEmitter> documentation
    ("The SpecialGluonEmitter class implements the standard classical "
     "gluon emission from a colour dipole between special partons "
     "(junctions, remnants or coloured resonance products).");

}

