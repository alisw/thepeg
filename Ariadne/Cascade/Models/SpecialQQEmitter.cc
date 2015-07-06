// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SpecialQQEmitter class.
//

#include "SpecialQQEmission.h"
#include "SpecialQQEmitter.h"
#include "FSQQEmitter.h"
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


#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Ariadne5;

SpecialQQEmitter::SpecialQQEmitter() {}

SpecialQQEmitter::~SpecialQQEmitter() {}

IBPtr SpecialQQEmitter::clone() const {
  return new_ptr(*this);
}

IBPtr SpecialQQEmitter::fullclone() const {
  return new_ptr(*this);
}


bool SpecialQQEmitter::canHandle(const DipoleBase & e) const {
  tcQCDPtr d = dynamic_ptr_cast<tcQCDPtr>(&e);
  if ( !d ) return false;
  bool ispec = d->iPart()->special(d->oPart());
  bool ospec = d->oPart()->special(d->iPart());
  if ( d->iPart()->isG() && !ispec && ospec ) return true;
  if ( d->oPart()->isG() && !ospec && ispec ) return true;
  return false;
}

bool SpecialQQEmitter::touched(const DipoleBase & dipole) const {
  if ( dipole.touched() ) return true;
  tcQCDPtr d = dynamic_cast<const QCDDipole *>(&dipole);
  return d->iPart()->touched() || d->oPart()->touched() ||
    d->iPart()->touchedNeighbours(d) || d->oPart()->touchedNeighbours(d);
}  

bool SpecialQQEmitter::
overrides(const EmitterBase & em, DipoleBase &) const {
  if ( typeid(em) == typeid(FSQQEmitter) ) return true;
  return false;
}

EmPtr SpecialQQEmitter::
generate(const DipoleBase & dipole, Energy rhomin, Energy rhomax) const {
  const QCDDipole & d = dynamic_cast<const QCDDipole &>(dipole);

  EmPtr esel;
  tParPtr ip = d.iPart();
  tParPtr op = d.oPart();
  bool ijunc = false;
  //  bool ojunc = false;
  if ( tJunctionPtr j = dynamic_ptr_cast<tJunctionPtr>(ip) ) {
    ijunc = true;
    ip = j->getRandomRecoiler(&d);
  }
  if ( tJunctionPtr j = dynamic_ptr_cast<tJunctionPtr>(op) ) {
    if ( ijunc ) return EmPtr(); // ignore dipoles between two junctions.
    op = j->getRandomRecoiler(&d);
    //    ojunc = true;
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
      << "SpecialQQEmitter '" << name() << "' could not find a "
      << "ColourResonanceModel object assigned to the AriadneHandler class."
      << Exception::runerror;

  // Now generate special emissions for the special partons (unless
  // they were originally junctions) and setup pseudo partons for the
  // main gluon emission.
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

  Energy2 S = (pip.p + pop.p).m2();
  if ( S <= sqr(rhomin*2.0) ) return EmPtr();
  double C1 = 0.0;
  double C3 = 0.0;
  if ( d.iPart()->isG() && !d.iPart()->special(d.oPart()) )
    C1 = 1.0/(4.0*Constants::pi*(1.0 + S/d.prev()->sdip()));
  if ( d.oPart()->isG() && !d.oPart()->special(d.iPart()) )
    C3 = 1.0/(4.0*Constants::pi*(1.0 + S/d.next()->sdip()));
  if ( C1 + C3 <= 0.0 ) return EmPtr();

  tParPtr gluon = C1 > 0.0? d.iPart(): d.oPart();
  tRemParPtr orem = C1 > 0.0? remop: remip;
  PseudoParton pother = C1 > 0.0? pop: pip;
  tParPtr other = C1 > 0.0 ? op: ip;

  // We will pretend parton 3 is the gluon even if it isn't. Hence
  // parton 1 may have a mass or not.
  double y1 = 0.0;
  if ( C1 == 0.0 ) y1 = max(pip.p.mass2()/S, 0.0);
  if ( C3 == 0.0 ) y1 = max(pop.p.mass2()/S, 0.0);

  double sy1 = sqrt(y1);
  Energy W = sqrt(S);

  rhomax = min(rhomax, W/2.0);
  if ( rhomax <= rhomin ) return EmPtr();

  const int nfl = Current<AriadneHandler>()->nFlav();
  SpecialQQEmission e(*this, dipole, 0, ZERO, S, pip, pop);

  for ( int ifl = 1; ifl <= nfl; ++ifl ) {

    tcPDPtr q = generator()->getParticleData(ifl);
    if ( !q ) continue;

    e.mq = q->generateMass();
    e.ifl = ifl;

    double syq = e.mq/W;

    Energy rhom = rhomax;

    if ( rhom <= rhomin ) continue;

    double yint = 2.0*acosh(0.5*W/rhomin);

    double CW = preweight(e);
    C1 *= CW;
    C3 *= CW;

    while (true) {
      if ( rhom <= rhomin ) break;
      double weight =
	FSQQEmitter::qqbarEngine(e, rhomin, rhom, W, C1 + C3, yint, sy1, syq);
      if ( weight < 0.0 ) break;
      rhom = e.rho;
      if ( weight == 0.0 ) continue;
      weight *= reweight(e);

      try {

	if ( orem ) {
	  e.genmom = getMomenta(S, e.x1, e.x3, pother.p.mass(), e.mq, e.mq,
				true, false, UseRandom::rnd(2.0*Constants::pi),
				pother.p, gluon->momentum());
	  weight *=
	    (e.wrem1 = remmod.reweightFS(orem, e.rho,
					 e.genmom.second, e.genmom.first));
	}
      } catch ( ImpossibleKinematics ) {
	continue;
      }

      if ( weight < UseRandom::rnd() ) continue;
  
      // Save information in the Emission object.
      e.ymax = log(W/Current<AriadneHandler>()->pTCut());
      e.ifl = rndbool(C1, C3)? -ifl: ifl;
      e.radiators.push_back(ip);
      e.radiators.push_back(op);
      e.colourParent = ip;
      e.mainParent = e.antiColourParent = op;
      e.pold = make_pair(pop.p, pip.p);
      e.od = d.next();
      if ( e.ifl < 0 ) {
	swap(e.radiators[0], e.radiators[1]);
	swap(e.pold.first, e.pold.second);
	e.mainParent = d.iPart();
	e.y = -e.y;
	e.yo = -e.yo;
	e.od = d.prev();
      }

      rhomin = rhom;
      esel = new_ptr(e);
      break;

    }
  }
  return esel;

}

bool SpecialQQEmitter::
perform(const Emission & emission) const {
  tColourResonanceModelPtr resmod =
    Current<AriadneHandler>()->colourResonanceModel();
  QCDDipole & d = dynamic_cast<QCDDipole &>(*emission.dipole);
  const SpecialQQEmission & e =
    dynamic_cast<const SpecialQQEmission &>(emission);
  tParPtr ip = d.iPart();
  tParPtr op = d.oPart();

  tParPtr gluon = e.ifl > 0? d.oPart(): d.iPart();
  PseudoParton pother = e.ifl > 0? e.pip: e.pop;
  tParPtr other = e.ifl > 0 ? ip: op;

  e.pold = make_pair(gluon->momentum(), pother.p);

  tRemParPtr orem = dynamic_ptr_cast<tRemParPtr>(other);
  if ( orem && orem->hard() ) orem = tRemParPtr();

  tResParPtr ores = dynamic_ptr_cast<tResParPtr>(other);
  // If none of the pseudo-partons were remnants we cannot be sure the
  // momenta of the pseudo-partons hasn't changed. Also the new
  // momenta were not pre-generated.
  if ( !orem ) {
    pother = PseudoParton(other->isG(), !other->isG(), PseudoParton::normal,
		       other->momentum(), other);
    if ( ores ) pother = resmod->getPseudoParton(ores);
    e.pold = make_pair(gluon->momentum(), pother.p);
    try {
      e.genmom = getMomenta(d.sdip(), e.x1, e.x3,
			    other->momentum().mass(), e.mq, e.mq,
			    true, false, rnd(2.0*Constants::pi),
			    pother.p, gluon->momentum());
    }
    catch ( ImpossibleKinematics ) {
      return false;
    }
  }

  pair<tParPtr,tParPtr> qqbar = FSQQEmitter::splitGluon(d, gluon, e.ifl);

  ParPtr q = qqbar.first;
  ParPtr qbar = qqbar.second;

  if ( orem ) orem->setMomentum(e.genmom.first);
  else if ( ores ) ores->setMomentum(e.genmom.first);
  else other->momentum() = e.genmom.first;

  if ( e.ifl < 0 ) swap(e.genmom.second, e.genmom.third);

  qbar->momentum() = e.genmom.second;
  q->momentum() = e.genmom.third;
  e.partons.push_back(q);
  e.partons.push_back(qbar);

  SpecialQQEmission eo(e);
  swap(eo.y, eo.yo);
  EmPtr peo = new_ptr(eo);
  if ( e.ifl > 0.0 ) {
    qbar->emission(&e);
    q->emission(peo);
  } else {
    q->emission(&e);
    qbar->emission(peo);
  }
  return true;

}

void SpecialQQEmitter::revert(const Emission & emission) const {
  tColourResonanceModelPtr resmod =
    Current<AriadneHandler>()->colourResonanceModel();
  QCDDipole & d = dynamic_cast<QCDDipole &>(*emission.dipole);
  const SpecialQQEmission & e =
    dynamic_cast<const SpecialQQEmission &>(emission);
  tParPtr gluon =
    FSQQEmitter::fuseQQBar(d, e.partons[0], e.partons[1], e.od, e.mainParent);

  gluon->setMomentum(e.pold.first);

  tParPtr other = e.ifl > 0? e.colourParent: e.antiColourParent;
  PseudoParton pother = e.ifl > 0? e.pip: e.pop;
  pother.realParton->setMomentum(e.pold.second);
  d.state()->untouchHadronicState();

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void SpecialQQEmitter::persistentOutput(PersistentOStream &) const {}

void SpecialQQEmitter::persistentInput(PersistentIStream &, int) {}

DescribeClass<SpecialQQEmitter,FSQQEmitter>
describeAriadne5SpecialQQEmitter("Ariadne5::SpecialQQEmitter",
				    "libAriadne5.so");

void SpecialQQEmitter::Init() {

  static ClassDocumentation<SpecialQQEmitter> documentation
    ("The SpecialQQEmitter class implements the standard classical "
     "gluon emission from a colour dipole between special partons "
     "(junctions, remnants or coloured resonance products).");

}

