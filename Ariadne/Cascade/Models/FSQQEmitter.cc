// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FSQQEmitter class.
//

#include "FSQQEmitter.h"
#include "FSQQEmission.h"
#include "Ariadne/Cascade/QCDDipole.h"
#include "Ariadne/Cascade/DipoleState.h"
#include "Ariadne/Cascade/RemnantParton.h"
#include "Ariadne/Cascade/AriadneHandler.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Utilities/UtilityBase.h"
#include "ThePEG/Utilities/Current.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Repository/CurrentGenerator.h"


#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Ariadne5;

FSQQEmitter::FSQQEmitter() {}

FSQQEmitter::~FSQQEmitter() {}

IBPtr FSQQEmitter::clone() const {
  return new_ptr(*this);
}

IBPtr FSQQEmitter::fullclone() const {
  return new_ptr(*this);
}


bool FSQQEmitter::canHandle(const DipoleBase & dipole) const {
  tcQCDPtr d = dynamic_cast<const QCDDipole *>(&dipole);
  if ( !d ) return false;
  if ( d->iPart()->special(d->oPart()) || d->oPart()->special(d->iPart()) )
    return false;
  return d->iPart()->isG() || d->oPart()->isG();
}

bool FSQQEmitter::overrides(const EmitterBase &, DipoleBase &) const {
  return false;
}

bool FSQQEmitter::touched(const DipoleBase & dipole) const {
  if ( EmitterBase::touched(dipole) ) return true;
  if ( tcQCDPtr d = dynamic_ptr_cast<tcQCDPtr>(&dipole) ) {
    if ( d->next() && d->next()->touched() ) return true;
    if ( d->prev() && d->prev()->touched() ) return true;
  }
  return false;
}

EmPtr FSQQEmitter::
generate(const DipoleBase & dipole, Energy rhomin, Energy rhomax) const {
  // *** ATTENTION *** investigate why these emissions are rejected by
  // gluon pt veto so often.

  FSQQEmPtr esel;

  const QCDDipole & d = dynamic_cast<const QCDDipole &>(dipole);
  tcParPtr ip = d.iPart();
  tcParPtr op = d.oPart();

  Energy2 S = d.sdip();
  if ( S <= sqr(2.0*rhomin) ) return EmPtr();
  double C1 = 0.0;
  double C3 = 0.0;
  if ( ip->isG() ) C1 = 1.0/(4.0*Constants::pi*(1.0 + S/d.prev()->sdip()));
  if ( op->isG() ) C3 = 1.0/(4.0*Constants::pi*(1.0 + S/d.next()->sdip()));
  if ( C1 + C3 <= 0.0 ) return EmPtr();

  // We will pretend parton 3 is the gluon even if it isn't. Hence
  // parton 1 may have a mass or not.
  double y1 = 0.0;
  if ( !ip->isG() ) y1 = max(ip->momentum().mass2()/S, 0.0);
  if ( !op->isG() ) y1 = max(op->momentum().mass2()/S, 0.0);

  double sy1 = sqrt(y1);
  Energy W = sqrt(S);

  rhomax = min(rhomax, W/2.0);
  if ( rhomax <= rhomin ) return EmPtr();

  const int nfl = Current<AriadneHandler>()->nFlav();
  FSQQEmission e(*this, dipole, 0, ZERO);

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
      double weight = qqbarEngine(e, rhomin, rhom, W, C1 + C3, yint, sy1, syq);
      if ( weight < 0.0 ) break;
      rhom = e.rho;
      if ( weight == 0.0 ) continue;
      weight *= reweight(e);

      // *** ATTENTION *** Add other weights here.
 
      if ( weight < rnd() ) continue;

      e.ymax = log(W/Current<AriadneHandler>()->pTCut());
      e.ifl = rndbool(C1, C3)? -ifl: ifl;
      e.radiators.push_back(d.iPart());
      e.radiators.push_back(d.oPart());
      e.colourParent = d.iPart();
      e.mainParent = e.antiColourParent = d.oPart();
      e.pold = make_pair(d.oPart()->momentum(), d.iPart()->momentum());
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

double FSQQEmitter::
qqbarEngine(FSQQEmission & e, Energy rhomin, Energy rhomax,
	    Energy W, double C, double yint, double sy1, double syq) {
  rhomax = sqrt(rndsud(C*yint, sqr(rhomax), sqr(rhomin)));
  if ( rhomax <= rhomin ) return -1.0;
  e.rho = rhomax;
  double ymax = acosh(0.5*W/rhomax);
  e.y = 2.0*ymax*UseRandom::rnd() - ymax;
  e.x1 = 1.0 - rhomax*exp( e.y)/W + sqr(sy1) - sqr(2.0*syq);
  e.x3 = 1.0 - rhomax*exp(-e.y)/W + sqr(syq) - sqr(sy1 + syq);

  if ( !check(e.x1, e.x3, sqr(sy1), sqr(syq), sqr(syq)) ) return 0.0;

  double x2 = 2.0 - e.x1 - e.x3;
  e.yo = 0.5*log((1.0 - e.x1 + sqr(sy1) - sqr(2.0*syq))/
		 (1.0 - x2 + sqr(syq) - sqr(sy1 + syq)));


  return (sqr(1.0 - e.x3 + sqr(syq)) + sqr(1.0 - x2 + sqr(syq)))*
    (sqr(rhomax/W)/(1.0 - e.x1 + sqr(sy1)))*2.0*ymax/yint;


}

bool FSQQEmitter::
perform(const Emission & emission) const {
  QCDDipole & d = dynamic_cast<QCDDipole &>(*emission.dipole);
  const FSQQEmission & e = dynamic_cast<const FSQQEmission &>(emission);

  tParPtr gluon = e.ifl > 0? d.oPart(): d.iPart();
  tParPtr other = e.ifl > 0? d.iPart(): d.oPart();
  e.pold = make_pair(gluon->momentum(), other->momentum());
  try {
    e.genmom = getMomenta(d.sdip(), e.x1, e.x3,
			  other->momentum().mass(), e.mq, e.mq,
			  true, false, rnd(2.0*Constants::pi),
			  other->momentum(), gluon->momentum());
  }
  catch ( ImpossibleKinematics ) {
    return false;
  }

  pair<tParPtr,tParPtr> qqbar = splitGluon(d, gluon, e.ifl);

  ParPtr q = qqbar.first;
  ParPtr qbar = qqbar.second;

  other->momentum() = e.genmom.first;
  if ( e.ifl < 0 ) swap(e.genmom.second, e.genmom.third);

  qbar->momentum() = e.genmom.second;
  q->momentum() = e.genmom.third;
  e.partons.push_back(q);
  e.partons.push_back(qbar);

  FSQQEmPtr eo = new_ptr(e);
  swap(eo->y, eo->yo);
  if ( e.ifl > 0.0 ) {
    qbar->emission(&e);
    q->emission(eo);
  } else {
    q->emission(&e);
    qbar->emission(eo);
  }
  return true;

}

void FSQQEmitter::revert(const Emission & emission) const {
  QCDDipole & d = dynamic_cast<QCDDipole &>(*emission.dipole);
  const FSQQEmission & e = dynamic_cast<const FSQQEmission &>(emission);

  tParPtr gluon = fuseQQBar(d, e.partons[0], e.partons[1], e.od, e.mainParent);

  tParPtr other = e.ifl > 0? e.colourParent: e.antiColourParent;

  gluon->setMomentum(e.pold.first);
  other->setMomentum(e.pold.second);

}

tParPtr FSQQEmitter::
fuseQQBar(QCDDipole & d, tParPtr q, tParPtr qbar, tQCDPtr od, tParPtr gluon) {
  // If we do not have an old gluon, create one.
  if ( gluon ) {
    d.state()->addHadronicFS(gluon);
    gluon->untouch();
  } else {
    gluon =
      d.state()->create(CurrentGenerator()->getParticleData(ParticleID::g), q);
    gluon->setVertex(q->vertex() + qbar->vertex());
  }
  if ( d.iPart() == q ) {
    od->next(&d);
    d.prev(od);
    d.iPart(gluon);
    od->oPart(gluon);
    d.oPart()->untouch();
  } else {
    od->prev(&d);
    d.next(od);
    d.oPart(gluon);
    od->iPart(gluon);
    d.iPart()->untouch();
  }
  gluon->origOCol(q->origOCol());
  gluon->origICol(qbar->origICol());
  if ( d.prev() ) d.prev()->untouch();
  if ( d.next() ) d.next()->untouch();
  d.touch();
  d.state()->forgetParton(q);
  d.state()->forgetParton(qbar);

  return gluon;

}

pair<tParPtr,tParPtr>
FSQQEmitter::splitGluon(QCDDipole & d, tParPtr g, int ifl) {
  pair<tParPtr,tParPtr> ret;
  ParPtr q = ret.first = 
    d.state()->create(CurrentGenerator()->getParticleData(abs(ifl)), g);
  ParPtr qbar = ret.second = 
    d.state()->create(CurrentGenerator()->getParticleData(-abs(ifl)), g);
  if ( d.prev() ) d.prev()->touch();
  if ( d.next() ) d.next()->touch();
  if ( g == d.iPart() ) {
    d.iPart(q);
    d.prev()->oPart(qbar);
    d.prev()->next(tQCDPtr());
    d.prev(tQCDPtr());
  } else {
    d.oPart(qbar);
    d.next()->iPart(q);
    d.next()->prev(tQCDPtr());
    d.next(tQCDPtr());
  }
  d.touch();
  d.iPart()->touch();
  d.oPart()->touch();
  d.state()->remove(g);

  q->origOCol(g->origOCol());
  g->origOCol(tColinePtr());
  qbar->origICol(g->origICol());
  g->origICol(tColinePtr());

  q->setVertex(g->vertex());
  qbar->setVertex(g->vertex());

  return ret;

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FSQQEmitter::persistentOutput(PersistentOStream &) const {}

void FSQQEmitter::persistentInput(PersistentIStream &, int) {}

DescribeClass<FSQQEmitter,EmitterBase>
describeAriadne5FSQQEmitter("Ariadne5::FSQQEmitter", "libAriadne5.so");

void FSQQEmitter::Init() {

  static ClassDocumentation<FSQQEmitter> documentation
    ("The FSQQEmitter class implements the final-state splitting of gluons "
     "into a q-qbar pair.");

}

