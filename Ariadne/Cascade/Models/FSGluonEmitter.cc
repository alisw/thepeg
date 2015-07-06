// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FSGluonEmitter class.
//

#include "FSGluonEmission.h"
#include "FSGluonEmitter.h"
#include "Ariadne/Cascade/QCDDipole.h"
#include "Ariadne/Cascade/DipoleState.h"
#include "Ariadne/Cascade/AriadneHandler.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Utilities/UtilityBase.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Repository/CurrentGenerator.h"


#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Ariadne5;

FSGluonEmitter::FSGluonEmitter() {}

FSGluonEmitter::~FSGluonEmitter() {}

IBPtr FSGluonEmitter::clone() const {
  return new_ptr(*this);
}

IBPtr FSGluonEmitter::fullclone() const {
  return new_ptr(*this);
}


bool FSGluonEmitter::canHandle(const DipoleBase & e) const {
  tcQCDPtr d = dynamic_ptr_cast<tcQCDPtr>(&e);
  if ( !d ) return false;
  if ( d->iPart()->special(d->oPart()) || d->oPart()->special(d->iPart()) )
    return false;
  return true;
}

bool FSGluonEmitter::overrides(const EmitterBase &, DipoleBase &) const {
  return false;
}

double FSGluonEmitter::
gluonEngine(FSGluonEmission & e, Energy rhomin, Energy rhomax,
	    Energy W, double C, double yint, int n1, int n3) {
  rhomax = sqrt(rndsud(2.0*C*yint, sqr(rhomax), sqr(rhomin)));
  if ( rhomax <= rhomin ) return -1.0;
  e.rho = rhomax;
  double ymax = acosh(0.5*W/rhomax);
  e.y = 2.0*ymax*UseRandom::rnd() - ymax;
  e.x1 = 1.0 - rhomax*exp( e.y)/W + e.y1 - e.y3;
  e.x3 = 1.0 - rhomax*exp(-e.y)/W + e.y3 - e.y1;

  if ( !check(e.x1, e.x3, e.y1, 0.0, e.y3) ) return 0.0;

  double weight = (pow(e.x1, n1) + pow(e.x3, n3))/2.0;
  weight *= 2.0*ymax/yint;

  // Special handling of heavy quarks, the 'dead-cone' effect.
  double xm = exp(2.0*e.y);
  weight *= 1.0 - (e.y1*xm + e.y3/xm)/(e.x1 + e.x3 - 1.0);

  return weight;

}

EmPtr FSGluonEmitter::
generate(const DipoleBase & dipole, Energy rhomin, Energy rhomax) const {
  const QCDDipole & d = dynamic_cast<const QCDDipole &>(dipole);

  double C = (d.iPart()->isG() || d.oPart()->isG()? 3.0/4.0: 2.0/3.0)/
    Constants::pi;
  int n1 =  d.iPart()->isG()? 3: 2;
  int n3 =  d.oPart()->isG()? 3: 2;
  Energy2 S = d.sdip();
  if ( S <= sqr(2.0*rhomin) ) return EmPtr();
  Energy W = sqrt(S);
  double y1 = max(d.iPart()->momentum().mass2()/S, 0.0);
  double y3 = max(d.oPart()->momentum().mass2()/S, 0.0);
  FSGluonEmission e(*this, dipole, y1, y3);

  C *= preweight(e);

  rhomax = min(rhomax, W/2.0);
  double yint = 2.0*acosh(0.5*W/rhomin);

  while (true) {
    if ( rhomax <= rhomin ) return EmPtr();
    double weight = gluonEngine(e, rhomin, rhomax, W, C, yint, n1, n3);
    if ( weight < 0.0 ) return EmPtr();
    rhomax = e.rho;
    if ( weight == 0.0 ) continue;
    weight *= reweight(e);
    if ( weight > UseRandom::rnd() ) break;
  }

  // Save information in the Emission object.
  e.ymax = log(W/Current<AriadneHandler>()->pTCut());
  e.radiators.push_back(d.oPart());
  e.radiators.push_back(d.iPart());
  e.colourParent = d.iPart();
  e.mainParent = e.antiColourParent = d.oPart();
  e.pold = make_pair(d.iPart()->momentum(), d.oPart()->momentum());
  // The main parent is the one which looses most energy.
  if ( sqr(e.x3) > UseRandom::rnd()*(sqr(e.x1) + sqr(e.x3)) ) {
    e.mainParent = d.iPart();
    swap(e.radiators[0], e.radiators[1]);
  }

  return new_ptr(e);

}


bool FSGluonEmitter::
perform(const Emission & emission) const {
  QCDDipole & d = dynamic_cast<QCDDipole &>(*emission.dipole);
  const FSGluonEmission & e = dynamic_cast<const FSGluonEmission &>(emission);

  tParPtr ip = d.iPart();
  tParPtr op = d.oPart();
  e.pold = make_pair(ip->momentum(), op->momentum());
  try {
    e.genmom = getMomenta(d.sdip(), e.x1, e.x3,
		   ip->momentum().mass(), ZERO, op->momentum().mass(),
		   ip->isG(), op->isG(), rnd(2.0*Constants::pi),
		   ip->momentum(), op->momentum());
  }
  catch ( ImpossibleKinematics ) {
    return false;
  }

  LorentzRotation R = Utilities::getBoostFromCM(make_pair(ip->momentum(),
							  op->momentum()));
 

  ip->momentum() = e.genmom.first;
  op->momentum() = e.genmom.third;

  tParPtr g = insertGluon(d, e.mainParent).first;

  g->momentum() = e.genmom.second;
  g->setVertex(R.inverse()*((sqr(e.x3)*(R*ip->vertex())
			     + sqr(e.x1)*(R*op->vertex()))/
			    (sqr(e.x3) + sqr(e.x1))));

  g->emission(&e);
  e.partons.push_back(g);

  return true;
}

pair<ParPtr, QCDPtr> FSGluonEmitter::
insertGluon(QCDDipole & d, tParPtr p)  {
  pair<ParPtr, QCDPtr> ret;
  tParPtr g = ret.first = 
    d.state()->create(CurrentGenerator()->getParticleData(ParticleID::g), p);
  tQCDPtr newd = ret.second = d.state()->create<QCDDipole>();
  newd->colourIndex(d.colourIndex());
  if ( !g || !newd )
    Throw<Exception>()
      << "Ariadne could not create new dipoles or partons."
      << Exception::abortnow;

  tParPtr ip = d.iPart();
  tParPtr op = d.oPart();
  ip->touch();
  op->touch();
  if ( d.next() ) d.next()->touch();
  if ( d.prev() ) d.prev()->touch();
  d.touch();

  newd->iPart(g);
  newd->oPart(op);
  d.oPart(g);
  if ( d.next() ) d.next()->prev(newd);
  newd->next(d.next());
  d.next(newd);
  newd->prev(&d);
  if ( op == p ) {
    newd->generateColourIndex();
    g->origICol(op->origICol());
    op->origICol(tColinePtr());
  } else {
    d.generateColourIndex();
    g->origOCol(ip->origOCol());
    ip->origOCol(tColinePtr());
  }

  return ret;
}

void FSGluonEmitter::removeGluon(QCDDipole & d, tParPtr g, tParPtr p) {
  tQCDPtr newd = d.next();
  if ( p == d.iPart() ) {
    d.colourIndex(newd->colourIndex());
    d.iPart()->origOCol(g->origOCol());
  } else {
    newd->oPart()->origICol(g->origICol());
  }
  d.oPart(newd->oPart());
  d.next(newd->next());
  if ( d.next() ) d.next()->prev(&d);
  d.iPart()->untouch();
  d.oPart()->untouch();
  if ( d.next() ) d.next()->untouch();
  if ( d.prev() ) d.prev()->untouch();

  d.state()->forgetParton(g);
  d.state()->forgetDipole(newd);
  
    
}

void FSGluonEmitter::revert(const Emission & emission) const {
  QCDDipole & d = dynamic_cast<QCDDipole &>(*emission.dipole);
  const FSGluonEmission & e = dynamic_cast<const FSGluonEmission &>(emission);

  removeGluon(d, d.oPart(), e.mainParent);

  d.iPart()->setMomentum(e.pold.first);
  d.oPart()->setMomentum(e.pold.second);

}




// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FSGluonEmitter::persistentOutput(PersistentOStream &) const {}

void FSGluonEmitter::persistentInput(PersistentIStream &, int) {}

DescribeClass<FSGluonEmitter,EmitterBase>
describeAriadne5FSGluonEmitter("Ariadne5::FSGluonEmitter", "libAriadne5.so");

void FSGluonEmitter::Init() {

  static ClassDocumentation<FSGluonEmitter> documentation
    ("The FSGluonEmitter class implements the standard classical gluon "
     "emission from a final-state colour dipole.");

}

