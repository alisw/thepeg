// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Dipole class.
//

#include "ExtendedDipole.h"
#include "Parton.h"
#include "String.h"
#include "DipoleState.h"
#include "CascadeHandler.h"
#include "MECorrBase.h"
#include "ThePEG/PDT/EnumParticles.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "Dipole.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Utilities/UtilityBase.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

using namespace Ariadne;

ClonePtr Dipole::clone() const {
  return new_ptr(*this);
}

void Dipole::init(tParPtr ip, tParPtr op, bool respectScale) {
  iPart(ip);
  oPart(op);
  touch();
  reset(sdip()/4.0);
  if ( !respectScale ) return;
  Energy2 maxscale = sdip()/4.0;
  if ( ip->orig() ) maxscale = min(maxscale, ip->orig()->scale());
  if ( op->orig() ) maxscale = min(maxscale, op->orig()->scale());
  maxScale(maxscale);
}

void Dipole::reset(Energy2 pt2max) {
  setup(pt2max);
  lastPT2(-1.0*GeV2);
  resonance(tcPDPtr());
}

Energy2 Dipole::generate(Energy2 pt2min, Energy2 pt2max) {
  if ( !iPart()->touched() && !oPart()->touched() && !touched() ) {
    return lastPT2();
  }
  reset(pt2max);
  genG(pt2min, pt2max);
  genQ(max(lastPT2(), pt2min), pt2max);
  return lastPT2();
}

void Dipole::genQ(Energy2 pt2min, Energy2 pt2max) {
  Energy2 S = sdip();
  double C1 = 0.0;
  if ( iPart()->isG() ) C1 = 1.0/(4.0*Constants::pi*(1.0 + S/prev()->sdip()));
  double C3 = 0.0;
  if ( oPart()->isG() ) C3 = 1.0/(4.0*Constants::pi*(1.0 + S/next()->sdip()));
  if ( C1 + C3 <= 0.0 ) return;

  // We will pretend parton 3 is the gluon even if it isn't. Hence
  // parton 1 may have a mass or not.
  double y1 = 0.0;
  if ( !iPart()->isG() ) y1 = max(iPart()->momentum().mass2()/S, 0.0);
  if ( !oPart()->isG() ) y1 = max(oPart()->momentum().mass2()/S, 0.0);
  double sy1 = sqrt(y1);
  Energy W = sqrt(S);
  pt2max = min(pt2max, S/4.0);
  if ( pt2max <= pt2min ) return;

  //Vector containing x1, x3 and the quark mass in GeV.
  vector<double> var(3);

  for ( int ifl = 1; ifl <= handler()->nFlav(); ++ifl ) {

    tcPDPtr q = handler()->generator()->getParticleData(ifl);
    if ( !q ) continue;
    Energy mq = q->generateMass();
    double syq = mq/W;
    double yq = sqr(syq);
    var[2] = mq/GeV;

    Energy2 pt2m = pt2max;
    if ( pt2m <= pt2min ) continue;
    double yint = 2.0*acosh(sqrt(S/pt2min)/2.0);

    double CW = preweight(this, state(), ifl, FSGtoQQ());

    while ( true ) {
      pt2m = rndsud((C1 + C3)*yint*CW, pt2m, pt2min);
      if ( pt2m <= pt2min ) break;
      double ymax = acosh(sqrt(S/pt2m)/2.0);
      double y = 2.0*ymax*handler()->rnd() - ymax;
      var[0] = 1.0 - sqrt(pt2m/S)*exp( y) + y1 - 4.0*yq;
      var[1] = 1.0 - sqrt(pt2m/S)*exp(-y) + yq - sqr(sy1 + syq);
      double x2 = 2.0 - var[0] - var[1];

      if ( !check(var[0], var[1], y1, yq, yq) ) continue;

      double weight = (sqr(1.0 - var[1] + yq) + sqr(1.0 - x2 + yq))*
	((pt2m/S)/(1.0 - var[0] + y1))*2.0*ymax/yint;

      // *** ATTENTION *** Add other weights here.

      weight *= reweight(this, state(), ifl, pt2m, var, FSGtoQQ());

      if ( weight < handler()->rnd() ) continue;
      
      lastPT2(pt2m);
      lastEType(FSGtoQQ());
      flav = handler()->rndbool(C1, C3)? -ifl: ifl;
      genVar = var;
      pt2min = pt2m;
      break;
    }
  }

}

void Dipole::genG(Energy2 pt2min, Energy2 pt2max) {
  double C = (iPart()->isG() || oPart()->isG()? 3.0/4.0: 2.0/3.0)/
    Constants::pi;
  C *= preweight(this, state(), ParticleID::g, FSGluon());
  int n1 =  iPart()->isG()? 3: 2;
  int n3 =  oPart()->isG()? 3: 2;
  Energy2 S = sdip();
  double y1 = max(iPart()->momentum().mass2()/S, 0.0);
  double y3 = max(oPart()->momentum().mass2()/S, 0.0);
  pt2max = min(pt2max, S/4.0);
  if ( pt2max <= pt2min ) return;
  double yint = 2.0*acosh(sqrt(S/pt2min)/2.0);

  //Vector containing x1 and x3
  vector<double> var(2);
  while (true) {
    pt2max = rndsud(2.0*C*yint, pt2max, pt2min);
    if ( pt2max <= pt2min ) return;
    double ymax = acosh(sqrt(S/pt2max)/2.0);
    double y = 2.0*ymax*handler()->rnd() - ymax;
    var[0] = 1.0 - sqrt(pt2max/S)*exp( y) + y1 - y3;
    var[1] = 1.0 - sqrt(pt2max/S)*exp(-y) + y3 - y1;

    if ( !check(var[0], var[1], y1, 0.0, y3) ) continue;

    double weight = (pow(var[0], n1) + pow(var[1], n3))/2.0;
    weight *= 2.0*ymax/yint;

    weight *= reweight(this, state(), ParticleID::g, pt2max, var, FSGluon());

    //Special handling of heavy quarks, the 'dead-cone' effect.
    double xm = exp(2.0*y);
    weight *= 1.0 - (y1*xm + y3/xm)/(var[0] + var[1] - 1.0);

    if ( weight < handler()->rnd() ) continue;

    lastPT2(pt2max);
    lastEType(FSGluon());
    flav = ParticleID::g;
    genVar.swap(var);
    return;
  }
}

tParPtr Dipole::perform() {
  if ( !flav ) return tParPtr();
  if ( flav == ParticleID::g )  return performG();
  else return  performQ();
}

tParPtr Dipole::performQ() {
  //genVar contains x1 and x3
  tParPtr gluon = flav > 0? oPart(): iPart();
  tParPtr other = flav > 0? iPart(): oPart();
  Lorentz5Momentum p1(other->momentum().mass());
  Lorentz5Momentum p2(genVar[2]*GeV);
  Lorentz5Momentum p3(genVar[2]*GeV);

  try {
    SimplePhaseSpace::CMS(p1, p2, p3, sdip(), genVar[0], genVar[1], 0.0, 0.0, 0.0);
  }
  catch ( ImpossibleKinematics ) {
    return tParPtr();
  }

  LorentzRotation R;
  R.rotateZ(handler()->rnd(2.0*Constants::pi));
  R.transform(Utilities::getBoostFromCM(make_pair(other->momentum(),
						  gluon->momentum())));
  p1.transform(R);
  p2.transform(R);
  p3.transform(R);

  other->momentum() = p1;

  return splitGluon(p2, p3);

}

tParPtr Dipole::performG() {
  //genVar contains x1 and x3
  Lorentz5Momentum p1(iPart()->momentum().mass());
  Lorentz5Momentum p2;
  Lorentz5Momentum p3(oPart()->momentum().mass());

  try {
    SimplePhaseSpace::CMS(p1, p2, p3, sdip(), genVar[0], genVar[1], 0.0, 0.0, 0.0);
  }
  catch ( ImpossibleKinematics ) {
    return tParPtr();
  }

  double Psi = Constants::pi - p3.theta();
  double beta = 0.0;
  if ( iPart()->isG() && oPart()->isG() )
    beta = Psi*sqr(genVar[1])/(sqr(genVar[0]) + sqr(genVar[1])); // minimize pt
  else if ( oPart()->isG() )
    beta = Psi;
  else if ( !iPart()->isG() && sqr(genVar[1]) > handler()->rnd()*
      (sqr(genVar[0]) + sqr(genVar[1])) )
    beta = Psi;

  LorentzRotation R;
  R.rotateY(-beta);
  R.rotateZ(handler()->rnd(2.0*Constants::pi));
  R.transform(Utilities::getBoostFromCM(make_pair(iPart()->momentum(),
						  oPart()->momentum())));
  p1.transform(R);
  p2.transform(R);
  p3.transform(R);

  iPart()->momentum() = p1;
  oPart()->momentum() = p3;
  return emitGluon(p2);
}

tParPtr Dipole::splitGluon(const Lorentz5Momentum & p2,
			const Lorentz5Momentum & p3) {
  ParPtr q =
    state()->create(handler()->generator()->getParticleData(abs(flav)));
  ParPtr qbar =
    state()->create(handler()->generator()->getParticleData(-abs(flav)));
  touch();
  if ( flav > 0 ) {
    q->parents(make_pair(oPart(), iPart()));
    qbar->parents(make_pair(oPart(), iPart()));
    q->momentum() = p2;
    qbar->momentum() = p3;
    next()->touch();
    iPart()->touch();
    oPart()->string()->split(oPart(), q, qbar, true);
    return qbar;
  } else {
    q->parents(make_pair(iPart(), oPart()));
    qbar->parents(make_pair(iPart(), oPart()));
    qbar->momentum() = p2;
    q->momentum() = p3;
    prev()->touch();
    oPart()->touch();
    iPart()->string()->split(iPart(), q, qbar, true);
    return q;
  }
}

ParPtr Dipole::emitGluon(const Lorentz5Momentum & pg) {
  ParPtr g =
    state()->create(handler()->generator()->getParticleData(ParticleID::g));
  DipPtr d = state()->create<Dipole>();
  if ( !g || !d ) return g;
  g->momentum() = pg;
  iPart()->touch();
  oPart()->touch();
  touch();
  g->parents(make_pair(iPart(), oPart()));
  oPart()->iDip(d);
  d->init(g, oPart());
  g->oDip(d);
  g->iDip(this);
  oPart(g);

  g->string(iPart()->string());
  d->generateColourIndex();

  return g;
}

double Dipole::emissionProbability(){
  Energy2 S = sdip();
  double y1 = max(iPart()->momentum().mass2()/S, 0.0);
  double y3 = max(oPart()->momentum().mass2()/S, 0.0);
  double scale = sqr(handler()->pTCut()) / lastPT2();
  double me = 1.0;
  if ( flav == ParticleID::g )
    me *= preweight(this, state(), flav, FSGluon()) *
      reweight(this, state(), flav, lastPT2(), genVar, FSGluon());
  else
    me *= preweight(this, state(), flav, FSGtoQQ()) *
      reweight(this, state(), flav, lastPT2(), genVar, FSGtoQQ());

  if(flav == ParticleID::g){
    double C = (iPart()->isG() || oPart()->isG()? 3.0/4.0: 2.0/3.0)/
      Constants::pi;
    int n1 =  iPart()->isG()? 3: 2;
    int n3 =  oPart()->isG()? 3: 2;
    //dead cone
    double xm = (1 - genVar[0] + y1 - y3) / (1 - genVar[1] + y3 - y1);
    double dcw = max(0.0, 1.0 - (y1*xm + y3/xm)/(genVar[0] + genVar[1] - 1.0));
    return C*dcw*(pow(genVar[0], n1) + pow(genVar[1], n3)) * scale * me;
  }
  else{
    double C = flav > 0 ? S/(4.0*Constants::pi*(S + next()->sdip())):
      S/(4.0*Constants::pi*(S + prev()->sdip()));
    double yp = max(y1,y3);
    double yq = genVar[2]*GeV2 / S;
    double x2 = 2.0 - genVar[0] - genVar[1];
    return C*(sqr(1.0 - genVar[1] + yq) + sqr(1.0 - x2 + yq))/
      (1.0 - genVar[0] + yp) * lastPT2()/S * scale * me;
  }
}

double Dipole::coupling(){
  setup(state()->constructedPT2());
  switch(handler()->runningCoupling()){
  case CascadeHandler::externalRunning:
    return handler()->standardModel()->alphaS(lastPT2());
  case CascadeHandler::noRunning:
    return alpha0();
  case CascadeHandler::simpleRunning:
    return -alpha0()/log(lambdaQCD2()/lastPT2());
  case CascadeHandler::internalRunning:
    return handler()->internalAlphaS()->value(lastPT2(), handler()->SM());
  }
  return 0.0;
}

Emitter::DipoleStateVector Dipole::constructStep(){
  DipoleStateVector vec;
  if(oPart()->isG() && ! dynamic_ptr_cast<tExDipPtr> (next()) ){
    TranslationMap trans;
    DipoleStatePtr dipstate = state()->fullclone(trans);
    tDipPtr newdip = trans.translate(tDipPtr(this));
    dipstate->selected(newdip);
    if(newdip->undoOG()){
      vec.push_back(dipstate);
    }
  }

  if(!oPart()->isG() || !iPart()->isG()){
    const HardSubSys::PartonSet partset = state()->hardSubSys().active();
    for(HardSubSys::PartonSet::const_iterator it = partset.begin();
        it != partset.end(); ++it ){
      if( (iPart()->data().id() == -(*it)->data().id() && 
            iPart()->string() != (*it)->string()) ||
          (oPart()->data().id() == -(*it)->data().id() && 
           oPart()->string() != (*it)->string()) ){
        TranslationMap trans;
        DipoleStatePtr dipstate = state()->fullclone(trans);
        tDipPtr newdip = trans.translate(tDipPtr(this));
        tParPtr q = trans.translate(*it);
        dipstate->selected(newdip);
        if(newdip->undoQ(q)){
          vec.push_back(dipstate);
        }
      }
    }
  }
  return vec;
}

bool Dipole::undoOG(){
  genVar.clear();
  tDipPtr dip = next();
  if(!(dip && oPart()->isG())){
    return false;
  }
  tParPtr nPart = dip->oPart();

  LorentzMomentum p1 = iPart()->momentum();
  LorentzMomentum p2 = oPart()->momentum();
  LorentzMomentum p3 = nPart->momentum();

  LorentzRotation R = Utilities::boostToCM(makeTriplet(&p1, &p2, &p3));

  Energy2 S = (p1 + p2 + p3).m2();
  Energy W = sqrt(S);
  genVar.push_back(2.0*p1.e()/W);
  genVar.push_back(2.0*p3.e()/W);
  flav = ParticleID::g;

  double Psi = Constants::pi - p3.theta();
  double beta = 0.0;
  if ( (iPart()->isG() && nPart->isG()) || (!iPart()->isG() && !nPart->isG()) )
    beta = Psi*sqr(genVar[1])/(sqr(genVar[0]) + sqr(genVar[1])); // minimize pt
  else if ( nPart->isG() )
    beta = Psi;

  R.rotateY(-beta);
  R.invert();

  Lorentz5Momentum p1n(iPart()->momentum().mass());
  Lorentz5Momentum p3n(nPart->momentum().mass());

  try{
    SimplePhaseSpace::CMS(p1n, p3n, S, 1.0, 0.0);
  }
  catch ( ImpossibleKinematics ) {
    return false;
  }

  lastPT2(S*(1 - genVar[0] + (p1n.mass2() - p3n.mass2())/S)*
	  (1 - genVar[1] + (p3n.mass2() - p1n.mass2())/S));

  p1n.transform(R);
  p3n.transform(R);

  state()->removeEmitter(tEmiPtr(dip));
  state()->remove(oPart());

  oPart(nPart);
  oPart()->iDip(this);

  iPart()->momentum() = p1n;
  oPart()->momentum() = p3n;

  return true;
}

bool Dipole::undoQ(tParPtr q2){
  genVar.clear();
  tParPtr q1;
  tParPtr other;
  if(iPart()->data().id() == -q2->data().id()){
    q1 = iPart();
    other = oPart();
  }
  else if(oPart()->data().id() == -q2->data().id()){
    q1 = oPart();
    other = iPart();
  }
  else{
    return false;
  }

  LorentzMomentum p1 = other->momentum();
  LorentzMomentum p2 = q1->momentum();
  LorentzMomentum p3 = q2->momentum();

  LorentzRotation R = Utilities::getBoostFromCM(makeTriplet(&p1, &p2, &p3));

  Energy2 S = (p1 + p2 + p3).m2();
  flav = q1->data().id();
  Energy qmass = (q1->momentum().mass()+q2->momentum().mass())/2;
  genVar.push_back(1.0 - ((p2 + p3).m2() - other->momentum().mass2())/S);
  genVar.push_back(1.0 - ((p2 + p1).m2() - sqr(qmass))/S);
  genVar.push_back(qmass/GeV);
  lastPT2(S*(1 - genVar[0] + (other->momentum().mass2() - 4*sqr(qmass))/S)*
	  (1 - genVar[1] + (sqr(qmass) - sqr(other->momentum().mass()+qmass))/S));

  Lorentz5Momentum p1n(other->momentum().mass());
  Lorentz5Momentum p3n;

  try{
    SimplePhaseSpace::CMS(p1n, p3n, S, 1.0, 0.0);
  }
  catch ( ImpossibleKinematics ) {
    return tParPtr();
  }

  p1n.transform(R);
  p3n.transform(R);

  tParPtr g =
    state()->create(handler()->generator()->getParticleData(ParticleID::g));
  if(!g){
    return g;
  }

  q1->string()->join(q2, g);

  state()->remove(q1);
  state()->remove(q2);

  other->momentum() = p1n;
  g->momentum() = p3n;

  return true;
}

tDipPtr Dipole::prev() const {
  return iPart()? iPart()->iDip(): tDipPtr();
}

tDipPtr Dipole::next() const {
  return oPart()? oPart()->oDip(): tDipPtr();
}

Energy2 Dipole::sdip() const {
  return (iPart()->momentum() + oPart()->momentum()).m2();
}

void Dipole::generateColourIndex() {
  int ncol = handler()->nCol();
  if ( ncol < 3 ) return;
  int idx1 = prev()? prev()->colourIndex(): 0;
  int idx3 = next()? next()->colourIndex(): 0;
  int strndx = 0;
  if ( idx1 ) strndx = idx1/div;
  else if ( idx3 ) strndx = idx3/div;
  do {
    colourIndex(handler()->irnd(ncol) + 1 + strndx*div);
  } while ( colourIndex() == idx1 || colourIndex() == idx3 );
}

void Dipole::fillReferences(CloneSet & cset) const {
  CascadeBase::fillReferences(cset);
  cset.insert(oPart());
  cset.insert(iPart());
}

void Dipole::rebind(const TranslationMap & trans) {
  CascadeBase::rebind(trans);
  theOPart = trans.translate(oPart());
  theIPart = trans.translate(iPart());
}

void Dipole::persistentOutput(PersistentOStream & os) const {
  os << theIPart << theOPart << theColourIndex << theResonance
     << flav;
}

void Dipole::persistentInput(PersistentIStream & is, int) {
  is >> theIPart >> theOPart >> theColourIndex >> theResonance
     >> flav;
}

ClassDescription<Dipole> Dipole::initDipole;
// Definition of the static class description member.

void Dipole::Init() {}

void Dipole::debugme() const {
  Emitter::debugme();
  cerr << "D"
       << setw(3) << state()->index(this)
       << setw(3) << state()->index(iPart())
       << setw(3) << state()->index(oPart())
       << setw(10) << setprecision(3) << sqrt(max(lastPT2()/GeV2, 0.0))
       << (touched()? " *": "  ")
       << (state()->selected() == this? "<-": "  ");
}

bool Dipole::checkIntegrety() {
  return Emitter::checkIntegrety() &&
    iPart() &&
    oPart() &&
    this == iPart()->oDip() &&
    this == oPart()->iDip() &&
    iPart()->string() == oPart()->string();
}

