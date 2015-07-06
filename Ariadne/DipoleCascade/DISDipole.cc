// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DISDipole class.
//

#include "DISDipole.h"
#include "HardRemnant.h"
#include "MECorrBase.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Utilities/UtilityBase.h"
#include "ThePEG/Repository/EventGenerator.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "DISDipole.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Ariadne;

using namespace ThePEG;

ClonePtr DISDipole::clone() const {
  return new_ptr(*this);
}

void DISDipole::init(tParPtr ip, tParPtr op, bool respectScale) {
  ExtendedDipole::init(ip, op, respectScale);
  hardRem(dynamic_ptr_cast<tHardRemPtr>(iSoftRem() ? op : ip));
  softRem(iSoftRem() ? iSoftRem() : oSoftRem());
}

Energy2 DISDipole::generate(Energy2 pt2min, Energy2 pt2max) {
  if(state()->hardSubSys().touched()){
    touch();
  }
  if ( !iPart()->touched() && !oPart()->touched() && !touched() ) {
    return lastPT2();
  }
  Dipole::generate(pt2min, pt2max);
  genDISG(max(lastPT2(), pt2min), pt2max);
  //genBGF(max(lastPT2(), pt2min), pt2max);
  return lastPT2();
}

tParPtr DISDipole::perform() {

  tParPtr emitted = ExtendedDipole::perform();
  if ( !emitted ) return emitted;

  //replace DISDip with ExtDip.
  ExDipPtr d = state()->create<ExtendedDipole>();
  d->init(iPart(), oPart());
  iPart()->oDip(d);
  oPart()->iDip(d);
  d->touch();
  state()->removeEmitter(this);
  return emitted;

}

void DISDipole::genQ(Energy2 pt2min, Energy2 pt2max) {}

void DISDipole::genG(Energy2 pt2min, Energy2 pt2max) {
  double C = preweight(this, state(), ParticleID::g, FSGluon())*2.0/(3.0*Constants::pi);
  Energy2 S = sdip();
  Energy2 Q2 = hardRem()->Q2();
  double nu = y();
  double y3 = max(hardRem()->momentum().mass2()/S, 0.0);
  pt2max = min(pt2max, S/4.0);
  if ( pt2max <= pt2min ) return;
  double yint = 2.0*acosh(sqrt(S/pt2min)/2.0);

  //Vector containing x1, x3 and the rotation angle.
  vector<double> var(3);
  var[2] = handler()->rnd(2.0*Constants::pi);

  while (true) {
    pt2max = rndsud(2.0*C*yint, pt2max, pt2min);
    if ( pt2max <= pt2min ) return;
    double ymax = acosh(sqrt(S/pt2max)/2.0);
    double y = 2.0*ymax*handler()->rnd() - ymax;
    var[0] = 1.0 - sqrt(pt2max/S)*exp( y) - y3;
    var[1] = 1.0 - sqrt(pt2max/S)*exp(-y) + y3;

    if ( !check(var[0], var[1], 0.0, 0.0, y3) ) continue;

    //Check if the gluon has smaller light cone momenta than the hard
    //sub system.
    Lorentz5Momentum p1;
    Lorentz5Momentum p2;
    Lorentz5Momentum p3(hardRem()->momentum().mass());
    try {
      SimplePhaseSpace::CMS(p1, p2, p3, S, var[0], var[1], 0.0, 0.0, 0.0);
    }
    catch ( ImpossibleKinematics ) {
      continue;
    }
    //Check so that the gluon is not too far from the hard remnant.
    if(p2.plus() > hardRem()->Q2() / sqrt(S)){
      continue;
    }

    double z = Q2 / (S * (1 - var[1]) + Q2);
    double xi = 1 - (1 - var[0]) / var[1];
    //matrix element multiplied by jacobi determinant over 2z.
    double weight = ( sqr(z) + sqr(xi) + (2 + xi*z * (2 + 8*(1.0-nu)
          / (1 + sqr(1.0 - nu) ) ) ) * (1.0 - xi) * (1.0 - z) ) / 2.0;
    weight *= 2.0*ymax/yint;

    weight *= reweight(this, state(), ParticleID::g, pt2max, var, FSGluon());

    if ( weight < handler()->rnd() ) continue;

    lastPT2(pt2max);
    flav = ParticleID::g;
    genVar.swap(var);
    return;
  }
}

void DISDipole::genDISG(Energy2 pt2min, Energy2 pt2max) {
  double C = preweight(this, state(), ParticleID::g, RFGluon())*2.0/(3.0*Constants::pi);
  Energy2 S = (softRem()->momentum() +
      hardRem()->momentum()).m2();
  Energy2 mq2 = hardRem()->momentum().mass2();
  Energy W = sqrt(S);
  Energy2 Q2 = hardRem()->Q2();
  double nu = y();
  pt2max = min(pt2max, S/4.0);
  if ( pt2max <= pt2min ) return;
  double yint = 2.0*acosh(sqrt(S/pt2min)/2.0);

  //Vector containing rapidity.
  vector<double> var(1);

  while (true) {
    pt2max = rndsud(2.0*C*yint, pt2max, pt2min);
    if ( pt2max <= pt2min ) return;
    double ymax = acosh(sqrt(S/pt2max)/2.0);
    var[0] = 2.0*ymax*handler()->rnd() - ymax;

    //double aa1 = pow(sqr(softRem()->mu())/pt2max, softRem()->alpha()/2);
    //double aa3 = pow(sqr(hardRem()->mu())/pt2max, hardRem()->alpha()/2);

    //Calculations below are done assuming that the soft remnant is in
    //positive y direction.

    //double extrat = (S/sqr(exp(-var[0])/aa3 + exp(var[0])/aa1))/pt2max;
    double extrat = MaxPT(S, var[0], softRem()->mu(), hardRem()->mu(),
        softRem()->alpha(), hardRem()->alpha()) / sqrt(pt2max);
    //double extrat = sqr(softRem()->mu())/pt2max;

    double beta = handler()->beta();
    if ( beta < 0.0 &&  extrat < 1.0 ) continue;

    Energy pt = sqrt(pt2max);

    //Check so that the gluon is far enough from the hard remnant.
    if ( pt*exp(var[0]) < Q2 / W ) continue;
    
    //Check kinematics
    if ( pt*exp(-var[0]) >= W ) continue;

    //S of the hard system after the gluon is emitted.
    Energy ppq = (hardRem()->momentum().mass2() + pt2max) / 
      (W - pt*exp(-var[0]));
    Energy2 sh = W*(ppq + pt*exp(var[0]));

    //Check if there is enough phase space availible.
    if( sh >= S ) continue;

    //Check so that en inviraint mass of the gluon and the remnant is
    //larger than 4*mq^2, where mq is the mass of the extracted quark.
    if( (W - ppq) * pt * exp(-var[0]) - pt2max <
      4 * sqr(softRem()->data().generateMass()) ) {
      continue;
    }

    double z = (hardRem()->momentum().mass2() + Q2) / (sh + Q2);
    double xi = 1.0 - sqrt(pt2max / S) * exp( -var[0] );

    //matrix element multiplied by (1-z)(1-xi)/2
    double weight = ( sqr(z) + sqr(xi) + (2 + xi*z * (2 + 8*(1.0-nu)
          / (1 + sqr(1.0 - nu) ) ) ) * (1.0 - xi) * (1.0 - z) )/2.0;
    //Jacobi determinant divided by z(1-z)(1-xi)
    weight *= S / (S + mq2 * exp(-2*var[0]));
    weight *= 2.0*ymax/yint;

    // Probability to have emissions above the soft cutoff.
    if ( beta >= 0.0 ) {
      weight /= (1.0 + pow(extrat, -beta));
    }

    weight *= reweight(this, state(), ParticleID::g, pt2max, var, RFGluon());

    if ( weight < handler()->rnd() ) continue;

    lastPT2(pt2max);
    flav = ParticleID::g;
    genVar.swap(var);
    lastEType(RFGluon());
    return;
  }
}

void DISDipole::genBGF(Energy2 pt2min, Energy2 pt2max){
  if( softRem()->isG() || !softRem()->data().coloured() ){
    return;
  }

  HardSubSys & hardsys = state()->hardSubSys();
  Energy mq = hardsys.momentum().mass();
  Energy2 mq2 = sqr(mq);
  double x = softRem()->partonx();
  Energy2 Q2 = hardsys.Q2();
  double nu = y();

  //Calculate sdip and maximum possible pt2
  Energy2 S = (softRem()->momentum() + hardsys.momentum()).m2();
  //Use mt as the scale of the process.
  Energy2 mt2km = S / 4.0;
  //Lower mt2max to the maximum possible kinimatically
  Energy2 mt2max = min(pt2max, mt2km);
  Energy2 mt2min = max(pt2min, sqr(handler()->pTCut()) + mq2);
  if(mt2max <= mt2min){
    return;
  }

  int ifl = -softRem()->data().id();
  tcPDPtr q = handler()->generator()->getParticleData(ifl);

  //Use the geometric average of max and min as scale for the pdf.
  //If mt2max is smaller use that value.
  Energy2 mt2z = min( mt2max, sqrt( mt2km * mt2min ) );
  Energy2 mt2mean = sqrt( mt2z * mt2min );
  //Calculate zmax for pt2z and pt2min
  double zmax = (mq2 + Q2)/ (4.0 * mt2min + Q2);
  double zmax1 = (mq2 + Q2) / (4.0 * mt2mean + Q2);
  double zmax2 = (mq2 + Q2) / (4.0 * mt2z + Q2);
  tcPDPtr g = handler()->generator()->getParticleData(ParticleID::g);

  //multiply pdf by a safty factor of two
  //use zero if zmax is less than x
  double maxpdf = 2.0 * max(
      PDFRatio(softRem(), g, softRem()->dataPtr(), mt2z, x, zmax2),
      max( 
        PDFRatio(softRem(), g, softRem()->dataPtr(), mt2mean, x, zmax1),
        PDFRatio(softRem(), g, softRem()->dataPtr(), mt2min, x, zmax) ) );
  if( maxpdf <= 0.0){
    return;
  }

  //xi is the fraction of the positive light cone momentum taken by the
  //hard sub system. x < xi < 1
  double C = preweight(this, state(), ifl, ISGtoQQ())/(4.0*Constants::pi);
  double yint = 2.0*acosh(sqrt(S/(pt2min+mq2))/2.0);

  //Vector containing xi, z and the quark mass in GeV.
  vector<double> var(3);
  var[2] = mq/GeV;

  while ( true ) {
    mt2max = rndsud(2.0 * C * maxpdf * yint, mt2max, mt2min);
    if ( mt2max <= mt2min ) break;

    double ymax = acosh(sqrt(S/mt2max)/2.0);
    double y = 2.0*ymax*handler()->rnd() - ymax;

    Energy2 sh = mt2max*sqr(2*cosh(y));
    //check if y and pt2 is kinematically possible
    if( sh > S ){
      continue;
    }

    var[1] = (mq2 + Q2)/(sh + Q2);
    var[0] = 1.0 / (1.0 + exp(-2*y));

    double weight = (sqr(var[1]) + sqr(1 - var[1])) *
      (sqr(var[0]) + sqr(1 - var[0])) / (var[0] * (1 - var[0])) +
      16 * var[1] * (1 - var[1]) * (1 - nu) / (1 + sqr(1 - nu));
    weight *= 
      PDFRatio(softRem(), g, softRem()->dataPtr(), mt2max, x, var[1])
      / maxpdf;
    weight *= 2.0*ymax/yint;

    //multiply by jacobi determinant / 2.0
    weight *= var[1] * mt2max / (sh + Q2);

    weight *= reweight(this, state(), ifl, mt2max - mq2, var, ISGtoQQ());

    if ( weight > 1.0 ) handler()->generator()->logWarning(
        WeightException() << "Ariadne::DISDipole::genBGF "
        << "failed to overestimate the PDF ratio. "
        << "If this hapens too often you should contact the authors."
        << Exception::warning);			   

    if ( weight < handler()->rnd() ) continue;
      
    lastPT2(mt2max);
    flav = ifl;
    genVar.swap(var);
    //set the emission type
    lastEType(ISGtoQQ());
    break;
  }
}

tParPtr DISDipole::performInitialQ(){
  tSoftRemPtr rem;
  if(flav < 0){
    rem = iSoftRem();
  }
  else{
    rem = oSoftRem();
  }
  if(!rem){
    return tParPtr();
  }

  Lorentz5Momentum pr = rem->momentum();
  Lorentz5Momentum ph = state()->hardSubSys().momentum();

  //boost to cm
  LorentzRotation rcm = Utilities::boostToCM(make_pair(&pr, &ph));

  Energy mq = genVar[2]*GeV;
  double xi = genVar[0];
  double z = genVar[1];
  Energy2 Q2 = state()->hardSubSys().Q2();
  Energy2 mh2 = ph.mass2();
  double phi = 2.0 * Constants::pi * handler()->rnd();
  //lastPT2 returns mt2.
  Energy pt = sqrt(lastPT2() - sqr(mq));
  Transverse<Energy> tpt( pt * cos(phi), pt * sin(phi) );

  //pplus of q plus hard sub sys.
  Energy pplus = (1.0/z * (mh2 + Q2) - Q2) / ph.minus();
  Energy pminus = ph.minus();

  Energy dpp = pplus - ph.plus();
  if(dpp > pr.plus()){
    return tParPtr();
  }
  pr = lightCone(pr.plus() - dpp, pr.minus(), pr.x(), pr.y());

  Lorentz5Momentum pq = lightCone((1-xi) * pplus,
      lastPT2() / ((1-xi) * pplus), tpt);
  ph = lightCone(pplus - pq.plus(), pminus - pq.minus(), -tpt);
  pq.setMass(mq);

  //boost to lab frame
  rcm.invert();
  ph.transform(rcm);
  pq.transform(rcm);
  pr.transform(rcm);

  state()->hardSubSys().setMomentum(ph, pr.z() > 0.0*GeV);
  state()->hardSubSys().touch();

  rem->momentum() = pr;
  rem->touch();

  //create quark
  ParPtr q = state()->create(handler()->generator()->getParticleData(flav));
  q->momentum() = pq;
  q->parents( make_pair(iPart(), oPart()) );
  q->touch();

  return q;
}

double DISDipole::y(){
  if(theY < 0.0){
    LorentzMomentum P = softRem()->parentMomentum();
    LorentzMomentum q = state()->totalMomentum() - P;
    tcPPair leptons = state()->scatteredLeptons();
    LorentzMomentum l = leptons.first ? leptons.first->momentum() :
      leptons.second->momentum();
    theY = P.dot(q) / P.dot(l + q);
  }
  return theY;
}

DISDipole::~DISDipole() {}

void DISDipole::fillReferences(CloneSet & cset) const {
  ExtendedDipole::fillReferences(cset);
  cset.insert(softRem());
  cset.insert(hardRem());
}

void DISDipole::rebind(const TranslationMap & trans) {
  ExtendedDipole::rebind(trans);
  theSoftRem = trans.translate(softRem());
  theHardRem = trans.translate(hardRem());
}

void DISDipole::persistentOutput(PersistentOStream & os) const {
  os << theHardRem << theSoftRem;
}

void DISDipole::persistentInput(PersistentIStream & is, int) {
  is >> theHardRem >> theSoftRem;
}

void DISDipole::debugme() const {
  ExtendedDipole::debugme();
  cerr << "dis";
}

ClassDescription<DISDipole> DISDipole::initDISDipole;
// Definition of the static class description member.

void DISDipole::Init() {}

