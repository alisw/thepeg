// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ExtendedDipole class.
//

#include "ExtendedDipole.h"
#include "DISDipole.h"
#include "Ariadne/DipoleCascade/SoftRemnant.h"
#include "Ariadne/DipoleCascade/MECorrBase.h"
#include "Ariadne/DipoleCascade/HardRemnant.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Utilities/UtilityBase.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ExtendedDipole.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Ariadne;

ClonePtr ExtendedDipole::clone() const {
  return new_ptr(*this);
}

void ExtendedDipole::init(tParPtr ip, tParPtr op, bool respectScale) {
  Dipole::init(ip, op, respectScale);
  iSoftRem(dynamic_ptr_cast<tSoftRemPtr>(ip));
  oSoftRem(dynamic_ptr_cast<tSoftRemPtr>(op));
  iHardRem(dynamic_ptr_cast<tHardRemPtr>(ip));
  oHardRem(dynamic_ptr_cast<tHardRemPtr>(op));
}

Energy2 ExtendedDipole::generate(Energy2 pt2min, Energy2 pt2max) {
  if(state()->hardSubSys().touched()){
    touch();
  }
  if ( !iPart()->touched() && !oPart()->touched() && !touched() ) {
    return lastPT2();
  }
  Dipole::generate(pt2min, pt2max);
  genHRemG(max(lastPT2(), pt2min), pt2max);
  genIniG(max(lastPT2(), pt2min), pt2max);
  //genRemG(max(lastPT2(), pt2min), pt2max);
  //genSeaQ(max(lastPT2(), pt2min), pt2max);
  //genQtoGQ(max(lastPT2(), pt2min), pt2max);
  return lastPT2();
}

tParPtr ExtendedDipole::perform() {
  state()->hardSubSys().touch();
  if ( isType(FSGluon()) || isType(FSGtoQQ()) ) return Dipole::perform();
  if ( isType(RFGluon()) ) return performIniG();
  if ( isType(RRGluon()) ) return performRemG();
  if ( isType(HRGluon()) ) return performHRemG();
  if ( isType(ISGtoQQ()) ) return performSeaQ();
  if ( isType(ISQtoGQ()) ) return performQtoGQ();
  return tParPtr();
}

double ExtendedDipole::
newa(double xt, double y, double a0, double beta, double R) {
  double z0 = sqrt(a0);
  double yp = y - log(z0);
  double r = xt*2.0*cosh(yp)/z0;
  double rb = pow(r, beta);
  double zb = pow(z0, -beta);
  double z =
    z0*pow((rb*R + (1.0 + rb - R)*zb)/(1.0 + rb*(1.0 - R) + R*zb),1.0/beta);
  return sqr(z);
}

void ExtendedDipole::genG(Energy2 pt2min, Energy2 pt2max) {
  if(iSoftRem() && oSoftRem()){
    return;
  }
  if( !(iSoftRem() || oSoftRem()) ){
    return;
  }

  tSoftRemPtr rem = iSoftRem()? iSoftRem(): oSoftRem();
  tParPtr part = iSoftRem()? oPart(): iPart();

  double C = (part->isG()? 3.0/4.0: 2.0/3.0)/
    Constants::pi;
  C *= preweight(this, state(), ParticleID::g, FSGluon());
  int n1 =  rem->isG()? 3: 2;
  int n3 =  part->isG()? 3: 2;
  Energy2 S = sdip();
  double y3 = max(part->momentum().mass2()/S, 0.0);
  pt2max = min(pt2max, S/4.0);
  if ( pt2max <= pt2min ) return;
  double yint = 2.0*acosh(sqrt(S/pt2min)/2.0);

  //Vector containing x1, x3 and the rotation angle.
  vector<double> var(3);

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
    Lorentz5Momentum p3(part->momentum().mass());
    try {
      SimplePhaseSpace::CMS(p1, p2, p3, S, var[0], var[1], 0.0, 0.0, 0.0);
    }
    catch ( ImpossibleKinematics ) {
      continue;
    }
    var[2] = handler()->rnd(2.0*Constants::pi);
    p2.rotateZ(var[2]);
    p2.transform(Utilities::getBoostFromCM(make_pair(rem->momentum(),
          part->momentum())));
    if( p2.plus() > state()->hardSubSys().momentum().plus() ||
        p2.minus() > state()->hardSubSys().momentum().minus() ){
      continue;
    }

    double weight = (pow(var[0], n1) + pow(var[1], n3))/2.0;
    weight *= 2.0*ymax/yint;

    weight *= reweight(this, state(), ParticleID::g, pt2max, var, FSGluon());

    // Special handling of heavy quarks, the 'dead-cone' effect.
    double xm = exp(2.0*y);
    weight *= 1.0 - (y3/xm)/(var[0] + var[1] - 1.0);

    if ( weight < handler()->rnd() ) continue;

    lastPT2(pt2max);
    flav = ParticleID::g;
    genVar.swap(var);
    return;
  }
}

void ExtendedDipole::genQ(Energy2 pt2min, Energy2 pt2max) {
  if(iSoftRem() && oSoftRem()){
    return;
  }

  tRemParPtr rem;
  bool irem;
  if(iSoftRem() || oSoftRem()){
    irem = iSoftRem();
    rem = irem? iSoftRem(): oSoftRem();
  }
  else if(iHardRem() || oHardRem()){
    irem = iHardRem();
    rem = irem? iHardRem(): oHardRem();
  }
  else{
    return;
  }
  tParPtr part = irem? oPart(): iPart();
  if(! part->isG()){
    return;
  }
  Energy2 S = sdip();
  double C = irem ? S/(4.0*Constants::pi*(S + next()->sdip())) :
    S/(4.0*Constants::pi*(S + prev()->sdip()));

  // We will pretend parton 3 is the gluon even if it isn't. Hence
  // parton 1 may have a mass or not.
  Energy W = sqrt(S);
  pt2max = min(pt2max, S/4.0);
  if ( pt2max <= pt2min ) return;

  //Vector containing x1, x3 and the quark mass in GeV
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
      pt2m = rndsud(C*yint*CW, pt2m, pt2min);
      if ( pt2m <= pt2min ) break;
      double ymax = acosh(sqrt(S/pt2m)/2.0);
      double y = 2.0*ymax*handler()->rnd() - ymax;

      double extrat = MaxPT(S, y, rem->mu(), rem->alpha()) / sqrt(pt2m);

      double beta = handler()->beta();
      if ( beta < 0.0 &&  extrat < 1.0 ) continue;

      var[0] = 1.0 - sqrt(pt2m/S)*exp( y) - 4.0*yq;
      var[1] = 1.0 - sqrt(pt2m/S)*exp(-y);
      double xx2 = 2.0 - var[0] - var[1];

      if ( !check(var[0], var[1], 0.0, yq, yq) ) continue;

      double weight = (sqr(1.0 - var[1] + yq) + sqr(1.0 - xx2 + yq))*
	((pt2m/S)/(1.0 - var[1]))*2.0*ymax/yint;

      // Probability to have emissions above the soft cutoff.
      if ( beta >= 0.0 ) weight /= (1.0 + pow(extrat, -beta));

      // *** ATTENTION *** Add other weights here.

      weight *= reweight(this, state(), ifl, pt2m, var, FSGtoQQ());

      if ( weight < handler()->rnd() ) continue;
      
      lastPT2(pt2m);
      flav = irem ? ifl: -ifl;
      genVar = var;
      pt2min = pt2m;
      break;
    }
  }
}

tParPtr ExtendedDipole::performG() {
  tSoftRemPtr rem = iSoftRem()? iSoftRem(): oSoftRem();
  tParPtr part = iSoftRem()? oPart(): iPart();

  Lorentz5Momentum p1;
  Lorentz5Momentum p2;
  Lorentz5Momentum p3(part->momentum().mass());
  try {
    SimplePhaseSpace::CMS(p1, p2, p3, sdip(), genVar[0], genVar[1], 0.0, 0.0, 0.0);
  }
  catch ( ImpossibleKinematics ) {
    return tParPtr();
  }

  LorentzRotation R;
  R.rotateZ(genVar[2]);
  R.transform(Utilities::getBoostFromCM(make_pair(rem->momentum(),
          part->momentum())));

  p1.transform(R);
  p2.transform(R);
  p3.transform(R);

  rem->momentum() = p1;
  part->momentum() = p3;
  return emitGluon(p2);
}

void ExtendedDipole::genHRemG(Energy2 pt2min, Energy2 pt2max) {
  if(iSoftRem() || oSoftRem()){
    return;
  }
  double C = (iPart()->isG() || oPart()->isG()? 3.0/4.0: 2.0/3.0)/
    Constants::pi;
  C *= preweight(this, state(), ParticleID::g, HRGluon());
  tHardRemPtr hr = iHardRem() ? iHardRem() : oHardRem();
  C *= (hr->Q2() + hr->momentum().mt2())/hr->Q2();

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

    weight *= reweight(this, state(), ParticleID::g, pt2max, var, HRGluon());

    //Special handling of heavy quarks, the 'dead-cone' effect.
    double xm = exp(2.0*y);
    weight *= 1.0 - (y1*xm + y3/xm)/(var[0] + var[1] - 1.0);

    if ( weight < handler()->rnd() ) continue;

    lastPT2(pt2max);
    flav = ParticleID::g;
    genVar.swap(var);
    lastEType(HRGluon());
    return;
  }
}

tParPtr ExtendedDipole::performHRemG() {
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

  //if ( hr && hr == iPart() &&
       //hr->Q2()/(hr->Q2() + 4.0*p1.mt2()) > handler()->rnd() )
    //hr = tHardRemPtr();
  //if ( hr && hr == oPart() &&
       //hr->Q2()/(hr->Q2() + 4.0*p3.mt2()) > handler()->rnd() )
    //hr = tHardRemPtr();
  if ( ( iHardRem() && iHardRem()->Q2() /
	 (iHardRem()->Q2() + p1.mt2()) > handler()->rnd() ) ||
       ( oHardRem() && oHardRem()->Q2() /
	 (oHardRem()->Q2() + p3.mt2()) > handler()->rnd() ) ) {
    return tParPtr();
  }

  iPart()->momentum() = p1;
  oPart()->momentum() = p3;
  return emitGluon(p2);

}

void ExtendedDipole::genIniG(Energy2 pt2min, Energy2 pt2max) {
  if(iSoftRem() && oSoftRem()){
    return;
  }
  if( !(iSoftRem() || oSoftRem()) ){
    return;
  }

  tSoftRemPtr rem = iSoftRem()? iSoftRem(): oSoftRem();

  double C = 2.0/3.0/Constants::pi;
  C *= preweight(this, state(), ParticleID::g, RFGluon());
  Energy2 S = (rem->momentum() + state()->hardSubSys().momentum()).m2();
  Energy W = sqrt(S);
  Energy2 mh2 = state()->hardSubSys().momentum().mass2();
  Energy2 Q2 = state()->hardSubSys().Q2();
  //-------------------------------------------------------------
    //C *= 1 + rem->momentum().perp2() / sqr(rem->mu());
    C *= 1 + rem->momentum().perp2() / (0.36*GeV2);
  //-------------------------------------------------------------
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

    //Calculations below are done assuming that the remnant is in
    //positive y direction.

    //double extrat = sqr(rem->mu())/pt2max;
    double extrat = MaxPT(S, var[0], rem->mu(), rem->alpha()) /
      sqrt(pt2max);

    double beta = handler()->beta();
    if ( beta < 0.0 &&  extrat < 1.0 ) continue;

    Energy pt = sqrt(pt2max);

    //Check so that the gluon has a higher p+ than the hard sub sys
    if ( pt*exp(var[0]) < mh2/W ) continue;
    
    //Check kinematics
    if ( pt*exp(-var[0]) >= W ) continue;

    //S of the hard system after the gluon is emitted.
    Energy pphs = (mh2 + pt2max)/(W - pt*exp(-var[0]));
    Energy2 sh = W*(pphs + pt*exp(var[0]));
    double z = (mh2 + Q2)/(sh + Q2);

    //Check if there is enough phase space availible.
    if( sh >= S ) continue;

    //Check so that en inviraint mass of the gluon and the remnant is
    //larger than 4*mq^2, where mq is the mass of the extracted quark.
    if( (W - pphs) * pt * exp(-var[0]) - pt2max <
      4 * sqr(rem->data().generateMass()) ) {
      continue;
    }

    double weight = ( 1.0 + sqr(z) )/2.0;
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

tParPtr ExtendedDipole::performIniG() {
  tSoftRemPtr rem = iSoftRem()? iSoftRem(): oSoftRem();
  Energy pt = sqrt(lastPT2());
  double y = genVar[0];
  Lorentz5Momentum ph = state()->hardSubSys().momentum();
  Energy mh = ph.mass();
  Lorentz5Momentum pr = rem->momentum();
  Energy2 S = (pr+ph).m2();
  Energy W = sqrt(S);
  //LorentzRotation rcm = Utilities::getBoostFromCM(make_pair(pr,ph));

  double phi = handler()->rnd(2.0*Constants::pi);
  Lorentz5Momentum pg = ptRapidity(0.0*GeV, pt, y, phi);
  pg.setMass(0.0*GeV);
  //ph = lightCone((sqr(mh) + lastPT2()) / (sqrt(S) - pg.minus()), 
      //sqrt(S) - pg.minus(), pt*cos(phi + Constants::pi), 
      //pt*sin(phi + Constants::pi));
  //ph.setMass(mh);
  //pr = lightCone(sqrt(S) - pg.plus() - ph.plus(), 0.0*GeV, 0.0*GeV, 0.0*GeV);
  //pr.setMass(0.0*GeV);
  //-------------------------------------------------------------
  LorentzRotation rcm;// = Utilities::getBoostFromCM(make_pair(pr,ph));
  rcm.boostZ((pr+ph).boostVector().z());
  Transverse <Energy> tpt( pt * cos(phi), pt * sin(phi) );
  /*double aa = pow(sqr(rem->mu())/lastPT2(), rem->alpha()/2);
  double amin = pt*exp(y)/W;
  double beta = handler()->beta();
  double r = sqrt( amin / (aa * pow( (1 + pow(amin/aa, beta)) / 
            handler()->rnd() - 1.0 , 1.0/beta ) ) );*/
  /*double r = 1.0;
  double ymax = acosh(sqrt(S/lastPT2())/2.0);
  if(lastPT2() > 4*pr.perp2()){
    r = (1 + y/ymax) / 2;
  }*/
  /*Energy2 pt2cut = sqr(handler()->pTCut());
  Energy2 Q2 = state()->hardSubSys().Q2();
  if( y < log(sqr(rem->mu())/pt2cut) - log(Q2/pt2cut) &&
      lastPT2() > pr.perp2() ){
    r = 0.0;
  }*/
  double r = 0.0;
  Transverse <Energy> remPt(pr.x() - r*tpt.x(), pr.y() - r*tpt.y());
  Transverse <Energy> hPt(ph.x() - (1-r)*tpt.x(), ph.y() -
      (1-r)*tpt.y());
  Energy2 mthr = (W - pg.minus()) * (W - pg.plus());
  Energy pz;
  try{
    pz = SimplePhaseSpace::getMagnitude(mthr,
        sqrt(sqr(mh) + hPt.pt2()), remPt.pt());
  }
  catch ( ImpossibleKinematics ) {
    return tParPtr();
  }
  Energy er = sqrt( sqr(pz) + remPt.pt2() );
  pr = lightCone(er + pz, er - pz, remPt);
  pr.setMass(0.0*GeV);
  pr.boost(0.0, 0.0, -pg.z()/sqrt(mthr + sqr(pg.z())));
  ph = lightCone(W - pg.plus() - pr.plus(), W - pg.minus() - pr.minus(), 
      hPt);
  ph.setMass(mh);
  //-------------------------------------------------------------

  ph.transform(rcm);
  pr.transform(rcm);
  pg.transform(rcm);

  state()->hardSubSys().setMomentum(ph, pr.z() > 0.0*GeV);
  rem->momentum() = pr;

  tParPtr emitted = emitGluon(pg);

  if ( !emitted ) return emitted;

  //-------------------------------------------------------------
  if(!(iHardRem() || oHardRem()) && 
      0.36*GeV2 / (rem->momentum().perp2() + 0.36*GeV2) <
        handler()->rnd()){
      /*sqr(rem->mu()) / (rem->momentum().perp2() + sqr(rem->mu())) <
        handler()->rnd()){*/
    return tParPtr();
  }
  //-------------------------------------------------------------
  return emitted;
}

void ExtendedDipole::genRemG(Energy2 pt2min, Energy2 pt2max) {
  if(!iSoftRem() || !oSoftRem()){
    return;
  }

  double C = (iPart()->isG() || oPart()->isG()? 3.0/4.0: 2.0/3.0)/Constants::pi;
  C *= preweight(this, state(), ParticleID::g, RRGluon());
  int n1 =  iPart()->isG()? 3: 2;
  int n3 =  oPart()->isG()? 3: 2;
  Energy2 S = sdip();
  pt2max = min(pt2max, S/4.0);
  if ( pt2max <= pt2min ) return;
  double yint = 2.0*acosh(sqrt(S/pt2min)/2.0);
  double yh = state()->hardSubSys().momentum().rapidity() - 
    (iSoftRem()->momentum() + oSoftRem()->momentum()).rapidity();
  Energy mh = state()->hardSubSys().momentum().mass();

  //Vector containing rapidity.
  vector<double> var(1);

  while (true) {
    pt2max = rndsud(2.0*C*yint, pt2max, pt2min);
    if ( pt2max <= pt2min ) return;
    double ymax = acosh(sqrt(S/pt2max)/2.0);
    var[0] = 2.0*ymax*handler()->rnd() - ymax;


    double extrat = MaxPT(S, var[0], iSoftRem()->mu(), 
        oSoftRem()->mu(), iSoftRem()->alpha(), oSoftRem()->alpha()) /
      sqrt(pt2max);

    double beta = handler()->beta();
    if ( beta < 0.0 &&  extrat < 1.0 ) continue;

    double x1 = 1.0 - sqrt(pt2max/S)*exp( var[0]);
    double x3 = 1.0 - sqrt(pt2max/S)*exp(-var[0]);

    //Check if the emission is kinimatically possible.
    //yh is the rapidity of the hard sub sys in the dipole rest frame.
    //Check if pplus and pminus for the hard sub sys and the gluon
    //is larger than sqrt(S).
    if( (sqrt(pt2max+sqr(mh))-mh)*exp(yh) + sqrt(pt2max)*exp(var[0]) + 
        10.0*GeV > sqrt(S) ||
        (sqrt(pt2max+sqr(mh))-mh)*exp(-yh) + sqrt(pt2max)*exp(-var[0]) +
        10.0*GeV > sqrt(S) ){
      continue;
    }
    //Check so that the invariant mass of the gluon and any of the
    //hardon remnants is larger than 4*mq^2, where mq is the mass of the
    //extracted quark.
    int posdir = oSoftRem()->momentum().z() > 0.0*GeV ? 1 : -1;
    if( ( sqrt(S) + (mh - sqrt(sqr(mh) + pt2max)) * exp(posdir * yh) ) *
        sqrt(pt2max) * exp(-posdir * var[0]) - pt2max < 
        4 * sqr(oSoftRem()->data().generateMass()) ||
        ( sqrt(S) + (mh - sqrt(sqr(mh) + pt2max)) * exp(-posdir * yh) ) *
        sqrt(pt2max) * exp(posdir * var[0]) - pt2max < 
        4 * sqr(iSoftRem()->data().generateMass()) ){
      continue;
    }

    double weight = (pow(x1, n1) + pow(x3, n3))/2.0;
    weight *= 2.0*ymax/yint;

    // Probability to have emissions above the soft cutoff.
    if ( beta >= 0.0 ) {
      weight /= (1.0 + pow(extrat, -beta));
    }

    weight *= reweight(this, state(), ParticleID::g, pt2max, var, RRGluon());

    if ( weight < handler()->rnd() ) continue;

    lastPT2(pt2max);
    flav = ParticleID::g;
    genVar.swap(var);
    lastEType(RRGluon());
    return;
  }
}

tParPtr ExtendedDipole::performRemG() {
  if(!iSoftRem() || !oSoftRem()){
    return tParPtr();
  }

  Energy pt = sqrt(lastPT2());
  double yg = genVar[0] + (iSoftRem()->momentum() +
      oSoftRem()->momentum()).rapidity();
  Lorentz5Momentum ph = state()->hardSubSys().momentum();
  double yh = ph.rapidity();

  double phi = handler()->rnd(2.0*Constants::pi);
  Lorentz5Momentum pg = ptRapidity(0.0*GeV, pt, yg, phi);
  pg.setMass(0.0*GeV);

  LorentzMomentum phn = ptRapidity(ph.mass(), pt, yh, phi + Constants::pi);
  state()->hardSubSys().setMomentum(phn, yg > yh);

  //Calculate change in pplus and pminus
  Energy dpp = pg.plus() + phn.plus() - ph.plus();
  Energy dpm = pg.minus() + phn.minus() - ph.minus();

  Lorentz5Momentum &pi = iSoftRem()->momentum();
  Lorentz5Momentum &po = oSoftRem()->momentum();
  if(pi.z() > 0.0*GeV){
    pi = lightCone(pi.plus() - dpp, pi.minus(), pi.x(), pi.y());
    po = lightCone(po.plus(), po.minus() - dpm, po.x(), po.y());
  }
  else{
    pi = lightCone(pi.plus(), pi.minus() - dpm, pi.x(), pi.y());
    po = lightCone(po.plus() - dpp, po.minus(), po.x(), po.y());
  }

  return emitGluon(pg);
  
}

ParPtr ExtendedDipole::emitGluon(const Lorentz5Momentum & pg) {
  ParPtr g =
    state()->create(handler()->generator()->getParticleData(ParticleID::g));
  if ( !g ) return g;
  g->momentum() = pg;
  g->parents(make_pair(iPart(), oPart()));
  g->string(iPart()->string());

  iPart()->touch();
  oPart()->touch();

  tDipPtr d;
  if(iSoftRem() || iHardRem()){
    if(oSoftRem() || oHardRem()){
      d = state()->create<ExtendedDipole>();
    }
    else{
      d = state()->create<Dipole>();
    }
    oPart()->iDip(d);
    d->init(g, oPart());
    g->oDip(d);
    g->iDip(this);
    oPart(g);
    oSoftRem(tSoftRemPtr());
  }
  else{
    d = state()->create<Dipole>();
    iPart()->oDip(d);
    d->init(iPart(), g);
    g->iDip(d);
    g->oDip(this);
    iPart(g);
    iSoftRem(tSoftRemPtr());
  }

  d->generateColourIndex();

  touch();

  return g;
}

void ExtendedDipole::genSeaQ(Energy2 pt2min, Energy2 pt2max){
  if ( pt2max <= pt2min ){
    return;
  }

  if(iSoftRem()){
    genSeaQ(iSoftRem(), pt2min, pt2max);
  }

  if(oSoftRem()){
    genSeaQ(oSoftRem(), max(pt2min, lastPT2()), pt2max);
  }
}

void ExtendedDipole::genSeaQ(tSoftRemPtr par,Energy2 pt2min, Energy2 pt2max){
  if( par->isG() || !par->data().coloured() ){
    return;
  }

  HardSubSys & hardsys = state()->hardSubSys();
  Energy2 mh2 = hardsys.momentum().mass2();
  double x = par->partonx();

  int ifl = -par->data().id();
  tcPDPtr q = handler()->generator()->getParticleData(ifl);
  Energy mq = q->generateMass();
  Energy2 mq2 = sqr(mq);

  //Calculate sdip and maximum possible pt2 assuming mq = 0
  Energy2 s = (par->momentum() + hardsys.momentum()).m2();
  Energy2 mt2km = sqr(s - mh2)/(4.0 * s) + mq2;
  //Lower pt2max to the maximum possible kinimatically
  Energy2 mt2max = min(pt2max, mt2km);
  Energy2 mt2min = max(pt2min, sqr(handler()->pTCut()) + mq2);
  if(mt2max <= mt2min){
    return;
  }

  //Use the geometric average of max and min as scale for the pdf.
  //If pt2max is smaller use that value.
  Energy2 mt2z = min( mt2max, sqrt( mt2km * mt2min ) );
  Energy2 mt2mean = sqrt( mt2z * mt2min );
  //Calculate zmax for pt2z and pt2min
  double zmax = mh2 / sqr( sqrt(mt2min - mq2 + mh2) + sqrt(mt2min) );
  double zmax1 = mh2 / sqr( sqrt(mt2mean - mq2 + mh2) + sqrt(mt2mean) );
  double zmax2 = mh2 / sqr( sqrt(mt2z - mq2 + mh2) + sqrt(mt2z) );
  tcPDPtr g = handler()->generator()->getParticleData(ParticleID::g);

  //multiply pdf by a safty factor of two
  //use zero if zmax is less than x
  double maxpdf = 2.0 * max(
      PDFRatio(par, g, par->dataPtr(), mt2z, x, zmax2),
      max( PDFRatio(par, g, par->dataPtr(), mt2mean, x, zmax1),
      PDFRatio(par, g, par->dataPtr(), mt2min, x, zmax) ) );
  if( maxpdf <= 0.0){
    return;
  }

  //xi is the fraction of the positive light cone momentum taken by the
  //hard sub system. x < xi < 1
  double C = 1.0/(4.0*Constants::pi) * (1.0 - x);
  C *= preweight(this, state(), ifl, ISGtoQQ());

  //Vector containing xi, z and the quark mass in GeV.
  vector<double> var(3);
  var[2] = mq/GeV;

  while ( true ) {
    mt2max = rndsud(C * maxpdf, mt2max, mt2min);
    if ( mt2max <= mt2min ) break;

    var[0] = x + (1.0 - x) * handler()->rnd();

    //check if xi is kinematically possible
    if( (mt2max + mh2 - mq2)/var[0] + mt2max/(1-var[0]) > mh2/x ){
      continue;
    }
    var[1] = mh2 * var[0] * (1-var[0]) /
      (mh2*(1-var[0]) + mq2*(var[0]-1) + mt2max);

    double weight = (sqr(var[1]) + sqr(1 - var[1])) / maxpdf * 
      PDFRatio(par, g, par->dataPtr(), mt2max+mq2, x, var[1]);

    //multiply by jacobi determinant
    weight *= mt2max / ( (var[0]-1) * mq2 + mt2max ) * var[1] / var[0];

    weight *=  reweight(this, state(), ifl, mt2max - mq2, var, ISGtoQQ());

    if ( weight > 1.0 ) handler()->generator()->logWarning(
        WeightException() << "Ariadne::ExtendedDipole::genSeaQ "
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

void ExtendedDipole::genQtoGQ(Energy2 pt2min, Energy2 pt2max){
  if ( pt2max <= pt2min ){
    return;
  }

  if(iSoftRem()){
    genQtoGQ(iSoftRem(), pt2min, pt2max);
  }

  if(oSoftRem()){
    genQtoGQ(oSoftRem(), max(pt2min, lastPT2()), pt2max);
  }
}

void ExtendedDipole::genQtoGQ(tSoftRemPtr par,Energy2 pt2min, Energy2 pt2max){
  if( !par->isG() ){
    return;
  }

  HardSubSys & hardsys = state()->hardSubSys();
  Energy2 mh2 = hardsys.momentum().mass2();
  double x = par->partonx();

  //Calculate sdip and maximum possible pt2 assuming mq = 0
  Energy2 s = (par->momentum() + hardsys.momentum()).m2();
  Energy2 pt2km = sqr(s - mh2)/(4.0 * s);
  //Lower pt2max to the maximum possible kinimatically
  if(pt2max > pt2km){
    pt2max = pt2km;
  }
  if(pt2max <= pt2min){
    return;
  }

  //Use the geometric average of max and min as scale for the pdf.
  //If pt2max is smaller use that value.
  Energy2 pt2z = min( pt2max, sqrt( pt2km * pt2min ) );
  //Calculate zmax for pt2z and pt2min
  double zmax = mh2 / sqr( sqrt(pt2min + mh2) + sqrt(pt2min) );
  double zmax2 = mh2 / sqr( sqrt(pt2z + mh2) + sqrt(pt2z) );
  tcPDPtr g = handler()->generator()->getParticleData(ParticleID::g);

  //Vector containing xi, z and the quark mass in GeV.
  vector<double> var(3);
  //xi is the fraction of the positive light cone momentum taken by the
  //hard sub system. x < xi < 1

  for ( int ifl = 1; ifl <= handler()->nFlav(); ++ifl ) {
    tcPDPtr q;
    if(par == iSoftRem()){
      q = handler()->generator()->getParticleData(-ifl);
    }
    else{
      q = handler()->generator()->getParticleData(ifl);
    }
    if(!q){
      continue;
    }

    Energy mq = q->generateMass();
    Energy2 mq2 = sqr(mq);
    var[2] = mq/GeV;

    Energy2 mt2m = min(pt2max, pt2km + mq2);
    Energy2 mt2min = max(pt2min, sqr(handler()->pTCut()) + mq2);

    //multiply pdf by a safty factor of two
    //use zero if zmax is less than x
    //sample the pdf at 0.1 as well to avoid problems
    //with valance quarks.
    double maxpdf = 2.0 * max(
        PDFRatio(par, q, g, pt2z+mq2, x, zmax2),
        max(PDFRatio(par, q, g, mt2min, x, zmax),
          (x > 0.1 || x/zmax > 0.1) ? 0.0 :
          PDFRatio(par, q, g, mt2min, x, x*10.0) ));
    if( maxpdf <= 0.0){
      continue;
    }

    double C = 4.0/(3.0*Constants::pi) * log(1.0 / x);
    C *= preweight(this, state(), ifl, ISQtoGQ());

    while ( true ) {
      mt2m = rndsud(C*maxpdf, mt2m, mt2min);
      if ( mt2m <= mt2min ) break;

      //generate xi according to dxi/xi
      var[0] = pow(x, handler()->rnd());

      //check if xi is kinematically possible
      if( (mt2m + mh2 - mq2)/var[0] + mt2m/(1-var[0]) > mh2/x ){
        continue;
      }
      var[1] = mh2 * var[0] * (1-var[0]) /
        (mh2*(1-var[0]) + mq2*(var[0]-1) + mt2m);

      double weight = (1.0 + sqr(1.0 - var[1])) / 2.0;
      weight *= PDFRatio(par, q, g, mt2m, x, var[0])/maxpdf;

      //multiply by jacobi determinant * z
      weight *= mt2m / ( mq2 * (var[0] - 1) + mt2m );

        weight *= reweight(this, state(), ifl, mt2m - mq2, var, ISQtoGQ());

      if ( weight > 1.0 ) handler()->generator()->logWarning(
          WeightException() << "Ariadne::ExtendedDipole::genQtoGQ "
          << "failed to overestimate the PDF ratio. "
          << "If this hapens too often you should contact the authors."
          << Exception::warning);			   

      if ( weight < handler()->rnd() ) continue;

      lastPT2(mt2m);
      flav = q->id();
      genVar = var;
      mt2min = mt2m;
      //set the emission type
      lastEType(ISQtoGQ());
      break;
    }
  }
}

tParPtr ExtendedDipole::performSeaQ(){
  tParPtr q = performInitialQ();
  if ( !q ) return q;

  tSoftRemPtr rem;
  if(flav < 0){
    rem = iSoftRem();
  }
  else{
    rem = oSoftRem();
  }

  rem->data(handler()->generator()->getParticleData(ParticleID::g));
  q->string( rem->string() );
  q->string()->moveEndpoint(rem, q);

  ExDipPtr d = state()->create<ExtendedDipole>();
  if(flav < 0){
    d->init(q, rem);
    q->oDip(d);
    rem->iDip(d);
  }
  else{
    d->init(rem, q);
    q->iDip(d);
    rem->oDip(d);
  }
  d->touch();

  touch();

  return q;
}

tParPtr ExtendedDipole::performQtoGQ(){
  tParPtr q = performInitialQ();
  if ( !q ) return q;

  //Change the flavor and split the sting.
  //Also replace the exteded dipole with an ordinary dipole
  //if appropriate.
  if(flav < 0){
    iSoftRem()->data(handler()->generator()->getParticleData(flav));
    iSoftRem()->string()->split(iSoftRem(), iSoftRem(), q);
    iSoftRem()->oDip(tDipPtr());
    if(oSoftRem()){
      iSoftRem(tSoftRemPtr());
    }
    else{
      DipPtr d = state()->create<Dipole>();
      d->init(q, oPart());
      q->oDip(d);
      oPart()->iDip(d);
      d->touch();
      state()->removeEmitter(this);
    }
  }
  else{
    oSoftRem()->data(handler()->generator()->getParticleData(flav));
    oSoftRem()->string()->split(oSoftRem(), q, oSoftRem());
    oSoftRem()->iDip(tDipPtr());
    if(iSoftRem()){
      oSoftRem(tSoftRemPtr());
    }
    else{
      DipPtr d = state()->create<Dipole>();
      d->init(iPart(), q);
      q->iDip(d);
      iPart()->oDip(d);
      d->touch();
      state()->removeEmitter(this);
    }
  }

  touch();
  return q;
}

tParPtr ExtendedDipole::performInitialQ(){
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

  double phi = 2.0 * Constants::pi * handler()->rnd();
  Energy mq = genVar[2]*GeV;
  Energy pt = sqrt(lastPT2() - sqr(mq));
  Transverse<Energy> tpt( pt * cos(phi), pt * sin(phi) );
  double xi = genVar[0];
  double z = genVar[1];

  Energy dpp = ph.plus() * (1.0 / z - 1.0);
  if(dpp > pr.plus()){
    return tParPtr();
  }
  pr = lightCone(pr.plus() - dpp, pr.minus(), pr.x(), pr.y());

  //pplus of q plus hard sub sys.
  Energy pplus = ph.plus()/z;
  Energy pminus = ph.minus();

  Lorentz5Momentum pq = lightCone((1-xi) * pplus,
      lastPT2()/((1-xi)*pplus), tpt);
  ph = lightCone(pplus - pq.plus(), pminus - pq.minus(), -tpt);
  pq.setMass(mq);

  //boost to lab frame
  rcm.invert();
  ph.transform(rcm);
  pq.transform(rcm);
  pr.transform(rcm);

  state()->hardSubSys().setMomentum(ph, pr.z() > 0.0*GeV);

  rem->momentum() = pr;
  rem->touch();

  //create quark
  ParPtr q = state()->create(handler()->generator()->getParticleData(flav));
  q->momentum() = pq;
  q->parents( make_pair(iPart(), oPart()) );
  q->touch();

  return q;
}

double ExtendedDipole::PDFRatio(tSoftRemPtr rem, tcPDPtr f1, tcPDPtr f2, 
    Energy2 pt2, double x, double z){
  if(x <= 0.0 || x >= 1.0 || z <= x || z >= 1.0){
    return 0.0;
  }

  double nom = rem->pdf().xfx(f1, pt2, x/z); 
  double denom = rem->pdf().xfx(f2, pt2, x);
  if(denom <= 0.0){
    return 0.0;
  }

  return nom/denom;
}

Energy ExtendedDipole::MaxPT(Energy2 S, double y, Energy mu1,
    Energy mu3, double alpha1, double alpha3){
  Energy W = sqrt(S);
  Energy pt = W*min( pow(pow(mu1/W, alpha1)/exp(y), 1/(1+alpha1)),
		      pow(pow(mu3/W, alpha3)/exp(-y), 1/(1+alpha3)));
  Energy ptm;
  do{
    ptm = pt;
    double b1 = exp(y) * pow(pt/mu1, alpha1);
    double b3 = exp(-y) * pow(pt/mu3, alpha3);
    pt = pt * W * ( (1 + alpha1) * b1 + (1 + alpha3) * b3 ) / 
      ( pt * sqr(b1 + b3) +  W * ( alpha1 * b1 + alpha3 * b3) );
  }
  while(abs(pt/ptm-1) > 1e-3);
  return pt;
}

Energy ExtendedDipole::MaxPT(Energy2 S, double y, Energy mu,
    double alpha){
  Energy W = sqrt(S);
  Energy pt = W*min( 1.0/exp(-y),
		     pow(pow(mu/W, alpha)/exp(y), 1/(1+alpha)));
  Energy ptm;
  do{
    ptm = pt;
    double b = exp(y) * pow(pt/mu, alpha);
    pt = pt * W * (exp(-y) + (1+alpha) * b) /
      (pt * sqr( exp(-y) + b) + alpha * b * W);
  }
  while(abs(pt/ptm-1) > 1e-3);
  return pt;
}

double ExtendedDipole::emissionProbability(){
  Energy2 S = sdip();
  double scale = sqr(handler()->pTCut()) / lastPT2();
  tSoftRemPtr rem = iSoftRem()? iSoftRem(): oSoftRem();
  tParPtr part = iSoftRem()? oPart(): iPart();
  double rw = preweight(this, state(), flav, *lastEType()) *
    reweight(this, state(), flav, lastPT2(), genVar, *lastEType());



  if ( isType(FSGluon()) ) {     
    double C = (part->isG()? 3.0/4.0: 2.0/3.0)/
      Constants::pi;
    int n1 =  rem->isG()? 3: 2;
    int n3 =  part->isG()? 3: 2;
    double y3 = max(part->momentum().mass2()/S, 0.0);
    //dead cone
    double xm = (1 - genVar[0] - y3) / (1 - genVar[1] + y3);
    double dcw = max(0.0, 1.0 - (y3/xm)/(genVar[0] + genVar[1] - 1.0));
    return C * dcw * (pow(genVar[0], n1) + pow(genVar[1], n3)) * scale * rw;
  }
  else if ( isType(FSGtoQQ()) ) {
    double C = flav > 0 ? S/(4.0*Constants::pi*(S + next()->sdip())):
      S/(4.0*Constants::pi*(S + prev()->sdip()));
    double yq = genVar[2]*GeV2 / S;
    double x2 = 2.0 - genVar[0] - genVar[1];
    return C*(sqr(1.0 - genVar[1] + yq) + sqr(1.0 - x2 + yq))/
      (1.0 - genVar[0]) * lastPT2()/S * scale * rw;
  }
  else if ( isType(RRGluon()) ) {
    double C = (iPart()->isG() || oPart()->isG()? 3.0/4.0: 2.0/3.0)/Constants::pi;
    int n1 =  iPart()->isG()? 3: 2;
    int n3 =  oPart()->isG()? 3: 2;
    double xt2 = lastPT2() / S;
    double x1 = 1.0 - sqrt(xt2)*exp( genVar[0]);
    double x3 = 1.0 - sqrt(xt2)*exp(-genVar[0]);
    return C*(pow(x1, n1) + pow(x3, n3)) * scale * rw;
  }
  else if ( isType(ISGtoQQ()) ) {
    double jacobi = lastPT2()/(genVar[0]*sqr(genVar[2]*GeV) + lastPT2())
      *genVar[1]/genVar[0];
    return 1.0/(4.0*Constants::pi)*(sqr(genVar[1]) + sqr(1 - genVar[1]))*jacobi*scale*rw;
  }
  else if ( isType(ISQtoGQ()) ) {
    double jacobi = lastPT2()/(genVar[0]*sqr(genVar[2]*GeV) + lastPT2())
      *genVar[1]/genVar[0];
    return 2.0/(3.0*Constants::pi)*(1.0 + sqr(1.0 - genVar[1]))/genVar[1]*jacobi*scale*rw;
  }
  else if ( isType(RFGluon()) ) {
    double C = 2.0/3.0/Constants::pi;
    Energy2 mh2 = state()->hardSubSys().momentum().mass2();
    Energy W = sqrt(S);
    Energy pphs = (mh2 + lastPT2())/
      (W - sqrt(lastPT2())*exp(-genVar[0]));
    Energy2 sh = W*(pphs + sqrt(lastPT2())*exp(genVar[0]));
    Energy2 Q2 = state()->hardSubSys().Q2();
    double z = (mh2 + Q2)/(sh + Q2);
    return C * (1.0 + sqr(z)) * scale * rw;
  }
  else {
    return 0.0;
  }
}

double ExtendedDipole::PDFRatio(){
  if ( isType(ISGtoQQ()) || isType(ISQtoGQ()) ){
    tSoftRemPtr rem;
    if(flav < 0){
      rem = iSoftRem();
    }
    else{
      rem = oSoftRem();
    }
    tcPDPtr q = handler()->generator()->getParticleData(flav);
    tcPDPtr g = handler()->generator()->getParticleData(ParticleID::g);
    if ( isType(ISGtoQQ()) )
      return PDFRatio(rem, g, q, lastPT2(), genVar[1], genVar[0]);
    else
      return PDFRatio(rem, q, g, lastPT2(), genVar[1], genVar[0]);
  }

  double extrat = 0.0;
  tSoftRemPtr rem = iSoftRem()? iSoftRem(): oSoftRem();
  if ( isType(FSGluon()) ) {
    return 1.0;
  }
  else if ( isType(FSGtoQQ()) ) {
    double yq = sqr(genVar[2]*GeV)/sdip();
    double aa = pow(sqr(rem->mu())/lastPT2(), rem->alpha());
    double expy = sqrt((1.0 - genVar[0] - 4*yq)/(1.0 - genVar[1]));
    extrat = sdip()*sqr(1.0/expy + expy/aa)/lastPT2();
  }
  else if ( isType(RRGluon()) ) {
    double aa1 = pow(sqr(iSoftRem()->mu())/lastPT2(),
		     iSoftRem()->alpha()/2);
    double aa3 = pow(sqr(oSoftRem()->mu())/lastPT2(),
		     oSoftRem()->alpha()/2);
    extrat = (sdip()/sqr(exp(-genVar[0])/aa3 +
			 exp(genVar[0])/aa1)) / lastPT2();
  }
  else if ( isType(RFGluon()) ) {
    double aa = pow(sqr(rem->mu())/lastPT2(), rem->alpha()/2);
    extrat = (sdip()/sqr(exp(-genVar[0]) +
			 exp(genVar[0])/aa)) / lastPT2();
  }
  else {
    return 0.0;
  }

  double beta = handler()->beta();
  if ( beta < 0.0 &&  extrat < 1.0 ){
    return extrat < 1.0 ? 0.0 : 1.0;
  }
  else{
    return 1.0 / (1.0 + pow(extrat, -beta / 2.0));
  }
}


Emitter::DipoleStateVector ExtendedDipole::constructStep(){
  DipoleStateVector vec;

  //undo gluon
  if( ( iSoftRem() && !oSoftRem() && oPart()->isG() ) ||
      ( oSoftRem() && !iSoftRem() && iPart()->isG() ) ) {
    //undo a gluon emission from a remnant-remnant dipole
    if( iSoftRem() && dynamic_ptr_cast<tSoftRemPtr>(next()->oPart()) ){
      TranslationMap trans;
      DipoleStatePtr dipstate = state()->fullclone(trans);
      tExDipPtr newdip = trans.translate(tExDipPtr(this));
      if(newdip->undoRemOG()){
        vec.push_back(dipstate);
      }
    }
    //undo a gluon emission from a non remnant-remnant dipole
    else if( iSoftRem() || !dynamic_ptr_cast<tSoftRemPtr>(prev()->iPart()) ){
      TranslationMap trans1;
      DipoleStatePtr dipstate1 = state()->fullclone(trans1);
      tExDipPtr newdip1 = trans1.translate(tExDipPtr(this));
      if(newdip1->undoIniG()){
        vec.push_back(dipstate1);
      }
      TranslationMap trans2;
      DipoleStatePtr dipstate2 = state()->fullclone(trans2);
      tExDipPtr newdip2 = trans2.translate(tExDipPtr(this));
      if(newdip2->undoG()){
        vec.push_back(dipstate2);
      }
    }
  }

  //undo initial state sea quark
  if( ( iSoftRem() && iSoftRem()->isG() && !oSoftRem() && !oPart()->isG() ) ||
      ( oSoftRem() && oSoftRem()->isG() && !iSoftRem() && !iPart()->isG() )){
    TranslationMap trans;
    DipoleStatePtr dipstate = state()->fullclone(trans);
    tExDipPtr newdip = trans.translate(tExDipPtr(this));
    if(newdip->undoSeaQ()){
      vec.push_back(dipstate);
    }
  }

  const HardSubSys::PartonSet partset = state()->hardSubSys().active();

  for(HardSubSys::PartonSet::const_iterator it = partset.begin();
      it != partset.end(); ++it ){
    if((*it)->isG()){
      continue;
    }

    //undo initial state QtoGQ
    if( ( iSoftRem() && iSoftRem()->data().id() == (*it)->data().id() && 
	  iSoftRem()->string() != (*it)->string() ) ||
        ( oSoftRem() && oSoftRem()->data().id() == (*it)->data().id() && 
	  oSoftRem()->string() != (*it)->string() ) ) {
      TranslationMap trans;
      DipoleStatePtr dipstate = state()->fullclone(trans);
      tExDipPtr newdip = trans.translate(tExDipPtr(this));
      tParPtr q = trans.translate(*it);
      if(newdip->undoQtoGQ(q)){
        vec.push_back(dipstate);
      }
    }

    //undo final state g->qqbar
    if( ( !iSoftRem() && iPart()->data().id() == -(*it)->data().id() && 
	  iPart()->string() != (*it)->string() ) ||
        ( !oSoftRem() && oPart()->data().id() == -(*it)->data().id() && 
	  oPart()->string() != (*it)->string() ) ) {
      TranslationMap trans;
      DipoleStatePtr dipstate = state()->fullclone(trans);
      tExDipPtr newdip = trans.translate(tExDipPtr(this));
      tParPtr q = trans.translate(*it);
      dipstate->selected(newdip);
      if(newdip->undoQ(q)){
        vec.push_back(dipstate);
      }
    }
  }
  return vec;
}

bool ExtendedDipole::undoRemOG(){
  genVar.clear();
  Lorentz5Momentum &pg = oPart()->momentum();
  if ( !removeGluon(oPart()) ) {
    return false;
  }
  Lorentz5Momentum ph = state()->hardSubSys().momentum();
  LorentzMomentum phn = ptRapidity(ph.mass(), 0.0*GeV, ph.rapidity(), 0.0);

  state()->hardSubSys().setMomentum(phn, pg.rapidity() > ph.rapidity());

  Energy dpp = pg.plus() + ph.plus() - phn.plus();
  Energy dpm = pg.minus() + ph.minus() - phn.minus();

  Lorentz5Momentum & pi = iSoftRem()->momentum();
  Lorentz5Momentum & po = oSoftRem()->momentum();
  if ( pi.z() > 0.0*GeV ) {
    pi = lightCone(pi.plus() + dpp, pi.minus(), pi.x(), pi.y());
    po = lightCone(po.plus(), po.minus() + dpm, po.x(), po.y());
  }
  else{
    pi = lightCone(pi.plus(), pi.minus() + dpm, pi.x(), pi.y());
    po = lightCone(po.plus() + dpp, po.minus(), po.x(), po.y());
  }

  flav = ParticleID::g;
  genVar.push_back(pg.rapidity() - (pi + po).rapidity());
  lastPT2(pg.perp2());
  lastEType(RRGluon());
  state()->selected(this);
  return true;
}

bool ExtendedDipole::undoIniG(){
  genVar.clear();
  tSoftRemPtr rem = iSoftRem() ? iSoftRem() : oSoftRem();
  tParPtr g = iSoftRem() ? oPart() : iPart();

  Energy2 sh = state()->hardSubSys().momentum().mass2();
  Lorentz5Momentum pg = g->momentum();
  if(!removeGluon(g)){
    return false;
  }

  Lorentz5Momentum ph = state()->hardSubSys().momentum();
  Energy2 mh2 = ph.mass2();
  Lorentz5Momentum pr = rem->momentum();
  Energy2 S = (ph + pg + pr).m2();
  Energy W = sqrt(S);
  LorentzRotation rcm = Utilities::boostToCM(makeTriplet(&pr, &pg, &ph));

  //Save gluon rapidity and z
  genVar.push_back(pg.rapidity());
  genVar.push_back(mh2 / sh);

  //Check so that the gluon has a higher p+ than the hard sub sys
  if ( pg.plus() < mh2/W ){
    return false;
  }

  ph = lightCone(mh2 / W, W, 0.0*GeV, 0.0*GeV);
  pr = lightCone(W - ph.plus(), 0.0*GeV, 0.0*GeV, 0.0*GeV);

  rcm.invert();
  ph.transform(rcm);
  pr.transform(rcm);
  state()->hardSubSys().setMomentum(ph, pr.z() > 0.0*GeV);
  rem->momentum() = pr;

  flav = ParticleID::g;
  lastPT2(pg.perp2());
  lastEType(RFGluon());
  state()->selected(this);
  return true;
}

bool ExtendedDipole::undoG(){
  genVar.clear();
  tSoftRemPtr rem = iSoftRem()? iSoftRem(): oSoftRem();
  tParPtr g = iSoftRem()? oPart(): iPart();
  tParPtr part = iSoftRem() ? oPart()->next(): iPart()->prev();

  lastPT2(g->invPT2());
  Lorentz5Momentum pg = g->momentum();
  //Check if the gluon has smaller light cone momenta than the hard
  //sub system.
  if( pg.plus() > state()->hardSubSys().momentum().plus() ||
      pg.minus() > state()->hardSubSys().momentum().minus() ){
    return false;
  }
  if(!removeGluon(g)){
    return false;
  }

  Lorentz5Momentum pp = part->momentum();
  Lorentz5Momentum pr = rem->momentum();
  Energy2 S = (pp + pg + pr).m2();
  Energy W = sqrt(S);
  LorentzRotation rcm = Utilities::boostToCM(makeTriplet(&pr, &pg, &pp));

  //Save x1 and x3
  genVar.push_back(2.0 * pr.e()/W);
  genVar.push_back(2.0 * pp.e()/W);

  pp = lightCone(pp.mass2() / W, W, 0.0*GeV, 0.0*GeV);
  pr = lightCone(W - pp.plus(), 0.0*GeV, 0.0*GeV, 0.0*GeV);

  rcm.invert();
  pp.transform(rcm);
  pr.transform(rcm);
  part->momentum() = pp;
  rem->momentum() = pr;

  flav = ParticleID::g;
  lastEType(FSGluon());
  state()->selected(this);
  return true;
}

bool ExtendedDipole::removeGluon(tParPtr g){
  state()->hardSubSys().remove(g);
  //If one of the dipoles is extended keep it.
  if( !dynamic_ptr_cast<tExDipPtr>(g->iDip()) && 
      dynamic_ptr_cast<tExDipPtr>(g->oDip()) ){
    g->oDip()->iPart(g->iDip()->iPart());
    g->prev()->oDip(g->oDip());
    state()->removeEmitter(g->iDip());
  }
  else{
    g->iDip()->oPart(g->oDip()->oPart());
    tExDipPtr od = dynamic_ptr_cast<tExDipPtr>(g->oDip());
    tExDipPtr id = dynamic_ptr_cast<tExDipPtr>(g->iDip());
    if(id && od){
      id->oSoftRem(od->oSoftRem());
    }
    g->next()->iDip(g->iDip());
    state()->removeEmitter(g->oDip());
  }
  return true;
}

bool ExtendedDipole::undoSeaQ(){
  tSoftRemPtr rem;
  tParPtr q;
  if(iSoftRem()){
    rem = iSoftRem();
    q = oPart();
    state()->selected(prev());
  }
  else{
    rem = oSoftRem();
    q = iPart();
    state()->selected(next());
  }

  flav = q->data().id();

  state()->removeEmitter(tExDipPtr(this));
  if(iSoftRem()){
    iSoftRem()->oDip(tDipPtr());
  }
  else{
    oSoftRem()->iDip(tDipPtr());
  }
  rem->string()->moveEndpoint(q, rem);
  rem->data(handler()->generator()->getParticleData(
        -q->data().id()));

  tExDipPtr exdip = dynamic_ptr_cast<tExDipPtr>(state()->selected());
  if(exdip){
    exdip->lastEType(ISGtoQQ());
    exdip->lastPT2(q->momentum().perp2());
    if(exdip->undoInitialQ(rem, q)){
      return true;
    }
  }
  return false;
}

bool ExtendedDipole::undoQtoGQ(tParPtr q){
  tSoftRemPtr rem;
  tDipPtr dip;
  if(iSoftRem() && iSoftRem()->data().id() == q->data().id()){
    rem = iSoftRem();
    dip = q->iDip();
  }
  else if(oSoftRem() && oSoftRem()->data().id() == q->data().id()){
    rem = oSoftRem();
    dip = q->oDip();
  }
  else{
    return false;
  }

  flav = q->data().id();
  rem->data(handler()->generator()->getParticleData(ParticleID::g));
  rem->string()->join(q, rem);

  //If the dipole that was connected to the quark was an ordinary dipole
  //replace it with an extended dipole.
  if(! dynamic_ptr_cast<tExDipPtr>(dip) ){
    tExDipPtr exdip = state()->create<ExtendedDipole>();
    exdip->init(dip->iPart(), dip->oPart());
    exdip->iPart()->oDip(exdip);
    exdip->oPart()->iDip(exdip);
    state()->removeEmitter(dip);
  }

  if(rem == iSoftRem()){
    if(flav < 0){
      state()->selected(this);
    }
    else{
      state()->selected(prev());
    }
  }
  else{
    if(flav < 0){
      state()->selected(next());
    }
    else{
      state()->selected(this);
    }
  }

  tExDipPtr exdip = dynamic_ptr_cast<tExDipPtr>(state()->selected());
  if(exdip){
    exdip->lastEType(ISQtoGQ());
    exdip->lastPT2(q->momentum().perp2());
    if(exdip->undoInitialQ(rem, q)){
      return true;
    }
  }
  return false;
}

bool ExtendedDipole::undoInitialQ(tSoftRemPtr rem, tParPtr q){
  genVar.clear();
  Energy2 sh = state()->hardSubSys().momentum().mass2();

  state()->hardSubSys().remove(q);
  Lorentz5Momentum ph = state()->hardSubSys().momentum();
  Lorentz5Momentum pq = q->momentum();
  Lorentz5Momentum pr = rem->momentum();
  LorentzRotation rcm = Utilities::boostToCM(makeTriplet(&pr, &pq, &ph));

  genVar.push_back( ph.plus()/(ph.plus() + pq.plus()) );
  genVar.push_back( ph.mass2()/sh );

  Energy W = (ph + pq + pr).m();
  ph = lightCone(ph.mass2()/W, W, 0.0*GeV, 0.0*GeV);
  pr = lightCone(W - ph.plus(), 0.0*GeV, 0.0*GeV, 0.0*GeV);

  rcm.invert();
  ph.transform(rcm);
  pr.transform(rcm);

  rem->momentum() = pr;
  state()->hardSubSys().setMomentum(ph, pr.z() > 0.0*GeV);
  genVar.push_back(pq.mass() / GeV);

  return true;
}

void ExtendedDipole::fillReferences(CloneSet & cset) const {
  Dipole::fillReferences(cset);
  cset.insert(oSoftRem());
  cset.insert(iSoftRem());
  cset.insert(oHardRem());
  cset.insert(iHardRem());
}

void ExtendedDipole::rebind(const TranslationMap & trans) {
  Dipole::rebind(trans);
  theOSoftRem = trans.translate(oSoftRem());
  theISoftRem = trans.translate(iSoftRem());
  theOHardRem = trans.translate(oHardRem());
  theIHardRem = trans.translate(iHardRem());
}

void ExtendedDipole::persistentOutput(PersistentOStream & os) const {
  os << theISoftRem << theOSoftRem << theIHardRem << theOHardRem;
}

void ExtendedDipole::persistentInput(PersistentIStream & is, int) {
  is >> theISoftRem >> theOSoftRem >> theIHardRem >> theOHardRem;
}

ClassDescription<ExtendedDipole> ExtendedDipole::initExtendedDipole;
// Definition of the static class description member.

void ExtendedDipole::Init() {}

void ExtendedDipole::debugme() const {
  Dipole::debugme();
  cerr << "ex";
}

bool ExtendedDipole::checkIntegrety(){
  return Dipole::checkIntegrety() &&
    (iSoftRem() || oSoftRem() || iHardRem() || oHardRem()) &&
    dynamic_ptr_cast<tSoftRemPtr>(iPart()) == iSoftRem() &&
    dynamic_ptr_cast<tSoftRemPtr>(oPart()) == oSoftRem() &&
    dynamic_ptr_cast<tHardRemPtr>(iPart()) == iHardRem() &&
    dynamic_ptr_cast<tHardRemPtr>(oPart()) == oHardRem();
}

