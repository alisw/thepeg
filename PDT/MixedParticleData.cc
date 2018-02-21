// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MixedParticleData class.
//

#include "MixedParticleData.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Utilities/Debug.h"
#include "ThePEG/PDT/DecayMode.h"

using namespace ThePEG;

void MixedParticleData::persistentOutput(PersistentOStream & os) const {
  os << ounit(_deltam,GeV) << ounit(_deltagamma,GeV) << _pqmag << _pqphase
     << _pq << _zmag << _zphase << _z << _x << _y << _prob;
}

void MixedParticleData::persistentInput(PersistentIStream & is, int) {
  is >> iunit(_deltam,GeV) >> iunit(_deltagamma,GeV) >> _pqmag >> _pqphase
     >> _pq >> _zmag >> _zphase >> _z >> _x >> _y >> _prob;
}

ClassDescription<MixedParticleData> MixedParticleData::initMixedParticleData;
// Definition of the static class description member.

void MixedParticleData::Init() {

  static ClassDocumentation<MixedParticleData> documentation
    ("The MixedParticleData class provides storage of the particle data"
     " for particles which undergo mixing.");

  static Parameter<MixedParticleData,Energy> interfaceDeltaM
    ("DeltaM",
     "The mass difference",
     &MixedParticleData::_deltam, GeV, 0.0*GeV, 0.0*GeV, 1.*GeV,
     false, false, Interface::limited,
     &MixedParticleData::setDeltaM, 0, 0, 0, 0);

  static Parameter<MixedParticleData,Energy> interfaceDeltaGamma
    ("DeltaGamma",
     "The width difference",
     &MixedParticleData::_deltagamma, GeV, 0.0*GeV, 0.0*GeV, 1.0*GeV,
     false, false, Interface::limited,
     &MixedParticleData::setDeltaGamma, 0, 0, 0, 0);

  static Parameter<MixedParticleData,double> interfacePQMagnitude
    ("PQMagnitude",
     "The value of |p/q|",
     &MixedParticleData::_pqmag, 1.0, 0.0, 10.0,
     false, false, Interface::limited,
     &MixedParticleData::setPQMagnitude, 0, 0, 0, 0);

  static Parameter<MixedParticleData,double> interfacePQPhase
    ("PQPhase",
     "The phase of p/q",
     &MixedParticleData::_pqmag, 0.0, 0.0, 2.*Constants::pi,
     false, false, Interface::limited,
     &MixedParticleData::setPQPhase, 0, 0, 0, 0);

  static Parameter<MixedParticleData,double> interfaceZMagnitude
    ("ZMagnitude",
     "The value of |z|",
     &MixedParticleData::_zmag, 0.0, 0.0, 1.0,
     false, false, Interface::limited,
     &MixedParticleData::setZMagnitude, 0, 0, 0, 0);

  static Parameter<MixedParticleData,double> interfaceZPhase
    ("ZPhase",
     "The phase of z",
     &MixedParticleData::_zmag, 0.0, 0.0, 2.*Constants::pi,
     false, false, Interface::limited,
     &MixedParticleData::setZPhase, 0, 0, 0, 0);

}

MixedParticleData::MixedParticleData(long newId, string newPDGName)
  : ParticleData(newId, newPDGName), _deltam(0.*GeV), _deltagamma(0.*GeV),
    _pqmag(1.), _pqphase(0.), _pq(1.,0.), _zmag(0.), _zphase(0.), _z(0.),
    _x(0.), _y(0.), _prob(make_pair(1.,0.)) {}

PDPtr MixedParticleData::
Create(long newId, string newPDGName) {
  return new_ptr(MixedParticleData(newId, newPDGName));
}

PDPair MixedParticleData::
Create(long newId, string newPDGName, string newAntiPDGName) {
  PDPair pap;
  pap.first = new_ptr(MixedParticleData(newId, newPDGName));
  pap.second = new_ptr(MixedParticleData(-newId, newAntiPDGName));
  antiSetup(pap);
  return pap;
}

PDPtr MixedParticleData::pdclone() const {
  return new_ptr(*this);
}

void MixedParticleData::setDeltaM(Energy m) {
  _deltam = m;
  MixedParticleData * apd =
    dynamic_cast<MixedParticleData*>(CC().operator->());
  if ( synchronized() && apd ) apd->_deltam = m;
}

void MixedParticleData::setDeltaGamma(Energy m) {
  _deltagamma = m;
  MixedParticleData * apd =
    dynamic_cast<MixedParticleData*>(CC().operator->());
  if ( synchronized() && apd ) apd->_deltagamma = m;
}

void MixedParticleData::setPQMagnitude(double m) {
  _pqmag = m;
  MixedParticleData * apd =
    dynamic_cast<MixedParticleData*>(CC().operator->());
  if ( synchronized() && apd ) apd->_pqmag = m;
}

void MixedParticleData::setPQPhase(double m) {
  _pqphase = m;
  MixedParticleData * apd =
    dynamic_cast<MixedParticleData*>(CC().operator->());
  if ( synchronized() && apd ) apd->_pqphase = m;
}

void MixedParticleData::setZMagnitude(double m) {
  _zmag = m;
  MixedParticleData * apd =
    dynamic_cast<MixedParticleData*>(CC().operator->());
  if ( synchronized() && apd ) apd->_zmag = m;
}

void MixedParticleData::setZPhase(double m) {
  _zphase = m;
  MixedParticleData * apd =
    dynamic_cast<MixedParticleData*>(CC().operator->());
  if ( synchronized() && apd ) apd->_zphase = m;
}

void MixedParticleData::doinit() {
  ParticleData::doinit();
  // calculate the complex parameters from the magnitudes and phases
  // and x and y from massive parameters
  // p/q
  _pq = _pqmag*Complex(cos(_pqphase),sin(_pqphase));
  // z
  _z  = _zmag *Complex(cos(_zphase ),sin(_zphase ));
  // x
  _x =     _deltam    /width();
  // y
  _y = 0.5*_deltagamma/width();
  // probabilities
  double zr = _z.real(), zi = _z.imag();
  double root = sqrt( (1 - 2 * zr * zr + 2 * zi * zi +  pow( zr,  4) 
		       + 2 * zr * zr * zi * zi +  pow( zi,  4)));
  double x2=sqr(_x),y2=sqr(_y),modqp=1./sqr(abs(_pq)),z2(sqr(zr)+sqr(zi));
  double mixprob = id()>0 ?
    -modqp*root*(x2+y2)/(2*zr*_y*(1+x2) - modqp*root*(x2+y2) 
			 - (1.+z2)*x2 - 2*zi*_x*(1-y2) +  y2*(1-z2)  - 2) :
    root*(x2+y2)/(2*modqp*zr*(1+x2)*_y+root*(x2+y2)+modqp*(1.+z2)*x2 
		  - 2*modqp*zi*_x*(1-y2)-modqp*y2*(1-z2)+2*modqp);
  _prob= make_pair(1.-mixprob,mixprob);
  if( Debug::level > 1 ) {
    generator()->log() << "Parameters for the mixing of " << PDGName() << " and " 
		       << CC()->PDGName() << "\n";
    generator()->log() << "x = " << _x << "\t y = " << _y << "\n";
    generator()->log() << "Integrated mixing probability = " << mixprob << "\n";
  }
}

pair<bool,Length> MixedParticleData::generateLifeTime() const {
  // first decide if mixes
  bool mix = UseRandom::rndbool(_prob.second);
  double wgt;
  Length ct;
  double zi(_z.imag()),zr(_z.real()),zabs(sqr(zi)+sqr(zr)),
    root(1.-zabs);
  Length ctau = hbarc/(width()-0.5*abs(_deltagamma));
  do {
    ct = UseRandom::rndExp(ctau);
    double gt = ct/cTau();
    if(id()>0) {
      if(!mix) {
	wgt = 0.5*(1.+zabs)*cosh(_y*gt)+0.5*(1.-zabs)*cos(_x*gt)
	  -zr*sinh(_y*gt)+zi*sin(_x*gt);
      }
      else {
	wgt = 0.5*root/sqr(abs(_pq))*(cosh(_y*gt)-cos(_x*gt));
      }
    }
    else {
      if(!mix) {
	wgt = 0.5*(1.+zabs)*cosh(_y*gt)+0.5*(1.-zabs)*cos(_x*gt)
	  +zr*sinh(_y*gt)-zi*sin(_x*gt);
      }
      else {
	wgt = 0.5*root*sqr(abs(_pq))*(cosh(_y*gt)-cos(_x*gt));
      }
    }
    wgt *= exp(-gt+ct/ctau);
  }
  while(UseRandom::rnd()>wgt);
  return make_pair(mix,ct);
}

pair<Complex,Complex> MixedParticleData::mixingAmplitudes(Length ct,bool part) const {
  double gt = ct/cTau();
  Complex ep = exp(Complex(-0.5*_y,-0.5*_x)*gt);
  Complex gp = 0.5*(ep+1./ep), gm = 0.5*(ep-1./ep);
  pair<Complex,Complex> output;
  if(part) {
    output.first  = gp + _z*gm;
    output.second = -sqrt(1.-sqr(_z))/_pq*gm;
  }
  else {
    output.first  = gp - _z*gm;
    output.second = -sqrt(1.-sqr(_z))*_pq*gm; 
  }
  return output;
}
