// -*- C++ -*-
//
// Rank3TensorWaveFunction.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2019 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Rank3TensorWaveFunction class.
//
// Author: Peter Richardson
//
#include "Rank3TensorWaveFunction.h"

using namespace ThePEG;
using namespace ThePEG::Helicity;

// calculate the actual wavefunction
void Rank3TensorWaveFunction::calculateWaveFunction(unsigned int ihel) {
  int jhel=ihel-3;
  assert(direction()!=intermediate);
  // check for a valid helicty combination
  assert( (jhel<=3 && jhel>=-3   && mass() >ZERO) || 
	 ((jhel==3 || jhel==-3)  && mass()==ZERO));
  // extract the momentum components
  double fact = direction()==outgoing ? -1. : 1;
  Energy ppx=fact*px(),ppy=fact*py(),ppz=fact*pz(),pee=fact*e(),pmm=mass();
  // calculate some kinematic quantites;
  Energy2 pt2 = sqr(ppx)+sqr(ppy);
  Energy pabs = sqrt(pt2+sqr(ppz));
  Energy pt = sqrt(pt2);
  // polarization vectors
  complex<double> epsp[4],epsm[4],eps0[4];
  // + helicity vector if needed
  if(jhel>=-1) {
    // calculate the overall phase
    complex<double> phase = pt==ZERO ? sqrt(0.5) : sqrt(0.5)*complex<double>(ppx/pt,-fact*ppy/pt);
    // first the no pt case
    if(pt==ZERO) {
      double sgnz = ppz<ZERO ? -1. : 1.;
      epsp[0]=-phase;
      epsp[1]= sgnz*phase*complex<double>(0,-fact);
      epsp[2]=0.;
      epsp[3]=0.;
    }
    else {
      InvEnergy opabs=1./pabs;
      InvEnergy opt  =1./pt;
      epsp[0]=phase*complex<double>(-ppz*ppx*opabs*opt,
				    fact*ppy*opt);
      epsp[1]=phase*complex<double>(-ppz*ppy*opabs*opt,
				    -fact*ppx*opt);
      epsp[2]=Complex(pt*opabs*phase);
      epsp[3]=0.;
    }
  }
  // - helicity vector if needed
  if(jhel<=1) {
    // calculate the overall phase
    complex<double> phase = pt==ZERO ? sqrt(0.5) : sqrt(0.5)*complex<double>(ppx/pt,fact*ppy/pt);
    // first the no pt case
    if(pt==ZERO) {
      double sgnz = ppz<ZERO ? -1. : 1.;
      epsm[0]= phase;
      epsm[1]= sgnz*phase*complex<double>(0,-fact);
      epsm[2]=0.;
      epsm[3]=0.;
    }
    else {
      InvEnergy opabs=1./pabs;
      InvEnergy opt  =1./pt;
      epsm[0]=phase*complex<double>(ppz*ppx*opabs*opt,
				    fact*ppy*opt);
      epsm[1]=phase*complex<double>(ppz*ppy*opabs*opt,
				    -fact*ppx*opt);
      epsm[2]=-Complex(pt*opabs*phase);
      epsm[3]=0.;
    }
  }
  // 0 helicity vector if needed
  if(abs(jhel)!=3) {
    if(pabs==ZERO) {
      eps0[0] = 0.;
      eps0[1] = 0.;
      eps0[2] = 1.;
      eps0[3] = 0.;
    }
    else {
      InvEnergy empabs=pee/pmm/pabs;
      eps0[0] = empabs*ppx;
      eps0[1] = empabs*ppy;
      eps0[2] = empabs*ppz;
      eps0[3] = pabs/pmm;
    }
  }
  // put the polarization vectors together to get the wavefunction
  double ort;
  switch (jhel) {
  case 3:
    for(int ix=0;ix<4;++ix)
      for(int iy=0;iy<4;++iy)
	for(int iz=0;iz<4;++iz) _wf(ix,iy,iz)=epsp[ix]*epsp[iy]*epsp[iz];
    break;
  case 2:
    ort = sqrt(1./3.);
    for(int ix=0;ix<4;++ix)
      for(int iy=0;iy<4;++iy)
	for(int iz=0;iz<4;++iz)
	  _wf(ix,iy,iz)= ort*(epsp[ix]*epsp[iy]*eps0[iz] +
			      epsp[ix]*eps0[iy]*epsp[iz] +
			      eps0[ix]*epsp[iy]*epsp[iz]);
    break;
  case 1:
    ort = sqrt(1./15.);
    for(int ix=0;ix<4;++ix)
      for(int iy=0;iy<4;++iy)
	for(int iz=0;iz<4;++iz)
	  _wf(ix,iy,iz)= ort*(epsp[ix]*epsp[iy]*epsm[iz] +
			      epsp[ix]*epsm[iy]*epsp[iz] +
			      epsm[ix]*epsp[iy]*epsp[iz]
			      +2.*(eps0[ix]*eps0[iy]*epsp[iz] +
				   eps0[ix]*epsp[iy]*eps0[iz] +
				   epsp[ix]*eps0[iy]*eps0[iz]));
    break;
  case 0:
    ort = sqrt(0.1);
    for(int ix=0;ix<4;++ix)
      for(int iy=0;iy<4;++iy)
	for(int iz=0;iz<4;++iz)
	  _wf(ix,iy,iz)= ort*(epsp[ix]*epsm[iy]*eps0[iz] +
			      epsp[ix]*eps0[iy]*epsm[iz] +
			      eps0[ix]*epsp[iy]*epsm[iz] +
			      eps0[ix]*epsm[iy]*epsp[iz] +
			      epsm[ix]*epsp[iy]*eps0[iz] +
			      epsm[ix]*eps0[iy]*epsp[iz] +
			      2.*eps0[ix]*eps0[iy]*eps0[iz]);
    break;
  case -1:
    ort = sqrt(1./15.);
    for(int ix=0;ix<4;++ix)
      for(int iy=0;iy<4;++iy)
	for(int iz=0;iz<4;++iz)
	  _wf(ix,iy,iz)= ort*(epsm[ix]*epsm[iy]*epsp[iz] +
			      epsm[ix]*epsp[iy]*epsm[iz] +
			      epsp[ix]*epsm[iy]*epsm[iz]
			      +2.*(eps0[ix]*eps0[iy]*epsm[iz] +
				   eps0[ix]*epsm[iy]*eps0[iz] +
				   epsm[ix]*eps0[iy]*eps0[iz]));
    break;
  case -2:
    ort = sqrt(1./3.);
    for(int ix=0;ix<4;++ix)
      for(int iy=0;iy<4;++iy)
	for(int iz=0;iz<4;++iz)
	  _wf(ix,iy,iz)= ort*(epsm[ix]*epsm[iy]*eps0[iz] +
			      epsm[ix]*eps0[iy]*epsm[iz] +
			      eps0[ix]*epsm[iy]*epsm[iz]);
    break;;
  case -3:
    for(int ix=0;ix<4;++ix)
      for(int iy=0;iy<4;++iy)
	for(int iz=0;iz<4;++iz) _wf(ix,iy,iz)=epsm[ix]*epsm[iy]*epsm[iz];
    break;
  default:
    ThePEG::Helicity::HelicityConsistencyError() 
      << "Invalid Helicity = " << jhel << " requested for Rank3Tensor" 
      << Exception::abortnow;
    break;
  }
}

void Rank3TensorWaveFunction::
calculateWaveFunctions(vector<LorentzRank3Tensor<double> > & waves,
		       tPPtr particle,Direction dir, bool massless) {
  tRank3TensorSpinPtr inspin = !particle->spinInfo() ? tRank3TensorSpinPtr() :
    dynamic_ptr_cast<tRank3TensorSpinPtr>(particle->spinInfo());
  waves.resize(7);
  if(inspin) {
    if(dir==outgoing) {
      for(unsigned int ix=0;ix<7;++ix)
	waves[ix]=inspin->getProductionBasisState(ix);
    }
    else {
      inspin->decay();
      for(unsigned int ix=0;ix<7;++ix)
	waves[ix]=inspin->getDecayBasisState(ix);
    }
  }
  else {
    assert(!particle->spinInfo());
    Rank3TensorWaveFunction wave(particle->momentum(),particle->dataPtr(),0,dir);
    for(unsigned int ix=0;ix<7;++ix) {
      if(massless&&ix>0&&ix<7) {
	waves[ix] = LorentzRank3Tensor<double>();
      }
      else {
	if(ix!=0) wave.reset(ix);
	waves[ix] = wave.wave();
      }
    }
  }
}

void Rank3TensorWaveFunction::
calculateWaveFunctions(vector<Rank3TensorWaveFunction> & waves,
		       tPPtr particle, Direction dir, bool massless) {
  tRank3TensorSpinPtr inspin = !particle->spinInfo() ? tRank3TensorSpinPtr() :
    dynamic_ptr_cast<tRank3TensorSpinPtr>(particle->spinInfo());
  waves.resize(7);
  if(inspin) {
    if(dir==outgoing) {
      for(unsigned int ix=0;ix<7;++ix)
	waves[ix]=Rank3TensorWaveFunction(particle->momentum(),
				     particle->dataPtr(),
				     inspin->getProductionBasisState(ix),dir);
    }
    else {
      inspin->decay();
      for(unsigned int ix=0;ix<7;++ix)
	waves[ix]=Rank3TensorWaveFunction(particle->momentum(),
				     particle->dataPtr(),
				     inspin->getDecayBasisState(ix),dir);
    }
  }
  else {
    assert(!particle->spinInfo());
    Rank3TensorWaveFunction wave(particle->momentum(),particle->dataPtr(),0,dir);
    for(unsigned int ix=0;ix<7;++ix) {
      if(massless&&ix>0&&ix<7) {
	waves[ix] = Rank3TensorWaveFunction(particle->momentum(),particle->dataPtr(),dir);
      }
      else {
	if(ix!=0) wave.reset(ix);
	waves[ix] = wave;
      }
    }
  }
}

void  Rank3TensorWaveFunction::
calculateWaveFunctions(vector<LorentzRank3Tensor<double> > & waves,
		       RhoDMatrix & rho,
		       tPPtr particle,Direction dir,bool massless) {
  tRank3TensorSpinPtr inspin = !particle->spinInfo() ? tRank3TensorSpinPtr() :
    dynamic_ptr_cast<tRank3TensorSpinPtr>(particle->spinInfo());
  waves.resize(7);
  if(inspin) {
    if(dir==outgoing) {
      for(unsigned int ix=0;ix<7;++ix)
	waves[ix]=inspin->getProductionBasisState(ix);
      rho = RhoDMatrix(PDT::Spin3);
    }
    else {
      inspin->decay();
      for(unsigned int ix=0;ix<7;++ix)
	waves[ix]=inspin->getDecayBasisState(ix);
      rho = inspin->rhoMatrix();
    }
  }
  else {
    assert(!particle->spinInfo());
    Rank3TensorWaveFunction wave(particle->momentum(),particle->dataPtr(),0,dir);
    for(unsigned int ix=0;ix<7;++ix) {
      if(massless&&ix>0&&ix<7) {
	waves[ix] = LorentzRank3Tensor<double>();
      }
      else {
	if(ix!=0) wave.reset(ix);
	waves[ix] = wave.wave();
      }
    }
    rho = RhoDMatrix(PDT::Spin3);
  }
}

void  Rank3TensorWaveFunction::
calculateWaveFunctions(vector<Rank3TensorWaveFunction> & waves,
		       RhoDMatrix & rho,
		       tPPtr particle, Direction dir, bool massless) {
  tRank3TensorSpinPtr inspin = !particle->spinInfo() ? tRank3TensorSpinPtr() :
    dynamic_ptr_cast<tRank3TensorSpinPtr>(particle->spinInfo());
  waves.resize(7);
  if(inspin) {
    if(dir==outgoing) {
      for(unsigned int ix=0;ix<7;++ix)
	waves[ix]=Rank3TensorWaveFunction(particle->momentum(),
				     particle->dataPtr(),
				     inspin->getProductionBasisState(ix),dir);
      rho = RhoDMatrix(PDT::Spin3);
    }
    else {
      inspin->decay();
      for(unsigned int ix=0;ix<7;++ix)
	waves[ix]=Rank3TensorWaveFunction(particle->momentum(),
				     particle->dataPtr(),
				     inspin->getDecayBasisState(ix),dir);
      rho = inspin->rhoMatrix();
    }
  }
  else {
    assert(!particle->spinInfo());
    Rank3TensorWaveFunction wave(particle->momentum(),particle->dataPtr(),0,dir);
    for(unsigned int ix=0;ix<7;++ix) {
      if(massless&&ix>0&&ix<7) {
	waves[ix] = Rank3TensorWaveFunction(particle->momentum(),particle->dataPtr(),dir);
      }
      else {
	if(ix!=0) wave.reset(ix);
	waves[ix] = wave;
      }
    }
    rho = RhoDMatrix(PDT::Spin3);
  }
}

void  Rank3TensorWaveFunction::
constructSpinInfo(const vector<LorentzRank3Tensor<double> > & waves,
		  tPPtr part,Direction dir, bool time,bool ) {
  assert(waves.size()==7);
  tRank3TensorSpinPtr inspin = !part->spinInfo() ? tRank3TensorSpinPtr() :
    dynamic_ptr_cast<tRank3TensorSpinPtr>(part->spinInfo());
  if(inspin) {
    for(unsigned int ix=0;ix<7;++ix)
      if(dir==outgoing) inspin->setBasisState(ix,waves[ix]);
      else              inspin->setDecayState(ix,waves[ix]);
  }
  else {
    Rank3TensorSpinPtr temp = new_ptr(Rank3TensorSpinInfo(part->momentum(),time));
    part->spinInfo(temp);
    for(unsigned int ix=0;ix<7;++ix)
      if(dir==outgoing) temp->setBasisState(ix,waves[ix]);
      else              temp->setDecayState(ix,waves[ix]);
  }
}

void  Rank3TensorWaveFunction::
constructSpinInfo(const vector<Rank3TensorWaveFunction> & waves,
		  tPPtr part,Direction dir, bool time,bool ) {
  assert(waves.size()==7);
  tRank3TensorSpinPtr inspin = !part->spinInfo() ? tRank3TensorSpinPtr() :
    dynamic_ptr_cast<tRank3TensorSpinPtr>(part->spinInfo());
  if(inspin) {
    for(unsigned int ix=0;ix<7;++ix)
      if(dir==outgoing) inspin->setBasisState(ix,waves[ix].wave());
      else              inspin->setDecayState(ix,waves[ix].wave());
  }
  else {
    Rank3TensorSpinPtr temp = new_ptr(Rank3TensorSpinInfo(part->momentum(),time));
    part->spinInfo(temp);
    for(unsigned int ix=0;ix<7;++ix)
      if(dir==outgoing) temp->setBasisState(ix,waves[ix].wave());
      else              temp->setDecayState(ix,waves[ix].wave());
  }
}
