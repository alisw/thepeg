// -*- C++ -*-
//
// TensorWaveFunction.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2017 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TensorWaveFunction class.
//
// Author: Peter Richardson
//
#include "TensorWaveFunction.h"

using namespace ThePEG;
using namespace ThePEG::Helicity;

// calculate the actual wavefunction
void TensorWaveFunction::calculateWaveFunction(unsigned int ihel, TensorPhase tphase) {
  int jhel=ihel-2;
  assert(direction()!=intermediate);
  // check for a valid helicty combination
  assert((jhel<=2 && jhel>=-2   && mass() >ZERO) || 
	 ((jhel==2 || jhel==-2) && mass()==ZERO));
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
  if(jhel>=0) {
    // calculate the overall phase
    complex<double>phase;
    if(tphase==tensor_phase) {
      phase = pt==ZERO ? 1. : complex<double>(ppx/pt,-fact*ppy/pt);
    }
    else phase = 1.; 
    phase = phase*sqrt(0.5);
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
      epsp[2]=pt*opabs*phase;
      epsp[3]=0.;
    }
  }
  // - helicity vector if needed
  if(jhel<=0) {
    // calculate the overall phase
    complex<double> phase;
    if(tphase==tensor_phase) {
      phase = pt==ZERO ? 1. : complex<double>(ppx/pt,fact*ppy/pt);
    }
    else phase = 1.;
    phase = phase*sqrt(0.5);
    // first the no pt case
    if(pt==ZERO) {
      double sgnz;
      if(ppz<ZERO){sgnz=-1.;}
      else{sgnz=1.;}
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
      epsm[2]=-pt*opabs*phase;
      epsm[3]=0.;
    }
  }
  // 0 helicity vector if needed
  if(jhel<=1 && jhel>=-1) {
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
  case 2:
    for(int ix=0;ix<4;++ix)
      for(int iy=0;iy<4;++iy) _wf(ix,iy)=epsp[ix]*epsp[iy];
    break;
  case 1:
    ort = sqrt(0.5);
    for(int ix=0;ix<4;++ix)
      for(int iy=0;iy<4;++iy) _wf(ix,iy)=ort*( epsp[ix]*eps0[iy]+
					       eps0[ix]*epsp[iy]);
    break;
  case 0:
    ort = 1./sqrt(6.);
    for(int ix=0;ix<4;++ix)
      for(int iy=0;iy<4;++iy) _wf(ix,iy)=ort*(    epsp[ix]*epsm[iy]
					       +  epsm[ix]*epsp[iy]
					      +2.*eps0[ix]*eps0[iy]);
    break;
  case -1:
    ort = 1./sqrt(2.);
    for(int ix=0;ix<4;++ix)
      for(int iy=0;iy<4;++iy) _wf(ix,iy)=ort*( epsm[ix]*eps0[iy]+
					       eps0[ix]*epsm[iy]);
    break;
  case -2:
    for(int ix=0;ix<4;++ix)
      for(int iy=0;iy<4;++iy) _wf(ix,iy)=epsm[ix]*epsm[iy];
    break;
  default:
    ThePEG::Helicity::HelicityConsistencyError() 
      << "Invalid Helicity = " << jhel << " requested for Tensor" 
      << Exception::abortnow;
    break;
  }
}

void TensorWaveFunction::
calculateWaveFunctions(vector<LorentzTensor<double> > & waves,
		       tPPtr particle,Direction dir, bool massless,
		       TensorPhase phase) {
  tTensorSpinPtr inspin = !particle->spinInfo() ? tTensorSpinPtr() :
    dynamic_ptr_cast<tTensorSpinPtr>(particle->spinInfo());
  waves.resize(5);
  if(inspin) {
    if(dir==outgoing) {
      for(unsigned int ix=0;ix<5;++ix)
	waves[ix]=inspin->getProductionBasisState(ix);
    }
    else {
      inspin->decay();
      for(unsigned int ix=0;ix<5;++ix)
	waves[ix]=inspin->getDecayBasisState(ix);
    }
  }
  else {
    assert(!particle->spinInfo());
    TensorWaveFunction wave(particle->momentum(),particle->dataPtr(),0,
			    dir,phase);
    for(unsigned int ix=0;ix<5;++ix) {
      if(massless&&ix>0&&ix<5) {
	waves[ix] = LorentzTensor<double>();
      }
      else {
	if(ix!=0) wave.reset(ix);
	waves[ix] = wave.wave();
      }
    }
  }
}

void TensorWaveFunction::
calculateWaveFunctions(vector<TensorWaveFunction> & waves,
		       tPPtr particle, Direction dir, bool massless,
		       TensorPhase phase) {
  tTensorSpinPtr inspin = !particle->spinInfo() ? tTensorSpinPtr() :
    dynamic_ptr_cast<tTensorSpinPtr>(particle->spinInfo());
  waves.resize(5);
  if(inspin) {
    if(dir==outgoing) {
      for(unsigned int ix=0;ix<5;++ix)
	waves[ix]=TensorWaveFunction(particle->momentum(),
				     particle->dataPtr(),
				     inspin->getProductionBasisState(ix),dir);
    }
    else {
      inspin->decay();
      for(unsigned int ix=0;ix<5;++ix)
	waves[ix]=TensorWaveFunction(particle->momentum(),
				     particle->dataPtr(),
				     inspin->getDecayBasisState(ix),dir);
    }
  }
  else {
    assert(!particle->spinInfo());
    TensorWaveFunction wave(particle->momentum(),particle->dataPtr(),0,
			    dir,phase);
    for(unsigned int ix=0;ix<5;++ix) {
      if(massless&&ix>0&&ix<5) {
	waves[ix] = TensorWaveFunction(particle->momentum(),particle->dataPtr(),dir);
      }
      else {
	if(ix!=0) wave.reset(ix);
	waves[ix] = wave;
      }
    }
  }
}

void  TensorWaveFunction::
calculateWaveFunctions(vector<LorentzTensor<double> > & waves,
		       RhoDMatrix & rho,
		       tPPtr particle,Direction dir,bool massless,
		       TensorPhase phase) {
  tTensorSpinPtr inspin = !particle->spinInfo() ? tTensorSpinPtr() :
    dynamic_ptr_cast<tTensorSpinPtr>(particle->spinInfo());
  waves.resize(5);
  if(inspin) {
    if(dir==outgoing) {
      for(unsigned int ix=0;ix<5;++ix)
	waves[ix]=inspin->getProductionBasisState(ix);
      rho = RhoDMatrix(PDT::Spin2);
    }
    else {
      inspin->decay();
      for(unsigned int ix=0;ix<5;++ix)
	waves[ix]=inspin->getDecayBasisState(ix);
      rho = inspin->rhoMatrix();
    }
  }
  else {
    assert(!particle->spinInfo());
    TensorWaveFunction wave(particle->momentum(),particle->dataPtr(),0,
			    dir,phase);
    for(unsigned int ix=0;ix<5;++ix) {
      if(massless&&ix>0&&ix<5) {
	waves[ix] = LorentzTensor<double>();
      }
      else {
	if(ix!=0) wave.reset(ix);
	waves[ix] = wave.wave();
      }
    }
    rho = RhoDMatrix(PDT::Spin2);
  }
}

void  TensorWaveFunction::
calculateWaveFunctions(vector<TensorWaveFunction> & waves,
		       RhoDMatrix & rho,
		       tPPtr particle, Direction dir, bool massless,
		       TensorPhase phase) {
  tTensorSpinPtr inspin = !particle->spinInfo() ? tTensorSpinPtr() :
    dynamic_ptr_cast<tTensorSpinPtr>(particle->spinInfo());
  waves.resize(5);
  if(inspin) {
    if(dir==outgoing) {
      for(unsigned int ix=0;ix<5;++ix)
	waves[ix]=TensorWaveFunction(particle->momentum(),
				     particle->dataPtr(),
				     inspin->getProductionBasisState(ix),dir);
      rho = RhoDMatrix(PDT::Spin2);
    }
    else {
      inspin->decay();
      for(unsigned int ix=0;ix<5;++ix)
	waves[ix]=TensorWaveFunction(particle->momentum(),
				     particle->dataPtr(),
				     inspin->getDecayBasisState(ix),dir);
      rho = inspin->rhoMatrix();
    }
  }
  else {
    assert(!particle->spinInfo());
    TensorWaveFunction wave(particle->momentum(),particle->dataPtr(),0,
			    dir,phase);
    for(unsigned int ix=0;ix<5;++ix) {
      if(massless&&ix>0&&ix<5) {
	waves[ix] = TensorWaveFunction(particle->momentum(),particle->dataPtr(),dir);
      }
      else {
	if(ix!=0) wave.reset(ix);
	waves[ix] = wave;
      }
    }
    rho = RhoDMatrix(PDT::Spin2);
  }
}

void  TensorWaveFunction::
constructSpinInfo(const vector<LorentzTensor<double> > & waves,
		  tPPtr part,Direction dir, bool time,bool ) {
  assert(waves.size()==5);
  tTensorSpinPtr inspin = !part->spinInfo() ? tTensorSpinPtr() :
    dynamic_ptr_cast<tTensorSpinPtr>(part->spinInfo());
  if(inspin) {
    for(unsigned int ix=0;ix<5;++ix)
      if(dir==outgoing) inspin->setBasisState(ix,waves[ix]);
      else              inspin->setDecayState(ix,waves[ix]);
  }
  else {
    TensorSpinPtr temp = new_ptr(TensorSpinInfo(part->momentum(),time));
    part->spinInfo(temp);
    for(unsigned int ix=0;ix<5;++ix)
      if(dir==outgoing) temp->setBasisState(ix,waves[ix]);
      else              temp->setDecayState(ix,waves[ix]);
  }
}

void  TensorWaveFunction::
constructSpinInfo(const vector<TensorWaveFunction> & waves,
		  tPPtr part,Direction dir, bool time,bool ) {
  assert(waves.size()==5);
  tTensorSpinPtr inspin = !part->spinInfo() ? tTensorSpinPtr() :
    dynamic_ptr_cast<tTensorSpinPtr>(part->spinInfo());
  if(inspin) {
    for(unsigned int ix=0;ix<5;++ix)
      if(dir==outgoing) inspin->setBasisState(ix,waves[ix].wave());
      else              inspin->setDecayState(ix,waves[ix].wave());
  }
  else {
    TensorSpinPtr temp = new_ptr(TensorSpinInfo(part->momentum(),time));
    part->spinInfo(temp);
    for(unsigned int ix=0;ix<5;++ix)
      if(dir==outgoing) temp->setBasisState(ix,waves[ix].wave());
      else              temp->setDecayState(ix,waves[ix].wave());
  }
}
