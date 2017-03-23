// -*- C++ -*-
//
// SpinorBarWaveFunction.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2011 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SpinorBarWaveFunction class.
//
// Author: Peter Richardson
//

#include "SpinorBarWaveFunction.h"
#include "SpinorWaveFunction.h"

using namespace ThePEG;
  using namespace Helicity;

// calculate the Wavefunction
void SpinorBarWaveFunction::calculateWaveFunction(unsigned int ihel) {
  Direction dir=direction();
  if(dir==intermediate) ThePEG::Helicity::HelicityConsistencyError() 
    << "In SpinorBarWaveFunction::calcluateWaveFunction "
    << "particle must be incoming or outgoing not intermediate" 
    << Exception::abortnow;
  // check ihelicity is O.K.
  if(ihel>1)  ThePEG::Helicity::HelicityConsistencyError() 
    << "Invalid Helicity = " << ihel << " requested for SpinorBar" 
    << Exception::abortnow;
  // extract the momentum components
  double fact = dir==incoming ? 1. : -1.;
  Energy ppx=fact*px(),ppy=fact*py(),ppz=fact*pz(),pee=fact*e(),pmm=mass();
  // define and calculate some kinematic quantities
  Energy2 ptran2  = ppx*ppx+ppy*ppy;
  Energy pabs   = sqrt(ptran2+ppz*ppz);
  Energy ptran  = sqrt(ptran2);
  // first need to evalulate the 2-component helicity spinors 
  // this is the same regardless of which definition of the spinors
  // we are using
  Complex hel_wf[2];
  // compute the + spinor for + helicty particles and - helicity antiparticles
  if((dir==outgoing && ihel== 1) || (dir==incoming && ihel==0)) {
    // no transverse momentum 
    if(ptran==ZERO) {
      if(ppz>=ZERO) {
	hel_wf[0] = 1;
	hel_wf[1] = 0;
      }
      else {
	hel_wf[0] = 0;
	hel_wf[1] = 1;
      }
    }
    else {
      InvSqrtEnergy denominator = 1./sqrt(2.*pabs);
      SqrtEnergy rtppluspz = (ppz>=ZERO) ? sqrt(pabs+ppz) : ptran/sqrt(pabs-ppz); 
      hel_wf[0] = denominator*rtppluspz;
      hel_wf[1] = denominator/rtppluspz*complex<Energy>(ppx,-ppy);
    }
  }
  // compute the - spinor for - helicty particles and + helicity antiparticles
  else {
    // no transverse momentum
    if(ptran==ZERO) {
      if(ppz>=ZERO) {
	hel_wf[0] = 0;
	hel_wf[1] = 1;
      }
      // transverse momentum 
      else {
	hel_wf[0] = -1;
	hel_wf[1] =  0;
      }
    }
    else {
      InvSqrtEnergy denominator = 1./sqrt(2.*pabs);
      SqrtEnergy rtppluspz = (ppz>=ZERO) ? sqrt(pabs+ppz) : ptran/sqrt(pabs-ppz);
      hel_wf[0] = denominator/rtppluspz*complex<Energy>(-ppx,-ppy);
      hel_wf[1] = denominator*rtppluspz;
    }
  }
  SqrtEnergy upper, lower;
  SqrtEnergy eplusp  = sqrt(max(pee+pabs,ZERO));
  SqrtEnergy eminusp = ( pmm!=ZERO ) ? pmm/eplusp : ZERO;
  // set up the coefficients for the different cases
  if(dir==outgoing) {
    if(ihel==1) {
      upper = eplusp;
      lower = eminusp;
    }
    else {
      upper = eminusp;
      lower = eplusp;
    }
  }
  else {
    if(ihel==1) {
	upper = eminusp;
	lower = -eplusp;
    }
    else {
      upper =-eplusp;
      lower = eminusp;
    }
  }
  // now finally we can construct the spinors
  _wf = LorentzSpinorBar<double>((dir==incoming) ? v_spinortype : u_spinortype);
  _wf[0] = upper*hel_wf[0]*UnitRemoval::InvSqrtE;
  _wf[1] = upper*hel_wf[1]*UnitRemoval::InvSqrtE;
  _wf[2] = lower*hel_wf[0]*UnitRemoval::InvSqrtE;
  _wf[3] = lower*hel_wf[1]*UnitRemoval::InvSqrtE;
}

void SpinorBarWaveFunction::conjugate() {
  _wf=_wf.conjugate();
}

SpinorWaveFunction SpinorBarWaveFunction::bar() {
  Lorentz5Momentum p = momentum();
  if(direction()==outgoing) p *= -1.;
  tcPDPtr ptemp = particle();
  if(direction()==incoming&&particle()->CC())
    ptemp = particle()->CC();
  return SpinorWaveFunction(p,ptemp,_wf.bar(),direction());
}	  

void SpinorBarWaveFunction::
calculateWaveFunctions(vector<LorentzSpinorBar<SqrtEnergy> > & waves,
		       tPPtr particle,Direction dir) {
  tFermionSpinPtr inspin = !particle->spinInfo() ? tFermionSpinPtr() :
    dynamic_ptr_cast<tFermionSpinPtr>(particle->spinInfo());
  waves.resize(2);
  // spin info object exists
  if(inspin) {
    if(dir==outgoing) {
      waves[0] = inspin->getProductionBasisState(0).bar();
      waves[1] = inspin->getProductionBasisState(1).bar();
    }
    else {
      inspin->decay();
      if( (particle->id()>0&&inspin->getDecayBasisState(0).Type()!=u_spinortype) || 
	  (particle->id()<0&&inspin->getDecayBasisState(0).Type()!=v_spinortype)) {
	waves[0] = inspin->getDecayBasisState(0).conjugate().bar();
	waves[1] = inspin->getDecayBasisState(1).conjugate().bar();
      }
      else {
	waves[0] = inspin->getDecayBasisState(0).bar();
	waves[1] = inspin->getDecayBasisState(1).bar();
      }
    }
  }
  // do the calculation
  else {
    assert(!particle->spinInfo());
    SpinorBarWaveFunction wave(particle->momentum(),particle->dataPtr(),dir);
    for(unsigned int ix=0;ix<2;++ix) {
      wave.reset(ix);
      waves[ix] = wave.dimensionedWave();
    }
  }
}

void SpinorBarWaveFunction::
calculateWaveFunctions(vector<SpinorBarWaveFunction> & waves,
		       tPPtr particle,Direction dir) {
  tFermionSpinPtr inspin = !particle->spinInfo() ? tFermionSpinPtr() :
    dynamic_ptr_cast<tFermionSpinPtr>(particle->spinInfo());
  waves.resize(2);
  // spin info object exists
  if(inspin) {
    if(dir==outgoing) {
      for(unsigned int ix=0;ix<2;++ix)
	waves[ix] = SpinorBarWaveFunction(particle,
					  inspin->getProductionBasisState(ix).bar(),
					  dir);
    }
    else {
      inspin->decay();
      if((particle->id()>0&&inspin->getDecayBasisState(0).Type()!=u_spinortype) ||
	 (particle->id()<0&&inspin->getDecayBasisState(0).Type()!=v_spinortype)) {
	for(unsigned int ix=0;ix<2;++ix)
	  waves[ix] = SpinorBarWaveFunction(particle,
					    inspin->getDecayBasisState(ix).conjugate().bar(),dir);
      }
      else {
	for(unsigned int ix=0;ix<2;++ix)
	  waves[ix] = SpinorBarWaveFunction(particle,
					    inspin->getDecayBasisState(ix).bar(),dir);
      }
    }
  }
  // do the calculation
  else {
    assert(!particle->spinInfo());
    SpinorBarWaveFunction wave(particle->momentum(),particle->dataPtr(),dir);
    for(unsigned int ix=0;ix<2;++ix) {
      wave.reset(ix);
      waves[ix] = wave;
    }
  }
}

void SpinorBarWaveFunction::
calculateWaveFunctions(vector<LorentzSpinorBar<SqrtEnergy> > & waves,
		       RhoDMatrix & rho,
		       tPPtr particle,Direction dir) {
  tFermionSpinPtr inspin = !particle->spinInfo() ? tFermionSpinPtr() :
    dynamic_ptr_cast<tFermionSpinPtr>(particle->spinInfo());
  waves.resize(2);
  // spin info object exists
  if(inspin) {
    if(dir==outgoing) {
      waves[0] = inspin->getProductionBasisState(0).bar();
      waves[1] = inspin->getProductionBasisState(1).bar();
      rho = RhoDMatrix(PDT::Spin1Half);
    }
    else {
      inspin->decay();
      if((particle->id()>0&&inspin->getDecayBasisState(0).Type()!=u_spinortype) ||
	 (particle->id()<0&&inspin->getDecayBasisState(0).Type()!=v_spinortype)) {
	waves[0] = inspin->getDecayBasisState(0).conjugate().bar();
	waves[1] = inspin->getDecayBasisState(1).conjugate().bar();
      }
      else {
	waves[0] = inspin->getDecayBasisState(0).bar();
	waves[1] = inspin->getDecayBasisState(1).bar();
      }
      rho = inspin->rhoMatrix();
    }
  }
  // do the calculation
  else {
    assert(!particle->spinInfo());
    SpinorBarWaveFunction wave(particle->momentum(),particle->dataPtr(),dir);
    for(unsigned int ix=0;ix<2;++ix) {
      wave.reset(ix);
      waves[ix] = wave.dimensionedWave();
    }
    rho = RhoDMatrix(PDT::Spin1Half);
  }
}

void SpinorBarWaveFunction::
calculateWaveFunctions(vector<SpinorBarWaveFunction> & waves,
		       RhoDMatrix & rho,
		       tPPtr particle,Direction dir) {
  tFermionSpinPtr inspin = !particle->spinInfo() ? tFermionSpinPtr() :
    dynamic_ptr_cast<tFermionSpinPtr>(particle->spinInfo());
  waves.resize(2);
  // spin info object exists
  if(inspin) {
    if(dir==outgoing) {
      for(unsigned int ix=0;ix<2;++ix)
	waves[ix] = SpinorBarWaveFunction(particle,
					  inspin->getProductionBasisState(ix).bar(),
					  dir);
      rho = RhoDMatrix(PDT::Spin1Half);
    }
    else {
      inspin->decay();
      if((particle->id()>0&&inspin->getDecayBasisState(0).Type()!=u_spinortype) ||
	 (particle->id()<0&&inspin->getDecayBasisState(0).Type()!=v_spinortype)) {
	for(unsigned int ix=0;ix<2;++ix)
	  waves[ix] = SpinorBarWaveFunction(particle,
					    inspin->getDecayBasisState(ix).conjugate().bar(),dir);
      }
      else {
	for(unsigned int ix=0;ix<2;++ix)
	  waves[ix] = SpinorBarWaveFunction(particle,
					    inspin->getDecayBasisState(ix).bar(),dir);
      }
      rho = inspin->rhoMatrix();
    }
  }
  // do the calculation
  else {
    assert(!particle->spinInfo());
    SpinorBarWaveFunction wave(particle->momentum(),particle->dataPtr(),dir);
    for(unsigned int ix=0;ix<2;++ix) {
      wave.reset(ix);
      waves[ix] = wave;
    }
    rho = RhoDMatrix(PDT::Spin1Half);
  }
}

void SpinorBarWaveFunction::
constructSpinInfo(const vector<LorentzSpinorBar<SqrtEnergy> > & waves,
		  tPPtr part,Direction dir, bool time) {
  assert(waves.size()==2);
  tFermionSpinPtr inspin = !part->spinInfo() ? tFermionSpinPtr() :
    dynamic_ptr_cast<tFermionSpinPtr>(part->spinInfo());
  if(inspin) {
    for(unsigned int ix=0;ix<2;++ix) {
      if(( dir == outgoing &&  time) || 
	 ( dir == incoming && !time))
	inspin->setBasisState(ix,waves[ix].bar());
      else
	inspin->setDecayState(ix,waves[ix].bar());
    }
  }
  else {
    FermionSpinPtr temp = new_ptr(FermionSpinInfo(part->momentum(),time));
    part->spinInfo(temp);
    for(unsigned int ix=0;ix<2;++ix) {
      if(( dir == outgoing &&  time) || 
	 ( dir == incoming && !time))
	temp->setBasisState(ix,waves[ix].bar());
      else
	temp->setDecayState(ix,waves[ix].bar());
    }
  }
}

void SpinorBarWaveFunction::
constructSpinInfo(const vector<SpinorBarWaveFunction> & waves,
		  tPPtr part,Direction dir, bool time) {
  assert(waves.size()==2);
  tFermionSpinPtr inspin = !part->spinInfo() ? tFermionSpinPtr() :
    dynamic_ptr_cast<tFermionSpinPtr>(part->spinInfo());
  if(inspin) {
    for(unsigned int ix=0;ix<2;++ix)
      if (dir==outgoing) inspin->setBasisState(ix,waves[ix].dimensionedWf().bar());
      else               inspin->setDecayState(ix,waves[ix].dimensionedWf().bar());
  }
  else {
    FermionSpinPtr temp = new_ptr(FermionSpinInfo(part->momentum(),time));
    part->spinInfo(temp);
    for(unsigned int ix=0;ix<2;++ix) {
      if(( dir == outgoing &&  time) || 
	 ( dir == incoming && !time))
	temp->setBasisState(ix,waves[ix].dimensionedWf().bar());
      else
	temp->setDecayState(ix,waves[ix].dimensionedWf().bar());
    }
  }
}
