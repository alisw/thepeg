// -*- C++ -*-
//
// SpinorBarWaveFunction.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2019 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
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
#include "ThePEG/Helicity/HelicityFunctions.h"

using namespace ThePEG;
  using namespace Helicity;

// calculate the Wavefunction
void SpinorBarWaveFunction::calculateWaveFunction(unsigned int ihel) {
  _wf = HelicityFunctions::spinorBar(momentum(),ihel,direction());
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
      if( (particle->id()>0&&inspin->getDecayBasisState(0).Type()!=SpinorType::u) || 
	  (particle->id()<0&&inspin->getDecayBasisState(0).Type()!=SpinorType::v)) {
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
      if((particle->id()>0&&inspin->getDecayBasisState(0).Type()!=SpinorType::u) ||
	 (particle->id()<0&&inspin->getDecayBasisState(0).Type()!=SpinorType::v)) {
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
    calculateWaveFunctions(waves,particle->momentum(),particle->dataPtr(),dir);
  }
}

void SpinorBarWaveFunction::calculateWaveFunctions(vector<SpinorBarWaveFunction> & waves,
						   const Lorentz5Momentum & momentum,
						   tcPDPtr parton,Direction dir) {
  waves.resize(2);
  SpinorBarWaveFunction wave(momentum,parton,dir);
  for(unsigned int ix=0;ix<2;++ix) {
    wave.reset(ix);
    waves[ix] = wave;
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
      if((particle->id()>0&&inspin->getDecayBasisState(0).Type()!=SpinorType::u) ||
	 (particle->id()<0&&inspin->getDecayBasisState(0).Type()!=SpinorType::v)) {
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
      if((particle->id()>0&&inspin->getDecayBasisState(0).Type()!=SpinorType::u) ||
	 (particle->id()<0&&inspin->getDecayBasisState(0).Type()!=SpinorType::v)) {
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
      if(( dir == outgoing &&  time) ||
	 ( dir == incoming && !time))
	inspin->setBasisState(ix,waves[ix].dimensionedWf().bar());
      else
	inspin->setDecayState(ix,waves[ix].dimensionedWf().bar());
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
