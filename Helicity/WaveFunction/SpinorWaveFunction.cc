// -*- C++ -*-
//
// SpinorWaveFunction.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2019 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SpinorWaveFunction class.
//
// Author: Peter Richardson
//

#include "SpinorWaveFunction.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "SpinorBarWaveFunction.h"
#include "ThePEG/Helicity/HelicityFunctions.h"

using namespace ThePEG;
using namespace Helicity;

// calculate the Wavefunction
void SpinorWaveFunction::calculateWaveFunction(unsigned int ihel) {
  _wf = HelicityFunctions::spinor(momentum(),ihel,direction());
}

SpinorBarWaveFunction SpinorWaveFunction::bar() {
  Lorentz5Momentum p = momentum();
  if(direction()==outgoing) p *= -1.;
  tcPDPtr ptemp = particle();
  if(direction()==incoming&&particle()->CC())
    ptemp = particle()->CC();
  return SpinorBarWaveFunction(p,ptemp,_wf.bar(),direction());
}	  

void SpinorWaveFunction::
calculateWaveFunctions(vector<LorentzSpinor<SqrtEnergy> > & waves,
		       tPPtr particle,Direction dir) {
  tFermionSpinPtr inspin = !particle->spinInfo() ? tFermionSpinPtr() :
    dynamic_ptr_cast<tFermionSpinPtr>(particle->spinInfo());
  waves.resize(2);
  // spin info object exists
  if(inspin) {
    if(dir==outgoing) {
      waves[0] = inspin->getProductionBasisState(0);
      waves[1] = inspin->getProductionBasisState(1);
    }
    else {
      inspin->decay();
      if( (particle->id()>0&&inspin->getDecayBasisState(0).Type()!=SpinorType::u) || 
	  (particle->id()<0&&inspin->getDecayBasisState(0).Type()!=SpinorType::v)) {
	waves[0] = inspin->getDecayBasisState(0).conjugate();
	waves[1] = inspin->getDecayBasisState(1).conjugate();
      }
      else {
	waves[0] = inspin->getDecayBasisState(0);
	waves[1] = inspin->getDecayBasisState(1);
      }
    }
  }
  // do the calculation
  else {
    assert(!particle->spinInfo());
    SpinorWaveFunction wave(particle->momentum(),particle->dataPtr(),dir);
    for(unsigned int ix=0;ix<2;++ix) {
      wave.reset(ix);
      waves[ix] = wave.dimensionedWave();
    }
  }
}

void SpinorWaveFunction::
calculateWaveFunctions(vector<SpinorWaveFunction> & waves,
		       tPPtr particle,Direction dir) {
  tFermionSpinPtr inspin = !particle->spinInfo() ? tFermionSpinPtr() :
    dynamic_ptr_cast<tFermionSpinPtr>(particle->spinInfo());
  waves.resize(2);
  // spin info object exists
  if(inspin) {
    if(dir==outgoing) {
      for(unsigned int ix=0;ix<2;++ix)
	waves[ix] = SpinorWaveFunction(particle,
				       inspin->getProductionBasisState(ix),dir);
    }
    else {
      inspin->decay();
      if((particle->id()>0&&inspin->getDecayBasisState(0).Type()!=SpinorType::u) ||
	 (particle->id()<0&&inspin->getDecayBasisState(0).Type()!=SpinorType::v)) {
	for(unsigned int ix=0;ix<2;++ix)
	  waves[ix] = SpinorWaveFunction(particle,
					 inspin->getDecayBasisState(ix).conjugate(),dir);
      }
      else {
	for(unsigned int ix=0;ix<2;++ix)
	  waves[ix] = SpinorWaveFunction(particle,
					 inspin->getDecayBasisState(ix),dir);
      }
    }
  }
  // do the calculation
  else {
    assert(!particle->spinInfo());
    calculateWaveFunctions(waves,particle->momentum(),particle->dataPtr(),dir);
  }
}
void SpinorWaveFunction::
calculateWaveFunctions(vector<SpinorWaveFunction> & waves,
		       const Lorentz5Momentum & momentum,
		       tcPDPtr parton,Direction dir) {
  waves.resize(2);
  SpinorWaveFunction wave(momentum,parton,dir);
  for(unsigned int ix=0;ix<2;++ix) {
    wave.reset(ix);
    waves[ix] = wave;
  }
}

void SpinorWaveFunction::
calculateWaveFunctions(vector<LorentzSpinor<SqrtEnergy> > & waves,
		       RhoDMatrix & rho,
		       tPPtr particle,Direction dir) {
  tFermionSpinPtr inspin = !particle->spinInfo() ? tFermionSpinPtr() :
    dynamic_ptr_cast<tFermionSpinPtr>(particle->spinInfo());
  waves.resize(2);
  // spin info object exists
  if(inspin) {
    if(dir==outgoing) {
      waves[0] = inspin->getProductionBasisState(0);
      waves[1] = inspin->getProductionBasisState(1);
      rho = RhoDMatrix(PDT::Spin1Half);
    }
    else {
      inspin->decay();
      if((particle->id()>0&&inspin->getDecayBasisState(0).Type()!=SpinorType::u) ||
	 (particle->id()<0&&inspin->getDecayBasisState(0).Type()!=SpinorType::v)) {
	waves[0] = inspin->getDecayBasisState(0).conjugate();
	waves[1] = inspin->getDecayBasisState(1).conjugate();
      }
      else {
	waves[0] = inspin->getDecayBasisState(0);
	waves[1] = inspin->getDecayBasisState(1);
      }
      rho = inspin->rhoMatrix();
    }
  }
  // do the calculation
  else {
    assert(!particle->spinInfo());
    SpinorWaveFunction wave(particle->momentum(),particle->dataPtr(),dir);
    for(unsigned int ix=0;ix<2;++ix) {
      wave.reset(ix);
      waves[ix] = wave.dimensionedWave();
    }
    rho = RhoDMatrix(PDT::Spin1Half);
  }
}

void SpinorWaveFunction::
calculateWaveFunctions(vector<SpinorWaveFunction> & waves,
		       RhoDMatrix & rho,
		       tPPtr particle,Direction dir) {
  tFermionSpinPtr inspin = !particle->spinInfo() ? tFermionSpinPtr() :
    dynamic_ptr_cast<tFermionSpinPtr>(particle->spinInfo());
  waves.resize(2);
  // spin info object exists
  if(inspin) {
    if(dir==outgoing) {
      for(unsigned int ix=0;ix<2;++ix)
	waves[ix] = SpinorWaveFunction(particle,
				       inspin->getProductionBasisState(ix),dir);
      rho = RhoDMatrix(PDT::Spin1Half);
    }
    else {
      inspin->decay();
      if((particle->id()>0&&inspin->getDecayBasisState(0).Type()!=SpinorType::u) ||
	 (particle->id()<0&&inspin->getDecayBasisState(0).Type()!=SpinorType::v)) {
	for(unsigned int ix=0;ix<2;++ix)
	  waves[ix] = SpinorWaveFunction(particle,
					 inspin->getDecayBasisState(ix).conjugate(),dir);
      }
      else {
	for(unsigned int ix=0;ix<2;++ix)
	  waves[ix] = SpinorWaveFunction(particle,
					 inspin->getDecayBasisState(ix),dir);
      }
      rho = inspin->rhoMatrix();
    }
  }
  // do the calculation
  else {
    assert(!particle->spinInfo());
    SpinorWaveFunction wave(particle->momentum(),particle->dataPtr(),dir);
    for(unsigned int ix=0;ix<2;++ix) {
      wave.reset(ix);
      waves[ix] = wave;
    }
    rho = RhoDMatrix(PDT::Spin1Half);
  }
}

void SpinorWaveFunction::
constructSpinInfo(const vector<LorentzSpinor<SqrtEnergy> > & waves,
		  tPPtr particle,Direction dir,bool time) {
  assert(waves.size()==2);
  tFermionSpinPtr inspin = !particle->spinInfo() ? tFermionSpinPtr() :
    dynamic_ptr_cast<tFermionSpinPtr>(particle->spinInfo());
  if(inspin) {
    for(unsigned int ix=0;ix<2;++ix) {
      if(( dir == outgoing &&  time) || 
	 ( dir == incoming && !time))
	inspin->setBasisState(ix,waves[ix]);
      else
	inspin->setDecayState(ix,waves[ix]);
    }
  }
  else {
    FermionSpinPtr temp = new_ptr(FermionSpinInfo(particle->momentum(),time));
    particle->spinInfo(temp);
    for(unsigned int ix=0;ix<2;++ix) {
      if(( dir == outgoing &&  time) || 
	 ( dir == incoming && !time))
	temp->setBasisState(ix,waves[ix]);
      else
	temp->setDecayState(ix,waves[ix]);
    }
  }
}

void SpinorWaveFunction::
constructSpinInfo(const vector<SpinorWaveFunction> & waves,
		  tPPtr particle,Direction dir,bool time) {
  assert(waves.size()==2);
  tFermionSpinPtr inspin = !particle->spinInfo() ? tFermionSpinPtr() :
    dynamic_ptr_cast<tFermionSpinPtr>(particle->spinInfo());
  if(inspin) {
    for(unsigned int ix=0;ix<2;++ix) {
      if(( dir == outgoing &&  time) || 
	 ( dir == incoming && !time)) 
	inspin->setBasisState(ix,waves[ix].dimensionedWf());
      else
	inspin->setDecayState(ix,waves[ix].dimensionedWf());
    }
  }
  else {
    FermionSpinPtr temp = new_ptr(FermionSpinInfo(particle->momentum(),time));
    particle->spinInfo(temp);
    for(unsigned int ix=0;ix<2;++ix) {
      if(( dir == outgoing &&  time) || 
	 ( dir == incoming && !time))
	temp->setBasisState(ix,waves[ix].dimensionedWf());
      else
	temp->setDecayState(ix,waves[ix].dimensionedWf());
    }
  }
}
