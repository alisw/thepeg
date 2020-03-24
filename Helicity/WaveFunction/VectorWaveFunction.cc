// -*- C++ -*-
//
// VectorWaveFunction.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2019 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VectorWaveFunction class.
//
// Author: Peter Richardson
//

#include "VectorWaveFunction.h"
#include "ThePEG/Helicity/HelicityFunctions.h"

using namespace ThePEG;
using namespace ThePEG::Helicity;

// calculate the Wavefunction
void VectorWaveFunction::calculateWaveFunction(unsigned int ihel,VectorPhase vphase) {
  _wf = HelicityFunctions::polarizationVector(momentum(),ihel,direction(),vphase);
}


void VectorWaveFunction::
calculateWaveFunctions(vector<LorentzPolarizationVector> & waves,
		       tPPtr particle,Direction dir,bool massless,
		       VectorPhase phase) {
  tVectorSpinPtr inspin = !particle->spinInfo() ? tVectorSpinPtr() :
    dynamic_ptr_cast<tVectorSpinPtr>(particle->spinInfo());
  waves.resize(3);
  if(inspin) {
    if(dir==outgoing) {
      for(unsigned int ix=0;ix<3;++ix)
	waves[ix]=inspin->getProductionBasisState(ix);
    }
    else {
      inspin->decay();
      for(unsigned int ix=0;ix<3;++ix)
	waves[ix]=inspin->getDecayBasisState(ix);
    }
  }
  else {
    assert(!particle->spinInfo());
    VectorWaveFunction wave(particle->momentum(),particle->dataPtr(),0,
			    dir,phase);
    for(unsigned int ix=0;ix<3;++ix) {
      if(massless&&ix==1) {
	waves[ix] = LorentzPolarizationVector();
      }
      else {
	if(ix!=0) wave.reset(ix);
	waves[ix] = wave.wave();
      }
    }
  }
}

void VectorWaveFunction::
calculateWaveFunctions(vector<VectorWaveFunction> & waves,
		       tPPtr particle, Direction dir, bool massless,
		       VectorPhase phase) {
  tVectorSpinPtr inspin = !particle->spinInfo() ? tVectorSpinPtr() :
    dynamic_ptr_cast<tVectorSpinPtr>(particle->spinInfo());
  waves.resize(3);
  if(inspin) {
    if(dir==outgoing) {
      for(unsigned int ix=0;ix<3;++ix)
	waves[ix]=VectorWaveFunction(particle->momentum(),
				     particle->dataPtr(),
				     inspin->getProductionBasisState(ix),dir);
    }
    else {
      inspin->decay();
      for(unsigned int ix=0;ix<3;++ix)
	waves[ix]=VectorWaveFunction(particle->momentum(),
				     particle->dataPtr(),
				     inspin->getDecayBasisState(ix),dir);
    }
  }
  else {
    assert(!particle->spinInfo());
    calculateWaveFunctions(waves,particle->momentum(),particle->dataPtr(),
			   dir,massless,phase);
  }
}

void VectorWaveFunction::
calculateWaveFunctions(vector<VectorWaveFunction> & waves,
		       const Lorentz5Momentum & momentum, tcPDPtr parton, 
		       Direction dir, bool massless,
		       VectorPhase phase) {
  waves.resize(3);
  VectorWaveFunction wave(momentum,parton,0,dir,phase);
  for(unsigned int ix=0;ix<3;++ix) {
    if(massless&&ix==1) {
      waves[ix] = VectorWaveFunction(momentum,parton,dir);
    }
    else {
      if(ix!=0) wave.reset(ix);
      waves[ix] = wave;
    }
  }
}

void VectorWaveFunction::
calculateWaveFunctions(vector<LorentzPolarizationVector> & waves,
		       RhoDMatrix & rho,tPPtr particle,
		       Direction dir,bool massless,
		       VectorPhase phase) {
  tVectorSpinPtr inspin = !particle->spinInfo() ? tVectorSpinPtr() :
    dynamic_ptr_cast<tVectorSpinPtr>(particle->spinInfo());
  waves.resize(3);
  if(inspin) {
    if(dir==outgoing) {
      for(unsigned int ix=0;ix<3;++ix)
	waves[ix]=inspin->getProductionBasisState(ix);
      rho = RhoDMatrix(PDT::Spin1);
    }
    else {
      inspin->decay();
      for(unsigned int ix=0;ix<3;++ix)
	waves[ix]=inspin->getDecayBasisState(ix);
      rho = inspin->rhoMatrix();
    }
  }
  else {
    assert(!particle->spinInfo());
    VectorWaveFunction wave(particle->momentum(),particle->dataPtr(),0,
			    dir,phase);
    for(unsigned int ix=0;ix<3;++ix) {
      if(massless&&ix==1) {
	waves[ix] = LorentzPolarizationVector();
      }
      else {
	if(ix!=0) wave.reset(ix);
	waves[ix] = wave.wave();
      }
    }
    rho = RhoDMatrix(PDT::Spin1);
  }
}

void VectorWaveFunction::
calculateWaveFunctions(vector<VectorWaveFunction> & waves,
		       RhoDMatrix & rho,
		       tPPtr particle, Direction dir,bool massless,
		       VectorPhase phase) {
  tVectorSpinPtr inspin = !particle->spinInfo() ? tVectorSpinPtr() :
    dynamic_ptr_cast<tVectorSpinPtr>(particle->spinInfo());
  waves.resize(3);
  if(inspin) {
    if(dir==outgoing) {
      for(unsigned int ix=0;ix<3;++ix)
	waves[ix]=VectorWaveFunction(particle->momentum(),
				     particle->dataPtr(),
				     inspin->getProductionBasisState(ix),dir);
      rho = RhoDMatrix(PDT::Spin1);
    }
    else {
      inspin->decay();
      for(unsigned int ix=0;ix<3;++ix)
	waves[ix]=VectorWaveFunction(particle->momentum(),
				     particle->dataPtr(),
				     inspin->getDecayBasisState(ix),dir);
      rho = inspin->rhoMatrix();
    }
  }
  else {
    assert(!particle->spinInfo());
    VectorWaveFunction wave(particle->momentum(),particle->dataPtr(),0,
			    dir,phase);
    for(unsigned int ix=0;ix<3;++ix) {
      if(massless&&ix==1) {
	waves[ix] = VectorWaveFunction(particle->momentum(),
				       particle->dataPtr(),dir);
      }
      else {
	if(ix!=0) wave.reset(ix);
	waves[ix] = wave;
      }
    }
    rho = RhoDMatrix(PDT::Spin1);
  }
}

void VectorWaveFunction::
constructSpinInfo(const vector<LorentzPolarizationVector> & waves,
		  tPPtr part,Direction dir, bool time,bool ) {
  assert(waves.size()==3);
  tVectorSpinPtr inspin = !part->spinInfo() ? tVectorSpinPtr() :
    dynamic_ptr_cast<tVectorSpinPtr>(part->spinInfo());
  if(inspin) {
    for(unsigned int ix=0;ix<3;++ix)
      if(( dir == outgoing &&  time) ||
	 ( dir == incoming && !time))
	inspin->setBasisState(ix,waves[ix]);
      else
	inspin->setDecayState(ix,waves[ix]);
  }
  else {
    VectorSpinPtr temp = new_ptr(VectorSpinInfo(part->momentum(),time));
    part->spinInfo(temp);
    for(unsigned int ix=0;ix<3;++ix)
      if(( dir == outgoing &&  time) ||
	 ( dir == incoming && !time))
	temp->setBasisState(ix,waves[ix]);
      else
	temp->setDecayState(ix,waves[ix]);
  }
}

void VectorWaveFunction::
constructSpinInfo(const vector<VectorWaveFunction> & waves,
		  tPPtr part,Direction dir, bool time,bool ) {
  assert(waves.size()==3);
  tVectorSpinPtr inspin = !part->spinInfo() ? tVectorSpinPtr() :
    dynamic_ptr_cast<tVectorSpinPtr>(part->spinInfo());
  if(inspin) {
    for(unsigned int ix=0;ix<3;++ix)
      if(( dir == outgoing &&  time) ||
	 ( dir == incoming && !time))
	inspin->setBasisState(ix,waves[ix].wave());
      else
	inspin->setDecayState(ix,waves[ix].wave());
  }
  else {
    VectorSpinPtr temp = new_ptr(VectorSpinInfo(part->momentum(),time));
    part->spinInfo(temp);
    for(unsigned int ix=0;ix<3;++ix)
      if(( dir == outgoing &&  time) ||
	 ( dir == incoming && !time))
	temp->setBasisState(ix,waves[ix].wave());
      else
	temp->setDecayState(ix,waves[ix].wave());
  }
}
