// -*- C++ -*-
//
// VectorWaveFunction.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2017 Peter Richardson, Leif Lonnblad
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

using namespace ThePEG;
using namespace ThePEG::Helicity;

// calculate the Wavefunction
void VectorWaveFunction::calculateWaveFunction(unsigned int ihel,VectorPhase vphase) {
  Direction dir=direction();
  if(dir==intermediate)
    throw ThePEG::Helicity::HelicityConsistencyError() 
      << "In VectorWaveFunction::calcluateWaveFunction "
      << "particle must be incoming or outgoing not intermediate" 
      << Exception::abortnow;
  // check a valid helicity combination
  if(ihel==0 || ihel==2||(ihel==1&&mass()>ZERO)) {
    int jhel=ihel-1;
    // extract the momentum components
    double fact=-1.; if(dir==incoming){fact=1.;}
    Energy ppx=fact*px(),ppy=fact*py(),ppz=fact*pz(),pee=fact*e(),pmm=mass();
    // calculate some kinematic quantites;
    Energy2 pt2 = ppx*ppx+ppy*ppy;
    Energy pabs = sqrt(pt2+ppz*ppz);
    Energy pt = sqrt(pt2);
    // overall phase of the vector
    Complex phase;
    if(vphase==vector_phase) {
      if(pt==ZERO || ihel==1) phase = 1.;
      else if(ihel==0)            phase = Complex(ppx/pt,-fact*ppy/pt);
      else                        phase = Complex(ppx/pt, fact*ppy/pt);
    }
    else                          phase = 1.;
    if(ihel!=1) phase = phase/sqrt(2.);
    // first the +/-1 helicity states
    if(ihel!=1) {
      // first the no pt case
      if(pt==ZERO) {
	double sgnz;
	sgnz = ppz<ZERO ? -1. : 1.;
	_wf.setX(-complex<double>(jhel)*phase);
	_wf.setY(sgnz*phase*complex<double>(0,-fact));
	_wf.setZ(0.);
	_wf.setT(0.);
      }
      else {
	InvEnergy opabs=1./pabs;
	InvEnergy opt  =1./pt;
	_wf.setX(phase*complex<double>(-jhel*ppz*ppx*opabs*opt, fact*ppy*opt));
	_wf.setY(phase*complex<double>(-jhel*ppz*ppy*opabs*opt,-fact*ppx*opt));
	_wf.setZ(double(jhel*pt*opabs)*phase);
	_wf.setT(0.);
      }
    }
    // 0 component for massive vectors
    else {
      if(pabs==ZERO) {
	_wf.setX(0.);
	_wf.setY(0.);
	_wf.setZ(1.);
	_wf.setT(0.);
      }
      else {
	InvEnergy empabs=pee/pmm/pabs;
	_wf.setX(double(empabs*ppx));
	_wf.setY(double(empabs*ppy));
	_wf.setZ(double(empabs*ppz));
	_wf.setT(double(pabs/pmm));
      }
    }
  }
  // special return the momentum as a check of gauge invariance
  else if(ihel==10) {
    _wf.setX(double(px()/MeV));
    _wf.setY(double(py()/MeV));
    _wf.setZ(double(pz()/MeV));
    _wf.setT(double(e()/MeV));
  }
  // issue warning and return zero
  else {
    ThePEG::Helicity::HelicityConsistencyError() 
      << "Invalid Helicity = " << ihel << " requested for Vector " 
      << particle()->PDGName() << Exception::abortnow;
    _wf.setX(0.);_wf.setY(0.);_wf.setZ(0.);_wf.setT(0.);
  }
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
