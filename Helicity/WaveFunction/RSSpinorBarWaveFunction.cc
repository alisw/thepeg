// -*- C++ -*-
//
// RSSpinorBarWaveFunction.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2019 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RSSpinorBarWaveFunction class.
//

#include "RSSpinorBarWaveFunction.h"
#include "ThePEG/Helicity/HelicityFunctions.h"

using namespace ThePEG;
using namespace ThePEG::Helicity;

// calculate the Wavefunction
void RSSpinorBarWaveFunction::calculateWaveFunction(unsigned int ihel) {
  // if rest frame, make sure speical case is used
  Lorentz5Momentum ptemp=momentum();
  double pr = ptemp.vect().mag2()/sqr(momentum().e());
  if(pr<1e-20) {
    ptemp.setX(ZERO);
    ptemp.setY(ZERO);
    ptemp.setZ(ZERO);
  }
  assert(direction()!=intermediate);
  assert(ihel<=3);
  // only two valid helicities in massless case
  assert( mass()>ZERO || (ihel == 0 || ihel == 3 ) );
  // new direct calculation
  _wf = LorentzRSSpinorBar<double>(direction()==outgoing ? SpinorType::u : SpinorType::v);
  // compute the spinors
  std::array<LorentzSpinorBar<double>,2> spin;
  if(ihel!=0) spin[0] = HelicityFunctions::spinorBar(ptemp,1,direction());
  if(ihel!=3) spin[1] = HelicityFunctions::spinorBar(ptemp,0,direction());
  // compute the polarization vectors to construct the RS spinor
  std::array<LorentzPolarizationVector,3> eps;
  if(ihel>=2)
    eps[0] = HelicityFunctions::polarizationVector(ptemp,2,direction(),vector_phase);
  else
    eps[2] = HelicityFunctions::polarizationVector(ptemp,0,direction(),vector_phase);
  if(mass()!=ZERO&&ihel!=0&&ihel!=3)
    eps[1] =  HelicityFunctions::polarizationVector(ptemp,1,direction());
  
  // now we can put the bits together to compute the RS spinor
  double or3(sqrt(1./3.)),tor3(sqrt(2./3.));
  if(ihel==3) {
    for(unsigned int iy=0;iy<4;++iy) {
      _wf(0,iy) = eps[0].x()*spin[0][iy];
      _wf(1,iy) = eps[0].y()*spin[0][iy];
      _wf(2,iy) = eps[0].z()*spin[0][iy];
      _wf(3,iy) = eps[0].t()*spin[0][iy];
    }
  }
  else if(ihel==2) {
    for(unsigned int iy=0;iy<4;++iy) {
      _wf(0,iy) = or3*eps[0].x()*spin[1][iy]+tor3*eps[1].x()*spin[0][iy];
      _wf(1,iy) = or3*eps[0].y()*spin[1][iy]+tor3*eps[1].y()*spin[0][iy];
      _wf(2,iy) = or3*eps[0].z()*spin[1][iy]+tor3*eps[1].z()*spin[0][iy];
      _wf(3,iy) = or3*eps[0].t()*spin[1][iy]+tor3*eps[1].t()*spin[0][iy];
    }
  }
  else if(ihel==1) {
    for(unsigned int iy=0;iy<4;++iy) {
      _wf(0,iy) = or3*eps[2].x()*spin[0][iy]+tor3*eps[1].x()*spin[1][iy];
      _wf(1,iy) = or3*eps[2].y()*spin[0][iy]+tor3*eps[1].y()*spin[1][iy];
      _wf(2,iy) = or3*eps[2].z()*spin[0][iy]+tor3*eps[1].z()*spin[1][iy];
      _wf(3,iy) = or3*eps[2].t()*spin[0][iy]+tor3*eps[1].t()*spin[1][iy];
    }
  }
  else if(ihel==0) {
    for(unsigned int iy=0;iy<4;++iy) {
      _wf(0,iy) = eps[2].x()*spin[1][iy];
      _wf(1,iy) = eps[2].y()*spin[1][iy];
      _wf(2,iy) = eps[2].z()*spin[1][iy];
      _wf(3,iy) = eps[2].t()*spin[1][iy];
    }
  }
  // this makes the phase choice the same as madgraph, useful for debugging only
  // Energy pt = ptemp.perp();
  // double fact = direction()==incoming ? 1. : -1.;
  // Complex emphi(fact*ptemp.x()/pt,-fact*ptemp.y()/pt);
  // if(ihel==3) {
  //   for(unsigned int ix=0;ix<4;++ix)
  //     for(unsigned int iy=0;iy<4;++iy)
  // 	_wf(ix,iy) /= emphi;
  // }
  // else if(ihel==1) {
  //   for(unsigned int ix=0;ix<4;++ix)
  //     for(unsigned int iy=0;iy<4;++iy)
  // 	_wf(ix,iy) *= emphi;
  // }
  // else if(ihel==0) {
  //   for(unsigned int ix=0;ix<4;++ix)
  //     for(unsigned int iy=0;iy<4;++iy)
  // 	_wf(ix,iy) *= sqr(emphi);
  // }
}

void RSSpinorBarWaveFunction::
calculateWaveFunctions(vector<LorentzRSSpinorBar<SqrtEnergy> > & waves,
		       tPPtr particle,Direction dir) {
  tRSFermionSpinPtr inspin = !particle->spinInfo() ? tRSFermionSpinPtr() :
    dynamic_ptr_cast<tRSFermionSpinPtr>(particle->spinInfo());
  waves.resize(4);
  // spin info object exists
  if(inspin) {
    if(dir==outgoing) {
      for(unsigned int ix=0;ix<4;++ix)
	waves[ix] = inspin->getProductionBasisState(ix).bar();
    }
    else {
      inspin->decay();
      for(unsigned int ix=0;ix<4;++ix)
	waves[ix] = inspin->getDecayBasisState(ix).bar();
    }
  }
  // do the calculation
  else {
    assert(!particle->spinInfo());
    RSSpinorBarWaveFunction wave(particle->momentum(),particle->dataPtr(),dir);
    for(unsigned int ix=0;ix<4;++ix) {
      wave.reset(ix);
      waves[ix] = wave.dimensionedWf();
    }
  }
}

void RSSpinorBarWaveFunction::
calculateWaveFunctions(vector<RSSpinorBarWaveFunction> & waves,
		       tPPtr particle,Direction dir) {
  tRSFermionSpinPtr inspin = !particle->spinInfo() ? tRSFermionSpinPtr() :
    dynamic_ptr_cast<tRSFermionSpinPtr>(particle->spinInfo());
  waves.resize(4);
  // spin info object exists
  if(inspin) {
    if(dir==outgoing) {
      for(unsigned int ix=0;ix<4;++ix)
	waves[ix] = RSSpinorBarWaveFunction(particle,
					    inspin->getProductionBasisState(ix).bar(),
					    dir);
    }
    else {
      inspin->decay();
      for(unsigned int ix=0;ix<4;++ix)
	waves[ix] = RSSpinorBarWaveFunction(particle,
					    inspin->getDecayBasisState(ix).bar(),
					    dir);
    }
  }
  // do the calculation
  else {
    assert(!particle->spinInfo());
    calculateWaveFunctions(waves,particle->momentum(),particle->dataPtr(),dir);
    }
}

void RSSpinorBarWaveFunction::
calculateWaveFunctions(vector<RSSpinorBarWaveFunction> & waves,
		       const Lorentz5Momentum & momentum,
		       tcPDPtr parton,Direction dir) {
  waves.resize(4);
  RSSpinorBarWaveFunction wave(momentum,parton,dir);
  for(unsigned int ix=0;ix<4;++ix) {
    wave.reset(ix);
    waves[ix] = wave;
  }
}

void RSSpinorBarWaveFunction::
calculateWaveFunctions(vector<LorentzRSSpinorBar<SqrtEnergy> > & waves,
		       RhoDMatrix & rho,
		       tPPtr particle,Direction dir) {
  tRSFermionSpinPtr inspin = !particle->spinInfo() ? tRSFermionSpinPtr() :
    dynamic_ptr_cast<tRSFermionSpinPtr>(particle->spinInfo());
  waves.resize(4);
  // spin info object exists
  if(inspin) {
    if(dir==outgoing) {
      for(unsigned int ix=0;ix<4;++ix)
	waves[ix] = inspin->getProductionBasisState(ix).bar();
      rho = RhoDMatrix(PDT::Spin3Half);
    }
    else {
      inspin->decay();
      for(unsigned int ix=0;ix<4;++ix)
	waves[ix] = inspin->getDecayBasisState(ix).bar();
      rho = inspin->rhoMatrix();
    }
  }
  // do the calculation
  else {
    assert(!particle->spinInfo());
    RSSpinorBarWaveFunction wave(particle->momentum(),particle->dataPtr(),dir);
    for(unsigned int ix=0;ix<4;++ix) {
      wave.reset(ix);
      waves[ix] = wave.dimensionedWf();
    }
    rho = RhoDMatrix(PDT::Spin3Half);
  }
}

void RSSpinorBarWaveFunction::
calculateWaveFunctions(vector<RSSpinorBarWaveFunction> & waves,
		       RhoDMatrix & rho,
		       tPPtr particle,Direction dir) {
  tRSFermionSpinPtr inspin = !particle->spinInfo() ? tRSFermionSpinPtr() :
    dynamic_ptr_cast<tRSFermionSpinPtr>(particle->spinInfo());
  waves.resize(4);
  // spin info object exists
  if(inspin) {
    if(dir==outgoing) {
      for(unsigned int ix=0;ix<4;++ix)
	waves[ix] = RSSpinorBarWaveFunction(particle,
					    inspin->getProductionBasisState(ix).bar(),
					    dir);
      rho = RhoDMatrix(PDT::Spin3Half);
    }
    else {
      inspin->decay();
      for(unsigned int ix=0;ix<4;++ix)
	waves[ix] = RSSpinorBarWaveFunction(particle,
					    inspin->getDecayBasisState(ix).bar(),
					    dir);
      rho = inspin->rhoMatrix();
    }
  }
  // do the calculation
  else {
    assert(!particle->spinInfo());
    RSSpinorBarWaveFunction wave(particle->momentum(),particle->dataPtr(),dir);
    for(unsigned int ix=0;ix<4;++ix) {
      wave.reset(ix);
      waves[ix] = wave;
    }
    rho = RhoDMatrix(PDT::Spin3Half);
  }
}

void RSSpinorBarWaveFunction::
constructSpinInfo(const vector<LorentzRSSpinorBar<SqrtEnergy> > & waves,
		  tPPtr particle,Direction dir, bool time) {
  assert(waves.size()==4);
  tRSFermionSpinPtr inspin = !particle->spinInfo() ? tRSFermionSpinPtr() :
    dynamic_ptr_cast<tRSFermionSpinPtr>(particle->spinInfo());
  if(inspin) {
    for(unsigned int ix=0;ix<4;++ix)
      if(dir==outgoing) inspin->setBasisState(ix,waves[ix].bar());
      else              inspin->setDecayState(ix,waves[ix].bar());
  }
  else {
    RSFermionSpinPtr temp = new_ptr(RSFermionSpinInfo(particle->momentum(),time));
    particle->spinInfo(temp);
    for(unsigned int ix=0;ix<4;++ix)
      if(dir==outgoing) temp->setBasisState(ix,waves[ix].bar());
      else              temp->setDecayState(ix,waves[ix].bar());
  }
}

void RSSpinorBarWaveFunction::
constructSpinInfo(const vector<RSSpinorBarWaveFunction> & waves,
		  tPPtr part,Direction dir, bool time) {
  assert(waves.size()==4);
  tRSFermionSpinPtr inspin = !part->spinInfo() ? tRSFermionSpinPtr() :
    dynamic_ptr_cast<tRSFermionSpinPtr>(part->spinInfo());
  if(inspin) {
    for(unsigned int ix=0;ix<4;++ix)
      if (dir==outgoing) inspin->setBasisState(ix,waves[ix].dimensionedWf().bar());
      else               inspin->setDecayState(ix,waves[ix].dimensionedWf().bar());
  }
  else {
    RSFermionSpinPtr temp = new_ptr(RSFermionSpinInfo(part->momentum(),time));
    part->spinInfo(temp);
    for(unsigned int ix=0;ix<4;++ix)
      if(dir==outgoing) temp->setBasisState(ix,waves[ix].dimensionedWf().bar());
      else              temp->setDecayState(ix,waves[ix].dimensionedWf().bar());
  }
}
