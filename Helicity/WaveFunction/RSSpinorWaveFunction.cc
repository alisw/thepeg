// -*- C++ -*-
//
// RSSpinorWaveFunction.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2019 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RSSpinorWaveFunction class.
//

#include "RSSpinorWaveFunction.h"
#include "ThePEG/Helicity/HelicityFunctions.h"

using namespace ThePEG;
using namespace ThePEG::Helicity;

// calculate the Wavefunction
void RSSpinorWaveFunction::calculateWaveFunction(unsigned int ihel) {
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
  _wf = LorentzRSSpinor<double>(direction()==outgoing ? SpinorType::v : SpinorType::u);
  // compute the spinors
  std::array<LorentzSpinor<double>,2> spin;
  if(ihel!=0) spin[0] = HelicityFunctions::spinor(ptemp,1,direction());
  if(ihel!=3) spin[1] = HelicityFunctions::spinor(ptemp,0,direction());
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

void RSSpinorWaveFunction::
calculateWaveFunctions(vector<LorentzRSSpinor<SqrtEnergy> > & waves,
		       tPPtr particle,Direction dir) {
  tRSFermionSpinPtr inspin = !particle->spinInfo() ? tRSFermionSpinPtr() :
    dynamic_ptr_cast<tRSFermionSpinPtr>(particle->spinInfo());
  waves.resize(4);
  // spin info object exists
  if(inspin) {
    if(dir==outgoing) {
      for(unsigned int ix=0;ix<4;++ix)
	waves[ix] = inspin->getProductionBasisState(ix);
    }
    else {
      inspin->decay();
      for(unsigned int ix=0;ix<4;++ix)
	waves[ix] = inspin->getDecayBasisState(ix);
    }
  }
  // do the calculation
  else {
    assert(!particle->spinInfo());
    RSSpinorWaveFunction wave(particle->momentum(),particle->dataPtr(),dir);
    for(unsigned int ix=0;ix<4;++ix) {
      wave.reset(ix);
      waves[ix] = wave.dimensionedWf();
    }
  }
}

void RSSpinorWaveFunction::
calculateWaveFunctions(vector<RSSpinorWaveFunction> & waves,
		       tPPtr particle,Direction dir) {
  tRSFermionSpinPtr inspin = !particle->spinInfo() ? tRSFermionSpinPtr() :
    dynamic_ptr_cast<tRSFermionSpinPtr>(particle->spinInfo());
  waves.resize(4);
  // spin info object exists
  if(inspin) {
    if(dir==outgoing) {
      for(unsigned int ix=0;ix<4;++ix)
	waves[ix] = RSSpinorWaveFunction(particle,
					 inspin->getProductionBasisState(ix),dir);
    }
    else {
      inspin->decay();
      for(unsigned int ix=0;ix<4;++ix)
	waves[ix] = RSSpinorWaveFunction(particle,
					 inspin->getDecayBasisState(ix),dir);
    }
  }
  // do the calculation
  else {
    assert(!particle->spinInfo());
    calculateWaveFunctions(waves,particle->momentum(),particle->dataPtr(),dir);
  }
}

void RSSpinorWaveFunction::
calculateWaveFunctions(vector<RSSpinorWaveFunction> & waves,
		       const Lorentz5Momentum & momentum,
		       tcPDPtr parton,Direction dir) {
  waves.resize(4);
  RSSpinorWaveFunction wave(momentum,parton,dir);
  for(unsigned int ix=0;ix<4;++ix) {
    wave.reset(ix);
    waves[ix] = wave;
  }
}

void RSSpinorWaveFunction::
calculateWaveFunctions(vector<LorentzRSSpinor<SqrtEnergy> > & waves,
		       RhoDMatrix & rho,
		       tPPtr particle,Direction dir) {
  tRSFermionSpinPtr inspin = !particle->spinInfo() ? tRSFermionSpinPtr() :
    dynamic_ptr_cast<tRSFermionSpinPtr>(particle->spinInfo());
  waves.resize(4);
  // spin info object exists
  if(inspin) {
    if(dir==outgoing) {
      for(unsigned int ix=0;ix<4;++ix)
	waves[ix] = inspin->getProductionBasisState(ix);
      rho = RhoDMatrix(PDT::Spin3Half);
    }
    else {
      inspin->decay();
      for(unsigned int ix=0;ix<4;++ix)
	waves[ix] = inspin->getDecayBasisState(ix);
      rho = inspin->rhoMatrix();
    }
  }
  // do the calculation
  else {
    assert(!particle->spinInfo());
    RSSpinorWaveFunction wave(particle->momentum(),particle->dataPtr(),dir);
    for(unsigned int ix=0;ix<4;++ix) {
      wave.reset(ix);
      waves[ix] = wave.dimensionedWf();
    }
    rho = RhoDMatrix(PDT::Spin3Half);
  }
}

void RSSpinorWaveFunction::
calculateWaveFunctions(vector<RSSpinorWaveFunction> & waves,
		       RhoDMatrix & rho,
		       tPPtr particle,Direction dir) {
  tRSFermionSpinPtr inspin = !particle->spinInfo() ? tRSFermionSpinPtr() :
    dynamic_ptr_cast<tRSFermionSpinPtr>(particle->spinInfo());
  waves.resize(4);
  // spin info object exists
  if(inspin) {
    if(dir==outgoing) {
      for(unsigned int ix=0;ix<4;++ix)
	waves[ix] = RSSpinorWaveFunction(particle,
					 inspin->getProductionBasisState(ix),dir);
      rho = RhoDMatrix(PDT::Spin3Half);
    }
    else {
      inspin->decay();
      for(unsigned int ix=0;ix<4;++ix)
	waves[ix] = RSSpinorWaveFunction(particle,
					 inspin->getDecayBasisState(ix),dir);
      rho = inspin->rhoMatrix();
    }
  }
  // do the calculation
  else {
    assert(!particle->spinInfo());
    RSSpinorWaveFunction wave(particle->momentum(),particle->dataPtr(),dir);
    for(unsigned int ix=0;ix<4;++ix) {
      wave.reset(ix);
      waves[ix] = wave;
    }
    rho = RhoDMatrix(PDT::Spin3Half);
  }
}

void RSSpinorWaveFunction::
constructSpinInfo(const vector<LorentzRSSpinor<SqrtEnergy> > & waves,
		  tPPtr particle,Direction dir,bool time) {
  assert(waves.size()==4);
  tRSFermionSpinPtr inspin = !particle->spinInfo() ? tRSFermionSpinPtr() :
    dynamic_ptr_cast<tRSFermionSpinPtr>(particle->spinInfo());
  if(inspin) {
    for(unsigned int ix=0;ix<4;++ix)
      if(dir==outgoing) inspin->setBasisState(ix,waves[ix]);
      else              inspin->setDecayState(ix,waves[ix]);
  }
  else {
    RSFermionSpinPtr temp = new_ptr(RSFermionSpinInfo(particle->momentum(),time));
    particle->spinInfo(temp);
    for(unsigned int ix=0;ix<4;++ix)
      if(dir==outgoing) temp->setBasisState(ix,waves[ix]);
      else              temp->setDecayState(ix,waves[ix]);
  }
}

void RSSpinorWaveFunction::
constructSpinInfo(const vector<RSSpinorWaveFunction> & waves,
		  tPPtr particle,Direction dir,bool time) {
  assert(waves.size()==4);
  tRSFermionSpinPtr inspin = !particle->spinInfo() ? tRSFermionSpinPtr() :
    dynamic_ptr_cast<tRSFermionSpinPtr>(particle->spinInfo());
  if(inspin) {
    for(unsigned int ix=0;ix<4;++ix)
      if(dir==outgoing) inspin->setBasisState(ix,waves[ix].dimensionedWf());
      else              inspin->setDecayState(ix,waves[ix].dimensionedWf());
  }
  else {
    RSFermionSpinPtr temp = new_ptr(RSFermionSpinInfo(particle->momentum(),time));
    particle->spinInfo(temp);
    for(unsigned int ix=0;ix<4;++ix)
      if(dir==outgoing) temp->setBasisState(ix,waves[ix].dimensionedWf());
      else              temp->setDecayState(ix,waves[ix].dimensionedWf());
  }
}
