// -*- C++ -*-
//
// SpinorWaveFunction.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2017 Peter Richardson, Leif Lonnblad
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

using namespace ThePEG;
using namespace Helicity;

// calculate the Wavefunction
void SpinorWaveFunction::calculateWaveFunction(unsigned int ihel) {
  // check helicity is O.K.
  Direction dir = direction();
  if(dir==intermediate) throw ThePEG::Helicity::HelicityConsistencyError() 
    << "In SpinorWaveFunction::calcluateWaveFunction "
    << "particle must be incoming or outgoing not intermediate" 
    << Exception::abortnow;
  if(ihel>1) throw ThePEG::Helicity::HelicityConsistencyError() 
    << "Invalid Helicity = " << ihel << " requested for Spinor" 
    << Exception::abortnow;
  // extract the momentum components
  double fact=-1.; if(dir==incoming){fact=1.;}
  Energy ppx=fact*px(),ppy=fact*py(),ppz=fact*pz(),pee=fact*e(),pmm=mass();
  // define and calculate some kinematic quantities
  Energy2 ptran2  = ppx*ppx+ppy*ppy;
  Energy pabs   = sqrt(ptran2+ppz*ppz);
  Energy ptran  = sqrt(ptran2);
  // first need to evalulate the 2-component helicity spinors 
  // this is the same regardless of which definition of the spinors
  // we are using
  complex <double> hel_wf[2];
  // compute the + spinor for + helicty particles and - helicity antiparticles
  if((dir==incoming && ihel== 1) || (dir==outgoing && ihel==0)) {
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
      hel_wf[1] = Complex(denominator/rtppluspz*complex<Energy>(ppx,ppy));
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
      hel_wf[0] = Complex(denominator/rtppluspz*complex<Energy>(-ppx,ppy));
      hel_wf[1] = denominator*rtppluspz;
    }
  }

  SqrtEnergy upper,lower;
  SqrtEnergy eplusp  = sqrt(max(pee+pabs,ZERO));
  SqrtEnergy eminusp = ( pmm != ZERO ) ? pmm/eplusp : ZERO;
  // set up the coefficients for the different cases
  if(dir==incoming) {
    if(ihel==1) {
      upper = eminusp;
      lower = eplusp;
    }
    else {
      upper = eplusp;
      lower = eminusp;
    }
  }
  else {
    if(ihel==1) {
      upper = -eplusp;
      lower = eminusp;
    }
    else {
      upper = eminusp;
      lower =-eplusp;
    }
  }
  // now finally we can construct the spinors
  _wf = LorentzSpinor<double>( (dir==incoming) ? SpinorType::u : SpinorType::v);
  _wf[0] = Complex(upper*hel_wf[0]*UnitRemoval::InvSqrtE);
  _wf[1] = Complex(upper*hel_wf[1]*UnitRemoval::InvSqrtE);
  _wf[2] = Complex(lower*hel_wf[0]*UnitRemoval::InvSqrtE);
  _wf[3] = Complex(lower*hel_wf[1]*UnitRemoval::InvSqrtE);
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
    SpinorWaveFunction wave(particle->momentum(),particle->dataPtr(),dir);
    for(unsigned int ix=0;ix<2;++ix) {
      wave.reset(ix);
      waves[ix] = wave;
    }
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
