// -*- C++ -*-
//
// RSSpinorWaveFunction.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2011 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RSSpinorWaveFunction class.
//

#include "RSSpinorWaveFunction.h"

using namespace ThePEG;
using namespace ThePEG::Helicity;

// calculate the Wavefunction
void RSSpinorWaveFunction::calculateWaveFunction(unsigned int ihel) {
  LorentzRSSpinor<double> news;
  if(direction()==incoming) 
    news=LorentzRSSpinor<double>(u_spinortype);
  else 
    news=LorentzRSSpinor<double>(v_spinortype);
  unsigned int ix,iy;
  // check helicity and type
  assert(direction()!=intermediate);
  assert(ihel<=3);
  // massive
  // only two valid helicities in massless case
  assert( mass()>ZERO || ( ihel == 0 || ihel==3 ) );
  // extract the momentum components
  // compute the normal spinors to construct the RS spinor
  Complex hel_wf[2][2];
  if(direction()==incoming) {
    // the + spinor
    hel_wf[0][0] = 1.;
    hel_wf[0][1] = 0.;
    // the - spinor
    hel_wf[1][0] = 0.;
    hel_wf[1][1] = 1.;
  }
  else {
    // the + spinor
    hel_wf[0][0] = 0.;
    hel_wf[0][1] = 1.;
    // the - spinor
    hel_wf[1][0] = 1.;
    hel_wf[1][1] = 0.;
  }
  // prefactors
  double fact = direction()==incoming ? 1. : -1.;
  Energy pmm=mass(),pee=fact*e();
  Energy pabs   = sqrt(sqr(px())+sqr(py())+sqr(pz()));
  SqrtEnergy eplusp = sqrt(pee+pabs);
  SqrtEnergy eminusp = ( pmm == ZERO ) ? ZERO : pmm/eplusp;
  SqrtEnergy upper[2],lower[2];
  if(direction()==incoming) {
    upper[0] = eminusp;
    lower[0] = eplusp ;
    upper[1] = eplusp ;
    lower[1] = eminusp;
  }
  else {
    upper[0] =-eplusp ;
    lower[0] = eminusp;
    upper[1] = eminusp;
    lower[1] =-eplusp ;
  }
  // now construct the spinors
  complex<SqrtEnergy> spinor[2][4];
  for(ix=0;ix<2;++ix) {
    spinor[ix][0] = upper[ix]*hel_wf[ix][0];
    spinor[ix][1] = upper[ix]*hel_wf[ix][1];
    spinor[ix][2] = lower[ix]*hel_wf[ix][0];
    spinor[ix][3] = lower[ix]*hel_wf[ix][1];
  }
  // compute the polarization vectors to construct the RS spinor
  Complex vec[3][4],ii(0.,1.);
  double ort = sqrt(0.5);
  double r1 = ( pmm == ZERO ) ? 0. : double(pee /pmm);
  double r2 = ( pmm == ZERO ) ? 0. : double(pabs/pmm);
  if(direction()==incoming) {
    vec[0][0] =-ort;
    vec[0][1] =-ort*ii;
    vec[0][2] = 0.;
    vec[0][3] = 0.;
    vec[1][0] = 0.;
    vec[1][1] = 0.;
    vec[1][2] = r1;
    vec[1][3] = r2;
    vec[2][0] = ort;
    vec[2][1] =-ort*ii;
    vec[2][2] = 0.;
    vec[2][3] = 0.;
  }
  else {
    vec[0][0] = ort;
    vec[0][1] =-ort*ii;
    vec[0][2] = 0.;
    vec[0][3] = 0.;
    vec[1][0] = 0.;
    vec[1][1] = 0.;
    vec[1][2] =-r1;
    vec[1][3] =-r2;
    vec[2][0] =-ort;
    vec[2][1] =-ort*ii;
    vec[2][2] = 0.;
    vec[2][3] = 0.;
  }
  // now we can put the bits together to compute the RS spinor
  double or3(sqrt(1./3.)),tor3(sqrt(2./3.));
  if(ihel==3) {
    for(ix=0;ix<4;++ix)
      for(iy=0;iy<4;++iy)
	news(ix,iy)=UnitRemoval::InvSqrtE*vec[0][ix]*spinor[0][iy];
  }
  else if(ihel==2) {
    for(ix=0;ix<4;++ix)
      for(iy=0;iy<4;++iy)
	news(ix,iy)=UnitRemoval::InvSqrtE*
	  (or3*vec[0][ix]*spinor[1][iy]+tor3*vec[1][ix]*spinor[0][iy]);
  }
  else if(ihel==1) {
    for(ix=0;ix<4;++ix)
      for(iy=0;iy<4;++iy)
	news(ix,iy)=UnitRemoval::InvSqrtE*
	  (or3*vec[2][ix]*spinor[0][iy]+tor3*vec[1][ix]*spinor[1][iy]);
  }
  else if(ihel==0) {
    for(ix=0;ix<4;++ix)
      for(iy=0;iy<4;++iy)
	news(ix,iy)=UnitRemoval::InvSqrtE*vec[2][ix]*spinor[1][iy];
  }
  // spinor is currently along the z axis, rotate so in right direction
  if(pabs/pmm>1e-8) {
    Axis axis;
    axis.setX(fact*momentum().x()/pabs);
    axis.setY(fact*momentum().y()/pabs);
    axis.setZ(fact*momentum().z()/pabs);
    double sinth(sqrt(sqr(axis.x())+sqr(axis.y())));
    if(sinth>1e-8) {
      LorentzRotation rot;
      rot.setRotate(acos(axis.z()),Axis(-axis.y()/sinth,axis.x()/sinth,0.));
      _wf= news.transform(rot); 
    }
    else if (axis.z()<0.) {
      LorentzRotation rot;
      rot.setRotateX(Constants::pi);
      _wf= news.transform(rot);
    }
    else {
      _wf = news;
    }    
  }
  else {
    _wf=news;
  }
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
    RSSpinorWaveFunction wave(particle->momentum(),particle->dataPtr(),dir);
    for(unsigned int ix=0;ix<4;++ix) {
      wave.reset(ix);
      waves[ix] = wave;
    }
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
