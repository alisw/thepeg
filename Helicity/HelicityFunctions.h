// -*- C++ -*-
//
// HelicityFunctions.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2019 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_HelicityFunctions_H
#define ThePEG_HelicityFunctions_H
//
// This is the declaration of the HelicityFunctions header for common functions
// used in helicity calculations to avoid duplication of code

#include "ThePEG/Vectors/LorentzVector.h"
#include "LorentzSpinor.h"
#include "LorentzSpinorBar.h"

namespace ThePEG {
namespace Helicity {
namespace HelicityFunctions {

inline LorentzPolarizationVector polarizationVector(const Lorentz5Momentum & p,
						    unsigned int ihel,
						    Direction dir,
						    VectorPhase vphase=default_vector_phase) {
  // check the direction
  assert(dir!=intermediate);
  // special helicity combination for gauge invariance tests
  if(ihel==10) return p*UnitRemoval::InvE;
  // check a valid helicity combination
  assert(ihel==0 || ihel == 2 || ((ihel==1 || ihel==3) && p.mass()>ZERO ));
  // convert the helicitty from 0,1,2 to -1,0,1
  int jhel=ihel-1;
  // extract the momentum components
  double fact = dir==outgoing ? -1 : 1;
  Energy ppx=fact*p.x(),ppy=fact*p.y(),ppz=fact*p.z(),pee=fact*p.e(),pmm=p.mass();
  // calculate some kinematic quantites
  Energy2 pt2 = ppx*ppx+ppy*ppy;
  Energy pabs = sqrt(pt2+ppz*ppz);
  Energy pt = sqrt(pt2);
  // zero subtracted
  if(ihel==3) {
    InvEnergy pre = pmm/pabs/(pee+pabs);
    return LorentzPolarizationVector(double(pre*ppx),double(pre*ppy),double(pre*ppz),-double(pre*pabs));
  }
  // overall phase of the vector
  Complex phase(1.);
  if(vphase==vector_phase) {
    if(pt==ZERO || ihel==1) phase = 1.;
    else if(ihel==0)        phase = Complex(ppx/pt,-fact*ppy/pt);
    else                    phase = Complex(ppx/pt, fact*ppy/pt);
  }
  if(ihel!=1) phase = phase/sqrt(2.);
  // first the +/-1 helicity states
  if(ihel!=1) {
    // first the zero pt case
    if(pt==ZERO) {
      double sgnz = ppz<ZERO ? -1. : 1.;
      return LorentzPolarizationVector(-complex<double>(jhel)*phase,
				       sgnz*phase*complex<double>(0,-fact),
				       0.,0.);
    }
    else {
      InvEnergy opabs=1./pabs;
      InvEnergy opt  =1./pt;
      return LorentzPolarizationVector(phase*complex<double>(-jhel*ppz*ppx*opabs*opt, fact*ppy*opt),
				       phase*complex<double>(-jhel*ppz*ppy*opabs*opt,-fact*ppx*opt),
				       double(jhel*pt*opabs)*phase,0.);
    }
  }
  // 0 component for massive vectors
  else {
    if(pabs==ZERO) {
      return LorentzPolarizationVector(0.,0.,1.,0.);
    }
    else {
      InvEnergy empabs=pee/pmm/pabs;
      return LorentzPolarizationVector(double(empabs*ppx),double(empabs*ppy),
				       double(empabs*ppz),double(pabs/pmm));
    }
  }
}


inline LorentzSpinor<SqrtEnergy> dimensionedSpinor(const Lorentz5Momentum & p,
						   unsigned int ihel,
						   Direction dir) {
  // check direction and helicity
  assert(dir!=intermediate);
  assert(ihel<=1);
  // extract the momentum components
  double fact = dir==incoming ? 1 : -1.;
  Energy ppx=fact*p.x(),ppy=fact*p.y(),ppz=fact*p.z(),pee=fact*p.e(),pmm=p.mass();
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
  return LorentzSpinor<SqrtEnergy>(upper*hel_wf[0],upper*hel_wf[1],
				   lower*hel_wf[0],lower*hel_wf[1],
				   (dir==incoming) ? SpinorType::u : SpinorType::v);
}

inline LorentzSpinor<double> spinor(const Lorentz5Momentum & p,
				    unsigned int ihel,
				    Direction dir) {
  LorentzSpinor<SqrtEnergy> temp = dimensionedSpinor(p,ihel,dir);
  return LorentzSpinor<double>(temp.s1()*UnitRemoval::InvSqrtE,
			       temp.s2()*UnitRemoval::InvSqrtE,
			       temp.s3()*UnitRemoval::InvSqrtE,
			       temp.s4()*UnitRemoval::InvSqrtE,temp.Type());
}

inline LorentzSpinorBar<SqrtEnergy> dimensionedSpinorBar(const Lorentz5Momentum & p,
							 unsigned int ihel,
							 Direction dir) {
  // check direction and helicity
  assert(dir!=intermediate);
  assert(ihel<=1);
  // extract the momentum components
  double fact = dir==incoming ? 1. : -1.;
  Energy ppx=fact*p.x(),ppy=fact*p.y(),ppz=fact*p.z(),pee=fact*p.e(),pmm=p.mass();
  // define and calculate some kinematic quantities
  Energy2 ptran2  = ppx*ppx+ppy*ppy;
  Energy pabs   = sqrt(ptran2+ppz*ppz);
  Energy ptran  = sqrt(ptran2);
  // first need to evalulate the 2-component helicity spinors
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
      hel_wf[1] = Complex(denominator/rtppluspz*complex<Energy>(ppx,-ppy));
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
      hel_wf[0] = Complex(denominator/rtppluspz*complex<Energy>(-ppx,-ppy));
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
  return LorentzSpinorBar<SqrtEnergy>(upper*hel_wf[0],
				      upper*hel_wf[1],
				      lower*hel_wf[0],
				      lower*hel_wf[1],
				      (dir==incoming) ? SpinorType::v : SpinorType::u);
}

inline LorentzSpinorBar<double> spinorBar(const Lorentz5Momentum & p,
					  unsigned int ihel,
					  Direction dir) {
  LorentzSpinorBar<SqrtEnergy> temp = dimensionedSpinorBar(p,ihel,dir);
  return LorentzSpinorBar<double>(temp.s1()*UnitRemoval::InvSqrtE,
  				  temp.s2()*UnitRemoval::InvSqrtE,
  				  temp.s3()*UnitRemoval::InvSqrtE,
  				  temp.s4()*UnitRemoval::InvSqrtE,temp.Type());
}
}
}
}

#endif /* ThePEG_HelicityFunctions_H */
