// -*- C++ -*-
//
// LorentzRSSpinor.tcc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2017 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LorentzRSSpinor class.
//
// Author: Peter Richardson
//

#include "LorentzRSSpinor.h"
#include "LorentzRSSpinorBar.h"

using namespace ThePEG;
using namespace ThePEG::Helicity;

// return the barred spinor
template<typename Value> LorentzRSSpinorBar<Value>
LorentzRSSpinor<Value>::bar() const {
  complex<Value> out[4][4];
  unsigned int ix;
  // HELAS
  for(ix=0;ix<4;++ix) {
    out[ix][0] = conj(_spin[ix][2]);
    out[ix][1] = conj(_spin[ix][3]);
    out[ix][2] = conj(_spin[ix][0]);
    out[ix][3] = conj(_spin[ix][1]);
  }
  return LorentzRSSpinorBar<Value>(out[0][0],out[0][1],out[0][2],out[0][3],
				   out[1][0],out[1][1],out[1][2],out[1][3],
				   out[2][0],out[2][1],out[2][2],out[2][3],
				   out[3][0],out[3][1],out[3][2],out[3][3],_type);
}

// boost the spinor
template<typename Value> LorentzRSSpinor<Value> & 
LorentzRSSpinor<Value>::boost(double bx,double by,double bz) {
  // work out beta and chi
  double b2(bx*bx+by*by+bz*bz),beta(sqrt(b2)),chi(atanh(beta));
  double sinhchi(sinh(0.5*chi)/beta),coshchi(cosh(0.5*chi));
  double gamma = 1.0/sqrt(1.0-b2);
  double gmmone = b2 >0 ? (gamma-1.)/b2 : 0.0;
  double bvec[3]={bx,by,bz};
  unsigned int ix,iy,ixa,iya;
  // vector boost matrix
  double boostv[4][4];
  for(ix=0;ix<3;++ix)
    {
      for(iy=0;iy<3;++iy){boostv[ix][iy]=bvec[ix]*bvec[iy]*gmmone;}
      boostv[ix][ix]+=1;
      boostv[ix][3]=gamma*bvec[ix];
      boostv[3][ix]=boostv[ix][3];
    }
  boostv[3][3]=gamma;
  // spinor boost matrix
  Complex boosts[4][4],ii(0.,1.),nxminy(bx-ii*by),nxpiny(bx+ii*by);
  boosts[0][0] = coshchi-sinhchi*bz;
  boosts[0][1] = -sinhchi*nxminy;
  boosts[0][2] = 0.;
  boosts[0][3] = 0.;
  boosts[1][0] = -sinhchi*nxpiny;
  boosts[1][1] = coshchi+sinhchi*bz;
  boosts[1][2] = 0.;
  boosts[1][3] = 0.;
  boosts[2][0] = 0.;
  boosts[2][1] = 0.;
  boosts[2][2] = coshchi+sinhchi*bz;
  boosts[2][3] = +sinhchi*nxminy;
  boosts[3][0] = 0.;
  boosts[3][1] = 0.;
  boosts[3][2] = +sinhchi*nxpiny;
  boosts[3][3] = coshchi-sinhchi*bz;
  Complex out[4][4];
  // apply the boost
  for(ix=0;ix<4;++ix)
    {
      for(iy=0;iy<4;++iy)
	{
	  out[ix][iy]=0.;
	  for(ixa=0;ixa<4;++ixa)
	    {
	      for(iya=0;iya<4;++iya)
		{out[ix][iy]+=boostv[ix][ixa]*boosts[iy][iya]*_spin[ixa][iya];}
	    }
	}
    }
  *this=LorentzRSSpinor<Value>(out[0][0],out[0][1],out[0][2],out[0][3],
			out[1][0],out[1][1],out[1][2],out[1][3],
			out[2][0],out[2][1],out[2][2],out[2][3],
			out[3][0],out[3][1],out[3][2],out[3][3],_type);
  return *this;
}

// boost the spinor
template<typename Value> LorentzRSSpinor<Value> & LorentzRSSpinor<Value>::boost(const Boost & boostvec)
{
  const double beta = boostvec.mag(),b2=beta*beta;
  const double bx=boostvec.x(),by=boostvec.y(),bz=boostvec.z();
  const double boostvec1[3] = {bx, by, bz};
  const double gamma = 1.0/sqrt(1.0-b2);
  const double gmmone = b2 >0 ? (gamma-1.)/b2 : 0.0;
  const double chi = atanh(beta);
  const double sinhchi = sinh(0.5*chi)/beta, coshchi = cosh(0.5*chi);
  complex<Value> out[4][4];
  Complex ii(0.,1.);
  const Complex nxminy=bx-ii*by;
  const Complex nxpiny=bx+ii*by;
  unsigned int ix,iy,ixa,iya;
  // vector boost matrix
  double boostv[4][4];
  for(ix=0;ix<3;++ix)
    {
      for(iy=0;iy<3;++iy){boostv[ix][iy]=boostvec1[ix]*boostvec1[iy]*gmmone;}
      boostv[ix][ix]+=1;
      boostv[ix][3]=gamma*boostvec1[ix];
      boostv[3][ix]=boostv[ix][3];
    }
  boostv[3][3]=gamma;
  // spinor boost matrix
  Complex boosts[4][4];
  boosts[0][0] = coshchi-sinhchi*bz;
  boosts[0][1] = -sinhchi*nxminy;
  boosts[0][2] = 0.;
  boosts[0][3] = 0.;
  boosts[1][0] = -sinhchi*nxpiny;
  boosts[1][1] = coshchi+sinhchi*bz;
  boosts[1][2] = 0.;
  boosts[1][3] = 0.;
  boosts[2][0] = 0.;
  boosts[2][1] = 0.;
  boosts[2][2] = coshchi+sinhchi*bz;
  boosts[2][3] = +sinhchi*nxminy;
  boosts[3][0] = 0.;
  boosts[3][1] = 0.;
  boosts[3][2] = +sinhchi*nxpiny;
  boosts[3][3] = coshchi-sinhchi*bz;
  // apply the boost
  for(ix=0;ix<4;++ix)
    {
      for(iy=0;iy<4;++iy)
	{
	  out[ix][iy]=complex<Value>();
	  for(ixa=0;ixa<4;++ixa)
	    {
	      for(iya=0;iya<4;++iya)
		{out[ix][iy]+=boostv[ix][ixa]*boosts[iy][iya]*_spin[ixa][iya];}
	    }
	}
    }
  *this= LorentzRSSpinor<Value>(out[0][0],out[0][1],out[0][2],out[0][3],
			 out[1][0],out[1][1],out[1][2],out[1][3],
			 out[2][0],out[2][1],out[2][2],out[2][3],
			 out[3][0],out[3][1],out[3][2],out[3][3],_type);
  return *this;
}

//general lorentz tranformation
template<typename Value> LorentzRSSpinor<Value> & LorentzRSSpinor<Value>::transform(const LorentzRotation & r)
{
  unsigned int ix,iy,ixa,iya;
  LorentzRSSpinor<Value> out;
  for(ix=0;ix<4;++ix)
    {
      for(iy=0;iy<4;++iy)
	{
	  out(ix,iy)=complex<Value>();
	  for(ixa=0;ixa<4;++ixa)
	    {
	      for(iya=0;iya<4;++iya)
		{out(ix,iy)+=r.one()(ix,ixa)*r.half()(iy,iya)*_spin[ixa][iya];}
	    }
	}
    }
  *this=out;
  return *this;
}
