// -*- C++ -*-
//
// LorentzSpinor.tcc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2019 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined member
// functions of the LorentzSpinor class.
//
// Author: Peter Richardson
//

#include "LorentzSpinor.h"
#include "LorentzSpinorBar.h"

using namespace ThePEG;
using namespace ThePEG::Helicity;

// return the barred spinor
template <typename Value>
LorentzSpinorBar<Value> LorentzSpinor<Value>::bar() const
{
  complex<Value> output[4];
  // HELAS
  output[0] = conj(_spin[2]);
  output[1] = conj(_spin[3]);
  output[2] = conj(_spin[0]);
  output[3] = conj(_spin[1]);
  return LorentzSpinorBar<Value>(output[0],output[1],output[2],output[3],_type);
}

// boost the spinor
template <typename Value>
LorentzSpinor<Value> & LorentzSpinor<Value>::boost(double bx,double by,double bz)
{
  // work out beta and chi
  double beta=sqrt(bx*bx+by*by+bz*bz);
  double chi = atanh(beta);
  double sinhchi = sinh(0.5*chi)/beta, coshchi = cosh(0.5*chi);
  // calculate the new spinor
  complex<Value> out[4];
  Complex ii(0.,1.);
  Complex nxminy=bx-ii*by;
  Complex nxpiny=bx+ii*by;
  out[0] = coshchi*_spin[0]+sinhchi*(-bz*_spin[0]-nxminy*_spin[1]);
  out[1] = coshchi*_spin[1]+sinhchi*(+bz*_spin[1]-nxpiny*_spin[0]);
  out[2] = coshchi*_spin[2]+sinhchi*(+bz*_spin[2]+nxminy*_spin[3]);
  out[3] = coshchi*_spin[3]+sinhchi*(-bz*_spin[3]+nxpiny*_spin[2]);
  for(unsigned int ix=0;ix<4;++ix) _spin[ix]=out[ix];
  return *this;
}

// boost the spinor
template <typename Value>
LorentzSpinor<Value> & LorentzSpinor<Value>::boost(const Boost & boostv)
{
  double beta = boostv.mag();
  double bx=boostv.x(),by=boostv.y(),bz=boostv.z();
  double chi = atanh(beta);
  double sinhchi = sinh(0.5*chi)/beta, coshchi = cosh(0.5*chi);
  complex<Value> out[4];
  Complex ii(0.,1.);
  Complex nxminy=bx-ii*by;
  Complex nxpiny=bx+ii*by;
  out[0] = coshchi*_spin[0]+sinhchi*(-bz*_spin[0]-nxminy*_spin[1]);
  out[1] = coshchi*_spin[1]+sinhchi*(+bz*_spin[1]-nxpiny*_spin[0]);
  out[2] = coshchi*_spin[2]+sinhchi*(+bz*_spin[2]+nxminy*_spin[3]);
  out[3] = coshchi*_spin[3]+sinhchi*(-bz*_spin[3]+nxpiny*_spin[2]);
  for(unsigned int ix=0;ix<4;++ix) _spin[ix]=out[ix];
  return *this;
}

// general transform
template <typename Value>
LorentzSpinor<Value> & LorentzSpinor<Value>::
transform(const SpinHalfLorentzRotation & r) {
  unsigned int ix,iy;
  complex<Value> out[4];
  for(ix=0;ix<4;++ix) {
    out[ix]=complex<Value>();
    for(iy=0;iy<4;++iy) out[ix]+=r(ix,iy)*_spin[iy];
  }
  for(ix=0;ix<4;++ix) _spin[ix]=out[ix];
  return *this;
}


// conjugation
template <typename Value>
LorentzSpinor<Value> LorentzSpinor<Value>::conjugate() const {
  SpinorType new_type;
  switch(_type) {
  case SpinorType::u:
    new_type=SpinorType::v;
    break;
  case SpinorType::v:
    new_type=SpinorType::u;
    break;
  case SpinorType::unknown:
  default:
    new_type=SpinorType::unknown;
    break;
  }
  return LorentzSpinor<Value>( conj(_spin[3]),
			      -conj(_spin[2]),
			      -conj(_spin[1]),
			      +conj(_spin[0]),
			       new_type);
}
