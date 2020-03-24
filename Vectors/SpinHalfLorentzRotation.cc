// -*- C++ -*-
//
// SpinHalfLorentzRotation.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SpinHalfLorentzRotation class.
//

#include "SpinHalfLorentzRotation.h"

using namespace ThePEG;

// default constructor
SpinHalfLorentzRotation::SpinHalfLorentzRotation() 
: _mx{{ {{1., 0., 0., 0.}},
        {{0., 1., 0., 0.}},
        {{0., 0., 1., 0.}},
        {{0., 0., 0., 1.}} }}
{}

// constructor giving the components of a Lorentz boost
SpinHalfLorentzRotation::
SpinHalfLorentzRotation(double bx, double by, double bz, double gamma) 
: _mx{}
{
  setBoost(bx,by,bz,gamma);
}

// constructor with boost vector
SpinHalfLorentzRotation::
SpinHalfLorentzRotation (const Boost & b, double gamma)
: _mx{}
{
  setBoost(b,gamma);
}

// protected set all elements constructor
SpinHalfLorentzRotation::
SpinHalfLorentzRotation(Complex c1c1,Complex c1c2,Complex c1c3,Complex c1c4,
			Complex c2c1,Complex c2c2,Complex c2c3,Complex c2c4,
			Complex c3c1,Complex c3c2,Complex c3c3,Complex c3c4,
			Complex c4c1,Complex c4c2,Complex c4c3,Complex c4c4) 
: _mx{{ {{c1c1, c1c2, c1c3, c1c4}},
        {{c2c1, c2c2, c2c3, c2c4}},
        {{c3c1, c3c2, c3c3, c3c4}},
        {{c4c1, c4c2, c4c3, c4c4}} }}
{}

// check for identity matrix
bool SpinHalfLorentzRotation::isIdentity() const {
  return _mx == MatrixT{{ {{1., 0., 0., 0.}},
                          {{0., 1., 0., 0.}},
                          {{0., 0., 1., 0.}},
                          {{0., 0., 0., 1.}} }};
}

// inverse ( inverse is gamma0 S dagger gamma 0)
SpinHalfLorentzRotation SpinHalfLorentzRotation::inverse() const {
  return SpinHalfLorentzRotation
    {conj(_mx[2][2]),conj(_mx[3][2]),conj(_mx[0][2]),conj(_mx[1][2]),
     conj(_mx[2][3]),conj(_mx[3][3]),conj(_mx[0][3]),conj(_mx[1][3]),
     conj(_mx[2][0]),conj(_mx[3][0]),conj(_mx[0][0]),conj(_mx[1][0]),
     conj(_mx[2][1]),conj(_mx[3][1]),conj(_mx[0][1]),conj(_mx[1][1])};
}

// specify the components of a lorentz boost
SpinHalfLorentzRotation & SpinHalfLorentzRotation::
setBoost (double bx, double by, double bz, double) {
  // work out beta and chi
  const double eps=1e-10;
  double beta(sqrt(bx*bx+by*by+bz*bz));
  double chi(atanh(beta)), chc(cosh(0.5*chi)), shc(0.5);
  if ( beta > eps ) 
    shc=sinh(0.5*chi)/beta;
  Complex ii(0.,1.),nxminy(bx-ii*by),nxplny(bx+ii*by);
  _mx[0][0]= chc-shc*bz;
  _mx[0][1]=-shc*nxminy;
  _mx[0][2]= 0.        ;
  _mx[0][3]= 0.        ;

  _mx[1][0]=-shc*nxplny;
  _mx[1][1]= chc+shc*bz;
  _mx[1][2]= 0.        ;
  _mx[1][3]= 0.        ;

  _mx[2][0]= 0.        ;
  _mx[2][1]= 0.        ;
  _mx[2][2]= chc+shc*bz;
  _mx[2][3]=+shc*nxminy;
  
  _mx[3][0]= 0.        ;
  _mx[3][1]= 0.        ;
  _mx[3][2]=+shc*nxplny;
  _mx[3][3]= chc-shc*bz;
  return *this;
}

// specify a boost vector
SpinHalfLorentzRotation & SpinHalfLorentzRotation::setBoost (const Boost & b,double) {
  // work out beta and chi
  const double eps=1e-10;
  double bx(b.x()),by(b.y()),bz(b.z()),beta(b.mag()),chi(atanh(beta)),
    chc(cosh(0.5*chi)),shc(0.5);
  if(beta>eps){shc=sinh(0.5*chi)/beta;}
  Complex ii(0.,1.),nxminy(bx-ii*by),nxplny(bx+ii*by);
  _mx[0][0]= chc-shc*bz;_mx[0][1]=-shc*nxminy;_mx[0][2]= 0.        ;_mx[0][3]= 0.        ;
  _mx[1][0]=-shc*nxplny;_mx[1][1]= chc+shc*bz;_mx[1][2]= 0.        ;_mx[1][3]= 0.        ;
  _mx[2][0]= 0.        ;_mx[2][1]= 0.        ;_mx[2][2]= chc+shc*bz;_mx[2][3]=+shc*nxminy;
  _mx[3][0]= 0.        ;_mx[3][1]= 0.        ;_mx[3][2]=+shc*nxplny;_mx[3][3]= chc-shc*bz;
  return *this;
}

// lorentz boost in x direction
SpinHalfLorentzRotation & SpinHalfLorentzRotation::setBoostX (double & bx) {
  // work out beta and chi
  double chi(atanh(bx)),shc(sinh(0.5*chi)),chc(cosh(0.5*chi));
  _mx[0][0]= chc;_mx[0][1]=-shc;_mx[0][2]= 0. ;_mx[0][3]= 0. ;
  _mx[1][0]=-shc;_mx[1][1]= chc;_mx[1][2]= 0. ;_mx[1][3]= 0. ;
  _mx[2][0]= 0  ;_mx[2][1]= 0. ;_mx[2][2]= chc;_mx[2][3]=+shc;
  _mx[3][0]= 0  ;_mx[3][1]= 0. ;_mx[3][2]=+shc;_mx[3][3]= chc;
  return *this;
}

// lorentz boost in y direction
SpinHalfLorentzRotation & SpinHalfLorentzRotation::setBoostY (double & by)
{
  // work out beta and chi
  double chi(atanh(by)),chc(cosh(0.5*chi));
  Complex shc(0.,sinh(0.5*chi));
  _mx[0][0]= chc;_mx[0][1]= shc;_mx[0][2]= 0. ;_mx[0][3]= 0. ;
  _mx[1][0]=-shc;_mx[1][1]= chc;_mx[1][2]= 0. ;_mx[1][3]= 0  ;
  _mx[2][0]= 0. ;_mx[2][1]= 0. ;_mx[2][2]= chc;_mx[2][3]=-shc;
  _mx[3][0]= 0. ;_mx[3][1]= 0. ;_mx[3][2]=+shc;_mx[3][3]= chc;
  return *this;
}

// lorentz boost in z direction 
SpinHalfLorentzRotation & SpinHalfLorentzRotation::setBoostZ (double & bz)
{
  // work out beta and chi
  double chi(atanh(bz)),shc(sinh(0.5*chi)),chc(cosh(0.5*chi));
  _mx[0][0]= chc-shc;_mx[0][1]= 0.     ;_mx[0][2]= 0.        ;_mx[0][3]= 0.     ;
  _mx[1][0]= 0.     ;_mx[1][1]= chc+shc;_mx[1][2]= 0.        ;_mx[1][3]= 0.     ;
  _mx[2][0]= 0.     ;_mx[2][1]= 0.     ;_mx[2][2]= chc+shc   ;_mx[2][3]= 0.     ;
  _mx[3][0]= 0.     ;_mx[3][1]= 0.     ;_mx[3][2]=+0.        ;_mx[3][3]= chc-shc;
  return *this;
}

// Pure boost along the x-axis; equivalent to LT = BoostX(beta) * LT
SpinHalfLorentzRotation & SpinHalfLorentzRotation::boostX(double bx)
{
  double chi(atanh(bx)),shc(sinh(0.5*chi)),chc(cosh(0.5*chi));
  MatrixT temp;
  for(size_t ix=0;ix<4;++ix)
    {
      temp[0][ix]= chc*_mx[0][ix]-shc*_mx[1][ix];
      temp[1][ix]= chc*_mx[1][ix]-shc*_mx[0][ix];
      temp[2][ix]= chc*_mx[2][ix]+shc*_mx[3][ix];
      temp[3][ix]= chc*_mx[3][ix]+shc*_mx[2][ix];
    }
  _mx = temp;
  return *this;
}
// Pure boost along the y-axis; equivalent to LT = BoostY(beta) * LT
SpinHalfLorentzRotation & SpinHalfLorentzRotation::boostY(double by)
{
  double chi(atanh(by)),chc(cosh(0.5*chi));
  Complex shc(0.,sinh(0.5*chi));
  MatrixT temp;
  for(size_t ix=0;ix<4;++ix)
    {
      temp[0][ix]= chc*_mx[0][ix]+shc*_mx[1][ix];
      temp[1][ix]= chc*_mx[1][ix]-shc*_mx[0][ix];
      temp[2][ix]= chc*_mx[2][ix]-shc*_mx[3][ix];
      temp[3][ix]= chc*_mx[3][ix]+shc*_mx[2][ix];
    }
  _mx = temp;
  return *this;
}

// Pure boost along the z-axis; equivalent to LT = BoostX(beta) * LT
SpinHalfLorentzRotation & SpinHalfLorentzRotation::boostZ(double bz)
{
  double chi(atanh(bz)),shc(sinh(0.5*chi)),chc(cosh(0.5*chi));
  MatrixT temp;
  for(size_t ix=0;ix<4;++ix)
    {
      temp[0][ix]=(chc-shc)*_mx[0][ix];
      temp[1][ix]=(chc+shc)*_mx[1][ix];
      temp[2][ix]=(chc+shc)*_mx[2][ix];
      temp[3][ix]=( chc-shc)*_mx[3][ix];
    }
  _mx = temp;
  return *this;
}

// General boost equivalent to LT = Boost(bx,by,bz) * LT
SpinHalfLorentzRotation & SpinHalfLorentzRotation::boost(double bx, double by, double bz, double gamma) {
  // calculation of gamma and beta
  double b2(bx*bx+by*by+bz*bz);
  if(gamma<1.) gamma = 1./sqrt(1.-b2);
  double beta = sqrt(b2);
  // work out beta and chi
  const double eps=1e-8;
  double chc = sqrt(0.5*(1+gamma));
  double shc = beta>eps ? sqrt(0.5*(gamma-1))/beta : 0.5+b2*(0.1875+0.12109375*b2);
  Complex ii(0.,1.),nxminy(bx-ii*by),nxplny(bx+ii*by);
  MatrixT temp;
  for(size_t ix=0;ix<4;++ix) {
    temp[0][ix]= (chc-shc*bz)*_mx[0][ix]-shc*nxminy  *_mx[1][ix];
    temp[1][ix]=-shc*nxplny  *_mx[0][ix]+(chc+shc*bz)*_mx[1][ix];
    temp[2][ix]= (chc+shc*bz)*_mx[2][ix]+shc*nxminy  *_mx[3][ix];
    temp[3][ix]= shc*nxplny  *_mx[2][ix]+(chc-shc*bz)*_mx[3][ix];
  }
  _mx = temp;
  return *this;
}

// General boost equivalent to LT = Boost(bv) * LT
SpinHalfLorentzRotation & SpinHalfLorentzRotation::boost(const Boost & b, double gamma) {
  // calculation of gamma and beta
  double b2(b.mag2());
  if(gamma<1.) gamma = 1./sqrt(1.-b2);
  double beta = sqrt(b2);
  // work out chi
  const double eps=1e-8;
  double bx(b.x()),by(b.y()),bz(b.z());
  double chc = sqrt(0.5*(1+gamma));
  double shc = beta>eps ? sqrt(0.5*(gamma-1))/beta : 0.5+b2*(0.1875+0.12109375*b2);
  Complex ii(0.,1.),nxminy(bx-ii*by),nxplny(bx+ii*by);
  MatrixT temp;
  for(size_t ix=0;ix<4;++ix) {
    temp[0][ix]= (chc-shc*bz)*_mx[0][ix]-shc*nxminy  *_mx[1][ix];
    temp[1][ix]=-shc*nxplny  *_mx[0][ix]+(chc+shc*bz)*_mx[1][ix];
    temp[2][ix]= (chc+shc*bz)*_mx[2][ix]+shc*nxminy  *_mx[3][ix];
    temp[3][ix]= shc*nxplny  *_mx[2][ix]+(chc-shc*bz)*_mx[3][ix];
  }
  _mx = temp;
  return *this;
}

std::ostream & SpinHalfLorentzRotation::print( std::ostream & os ) const {
  os << "\n   [ ( " <<
    std::setw(14) << std::setprecision(6) << s1s1() << "   " <<
    std::setw(14) << std::setprecision(6) << s1s2() << "   " <<
    std::setw(14) << std::setprecision(6) << s1s3() << "   " <<
    std::setw(14) << std::setprecision(6) << s1s4() << ")\n"
     << "     ( " <<
    std::setw(14) << std::setprecision(6) << s2s1() << "   " <<
    std::setw(14) << std::setprecision(6) << s2s2() << "   " <<
    std::setw(14) << std::setprecision(6) << s2s3() << "   " <<
    std::setw(14) << std::setprecision(6) << s2s4() << ")\n"
     << "     ( " <<
    std::setw(14) << std::setprecision(6) << s3s1() << "   " <<
    std::setw(14) << std::setprecision(6) << s3s2() << "   " <<
    std::setw(14) << std::setprecision(6) << s3s3() << "   " <<
    std::setw(14) << std::setprecision(6) << s3s4() << ")\n"
     << "     ( " <<
    std::setw(14) << std::setprecision(6) << s4s1() << "   " <<
    std::setw(14) << std::setprecision(6) << s4s2() << "   " <<
    std::setw(14) << std::setprecision(6) << s4s3() << "   " <<
    std::setw(14) << std::setprecision(6) << s4s4() << ") ]\n";
  return os;
}


// general rotation
SpinHalfLorentzRotation & 
SpinHalfLorentzRotation::setRotate(double phi, const Axis & axis) {
  double cp(cos(0.5*phi));
  // get the normalised components of the vector
  double amag(axis.mag()),ax(axis.x()/amag),ay(axis.y()/amag),az(axis.z()/amag);
  Complex ii(0.,1.),nxminy(ax-ii*ay),nxplny(ax+ii*ay),isp(0.,sin(0.5*phi));
  // rotatation matrix is the same in both conventions
  _mx[0][0]= cp-isp*az ;_mx[0][1]=-isp*nxminy;_mx[0][2]= 0.        ;_mx[0][3]= 0.        ;
  _mx[1][0]=-isp*nxplny;_mx[1][1]= cp+isp*az ;_mx[1][2]= 0.        ;_mx[1][3]= 0.        ;
  _mx[2][0]= 0.        ;_mx[2][1]= 0.        ;_mx[2][2]= cp-isp*az ;_mx[2][3]=-isp*nxminy;
  _mx[3][0]= 0.        ;_mx[3][1]= 0.        ;_mx[3][2]=-isp*nxplny;_mx[3][3]= cp+isp*az ;
  return *this;
}

// rotation about x
SpinHalfLorentzRotation & SpinHalfLorentzRotation::setRotateX(double& phi) {
  double cp(cos(0.5*phi));
  Complex isp(0.,sin(0.5*phi));
  // rotatation matrix is the same in both conventions
  _mx[0][0]= cp ;_mx[0][1]=-isp;_mx[0][2]= 0. ;_mx[0][3]= 0. ;
  _mx[1][0]=-isp;_mx[1][1]= cp ;_mx[1][2]= 0. ;_mx[1][3]= 0. ;
  _mx[2][0]= 0. ;_mx[2][1]= 0. ;_mx[2][2]= cp ;_mx[2][3]=-isp;
  _mx[3][0]= 0. ;_mx[3][1]= 0. ;_mx[3][2]=-isp;_mx[3][3]= cp ;
  return *this;
}

// rotation about y
SpinHalfLorentzRotation & SpinHalfLorentzRotation::setRotateY(double& phi) {
  double cp(cos(0.5*phi)),sp(sin(0.5*phi));
  // rotatation matrix is the same in both conventions
  _mx[0][0]= cp;_mx[0][1]=-sp;_mx[0][2]= 0.;_mx[0][3]= 0.;
  _mx[1][0]= sp;_mx[1][1]= cp;_mx[1][2]= 0.;_mx[1][3]= 0.;
  _mx[2][0]= 0.;_mx[2][1]= 0.;_mx[2][2]= cp;_mx[2][3]=-sp;
  _mx[3][0]= 0.;_mx[3][1]= 0.;_mx[3][2]= sp;_mx[3][3]= cp;
  return *this;
}

// rotation about z
SpinHalfLorentzRotation & SpinHalfLorentzRotation::setRotateZ(double& phi)
{
  double cp(cos(0.5*phi));
  Complex isp(0.,sin(0.5*phi));
  // rotatation matrix is the same in both conventions
  _mx[0][0]= cp-isp ;_mx[0][1]= 0.    ;_mx[0][2]= 0.    ;_mx[0][3]= 0.    ;
  _mx[1][0]= 0.     ;_mx[1][1]= cp+isp;_mx[1][2]= 0.    ;_mx[1][3]= 0.    ;
  _mx[2][0]= 0.     ;_mx[2][1]= 0.    ;_mx[2][2]= cp-isp;_mx[2][3]= 0.    ;
  _mx[3][0]= 0.     ;_mx[3][1]= 0.    ;_mx[3][2]= 0.    ;_mx[3][3]= cp+isp;
  return *this;
}


// product
SpinHalfLorentzRotation 
SpinHalfLorentzRotation::operator * (const SpinHalfLorentzRotation & lt) const {
  SpinHalfLorentzRotation temp{_mx};
  temp *= lt;
  return temp;
}

// multiply and assign
SpinHalfLorentzRotation & 
SpinHalfLorentzRotation::operator *= (const SpinHalfLorentzRotation & lt) {
  MatrixT temp;
  for(size_t ix=0;ix<4;++ix) {
    for(size_t iy=0;iy<4;++iy) {
      temp[ix][iy] = 0.0;
      for(size_t iz=0;iz<4;++iz)
        temp[ix][iy] += _mx[ix][iz] * lt._mx[iz][iy];
    }
  }
  _mx = temp;
  return *this;
}

// transform method
SpinHalfLorentzRotation & 
SpinHalfLorentzRotation::transform(const SpinHalfLorentzRotation & lt) {
  SpinHalfLorentzRotation temp{lt._mx};
  temp *= (*this);
  _mx = temp._mx;
  return *this;
}

// Rotation around the x-axis; equivalent to LT = RotationX(delta) * LT
SpinHalfLorentzRotation & SpinHalfLorentzRotation::rotateX(double phi) {
  double cp(cos(0.5*phi));
  Complex isp(0.,sin(0.5*phi));
  MatrixT temp;
  for(size_t ix=0;ix<4;++ix) {
    temp[0][ix]=  cp*_mx[0][ix]-isp*_mx[1][ix];
    temp[1][ix]=-isp*_mx[0][ix]+ cp*_mx[1][ix];
    temp[2][ix]=  cp*_mx[2][ix]-isp*_mx[3][ix];
    temp[3][ix]=-isp*_mx[2][ix]+ cp*_mx[3][ix];
  }
  _mx = temp;
  return *this;
}

// Rotation around the y-axis; equivalent to LT = RotationY(delta) * LT
SpinHalfLorentzRotation & SpinHalfLorentzRotation::rotateY(double phi) {
  double cp(cos(0.5*phi)),sp(sin(0.5*phi));
  MatrixT temp;
  for(size_t ix=0;ix<4;++ix) {
    temp[0][ix]= cp*_mx[0][ix]-sp*_mx[1][ix];
    temp[1][ix]= sp*_mx[0][ix]+cp*_mx[1][ix];
    temp[2][ix]= cp*_mx[2][ix]-sp*_mx[3][ix];
    temp[3][ix]= sp*_mx[2][ix]+cp*_mx[3][ix];
  }
  _mx = temp;
  return *this;
}

// Rotation around the z-axis; equivalent to LT = RotationZ(delta) * LT
SpinHalfLorentzRotation & SpinHalfLorentzRotation::rotateZ(double phi) {
  double cp(cos(0.5*phi));
  Complex isp(0.,sin(0.5*phi));
  MatrixT temp;
  for(size_t ix=0;ix<4;++ix) {
    temp[0][ix]= (cp-isp)*_mx[0][ix];
    temp[1][ix]= (cp+isp)*_mx[1][ix];
    temp[2][ix]= (cp-isp)*_mx[2][ix];
    temp[3][ix]= (cp+isp)*_mx[3][ix];
  }
  _mx = temp;
  return *this;
}

// Rotation around specified vector - LT = Rotation(delta,axis)*LT
SpinHalfLorentzRotation & 
SpinHalfLorentzRotation::rotate(double phi, const Axis & axis) {
  double cp(cos(0.5*phi)),amag(axis.mag()),
    ax(axis.x()/amag),ay(axis.y()/amag),az(axis.z()/amag);
  Complex ii(0.,1.),nxminy(ax-ii*ay),nxplny(ax+ii*ay),isp(0.,sin(0.5*phi));
  MatrixT temp;
  for(size_t ix=0;ix<4;++ix) {
    temp[0][ix]= (cp-isp*az)*_mx[0][ix]-isp*nxminy *_mx[1][ix];
    temp[1][ix]=-isp*nxplny *_mx[0][ix]+(cp+isp*az)*_mx[1][ix];
    temp[2][ix]= (cp-isp*az)*_mx[2][ix]-isp*nxminy *_mx[3][ix];
    temp[3][ix]=-isp*nxplny *_mx[2][ix]+(cp+isp*az)*_mx[3][ix];
  }
  _mx = temp;
  return *this;
}
