// -*- C++ -*-
//
// SpinOneLorentzRotation.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LorentzRotation class.
//

#include "SpinOneLorentzRotation.h"

using namespace ThePEG;
SpinOneLorentzRotation::
SpinOneLorentzRotation(double xx, double xy, double xz, double xt,
		       double yx, double yy, double yz, double yt,
		       double zx, double zy, double zz, double zt,
		       double tx, double ty, double tz, double tt)
{
  xx_() = xx; xy_() = xy; xz_() = xz; xt_() = xt;
  yx_() = yx; yy_() = yy; yz_() = yz; yt_() = yt;
  zx_() = zx; zy_() = zy; zz_() = zz; zt_() = zt;
  tx_() = tx; ty_() = ty; tz_() = tz; tt_() = tt;
}

bool SpinOneLorentzRotation::isIdentity() const {
  return 
       1.0 == xx() && 0.0 == xy() && 0.0 == xz() && 0.0 == xt()
    && 0.0 == yx() && 1.0 == yy() && 0.0 == yz() && 0.0 == yt()
    && 0.0 == zx() && 0.0 == zy() && 1.0 == zz() && 0.0 == zt()
    && 0.0 == tx() && 0.0 == ty() && 0.0 == tz() && 1.0 == tt();
}

SpinOneLorentzRotation SpinOneLorentzRotation::inverse() const {
  return SpinOneLorentzRotation( xx(), yx(), zx(),-tx(),
				 xy(), yy(), zy(),-ty(), 
				 xz(), yz(), zz(),-tz(),
				-xt(),-yt(),-zt(), tt());
}

// output operator
std::ostream & SpinOneLorentzRotation::print( std::ostream &  os) const {  
  os << "\n   [ ( " <<
    std::setw(14) << std::setprecision(6) << xx() << "   " <<
    std::setw(14) << std::setprecision(6) << xy() << "   " <<
    std::setw(14) << std::setprecision(6) << xz() << "   " <<
    std::setw(14) << std::setprecision(6) << xt() << ")\n"
     << "     ( " <<
    std::setw(14) << std::setprecision(6) << yx() << "   " <<
    std::setw(14) << std::setprecision(6) << yy() << "   " <<
    std::setw(14) << std::setprecision(6) << yz() << "   " <<
    std::setw(14) << std::setprecision(6) << yt() << ")\n"
     << "     ( " <<
    std::setw(14) << std::setprecision(6) << zx() << "   " <<
    std::setw(14) << std::setprecision(6) << zy() << "   " <<
    std::setw(14) << std::setprecision(6) << zz() << "   " <<
    std::setw(14) << std::setprecision(6) << zt() << ")\n"
     << "     ( " <<
    std::setw(14) << std::setprecision(6) << tx() << "   " <<
    std::setw(14) << std::setprecision(6) << ty() << "   " <<
    std::setw(14) << std::setprecision(6) << tz() << "   " <<
    std::setw(14) << std::setprecision(6) << tt() << ") ]\n";
  return os;
}

SpinOneLorentzRotation & SpinOneLorentzRotation::
setBoost (double bx, double by, double bz, double gamma) {
  double beta2 = sqr(bx) + sqr(by) + sqr(bz);
  if (beta2 >= 1.0) {
    throw Exception() << "Invalid boost (" << bx << ',' << by << ',' << bz 
		      << ") in SpinOneLorentzRotatio::setBoost" << Exception::eventerror;
  }
  if(gamma<0.) gamma = 1.0 / sqrt((1.-bz)*(1.+bz)-sqr(bx)-sqr(by));
  double bgamma = sqr(gamma) / (1.0 + gamma);
  xx_() = 1.0 + bgamma * bx * bx;
  yy_() = 1.0 + bgamma * by * by;
  zz_() = 1.0 + bgamma * bz * bz;
  xy_() = yx_() = bgamma * bx * by;
  xz_() = zx_() = bgamma * bx * bz;
  yz_() = zy_() = bgamma * by * bz;
  xt_() = tx_() = gamma * bx;
  yt_() = ty_() = gamma * by;
  zt_() = tz_() = gamma * bz;
  tt_() = gamma;
  return *this;
}

SpinOneLorentzRotation & 
SpinOneLorentzRotation::setRotate(double delta, const Axis & axis) {

  double sinDelta = sin(delta), cosDelta = cos(delta);
  double oneMinusCosDelta = 1.0 - cosDelta;

  Axis u = unitVector(axis);

  double uX = u.x();
  double uY = u.y();
  double uZ = u.z();

  double rxx = oneMinusCosDelta * uX * uX  +  cosDelta;
  double rxy = oneMinusCosDelta * uX * uY  -  sinDelta * uZ;
  double rxz = oneMinusCosDelta * uX * uZ  +  sinDelta * uY;

  double ryx = oneMinusCosDelta * uY * uX  +  sinDelta * uZ;
  double ryy = oneMinusCosDelta * uY * uY  +  cosDelta;
  double ryz = oneMinusCosDelta * uY * uZ  -  sinDelta * uX;

  double rzx = oneMinusCosDelta * uZ * uX  -  sinDelta * uY;
  double rzy = oneMinusCosDelta * uZ * uY  +  sinDelta * uX;
  double rzz = oneMinusCosDelta * uZ * uZ  +  cosDelta;

  xx_() = rxx; xy_() = rxy; xz_() = rxz; xt_() = 0.0;
  yx_() = ryx; yy_() = ryy; yz_() = ryz; yt_() = 0.0;
  zx_() = rzx; zy_() = rzy; zz_() = rzz; zt_() = 0.0;
  tx_() = 0.0; ty_() = 0.0; tz_() = 0.0; tt_() = 1.0;

  return  *this;
}

SpinOneLorentzRotation & 
SpinOneLorentzRotation::setRotateX(double delta) {
  double sinDelta = sin(delta), cosDelta = cos(delta);

  double ryy = cosDelta, ryz = -sinDelta;
  double rzy = sinDelta, rzz =  cosDelta;

  xx_() = 1.0; xy_() = 0.0; xz_() = 0.0; xt_() = 0.0;
  yx_() = 0.0; yy_() = ryy; yz_() = ryz; yt_() = 0.0;
  zx_() = 0.0; zy_() = rzy; zz_() = rzz; zt_() = 0.0;
  tx_() = 0.0; ty_() = 0.0; tz_() = 0.0; tt_() = 1.0;
  return  *this;
}

SpinOneLorentzRotation & 
SpinOneLorentzRotation::setRotateY(double delta) {
  double sinDelta = sin(delta), cosDelta = cos(delta);

  double rxx =  cosDelta, rxz =  sinDelta;
  double rzx = -sinDelta, rzz =  cosDelta;

  xx_() = rxx; xy_() = 0.0; xz_() = rxz; xt_() = 0.0;
  yx_() = 0.0; yy_() = 1.0; yz_() = 0.0; yt_() = 0.0;
  zx_() = rzx; zy_() = 0.0; zz_() = rzz; zt_() = 0.0;
  tx_() = 0.0; ty_() = 0.0; tz_() = 0.0; tt_() = 1.0;
  return  *this;
}

SpinOneLorentzRotation & 
SpinOneLorentzRotation::setRotateZ(double delta) {
  double sinDelta = sin(delta), cosDelta = cos(delta);

  double rxx = cosDelta, rxy = -sinDelta;
  double ryx = sinDelta, ryy =  cosDelta;

  xx_() = rxx; xy_() = rxy; xz_() = 0.0; xt_() = 0.0;
  yx_() = ryx; yy_() = ryy; yz_() = 0.0; yt_() = 0.0;
  zx_() = 0.0; zy_() = 0.0; zz_() = 1.0; zt_() = 0.0;
  tx_() = 0.0; ty_() = 0.0; tz_() = 0.0; tt_() = 1.0;
  return  *this;
}

SpinOneLorentzRotation 
SpinOneLorentzRotation::operator*(const SpinOneLorentzRotation & b) const {
  return SpinOneLorentzRotation
    (xx()*b.xx() + xy()*b.yx() + xz()*b.zx() + xt()*b.tx(),
     xx()*b.xy() + xy()*b.yy() + xz()*b.zy() + xt()*b.ty(),
     xx()*b.xz() + xy()*b.yz() + xz()*b.zz() + xt()*b.tz(),
     xx()*b.xt() + xy()*b.yt() + xz()*b.zt() + xt()*b.tt(),
     
     yx()*b.xx() + yy()*b.yx() + yz()*b.zx() + yt()*b.tx(),
     yx()*b.xy() + yy()*b.yy() + yz()*b.zy() + yt()*b.ty(),
     yx()*b.xz() + yy()*b.yz() + yz()*b.zz() + yt()*b.tz(),
     yx()*b.xt() + yy()*b.yt() + yz()*b.zt() + yt()*b.tt(),
     
     zx()*b.xx() + zy()*b.yx() + zz()*b.zx() + zt()*b.tx(),
     zx()*b.xy() + zy()*b.yy() + zz()*b.zy() + zt()*b.ty(),
     zx()*b.xz() + zy()*b.yz() + zz()*b.zz() + zt()*b.tz(),
     zx()*b.xt() + zy()*b.yt() + zz()*b.zt() + zt()*b.tt(),

     tx()*b.xx() + ty()*b.yx() + tz()*b.zx() + tt()*b.tx(),
     tx()*b.xy() + ty()*b.yy() + tz()*b.zy() + tt()*b.ty(),
     tx()*b.xz() + ty()*b.yz() + tz()*b.zz() + tt()*b.tz(),		
     tx()*b.xt() + ty()*b.yt() + tz()*b.zt() + tt()*b.tt());
}

