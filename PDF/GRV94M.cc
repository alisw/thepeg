// -*- C++ -*-
//
// GRV94M.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GRV94M class.
//

#include "GRV94M.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;

IBPtr GRV94M::clone() const {
  return new_ptr(*this);
}

IBPtr GRV94M::fullclone() const {
  return new_ptr(*this);
}

void GRV94M::setup(double l, Energy2 scale) const {
  GRVBase::setup(l, scale, mu2(), lam2());
}

double GRV94M::uv() const {
  return valens(1.304 + 0.863*S(),
		0.558 - 0.020*S(),
		0.183*S(),
		-0.113 + 0.283*S() - 0.321*S2(),
		6.843 - 5.089*S() + 2.647*S2() - 0.527*S3(),
		7.771 - 10.09*S() + 2.630*S2(),
		3.315 + 1.145*S() - 0.583*S2() + 0.154*S3());
}

double GRV94M::dv() const {
  return valens(0.102 - 0.017*S() + 0.005*S2(),
		0.270 - 0.019*S(),
		0.260,
		2.393 + 6.228*S() - 0.881*S2(),
		46.06 + 4.673*S() - 14.98*S2() + 1.331*S3(),
		17.83 - 53.47*S() + 21.24*S2(),
		4.081 + 0.976*S() - 0.485*S2() + 0.152*S3());
}

double GRV94M::del() const {
  return 0.5*valens(0.070 + 0.042*S() - 0.011*S2() + 0.004*S3(),
		    0.409 - 0.007*S(),
		    0.782 + 0.082*S(),
		    -29.65 + 26.49*S() + 5.429*S2(),
		    90.20 - 74.97*S() + 4.526*S2(),
		    0.0,
		    8.122 + 2.120*S() - 1.088*S2() + 0.231*S3());
}

double GRV94M::udb() const {
  return 0.5*lightsea(0.877,
		      0.561,
		      0.275,
		      0.0,
		      0.997,
		      3.210 - 1.866*S(),
		      7.300,
		      9.010 + 0.896*rootS() + 0.222*S2(),
		      3.077 + 1.446*S(),
		      3.173 - 2.445*rootS() + 2.207*S());
}

double GRV94M::sb() const {
  return heavysea(0.0,
		  0.756,
		  0.216,
		  1.690 + 0.650*rootS() - 0.922*S(),
		  -4.329 + 1.131*S(),
		  9.568 - 1.744*S(),
		  9.377 + 1.088*rootS() - 1.320*S() + 0.130*S2(),
		  3.031 + 1.639*S(),
		  5.837 + 0.815*S());
}

double GRV94M::cb() const {
  return heavysea(0.820,
		  0.98,
		  0.0,
		  -0.625 - 0.523*S(),
		  0.0,
		  1.896 + 1.616*S(),
		  4.12  + 0.683*S(),
		  4.36  + 1.328*S(),
		  0.677 + 0.679*S());
}

double GRV94M::bb() const {
  return heavysea(1.297,
		  0.99,
		  0.0,
		  - 0.193*S(),
		  0.0,
		  0.0,
		  3.447 + 0.927*S(),
		  4.68  + 1.259*S(),
		  1.892 + 2.199*S());
}

double GRV94M::gl() const {
  return lightsea(1.014,
		  1.738,
		  1.724 + 0.157*S(),
		  0.800 + 1.016*S(),
		  7.517 - 2.547*S(),
		  34.09 - 52.21*rootS() + 17.47*S(),
		  4.039 + 1.491*S(),
		  3.404 + 0.830*S(),
		  -1.112 + 3.438*S()  - 0.302*S2(),
		  3.256 - 0.436*S());
}

NoPIOClassDescription<GRV94M> GRV94M::initGRV94M;

void GRV94M::Init() {

  static ClassDocumentation<GRV94M> documentation
    ("Implements the GRV94M PDF parameterization.");

}

