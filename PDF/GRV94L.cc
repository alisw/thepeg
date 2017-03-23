// -*- C++ -*-
//
// GRV94L.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GRV94L class.
//

#include "GRV94L.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;

IBPtr GRV94L::clone() const {
  return new_ptr(*this);
}

IBPtr GRV94L::fullclone() const {
  return new_ptr(*this);
}

void GRV94L::setup(double l, Energy2 scale) const {
  GRVBase::setup(l, scale, mu2(), lam2());
}

double GRV94L::uv() const {
  return valens(2.284 + 0.802*S() + 0.055*S2(),
		0.590 - 0.024*S(),
		0.131 + 0.063*S(),
		-0.449 - 0.138*S() - 0.076*S2(),
		0.213 + 2.669*S() - 0.728*S2(),
		8.854 - 9.135*S() + 1.979*S2(),
		2.997 + 0.753*S() - 0.076*S2());
}

double GRV94L::dv() const {
  return valens(0.371 + 0.083*S() + 0.039*S2(),
		0.376,
		0.486 + 0.062*S(),
		-0.509 + 3.310*S() - 1.248*S2(),
		12.41 - 10.52*S() + 2.267*S2(),
		6.373 - 6.208*S() + 1.418*S2(),
		3.691 + 0.799*S() - 0.071*S2());
}

double GRV94L::del() const {
  return 0.5*valens(0.082 + 0.014*S() + 0.008*S2(),
		    0.409 - 0.005*S(),
		    0.799 + 0.071*S(),
		    -38.07 + 36.13*S() - 0.656*S2(),
		    90.31 - 74.15*S() + 7.645*S2(),
		    0.0,
		    7.486 + 1.217*S() - 0.159*S2());
}

double GRV94L::udb() const {
  return 0.5*lightsea(1.451,
		      0.271,
		      0.410 - 0.232*S(),
		      0.534 - 0.457*S(),
		      0.890 - 0.140*S(),
		      -0.981,
		      0.320 + 0.683*S(),
		      4.752 + 1.164*S() + 0.286*S2(),
		      4.119 + 1.713*S(),
		      0.682 + 2.978*S());
}

double GRV94L::sb() const {
  return heavysea(0.0,
		  0.914,
		  0.577,
		  1.798 - 0.596*S(),
		  -5.548 + 3.669*rootS() - 0.616*S(),
		  18.92 - 16.73*rootS() + 5.168*S(),
		  6.379 - 0.350*S()  + 0.142*S2(),
		  3.981 + 1.638*S(),
		  6.402);
}

double GRV94L::cb() const {
  return heavysea(0.888,
		  1.01,
		  0.37,
		  0.0,
		  0.0,
		  4.24  - 0.804*S(),
		  3.46  - 1.076*S(),
		  4.61  + 1.49 *S(),
		  2.555 + 1.961*S());
}

double GRV94L::bb() const {
  return heavysea(1.351,
		  1.00,
		  0.51,
		  0.0,
		  0.0,
		  1.848,
		  2.929 + 1.396*S(),
		  4.71  + 1.514*S(),
		  4.02  + 1.239*S());
}

double GRV94L::gl() const {
  return lightsea(0.524,
		  1.088,
		  1.742 - 0.930*S(),
		  - 0.399*S2(),
		  7.486 - 2.185*S(),
		  16.69 - 22.74*S()  + 5.779*S2(),
		  -25.59 + 29.71*S()  - 7.296*S2(),
		  2.792 + 2.215*S()  + 0.422*S2() - 0.104*S3(),
		  0.807 + 2.005*S(),
		  3.841 + 0.316*S());
}

NoPIOClassDescription<GRV94L> GRV94L::initGRV94L;

void GRV94L::Init() {

  static ClassDocumentation<GRV94L> documentation
    ("Implements the GRV94L PDF parameterization.");

}

