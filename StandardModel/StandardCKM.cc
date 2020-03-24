// -*- C++ -*-
//
// StandardCKM.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the StandardCKM class.
//

#include "StandardCKM.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;

IBPtr StandardCKM::clone() const {
  return new_ptr(*this);
}

IBPtr StandardCKM::fullclone() const {
  return new_ptr(*this);
}

vector< vector<double> > StandardCKM::getMatrix(unsigned int nFamilies) const {
  vector< vector<double> > ckm(nFamilies, vector<double>(nFamilies, 0.0));
  for ( unsigned int i = 0; i < nFamilies; ++i ) ckm[i][i] = 1.0;
  if ( nFamilies <= 1 ) return ckm;
  double s12 = sin(theta12);
  double c12 = cos(theta12);
  if ( nFamilies == 2 ) {
    ckm[0][0] = sqr(c12);
    ckm[0][1] = sqr(s12);
    ckm[1][0] = sqr(s12);
    ckm[1][1] = sqr(c12);
    return ckm;
  }
  double s13 = sin(theta13);
  double c13 = cos(theta13);
  double s23 = sin(theta23);
  double c23 = cos(theta23);
  double cd = cos(delta);
  ckm[0][0] = sqr(c12*c13);
  ckm[0][1] = sqr(s12*c13);
  ckm[0][2] = sqr(s13);
  ckm[1][0] = sqr(s12*c23)+sqr(c12*s23*s13)+2.0*s12*c23*c12*s23*s13*cd;
  ckm[1][1] = sqr(c12*c23)+sqr(s12*s23*s13)-2.0*c12*c23*s12*s23*s13*cd;
  ckm[1][2] = sqr(s23*c13);
  ckm[2][0] = sqr(s12*s23)+sqr(c12*c23*s13)-2.0*s12*s23*c12*c23*s13*cd;
  ckm[2][1] = sqr(c12*s23)+sqr(s12*c23*s13)+2.0*c12*s23*s12*c23*s13*cd;
  ckm[2][2] = sqr(c23*c13);
  return ckm;
}

void StandardCKM::persistentOutput(PersistentOStream & os) const {
  os << theta12 << theta13 << theta23 << delta;
}

void StandardCKM::persistentInput(PersistentIStream & is, int) {
  is >> theta12 >> theta13 >> theta23 >> delta;
}

ClassDescription<StandardCKM> StandardCKM::initStandardCKM;

void StandardCKM::Init() {

  static ClassDocumentation<StandardCKM> documentation
    ("Implements the standard parameterization of the CKM matrix in terms "
     "of three angles and a phase.");

  static Parameter<StandardCKM,double> interfaceTheta12
    ("theta_12",
     "The mixing angle between the first and second generation in the standard "
     "parameterization of the CKM matrix",
     &StandardCKM::theta12, 0.222357, 0.0, Constants::twopi, false, false, true);

  static Parameter<StandardCKM,double> interfaceTheta13
    ("theta_13",
     "The mixing angle between the first and third generation in the standard "
     "parameterization of the CKM matrix",
     &StandardCKM::theta13, 0.0003150, 0.0, Constants::twopi, false, false, true);

  static Parameter<StandardCKM,double> interfaceTheta23
    ("theta_23",
     "The mixing angle between the second and third generation in the standard "
     "parameterization of the CKM matrix",
     &StandardCKM::theta23, 0.039009, 0.0, Constants::twopi, false, false, true);

  static Parameter<StandardCKM,double> interfaceDelta
    ("delta",
     "The phase angle in the standard "
     "parameterization of the CKM matrix",
     &StandardCKM::delta, 1.35819, 0.0, Constants::twopi, false, false, true);

  interfaceTheta12.rank(10);
  interfaceTheta13.rank(9);
  interfaceTheta23.rank(8);
  interfaceDelta.rank(7);

}

