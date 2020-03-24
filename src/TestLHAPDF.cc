// -*- C++ -*-
//
// TestLHAPDF.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TestLHAPDF class.
//

#include "TestLHAPDF.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDF/PDFBase.h"

using namespace ThePEG;

NoPIOClassDescription<TestLHAPDF> TestLHAPDF::initTestLHAPDF;
// Definition of the static class description member.

bool approx(double x1, double x2, double eps = 1.0e-3) {
  eps *= abs(x1) + abs(x2);
  if ( eps == 0.0 ) return true;
  return abs(x1 - x2) < eps;
}

void TestLHAPDF::Init() {

  double PDF1P[]  = { 0.000000000,  0.254411259,  0.551724644,  0.679716202,
		      0.871195736,  0.872211343, 23.791350673,  0.919026004,
		      0.947097315,  0.679716202,  0.551724644,  0.254411259,
		      0.000000000 };
  double PDF1PB[] = { 0.000000000,  0.254411259,  0.551724644,  0.679716202,
		      0.947097315,  0.919026004, 23.791350673,  0.872211343,
		      0.871195736,  0.679716202,  0.551724644,  0.254411259,
		      0.000000000 };
  double PDF1N[]  = { 0.000000000,  0.254411259,  0.551724644,  0.679716202,
		      0.872211343,  0.871195736, 23.791350673,  0.947097315,
		      0.919026004,  0.679716202,  0.551724644,  0.254411259,
		      0.000000000 };
  double PDF1NB[] = { 0.000000000,  0.254411259,  0.551724644,  0.679716202,
		      0.919026004,  0.947097315, 23.791350673,  0.871195736,
		      0.872211343,  0.679716202,  0.551724644,  0.254411259,
		      0.000000000 };
  double PDFVP[]  = { 0.000000000,  0.000000000,  0.000000000,  0.000000000,
		      0.000000000,  0.000000000,  0.000000000,
		      0.919026004 - 0.872211343,  0.947097315 - 0.871195736,
		      0.000000000,  0.000000000,  0.000000000,
		      0.000000000 };
  double PDFVPB[] = { 0.000000000,  0.000000000,  0.000000000,  0.000000000,
		      0.947097315 - 0.871195736,  0.919026004 - 0.872211343,
		      0.000000000,  0.000000000,  0.000000000,  0.000000000,
		      0.000000000,  0.000000000,
		      0.000000000 };
  double PDFVN[]  = { 0.000000000,  0.000000000,  0.000000000,  0.000000000,
		      0.000000000,  0.000000000,  0.000000000,
		      0.947097315 - 0.871195736,  0.919026004 - 0.87221134,
		      0.000000000,  0.000000000,  0.000000000,
		      0.000000000 };
  double PDFVNB[] = { 0.000000000,  0.000000000,  0.000000000,  0.000000000,
		      0.919026004 - 0.872211343,  0.947097315 - 0.871195736,
		      0.000000000,  0.000000000,  0.000000000,  0.000000000,
		      0.000000000,  0.000000000,  0.00000000 };
  double PDF2N[]  = { 0.000000000,  0.277690000,  0.563260000,  0.756930000,
 	              0.894200000,  0.893470000, 26.288140000,  0.963870000, 
 	              0.945850000,  0.756930000,  0.563260000,  0.277690000, 
 	              0.000000000 };
  double PDF3NB[] = { 0.000000000,  0.257616025,  0.551677676,  0.551677676,
                      1.142526340,  1.092702640, 24.771891800,  1.031855390,
                      1.105400690,  0.551677676,  0.551677676,  0.257616025,
                      0.000000000 };
  double PDF4PB[] = { 0.000000000,  0.000000000,  0.551300347,  0.679214040,
		      0.946471989,  0.918402827, 23.770542200,  0.871592366,
		      0.870576794,  0.679214040,  0.551300347,  0.000000000,
		      0.000000000 };

  PDPtr proton = eventGenerator()->getParticleData(ParticleID::pplus);
  PDPtr neutron = eventGenerator()->getParticleData(ParticleID::n0);
  PDPtr pbar = eventGenerator()->getParticleData(ParticleID::pbarminus);
  PDPtr nbar = eventGenerator()->getParticleData(ParticleID::nbar0);
  PDPtr pi0 = eventGenerator()->getParticleData(ParticleID::pi0);
  PDVector partons;
  partons.push_back(eventGenerator()->getParticleData(ParticleID::tbar));
  partons.push_back(eventGenerator()->getParticleData(ParticleID::bbar));
  partons.push_back(eventGenerator()->getParticleData(ParticleID::cbar));
  partons.push_back(eventGenerator()->getParticleData(ParticleID::sbar));
  partons.push_back(eventGenerator()->getParticleData(ParticleID::ubar));
  partons.push_back(eventGenerator()->getParticleData(ParticleID::dbar));
  partons.push_back(eventGenerator()->getParticleData(ParticleID::g));
  partons.push_back(eventGenerator()->getParticleData(ParticleID::d));
  partons.push_back(eventGenerator()->getParticleData(ParticleID::u));
  partons.push_back(eventGenerator()->getParticleData(ParticleID::s));
  partons.push_back(eventGenerator()->getParticleData(ParticleID::c));
  partons.push_back(eventGenerator()->getParticleData(ParticleID::b));
  partons.push_back(eventGenerator()->getParticleData(ParticleID::t));

  PDFPtr pdf1 =
    eventGenerator()->getObject<PDFBase>("/LHAPDF/CTEQ6ll");
  assert(pdf1);
  for ( int i = 0; i < 13; ++i )
    if ( !approx(pdf1->xfx(proton, partons[i], 100.0*GeV2, 0.001), PDF1P[i]) )
      exit(1);
  for ( int i = 0; i < 13; ++i )
    if ( !approx(pdf1->xfx(pbar, partons[i], 100.0*GeV2, 0.001), PDF1PB[i]) )
      exit(1);
  for ( int i = 0; i < 13; ++i )
    if ( !approx(pdf1->xfx(neutron, partons[i], 100.0*GeV2, 0.001), PDF1N[i]) )
      exit(1);
  for ( int i = 0; i < 13; ++i )
    if ( !approx(pdf1->xfx(nbar, partons[i], 100.0*GeV2, 0.001), PDF1NB[i]) )
      exit(1);
  
  for ( int i = 0; i < 13; ++i )
    if ( !approx(pdf1->xfvx(proton, partons[i], 100.0*GeV2, 0.001), PDFVP[i]) )
      exit(1);
  for ( int i = 0; i < 13; ++i )
    if ( !approx(pdf1->xfvx(pbar, partons[i], 100.0*GeV2, 0.001), PDFVPB[i]) )
      exit(1);
  for ( int i = 0; i < 13; ++i )
    if ( !approx(pdf1->xfvx(neutron, partons[i], 100.0*GeV2, 0.001), PDFVN[i]) )
      exit(1);
  for ( int i = 0; i < 13; ++i )
    if ( !approx(pdf1->xfvx(nbar, partons[i], 100.0*GeV2, 0.001), PDFVNB[i]) )
      exit(1);

  PDFPtr pdf2 =
    eventGenerator()->getObject<PDFBase>("/LHAPDF/MRST2007");
  assert(pdf2);
  for ( int i = 0; i < 13; ++i )
    if ( !approx(pdf2->xfx(neutron, partons[i], 100.0*GeV2, 0.001), PDF2N[i]) )
      exit(1);

  PDFPtr pdf3 =
    eventGenerator()->getObject<PDFBase>("/LHAPDF/CT09MC1");
  assert(pdf3);
  for ( int i = 0; i < 13; ++i )
    if ( !approx(pdf3->xfx(nbar, partons[i], 100.0*GeV2, 0.001), PDF3NB[i]) )
      exit(1);

  PDFPtr pdf4 =
    eventGenerator()->getObject<PDFBase>("/LHAPDF/CTEQ6ll4");
  assert(pdf4);
  for ( int i = 0; i < 13; ++i )
    if ( !approx(pdf4->xfx(pbar, partons[i], 100.0*GeV2, 0.001), PDF4PB[i]) )
      exit(1);  
}
