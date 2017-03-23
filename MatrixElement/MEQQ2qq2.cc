// -*- C++ -*-
//
// MEQQ2qq2.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEQQ2qq class.
//

#include "MEQQ2qq2.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace ThePEG;

void MEQQ2qq::getDiagrams() const {
  tcPDPtr g = getParticleData(ParticleID::g);
  for ( int i = 1; i <= maxFlavour(); ++i ) {
    tcPDPtr q = getParticleData(i);
    tcPDPtr qb = q->CC();
    for ( int j = 1; j <= maxFlavour(); ++j ) {
      if ( i == j ) continue;
      tcPDPtr qp = getParticleData(j);
      tcPDPtr qbp = qp->CC();
      add(new_ptr((Tree2toNDiagram(2), q, qb, 1, g, 3, qp, 3, qbp, -1)));
    }
  }
}

double MEQQ2qq::me2() const {
  return comfac()*colA()*KfacA()*2.0/9.0;
}

Selector<const ColourLines *>
MEQQ2qq::colourGeometries(tcDiagPtr) const {
  Selector<const ColourLines *> sel;

  static ColourLines csST("1 3 4, -5 -3 -2");

  sel.insert(1.0, &csST);
  return sel;
}

Selector<MEQQ2qq::DiagramIndex>
MEQQ2qq::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) 
    if ( diags[i]->id() == -1 ) sel.insert(1.0, i);
  return sel;
}

NoPIOClassDescription<MEQQ2qq> MEQQ2qq::initMEQQ2qq;
// Definition of the static class description member.

void MEQQ2qq::Init() {

  static ClassDocumentation<MEQQ2qq> documentation
    ("The ThePEG::MEQQ2qq class describes the standard QCD "
     "\\f$q\\bar{q} \\rightarrow q' + \\bar{q}'\\f$ matrix element.");
}

