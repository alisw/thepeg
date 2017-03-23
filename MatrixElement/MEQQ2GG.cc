// -*- C++ -*-
//
// MEQQ2GG.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEQQ2GG class.
//

#include "MEQQ2GG.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace ThePEG;

IBPtr MEQQ2GG::clone() const {
  return new_ptr(*this);
}

IBPtr MEQQ2GG::fullclone() const {
  return new_ptr(*this);
}

void MEQQ2GG::getDiagrams() const {
  tcPDPtr g = getParticleData(ParticleID::g);
  for ( int i = 1; i <= maxFlavour(); ++i ) {
    tcPDPtr q = getParticleData(i);
    tcPDPtr qb = q->CC();
    add(new_ptr((Tree2toNDiagram(3), q, q, qb, 1, g, 2, g, -1)));
    add(new_ptr((Tree2toNDiagram(3), q, q, qb, 2, g, 1, g, -2)));
  }
}

double MEQQ2GG::me2() const {
  return comfac()*(colA() + colB())*KfacA()*16.0/27.0;
}

Selector<const ColourLines *>
MEQQ2GG::colourGeometries(tcDiagPtr diag) const {

  static ColourLines ctST("1 4, -4 2 5, -5 -3");
  static ColourLines ctSU("1 5, -5 2 4, -4 -3");

  Selector<const ColourLines *> sel;
  if ( diag->id() == -1 )
    sel.insert(1.0, &ctST);
  else
    sel.insert(1.0, &ctSU);
  return sel;
}

Selector<MEQQ2GG::DiagramIndex>
MEQQ2GG::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) 
    if ( diags[i]->id() == -1 ) sel.insert(colA(), i);
    else if ( diags[i]->id() == -2 )  sel.insert(colB(), i);
  return sel;
}

NoPIOClassDescription<MEQQ2GG> MEQQ2GG::initMEQQ2GG;
// Definition of the static class description member.

void MEQQ2GG::Init() {

  static ClassDocumentation<MEQQ2GG> documentation
    ("The ThePEG::MEQQ2GG class describes the standard QCD "
     "\\f$q\\bar{q} \\rightarrow gg\\f$ matrix element.");

}

