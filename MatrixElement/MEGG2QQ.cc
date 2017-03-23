// -*- C++ -*-
//
// MEGG2QQ.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEGG2QQ class.
//

#include "MEGG2QQ.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace ThePEG;

IBPtr MEGG2QQ::clone() const {
  return new_ptr(*this);
}

IBPtr MEGG2QQ::fullclone() const {
  return new_ptr(*this);
}

void MEGG2QQ::getDiagrams() const {
  tcPDPtr g = getParticleData(ParticleID::g);
  for ( int i = 1; i <= maxFlavour(); ++i ) {
    tcPDPtr q = getParticleData(i);
    tcPDPtr qb = q->CC();
    add(new_ptr((Tree2toNDiagram(3), g, q, g, 1, q, 2, qb, -1)));
    add(new_ptr((Tree2toNDiagram(3), g, qb, g, 2, q, 1, qb, -2)));
  }
}

double MEGG2QQ::me2() const {
  return comfac()*(colA() + colB())*KfacA()/12.0;
}

Selector<const ColourLines *>
MEGG2QQ::colourGeometries(tcDiagPtr diag) const {

  static ColourLines ctST("1 4, -5 -3, 3 2 -1");
  static ColourLines cuSU("1 -2 -3, 3 4, -5 -1");

  Selector<const ColourLines *> sel;
  if ( diag->id() == -1 )
    sel.insert(1.0, &ctST);
  else
    sel.insert(1.0, &cuSU);
  return sel;
}

Selector<MEGG2QQ::DiagramIndex>
MEGG2QQ::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) 
    if ( diags[i]->id() == -1 ) sel.insert(colA(), i);
    else if ( diags[i]->id() == -2 )  sel.insert(colB(), i);
  return sel;
}

NoPIOClassDescription<MEGG2QQ> MEGG2QQ::initMEGG2QQ;
// Definition of the static class description member.

void MEGG2QQ::Init() {

  static ClassDocumentation<MEGG2QQ> documentation
    ("The ThePEG::MEGG2QQ class describes the standard QCD "
     "\\f$gg\\rightarrow q\\bar{q}\\f$ matrix element.");

}

