// -*- C++ -*-
//
// MEGG2GG.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEGG2GG class.
//

#include "MEGG2GG.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace ThePEG;

IBPtr MEGG2GG::clone() const {
  return new_ptr(*this);
}

IBPtr MEGG2GG::fullclone() const {
  return new_ptr(*this);
}

void MEGG2GG::getDiagrams() const {
  tcPDPtr g = getParticleData(ParticleID::g);
  add(new_ptr((Tree2toNDiagram(3), g, g, g, 1, g, 2, g, -1)));
  add(new_ptr((Tree2toNDiagram(3), g, g, g, 2, g, 1, g, -2)));
  add(new_ptr((Tree2toNDiagram(2), g, g, 1, g, 3, g, 3, g, -3)));
}

double MEGG2GG::me2() const {
  return 9.0*comfac()*((colA1() + colA2() + colB1() + colB2())*KfacA() +
			 (colC1() + colC2())*Kfac())/16.0;
}

Selector<const ColourLines *>
MEGG2GG::colourGeometries(tcDiagPtr diag) const {
  static ColourLines ctST("1 -2 -3, 3 5, -5 2 4, -4 -1");
  static ColourLines ctTS("1 4, -4 -2 5, -5 -3, 3 2 -1");
  static ColourLines ctUT("1 -2 5, -5 -3, 3 2 4, -4 -1");
  static ColourLines ctTU("1 4, -4 -2 -3, 3 5, -5 2 -1");

  static ColourLines cuSU("1 -2 -3, 3 4, -4 2 5, -5 -1");
  static ColourLines cuUS("1 5, -5 -2 4, -4 -3, 3 2 -1");
  static ColourLines cuTU("1 -2 4, -4 -3, 3 2 5, -5 -1");
  static ColourLines cuUT("1 5, -5 -2 -3, 3 4, -4 2 -1");

  static ColourLines csTS("1 3 4, -4 5, -5 -3 -2, 2 -1");
  static ColourLines csST("1 -2, 2 3 5, -5 4, -4 -3 -1");
  static ColourLines csUS("1 3 5, -5 4, -4 -3 -2, 2 -1");
  static ColourLines csSU("1 -2, 2 3 4, -4 5, -5 -3 -1");

  
  Selector<const ColourLines *> sel;
  if ( diag->id() == -1 ) {
    sel.insert(colA1(), &ctST);
    sel.insert(colA1(), &ctTS);
    sel.insert(colC2(), &ctUT);
    sel.insert(colC2(), &ctTU);
  } else if ( diag->id() == -2 ) {
    sel.insert(colB2(), &cuSU);
    sel.insert(colB2(), &cuUS);
    sel.insert(colC1(), &cuTU);
    sel.insert(colC1(), &cuUT);
  } else {
    sel.insert(colA2(), &csST);
    sel.insert(colA2(), &csTS);
    sel.insert(colB1(), &csSU);
    sel.insert(colB1(), &csUS);
  }
  return sel;
}

Selector<MEGG2GG::DiagramIndex>
MEGG2GG::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i )
    if ( diags[i]->id() == -1 ) sel.insert(colA1() + colC2(), i);
    else if ( diags[i]->id() == -2 ) sel.insert(colC1() + colB2(), i);
    else sel.insert(colB1() + colA2(), i);
  return sel;
}

NoPIOClassDescription<MEGG2GG> MEGG2GG::initMEGG2GG;
// Definition of the static class description member.

void MEGG2GG::Init() {

  static ClassDocumentation<MEGG2GG> documentation
    ("The ThePEG::MEGG2GG class describes the standard QCD "
     "\\f$gg \\rightarrow gg\\f$ matrix element.");

}

