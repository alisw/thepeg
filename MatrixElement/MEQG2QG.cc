// -*- C++ -*-
//
// MEQG2QG.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEQG2QG class.
//

#include "MEQG2QG.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Utilities/Exception.h"

using namespace ThePEG;

IBPtr MEQG2QG::clone() const {
  return new_ptr(*this);
}

IBPtr MEQG2QG::fullclone() const {
  return new_ptr(*this);
}

void MEQG2QG::getDiagrams() const {
  tcPDPtr g = getParticleData(ParticleID::g);
  for ( int i = -maxFlavour(); i <= maxFlavour(); ++i ) {
    if ( i == 0 ) continue;
    tcPDPtr q = getParticleData(i);
    tcPDPtr qb = q->CC();
    add(new_ptr((Tree2toNDiagram(3), q, g, g, 1, q, 2, g, -1)));
    add(new_ptr((Tree2toNDiagram(3), q, q, g, 2, q, 1, g, -2)));
    add(new_ptr((Tree2toNDiagram(2), q, g, 1, q, 3, q, 3, g, -3)));
  }
}

double MEQG2QG::me2() const {
  return 2.0*comfac()*((colA1() + colA2())*KfacA() +
		       (colB1() + colB2())*Kfac())/9.0;
}

Selector<const ColourLines *>
MEQG2QG::colourGeometries(tcDiagPtr diag) const {

  static ColourLines ctST("1 -2 -3, 3 5, -5 2 4");
  static ColourLines ctTS("-4 -2 5, -5 -3, 3 2 -1");
  static ColourLines ctUT("1 -2 5, -5 -3, 3 4");
  static ColourLines ctTU("-4 -3, 3 5, -5 2 -1");
  static ColourLines cuTU("1 5, -5 2 -3, 3 4");
  static ColourLines cuUT("-4 -3, 3 -2 5, -1 -5");
  static ColourLines csST("1 -2, 2 3 5, -5 4");
  static ColourLines csTS("-4 5, -5 -3 -2, 2 -1");

  Selector<const ColourLines *> sel;
  int q = diag->partons()[0]->id();
  if ( diag->id() == -1 ) {
    if ( q > 0 ) {
      sel.insert(colA1(), &ctST);
      sel.insert(colB1(), &ctUT);
    } else {
      sel.insert(colA1(), &ctTS);
      sel.insert(colB1(), &ctTU);
    }
  } else if ( diag->id() == -2 ) {
    if ( q > 0 ) sel.insert(1.0, &cuTU);
    else sel.insert(1.0, &cuUT);
  } else {
    if ( q > 0 ) sel.insert(1.0, &csST);
    else sel.insert(1.0, &csTS);
  }
  return sel;
}

Selector<MEQG2QG::DiagramIndex>
MEQG2QG::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) 
    if ( diags[i]->id() == -1 ) sel.insert(colA1() + colB1(), i);
    else if ( diags[i]->id() == -2 )  sel.insert(colB2(), i);
    else sel.insert(colA2(), i);
  return sel;
}

NoPIOClassDescription<MEQG2QG> MEQG2QG::initMEQG2QG;
// Definition of the static class description member.

void MEQG2QG::Init() {

  static ClassDocumentation<MEQG2QG> documentation
    ("The ThePEG::MEQG2QG class describes the standard QCD "
     "\\f$qg\\rightarrow qg\\f$ matrix element.");

}

