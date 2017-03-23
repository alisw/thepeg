// -*- C++ -*-
//
// MEQq2Qq4.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEQq2Qq class.
//

#include "MEQq2Qq4.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace ThePEG;

void MEQq2Qq::getDiagrams() const {
  tcPDPtr g = getParticleData(ParticleID::g);
  for ( int i = 1; i <= maxFlavour(); ++i ) {
    tcPDPtr q = getParticleData(i);
    tcPDPtr qb = q->CC();
    for ( int j = 1; j <= maxFlavour(); ++j ) {
      if ( i == j ) continue;
      tcPDPtr qp = getParticleData(j);
      tcPDPtr qbp = qp->CC();
      add(new_ptr((Tree2toNDiagram(3), q, g, qp, 1, q, 2, qp, -1)));
      add(new_ptr((Tree2toNDiagram(3), qb, g, qp, 1, qb, 2, qp, -2)));
      add(new_ptr((Tree2toNDiagram(3), q, g, qbp, 1, q, 2, qbp, -3)));
      add(new_ptr((Tree2toNDiagram(3), qb, g, qbp, 1, qb, 2, qbp, -4)));
    }
  }
}

double MEQq2Qq::me2() const {
  return comfac()*colA()*KfacA()*2.0/9.0;
}

Selector<const ColourLines *>
MEQq2Qq::colourGeometries(tcDiagPtr diag) const {

  static ColourLines ctUT("1 -2 5, 3 2 4");
  static ColourLines ctST("3 2 -1, -4 -2 5");
  static ColourLines ctTS("1 -2 -3, -5 2 4");
  static ColourLines ctTU("-4 -2 -3, -5 2 -1");

  Selector<const ColourLines *> sel;
  if ( diag->id() == -1 )
    sel.insert(1.0, &ctUT);
  else if ( diag->id() == -2 )
    sel.insert(1.0, &ctST);
  else if ( diag->id() == -3 )
    sel.insert(1.0, &ctTS);
  else if ( diag->id() == -4 )
    sel.insert(1.0, &ctTU);
  return sel;
}

Selector<MEQq2Qq::DiagramIndex>
MEQq2Qq::diagrams(const DiagramVector & diags) const {
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) 
    if ( diags[i]->id() == -1 || diags[i]->id() == -2 ||
	 diags[i]->id() == -3 || diags[i]->id() == -4 ) sel.insert(1.0, i);
  return sel;
}

NoPIOClassDescription<MEQq2Qq> MEQq2Qq::initMEQq2Qq;
// Definition of the static class description member.

void MEQq2Qq::Init() {

  static ClassDocumentation<MEQq2Qq> documentation
    ("The ThePEG::MEQq2Qq class describes the standard QCD "
     "\\f$qq' \\rightarrow qq'\\f$ matrix element.");
}

