// -*- C++ -*-
//
// MENCDIS.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MENCDIS class.
//

#include "MENCDIS.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/StandardModel/StandardModelBase.h"
#include "ThePEG/Handlers/StandardXComb.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;

MENCDIS::MENCDIS()
  : mZ2(ZERO) {}

MENCDIS::MENCDIS(const MENCDIS & x)
  : ME2to2QCD(x), mZ2(x.mZ2) {}

MENCDIS::~MENCDIS() {}

unsigned int MENCDIS::orderInAlphaS() const {
  return 0;
}

unsigned int MENCDIS::orderInAlphaEW() const {
  return 2;
}

void MENCDIS::getDiagrams() const {
  tcPDPtr gamma = getParticleData(ParticleID::gamma);
  tcPDPtr Z0 = getParticleData(ParticleID::Z0);
  tcPDPtr ep = getParticleData(ParticleID::eplus);
  tcPDPtr em = getParticleData(ParticleID::eminus);
  for ( int i = -maxFlavour(); i <= maxFlavour(); ++i ) {
    if ( !i ) continue;
    tcPDPtr q = getParticleData(i);
    add(new_ptr((Tree2toNDiagram(3), q, gamma, em, 1, q, 2, em, -1)));
    add(new_ptr((Tree2toNDiagram(3), q, Z0, em, 1, q, 2, em, -2)));
    add(new_ptr((Tree2toNDiagram(3), q, gamma, ep, 1, q, 2, ep, -1)));
    add(new_ptr((Tree2toNDiagram(3), q, Z0, ep, 1, q, 2, ep, -2)));
  }
}

Energy2 MENCDIS::scale() const {
  return -tHat();
}

double MENCDIS::me2() const {
  double lastG = 0.0;
  double lastIntr = 0.0;
  double lastZ = 0.0;

  //pq is a vector in the same direction as the quark with zero mass.
  Lorentz5Momentum pq = meMomenta()[0];
  pq.setMass(ZERO);
  pq.rescaleEnergy();
  double y = 1.0 - pq.dot(meMomenta()[3]) / pq.dot(meMomenta()[1]);
  Energy4 F2Coeff = sqr(sHat()) * (1 + sqr(1-y));
  Energy4 F3Coeff = sqr(sHat()) * (1 - sqr(1-y));
  double C = 16 * SM().sin2ThetaW() * ( 1.0 - SM().sin2ThetaW() );
  if(mePartonData()[0]->id() < 0){
    F3Coeff = -F3Coeff;
  }
  if(mePartonData()[1]->id() < 0){
    F3Coeff = -F3Coeff;
  }
  if( abs(mePartonData()[0]->id())%2 == 0 ){
    lastG = F2Coeff * sqr(SM().eu()) / sqr(tHat());
    lastIntr = -2*SM().eu()*(F2Coeff*SM().ve()*SM().vu() + 
      2*F3Coeff*SM().ae()*SM().au()) / (-tHat() * (-tHat() + mZ2) * C);
    lastZ = ( F2Coeff * (sqr(SM().ae())+sqr(SM().ve())) * 
      (sqr(SM().au())+sqr(SM().vu())) + 4.0 * F3Coeff*SM().ve()*
      SM().ae()*SM().au()*SM().vu() ) / sqr((-tHat() + mZ2) * C);
  }
  else{
    lastG = F2Coeff * sqr(SM().ed()) / sqr(tHat());
    lastIntr = -2*SM().ed()*(F2Coeff*SM().ve()*SM().vd() +
      2*F3Coeff*SM().ae()*SM().ad()) / (-tHat() * (-tHat() + mZ2) * C);
    lastZ = ( F2Coeff * (sqr(SM().ae())+sqr(SM().ve())) * 
      (sqr(SM().ad())+sqr(SM().vd())) + 4.0 * F3Coeff*SM().ve()*
      SM().ae()*SM().ad()*SM().vd() ) / sqr((-tHat() + mZ2) * C);
  }

  DVector save;
  meInfo(save << lastG << lastZ);
  return (lastG + lastIntr + lastZ) * sqr(SM().alphaEM(scale())) *
    32.0 * sqr(Constants::pi);
}

Selector<MENCDIS::DiagramIndex>
MENCDIS::diagrams(const DiagramVector & diags) const {
  if ( lastXCombPtr() ) {
    lastG = meInfo()[0];
    lastZ = meInfo()[1];
  }
  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) {
    if ( diags[i]->id() == -1 ) sel.insert(lastG, i);
    else if ( diags[i]->id() == -2 ) sel.insert(lastZ, i);
  }
  return sel;
}

Selector<const ColourLines *>
MENCDIS::colourGeometries(tcDiagPtr diag) const {

  static ColourLines c("1 4");
  static ColourLines cb("-1 -4");

  Selector<const ColourLines *> sel;
  if ( diag->partons()[0]->id() > 0 )
    sel.insert(1.0, &c);
  else
    sel.insert(1.0, &cb);
  return sel;
}

IBPtr MENCDIS::clone() const {
  return new_ptr(*this);
}

IBPtr MENCDIS::fullclone() const {
  return new_ptr(*this);
}

void MENCDIS::doinit() {
  ME2to2QCD::doinit();
  tcPDPtr Z0 = getParticleData(ParticleID::Z0);
  mZ2 = sqr(Z0->mass());
}

void MENCDIS::persistentOutput(PersistentOStream & os) const {
  os << ounit(mZ2, GeV2) << lastG << lastZ;
}

void MENCDIS::persistentInput(PersistentIStream & is, int) {
  is >> iunit(mZ2, GeV2) >> lastG >> lastZ;
}

ClassDescription<MENCDIS> MENCDIS::initMENCDIS;

void MENCDIS::Init() {

  static ClassDocumentation<MENCDIS> documentation
    ("The ThePEG::MENCDIS class implements the full"
     "\\f$e^\\pm q \\rightarrow e^\\pm q\\f$ "
     "matrix element including the interference terms.");

}

