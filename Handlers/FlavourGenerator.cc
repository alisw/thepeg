// -*- C++ -*-
//
// FlavourGenerator.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FlavourGenerator class.
//

#include "FlavourGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/ParticleData.h"

using namespace ThePEG;

tcPDPair FlavourGenerator::alwaysGenerateHadron(tcPDPtr q) const {
  tcPDPair hq = generateHadron(q);
  if ( !hq.first || ! hq.second ) throw FlavourGeneratorException()
    << "Flavour generator '" << name() << "' was not able to generate a "
    << "hadron from the flavour " << q->PDGName() << "."
    << Exception::runerror;
  return hq;
}


tcPDPtr FlavourGenerator::getBaryon(tcPDPtr q1, tcPDPtr q2, tcPDPtr q3) const {
  return getBaryon(q1->id(), q2->id(), q3->id());
}

tcPDPtr FlavourGenerator::getBaryon(long iq1, long iq2, long iq3) const {
  return getBaryon(getParticleData(iq1),
		   getParticleData(iq2),
		   getParticleData(iq3));
}

tcPDPtr FlavourGenerator::
alwaysGetBaryon(tcPDPtr q1, tcPDPtr q2, tcPDPtr q3) const {
  tcPDPtr ret = getBaryon(q1, q2, q3);
  if ( !ret ) throw FlavourGeneratorException()
    << "Flavour generator '" << name() << "' was not able to get a "
    << "baryon from the flavours " << q1->PDGName() << "," << q2->PDGName()
    << " and " << q3->PDGName() << "." << Exception::runerror;
  return ret;
}

tcPDPtr FlavourGenerator::alwaysGetBaryon(long iq1, long iq2, long iq3) const {
  tcPDPtr ret = getBaryon(iq1, iq2, iq3);
  if ( !ret ) throw FlavourGeneratorException()
    << "Flavour generator '" << name() << "' was not able to get a "
    << "baryon from the flavours " << iq1 << "," << iq2 << " and " << iq3
    << "." << Exception::runerror;
  return ret;
}

tcPDPtr FlavourGenerator::getHadron(tcPDPtr q1, tcPDPtr q2) const {
  return getHadron(q1->id(), q2->id());
}

tcPDPtr FlavourGenerator::getHadron(long iq1, long iq2) const {
  return getHadron(getParticleData(iq1), getParticleData(iq2));
}

tcPDPtr FlavourGenerator::alwaysGetHadron(tcPDPtr q1, tcPDPtr q2) const {
  tcPDPtr ret = getHadron(q1, q2);
  if ( !ret ) throw FlavourGeneratorException()
    << "Flavour generator '" << name() << "' was not able to get a "
    << "hadron from the flavours " << q1->PDGName() << " and "
    << q2->PDGName() << "." << Exception::runerror;
  return ret;
}

tcPDPtr FlavourGenerator::alwaysGetHadron(long iq1, long iq2) const {
  tcPDPtr ret = getHadron(iq1, iq2);
  if ( !ret ) throw FlavourGeneratorException()
    << "Flavour generator '" << name() << "' was not able to get a "
    << "hadron from the flavours " << iq1 << " and " << iq2 << "."
    << Exception::runerror;
  return ret;
}

AbstractNoPIOClassDescription<FlavourGenerator>
 FlavourGenerator::initFlavourGenerator;

void FlavourGenerator::Init() {

  static ClassDocumentation<FlavourGenerator> documentation
    ("This class is used to generate hadron types using a given flavour "
     "content.");

}

