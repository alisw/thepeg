// -*- C++ -*-
//
// QuarksToHadronsDecayer.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the QuarksToHadronsDecayer class.
//

#include "QuarksToHadronsDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;

QuarksToHadronsDecayer::~QuarksToHadronsDecayer() {}

IBPtr QuarksToHadronsDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr QuarksToHadronsDecayer::fullclone() const {
  return new_ptr(*this);
}

bool QuarksToHadronsDecayer::accept(const DecayMode & dm) const {
  int col = 0;
  int acol = 0;
  if ( !dm.productMatchers().empty() ) {
    for ( MatcherMSet::const_iterator it = dm.productMatchers().begin();
	  it != dm.productMatchers().end(); ++it ) {
      if ( typeid(**it) == typeid(MatchLightQuark) ) ++col;
      else if ( typeid(**it) == typeid(MatchLightAntiQuark) ) ++acol;
      else return false;
    }
    if ( col != 1 || col != acol ) return false;
  }
  if ( dm.orderedProducts().size() + col + acol < 2 ||
       !dm.cascadeProducts().empty() || dm.wildProductMatcher() ) return false;
  for ( int i = 0, N = dm.orderedProducts().size(); i < N; ++i ) {
    if ( DiquarkMatcher::Check(*dm.orderedProducts()[i]) ) {
      if ( i + 1 != N ) return false;
      if ( dm.orderedProducts()[i]->id() < 0 ) ++col;
      else ++acol;
    }
    if ( QuarkMatcher::Check(*dm.orderedProducts()[i]) ) {
      if ( dm.orderedProducts()[i]->id() > 0 ) ++col;
      else ++acol;
    }
  }
  if ( acol != col || col < 1 || col > 2 ) return false;
  return true;
}

PVector QuarksToHadronsDecayer::decay(const DecayMode & dm,
				      const Particle & parent) const {
  PVector children;
  tcPDVector quarks;
  if ( !dm.productMatchers().empty() ) {
    tcPDPtr pd = getParticleData(flavourGenerator()->selectQuark());
    quarks.push_back(pd);
    quarks.push_back(pd->CC());
  }
  Energy summq = ZERO;
  Energy summp = ZERO;
  tPDVector prods = dm.orderedProducts();
  for ( int i = 0, N = prods.size(); i < N; ++i )
    if ( QuarkMatcher::Check(*prods[i]) || DiquarkMatcher::Check(*prods[i])) {
      quarks.push_back(prods[i]);
      summq += quarks.back()->mass();
    } else {
      children.push_back(prods[i]->produceParticle());
      summp += children.back()->mass();
    }

  Energy summh = ZERO;
  PVector hadrons;

  if ( !quarks.empty() ) do {

    hadrons = getHadrons(getN(parent.mass(), summq, quarks.size()), quarks);

    summh = ZERO;
    for ( int i = 0, N = hadrons.size(); i < N; ++i )
      summh += hadrons[i]->mass();

  } while ( hadrons.empty() || summp + summh >= parent.mass() );

  children.insert(children.end(), hadrons.begin(), hadrons.end());

  distribute(parent, children);

  finalBoost(parent, children);
  setScales(parent, children);

  return children;

}

int QuarksToHadronsDecayer::getN(Energy m0, Energy summq, int Nq) const {
  int Nh = fixedN();
  if ( Nh >= 2 ) return Nh;

  double c = c1()*log((m0 - summq)/c2()) + c3();
  if ( c < 0.0 ) return minN();
  while ( true ) {
    using namespace Constants;
    Nh = int(0.5 + double(Nq)/4.0 + c +
	     sqrt(-2.0*c*log(max(1.0e-10, rnd())))*sin(2.0*pi*rnd()));
    if ( Nh >= minN() ) return Nh;
  }
}

PVector QuarksToHadronsDecayer::
getHadrons(int Nh, tcPDVector quarks) const {
  PVector hadrons;
  Nh -= quarks.size()/2;
  while ( Nh-- > 0 ) {
    int i = irnd(quarks.size() - 1);
    tcPDPair hq = flavourGenerator()->alwaysGenerateHadron(quarks[i]);
    hadrons.push_back(hq.first->produceParticle());
    quarks[i] = hq.second;
  }
  if ( DiquarkMatcher::Check(*quarks[0]) && DiquarkMatcher::Check(*quarks[1]) )
    return PVector();
  tcPDPtr h = flavourGenerator()->alwaysGetHadron(quarks[0], quarks[1]);
  hadrons.push_back(h->produceParticle());
  if ( quarks.size() <= 2 ) return hadrons;
  if ( DiquarkMatcher::Check(*quarks[2]) && DiquarkMatcher::Check(*quarks[3]) )
    return PVector();
  h = flavourGenerator()->alwaysGetHadron(quarks[2], quarks[3]);
  hadrons.push_back(h->produceParticle());
  return hadrons;
}

void QuarksToHadronsDecayer::
distribute(const Particle & parent, PVector & children) const {
  do {
    try {
      SimplePhaseSpace::CMSn(children, parent.mass());
    }
    catch ( ImpossibleKinematics ) {
      children.clear();
      return;
    }
  } while ( reweight(parent, children) < rnd() );
}

double QuarksToHadronsDecayer::
reweight(const Particle &, const PVector &) const {
  return 1.0;
}

void QuarksToHadronsDecayer::persistentOutput(PersistentOStream & os) const {
  os << theFixedN << theMinN << theC1 << ounit(theC2,GeV) << theC3 << theFlavourGenerator;
}

void QuarksToHadronsDecayer::persistentInput(PersistentIStream & is, int) {
  is >> theFixedN >> theMinN >> theC1 >> iunit(theC2,GeV) >> theC3 >> theFlavourGenerator;
}

ClassDescription<QuarksToHadronsDecayer> QuarksToHadronsDecayer::initQuarksToHadronsDecayer;
// Definition of the static class description member.

void QuarksToHadronsDecayer::Init() {

  static ClassDocumentation<QuarksToHadronsDecayer> documentation
    ("This class decays particles to nq (2 or 4) quarks which then are "
     "decayes to hadrons according to phase space. The number of final "
     "hadrons can either be given by a fixed number or as a Gaussian "
     "multiplicity distribution centered around c+nq/4+c3 and a width "
     "sqrt(c), where c = c1 log((m - summ)/c2), m is the mass of the "
     "decaying particle, summ the sum of the quark masses and ci real "
     "parameters.");

  static Parameter<QuarksToHadronsDecayer,int> interfaceFixedN
    ("FixedN",
     "The fixed number of hadrons to be produced. If less than 2, the "
     "number is instead given by a gaussian multiplicity distribution.",
     &QuarksToHadronsDecayer::theFixedN, 0, 0, 10,
     true, false, true);

  static Parameter<QuarksToHadronsDecayer,int> interfaceMinN
    ("MinN",
     "The minimum hadrons to be produced.",
     &QuarksToHadronsDecayer::theMinN, 2, 2, 10,
     true, false, true);

  static Parameter<QuarksToHadronsDecayer,double> interfaceC1
    ("C1",
     "The c1 parameter of the gaussian multiplicity distribution centered "
     "around c1 log((m - summ)/c2) +c3.",
     &QuarksToHadronsDecayer::theC1, 4.5, 0.0, 10.0,
     true, false, true);

  static Parameter<QuarksToHadronsDecayer,Energy> interfaceC2
    ("C2",
     "The c2 parameter of the gaussian multiplicity distribution centered "
     "around c1 log((m - summ)/c2) +c3.",
     &QuarksToHadronsDecayer::theC2, GeV, 0.7*GeV, ZERO, 10.0*GeV,
     true, false, true);

  static Parameter<QuarksToHadronsDecayer,double> interfaceC3
    ("C3",
     "The c3 parameter of the gaussian multiplicity distribution centered "
     "around c1 log((m - summ)/c2) +c3.",
     &QuarksToHadronsDecayer::theC3, 0.0, 0.0, 10.0,
     true, false, true);

  static Reference<QuarksToHadronsDecayer,FlavourGenerator>
    interfaceFlavourGenerator
    ("FlavourGenerator",
     "The object in charge of generating hadrons spieces from given quark "
     "flavours.",
     &QuarksToHadronsDecayer::theFlavourGenerator,
     true, false, true, false, true);

  interfaceFixedN.rank(10);
  interfaceMinN.rank(9);
  interfaceFlavourGenerator.rank(8);
  interfaceMinN.setHasDefault(false);;

}

