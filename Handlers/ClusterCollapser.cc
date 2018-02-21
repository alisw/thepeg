// -*- C++ -*-
//
// ClusterCollapser.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ClusterCollapser class.
//

#include "ClusterCollapser.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/Step.h"
#include "ThePEG/Utilities/UtilityBase.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "ThePEG/Utilities/EnumIO.h"
#include "ThePEG/Utilities/Throw.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ClusterCollapser.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;

ClusterCollapser::~ClusterCollapser() {}

IBPtr ClusterCollapser::clone() const {
  return new_ptr(*this);
}

IBPtr ClusterCollapser::fullclone() const {
  return new_ptr(*this);
}

void ClusterCollapser::
handle(EventHandler &, const tPVector & tagged,
       const Hint &) {
  collapse(tagged, newStep());
}

vector<ColourSinglet> ClusterCollapser::
collapse(tPVector tagged, tStepPtr newstep) {
  vector<ColourSinglet> newTagged;
  SingletMap clusters = getSinglets(tagged);

  // Go through all clusters below the cut.
  while ( !clusters.empty() && clusters.begin()->first < cut() ) {

    SingletMap::iterator clit = clusters.begin();

    ColourSinglet & cl = clit->second;

    // If a cluster contains too many junctions, split them into
    // several singlets.
    while ( cl.nPieces() > 3 ) insert(clusters, cl.splitInternal());

    // If a cluster contains a junktion and a diquark, split the
    // diquark and make two simple strings.
    while ( diQuarkJunction(cl) )
      insert(clusters, splitDiQuarkJunction(cl, newstep));

    // First try to collapse into two particles.
    int ntry = nTry2();
    while ( ntry > 0 ) {
      if ( collapse2(newstep, cl) ) break;
      --ntry;
    }

    // If that didn't work collapse into one particle and shuffle some
    // energy to or from the other tagged particles.
    if ( ntry == 0 ) {

      // If this was a di-diquark cluster, split it into two.
      if ( diDiQuark(cl) ) {
	insert(clusters, splitDiDiQuark(cl, newstep));
	updateTagged(tagged);
      }
      collapse(newstep, cl, tagged);
    }

    updateTagged(tagged);

    // Remove the collapsed cluster.
    clusters.erase(clit);

    // Recalculate masses of the remaining clusters, insert them in a
    // temporary map and swap this map for the old map.
    multimap<Energy,ColourSinglet> newClusters;
    while ( !clusters.empty() ) {
      ColourSinglet & cl = clusters.begin()->second;
      for ( int i = 0, N = cl.partons().size(); i < N; ++i )
	cl.partons()[i] = cl.partons()[i]->final();
      insert(newClusters, cl);
      clusters.erase(clusters.begin());
    }
    clusters.swap(newClusters);

  }

  // Return the left-over clusters in a vector.
  newTagged.resize(clusters.size());
  for ( int is = 0, NS = newTagged.size(); is < NS; ++is ) {
    newTagged[is].swap(clusters.begin()->second);
    clusters.erase(clusters.begin());
  }
  return newTagged;

}

void ClusterCollapser::updateTagged(tPVector & tagged) const {
  tPVector::iterator it = tagged.begin();
  set<tPPtr> children;
  while ( it != tagged.end() ) {
    *it = (**it).final();
    if ( (**it).decayed() ) {
      children.insert((**it).children().begin(), (**it).children().end());
      it = tagged.erase(it);
    }
    else ++it;
  }
  tagged.insert(tagged.end(), children.begin(), children.end());

}

Energy ClusterCollapser::mass(const ColourSinglet & cl) {
  LorentzMomentum sump;
  Energy summ = ZERO;
  for ( int i = 0, N = cl.partons().size(); i < N; ++i ) {
    summ += cl.parton(i)->data().constituentMass();
    sump += cl.parton(i)->momentum();
  }
  return sump.m() - summ;
}

void ClusterCollapser::insert(SingletMap & mmap, const ColourSinglet & cl) {
  mmap.insert(make_pair(mass(cl), cl));
}

bool ClusterCollapser::diDiQuark(const ColourSinglet & cs) {
  return ( cs.nPieces() == 1 &&
	   DiquarkMatcher::Check(cs.piece(1).front()->data()) &&
	   DiquarkMatcher::Check(cs.piece(1).back()->data()) );
}

ColourSinglet ClusterCollapser::
splitDiDiQuark(ColourSinglet & cs, tStepPtr newStep) const {
  ColourSinglet ret;

  // Split the first diquark
  tcPPtr diq = cs.piece(1).front();
  PPair qq1;
  qq1.first = getParticle(diq->id()/1000);
  qq1.second = getParticle((diq->id()/100)%10);
  if ( qq1.first->mass() + qq1.second->mass() >= diq->mass() ) {
    // If the sum of the quarks masses is larger than the diuark mass,
    // set the new quark masses to zero.
    qq1.first->set5Momentum(Lorentz5Momentum());
    qq1.second->set5Momentum(Lorentz5Momentum());
  }

  // Distribut the quarks evenly in the cms of the diquark and add
  // them as children of the diquark to the new step.
  SimplePhaseSpace::CMS(sqr(diq->mass()), qq1.first, qq1.second);
  qq1.first->boost(diq->momentum().boostVector());
  qq1.second->boost(diq->momentum().boostVector());
  newStep->addDecayProduct(diq, qq1.first);
  newStep->addDecayProduct(diq, qq1.second);

  // Split the second diquark
  diq = cs.piece(1).back();
  PPair qq2;
  qq2.first = getParticle(diq->id()/1000);
  qq2.second = getParticle((diq->id()/100)%10);
  if ( qq2.first->mass() + qq2.second->mass() >= diq->mass() ) {
    // If the sum of the quarks masses is larger than the diuark mass,
    // set the new quark masses to zero.
    qq2.first->set5Momentum(Lorentz5Momentum());
    qq2.second->set5Momentum(Lorentz5Momentum());
  }
  
    // Distribut the quarks evenly in the cms of the diquark and add
  // them as children of the diquark to the new step.
  SimplePhaseSpace::CMS(sqr(diq->mass()), qq2.first, qq2.second);
  qq2.first->boost(diq->momentum().boostVector());
  qq2.second->boost(diq->momentum().boostVector());
  newStep->addDecayProduct(diq, qq2.first);
  newStep->addDecayProduct(diq, qq2.second);

  if ( rndbool() ) swap(qq1.first, qq1.second);

  return ret = cs.splitDiDiQuark(qq1, qq2);

}


bool ClusterCollapser::diQuarkJunction(const ColourSinglet & cs) {
  if ( cs.nPieces() < 3 ) return false;
  for ( int i = 1, N = cs.nPieces(); i <= N; ++i )
    if ( DiquarkMatcher::Check(cs.piece(i).front()->data()) ||
	 DiquarkMatcher::Check(cs.piece(i).back()->data()) ) return true;
  return false;
}

ColourSinglet ClusterCollapser::
splitDiQuarkJunction(ColourSinglet & cs, tStepPtr newStep) const {
  ColourSinglet ret;

  // Find the diquarks in the singlet and pick one randomly.
  vector<ColourSinglet::Index> diqs;
  for ( ColourSinglet::Index i = 1, N = cs.nPieces(); i <= N; ++i ) {
    if ( DiquarkMatcher::Check(cs.piece(i).front()->data()) )
      diqs.push_back(-i);
    if ( DiquarkMatcher::Check(cs.piece(i).back()->data()) )
      diqs.push_back(i);
  }
  if ( diqs.empty() ) return ret;
  ColourSinglet::Index seli = diqs[UseRandom::irnd(diqs.size())];
  tcPPtr diq = seli > 0? cs.piece(seli).back(): cs.piece(seli).front();

  // Create the to quarks
  PPair qq;
  qq.first = getParticle(diq->id()/1000);
  qq.second = getParticle((diq->id()/100)%10);
  if ( qq.first->mass() + qq.second->mass() >= diq->mass() ) {
    // If the sum of the quarks masses is larger than the diuark mass,
    // set the new quark masses to zero.
    qq.first->set5Momentum(Lorentz5Momentum());
    qq.second->set5Momentum(Lorentz5Momentum());
  }

  // Distribut the quarks evenly in the cms of the diquark and add
  // them as children of the diquark to the new step.
  SimplePhaseSpace::CMS(sqr(diq->mass()), qq.first, qq.second);
  qq.first->boost(diq->momentum().boostVector());
  qq.second->boost(diq->momentum().boostVector());
  newStep->addDecayProduct(diq, qq.first);
  newStep->addDecayProduct(diq, qq.second);

  ret = cs.splitDiQuarkJunction(seli, diq, qq);

  return ret;

}

ClusterCollapser::SingletMap
ClusterCollapser::getSinglets(const tPVector & pv) const {
  SingletMap ret;

  // Get initial singlets
  vector<ColourSinglet> clus = ColourSinglet::getSinglets(pv.begin(), pv.end());

  // Return the singlets ordered in mass.
  for ( int i = 0, N = clus.size(); i < N; ++i )
    if ( !clus[i].partons().empty() ) insert(ret, clus[i]);
  return ret;
}

tPVector ClusterCollapser::
getCompensators(Energy mh, const ColourSinglet & cs,
		const tPVector & tagged, tStepPtr newStep) const {
  tPVector ret;
  tcPVector comp;
  // First find the particles which are not a part of the collapsing
  // cluster.
  tParticleSet compset;
  for ( int i = 0, N = tagged.size(); i < N; ++i )
    if ( !member(cs.partons(), tagged[i]) ) compset.insert(tagged[i]);

  LorentzMomentum pcomp;
  LorentzMomentum pc = cs.momentum();
  // start by only looking at other strings.
  bool alsoSinglets = false;
  do {
    // Return an empty vector if no particles left.
    if ( compset.empty() ) break;

    // Now find the particle which is closest in phase space to the
    // cluster.
    Energy2 dist = Constants::MaxEnergy2;
    tParticleSet::iterator sel = compset.end();
    for ( tParticleSet::iterator it = compset.begin();
	  it != compset.end(); ++it ) {
      if ( !(**it).coloured() && !alsoSinglets ) continue;
      if ( -(pc - (**it).momentum()).m2() < dist ) {
	dist = -(pc - (**it).momentum()).m2();
	sel = it;
      }
    }

    if ( sel == compset.end() ) {
      if ( alsoSinglets ) break;
      else {
	alsoSinglets = true;
	continue;
      }
    }

    // Add to the temporary vector.
    comp.push_back(*sel);
    pcomp += (**sel).momentum();
    compset.erase(sel);

    // If there was not enough energy, find an additional compensator
    // particle. Also check that compensators have mass to avoid boost
    // problems.
  } while ( comp.empty() || (pc + pcomp).m() <= mh + pcomp.m() ||
	    ( comp.size() > 1 && pcomp.m2() <= ZERO ) );


  // If this didn't work, let's try to fix it by disregarding the
  // closest particle.
  tcPVector::size_type end = comp.size();
  while ( (pc + pcomp).m() <= mh + pcomp.m() ||
	  ( comp.size() > 1 && pcomp.m2() <= ZERO ) ) {
    if ( end == comp.size() ) {
      if ( comp.size() < 2 ) return ret;
      comp.erase(comp.begin());
      end = 1;
    } else
      ++end;
    pcomp = Utilities::sumMomentum(comp.begin(), comp.begin() + end);

  }


  // Now copy the compensators, add them to the new set and return them.
  ret.resize(end);
  for ( tcPVector::size_type i = 0; i < end; ++i )
    ret[i] = newStep->copyParticle(comp[i]);
  return ret;
}

void ClusterCollapser::
collapse(tStepPtr newStep, const ColourSinglet & cs,
	 const tPVector & tagged) const {

  // Produce the hadron to collapse into, set its momentum to the one
  // of the collapsing cluster and generate the mass to be set.
  tcPDPtr hd = getHadron(cs);
  LorentzMomentum pc = cs.momentum();
  PPtr h = hd->produceParticle(pc);
  Energy mh = hd->generateMass();

  // Select the partons to be used in momentum compensation
  tPVector comp = getCompensators(mh, cs, tagged, newStep);

  if ( comp.empty() ) {
    if ( errorlevel ) throw ClusterException(*this)
      << "Could not find particles to shuffle momentum."
      << errorlevel;
    h->set5Momentum(Lorentz5Momentum(pc, mh));
  } else {

    // Boost the hadron and the compensating particles into their cms
    // with the hadon along the z-axis.
    comp.push_back(h);
    LorentzRotation R =
      Utilities::boostToCM(comp.begin(), comp.end(), comp.end() - 1).inverse();

    // Give the hadron and the compensators their correct mass and
    // momentum and bost back
    LorentzMomentum pcomp =
      Utilities::sumMomentum(comp.begin(), comp.end() - 1);
    Energy2 s = (pcomp + h->momentum()).m2();

    try {    
      Energy pnew = SimplePhaseSpace::getMagnitude(s, mh, pcomp.m());
      h->set5Momentum(R*Lorentz5Momentum(ZERO, ZERO,
					 pnew, sqrt(sqr(pnew) + sqr(mh)), mh));
      
      comp.pop_back();
      R = R*LorentzRotation(0.0, 0.0, -(pcomp.e()*pcomp.z() +
					sqrt(sqr(pnew) + pcomp.m2())*pnew)/
			    (sqr(pnew) + sqr(pcomp.e())));
    }
    catch ( ... ) {
      Throw<Exception>()
	// *** TODO *** Check cicumstances for this.
	<< "Impossible kinematics found in Ariadne5::DipoleState::purgeGluon."
	<< Exception::eventerror;
    }    
    Utilities::transform(comp, R);
  }

  // Add the new particle to the step.
  if ( !newStep->
       addDecayProduct(cs.partons().begin(), cs.partons().end(), h) ) {
    throw ClusterException(*this)
      << "Could not add decay products to the new step" << Exception::abortnow;
  }

}

bool ClusterCollapser::
collapse2(tStepPtr newStep, const ColourSinglet & cs) const {

  // First get the two particles into which to decay the cluster.
  tcPDPair pdp = getHadrons(cs);
  PVector h(2);
  h[0] = pdp.first->produceParticle();
  h[1] = pdp.second->produceParticle();

  // Check the invariant mass of the cluster and return false if there
  // was not enough energy.
  LorentzMomentum pc = cs.momentum();
  Energy2 s = pc.m2();
  if ( sqr(h[0]->mass() + h[1]->mass()) >= s && diDiQuark(cs) ) {
    // In the special case of di-diquars we try to take the flavours
    // to create a meson pair instead. We begin by finding the quarks
    PDPair qq = make_pair(getParticleData(cs.piece(1).front()->id()/1000),
			 getParticleData((cs.piece(1).front()->id()/100)%10));
    PDPair aqq = make_pair(getParticleData(cs.piece(1).back()->id()/1000),
			  getParticleData((cs.piece(1).back()->id()/100)%10));
    if ( UseRandom::rndbool() ) swap(qq.first, qq.second);
    h[0] = flavGen->getHadron(qq.first, aqq.first)->produceParticle();
    h[1] = flavGen->getHadron(qq.second, aqq.second)->produceParticle();
  }
  if ( sqr(h[0]->mass() + h[1]->mass()) >= s ) return false;

  // Now set the momenta of the hadrons (distributed isotropically in
  // the cluster cm system).
  SimplePhaseSpace::CMS(s, h[0], h[1]);
  Utilities::transform(h, LorentzRotation(pc.boostVector()));

  // Add the hadrons as decay products of the partons in the
  // cluster. Returns false if it fails (although throwing an
  // exception may be more appropriate).
  if ( !newStep->addDecayProduct(cs.partons().begin(), cs.partons().end(),
				h.begin(), h.end()) ) {
    throw ClusterException(*this)
      << "Could not add decay products to the new step" << Exception::abortnow;
  }

  return true;

}

tcPDPair ClusterCollapser::getHadrons(const ColourSinglet & cs) const {
  tcPDPair ret;

  if ( cs.nPieces() == 3 ) {

    // This is a string with a junction. First find the tiplets.
    tcPDVector quarks = cs.getTripletData();
    if ( quarks.size() == 3 ) {
      // Three-quark junction. Select a quark to join into a meson
      // with a randomly chosed flavour.
      int i = UseRandom::irnd(3);
      tcPDPtr q = pickFlavour();
      if ( quarks[i]->iColour() == q->iColour() ) q = q->CC();
      ret.first = flavGen->getHadron(quarks[i], q);
      quarks[i] = q->CC();
      ret.second = flavGen->getBaryon(quarks[0], quarks[1], quarks[2]);
    }
    else if ( errorlevel )
      throw ClusterException(*this)
	<< "Too many diquarks in a junction string." << errorlevel;
  }

  else if ( cs.nPieces() != 1 ) {
    throw ClusterException(*this)
      << "Inconsistent number of string pieces in a cluster"
      << Exception::abortnow;
  }

  else if ( cs.piece(1).front()->data().iColour() == PDT::Colour8 ) {
    // This was a closed gluon loop. Create two random flavours.
    tcPDPtr q1 = pickFlavour();
    tcPDPtr q2 = pickFlavour();
    ret.first = flavGen->getHadron(q1, q2->CC());
    ret.second = flavGen->getHadron(q1->CC(), q2);
  }

  else {
    // This was a simple flat string. Pick a new flavour.
    tcPDPtr q = pickFlavour();
    if ( cs.piece(1).front()->data().iColour() == q->iColour() ) q = q->CC();
    ret.first = flavGen->getHadron(cs.piece(1).front()->dataPtr(), q);
    ret.second = flavGen->getHadron(cs.piece(1).back()->dataPtr(), q->CC());
  }

  if ( !ret.first || !ret.second ) throw ClusterException(*this)
    << "Could not generate hadrons from flavours " << cs.piece(1).front()->id()
    << " and " << cs.piece(1).back()->id() << "." << Exception::runerror;

  return ret;
}

tcPDPtr ClusterCollapser::pickFlavour() const {
  return getParticleData(2 + rndsign(1.0, 1.0, pStrange));
}

tcPDPtr ClusterCollapser::getHadron(const ColourSinglet & cs) const {
  if ( cs.nPieces() == 3 ) {

    // This is a string with a junction. First find the tiplets.
    tcPDVector quarks = cs.getTripletData();
    if ( quarks.size() == 3 ) {
      return flavGen->getBaryon(quarks[0], quarks[1], quarks[2]);
    }
    else if ( errorlevel )
      throw ClusterException(*this)
	<< "Too many diquarks in a junction string." << errorlevel;
  }

  else if ( cs.nPieces() != 1 ) {
    throw ClusterException(*this)
      << "Inconsistent number of string pieces in a cluster"
      << Exception::abortnow;
  }

  else if ( cs.piece(1).front()->data().iColour() == PDT::Colour8 ) {
    // This was a closed gluon loop. Create a random flavour.
    tcPDPtr q = pickFlavour();
    return flavGen->getHadron(q, q->CC());
  }

  // This was a simple flat string.
  return flavGen->getHadron(cs.piece(1).front()->dataPtr(),
			    cs.piece(1).back()->dataPtr());

}

  
void ClusterCollapser::persistentOutput(PersistentOStream & os) const {
  os << ounit(theEnergyCut, GeV) << theNTry2 << flavGen << oenum(errorlevel)
     << pStrange;
}

void ClusterCollapser::persistentInput(PersistentIStream & is, int) {
  is >> iunit(theEnergyCut, GeV) >> theNTry2 >> flavGen >> ienum(errorlevel)
     >> pStrange;
}

ClassDescription<ClusterCollapser> ClusterCollapser::initClusterCollapser;
// Definition of the static class description member.

void ClusterCollapser::Init() {

  static ClassDocumentation<ClusterCollapser> documentation
    ("The ThePEG::ClusterCollapser class can either be used as a "
     "preprocessor of a string fragmentation handler, or as a separate"
     "step handler to collapse small colour singlet systems of partons "
     "into one or two particles.");

  static Parameter<ClusterCollapser, Energy> interfaceEnergyCut
    ("EnergyCut",
     "If the invariant mass of a cluster, minus the constituent masses of its "
     "partons is below this cut (in GeV), it will be collapsed into one "
     "or two particles.",
     &ClusterCollapser::theEnergyCut, GeV, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Parameter<ClusterCollapser, int> interfaceNTry2
    ("NTry2",
     "The number of attempts to collapse a cluster into two particles, "
     "before it is collapsed into one particle.",
     &ClusterCollapser::theNTry2, 2, 0, 100, false, false, true);

  static Parameter<ClusterCollapser, double> interfacePStrange
    ("pStrange",
     "The relative probability to produce a s-sbar pair in a split as "
     "compared to a u-ubar or d-dbar pair.",
     &ClusterCollapser::pStrange, 1.0/3.0, 0.0, 2.0, false, false, true);

  static Switch<ClusterCollapser,Exception::Severity> interfaceLevel
    ("ErrorLevel",
     "What to do if a cluster could not be collapsed, or if momentum "
     "could not be conserved.",
     &ClusterCollapser::errorlevel, Exception::eventerror, true, false);
  static SwitchOption interfaceLevelNothing
    (interfaceLevel, "Nothing",
     "Do nothing, clusters may not collapse or momentum may not be conserved.",
     Exception::Severity(0));
  static SwitchOption interfaceLevelWarning
    (interfaceLevel, "Warning",
     "Report a warning, clusters may not collapse or momentum may not "
     "be conserved.", Exception::warning);
  static SwitchOption interfaceLevelEventError
    (interfaceLevel, "EventError",
     "Discard the whole event.", Exception::eventerror);
  static SwitchOption interfaceLevelRunError
    (interfaceLevel, "RunError",
     "End the run, printout the offending event.", Exception::runerror);
  static SwitchOption interfaceLevelAbort
    (interfaceLevel, "Abort",
     "Abort and dump core.", Exception::abortnow);

  static Reference<ClusterCollapser,FlavourGenerator> interfaceFlavGen
    ("FlavourGenerator",
     "The object used to combine quarks and diquarks into hadrons.",
     &ClusterCollapser::flavGen, true, false, true, false);
     
  interfaceEnergyCut.rank(10);
  interfacePStrange.rank(9);
  interfaceFlavGen.rank(8);
  interfaceNTry2.rank(7);
  interfaceLevel.rank(6);

}

