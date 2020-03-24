// -*- C++ -*-
//
// SimpleBaryonRemnantDecayer.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SimpleBaryonRemnantDecayer class.
//

#include "SimpleBaryonRemnantDecayer.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/PDT/DecayMode.h"
#include "ThePEG/PDT/RemnantData.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include "ThePEG/Utilities/UtilityBase.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Utilities/HoldFlag.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/PDF/BeamParticleData.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/EventRecord/Step.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;

SimpleBaryonRemnantDecayer::~SimpleBaryonRemnantDecayer() {}

IBPtr SimpleBaryonRemnantDecayer::clone() const {
  return new_ptr(*this);
}

IBPtr SimpleBaryonRemnantDecayer::fullclone() const {
  return new_ptr(*this);
}

bool SimpleBaryonRemnantDecayer::accept(const DecayMode &) const {
  return true;
}

bool SimpleBaryonRemnantDecayer::
canHandle(tcPDPtr particle, tcPDPtr parton) const {
  return BaryonMatcher::Check(*particle) &&
    StandardQCDPartonMatcher::Check(*parton);
}

bool SimpleBaryonRemnantDecayer::
checkExtract(tcPPtr parent, tcPPtr parton, const LorentzMomentum & pnew) const {
  return pnew.e() > parent->nominalMass() + parton->data().constituentMass();
}

ParticleVector SimpleBaryonRemnantDecayer::
decay(const DecayMode & dm, const Particle & p, Step & step) const {
  ParticleVector children;
  tcRemPPtr remnant = dynamic_ptr_cast<tcRemPPtr>(&p);
  if ( !remnant ) return children;
  tRemPDPtr rpd = data(remnant);
  PVector ex = extracted(remnant);
  tcPPtr particle = parent(remnant);
  if ( !particle || ex.empty() || !rpd ) return children;

  // We can't handle multiple extractions (yet)
  if ( ex.size() != 1 ) return children;

  tPPtr parton = ex[0];
  tPVector subsys = getSubSystem(particle, parton);
  tPVector subpart;
  LorentzMomentum pitot = Utilities::sumMomentum(subsys) + remnant->momentum();

  Energy2 s = ZERO;
  Energy2 shat = ZERO;
  LorentzMomentum psub;
  Energy minmass = particle->nominalMass() +
    2.0*parton->nominalMass() + margin();

  unsigned int closeskip = 0;

  while ( closeskip < subsys.size() ) {

    subpart = tPVector(subsys.begin() + closeskip, subsys.end());

    while ( !subpart.empty() ) {

      psub = Utilities::sumMomentum(subpart);

      s = max((remnant->momentum() + psub).m2(), ZERO);
      shat = max(psub.m2(), ZERO);
      if ( subpart.size() == 1 ) shat = max(subpart[0]->momentum().mass2(), ZERO);

      if ( sqrt(s) > sqrt(shat) + minmass ) break;

      subpart.pop_back();

    }

    if ( !subpart.empty() ) break;

    ++closeskip;

  }

  if ( subpart.empty() ) {
    if ( respectDISKinematics() ) {
      // If we couldn't find any way of shuffling momentum and this may
      // have been a DIS event, try to include the scattered lepton as
      // well.
      {
	HoldFlag<int> nodis(respectDIS, 0);
	children = decay(dm, p, step);
      }
      if ( respectDISKinematics() == 1 )
	Throw<NoDISRespect>()
	  << "The decay of the remnant '" << p.PDGName()
	  << "' changed the kinematics of a scattered lepton in a DIS event.\n"
	  << *step.collision()->event()
	  << Exception::warning;
      return children;
    } else {
      Throw<DecayFailed>()
	<< "Could not decay remnant '" << p.PDGName()
	<< "' since not enough energy-momentum was available."
	<< Exception::eventerror;
    }
  }
  const BaryonContent & bi = getBaryonInfo(particle->dataPtr());

  // Check if we are extracting a valence quark
  int iq = parton->id();
  vector<int> vflav = bi.flav;
  vector<int>::iterator v = find(vflav.begin(), vflav.end(), bi.sign*iq);
  double pval = 0.0;
  if (  v != vflav.end() ) {
    vflav.erase(v);
    tcPDFPtr pdf;
    const BeamParticleData * beamp =
      dynamic_cast<const BeamParticleData *>(&*particle->dataPtr());
    if ( beamp ) pdf = beamp->pdf();
    if ( pdf && !specialValence() ) {
      Energy2 scale = abs(parton->scale());
      double x = shat/s;
      pval = pdf->xfvx(particle->dataPtr(), parton->dataPtr(), scale, x)/
	pdf->xfx(particle->dataPtr(), parton->dataPtr(), scale, x);
    } else {
      pval = 1;
    }
  }

  Energy mr = ZERO;

  while ( true ) {

    children.clear();

    if ( rndbool(pval) ) {

      // A simple valence remnant.
      int idqr = 1000*max(vflav[0], vflav[1]) +	100*min(vflav[0], vflav[1]) + 3;
      if ( vflav[0] != vflav[1] && rndbool(0.25) ) idqr -= 2;
      children.push_back(getParticleData(bi.sign*idqr)->produceParticle());
      mr = children[0]->mass();

    } else {
      // We haven't extracted a valence so we first divide up the baryon
      // in a quark and a diquark.
      pair<int,int> r = bi.flavsel.select(UseRandom::current());
      int iqr = r.first*bi.sign;
      int idqr = r.second*bi.sign;

      if ( iq == ParticleID::g ) {
	children.push_back(getParticleData(iqr)->produceParticle());
	children.push_back(getParticleData(idqr)->produceParticle());
      } else if ( iq*iqr > 0 ) {
	children.push_back(getParticleData(idqr)->produceParticle());
	children.push_back(flavourGenerator().getHadron
			   (getParticleData(-iq),
			    getParticleData(iqr))->produceParticle());
      } else {
	children.push_back(getParticleData(iqr)->produceParticle());
	children.push_back(flavourGenerator().getHadron
			   (getParticleData(-iq),
			    getParticleData(idqr))->produceParticle());
      }
      
      TransverseMomentum ptr = pTGenerator()->generate();
      Energy2 mt02 = children[0]->momentum().mass2() + ptr.pt2();
      Energy2 mt12 = children[1]->momentum().mass2() + ptr.pt2();
      double z = zGenerator().generate(children[1]->dataPtr(),
				       children[0]->dataPtr(), mt12);
      mr = sqrt(mt02/(1.0 - z) + mt12/z);

      if ( sqrt(s) <= sqrt(shat) + mr ) continue;

      children[0]->set3Momentum(static_cast<const LorentzMomentum &>
				(lightCone((1.0 - z)*mr,
					   mt02/((1.0 - z)*mr), ptr)));
      children[1]->set3Momentum(static_cast<const LorentzMomentum &>
				(lightCone(z*mr, mt12/(z*mr), -ptr)));
    }

    break;

  }

  // Make copies of all final particles in the hard subsystem which
  // will take recoil.
  for ( unsigned int i = 0, N = subsys.size(); i < N; ++i ) {
    if ( subsys[i]->birthStep() != &step )
      subsys[i] = step.copyParticle(subsys[i]);
    if ( i >= closeskip && i - closeskip < subpart.size() )
      subpart[i - closeskip] = subsys[i];
  }

  // Boost part of the hard subsystem to give energy to the new
  // remnants.
  LorentzMomentum pr = remnant->momentum();
  LorentzRotation R = Utilities::getBoostToCM(make_pair(psub, pr));
  Energy pz = SimplePhaseSpace::getMagnitude(s, sqrt(shat), mr);
  if ( subpart.size() > 1 ) {
    LorentzRotation Rs(-(R*psub).boostVector());
    Rs.boost(0.0, 0.0, pz/sqrt(sqr(pz) + shat));
    Rs = Rs*R;
    R.invert();
    Rs = R*Rs;
    Utilities::transform(subpart, Rs);
  } else {
    subpart[0]->set5Momentum(
      Lorentz5Momentum(ZERO, ZERO, pz,
		       sqrt(sqr(pz) + shat), sqrt(shat)));
    R.invert();
    Utilities::transform(subpart, R);
  }

  LorentzRotation Rr(0.0, 0.0, -pz/sqrt(sqr(pz) + sqr(mr)));
  Rr = R*Rr;
  Utilities::transform(children, Rr);

  // Give the remnants and subsystem a transverse momentum by Lorentz
  // rotations.
  LorentzMomentum pr0 = ( pr = Utilities::sumMomentum(children) );
  LorentzMomentum psub0 = ( psub = Utilities::sumMomentum(subsys) );
  LorentzMomentum ksub = pr + psub - particle->momentum();
  R = Utilities::boostToCM(make_pair(&psub, &pr));
  TransverseMomentum kt;
  do {
    kt = pTGenerator()->generate();
  } while ( kt.pt() >= psub.z() );
  LorentzRotation Rtot = R;
  Rtot.rotateY(asin(kt.pt()/psub.z()));
  Rtot.rotateZ(kt.phi());
  Rtot = R.inverse()*Rtot;
  psub = Rtot*psub0;
  pr = Rtot*pr0;

  Utilities::transform(children, Rtot);
  if ( subsys.size() > 1 )
    Utilities::transform(subsys,
			 Utilities::getTransformToMomentum(psub0, psub, ksub));
  else
    subsys[0]->setMomentum(psub);
  // Make small z-boosts to correct 
  Utilities::transform(subsys, getZBoost(Utilities::sumMomentum(subsys), psub));

  LorentzMomentum pftot =
    Utilities::sumMomentum(subsys) + Utilities::sumMomentum(children);

  R = getZBoost(pftot, pitot);
  Utilities::transform(subsys, R);
  Utilities::transform(children, R);

  return children;
}

const SimpleBaryonRemnantDecayer::BaryonContent &
SimpleBaryonRemnantDecayer::getBaryonInfo(tcPDPtr baryon) const {
  map<tcPDPtr,BaryonContent>::iterator it = baryonmap.find(baryon);
  if ( it != baryonmap.end() ) return it->second;
  BaryonContent & bi = baryonmap[baryon];
  int pid = baryon->id();
  bi.sign = pid < 0? -1: 1;
  bi.flav = vector<int>(3);
  bi.flav[0] = (pid = abs(pid)/10)%10;
  bi.flav[1] = (pid /= 10)%10;
  bi.flav[2] = (pid /= 10)%10;
  bi.flavsel = VSelector< pair<int,int> >();
  for ( int iq1 = 0; iq1 < 3; ++iq1 ) {
    int iq2 = (iq1 + 1)%3;
    int iq3 = (iq2 + 1)%3;
    int idq = 1000*max(bi.flav[iq2], bi.flav[iq3]) +
      100*min(bi.flav[iq2], bi.flav[iq3]) + 3;
    bi.flavsel.insert(3.0, make_pair(bi.flav[iq1], idq));
    if ( bi.flav[iq2] == bi.flav[iq3] ) continue;
    bi.flavsel.insert(1.0, make_pair(bi.flav[iq1], idq - 2));
  }
  return bi;
}

bool SimpleBaryonRemnantDecayer::preInitialize() const {
  return RemnantDecayer::preInitialize() ||
    !theZGenerator || !theFlavourGenerator;
}

void SimpleBaryonRemnantDecayer::doinit() {
  RemnantDecayer::doinit();
  if ( !theZGenerator ) {
    theZGenerator = dynamic_ptr_cast<ZGPtr>
      (generator()->preinitCreate("ThePEG::SimpleZGenerator",
				  fullName() + "/ZGen",
				  "SimpleZGenerator.so"));
  }
  if ( !theFlavourGenerator ) {
    theFlavourGenerator = dynamic_ptr_cast<FlGPtr>
      (generator()->preinitCreate("ThePEG::SimpleFlavour",
				  fullName() + "/FlavGen",
				  "SimpleFlavour.so"));
  }
}


void SimpleBaryonRemnantDecayer::
persistentOutput(PersistentOStream & os) const {
  os << theZGenerator << theFlavourGenerator << ounit(theMargin,GeV)
     << useSpecialValence;
}

void SimpleBaryonRemnantDecayer::
persistentInput(PersistentIStream & is, int) {
  is >> theZGenerator >> theFlavourGenerator >> iunit(theMargin,GeV)
     >> useSpecialValence;
}

ClassDescription<SimpleBaryonRemnantDecayer> SimpleBaryonRemnantDecayer::initSimpleBaryonRemnantDecayer;
// Definition of the static class description member.

void SimpleBaryonRemnantDecayer::Init() {

  static ClassDocumentation<SimpleBaryonRemnantDecayer> documentation
    ("The SimpleBaryonRemnantDecayer class inherits from the RemnantDecayer "
     "class and is able to decay RemnantParticles produced by the "
     "SoftRemnantHandler class for the cases when a single parton has been "
     "extracted from a baryon.");

  static Reference<SimpleBaryonRemnantDecayer,ZGenerator> interfaceZGenerator
    ("ZGenerator",
     "The object responsible for generating momentum fractions in case "
     "of more than one remnant. If not set and the controlling EventGenerator "
     "has a default ZGenerator object, this will be used. Otherwise a "
     "SimpleZGenerator object created with default settings in the "
     "initialization will be used instead.",
     &SimpleBaryonRemnantDecayer::theZGenerator, true, false, true, true, true);

  static Reference<SimpleBaryonRemnantDecayer,FlavourGenerator>
    interfaceFlavourGenerator
    ("FlavourGenerator",
     "The object responsible for handling the flavour contents of a baryon. "
     "If not set and the controlling EventGenerator "
     "has a default FlavourGenerator object, this will be used. Otherwise a "
     "SimpleFlavour object created with default settings in the "
     "initialization will be used instead.",
     &SimpleBaryonRemnantDecayer::theFlavourGenerator,
     true, false, true, true, true);

  static Parameter<SimpleBaryonRemnantDecayer,Energy> interfaceMargin
    ("EnergyMargin",
     "The energy margin (in GeV) to be added to the sum of the parent and "
     "parton masses to determine if it is possible to construct the remnants "
     "with the given (upper limit of the) virtuality of the extracted parton.",
     &SimpleBaryonRemnantDecayer::theMargin, GeV, 1.0*GeV, ZERO, 10.0*GeV,
     false, false, true);

  static Switch<SimpleBaryonRemnantDecayer,bool> interfaceSpecialValence
    ("SpecialValence",
     "If true, an extracted valence quark will always give a di-quark remnant.",
     &SimpleBaryonRemnantDecayer::useSpecialValence, false, true, false);
  static SwitchOption interfaceSpecialValenceYes
    (interfaceSpecialValence,
     "Yes",
     "An extracted valence quark will always give a di-quark remnant.",
     true);
  static SwitchOption interfaceSpecialValenceNo
    (interfaceSpecialValence,
     "No",
     "An extracted valence flavour may be considered to be a sea-quark.",
     false);

  interfaceZGenerator.rank(8);
  interfaceFlavourGenerator.rank(7);

}

