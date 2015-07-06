// -*- C++ -*-
//
// FastJetFinder.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2007 Leif Lonnblad
// Copyright (C) 2009-2012 Simon Platzer
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FastJetFinder class.
//

#include "FastJetFinder.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Cuts/Cuts.h"
#include "fastjet/ClusterSequence.hh"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;

FastJetFinder::FastJetFinder() 
  : theDCut(ZERO), theConeRadius(0.7), 
    theVariant(kt), theMode(inclusive),
    theRecombination(recoE) {}

FastJetFinder::~FastJetFinder() {}

void FastJetFinder::doinit() {
  fastjet::ClusterSequence::set_fastjet_banner_stream(0);
}

void FastJetFinder::doinitrun() {
  fastjet::ClusterSequence::set_fastjet_banner_stream(& CurrentGenerator::log());
  // Force banner to be printed right now; suppresses later output.
  fastjet::ClusterSequence::print_banner();
}

IBPtr FastJetFinder::clone() const {
  return new_ptr(*this);
}

IBPtr FastJetFinder::fullclone() const {
  return new_ptr(*this);
}

bool FastJetFinder::cluster(tcPDVector & ptype, vector<LorentzMomentum> & p,
			  tcCutsPtr, tcPDPtr, tcPDPtr) const {
  if ( ptype.size() <= minOutgoing() ){
    return false;
  }

  tcPDVector::iterator di = ptype.begin();
  vector<LorentzMomentum>::iterator pi = p.begin();
  size_t index = 0;

  vector<fastjet::PseudoJet> recombinables;
  tcPDVector ptypeBuffer;
  vector<LorentzMomentum> pBuffer;
  for ( ; di != ptype.end(); ++di, ++pi, index++ ) {
    if ( !unresolvedMatcher()->check(**di) ) {
      ptypeBuffer.push_back(*di);
      pBuffer.push_back(*pi);
      continue;
    }
    recombinables.push_back(fastjet::PseudoJet( (*pi).x()/GeV,(*pi).y()/GeV,(*pi).z()/GeV, (*pi).t()/GeV ));
    recombinables.back().set_user_index(index);
  }

  fastjet::Strategy strategy = fastjet::Best;

  fastjet::RecombinationScheme recomb_scheme;
  if ( theRecombination == recoE )
    recomb_scheme = fastjet::E_scheme;
  else if ( theRecombination == recoPt)
    recomb_scheme = fastjet::pt_scheme;
  else assert(false);

  fastjet::JetAlgorithm jet_algorithm = fastjet::kt_algorithm;
  if ( theVariant == CA ) {
    jet_algorithm = fastjet::cambridge_algorithm;
  } else if ( theVariant == antiKt ) {
    jet_algorithm = fastjet::antikt_algorithm;
  } else if ( theVariant > 3 ) {
    jet_algorithm = fastjet::ee_genkt_algorithm;
  }

  fastjet::JetDefinition jet_def;
  if ( theVariant < 4 ) {
    jet_def = fastjet::JetDefinition(jet_algorithm, theConeRadius, recomb_scheme, strategy);
  } else {
    int power = 1;
    if ( theVariant == sphericalCA ) {
      power = 0;
    }
    if ( theVariant == sphericalAntiKt ) {
      power = -1;
    }
    jet_def = fastjet::JetDefinition(jet_algorithm, theConeRadius, power, recomb_scheme, strategy);
  }

  fastjet::ClusterSequence clust_seq(recombinables, jet_def);

  double dcut = 0.0;
  
  if ( theVariant != antiKt &&
       theVariant != sphericalAntiKt ) {
    dcut = theDCut/GeV2;
  } else {
    dcut = theDCut != ZERO ? GeV2/theDCut : ZERO;
  }

  vector<fastjet::PseudoJet> recoJets;
  if ( theMode == inclusive ) 
    recoJets  = clust_seq.inclusive_jets();
  else if ( theMode == exclusive )
    recoJets  = clust_seq.exclusive_jets(dcut);


  if ( recoJets.size() + pBuffer.size() == p.size() ){
    return false;
  }
  else {
    tcPDVector ptypeNew;
    vector<LorentzMomentum> pNew;
    for (vector<fastjet::PseudoJet>::const_iterator iter = recoJets.begin();
	 iter != recoJets.end(); iter++){
      ptypeNew.push_back(ptype[iter->constituents().begin()->user_index()]);
      pNew.push_back(LorentzMomentum(iter->px()*GeV, iter->py()*GeV, iter->pz()*GeV, iter->E()*GeV));
    }

    di = ptypeBuffer.begin();
    pi = pBuffer.begin();
    for (; di != ptypeBuffer.end(); di++, pi++){
      ptypeNew.push_back(*di);
      pNew.push_back(*pi);
    }
    ptype = ptypeNew;
    p = pNew;
    return true;
  }
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void FastJetFinder::persistentOutput(PersistentOStream & os) const {
  os << ounit(theDCut,GeV2) << theConeRadius << theVariant << theMode
     << theRecombination;
}

void FastJetFinder::persistentInput(PersistentIStream & is, int) {
  is >> iunit(theDCut,GeV2) >> theConeRadius >> theVariant >> theMode
     >> theRecombination;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<FastJetFinder,JetFinder>
  describeFastJetFinder("ThePEG::FastJetFinder", "FastJetFinder.so");

void FastJetFinder::Init() {

  static ClassDocumentation<FastJetFinder> documentation
    ("FastJetFinder implements the class of longitudinally invariant kt "
     "jet clustering algorithms, as relevant for cuts on the real "
     "emission contribution to a NLO calculation. Recombination is "
     "exclusively performed using the pt scheme.");


  static Parameter<FastJetFinder,Energy2> interfaceDCut
    ("DCut",
     "The distance cut, when acting exclusively. "
     "The inverse is taken for the anti-kt algorithm, "
     "while for the Cambridge/Aachen variant dCut/GeV2 is used.",
     &FastJetFinder::theDCut, GeV2, 0.0*GeV2, 0.0*GeV2, 0*GeV2,
     false, false, Interface::lowerlim);


  static Parameter<FastJetFinder,double> interfaceConeRadius
    ("ConeRadius",
     "The cone radius R used in inclusive mode.",
     &FastJetFinder::theConeRadius, 0.7, 0.0, 10.0,
     false, false, Interface::limited);

  static Switch<FastJetFinder,int> interfaceVariant
    ("Variant",
     "The variant to use.",
     &FastJetFinder::theVariant, kt, false, false);
  static SwitchOption interfaceVariantKt
    (interfaceVariant,
     "Kt",
     "Kt algorithm.",
     kt);
  static SwitchOption interfaceVariantCA
    (interfaceVariant,
     "CA",
     "Cambridge/Aachen algorithm.",
     CA);
  static SwitchOption interfaceVariantAntiKt
    (interfaceVariant,
     "AntiKt",
     "Anti kt algorithm.",
     antiKt);
  static SwitchOption interfaceVariantSphericalKt
    (interfaceVariant,
     "SphericalKt",
     "Spherical kt algorithm.",
     sphericalKt);
  static SwitchOption interfaceVariantSphericalCA
    (interfaceVariant,
     "SphericalCA",
     "Spherical Cambridge/Aachen algorithm.",
     sphericalCA);
  static SwitchOption interfaceVariantSphericalAntiKt
    (interfaceVariant,
     "SphericalAntiKt",
     "Spherical anti kt algorithm.",
     sphericalAntiKt);

  static Switch<FastJetFinder,int> interfaceMode
    ("Mode",
     "The mode to use.",
     &FastJetFinder::theMode, inclusive, false, false);
  static SwitchOption interfaceModeInclusive
    (interfaceMode,
     "Inclusive",
     "Find inclusive jets.",
     inclusive);
  static SwitchOption interfaceModeExclusive
    (interfaceMode,
     "Exclusive",
     "Find exclusive jets.",
     exclusive);

  static Switch<FastJetFinder,int> interfaceRecombination
    ("RecombinationScheme",
     "The recombination scheme to use.",
     &FastJetFinder::theRecombination, recoE, false, false);
  static SwitchOption interfaceRecombinationPt
    (interfaceRecombination,
     "Pt",
     "Add transverse momenta",
     recoPt);
  static SwitchOption interfaceRecombinationE
    (interfaceRecombination,
     "E",
     "Add the four-momenta",
     recoE);

}

