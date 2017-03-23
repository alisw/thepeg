// -*- C++ -*-
//
// PartonExtractor.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the PartonExtractor class.
//

#include "PartonExtractor.h"
#include "ThePEG/Handlers/XComb.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/PDF/NoPDF.h"
#include "ThePEG/PDF/RemnantHandler.h"
#include "ThePEG/PDF/BeamParticleData.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/Step.h"
#include "ThePEG/EventRecord/SubProcess.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Utilities/UtilityBase.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace ThePEG;

PartonExtractor::PartonExtractor()
  : theMaxTries(100), flatSHatY(false) {}

PartonExtractor::~PartonExtractor() {}

IBPtr PartonExtractor::clone() const {
  return new_ptr(*this);
}

IBPtr PartonExtractor::fullclone() const {
  return new_ptr(*this);
}

PartonPairVec PartonExtractor::
getPartons(Energy maxEnergy, const cPDPair & incoming,
	   const Cuts & kc) const {
  PartonPairVec result;
  PartonVector first;
  PDFCuts cuts1(kc, true, maxEnergy);
  PBPtr p1 =
    new_ptr(PartonBin(PDPtr(), PBPtr(), incoming.first, PDFPtr(), cuts1));
  addPartons(p1, cuts1,  theFirstPDF, first);
  PartonVector second;
  PDFCuts cuts2(kc, false, maxEnergy);
  PBPtr p2 =
    new_ptr(PartonBin(PDPtr(), PBPtr(), incoming.second, PDFPtr(), cuts2));
  addPartons(p2, cuts2, theSecondPDF, second);
  for ( PartonVector::iterator it1 = first.begin();
	it1 != first.end(); ++it1 )
    for ( PartonVector::iterator it2 = second.begin();
	it2 != second.end(); ++it2 )
      result.push_back(PBPair(*it1, *it2));

  // We add the original parton bins as well to avoid them being
  // deleted.
  result.push_back(PBPair(p1, p2));
  return result;
}

void PartonExtractor::
addPartons(tPBPtr incoming, const PDFCuts & cuts, tcPDFPtr pdf,
	   PartonVector & pbins) const {
  if(!pdf) pdf = getPDF(incoming->parton());
  if ( dynamic_ptr_cast<Ptr<NoPDF>::tcp>(pdf) ||
       incoming->parton() == incoming->particle() ) {
    pbins.push_back(incoming);
    return;
  }
  cPDVector partons = pdf->partons(incoming->parton());
  for ( int i = 0, N = partons.size(); i < N; ++i ) {
    PBPtr pb =
      new_ptr(PartonBin(incoming->parton(), incoming, partons[i], pdf, cuts));
    incoming->addOutgoing(pb);
    addPartons(pb, cuts, PDFPtr(), pbins);
  }

}

tcPDFPtr PartonExtractor::getPDF(tcPDPtr particle) const {
  for ( vector<PDFPtr>::const_iterator it = theSpecialDensities.begin();
	it != theSpecialDensities.end(); ++it )
    if ( (**it).canHandle(particle) ) return *it;
  Ptr<BeamParticleData>::tcp p = 
    dynamic_ptr_cast<Ptr<BeamParticleData>::tcp>(particle);
  if ( !p || !p->pdf() ) return noPDF();
  return p->pdf();
}

void PartonExtractor::select(tXCombPtr newXComb) {
  theLastXComb = newXComb;
}

tPBIPtr PartonExtractor::partonBinInstance(tcPPtr p) const {
  PartonBinInstanceMap::const_iterator it = partonBinInstances().find(p);
  return it == partonBinInstances().end()? PBIPtr(): it->second;
}

void PartonExtractor::
colourConnect(tPPtr particle, tPPtr parton, const tPVector & remnants) const {

  // Sorry cannot handle coloured resolved particles.
  if ( particle->coloured() ) throw RemColException(*this);

  // First connect the loose colour line from the extacted parton.
  if ( parton->hasColour() )
    findConnect(parton->colourLine(), parton, true,
		remnants.rbegin(), remnants.rend());

  // First connect the loose anti-colour line from the extacted parton.
  if ( parton->hasAntiColour() )
    findConnect(parton->antiColourLine(), parton, false,
		remnants.begin(), remnants.end());

  // Go through the rest of the remnants and create new colour lines
  // if needed. Go through it forwards and backwards to catch possible
  // inconsistencies.
  for ( tPVector::const_iterator it = remnants.begin();
	it != remnants.end(); ++it ) {
    if ( (**it).hasAntiColour() && !(**it).antiColourLine() )
      findConnect(ColourLine::create(*it, true), *it, false,
		it + 1, remnants.end());
  }
  for ( tPVector::const_reverse_iterator it = remnants.rbegin();
	it != remnants.rend(); ++it ) {
    if ( (**it).hasColour() && !(**it).colourLine() )
      findConnect(ColourLine::create(*it), *it, true, it + 1, remnants.rend());
  } 
}

Energy2 PartonExtractor::newScale() {
  return lastScale();
}

pair<int,int> PartonExtractor::nDims(const PBPair & pbins) {
  // if photon from a lepton or proton generate scale
  bool genscale[2]={false,false};
  for(unsigned int ix=0;ix<2;++ix) {
    PBPtr bin = ix==0 ? pbins.first : pbins.second;
    if (!bin || !bin->particle() || !bin->parton()) continue;
    if(bin->pdf()->partons(bin->particle()).size()==1 &&
       bin->particle()->id()!=bin->parton()->id())
      genscale[ix]=true;
  }
  return make_pair(pbins.first ->nDim(genscale[0]),
		   pbins.second->nDim(genscale[1]));
}

void PartonExtractor::prepare(const PBIPair & pbins) {
  partonBinInstances().clear();
  pbins.first->prepare();
  pbins.second->prepare();
}

void PartonExtractor::updatePartonBinInstances(const PBIPair & pbins) {
  partonBinInstances().clear();
  tPBIPtr current = pbins.first;
  while ( current->incoming() ) {
    partonBinInstances()[current->parton()] = current;
    current = current->incoming();
  }
  current = pbins.second;
  while ( current->incoming() ) {
    partonBinInstances()[current->parton()] = current;
    current = current->incoming();
  }
}

bool PartonExtractor::
generateL(const PBIPair & pbins, const double * r1, const double * r2) {
  Direction<0> dir(true);
  generateL(*pbins.first, r1);
  dir.reverse();
  generateL(*pbins.second, r2);
  if ( !flatSHatY || pbins.first->hasPoleIn1() || pbins.second->hasPoleIn1() )
    return true;

  Energy2 shmax = lastCuts().sHatMax();
  Energy2 shmin = lastCuts().sHatMin();
  Energy2 sh = shmin*pow(shmax/shmin, *r1);
  double ymax = lastCuts().yHatMax();
  double ymin = lastCuts().yHatMin();
  double km = log(shmax/shmin);
  ymax = min(ymax, log(lastCuts().x1Max()*sqrt(lastS()/sh)));
  ymin = max(ymin, -log(lastCuts().x2Max()*sqrt(lastS()/sh)));

  double y = ymin + (*r2)*(ymax - ymin);
  double l1 = 0.5*log(lastS()/sh) - y;
  double l2 = 0.5*log(lastS()/sh) + y;

  pbins.first->li(l1 - pbins.first->l() + pbins.first->li());
  pbins.first->l(l1);
  pbins.first->jacobian(km*(ymax - ymin));
  pbins.second->li(l2 - pbins.second->l() + pbins.second->li());
  pbins.second->l(l2);
  pbins.second->jacobian(1.0);
  return ( pbins.first->li() >= 0.0 && pbins.second->li() >= 0.0 );
}

Energy2 PartonExtractor::
generateSHat(Energy2, const PBIPair & pbins,
	     const double * r1, const double * r2,
	     bool haveMEPartons) {
  Direction<0> dir(true);
  if(pbins.first->bin()->pdfDim()<=1) pbins.first->scale(-lastScale());
  if ( !generate(*pbins.first, r1, lastSHat(),
		 pbins.first->getFirst()->parton()->momentum(),
		 haveMEPartons) )
    return -1.0*GeV2;
  dir.reverse();
  if(pbins.second->bin()->pdfDim()<=1) pbins.second->scale(-lastScale());
  if ( !generate(*pbins.second, r2, lastSHat(),
		 pbins.second->getFirst()->parton()->momentum(),
		 haveMEPartons) )
    return -1.0*GeV2;
  
  return (pbins.first->parton()->momentum() +
	  pbins.second->parton()->momentum()).m2();
}

void PartonExtractor::
generateL(PartonBinInstance & pb, const double * r) {
  if ( !pb.incoming() ) return;

  pb.parton(pb.partonData()->produceParticle(Lorentz5Momentum()));
  generateL(*pb.incoming(), r + pb.bin()->pdfDim() + pb.bin()->remDim());
  pb.particle(pb.incoming()->parton());

  if ( pb.li() >= 0 ) return;

  double jac = 1.0;
  if ( pb.bin()->pdfDim() )
    pb.li(pb.pdf()->flattenL(pb.particleData(), pb.partonData(),
			     pb.bin()->cuts(), *r++, jac));
  pb.scale(-1.0*GeV2);
  if ( pb.bin()->pdfDim() > 1 )
    pb.scale(pb.pdf()->flattenScale(pb.particleData(), pb.partonData(),
				    pb.bin()->cuts(), pb.li(), *r++, jac)
	     *pb.bin()->cuts().scaleMaxL(pb.li()));
  pb.jacobian(jac);
  pb.l(pb.incoming()->l() + pb.li());
}

bool PartonExtractor::
generate(PartonBinInstance & pb, const double * r,
	 Energy2 shat, const Lorentz5Momentum & first,
	 bool haveMEPartons) {
  if ( !pb.incoming() ) return true;
  if ( !generate(*pb.incoming(), r + pb.bin()->pdfDim() + pb.bin()->remDim(),
		 shat/pb.xi(), first) )
    return false;
  pb.remnantWeight(1.0);
  pb.parton()->setMomentum
    (pb.remnantHandler()->generate(pb, r + pb.bin()->pdfDim(), pb.scale(), shat,
				   pb.particle()->momentum(),haveMEPartons));
  if ( pb.remnantWeight() <= 0.0 ) return false;
  partonBinInstances()[pb.parton()] = &pb;
  return true;
}

void PartonExtractor::
constructRemnants(const PBIPair & pbins, tSubProPtr sub, tStepPtr step) const {
  partonBinInstances().clear();
  LorentzMomentum k1 = pbins.first->parton()->momentum();
  LorentzMomentum k2 = pbins.second->parton()->momentum();
  LorentzMomentum Ph = k1 + k2;
  LorentzMomentum Phold = Ph;
  LorentzRotation Rh = Utilities::getBoostToCM(make_pair(k1, k2));

  bool pickside = rndbool();
  if ( pickside && pbins.first->incoming() ) {
    Direction<0> dir(true);
    constructRemnants(*pbins.first, Ph, k2);
    construct(*pbins.first, step, false);
  }
  if ( pbins.second->incoming() ) {
    Direction<0> dir(false);
    constructRemnants(*pbins.second, Ph, pbins.first->parton()->momentum());
    construct(*pbins.second, step, false);
  }
  if ( (!pickside) && pbins.first->incoming() ) {
    Direction<0> dir(true);
    constructRemnants(*pbins.first, Ph, pbins.second->parton()->momentum());
    construct(*pbins.first, step, false);
  }
  //  LorentzRotation rot = Utilities::transformToMomentum(Phold, Ph);
  k1 = pbins.first->parton()->momentum();
  k2 = pbins.second->parton()->momentum();
  LorentzRotation rot = Utilities::getBoostFromCM(make_pair(k1, k2))*Rh;
  Utilities::transform(sub->outgoing(), rot);
  Utilities::transform(sub->intermediates(), rot);
  Ph = k1 + k2;
  if ( abs(Ph.m2() - Phold.m2())/Phold.m2() > 0.000001 )
    cerr << Ph.m2()/GeV2 << " was (" << Phold.m2()/GeV2 << ")" << endl;
}

void PartonExtractor::
constructRemnants(PartonBinInstance & pb, LorentzMomentum & Ph,
		  const LorentzMomentum & k) const {
  LorentzMomentum P = pb.particle()->momentum();
  DVector r = UseRandom::rndvec(pb.bin()->remDim());
  if ( r.empty() ) r.push_back(0.0);
  pb.parton()->setMomentum(pb.remnantHandler()->
			   generate(pb, &r[0], pb.scale(), Ph.m2(), P));
  if ( pb.remnantWeight() <= 0.0 ) throw Veto();
  pb.remnantHandler()->boostRemnants(pb);
  LorentzMomentum Pr = Utilities::sumMomentum(pb.remnants());
  transformRemnants(Ph, Pr, k, pb.particle()->momentum());
  pb.parton()->setMomentum(pb.particle()->momentum() - Pr);
  try {
    Utilities::setMomentum(pb.remnants().begin(),
			   pb.remnants().end(),
			   static_cast<const LorentzMomentum &>(Pr));
  }
  catch ( ThePEG::Exception & e) {
    throw e;
  }
  catch ( ThePEG::Veto ) {
    throw;
  }
  catch ( std::exception & e ) {
    throw Exception() << "Caught non-ThePEG exception " << e.what() << "in "
		      << "PartonExtractor::constructRemnants"
		      << Exception::eventerror;
  }
  partonBinInstances()[pb.parton()] = &pb;
  if ( !pb.incoming()->incoming() ) return;

  // We get here if we need to construct remnants recursively.
  LorentzMomentum Phnew = Ph + Pr;
  constructRemnants(*pb.incoming(), Phnew, k);
  LorentzRotation rot = Utilities::transformToMomentum(Ph + Pr, Phnew);
  Utilities::transform(pb.remnants(), rot);
  Ph.transform(rot);
}

LorentzRotation PartonExtractor::
boostRemnants(PBIPair & bins, LorentzMomentum k1, LorentzMomentum k2,
	 bool side1, bool side2) const {
  if ( !side1 && !side2 ) return LorentzRotation();

  LorentzMomentum P1 =
    bins.first? LorentzMomentum(bins.first->parton()->momentum()): k1;
  LorentzMomentum Pr1;
  if ( side1 ) {
    P1 = bins.first->particle()->momentum();
    Pr1 = Utilities::sumMomentum(bins.first->remnants());
  }
  LorentzMomentum P2 =
    bins.second? LorentzMomentum(bins.second->parton()->momentum()): k2;
  LorentzMomentum Pr2;
  if ( side2 ) {
    P2 = bins.second->particle()->momentum();
    Pr2 = Utilities::sumMomentum(bins.second->remnants());
  }

  LorentzRotation Rh = Utilities::getBoostToCM(make_pair(k1, k2));
  LorentzMomentum Ph = k1 + k2;
  //  LorentzMomentum Phold = Ph;

  bool otherside = rndbool();
  if ( otherside && side2 ){
    Direction<0> dir(false);
    transformRemnants(Ph, Pr2, k1, P2);
    k2 = P2 - Pr2;
  }
  if ( side1 ){
    Direction<0> dir(true);
    transformRemnants(Ph, Pr1, k2, P1);
    k1 = P1 - Pr1;
  }
  if ( side2 && !otherside ) {
    Direction<0> dir(false);
    transformRemnants(Ph, Pr2, k1, P2);
    k2 = P2 - Pr2;
  }

  if ( bins.first ) {
    if ( bins.first->remnants().size() == 1 )
      bins.first->remnants()[0]->setMomentum(Pr1);
    else
      Utilities::setMomentum(bins.first->remnants().begin(),
			     bins.first->remnants().end(),
			     static_cast<const LorentzMomentum &>(Pr1));
    bins.first->parton()->setMomentum(k1);
  }

  if ( bins.second ) {
    if ( bins.second->remnants().size() == 1 )
      bins.second->remnants()[0]->setMomentum(Pr2);
    else
      Utilities::setMomentum(bins.second->remnants().begin(),
			     bins.second->remnants().end(),
			     static_cast<const LorentzMomentum &>(Pr2));
    bins.second->parton()->setMomentum(k2);
  }

  Rh.transform(Utilities::getBoostFromCM(make_pair(k1, k2)));

  // LorentzMomentum phh = Rh*Phold;

  return Rh;

  //  return Utilities::transformToMomentum(Phold, Ph);

}

void PartonExtractor::
transformRemnants(LorentzMomentum & Ph, LorentzMomentum & Pr,
		  const LorentzMomentum & k, const LorentzMomentum & P) const {
  // don't do this for very soft remnants, as
  // we may run into numerical troubles; threshold
  // needs to become a parameter at some point
  if ( Pr.vect().mag2()/k.vect().mag2() < 1e-10 &&
       sqr(Pr.e()/k.e()) < 1e-10 )
    return;
  TransverseMomentum pt = Pr;
  try {
    if ( Direction<0>::pos() )
      SimplePhaseSpace::CMS(Pr, Ph, (P + k).m2(), 1.0, 0.0);
    else
      SimplePhaseSpace::CMS(Ph, Pr, (k + P).m2(), 1.0, 0.0);
    LorentzRotation rpt;
    if ( sqr(Pr.z()) > ZERO ) rpt.rotateY(asin(pt.pt()/Pr.z()));
    rpt.rotateZ(pt.phi());
    rpt = Direction<0>::pos()?
      Utilities::getBoostFromCM(make_pair(P, k))*rpt:
      Utilities::getBoostFromCM(make_pair(k, P))*rpt;
    Ph.transform(rpt);
    Pr.transform(rpt);
  } catch ( ImpossibleKinematics ) {}
}


double PartonExtractor::fullFn(const PBIPair & pbins, Energy2 scale,
			       pair<bool,bool> noLastPDF) {
  if(pbins.first->bin()->pdfDim()<=1) pbins.first->scale(scale);
  if(pbins.second->bin()->pdfDim()<=1) pbins.second->scale(scale);
  return fullFn(*pbins.first,noLastPDF.first)*fullFn(*pbins.second,noLastPDF.second);
}

double PartonExtractor::fullFn(const PartonBinInstance & pb,
			       bool noLastPDF) {
  if ( !pb.incoming() ) return 1.0;
  if (noLastPDF)
    return 
      fullFn(*pb.incoming(),false) * pb.jacobian() * 
      pb.remnantWeight() * exp(-pb.li());
  return fullFn(*pb.incoming(),false) * pb.jacobian() * pb.remnantWeight() *
    pb.pdf()->xfl(pb.particleData(), pb.partonData(), pb.scale(),
		  pb.li(), pb.incoming()->scale());
}

void PartonExtractor::
construct(const PBIPair & pbins, tStepPtr step) const {
  // if a long chain we need to break some mother/child relationships
  if(pbins.first->incoming()) {
    if(pbins.first->incoming()->incoming()) {
      if(!pbins.first->parton()->parents().empty()) {
	tParticleVector parents=pbins.first->parton()->parents();
	tPPtr parton = pbins.first->parton();
	for(unsigned int ix=0;ix<parents.size();++ix) parents[ix]->abandonChild(parton);
      }
    }
  }
  if(pbins.second->incoming()) {
    if(pbins.second->incoming()->incoming()) {
      if(!pbins.second->parton()->parents().empty()) {
	tParticleVector parents=pbins.second->parton()->parents();
	tPPtr parton = pbins.second->parton();
	for(unsigned int ix=0;ix<parents.size();++ix) parents[ix]->abandonChild(parton);
      }
    }
  }
  Direction<0> dir(true);
  construct(*pbins.first, step);
  dir.reverse();
  construct(*pbins.second, step);
}

void PartonExtractor::
construct(PartonBinInstance & pb, tStepPtr step, bool boost) const {
  if ( !pb.incoming() ) return;
  if ( boost ) pb.remnantHandler()->boostRemnants(pb);
  if ( pb.incoming()->incoming() ) {
    step->insertIntermediate(pb.particle(),pb.incoming()->particle(),pb.parton());
  }
  tPVector rem(pb.remnants().begin(), pb.remnants().end());
  if ( !step->addDecayProduct(pb.particle(), rem.begin(), rem.end(), false) )
    {}
  colourConnect(pb.particle(), pb.parton(), rem);
  construct(*pb.incoming(), step);
}

PBIPair PartonExtractor::newRemnants(tPPair oldp, tPPair newp, tStepPtr step) {
  PBIPair pb;
  Direction<0> dir(true);
  pb.first = newRemnants(partonBinInstance(oldp.first),
			 newp.first, newp.second->momentum());
  dir.reverse();
  pb.second = newRemnants(partonBinInstance(oldp.second),
			  newp.second, newp.first->momentum());
  addNewRemnants(partonBinInstance(oldp.first), pb.first, step);
  addNewRemnants(partonBinInstance(oldp.second), pb.second, step);
  return pb;
}

PBIPtr PartonExtractor::
newRemnants(tPBIPtr oldpb, tPPtr newp, const LorentzMomentum & k) {
  if ( ! oldpb || !oldpb->incoming() ) return oldpb;
  Energy2 shat = (k + newp->momentum()).m2();

  // Loop over all possible PartonBin sisters to find the one
  // corresponding to the newly extracted parton.
  const PartonBin::PBVector & sisters = oldpb->incoming()->bin()->outgoing();
  for ( int i = 0, N = sisters.size(); i < N; ++i )
    if ( sisters[i]->parton() == newp->dataPtr() ) {
      // Setup necessary info in new PartonBinInstance object.
      PBIPtr newpb = new_ptr(PartonBinInstance(sisters[i], oldpb->incoming()));
      newpb->particle(oldpb->particle());
      newpb->parton(newp);
      newpb->li(log(oldpb->particle()->momentum().dirPlus()/
		    newp->momentum().dirPlus()));
      newpb->l(oldpb->l() - oldpb->li() + newpb->li());
      Energy2 sc = -newp->scale();
      newpb->scale(newp->scale());
      if ( oldpb->incoming()->incoming() )
	sc = -newpb->particle()->momentum().m2();
      // Now we can construct the new remnants.
      newpb->remnantWeight(1.0);
      if ( !newpb->remnantHandler()->
	   recreateRemnants(*newpb, oldpb->parton(), newp, newpb->li(),
			    sc, shat, newpb->particle()->momentum()) )
	throw Veto();
      if ( newpb->remnantWeight() <= 0.0 ) throw Veto();
      return newpb;
    }
  throw Veto();
}

void PartonExtractor::
addNewRemnants(tPBIPtr oldpb, tPBIPtr newpb, tStepPtr step) {
  if ( oldpb == newpb ) return;
  if ( oldpb->parton() != newpb->parton() ) {
    step->removeDecayProduct(newpb->particle(), oldpb->parton());
    if ( !step->addDecayProduct(newpb->particle(), newpb->parton()) )
      throw Veto();
  }
  tPVector rem(newpb->remnants().begin(), newpb->remnants().end());
  colourConnect(newpb->particle(), newpb->parton(), rem);
  partonBinInstances()[newpb->parton()] = newpb;
  if ( !step->addDecayProduct(oldpb->remnants().begin(),
			      oldpb->remnants().end(),
 			      rem.begin(), rem.end()) )
    throw Veto();
}

void PartonExtractor::persistentOutput(PersistentOStream & os) const {
  os << theLastXComb << theSpecialDensities << theNoPDF << theMaxTries
     << flatSHatY << theFirstPDF << theSecondPDF;
}

void PartonExtractor::persistentInput(PersistentIStream & is, int) {
  is >> theLastXComb >> theSpecialDensities >> theNoPDF >> theMaxTries
     >> flatSHatY >> theFirstPDF >> theSecondPDF;
}

ClassDescription<PartonExtractor> PartonExtractor::initPartonExtractor;

void PartonExtractor::Init() {

  static ClassDocumentation<PartonExtractor> documentation
    ("There is no documentation for the ThePEG::PartonExtractor class");

  static RefVector<PartonExtractor,PDFBase> interfaceSpecialDensities
    ("SpecialDensities",
     "A list of parton density objects to be used for incoming particles "
     "overriding possible densities given for particles of the "
     "BeamParticleData class.",
     &PartonExtractor::theSpecialDensities, 0, false, false, true, false);

  static Reference<PartonExtractor,PDFBase> interfaceNoPDF
    ("NoPDF",
     "A fixed reference to a NoPDF object to be used for particles without "
     "substructure.",
     &PartonExtractor::theNoPDF, true, true, true, false);

  static Parameter<PartonExtractor,int> interfaceMaxTries
    ("MaxTries",
     "The maximum number of attempts allowed when trying to generate "
     "remnants.",
     &PartonExtractor::theMaxTries, 100, 1, 1000, false, false, true);

  static Switch<PartonExtractor,bool> interfaceFlatSHatY
    ("FlatSHatY",
     "The possibility to override the l-generation in the PDFs and generate "
     "a flat distribution in \\f$\\log(\\hat{s})\\f$ and \\f$y\\f$. This only "
     "applies if the parton densities do not have poles in \\f$x=1\\f$.",
     &PartonExtractor::flatSHatY, false, false, false);

  static SwitchOption interfaceFlatSHatY0
    (interfaceFlatSHatY,
     "Off", "Use the l-generation defined by the PDFs", false);

  static SwitchOption interfaceFlatSHatY1
    (interfaceFlatSHatY,
     "On", "Generate flat rapidity and \\f$\\log(\\hat{s})\\f$", true);

  static Reference<PartonExtractor,PDFBase> interfaceFirstPDF
    ("FirstPDF",
     "PDF to override the default PDF for the first beam particle",
     &PartonExtractor::theFirstPDF, false, false, true, true, false);

  static Reference<PartonExtractor,PDFBase> interfaceSecondPDF
    ("SecondPDF",
     "PDF to override the default PDF for the second beam particle",
     &PartonExtractor::theSecondPDF, false, false, true, true, false);

}

RemColException::RemColException(const PartonExtractor & pe) {
  theMessage << "Parton extractor '" << pe.name() << "' failed to connect "
	     << "the colours of the outgoing partons and the remnants.";
  severity(maybeabort);
}
  
void PartonExtractor::dofinish() {
  partonBinInstances().clear();
  HandlerBase::dofinish();
}
