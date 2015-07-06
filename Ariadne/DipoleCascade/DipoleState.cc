// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DipoleState class.
//

#include "DipoleState.h"
#include "Parton.h"
#include "HardRemnant.h"
#include "ExtendedDipole.h"
#include "DISDipole.h"
#include "EMDipole.h"
#include "String.h"
#include "MECorrBase.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/EventRecord/SubProcess.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/ColourSinglet.h"
#include "ThePEG/EventRecord/ColourLine.h"
#include "ThePEG/EventRecord/Step.h"
#include "ThePEG/Utilities/ObjectIndexer.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include "ThePEG/Utilities/EnumIO.h"
#include "ThePEG/Utilities/UtilityBase.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Utilities/MaxCmp.h"
#include "ThePEG/PDF/PartonExtractor.h"
#include "ThePEG/Handlers/XComb.h"
#include "ThePEG/Config/algorithm.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "DipoleState.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Ariadne;

DipoleState::~DipoleState() {}

ClonePtr DipoleState::clone() const {
  return new_ptr(*this);
}

tPPtr DipoleState::singleResonance(const PPair & inc) const {
  if ( inc.first->children().size() == 1 &&
       inc.second->children().size() == 1 &&
       inc.first->children()[0] ==
       inc.second->children()[0] )
    return inc.first->children()[0];
  return tPPtr();
}

DipoleState::ProcessType DipoleState::getType(const PPair & inc) const {
  if ( LeptonMatcher::Check(inc.first->data()) &&
       LeptonMatcher::Check(inc.second->data()) ) {
    if ( singleResonance(inc) ) return annihilation;
    return leptonlepton;
  }
  if ( LeptonMatcher::Check(inc.first->data()) ) {
    if ( inc.second->previous() || inc.second->parents().size() > 0 )
      return leptonhadron;
    return DIPSYlephad;
  }
  if ( LeptonMatcher::Check(inc.second->data()) ) {
    if ( inc.first->previous() || inc.first->parents().size() > 0 )
      return hadronlepton;
    return DIPSYhadlep;
  }
    if ( inc.second->previous() || inc.second->parents().size() > 0 ||
	 inc.first->previous() || inc.first->parents().size() > 0 )
      return hadronhadron;
    return DIPSYhadhad;
}

pair<LorentzMomentum,LorentzMomentum>
DipoleState::setIncoming(const PPair & inc) {

  pair<LorentzMomentum,LorentzMomentum> incmom;
  theRemnants = tRemPair();
  if ( pType() == annihilation || pType() == DIPSYlephad ||
       pType() == DIPSYhadlep || pType() == DIPSYhadhad ) {
    theParticles = inc;
    incmom = make_pair(inc.first->momentum(), inc.second->momentum());
  } else {
    if ( pType() == hadronlepton || pType() == hadronhadron ) {
      theParticles.first = inc.first->parents()[0];
      incmom.first =
	lightCone(particles().first->momentum().plus(), 0.0*GeV);
    } else {
      Energy2 virt = 0.0*GeV2;
      if ( inc.first->children().size() > 1 ) {
	for ( int i = 0, N = inc.first->children().size(); i < N; ++i ) {
	  if ( LeptonMatcher::Check(inc.first->children()[i]->data()) )
	    theScatteredLeptons.first = inc.first->children()[i];
	  if ( inc.first->children()[i]->momentum().m2() < virt ) {
	    theParticles.first = inc.first->children()[i];
	    virt = particles().first->momentum().m2();
	  }
	}
	if ( virt < 0.0*GeV2 ) incmom.first = particles().first->momentum();
      } else {
	tPPtr p = inc.first->children()[0];
	for ( int i = 0, N = p->parents().size(); i < N; ++i )
	  if ( p->parents()[i]->momentum().m2() < virt ) {
	    theParticles.first = p->parents()[i];
	    virt = particles().first->momentum().m2();
	  }
	if ( virt < 0.0*GeV2 ) incmom.first = -particles().first->momentum();
	theScatteredLeptons.first = p;
      }
    }
    if ( pType() == leptonhadron || pType() == hadronhadron ) {
      theParticles.second = inc.second->parents()[0];
      incmom.second =
	lightCone(0.0*GeV, particles().second->momentum().minus());
    } else {
      Energy2 virt = 0.0*GeV2;
      if ( inc.second->children().size() > 1 ) {
	for ( int i = 0, N = inc.second->children().size(); i < N; ++i ) {
	  if ( LeptonMatcher::Check(inc.second->children()[i]->data()) )
	    theScatteredLeptons.second = inc.second->children()[i];
	  if ( inc.second->children()[i]->momentum().m2() < virt ) {
	    theParticles.second = inc.second->children()[i];
	    virt = particles().second->momentum().m2();
	  }
	}
	if ( virt < 0.0*GeV2 ) incmom.second = particles().second->momentum();
      } else {
	tPPtr p = inc.second->children()[0];
	for ( int i = 0, N = p->parents().size(); i < N; ++i )
	  if ( p->parents()[i]->momentum().m2() < virt ) {
	    theParticles.second = p->parents()[i];
	    virt = particles().second->momentum().m2();
	  }
	if ( virt < 0.0*GeV2 ) incmom.second = -particles().second->momentum();
	theScatteredLeptons.second = p;
      }	
    }
  }

  return incmom;

}

LorentzRotation DipoleState::init(tSubProPtr sub) {

  thePType = getType(sub->incoming());
  pair<LorentzMomentum,LorentzMomentum> incmom = setIncoming(sub->incoming());

  LorentzRotation rot = Utilities::getBoostToCM(incmom);

  sub->transform(rot);

  if ( pType() == leptonlepton || pType() == leptonhadron ) {
    MaxCmp<double> maxrap;
    for ( int i = 0, N = sub->outgoing().size(); i < N; ++i )
      if ( !LeptonMatcher::Check(sub->outgoing()[i]->data()) &&
	   maxrap(sub->outgoing()[i]->momentum().rapidity()) )
	theScatteredQuarks.first = sub->outgoing()[i];
  }
  if ( pType() == leptonlepton || pType() == hadronlepton ) {
    MaxCmp<double> maxrap;
    for ( int i = 0, N = sub->outgoing().size(); i < N; ++i )
      if ( !LeptonMatcher::Check(sub->outgoing()[i]->data()) &&
	   maxrap(-sub->outgoing()[i]->momentum().rapidity()) )
	theScatteredQuarks.second = sub->outgoing()[i];
  }

  // If incoming hadrons, store information in the incoming parton
  // objects.
  if ( pType() == hadronlepton || pType() == hadronhadron ) {
    theRemnants.first = create<SoftRemnant>();
    remnants().first->orig(sub->incoming().first);
    remnants().first->mu(getSoftMu(particles().first));
    remnants().first->alpha(getSoftAlpha(particles().first));
    remnants().first->pdf(handler()->firstPDF());
    remnants().first->parentData(particles().first->dataPtr());
    remnants().first->parentMomentum(rot*incmom.first);
    remnants().first->momentum() =
      remnants().first->parentMomentum() - remnants().first->momentum();
  }
  if ( pType() == leptonhadron || pType() == hadronhadron ) {
    theRemnants.second = create<SoftRemnant>();
    remnants().second->orig(sub->incoming().second);
    remnants().second->mu(getSoftMu(particles().second));
    remnants().second->alpha(getSoftAlpha(particles().second));
    remnants().second->pdf(handler()->secondPDF());
    remnants().second->parentData(particles().second->dataPtr());
    remnants().second->parentMomentum(rot*incmom.second);
    remnants().second->momentum() =
      remnants().second->parentMomentum() - remnants().second->momentum();
  }

  for ( int i = 0, N = sub->intermediates().size(); i < N; ++i )
    if ( sub->intermediates()[i] != particles().first &&
	 sub->intermediates()[i] != particles().second )
      hardSubSys().addIntermediate(sub->intermediates()[i]);

  tcParticleSet coloured;
  for ( int i = 0, N = sub->outgoing().size(); i < N; ++i )
    if ( sub->outgoing()[i]->coloured() )
      coloured.insert(sub->outgoing()[i]);
    else if ( sub->outgoing()[i] != scatteredLeptons().first &&
	      sub->outgoing()[i] != scatteredLeptons().second )
      hardSubSys().add(tPPtr(sub->outgoing()[i]));

  PPair rems;
  if ( remnants().first ) {
    tcPPtr p = remnants().first->orig();
    rems.first = p->data().CC()->produceParticle();
    if ( p->colourLine() ) p->colourLine()->addAntiColoured(rems.first);
    if ( p->antiColourLine() ) p->antiColourLine()->addColoured(rems.first);
    coloured.insert(rems.first);
  }
  if ( remnants().second ) {
    tcPPtr p = remnants().second->orig();
    rems.second = p->data().CC()->produceParticle();
    if ( p->colourLine() ) p->colourLine()->addAntiColoured(rems.second);
    if ( p->antiColourLine() ) p->antiColourLine()->addColoured(rems.second);
    coloured.insert(rems.second);
  }

  vector<ColourSinglet> singlets =
    ColourSinglet::getSinglets(coloured.begin(), coloured.end());

  // Now extract all strings, partons and connect them with dipoles.
  for ( int i = 0, N = singlets.size(); i < N; ++i )
    createString(singlets[i], false, rems);

  // Finally rotate so that all incoming partons are collinear with
  // the incoming hadrons.
  if ( pType() == hadronhadron ) {
    ThreeVector<double> b = hardSubSys().momentum().boostVector();
    b = ThreeVector<double>(b.x(), b.y(), 0.0);
    hardSubSys().transform(LorentzRotation(-b));
    remnants().first->momentum() =
      lightCone(remnants().first->parentMomentum().plus() +
		remnants().second->parentMomentum().plus() -
		hardSubSys().momentum().plus(), 0.0*GeV);
    remnants().second->momentum() =
      lightCone(0.0*GeV, remnants().first->parentMomentum().minus() +
		remnants().second->parentMomentum().minus() -
		hardSubSys().momentum().minus());

  }

  findMECorrs();

  return rot;
}

void DipoleState::createString(const ColourSinglet & sing, bool respectScale,
			       tPPair rems) {
  bool extended = false;
  if ( sing.nPieces() > 1 ) throw JunctionException()
    << "The Ariadne::CascadeHandler '" << handler()->name()
    << "' cannot handle junction strings." << Exception::runerror;
  bool gluonring = false;
  if ( sing.partons().size() < 2 ) throw StringException()
    << "The Ariadne::CascadeHandler '" << handler()->name()
    << "' found a string with less than two partons." << Exception::runerror;
  if ( sing.partons()[0]->hasAntiColour() &&
       sing.partons()[0]->hasColour() ) gluonring = true;

  StrPtr str = create<String>();
  vector<tParPtr> parts(sing.partons().size());
  for ( int j = 0, M = sing.partons().size(); j < M; ++j ) {
    if ( rems.first && sing.partons()[j] == rems.first ) {
      parts[j] = remnants().first;
      extended = true;
    }
    else if ( rems.second && sing.partons()[j] == rems.second ) {
      parts[j] = remnants().second;
      extended = true;
    }
    else if ( ( pType() == leptonhadron || pType() == leptonlepton ) &&
	      sing.partons()[j] == scatteredQuarks().first ) {
      HardRemPtr r = create<HardRemnant>();
      parts[j] = r;
      parts[j]->orig(sing.partons()[j]);
      r->mu(getHardMu(particles().first));
      r->alpha(getHardAlpha(particles().first));
      r->setQ2(abs(particles().first->momentum().m2()));
      hardSubSys().add(parts[j]);
      hardSubSys().hardRemnant(r, true);
    }
    else if ( ( pType() == hadronlepton || pType() == leptonlepton ) &&
	      sing.partons()[j] == scatteredQuarks().second ) {
      HardRemPtr r = create<HardRemnant>();
      parts[j] = r;
      parts[j]->orig(sing.partons()[j]);
      r->mu(getHardMu(particles().second));
      r->alpha(getHardAlpha(particles().second));
      r->setQ2(abs(particles().second->momentum().m2()));
      hardSubSys().add(parts[j]);
      hardSubSys().hardRemnant(r, false);
    }
    else {
      parts[j] = create(sing.partons()[j]);
    }
    parts[j]->string(str);
  }
  if ( !gluonring && sing.partons()[0]->hasAntiColour() )
    std::reverse(parts.begin(), parts.end());
  str->endpoints(parts.back(), *parts.begin());
  //  strings.insert(str);


  // Selecting a resonance will encode the information needed for a
  // dipole to generate emissions according to leading-order matrix
  // elements.
  tcPDPtr resonance = tcPDPtr();
  if ( sing.partons().size() == 2 ) {

    // Check if this is a final-state dipole steming from the decay of
    // a colour-singlet resonance.
    if ( sing.partons()[0]->parents().size() == 1 &&
	 sing.partons()[1]->parents().size() == 1 &&
	 sing.partons()[0]->parents()[0]->children().size() == 2 &&
	 sing.partons()[0]->parents()[0] == sing.partons()[1]->parents()[0] )
      resonance = sing.partons()[0]->parents()[0]->dataPtr();

    // Check if this is a remnant--remnant dipole from a leading-order
    // Drell-Yan production of a colour-singlet resonance.
    if ( ( ( parts[0] == remnants().first &&
	     parts[1] == remnants().second ) ||
	   ( parts[1] == remnants().first &&
	     parts[0] == remnants().second ) ) &&
	 ( sing.partons()[0]->children().size() == 1 &&
	   sing.partons()[1]->children().size() == 1 &&
	   sing.partons()[0]->children()[0]->parents().size() == 2 &&
	   sing.partons()[0]->children()[0] ==
	   sing.partons()[1]->children()[0] ) )
      resonance = sing.partons()[0]->children()[0]->dataPtr();

    // Check if this was a leading-order DIS-like event.
    tSoftRemPtr rp = dynamic_ptr_cast<tSoftRemPtr>(parts[0]);
    if ( !rp ) rp = dynamic_ptr_cast<tSoftRemPtr>(parts[1]);
    if ( rp && rp->mu() < 0.0*GeV ) {
      if ( member(rp->orig()->parents(), particles().first ) )
	resonance = particles().first->dataPtr();
      if ( member(rp->orig()->parents(), particles().second ) )
	resonance = particles().second->dataPtr();
    }
  }
  for ( int j = 1, M = sing.partons().size(); j < M; ++j ) {
    DipPtr dip;
    if ( dynamic_ptr_cast<tSoftRemPtr>(parts[j]) ||
	 dynamic_ptr_cast<tSoftRemPtr>(parts[j - 1]) ) {
      if ( dynamic_ptr_cast<tHardRemPtr>(parts[j]) ||
	   dynamic_ptr_cast<tHardRemPtr>(parts[j - 1]) )
	dip = create<DISDipole>();
      else
	dip = create<ExtendedDipole>();
    } else
      dip = create<Dipole>();
    dip->init(parts[j], parts[j - 1], respectScale);
    parts[j - 1]->iDip(dip);
    parts[j]->oDip(dip);
    dip->generateColourIndex();
    dip->resonance(resonance);
  }
  if ( gluonring ) {
    DipPtr dip;
    if ( parts[0] == remnants().first ||
	 parts[0] == remnants().second ||
	 parts.back() == remnants().first ||
	 parts.back() == remnants().second )
      dip = create<ExtendedDipole>();
    else
      dip = create<Dipole>();
    dip->init(parts[0], parts.back(), respectScale);
    parts.back()->iDip(dip);
    parts[0]->oDip(dip);
    dip->generateColourIndex();
    dip->resonance(resonance);
  }
  else if ( !extended && handler()->photonEmissions() ) {
    if( parts[0]->data().charged() && parts.back()->data().charged() ) {
      EMDipPtr emdip = create<EMDipole>();
      emdip->init(parts[0], parts.back(), respectScale);
      str->EMDip(emdip);
    }
  }

}

bool DipoleState::init(const tPVector & out) {
  // First extract all singlets.
  vector<ColourSinglet> singlets =
    ColourSinglet::getSinglets(out.begin(), out.end());
  if ( singlets.empty() ) return false;

  // Now extract all strings, partons and connect them with dipoles.
  for ( int i = 0, N = singlets.size(); i < N; ++i )
    createString(singlets[i], true);

  thePType = unknown;

  findMECorrs();

  return true;
}

Energy2 DipoleState::sTot() const {
  LorentzMomentum ptot = hardSubSys().momentum();
  if ( remnants().first ) ptot += remnants().first->momentum();
  if ( remnants().second ) ptot += remnants().second->momentum();
  return ptot.m2();
}

Energy DipoleState::getSoftMu(tcPPtr p) const {
  // *** ATTENTION *** Here we may want to give separate mus for
  // different hadrons.
  if ( BaryonMatcher::Check(p->data()) ) return handler()->softMu();
  else if ( MesonMatcher::Check(p->data()) ) return handler()->softMu();
  else return handler()->softMu();
}

Energy DipoleState::getHardMu(tcPPtr p) const {
  // *** ATTENTION *** Maybe we need some options here.
  return abs(p->momentum().mass());
}

double DipoleState::getSoftAlpha(tcPPtr p) const {
  // *** ATTENTION *** Here we may want to give separate alphas for
  // different hadrons.
  if ( BaryonMatcher::Check(p->data()) ) return handler()->softAlpha();
  else if ( MesonMatcher::Check(p->data()) ) return handler()->softAlpha();
  else return handler()->softAlpha();
}

double DipoleState::getHardAlpha(tcPPtr p) const {
  // *** ATTENTION *** Maybe we need some options here.
  return handler()->hardAlpha();
}

Energy2 DipoleState::select(Energy2 pt2min, Energy2 pt2max) {
  theSelected = tEmiPtr();
  Energy2 pt2 = 0.0*GeV2;
  Energy2 pt2sel = pt2min;
  for ( EmitterSet::iterator it = emitters.begin();
	it != emitters.end(); ++it ) {
    Energy2 pt2m = pt2max;
    if ( (**it).maxScale() > 0.0*GeV2 ) pt2m = (**it).maxScale();
    if ( (pt2 = (**it).generate(pt2min, pt2m)) > pt2sel ) {
      theSelected = *it;
      pt2sel = pt2;
    }
  }
  for ( BaseSet::iterator it = objects.begin(); it != objects.end(); ++it)
    (**it).untouch();
  hardSubSys().untouch();
  return pt2sel;
}

bool DipoleState::perform() {
  if ( !selected() ) return false;

  //Save four momentum before the emission.
  LorentzMomentum ptot = totalMomentum();

  if ( !selected()->performEmission() ) {
    selected()->touch();
    return false;
  }
  ++theNe;
  selected()->removeMECorr();

  //Check conservation of four momentum.
  Energy2 s = ptot.m2();
  ptot -= totalMomentum();
  if( (ptot.vect().mag2() + sqr(ptot.e()) )/s > 1.0e-16){
        throw MomentumException() << "Ariadne::DipoleState::perform: "
        << "Emission did not conserve four momentum."
        << Exception::eventerror;
  }

  //Check if all gluons are above the cutoff.
  // if(!checkGluonPT()){
  //   return false;
  // }

  return true;
}

void DipoleState::
fill(const tPVector & initial, const LorentzRotation & rot, tStepPtr step) {

  PVector final;
  Energy2 scale = sqr(handler()->pTCut());
  if ( selected() ) scale = max(scale, selected()->lastPT2());
  HardSubSys::PartonSet partons = hardSubSys().active();
  partons.insert(hardSubSys().produced().begin(),
		 hardSubSys().produced().end());
  for ( HardSubSys::PartonSet::iterator it = partons.begin();
	it != partons.end(); ++it ) {
    final.push_back((**it).produceParticle(rot));
    final.back()->scale(scale);
    tcParticleSet parents;
    (**it).getOriginalParents(inserter(parents));
    step->addDecayProduct(parents.begin(), parents.end(), final.back(), false);
  }

  for ( HardSubSys::PartonSet::iterator it = partons.begin();
	it != partons.end(); ++it )
    if ( (**it).oDip() )
      ColourLine::create((**it).next()->particle(), (**it).particle());

}

LorentzMomentum DipoleState::totalMomentum() const {
  LorentzMomentum ptot = hardSubSys().momentum();
  if(remnants().first){
    ptot += remnants().first->momentum();
  }
  if(remnants().second){
    ptot += remnants().second->momentum();
  }
  return ptot;
}

bool DipoleState::checkGluonPT(){
  HardSubSys::PartonSet active = hardSubSys().active();
  for ( HardSubSys::PartonSet::const_iterator it = active.begin();
	it != active.end(); ++it ){
    if((*it)->data().id() == 21 && ( (*it)->touched() ||
          (*it)->prev()->touched() || (*it)->next()->touched() )){
      if((*it)->invPT2() < sqr(handler()->pTCut())){
        return false;
      }
    }
  }
  return true;
}

bool DipoleState::constructHistory(int steps){
  theHistory.clear();
  if(steps < 0){
    return false;
  }
  allHistories(steps);
  selected(tEmiPtr());
  if(hasOrderedHistory(steps)){
    removeUnorderedHistories(steps);
  }
  return selectHistory(steps);
}

void DipoleState::allHistories(int steps){
  theHistory.clear();
  if(steps < 1){
    return;
  }
  for ( EmitterSet::const_iterator it = emitters.begin();
	it != emitters.end(); ++it ){
    Emitter::DipoleStateVector vec = (*it)->constructStep();
    while(!vec.empty()){
      vec.back()->theNe--;
      vec.back()->allHistories(steps - 1);
      vec.back()->findMECorrs();
      double prob;
      if(steps == 1){
	prob = vec.back()->selected()->emissionProbability();
      }
      else{
	prob = vec.back()->selected()->emissionProbability() * 
	  vec.back()->theHistory.sum();
      }
      if(prob > 0.0){
	theHistory.insert(prob, vec.back());
      }
      vec.pop_back();
    }
  }
}

bool DipoleState::hasOrderedHistory(int steps){
  if(steps < 1){
    return true;
  }
  for ( DipoleStateSelector::const_iterator it = theHistory.begin();
	it != theHistory.end(); ++it ){
    if((!selected() || selected()->lastPT2() < (*it).second->selected()->lastPT2()) &&
       (*it).second->hasOrderedHistory(steps - 1)){
      return true;
    }
  }
  return false;
}

bool DipoleState::removeUnorderedHistories(int steps){
  if(steps < 1){
    return true;
  }
  if(theHistory.empty()){
    return false;
  }
  DipoleStateSelector s;
  double lastweight = 0.0;
  for ( DipoleStateSelector::const_iterator it = theHistory.begin();
	it != theHistory.end(); ++it ){
    double weight = (*it).first - lastweight;
    lastweight = (*it).first;
    if((!selected() || selected()->lastPT2() < (*it).second->selected()->lastPT2()) &&
       (*it).second->removeUnorderedHistories(steps - 1)){
      if(steps > 1){
        weight = (*it).second->selected()->emissionProbability() * 
          (*it).second->theHistory.sum();
      }
      s.insert(weight, (*it).second);
    }
  }
  theHistory.swap(s);
  return (!theHistory.empty());
}

bool DipoleState::sudakovVeto(int steps){
  if(steps < 0){
    return true;
  }

  Energy2 pt2max;
  if(steps == 0){
    pt2max = sTot() / 4.0;
  }
  else{
    if(!theSelectedHistory){
      return true;
    }
    pt2max = theSelectedHistory->selected()->lastPT2();
    if(theSelectedHistory->sudakovVeto(steps-1)){
      theSelectedHistory = DipoleStatePtr();
      return true;
    }
  } 

  theSelectedHistory = DipoleStatePtr();

  //Selected is null for the highest order state. 
  //No veto should be done here.
  if(!selected()){
    return false;
  }

  for ( EmitterSet::iterator it = emitters.begin();
	it != emitters.end(); ++it ){
    (*it)->touch();
  }

  Energy2 pt2min = selected()->lastPT2();
  if( select(pt2min, pt2max) > pt2min){
    return true;
  }
  return false;
}

bool DipoleState::selectHistory(int steps){
  if(steps < 1){
    return true;
  }
  if(theHistory.empty()){
    return false;
  }
  theSelectedHistory = theHistory.select(handler()->rnd());
  theHistory.clear();

  if(theSelectedHistory->selectHistory(steps - 1)){
    if(selected() && selected()->lastPT2() > 
        theSelectedHistory->selected()->lastPT2()){
      theSelectedHistory->selected()->lastPT2(selected()->lastPT2());
    }
    return true;
  }
  else{
    theSelectedHistory = DipoleStatePtr();
    return false;
  }
}

double DipoleState::couplingProduct(int steps){
  if(steps < 1){
    return 1.0;
  }
  if(theSelectedHistory){
    return theSelectedHistory->selected()->coupling()*
      theSelectedHistory->couplingProduct(steps - 1);
  }
  else{
    return 0.0;
  }
}

double DipoleState::PDFRatioProduct(int steps){
  if(steps < 1){
    return 1.0;
  }
  if(theSelectedHistory){
    return theSelectedHistory->selected()->PDFRatio()*
      theSelectedHistory->PDFRatioProduct(steps - 1);
  }
  else{
    return 0.0;
  }
}

void DipoleState::
fill(tSubProPtr sub, const LorentzRotation & rot, tStepPtr step) {

  PVector hard;

  // Determine the scale to be assigned to the produced particles.
  Energy2 scale = sqr(handler()->pTCut());
  if ( selected() ) scale = max(scale, selected()->lastPT2());

  // First add all changed colour-singlet particles from the hard
  // subsystem to the step
  LorentzRotation trot = rot*hardSubSys().totalRotation();
  const tPVector & inter = hardSubSys().intermediates();
  const tPVector & inits = hardSubSys().initial();
  for ( int i = 0, N = inter.size(); i < N; ++i ) {
    hard.push_back(inter[i]->data().produceParticle(trot*inter[i]->momentum()));
    step->setCopy(inter[i], hard.back());
  }
  for ( int i = 0, N = inits.size(); i < N; ++i ) {
    hard.push_back(inits[i]->data().produceParticle(trot*inits[i]->momentum()));
    step->setCopy(inits[i], hard.back());
  }
  for ( int i = 0, N = inter.size(); i < N; ++i )
    for ( int ip = 0, Np = inter[i]->parents().size(); ip < Np; ++ip )
      if ( inter[i]->parents()[ip]->next() )
	step->addDecayProduct(inter[i]->parents()[ip]->final(),
			      inter[i]->final());
  for ( int i = 0, N = inits.size(); i < N; ++i )
    for ( int ip = 0, Np = inits[i]->parents().size(); ip < Np; ++ip )
      if ( inits[i]->parents()[ip]->next() )
	step->addDecayProduct(inits[i]->parents()[ip]->final(),
			      inits[i]->final());

  // Now add all produced colour-singlet particles from the hard
  // subsystem to the step
  for ( HardSubSys::PartonSet::iterator it = hardSubSys().produced().begin();
	it != hardSubSys().produced().end(); ++it ) {
    tcParticleSet parents;
    (**it).getOriginalParents(inserter(parents));
    hard.push_back((**it).produceParticle(rot));
    step->addDecayProduct(parents.begin(), parents.end(), hard.back());
  }

  // For the active particles in the cascade we first have to
  // specially treat the incoming partons.
  PPair incoming = sub->incoming();
  PPair newincs = incoming;
  if ( pType() == hadronlepton || pType() == hadronhadron )
    newincs.first = remnants().first->produceParticle(rot);
  if ( pType() == leptonhadron || pType() == hadronhadron )
    newincs.second = remnants().second->produceParticle(rot);


  // Then we can create all final state coloured particles.
  for ( HardSubSys::PartonSet::iterator it = hardSubSys().active().begin();
	it != hardSubSys().active().end(); ++it ) {
    tcParticleSet parents;
    (**it).getOriginalParents(inserter(parents));
    if ( parents.size() == 1 ) continue;
    hard.push_back((**it).produceParticle(rot));
    hard.back()->scale(scale);
    step->addDecayProduct(parents.begin(), parents.end(), hard.back(), false);
  }
  for ( HardSubSys::PartonSet::iterator it = hardSubSys().active().begin();
	it != hardSubSys().active().end(); ++it ) {
    tcParticleSet parents;
    (**it).getOriginalParents(inserter(parents));
    if ( parents.size() != 1 ) continue;
    hard.push_back((**it).produceParticle(rot));
    hard.back()->scale(scale);
    tcPPtr prev = *parents.begin();
    step->setCopy(prev, hard.back());
    for ( int ip = 0, Np = prev->parents().size(); ip < Np; ++ip )
      if ( prev->parents()[ip]->next() )
	step->addDecayProduct(prev->parents()[ip]->final(), hard.back(), false);
  }

  // Connect all colours and take special care of the remnants.
  for ( HardSubSys::PartonSet::iterator it = hardSubSys().active().begin();
	it != hardSubSys().active().end(); ++it )
    if ( (**it).oDip() && (**it).next() != remnants().first &&
	 (**it).next() != remnants().second )
      ColourLine::create((**it).next()->particle(), (**it).particle());
  if ( remnants().first ) {
    if ( remnants().first->oDip() ) {
      ColinePtr cl = new_ptr(ColourLine());
      cl->addColoured(remnants().first->particle());
      cl->addColoured(remnants().first->next()->particle(),
		      dynamic_ptr_cast<tSoftRemPtr>(remnants().first->next()));
    }
    if ( remnants().first->iDip() ) {
      ColinePtr cl = new_ptr(ColourLine());
      cl->addAntiColoured(remnants().first->particle());
      cl->addColoured(remnants().first->prev()->particle(),
		      !dynamic_ptr_cast<tSoftRemPtr>(remnants().first->prev()));
    }
  }
  if ( remnants().second ) {
    if ( remnants().second->oDip() &&
	 !dynamic_ptr_cast<tSoftRemPtr>(remnants().second->next()) ) {
      ColinePtr cl = new_ptr(ColourLine());
      cl->addColoured(remnants().second->particle());
      cl->addColoured(remnants().second->next()->particle());
    }
    if ( remnants().second->iDip() &&
	 !dynamic_ptr_cast<tSoftRemPtr>(remnants().second->prev()) ) {
      ColinePtr cl = new_ptr(ColourLine());
      cl->addAntiColoured(remnants().second->particle());
      cl->addAntiColoured(remnants().second->prev()->particle());
    }
  }

  sub->transform(rot);

  // Finally we may have to 're-extract' the incoming particles and boost
  // the produced hard subsystem accordingly.
  if ( remnants().first || remnants().second ) {
    PartonExtractor & pex = *(handler()->lastExtractor());
    Lorentz5Momentum p1 = newincs.first->momentum();
    Lorentz5Momentum p2 = newincs.second->momentum();
    PBIPair newbins = pex.newRemnants(incoming, newincs, step);
    LorentzRotation tot =
      pex.boostRemnants(newbins, p1, p2,
			newbins.first && newbins.first->incoming(),
			newbins.second && newbins.second->incoming());

    // Begin Testing
//     LorentzMomentum k1 = newbins.first->parton()->momentum();
//     LorentzMomentum k2 = newbins.second->parton()->momentum();
//     LorentzMomentum pp = hard[0]->momentum();
//     double yi = (p1 + p2).rapidity();
//     double yf = (k1 + k2).rapidity();
//     double ptzi = pp.perp()/GeV;
//     double yzi = pp.rapidity();
//     double pth = (k1 + k2).perp()/GeV;

//     LorentzRotation boostback((k1 + k2).boostVector());
//     LorentzRotation boost(-(p1 + p2).boostVector());
//     k1.transform(boostback.inverse());
//     k2.transform(boostback.inverse());
//     p1.transform(boost);
//     pp.transform(boost);

//     pp.rotateZ(-p1.phi());
//     pp.rotateY(-p1.theta());
//     pp.rotateZ(p1.phi());

//     pp.rotateZ(-k1.phi());
//     pp.rotateY(k1.theta());
//     pp.rotateZ(k1.phi());

//     pp.transform(boostback);
    
//     double ptzf = pp.perp()/GeV;
//     double yzf = pp.rapidity();

    // End Testing

    Utilities::transform(hard.begin(), hard.end(), tot);

    // Begin Testing

//     LorentzMomentum ppp = hard[0]->momentum();

    // End Testing

    if ( remnants().first ) {
      newincs.first->addChild(incoming.first);
      step->addIntermediate(newincs.first);
    }
    if ( remnants().second ) {
      newincs.second->addChild(incoming.second);
      step->addIntermediate(newincs.second);
    }
  }

}

bool DipoleState::replaceRemnant(tSoftRemPtr old, tSoftRemPtr p) {
  if ( remnants().first == old )
    theRemnants.first = p;
  else if ( remnants().second == old )
    theRemnants.second = p;
  else
    return false;
  remove(tParPtr(old));
  return true;
}

tcMECPtr DipoleState::findMECorr(tcEmiPtr dipole) const {
  for ( int i = 0, N = handler()->MECorrectors().size(); i < N; ++i )
    if ( handler()->MECorrectors()[i]->canHandle(dipole, this) )
      return handler()->MECorrectors()[i];
  return tcMECPtr();
}

void DipoleState::findMECorrs() {
  for ( EmitterSet::iterator it = emitters.begin();
	it != emitters.end(); ++it )
    (**it).theMECorr = findMECorr(*it);
}

void DipoleState::rebind(const TranslationMap & trans) {
  BaseSet old;
  old.swap(objects);
  objects.clear();
  trans.translate(inserter(objects), old.begin(), old.end());
  EmitterSet eold;
  eold.swap(emitters);
  emitters.clear();
  trans.translate(inserter(emitters), eold.begin(), eold.end());
  StringSet sold;
  sold.swap(strings);
  strings.clear();
  trans.translate(inserter(strings), sold.begin(), sold.end());
  theSelected = trans.translate(theSelected);
  theRemnants.first = trans.translate(theRemnants.first);
  theRemnants.second = trans.translate(theRemnants.second);
  theHardSubSys.rebind(trans);
  emiindx.clear();
  parindx.clear();
  strindx.clear();
}


DipoleStatePtr DipoleState::fullclone() const {
  TranslationMap trans;
  return fullclone(trans);
}

DipoleStatePtr DipoleState::fullclone(TranslationMap & trans) const {
  DipoleStatePtr copy = dynamic_ptr_cast<DipoleStatePtr>(clone());
  vector<ClonePtr> copies;

  copies.push_back(trans[tcDipoleStatePtr(this)] = copy);
  for ( BaseSet::const_iterator it = objects.begin();
	it != objects.end(); ++it )
    copies.push_back(trans[*it] = (**it).clone());

  for ( int i = 0, N = copies.size(); i < N; ++i ) copies[i]->rebind(trans);

  return copy;
}

DipoleStatePtr DipoleState::preclone(TranslationMap & trans) const {
  DipoleStatePtr copy = dynamic_ptr_cast<DipoleStatePtr>(clone());
  trans[tcDipoleStatePtr(this)] = copy;
  for ( BaseSet::const_iterator it = objects.begin();
	it != objects.end(); ++it ) trans[*it] = (**it).clone();

  return copy;
}

void DipoleState::postclone(const TranslationMap & trans) const {
  for ( TranslationMap::const_iterator it = trans.map().begin();
	it != trans.map().end(); ++it ) it->second->rebind(trans);
}



ParPtr DipoleState::create(PDPtr type) {
  ParPtr par = new_ptr(Parton());
  objects.insert(objects.end(), par);
  par->handler(handler());
  par->state(this);
  par->data(type);
  hardSubSys().add(tParPtr(par));
  return par;
}

ParPtr DipoleState::create(tcPPtr p) {
  ParPtr par = new_ptr(Parton());
  objects.insert(objects.end(), par);
  par->handler(handler());
  par->state(this);
  par->orig(p);
  hardSubSys().add(tParPtr(par));
  return par;
}

void DipoleState::remove(tParPtr p) {
  hardSubSys().remove(p);
}

void DipoleState::removeEmitter(tEmiPtr p) {
  EmitterSet::iterator it = find(emitters, p);
  if ( it != emitters.end() ) emitters.erase(it);
}

void DipoleState::removeString(tStrPtr str) {
  StringSet::iterator it = find(strings, str);
  if ( it != strings.end() ) strings.erase(it);
}

void DipoleState::persistentOutput(PersistentOStream & os) const {
  os << emitters << strings << objects << theSelected << theRemnants;
  hardSubSys().persistentOutput(os);
  os << oenum(thePType) << theNe << theParticles << theScatteredLeptons
     << theScatteredQuarks;
}

void DipoleState::persistentInput(PersistentIStream & is, int version) {
  is >> emitters >> strings >> objects >> theSelected >> theRemnants;
  hardSubSys().persistentInput(is, version);
  is >> ienum(thePType) >> theNe >> theParticles >> theScatteredLeptons
     >> theScatteredQuarks;
}

ClassDescription<DipoleState> DipoleState::initDipoleState;
// Definition of the static class description member.

void DipoleState::Init() {}

void DipoleState::debugme() const {
  CascadeBase::debugme();
  cerr << "DipoleState:" << endl;
  LorentzMomentum sum;
  for ( StringSet::iterator it = strings.begin(); it != strings.end(); ++it ) {
    (**it).debug();
    tParPtr curr = (**it).endpoints().first;
    tParPtr last = (**it).endpoints().second;
    do {
      sum += curr->momentum();
      curr->debug();
      cerr << endl;
      if ( curr->oDip() ) {
	curr->oDip()->debug();
	cerr << endl;
      }
      if ( curr == last ) break;
      curr = curr->next();
    } while ( curr );
    cerr << endl;
  }

  if ( !hardSubSys().initial().empty() ||
       !hardSubSys().produced().empty() ) {
    cerr << "Colourless:" << endl;
    for ( int i = 0, N = hardSubSys().initial().size(); i < N; ++i ) {
      LorentzMomentum p = hardSubSys().totalRotation()*
	hardSubSys().initial()[i]->momentum();
      sum += p;
      cerr << setw(6) << hardSubSys().initial()[i]->id()
	   << setw(22) << p.x()/GeV
	   << setw(9) << p.y()/GeV
	   << setw(9) << p.z()/GeV
	   << setw(9) << p.e()/GeV
	   << setw(9) << p.m()/GeV << endl;
    }
    for ( HardSubSys::PartonSet::iterator it = hardSubSys().produced().begin();
	  it != hardSubSys().produced().end(); ++it ) {
      sum += (**it).momentum();
      cerr << setw(6) << (**it).data().id()
	   << setw(22) << (**it).momentum().x()/GeV
	   << setw(9) << (**it).momentum().y()/GeV
	   << setw(9) << (**it).momentum().z()/GeV
	   << setw(9) << (**it).momentum().e()/GeV
	   << setw(9) << (**it).momentum().m()/GeV << endl;
    }
    cerr << endl;
  }

  cerr << "Sum of momenta:" << setprecision(3)
       << setw(13) << ( abs(sum.x()) > MeV? sum.x(): 0.0*GeV )/GeV
       << setw(9) << ( abs(sum.y()) > MeV? sum.y(): 0.0*GeV )/GeV
       << setw(9) << ( abs(sum.z()) > MeV? sum.z(): 0.0*GeV )/GeV
       << setw(9) << sum.e()/GeV
       << setw(9) << sum.m()/GeV << endl
       << setprecision(12) << sum.e()/GeV << endl << endl;

  if(theSelectedHistory){
    cerr << "The Selected History:" << endl;
    theSelectedHistory->debug();
  }

  if(theHistory.size() > 0){
    cerr << "All Histories:" << endl;
    for(DipoleStateSelector::const_iterator it = theHistory.begin();
        it != theHistory.end(); it++){
      tEmiPtr sel = it->second->selected();
      cerr << "Emission probability: " << sel->emissionProbability() << endl;
      cerr << "genVar:";
      for(vector<double>::const_iterator vit =  sel->genVar.begin();
          vit != sel->genVar.end(); vit++){
        cerr << " " << *vit;
      }
      cerr << endl;

//       tExDipPtr edip = dynamic_ptr_cast<tExDipPtr>(sel);
//       if(edip){
//         cerr << "Emission type: " << edip->eClass() << endl;
//       }

      it->second->debug();
    }
  }
}

bool DipoleState::checkIntegrety(){
  for ( EmitterSet::iterator it = emitters.begin();
      it != emitters.end(); ++it ) {
    if( ! (*it)->checkIntegrety() ){
      return false;
    }
  }
  return true;
}

int DipoleState::index(tcCascadeBasePtr o) {
  emiindx(tcEmiPtr());
  parindx(tcParPtr());
  strindx(tcStrPtr());
  if ( dynamic_ptr_cast<tcEmiPtr>(o) )
    return emiindx(dynamic_ptr_cast<tcEmiPtr>(o));
  if ( dynamic_ptr_cast<tcParPtr>(o) )
    return parindx(dynamic_ptr_cast<tcParPtr>(o));
  if ( dynamic_ptr_cast<tcStrPtr>(o) )
    return strindx(dynamic_ptr_cast<tcStrPtr>(o));
  return 0;
}

