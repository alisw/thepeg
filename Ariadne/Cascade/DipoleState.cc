// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DipoleState class.
//

#include "DipoleState.h"
#include "Resonance.h"
#include "ResonanceParton.h"
#include "RemnantParton.h"
#include "AriadneHandler.h"
#include "EMDipole.h"
#include "QCDDipole.h"
#include "StateDipole.h"
#include "Junction.h"
#include "PartonTraits.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/SubProcess.h"
#include "ThePEG/EventRecord/Step.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/PDF/PartonExtractor.h"
#include "ThePEG/PDF/PartonBinInstance.h"
#include "ThePEG/Config/algorithm.h"
#include "ThePEG/Utilities/UtilityBase.h"
#include "ThePEG/Utilities/MaxCmp.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Utilities/UtilityBase.h"
#include "ThePEG/Utilities/DebugItem.h"
#include "ThePEG/Utilities/StringUtils.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "Models/FSGluonEmitter.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Ariadne5;

DipoleState::DipoleState(tPPair inc): theIncoming(inc), theNe(0) {
  dipindx(tcDBPtr());
  parindx(tcParPtr());
}

DipoleState::~DipoleState() {}

ClonePtr DipoleState::clone() const {
  return new_ptr(*this);
}

void DipoleState::setup(SubProcess & sub) {

  set<tPPtr> left(sub.outgoing().begin(), sub.outgoing().end());

  // First check if this is a DIS-like event.
  pair<tRemParPtr,tRemParPtr> leptons = handler()->findDISLeptons(sub, *this);
  pair<tRemParPtr,tRemParPtr> quarks =
    handler()->findDISQuarks(leptons, sub, *this);
  if ( quarks.first ) {
    quarks.first->setupHard(1);
    leptons.first->setupHard(1);
    theRemnants.first.push_back(quarks.first);
    theRemnants.first.push_back(leptons.first);
    addHadronicFS(quarks.first);
    addHardFS(leptons.first);
    left.erase(quarks.first->orig());
    left.erase(leptons.first->orig());
  }
  if ( quarks.second ) {
    quarks.second->setupHard(-1);
    leptons.second->setupHard(-1);
    theRemnants.second.push_back(quarks.second);
    theRemnants.second.push_back(leptons.second);
    addHadronicFS(quarks.second);
    addHardFS(leptons.second);
    left.erase(quarks.second->orig());
    left.erase(leptons.second->orig());
  }

  // Now setup the rest of the remnants.
  tPBIPtr pb = handler()->lastXComb().partonBinInstance(sub.incoming().first);
  while ( pb && pb->incoming() ) {
    tRemParPtr rem = create<RemnantParton>();
    rem->setup(*pb, 1);
    addFS(rem);
    theRemnants.first.push_back(rem);
    pb = pb->incoming();
  }
  pb = handler()->lastXComb().partonBinInstance(sub.incoming().second);
  while ( pb && pb->incoming() ) {
    tRemParPtr rem = create<RemnantParton>();
    rem->setup(*pb, -1);
    addFS(rem);
    theRemnants.second.push_back(rem);
    pb = pb->incoming();
  }

  // After that we should identify all resonances.
  tPVector resv = handler()->resonances(sub);
  theResonances.resize(resv.size());
  for ( int i = 0, N = resv.size(); i < N; ++i ) {
    tResPtr res = create<Resonance>();
    res->orig(resv[i]);
    res->decaySystem(i + 1);
    theResonances[i] = res;
  }
  for ( int i = 0, N = resv.size(); i < N; ++i ) {
    tResPtr res = theResonances[i];
    if ( resv[i]->parents().size() == 1 ) {
      tPVector::iterator pit = find(resv, resv[i]->parents()[0]);
      if ( pit != resv.end() ) 
	res->parentResonance(theResonances[pit - resv.begin()]);
    }
    // Children in the final state should be of type ResonanceParton.
    for ( int ic = 0, Nc = resv[i]->children().size(); ic < Nc; ++ic ) {
      set<tPPtr>::iterator pit = find(left, resv[i]->children()[ic]);
      if ( pit != left.end() ) {
	tResParPtr rp = create<ResonanceParton>();
	rp->orig(*pit);
	rp->setResonance(res);
	addHadronicFS(rp);
	left.erase(pit);
      }
    }
  }

  // If there are decayed partons left, they are not to be considered
  // resonances, so we remove them and add their kids instead.
  fixupDecays(left);

  // What is left is just ordinary partons.
  setup(left);

  // Save a pointer to the SubProcess
  subprocess = &sub;

}

bool DipoleState::addChildren(tPVector & v, tPPtr p) {
  if ( !p->decayed() && !p->next() ) return false;
  if ( p->next() ) {
    if ( !addChildren(v, p->next()) ) v.push_back(p->next());
  } else {
    for ( int i = 0, N = p->children().size(); i < N; ++i )
      if ( !addChildren(v, p->children()[i]) ) v.push_back(p->children()[i]);
  }
  return true;
}

void DipoleState::fixupDecays(set<tPPtr> & out) {
  tPVector removed;
  tPVector added;
  for ( set<tPPtr>::iterator it = out.begin(); it != out.end(); ++it )
    if ( addChildren(added, *it) ) removed.push_back(*it);
  for ( int i = 0, N = removed.size(); i < N; ++i ) out.erase(removed[i]);
  out.insert(added.begin(), added.end());
}


void DipoleState::setup(const set<tPPtr> & out) {
  subprocess = tSubProPtr();

  // Here we only have ordinary partons, life is simple.
  for ( set<tPPtr>::iterator pit = out.begin(); pit != out.end(); ++pit ) {
    tParPtr p = create<Parton>();
    p->orig(*pit);
    addHadronicFS(p);
  }

  // Now identify all QCD dipoles.
  vector<tQCDPtr> dips = handler()->findQCDDipoles(*this);
  active.insert(dips.begin(), dips.end());

  // And the EM dipoles.
  vector<tEMDipPtr> emdips = handler()->findEMDipoles(*this);
  active.insert(emdips.begin(), emdips.end());

  // Finally we add a StateDipole to be used by global models.
  
  active.insert(create<StateDipole>());

  sumTotalMomentum();

}

Energy DipoleState::select(Energy rhomin, Energy rhomax) {
  theSelected = EmPtr();
  //  Energy rho = ZERO;
  MaxCmp<Energy> rhosel;
  for ( set<tDBPtr>::iterator it = active.begin(); it != active.end(); ++it ) {
    pair<set<EmissionGenerator>::iterator, bool> ins =
      generators.insert(EmissionGenerator(*it));
    const EmissionGenerator & gen = *ins.first;
    if ( ins.second ) gen.init();
    else gen.reinit();
    if ( rhosel(gen.generate(rhomin, rhomax)) && gen.emission ) {
      if ( !gen.emission->geno ) gen.emission->geno = theNe + 1;
      theSelected = gen.emission;
      rhomin = gen.rho();
    }
  }
  for ( BaseSet::iterator it = objects.begin(); it != objects.end(); ++it)
    (**it).untouch();
  if ( selected() ) emissions.push_back(selected());
  untouch();
  return rhosel;
}

bool DipoleState::perform() {
  static DebugItem printsteps("Ariadne5::PrintSteps", 6);
  touch();
  if ( !selected() ) return false;

  if ( printsteps )
    cerr << "Ariadne5::PrintStep: Performing emission " << theNe + 1
	 << " of type " << typeid(*selected()).name() << endl;

  // Save four momentum before the emission.
  LorentzMomentum ptot = totalMomentum();

  if ( !selected()->perform() ) {
    selected()->dipole->touch();
    if ( printsteps )
      cerr << "Ariadne5::PrintStep: Emission of type "
	   << typeid(*selected()).name() << " failed" << endl;
    return false;
  }

  // Inform other partons that an emission was performed.
  for ( set<tParPtr>::const_iterator it = finalState().begin();
	it != finalState().end(); ++it )
    (**it).notify(*selected());

  ++theNe;

  selected()->emno = theNe;

  // Check conservation of four momentum.
  Energy2 s = ptot.m2();
  ptot -= sumTotalMomentum();
  if( (ptot.vect().mag2() + sqr(ptot.e()) )/s > 1.0e-16) 
    throw MomentumException()
      << "Ariadne5::DipoleState::perform: Emission of type "
      << StringUtils::typeName(typeid(*selected()))
      << " did not conserve four momentum." << Exception::eventerror;

  // Do additional consistency checks.
  return selected()->failsafe || 
    Current<AriadneHandler>()->checkState(*this, selected());
}

void DipoleState::fill(Step & step) const {
  set<tParPtr> done;
  // Determine the scale to be assigned to the produced particles.
  Energy2 scale = sqr(handler()->pTCut());
  if ( selected() ) scale = max(scale, sqr(selected()->rho));

  // First create all the ThePEG::Particles
  for ( set<tParPtr>::const_iterator it = finalState().begin();
	it != finalState().end(); ++it )
    if ( (**it).produceParticle() && (**it).particle() )
      (**it).particle()->scale(scale);

  // Then it's time to fix all colour lines.
  vector<ColinePtr> colines;
  for ( set<tDBPtr>::const_iterator it = active.begin();
	it != active.end(); ++it )
    if ( ColinePtr cl = (**it).colourLine() ) colines.push_back(cl);

  // Now collect all relevant particles (ie. ell except soft
  // remnants), and find all original final state partons that have actually
  // radiated.
  set<tPPtr> radiated;
  vector<PPtr> final;
  for ( set<tParPtr>::const_iterator it = hardFS().begin();
	it != hardFS().end(); ++it ) {
    if ( !(**it).orig() ) (**it).getOriginalParents(inserter(radiated));
    if ( tPPtr p = (**it).particle() ) final.push_back(p);
  }

  // Now we can fix all mother-daughter relationships and add new
  // particles to the step.
  for ( set<tParPtr>::const_iterator it = hardFS().begin();
	it != hardFS().end(); ++it )
    if ( (**it).orig() && !member(radiated, (**it).orig()) &&
	 (**it).adoptedOriginals.empty() )
      step.setCopy((**it).orig(), (**it).particle());
    else {
      set<tPPtr> parents;
      (**it).getOriginalParents(inserter(parents));
      step.addDecayProduct(parents.begin(), parents.end(),
			   (**it).particle(), false);
    }

  // Now add all intermediate resonances them to the step and fix up
  // the parenthood.
  for ( int i = 0, N = resonances().size(); i < N; ++i ) {
    final.push_back(resonances()[i]->produceParticle());
    step.setCopy(resonances()[i]->orig(), resonances()[i]->particle());
    if ( tColinePtr cl = resonances()[i]->orig()->colourLine() )
      cl->addColoured(resonances()[i]->particle());
    if ( tColinePtr cl = resonances()[i]->orig()->antiColourLine() )
      cl->addAntiColoured(resonances()[i]->particle());
  }
  for ( int i = 0, N = resonances().size(); i < N; ++i ) {
    tPPtr orig = resonances()[i]->orig();
    for ( int ic = 0, Nc = orig->children().size(); ic < Nc; ++ic ) {
      tPPtr child = orig->children()[ic];
      if ( child->next() )
	step.addDecayProduct(orig->final(), child->final());
      else if ( child->children().size() )
	step.addDecayProduct(orig->final(), child->children().begin(),
			      child->children().end(), false);
    }
  }


  if ( !subprocess ) return;

  // Now start with the innermost soft remnants and go outwards and
  // reextract the incoming particles.
  // First skip all hard remnants.
  vector<tRemParPtr>::const_iterator rit1 = remnants().first.begin();
  while ( rit1 != remnants().first.end() && (**rit1).hard() ) ++rit1;
  vector<tRemParPtr>::const_iterator rit2 = remnants().second.begin();
  while ( rit2 != remnants().second.end() && (**rit2).hard() ) ++rit2;
  PartonExtractor & pex = *(handler()->lastExtractor());
  while ( rit1 != remnants().first.end() || rit2 != remnants().second.end() ) {
    // If there are no soft remnants on one side oldinc=newinc is
    // taken to be the colliding particle and nothing will happen on
    // that side.
    PPair oldinc = subprocess->incoming();
    PPair newinc = subprocess->incoming();
    if ( rit1 != remnants().first.end() ) {
      oldinc.first = (**rit1).originalExtracted();
      newinc.first = (**rit1).extracted();
    }
    if ( rit2 != remnants().second.end() ) {
      oldinc.second = (**rit2).originalExtracted();
      newinc.second = (**rit2).extracted();
    }
    Lorentz5Momentum p1 = newinc.first->momentum();
    Lorentz5Momentum p2 = newinc.second->momentum();

    // Re-extract the remnants and calculate the boost needed to put
    // them on-shell.
    PBIPair newbins = pex.newRemnants(oldinc, newinc, &step);
    LorentzRotation tot =
      pex.boostRemnants(newbins, p1, p2,
			newbins.first && newbins.first->incoming(),
			newbins.second && newbins.second->incoming());

    // Perform the boost also for the final state.
    Utilities::transform(final.begin(), final.end(), tot);

    // Add the new remnants to the final state.
    if ( rit1 != remnants().first.end() ) {
      if ( newinc.first != oldinc.first ) {
	newinc.first->addChild(oldinc.first);
	step.addIntermediate(newinc.first);
      }
      final.insert(final.end(), newbins.first->remnants().begin(),
		   newbins.first->remnants().end());
      ++rit1;
    }
    if ( rit2 != remnants().second.end() ) {
      if ( newinc.second != oldinc.second ) {
	newinc.second->addChild(oldinc.second);
	step.addIntermediate(newinc.second);
      }
      final.insert(final.end(), newbins.second->remnants().begin(),
		   newbins.second->remnants().end());
      ++rit2;
    }
  }
}

const Lorentz5Momentum & DipoleState::sumTotalMomentum() {
  theTotalMomentum = Utilities::sumMomentum(finalState());
  theHardMomentum = Utilities::sumMomentum(hardFS());
  theHadronicMomentum = Utilities::sumMomentum(hadronicFS());
  theTotalMomentum.rescaleMass();
  theHardMomentum.rescaleMass();
  theHadronicMomentum.rescaleMass();
  return theTotalMomentum;
}

pair<double,int> DipoleState::lambdaMeasure(Energy2 scale) const {
  pair<double,int> lam (0.0, 0);
  for ( set<tDBPtr>::const_iterator it = activeDipoles().begin();
	it != activeDipoles().end(); ++it )
    if ( tQCDPtr d = dynamic_ptr_cast<tQCDPtr>(*it) ) {
      Energy2 m2 = d->sdip();
      if ( m2 > ZERO ) lam.first += log(m2/scale);
      ++lam.second;
    }
  return lam;
}

int DipoleState::crossings() const {
  int n = 0;
  for ( set<tDBPtr>::const_iterator it = activeDipoles().begin();
	it != activeDipoles().end(); ++it )
    if ( tQCDPtr d = dynamic_ptr_cast<tQCDPtr>(*it) )
      if ( d->iPart()->momentum().z()*d->oPart()->momentum().z() < ZERO ) ++n;
  return n;
}

int DipoleState::gluonsBelow() const {
  int n = 0;
  Energy2 pt2cut = sqr(Current<AriadneHandler>()->pTCut());
  for ( set<tDBPtr>::const_iterator it = activeDipoles().begin();
	it != activeDipoles().end(); ++it )
    if ( tQCDPtr d = dynamic_ptr_cast<tQCDPtr>(*it) )
      if ( d->next() && d->iPart() != d->next()->oPart() &&
	   EmitterBase::invPT2(d->iPart(), d->oPart(), d->next()->oPart())
	   < pt2cut ) ++n;
  return n;
}

double DipoleState::folding() const {
  double sumdy = 0.0;
  MaxCmp<double> ymax;
  MinCmp<double> ymin;
  for ( set<tDBPtr>::const_iterator it = activeDipoles().begin();
	it != activeDipoles().end(); ++it )
    if ( tQCDPtr d = dynamic_ptr_cast<tQCDPtr>(*it) ) {
      sumdy += abs(d->iPart()->momentum().rapidity() -
		   d->oPart()->momentum().rapidity());
      ymax(d->iPart()->momentum().rapidity());
      ymax(d->oPart()->momentum().rapidity());
      ymin(d->iPart()->momentum().rapidity());
      ymin(d->oPart()->momentum().rapidity());
    }      
  return sumdy/(ymax.value() - ymin.value());
}


void DipoleState::purgeGluons(Energy cut) {
  // In the first round we check all gluons. After that we only chck
  // those who where bleow in the previous round or those who were
  // affected of the previous reabsorption.
  set<tParPtr> maybelow;
  bool donefirst = false;

  while ( true ) { // Loop until all gluons are above cut.
                   // Choos the softest gluon in each round.
    MinCmp<Energy2,tQCDPtr> sel(sqr(cut));
    for ( set<tDBPtr>::const_iterator it = activeDipoles().begin();
	  it != activeDipoles().end(); ++it )
      if ( tQCDPtr d = dynamic_ptr_cast<tQCDPtr>(*it) ) {
	if ( !d->next() ) continue;
	if ( d->iPart() == d->next()->oPart() ) continue;
	Energy2 pt2 = EmitterBase::invPT2(d->iPart(), d->oPart(), d->next()->oPart());
	if ( d->prev() && d->prev()->colourIndex() == d->next()->colourIndex() && 
	     d->next()->next() && d->colourIndex() == d->next()->next()->colourIndex() ) continue;
	if ( donefirst && !member(maybelow, d->oPart()) ) continue;
	sel(pt2, d);
	if ( pt2 < sqr(cut) ) // If below, consider it also the next round
	  maybelow.insert(d->oPart());
	else
	  maybelow.erase(d->oPart());
      }
    if ( !sel.index() ) break;

    QCDDipole & d = *sel.index();
    purgeGluon(d);

    // Flag affected gluons so that they will be checked next round
    donefirst = true;
    if ( d.prev() ) maybelow.insert(d.prev()->iPart());
    maybelow.insert(d.iPart());
    maybelow.insert(d.oPart());
    if ( d.next() ) maybelow.insert(d.next()->oPart());

  }

}

void DipoleState::purgeGluon(QCDDipole & d) {
  ParPtr g = d.oPart();
  // First find the momenta of the connecting partons after removal.
  Lorentz5Momentum p1 = d.iPart()->momentum();
  Lorentz5Momentum p2 = d.oPart()->momentum();
  Lorentz5Momentum p3 = d.next()->oPart()->momentum();
  Energy2 S = (p1 + p2 + p3).m2();
  LorentzRotation R = Utilities::boostToCM(makeTriplet(&p1, &p2, &p3));
  Energy W = sqrt(S);
  double x1 = 2.0*p1.e()/W;
  double x3 = 2.0*p3.e()/W;
  bool nr1 = d.iPart()->isG() || dynamic_ptr_cast<tRemParPtr>(d.iPart());
  bool nr3 = d.next()->oPart()->isG() || dynamic_ptr_cast<tRemParPtr>(d.next()->oPart());
  double Psi = Constants::pi - p3.theta();
  double beta = 0.0;
  if ( ( nr1 && nr3 ) || (!nr1 && !nr3 ) )
    beta = Psi*sqr(x3)/(sqr(x1) + sqr(x3)); // minimize pt
  else if ( nr1 )
    beta = Psi;
  R.rotateY(-beta);
  R.invert();
  SimplePhaseSpace::CMS(p1, p3, S, 1.0, 0.0);
  p1.transform(R);
  p3.transform(R);

  // Check which dipole should be removed, the current or the
  // next. Take the smallest one if colourwise possible.
  bool remcur = d.sdip() < d.next()->sdip();
  if ( remcur && d.prev() &&
       d.prev()->colourIndex() == d.next()->colourIndex() ) remcur = false;
  if ( !remcur && d.next()->next() &&
       d.colourIndex() ==  d.next()->next()->colourIndex() ) remcur = true;

  // Do the actual removal set new momenta and mark affected parton as touched.
  FSGluonEmitter::removeGluon(d, g, (remcur? d.iPart(): d.oPart()));
  objects.insert(g); // Don't completely forget the gluon
 // We have to keep track of the original partons even if the
 // corresponding gluon was absorbed.
  (remcur? d.iPart(): d.oPart())->adoptOriginals(g);

  d.iPart()->setMomentum(p1);
  d.oPart()->setMomentum(p3);
  d.touch();
  d.iPart()->touch();
  d.oPart()->touch();
  if ( d.next() ) d.next()->touch();
  if ( d.prev() ) d.prev()->touch();
}

tParPtr DipoleState::create(tcPDPtr pd, tcParPtr parent, bool hfs) {
  tParPtr p = create<Parton>();
  p->data(pd);
  if ( parent ) p->origSystem(parent->system());
  if ( hfs ) addHadronicFS(p);
  return p;
}

void DipoleState::forgetParton(tParPtr p) {
  theHadronicFinalState.erase(p);
  theHardFinalState.erase(p);
  theFinalState.erase(p);
  objects.erase(p);
}

void DipoleState::forgetDipole(tDBPtr d) {
  active.erase(d);
  objects.erase(d);
}


DipoleStatePtr DipoleState::fullclone() const {
  TranslationMap trans;
  DipoleStatePtr copy = preclone(trans);
  copy->postclone(trans);
  return copy;
}

DipoleStatePtr DipoleState::preclone(TranslationMap & trans) const {
  CloneSet tocopy(objects.begin(), objects.end());
  tocopy.insert(tcDipoleStatePtr(this));
  CloneSet tocheck = tocopy;
  while ( ! tocheck.empty() ) {
    CloneSet additional;
    for ( CloneSet::const_iterator it = tocheck.begin();
	  it != tocheck.end(); ++it )
      if ( *it ) (**it).fillReferences(additional);
    tocheck.clear();
    for ( CloneSet::const_iterator it = additional.begin();
	  it != additional.end(); ++it )
      if ( tocopy.insert(*it).second )
	tocheck.insert(*it);
  }

  tocopy.erase(cDipoleStatePtr(this));
  DipoleStatePtr copy = dynamic_ptr_cast<DipoleStatePtr>(clone());
  trans[tcDipoleStatePtr(this)] = copy;
  for ( CloneSet::const_iterator it = tocopy.begin();
	it != tocopy.end(); ++it ) if ( *it ) trans[*it] = (**it).clone();
  return copy;
}

void DipoleState::postclone(const TranslationMap & trans) const {
  for ( TranslationMap::const_iterator it = trans.map().begin();
	it != trans.map().end(); ++it ) it->second->rebind(trans);
}

void DipoleState::fillReferences(CloneSet & cs) const {
  CascadeBase::fillReferences(cs);
  cs.insert(theSelected);
  for ( set<EmissionGenerator>::iterator it = generators.begin();
	it != generators.end(); ++it )
    if ( it->emission ) cs.insert(it->emission);
}

void DipoleState::rebind(const TranslationMap & trans) {
  BaseSet old;
  old.swap(objects);
  objects.clear();
  trans.translate(inserter(objects), old.begin(), old.end());
  set<tParPtr> oldfs;
  oldfs.swap(theFinalState);
  theFinalState.clear();
  trans.translate(inserter(theFinalState), oldfs.begin(), oldfs.end());
  oldfs.swap(theHardFinalState);
  theHardFinalState.clear();
  trans.translate(inserter(theHardFinalState), oldfs.begin(), oldfs.end());
  oldfs.swap(theHadronicFinalState);
  theHadronicFinalState.clear();
  trans.translate(inserter(theHadronicFinalState), oldfs.begin(), oldfs.end());
  vector<tRemParPtr> orem;
  orem.swap(theRemnants.first);
  theRemnants.first.clear();
  trans.translate(inserter(theRemnants.first), orem.begin(), orem.end());
  orem.swap(theRemnants.second);
  theRemnants.second.clear();
  trans.translate(inserter(theRemnants.second), orem.begin(), orem.end());
  vector<tResPtr> ores;
  ores.swap(theResonances);
  theResonances.clear();
  trans.translate(inserter(theResonances), ores.begin(), ores.end());
  set<tDBPtr> odip;
  odip.swap(active);
  active.clear();
  trans.translate(inserter(active), odip.begin(), odip.end());
  theSelected = trans.translate(theSelected);
  set<EmissionGenerator> gold;
  gold.swap(generators);
  for ( set<EmissionGenerator>::iterator it = gold.begin();
   	it != gold.end(); ++ it ) {
    EmissionGenerator gen(trans.translate(it->dipole));
    gen.emission = trans.translate(it->emission);
    generators.insert(gen);
  }
  dipindx.clear();
  parindx.clear();
    
}

void DipoleState::persistentOutput(PersistentOStream & os) const {
  os << subprocess << theIncoming << objects << theFinalState << theRemnants
     << theResonances
     << active << ounit(theTotalMomentum, GeV) << ounit(theHardMomentum, GeV)
     << ounit(theHadronicMomentum, GeV) << generators.size()
     << theNe << theSelected;
  for ( set<EmissionGenerator>::const_iterator it = generators.begin();
	it != generators.end(); ++it ) os << *it;
}

void DipoleState::persistentInput(PersistentIStream & is, int) {
  int gsize;
  is >> subprocess >> theIncoming >> objects >> theFinalState >> theRemnants
     >> theResonances
     >> active >> iunit(theTotalMomentum, GeV) >> iunit(theHardMomentum, GeV)
     >> iunit(theHadronicMomentum, GeV) >> gsize >>
    theNe >> theSelected;
  generators.clear();
  EmissionGenerator eg;
  while ( gsize-- ) {
    is >> eg;
    generators.insert(eg);
  }
}

void DipoleState::touchHadronicState() {
  for ( int i = 0, N = theRemnants.first.size(); i < N; ++i )
    if ( hadronicFS().find(theRemnants.first[i]) == hadronicFS().end() )
      theRemnants.first[i]->touch();
  for ( int i = 0, N = theRemnants.second.size(); i < N; ++i )
    if ( hadronicFS().find(theRemnants.second[i]) == hadronicFS().end() )
      theRemnants.second[i]->touch();
}

void DipoleState::untouchHadronicState() {
  for ( int i = 0, N = theRemnants.first.size(); i < N; ++i )
    if ( hadronicFS().find(theRemnants.first[i]) == hadronicFS().end() )
      theRemnants.first[i]->untouch();
  for ( int i = 0, N = theRemnants.second.size(); i < N; ++i )
    if ( hadronicFS().find(theRemnants.second[i]) == hadronicFS().end() )
      theRemnants.second[i]->untouch();
}

void DipoleState::transformHadronicState(const LorentzRotation & R) {
  UtilityBase::transform(theHadronicFinalState, R);
  touchHadronicState();
}


// The following static variable is needed for the type description
// system in ThePEG.
DescribeClass<DipoleState,CascadeBase>
  describeAriadne5DipoleState("Ariadne5::DipoleState", "libAriadne5.so");

void DipoleState::Init() {}

pair<tcQCDPtr,tcQCDPtr> DipoleState::StringEnds(tcQCDPtr d) {
  pair<tcQCDPtr,tcQCDPtr> ret = make_pair(d, d);
  while ( tcQCDPtr dn = ret.second->next() ) {
    if ( dn == d ) return ret;
    ret.second = dn;
  }

  while ( tcQCDPtr dp = ret.first->prev() ) {
    if ( dp == d ) return ret;
    ret.first = dp;
  }

  return ret;

}
    

void DipoleState::debugme() const {
  CascadeBase::debugme();
  cerr << "DipoleState:" << endl;
  LorentzMomentum sum;
  set<tParPtr> allpartons = finalState();
  set<tDBPtr> tdipoles = activeDipoles();
  set<tcDBPtr> dipoles(tdipoles.begin(), tdipoles.end());
  cerr << "> active dipoles:" << endl;
  while ( !dipoles.empty() ) {
    if ( tcQCDPtr d = dynamic_ptr_cast<tcQCDPtr>(*dipoles.begin()) ) {
      dipoles.erase(d);
      pair<tcQCDPtr,tcQCDPtr> range = StringEnds(d);
      tParPtr p1 = range.first->iPart();
      tParPtr p2 = range.second->oPart();
      if ( p1 == p2 )
	cerr << "  string loop:" << endl;
      else
	cerr << "  string: " << index(p1)
	     << " ... " << index(p2) << endl;
      while ( range.first ) {
	range.first->iPart()->debug();
	cerr << endl;
	sum += range.first->iPart()->momentum();
	allpartons.erase(range.first->iPart());
	range.first->debug();
	cerr << endl;
	dipoles.erase(range.first);
	if ( range.first == range.second )
	  range.first = tcQCDPtr();
	else
	  range.first = range.first->next();
      }
      if ( p1 != p2 ) {
	p2->debug();
	cerr << endl;
	sum += p2->momentum();
	allpartons.erase(p2);
      }
    } else {
      cerr << "  other dipole:" << endl;
      (**dipoles.begin()).debug();
      cerr << endl;
      dipoles.erase(dipoles.begin());
    }
  }

  cerr << "> rest of final state:" << endl;
  while ( !allpartons.empty() ) {
    (**allpartons.begin()).debug();
    cerr << endl;
    sum += (**allpartons.begin()).momentum();
    allpartons.erase(allpartons.begin());
  }

  cerr << "> sum of momenta:" << setprecision(3)
       << setw(9) << ( abs(sum.x()) > MeV? sum.x(): 0.0*GeV )/GeV
       << setw(9) << ( abs(sum.y()) > MeV? sum.y(): 0.0*GeV )/GeV
       << setw(9) << ( abs(sum.z()) > MeV? sum.z(): 0.0*GeV )/GeV
       << setw(9) << sum.e()/GeV
       << setw(9) << sum.m()/GeV << endl;

  if ( selected() ) {
    cerr << "> selected emission:" << endl;
    selected()->debug();
    cerr << "selected dipole: " << index(selected()->cdipole);
    cerr << endl;
  }
  
}

void DipoleState::debugEmissions() const {
  for ( unsigned i = 0; i < emissions.size(); ++i )
    if ( emissions[i] ) emissions[i]->debug();
}

bool DipoleState::checkIntegrity() {
  for ( set<tDBPtr>::iterator it = active.begin();
      it != active.end(); ++it ) {
    if( ! (*it)->checkIntegrity() ) {
      return false;
    }
  }
  return true;
}

int DipoleState::index(tcCascadeBasePtr o) const {
  dipindx(tcDBPtr());
  parindx(tcParPtr());
  if ( dynamic_ptr_cast<tcDBPtr>(o) )
    return dipindx(dynamic_ptr_cast<tcDBPtr>(o));
  if ( dynamic_ptr_cast<tcParPtr>(o) )
    return parindx(dynamic_ptr_cast<tcParPtr>(o));
  return 0;
}

SaveDipoleState::SaveDipoleState(tDipoleStatePtr state, bool force)
  : forced(force) {
  if ( force ||
       !( state->selected()->failsafe || state->selected()->reversible ) )
    backup = state->preclone(trans);
  else
    backup = state;
}

tDipoleStatePtr SaveDipoleState::revert() {
  static DebugItem printsteps("Ariadne5::PrintSteps", 6);
  if ( printsteps )
    cerr << "Ariadne5::PrintStep: Reverting emission " << backup->theNe + 1
	 << " of type " << typeid(*(backup->selected())).name() << endl;
  if ( forced ||
       !( backup->selected()->failsafe || backup->selected()->reversible ) )
    backup->postclone(trans);
  else if ( backup->selected()->failsafe )
    Throw<Exception>()
      << "The emitter model " << backup->selected()->model->fullName()
      << " has reported that an emission is failsafe, but the emission "
      << "failed anyway. Please contact the author to have this error "
      << "corrected." << Exception::runerror;
  else if ( backup->selected()->reversible ) {
  // Save four momentum before reverting.
    LorentzMomentum ptot = backup->totalMomentum();
    backup->selected()->revert();
  Energy2 s = ptot.m2();
  ptot -= backup->sumTotalMomentum();
  if( (ptot.vect().mag2() + sqr(ptot.e()) )/s > 1.0e-16) 
    throw DipoleState::MomentumException()
      << "Ariadne5::DipoleState::perform: "
      << "Emission did not conserve four momentum."
      << Exception::eventerror;
  }
  backup->untouch();
  return backup;
}

