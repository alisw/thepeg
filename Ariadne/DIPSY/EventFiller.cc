// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EventFiller class.
//

#include "DipoleEventHandler.h"
#include "EventFiller.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Utilities/Current.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/EventRecord/Step.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/SubProcess.h"
#include "RealPartonState.h"
#include "ParticleInfo.h"
#include "ThePEG/Utilities/Debug.h"
#include "ThePEG/Utilities/DebugItem.h"
#include "ThePEG/Utilities/UtilityBase.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Handlers/LuminosityFunction.h"

#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Reference.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "../Cascade/EmitterBase.h"

#include <iostream>
#include <fstream>

using namespace DIPSY;

EventFiller::EventFiller()
  : currentWeight(0.0), theRecoilScheme(0), theSingleMother(0),
    theDGLAPinPT(0), theEffectiveWeights(0), theFSSwingTime(0.0),
    theFSSwingTimeStep(0.1),
    theValenceChargeNormalisation(0), thePTCut(ZERO), theSoftRemove(1) {}
EventFiller::~EventFiller() {}

IBPtr EventFiller::clone() const {
  return new_ptr(*this);
}

IBPtr EventFiller::fullclone() const {
  return new_ptr(*this);
}

double EventFiller::fill(Step & step, DipoleEventHandler & eh, tPPair inc,
			 DipoleState & dl, DipoleState & dr,
			 const ImpactParameters & b) const {

  // Get a list of possible dipole-dipole interactions
  FList fl = eh.xSecFn().flist(dl, dr, b);

  //Sum over the (ununitarised) interaction probabilities 2*f_{ij}.
  //Note that the factor 2 is included in flist(), to give the correct
  //ND interaction probability.
  double sum = 0.0;
  for ( FList::iterator it = fl.begin(); it != fl.end(); ++it )
    sum += it->first.first;

  //Just check that they are not all 0.
  if ( sum == 0.0 ) {
    return 0.0;
  }

  //Non-diffractive interaction probability is 1 - exp(-Sum(2*f_{ij}))
  double weight = eh.xSecFn().unitarize(sum);

  //Combine the weights from all the sources.
  currentWeight =
    weight*sqr(hbarc)*dr.weight()*dl.weight()*b.weight()/eh.maxXSec();

  // Select the interactions which should be performed.
  pair<RealPartonStatePtr, RealPartonStatePtr> realStates = 
    selectInteractions(fl, b, eh.xSecFn());

  //If no interactions found, discard event.
  if ( realStates.first->interactions.empty() ) {
    return 0.0;
  }

  //Counts the number of participants in a HI event.
  //TODO: interface, or trigger on AA.
  countParticipants(dl, dr, b.bVec().pt());

  //Figure out which partons to keep, and which to remove.
  //Also sort out colour flow, momenta etc.
  vector<String> strings = extractStrings(dl, dr, realStates, b);
  DipoleState & finalState = dl;

  //Discard event if no final state partons were found.
  if ( strings.empty() || ! fillStep(step, inc, strings) ) {
    weight = 0.0;
    currentWeight = 0.0;
    return weight;
  }

  // The last thing we do is to fix up valens configurations.
  finalState.fixValence(step);

  return weight;
}

void EventFiller::countParticipants(const DipoleState & dl,
				    const DipoleState & dr, const InvEnergy b) const {
  //The initial dipoles. Note that these will in general have emitted
  //something, and thus no longer be active. They may have been
  //reconnected by a new dipole in the absorption process though.
  vector<DipolePtr> valenceDip = dl.initialDipoles();
  valenceDip.insert(valenceDip.end(),
		    dr.initialDipoles().begin(), dr.initialDipoles().end());

  //count the nonparticipating nucleons
  int untouchedNucleons = 0;
  for ( int i = 0; i < int(valenceDip.size()); i++ ) {
    //First find the three valence partons in the nucleon that
    //the dipole belongs to.
    //Two are connected to the dipole, easy.
    tPartonPtr p1 = valenceDip[i]->partons().first;
    tPartonPtr p2 = valenceDip[i]->partons().second;
    //The three dipoles in valenceDip are always after each other,
    //so the third parton can be found in the dipole before or
    //after (or both).
    tPartonPtr p3;
    if ( i != int(valenceDip.size())-1 &&
	 valenceDip[i+1]->partons().first == p2  )
      p3 = valenceDip[i+1]->partons().second;
    else if ( i != int(valenceDip.size())-1 &&
	      valenceDip[i+1]->partons().second == p1  )
      p3 = valenceDip[i+1]->partons().first;
    else if ( i != 0 && valenceDip[i-1]->partons().first == p2  )
      p3 = valenceDip[i-1]->partons().second;
    else if ( i != 0 && valenceDip[i-1]->partons().second == p1  )
      p3 = valenceDip[i-1]->partons().first;

    //If none of the valence partons have interacted
    //ie (have children that interacted), the nucleon is considered
    //to not have interacted. Note that with swings, a single interaction
    //can make several nucleons participating.
    if ( !(p1->interacted()) &&
	 !(p2->interacted()) &&
	 !( p3 && p3->interacted() )  ) {
      valenceDip[i]->participating(false);
      untouchedNucleons++;
    }
    else
      valenceDip[i]->participating(true);
  }

  //We will triple count, as we loop through all 3 dipoles in
  //each nucleon. Divide by 3 to make up for this.
  untouchedNucleons /= 3;
}

/*
 * This is the probability that a parton is interacting, given that
 * the state interacts, but the partons is not selected as primary.
 *
 * Prob of p interacting: p
 * Prob of p interacting, given that the state interacts: x = p/totalP
 * Prob of p being selected as primary, given state interacts:  y = p/sumP
 * Prob of p being selected as non-primary, given state
 * interacts: z
 * x = y + z ==> z = p/totalP - p/sumP 
 * Prob of p not being primary, given state interacts: A = (sumP-p)/sumP
 * And now what we are looking for:
 * Prob of p being selected as non-primary, given state interacts
 * and not selected as primary: B
 * z = A*B ==> B = z/A = p(1/totalP - 1/sumP)/((sumP-p)/sumP) =
 * = p(sumP/totalP - 1)/(sumP-p)
 */
double correctedProb(double totalP, double sumP, double p) {
  return p*(sumP/totalP - 1.0)/(sumP - p);
}

pair<RealPartonStatePtr, RealPartonStatePtr>
EventFiller::selectInteractions(const FList & fl, const ImpactParameters & b,
				const DipoleXSec & xSec) const {
  DipolePairVector interactions;

  double sumfij = 0.0; //the sum of the individual nonuntarised int probs
  double sumUP = 0.0;   //The sum of the individual unitarised int probs

  //Set up a selector, and a map sorted on pt of the interaction.
  Selector<FList::const_iterator, double> sel;
  DipolePairMap ordered;

  //Loop through all possible interactions and fill
  //@sel and @ordered.
  for ( FList::const_iterator it = fl.begin(); it != fl.end(); ++it ) {

    //If no interaction probability, don't bother.
    if ( it->first.second == 0.0 ) continue;

    //insert into the selector.
    sel.insert(it->first.second, it);

    sumfij += it->first.first;
    sumUP += it->first.second;

    //Calculate interaction recoils by calling the DipoleXSec object.
    DipoleXSec::InteractionRecoil rec =
      xSec.recoil(it->second.first->partons(), it->second.second->partons(), b);

    //Insert in @ordered, sorting on max pt recoil.
    double pt = max(max(rec.first.first.pt(), rec.first.second.pt()),
		    max(rec.second.first.pt(), rec.second.second.pt()))/GeV;
    ordered.insert(make_pair(pt, it));
  }

  //interaction probability (for at least one interaction)
  //of the two states.
  double totalP = Current<DipoleEventHandler>()->xSecFn().unitarize(sumfij);

  //create the real states.
  RealPartonStatePtr rrs = new_ptr(RealPartonState());
  RealPartonStatePtr lrs = new_ptr(RealPartonState());

  //Add the valence partons (as they are always real) and save.
  lrs->addValence(fl.begin()->second.first->dipoleState());
  rrs->addValence(fl.begin()->second.second->dipoleState());
  lrs->saveState();
  rrs->saveState();

  DipolePairMap potential;
  DipolePairMap failedPrims;
  bool found = false;
  int counter = 0;
  while ( !found ) {
    potential.clear();

    //select a first interaction (since there has to be at least one).
    FList::const_iterator prim = sel[UseRandom::rnd()];

    //Go through the other interactions and check if they are also
    //interacting. A modified probability is used, to make up for the
    //bias introduced in selecting a primary interaction.
    for ( FList::const_iterator it = fl.begin(); it != fl.end(); ++it )
      if ( it == prim || correctedProb(totalP,sumUP,it->first.second) 
	                   > UseRandom::rnd() ) {
	DipoleXSec::InteractionRecoil rec =
      xSec.recoil(it->second.first->partons(), it->second.second->partons(), b);
      double pt = max(max(rec.first.first.pt(), rec.first.second.pt()),
		      max(rec.second.first.pt(), rec.second.second.pt()))/GeV;

      //If the interaction passed the amplitude test, add it to the list
      //of potential interations.
	potential.insert(make_pair(pt, it));
      }

    //Keep track on how many times we have tried to find a consistent
    //set of interactions. Give up after 10 tries.
    counter++;
    if ( counter > 10 ) {
      overTenEvents += currentWeight;
      return make_pair(lrs, rrs);
    }

    //test all potetential interactions in order, skipping already failed dips.
    //failed first interactions inserted in failedprims,
    //and from sel and ordered if first in ordered.
    int i = 0;
    for ( DipolePairMap::iterator it = potential.begin(); it != potential.end(); ++it ) {
      i++;

      //Check first interaction, and already failed as first.
      //Then autofail without checking again.
      if ( failedPrims.find(it->first) != failedPrims.end() && !found ) {
	continue;
      }

      //Try to add the interaction.
      if (addInteraction(it->second, lrs, rrs, interactions, b, xSec)) {
	found = true;
      }
      else {
	//remember if it was first interaction and failed, to not
	//try it again later.
	if ( !found ) {
	  failedPrims.insert(*it);
	  if ( it == ordered.begin() ) {
	    ordered.erase(it);
	    sel.erase(it->second);
	  }
	}
      }
    }

    //if sel (or ordered) empty, give up event. no int can be first.
    if ( sel.empty() ) return make_pair(lrs, rrs);
  }

  //Make sure that the real state partons are on-shell.
  //May not actually be needed, not sure... TODO: check!
  for( DipolePairVector::iterator it = interactions.begin();
      it != interactions.end(); it++ ) {
    lrs->setOnShell((*it)->second.first);
    rrs->setOnShell((*it)->second.second);
  }
  for ( RealParton::RealPartonSet::iterator it = lrs->valence.begin();
	it != lrs->valence.end(); it++ )
    (*it)->setOnShell();
  for ( RealParton::RealPartonSet::iterator it = rrs->valence.begin();
	it != rrs->valence.end(); it++ )
    (*it)->setOnShell();

  //return the real states.
  return make_pair(lrs, rrs);
}



bool EventFiller::
addInteraction(FList::const_iterator inter, RealPartonStatePtr lrs,
	       RealPartonStatePtr rrs, DipolePairVector & inters,
	       const ImpactParameters & b, const DipoleXSec & xSec) const {
    rescatter += currentWeight;

  //With some settings, only some of the partons actually interact.
  //Check which here. Other tunes will just return 4x true.
  pair<pair<bool, bool>, pair<bool, bool> > doesInt =
    xSec.doesInt(inter->second.first->partons(),
		 inter->second.second->partons(), b);

  //Calculate the recoils from the interaction.
  DipoleXSec::InteractionRecoil recs =
    xSec.recoil(inter->second.first->partons(),
		inter->second.second->partons(), b, doesInt);

  //Add the interacting partons and their parents to the real state
  //and go through the evolution checking for ordering, local pt-max, etc.
  if ( !lrs->fullControlEvolution
       (inter->second.first, inter->second.second,
	doesInt.first.first, doesInt.first.second,
	recs.first.first.pt(), recs.first.second.pt()) ) {
    return false;
  }

  //Add the interaction to the other state, and check that as well.
  if ( !rrs->fullControlEvolution
       (inter->second.second, inter->second.first,
	doesInt.second.first, doesInt.second.second,
	recs.second.first.pt(), recs.second.second.pt()) ) {
    //If this second evolution failed, remove the interaction from
    //the first
    //state by reverting it to the previously saved state.
    lrs->revertToPrevious(inter->second.first);

    return false;
  }

  //If both evolutions could accomodate the interacting partons,
  //add the interaction recoil as well, 
  //and check that things still are ok.
  inters.push_back(inter);
  if ( !controlRecoils(inters, lrs, rrs, b, xSec, doesInt) ) {
    //If failed, remove the interaction from the evolutions.
    lrs->revertToPrevious(inter->second.first);
    rrs->revertToPrevious(inter->second.second);
    inters.pop_back();

    return false;
  }

  //If the interactions was ok, mark the dipoles as interacting
  //and save the real states.
  inter->second.first->interact(*inter->second.second);
  inter->second.second->interact(*inter->second.first);
  lrs->saveState();
  rrs->saveState();

  return true;
}

vector<EventFiller::String>
EventFiller::extractStrings(DipoleState & dl, DipoleState & dr,
			    pair<RealPartonStatePtr, RealPartonStatePtr> realStates, 
                            const ImpactParameters & b ) const {

  //Just rename a bit for convenience.
  DipoleStatePtr rightState = &dr;
  DipoleStatePtr leftState = &dl;
  RealPartonStatePtr lrs = realStates.first;
  RealPartonStatePtr rrs = realStates.second;


  //grow back the partons marked as virtuals.
  //This will be the ones without interacting children, but not the
  //ones removed during ordering/pt-max-fixing.
  removeVirtuals(leftState);
  removeVirtuals(rightState);

  //The unordered and not-ok pt-max were removed by pairing them
  //up with other partons that should stay. Now merge all the momentum
  //into a single parton, and flag the others as virtual.
  lrs->mergeVirtuals();
  rrs->mergeVirtuals();
  //Save to update the parton information from the realPartons.
  lrs->saveState();
  rrs->saveState();

  //A high-pt parton (so very localised in x_T) can not recoil all
  //of a low-pt (so smeared out in x_T) parent, but will rather shoot
  //out a recoiling part of the parent. Fix this by replacing recoiled
  //low-pt parents by a remnant, and a recoiler.
  lrs->addRecoilers();
  rrs->addRecoilers();
  //Save states to update partons.
  lrs->saveState();
  rrs->saveState();

  //balance p+/p- by moving the entire state. Was first done in the 
  //interaction (where the interacting partons had to provide the 
  //energy), but the modifications here changes kinematics so that
  //some more balancing may be needed. Should very rarely be
  //large changes.
  fixBoost(lrs, rrs);

  //Translate the right state to it's right place in x_T,
  //and mirror it in rapidity so that it actually comes from the
  //other side. Finally merge the two cascades into one state.
  rightState->translate(b);
  rightState->mirror(0.0);
  DipoleStatePtr finalState = leftState->merge(rightState);

  vector<pair<DipolePtr, DipolePtr> > swinging =
    Current<DipoleEventHandler>()->xSecFn().getColourExchanges(lrs, rrs);

  //swing the interacting dipoles.
  for ( int i = 0; i < int(swinging.size()); i++ ) {
    Current<DipoleEventHandler>()->
      xSecFn().reconnect(swinging[i].first, swinging[i].second);
  }

  //now remove the partons that were made virtual by mergeVirtuals.
  //This is so that it should be easier to reconnect the colour
  //through the interacting dipoles with the state from the other side.
  removeVirtuals(finalState);

 //do final state swings, aka pythias colour reconnection.
  double yrange = FSSwingTime();
  double ystep = FSSwingTimeStep();
  for ( double y = 0.0; y < yrange; y += ystep ) {
    finalState->swingFS(y, y + ystep);
  }

  //  finalState->lambdaMeasure(0.36*GeV2, histDipLength2, histDipMass2);

  //compensate for the 6 valence charges in the proton
  finalState->normaliseValenceCharge(theValenceChargeNormalisation);

  //make sure the non-particpating nucleons (if AA) have the
  //original colour flow.
  finalState->restoreNonparticipants();

  //try to find and fix any problems or inconsistencies that may have
  //popped up. Error messages are sent if anything is found.
  dodgeErrors(finalState);

  //If you want to save the initial state event to file (mainly AA).
  //TODO: interface!
  //TODO: stop event after this.
  //finalState->saveGluonsToFile(currentWeight);

  //If you want to make a fancy movie in mathematica of your event. :)
  // finalState->printForMovie(0, 1000);

  //Sort the remaining partons in strings.
  vector<String> ret = finalState->strings();

  return ret;
}

bool EventFiller::fillStep (Step & step, tPPair incoming,
			    const vector<String> & strings) const {
  Direction<0> dir(true);
  vector< vector<PPtr> > particles(strings.size());
  SubProPtr sub = new_ptr(SubProcess(incoming));
  for ( int is = 0, NS = strings.size(); is < NS; ++is ) {
    int N = strings[is].size();
    vector<bool> removed(N, false);
    if ( theSoftRemove && pT2Cut() > ZERO && N > 2 ) {
      bool changed = true;
      while ( changed ) {
	changed = false;
	for ( int i = 0; i < N; ++i ) {
	  if ( removed[i] ) continue;
	  if ( strings[is][i]->valence() && theSoftRemove == 2 ) continue;
	  if ( strings[is][i]->flavour() != ParticleID::g ) continue;
	  int i1 = (N + i - 1)%N;
	  while ( removed[i1] ) i1 = (N + i1 - 1)%N;
	  int i3 = (i + 1)%N;
	  while ( removed[i3] ) i3 = (i3 + 1)%N;
	  if ( i1 == i3 ) continue;
	  if ( invPT2(strings[is], i1, i, i3) < pT2Cut() ) {
	    if ( removeGluon(strings[is], i1, i, i3) ) {
	      removed[i] = changed = true;
	    } else {
	      return false;
	    }
	  }
	}
      }
    }	
    particles[is].resize(N);
    for ( int i = 0; i < N; ++i ) {
      if ( removed[i] ) continue;
      particles[is][i] = strings[is][i]->produceParticle();
      int ip = i - 1;
      while ( ip >= 0 && removed[ip] ) --ip;
      if ( ip >= 0 )
	particles[is][i]->colourConnect(particles[is][ip]);
      if ( i == N - 1 && strings[is][i]->flavour() == ParticleID::g ) {
	int i0 = 0;
	while ( i0 < i && removed[i0] ) ++i0;
	particles[is][i0]->colourConnect(particles[is][i]);
      }
      sub->addOutgoing(particles[is][i]);
    }
    if ( strings[is][0]->flavour() == ParticleID::g ) {
      int i0 = 0;
      int in = N - 1;
      while ( i0 < in && removed[i0] ) ++i0;
      while ( in > i0 && removed[in] ) --in;
      particles[is][i0]->colourConnect(particles[is][in]);
    }
  }

  step.addSubProcess(sub);

  return true;
}

void EventFiller::fixBoost(RealPartonStatePtr lrs, RealPartonStatePtr rrs) const {
  //access the CoM energy, through a kindof roundabout path...
  //TODO: must be a better way! Current?
  PartonPtr dummy = lrs->partons.begin()->first;
  Energy sqrtS = ((dummy->dipoles().first) ? dummy->dipoles().first:dummy->dipoles().second)->dipoleState().handler().lumiFn().maximumCMEnergy();

  //sum left and right p+ and p-.
  Energy leftPlus = ZERO;
  Energy leftMinus = ZERO;
  Energy rightPlus = ZERO;
  Energy rightMinus = ZERO;

  //loop and sum over the to-keep partons.
  for ( map<tPartonPtr,RealPartonPtr>::const_iterator
	  it =lrs->partons.begin();
	it != lrs->partons.end(); it++ ) {
    if ( it->second->keep == RealParton::NO ) continue;
    leftPlus += it->second->plus;
    leftMinus += it->second->minus;
  }
  for ( map<tPartonPtr,RealPartonPtr>::const_iterator
	  it = rrs->partons.begin();
	  it != rrs->partons.end(); it++ ) {
    if ( it->second->keep == RealParton::NO ) continue;
    rightPlus += it->second->minus;
    rightMinus += it->second->plus;
  }

  //Solve the 2:nd degree equation for how much energy has to be
  //transfered to set both states on shell.
  double A = (- rightPlus*rightMinus - sqr(sqrtS) + leftPlus*leftMinus)/(2.0*rightMinus*sqrtS);
  double B = rightPlus/rightMinus;
  double y1 = -1.0;
  double x1 = -1.0;
  if ( sqr(A) - B > 0.0 ) {
    //the factor to change right p- with.
    y1 = - A + sqrt(sqr(A) - B);
    // double y2 = - A - sqrt(sqr(A) - B);

    //The factor to change left p+ with.
    x1 = (y1*sqrtS - rightPlus)/(y1*leftPlus);
    // double x2 = (y2*sqrtS - rightPlus)/(y2*leftPlus);
  }

  //error handling
  if ( x1 < 0 || y1 < 0 ) {
    Throw<SpaceLikeGluons>()
      << "EventFiller::fixBoost gave negative or nan solution to boost equation, "
      << "will not balance momentum.\n"
      << "This is probably caused by a rouge gluon being far too virtual.\n"
      << Exception::warning;
    return;
  }

  //loop again, scaling the p+ and p- according to the solution above.
  for (   map<tPartonPtr,RealPartonPtr>::const_iterator it = lrs->partons.begin();
	  it != lrs->partons.end(); it++) {
    if ( it->second->keep == RealParton::NO ) continue;
    it->second->plus *= x1;
    it->second->minus /= x1;
  }
  for (   map<tPartonPtr,RealPartonPtr>::const_iterator it = rrs->partons.begin();
	  it != rrs->partons.end(); it++) {
    if ( it->second->keep == RealParton::NO ) continue;
    it->second->minus /= y1;
    it->second->plus *= y1;
  }

  //Save state to transfer the changes to the partons.
  lrs->saveState();
  rrs->saveState();
}

// void EventFiller::fixValence(Step & step, DipoleState & dl, DipoleState & dr) const {
//   list<PartonPtr> alll = dl.getPartons();
//   set<tcPartonPtr> vall;
//   for ( list<PartonPtr>::iterator pit = alll.begin(); pit != alll.end(); ++pit )
//     if ( (**pit).valence() ) vall.insert(*pit);
//   list<PartonPtr> allr = dr.getPartons();
//   set<tcPartonPtr> valr;
//   for ( list<PartonPtr>::iterator pit = allr.begin(); pit != allr.end(); ++pit )
//     if ( (**pit).valence() ) valr.insert(*pit);
//   vector<PPtr> lval;
//   vector<PPtr> rval;
//   for ( ParticleSet::iterator pit = step.particles().begin();
// 	pit != step.particles().end(); ++pit ) {
//     tcPartonPtr p = ParticleInfo::getParton(**pit);
//     if ( !p ) continue;
//     if ( member(vall, p) ) lval.push_back(*pit);
//     else if ( member(valr, p) ) rval.push_back(*pit);
//   }

//   if ( UseRandom::rndbool() ) {
//     dl.fixValence(step, lval);
//     dr.fixValence(step, rval);
//   } else {
//     dr.fixValence(step, rval);
//     dl.fixValence(step, lval);
//   }

// }

void EventFiller::dodgeErrors(DipoleStatePtr finalState) const {
  list<PartonPtr> partons = finalState->getPartons();
  for ( list<PartonPtr>::iterator it = partons.begin(); it != partons.end(); it++ ) {
    PartonPtr p = *it;
    //check for empty pointers
    if ( !p ) {
      Throw<SpaceLikeGluons>()
	<< "dodgeError found empty pointer from getPartons()! :o" << Exception::warning;
      continue;
    }
    //check for colour flow to itself
    //check for 0 monetum
    if ( p->pT().pt() == ZERO ) {
      Throw<SpaceLikeGluons>()
	<< "dodgeError found 0 pt gluon." << Exception::warning;
      p->pT(TransverseMomentum(UseRandom::rnd()*GeV,UseRandom::rnd()*GeV));
    }
    if ( p->plus() == ZERO || p->minus() == ZERO ) {
      Throw<SpaceLikeGluons>()
	<< "dodgeError found 0 lightcone momentum." << Exception::warning;
      p->minus(UseRandom::rnd()*GeV);
      p->plus(UseRandom::rnd()*GeV);
    }

    //check for nan momentum
    if ( isnan(p->pT().pt()/GeV) ) {
      Throw<SpaceLikeGluons>()
	<< "dodgeError found NAN pt gluon, fix pt." << Exception::warning;
      p->pT(TransverseMomentum(UseRandom::rnd()*GeV,UseRandom::rnd()*GeV));
    }
    if ( isnan(p->plus()/GeV) || isnan(p->minus()/GeV) ) {
      Throw<SpaceLikeGluons>()
	<< "dodgeError found NAN lightcone momentum, fix plus, minus." << Exception::warning;
      p->minus(UseRandom::rnd()*GeV);
      p->plus(UseRandom::rnd()*GeV);
    }

    //check for negative momentum
    if ( p->pT().pt() < ZERO ) {
      Throw<SpaceLikeGluons>()
	<< "dodgeError found negative pt gluon.... >_>, fix pt." << Exception::warning;
      p->pT(TransverseMomentum(UseRandom::rnd()*GeV,UseRandom::rnd()*GeV));
    }
    if ( p->plus() < ZERO || p->minus() < ZERO ) {
      Throw<SpaceLikeGluons>()
	<< "dodgeError found negative lightcone momentum, fix plus,minus." << Exception::warning;
      p->minus(UseRandom::rnd()*GeV);
      p->plus(UseRandom::rnd()*GeV);

      //check for off-shell momentum
      if ( sqr(p->plus()*p->minus() - p->pT().pt2() - sqr(p->mass())) > sqr(sqr(0.01*GeV)) ) {
	Throw<SpaceLikeGluons>()
	  << "dodgeError found off-shell parton, fix pt." << Exception::warning;
	if ( p->plus()*p->minus() < sqr(p->mass()) ) {
	  Throw<SpaceLikeGluons>()
	    << "dodgeError found insufficient energy for mass, fix plus,minus." << Exception::warning;
	  double ratio = 2.0*sqr(p->mass())/(p->plus()*p->minus()); //give a bit extra for pt
	  p->plus(p->plus()*sqrt(ratio));
	  p->minus(p->minus()*sqrt(ratio));
	}
	double mod = sqrt(p->plus()*p->minus() - sqr(p->mass()))/p->pT().pt();
	p->pT(p->pT()*mod);
      }
    }
  }
}

bool EventFiller::
controlRecoils(DipolePairVector & sel,
	       RealPartonStatePtr lrs, RealPartonStatePtr rrs,
	       const ImpactParameters & b, const DipoleXSec & xSec,
	       pair<pair<bool, bool>, pair<bool, bool> > doesInt) const {

  //Keep track on which partons are interacting.
  list<pair<bool, bool> >::const_iterator
    leftDoesInt = lrs->doesInts.begin();
  list<pair<bool, bool> >::const_iterator
    rightDoesInt = rrs->doesInts.begin();

  //Loop through the interaction, and
  for ( DipolePairVector::const_iterator inter = sel.begin();
	inter != sel.end(); inter++ ) {

    //Which partons are interacting in this dip-dip interaction
    pair<pair<bool, bool>, pair<bool, bool> >
      trueDoesInt = make_pair(*leftDoesInt, *rightDoesInt);
    //calculate recoils
    DipoleXSec::InteractionRecoil recoil =
      xSec.recoil((*inter)->second.first->partons(),
		  (*inter)->second.second->partons(), b, trueDoesInt);

    //call the DipoleXSec object to check kinematics and vetos
    //for the interaction. If things don't work, undo the changes.
    if ( !xSec.doInteraction(recoil, *inter, lrs, rrs, trueDoesInt, b) ) {
      //revert the total recoil between the states.
      lrs->totalRecoil -= recoil.first.first + recoil.first.second;
      rrs->totalRecoil -= recoil.second.first + recoil.second.second;
    }
  }
  return true;
}

void EventFiller::removeVirtuals(DipoleStatePtr state) const {

  list<PartonPtr> vP = state->getPartons();

  //loop through the partons, and remove the off shell ones.
  //this only does colour flow. RealPartonState does the kinematics
  //in mergeVirtuals.
  for ( list<PartonPtr>::iterator it = vP.begin(); it != vP.end(); it++) {
    tPartonPtr p = *it;
    //only handle off-shell parton.
    if ( !(p->onShell()) ) {
      //if there are only 2 or fewer on-shell partons in the colour chain,
      //it has to be swinged at some point, and better early than late,
      //as last-second forced swings can lead to silly colour connections.
      if ( p->nOnShellInChain() < 2 ) {
	//tell the absorber to find a swing anywhere in the colour chain
	if ( p->dipoles().first ) absorber()->swingLoop(p->dipoles().first, *state);
	else absorber()->swingLoop(p->dipoles().second, *state);
      }
      
      //once off-shell loops are handled, absorb the parton.
      absorber()->removeParton(p);
    }
  }
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void EventFiller::doinit() throw(InitException) {
  HandlerBase::doinit();
}

void EventFiller::dofinish() {}

void EventFiller::doinitrun() {
  HandlerBase::doinitrun();
}

double EventFiller::pTScale(DipoleState & state) const {
  return state.handler().emitter().pTScale();
}

Energy2 EventFiller::invPT2(const String & str, int i1, int i2, int i3) {
  LorentzMomentum p1 = str[i1]->momentum();
  LorentzMomentum p2 = str[i2]->momentum();
  LorentzMomentum p3 = str[i3]->momentum();
  Energy2 s = (p1 + p2 + p3).m2();
  if ( s < sqr(str[i1]->mass() + str[i3]->mass() ) )
    Throw<SpaceLikeGluons>()
      << "DIPSY produced space-like gluons. Three neighboring ones had "
      << "negative mass squared. This cannot be fixed at the moment. "
      << "Event discarded." << Exception::eventerror;
  Energy2 s12 = (p1 + p2).m2();
  Energy2 s23 = (p2 + p3).m2();
  if ( s12 < ZERO || s23 < ZERO ) return -1.0*GeV2;
  return s12*s23/s;
}

bool EventFiller::removeGluon(const String & str, int i1, int i2, int i3) {
  LorentzMomentum p1 = str[i1]->momentum();
  LorentzMomentum p2 = str[i2]->momentum();
  LorentzMomentum p3 = str[i3]->momentum();
  LorentzMomentum ptest = p1 + p2;
  if ( ptest.m2() < Constants::epsilon*1000.0*sqr(ptest.e()) ) {
    str[i1]->plus(ptest.plus());
    str[i1]->pT(TransverseMomentum(ptest.x(), ptest.y()));
    return true;
  }
  ptest = p3 + p2;
  if ( ptest.m2() < Constants::epsilon*1000.0*sqr(ptest.e()) ) {
    str[i3]->plus(ptest.plus());
    str[i3]->pT(TransverseMomentum(ptest.x(), ptest.y()));
    return true;
  }

  Energy2 S = (p1 + p2 + p3).m2();
  if ( S <= ZERO ) return false;
  LorentzRotation R = Utilities::boostToCM(makeTriplet(&p1, &p2, &p3));
  double Psi = Constants::pi - p3.theta();
  double beta = 0.0;
  Energy W = sqrt(S);
  double x1 = 2.0*p1.e()/W;
  double x3 = 2.0*p3.e()/W;
  bool g1 = str[i1]->flavour() == ParticleID::g;
  bool g3 = str[i3]->flavour() == ParticleID::g;
  if ( ( g1 && g3 ) || (!g1 && !g3 ) )
    beta = Psi*sqr(x3)/(sqr(x1) + sqr(x3)); // minimize pt
  else if ( g1 )
    beta = Psi;
  R.rotateY(-beta);
  R.invert();
  Lorentz5Momentum p1n(p1.m());
  Lorentz5Momentum p3n(p3.m());
  try {
    SimplePhaseSpace::CMS(p1n, p3n, S, 1.0, 0.0);
  } catch ( ImpossibleKinematics ) {
    return false;
  }
  p1n.transform(R);
  p3n.transform(R);
  str[i1]->plus(p1n.plus());
  str[i1]->pT(TransverseMomentum(p1n.x(), p1n.y()));
  str[i3]->plus(p3n.plus());
  str[i3]->pT(TransverseMomentum(p3n.x(), p3n.y()));
  return true;
}




void EventFiller::persistentOutput(PersistentOStream & os) const {
  os << theAbsorber << theRecoilScheme << theMode << theSingleMother
     << theDGLAPinPT << theEffectiveWeights << theFSSwingTime
     << theFSSwingTimeStep << theValenceChargeNormalisation
     << ounit(thePTCut, GeV) << theSoftRemove;
}

void EventFiller::persistentInput(PersistentIStream & is, int) {
  is >> theAbsorber >> theRecoilScheme >> theMode >> theSingleMother
     >> theDGLAPinPT >> theEffectiveWeights >> theFSSwingTime
     >> theFSSwingTimeStep >> theValenceChargeNormalisation
     >> iunit(thePTCut, GeV) >> theSoftRemove;
}

DescribeClass<EventFiller,HandlerBase>
describeDIPSYEventFiller("DIPSY::EventFiller", "libAriadne5.so libDIPSY.so");

void EventFiller::Init() {

  static ClassDocumentation<EventFiller> documentation
    ("The EventFiller class is able to produce an initial ThePEG::Step "
     "from two colliding DipoleStates.");


  static Reference<EventFiller,DipoleAbsorber> interfaceDipoleAbsorber
    ("DipoleAbsorber",
     "The object used to absorb non-interacting dipoles.",
     &EventFiller::theAbsorber, true, false, true, true, false);

  static Switch<EventFiller,int> interfaceMode
    ("Mode",
     "How the real state is found from the virtual cascade. Speed versus consistency.",
     &EventFiller::theMode, 0, true, false);
  static SwitchOption interfaceModeMode
    (interfaceMode,
     "Consistent",
     "The fully consistent version. Not currently reccomended.",
     0);
  static SwitchOption interfaceModeFast
    (interfaceMode,
     "Fast",
     "Checking only the new real partons introduced by each new interaction, rather than rechecking them all. Good for heavy ions. Not currently reccomended.",
     1);
  static SwitchOption interfaceModeSingle
    (interfaceMode,
     "SingleSweep",
     "Checking all partons, but sweeping just once in each direction, making it a lot faster. Good for heavy ions. May give an occasional unordered chain, but hopefully not too often. The recommended option.",
     2);
  static SwitchOption interfaceModeNonRecursive
    (interfaceMode,
     "NonRecursive",
     "Does all evo no matter what. Then removes biggest problem, and does it all over again. Never turns partons back on once switched off.",
     3);

  static Switch<EventFiller,int> interfaceEffectiveWeights
    ("EffectiveWeights",
     "How the p+ and pT recoils are distributed among the single partons in an effective parton.",
     &EventFiller::theEffectiveWeights, 0, true, false);
  static SwitchOption interfaceEffectiveWeightsPlusWeighted
    (interfaceEffectiveWeights,
     "PlusWeighted",
     "Weight pt and plus according to p+ of the individual partons.",
     0);
  static SwitchOption interfaceEffectiveWeightsPlusEvenWeighted
    (interfaceEffectiveWeights,
     "PlusEvenWeighted",
     "The plus is distibuted according to the plus of the partons, but the pt is shared evenly among the partons.",
     1);
  static SwitchOption interfaceEffectiveWeightsPlusSingleWeighted
    (interfaceEffectiveWeights,
     "PlusSingleWeighted",
     "The plus is distibuted according to the plus of the partons, but the pt is taken only by the colour connected parton",
     2);

  static Switch<EventFiller,int> interfaceRecoilScheme
    ("RecoilScheme",
     "How to give tread the positive light-cone momentum of a recoiling parton.",
     &EventFiller::theRecoilScheme, 0, true, false);
  static SwitchOption interfaceRecoilSchemePreservePlus
    (interfaceRecoilScheme,
     "PreservePlus",
     "blaha",
     0);
  static SwitchOption interfaceRecoilSchemeFixedY
    (interfaceRecoilScheme,
     "FixedY",
     "blahaaa",
     1);
  static SwitchOption interfaceRecoilSchemeFrameFan
    (interfaceRecoilScheme,
     "FrameFan",
     "dssklajhls",
     2);

  static Parameter<EventFiller,int> interfaceSingleMother
    ("SingleMother",
     "If an emission is regarded to come from a single parton rather than both "
     " partons in the emitting dipole.",
     &EventFiller::theSingleMother, 1, 1, 0, 0,
     true, false, Interface::lowerlim);

  static Parameter<EventFiller,int> interfaceDGLAPinPT
    ("DGLAPinPT",
     "If the DGLAP supression should be made in terms of pt rather than r. "
     " partons in the emitting dipole.",
     &EventFiller::theDGLAPinPT, 1, 1, 0, 0,
     true, false, Interface::lowerlim);


  static Parameter<EventFiller,double> interfaceFSSwingTimeStep
    ("FSSwingTimeStep",
     "How long time steps is are to be used for FS colour reconnections.",
     &EventFiller::theFSSwingTimeStep, 0.1, 0.0, 0,
     true, false, Interface::lowerlim);


  static Parameter<EventFiller,double> interfaceFSSwingTime
    ("FSSwingTime",
     "How long time is allowed for FS colour reconnections. 0 turns it off.",
     &EventFiller::theFSSwingTime, 0.0, 0.0, 0,
     true, false, Interface::lowerlim);

  static Switch<EventFiller,int> interfaceValenceChargeNormalisation
    ("ValenceChargeNormalisation",
     "How to treat the (too large) colour charge of the valence partons.",
     &EventFiller::theValenceChargeNormalisation, 0, true, false);
  static SwitchOption interfaceValenceChargeNormalisationNone
    (interfaceValenceChargeNormalisation,
     "None",
     "Dont do anything.",
     0);
  static SwitchOption interfaceValenceChargeNormalisationSwing
    (interfaceValenceChargeNormalisation,
     "Swing",
     "Swing some of the dipoles going to the valence partons",
     1);

  static Parameter<EventFiller,Energy> interfacePTCut
    ("PTCut",
     "The minimum invariant transverse momentum allowed for a gluon. "
     "Gluons below the cut will be removed from the final state. "
     "If zero, no gluons will be removed.",
     &EventFiller::thePTCut, GeV, 0.0*GeV, 0.0*GeV, 0*GeV,
     true, false, Interface::lowerlim);

  static Switch<EventFiller,int> interfaceSoftRemove
    ("SoftRemove",
     "Determines if gluons with invariant transverse momentum below "
     "<interface>PTCut</interface> should be reabsorbed.",
     &EventFiller::theSoftRemove, 1, true, false);
  static SwitchOption interfaceSoftRemoveOff
    (interfaceSoftRemove,
     "Off",
     "No gluons are absorbed",
     0);
  static SwitchOption interfaceSoftRemoveAll
    (interfaceSoftRemove,
     "All",
     "All soft gluons below the cut are absorbed.",
     1);
  static SwitchOption interfaceSoftRemoveNoValence
    (interfaceSoftRemove,
     "NoValence",
     "All except valence gluons are absorbed if below the cut.",
     2);

}





