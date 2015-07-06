// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RealPartonState class.
//

#include "RealPartonState.h"
#include "RealParton.h"
#include "Dipole.h"
#include "DipoleState.h"
#include "DipoleEventHandler.h"
#include "ThePEG/Utilities/Current.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/Debug.h"
#include "ThePEG/Utilities/DebugItem.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Repository/CurrentGenerator.h"

#include <iostream>

using namespace DIPSY;

RealPartonState::RealPartonState(): consistent(true), calls(0) {}

RealPartonState::~RealPartonState() {}

void RealPartonState::newInteraction(tDipolePtr intDip, tDipolePtr otherIntDip,
				     bool firstInt, bool secondInt,
				     Energy rec1, Energy rec2) {
  //The ancestors of the interacting partons, ie the partons that have
  //been involved in making the interacting dipole. 
  RealPartonSet ancestors;

  //diagnostics
  if ( Debug::level > 5 )
    cout << "adding interaction " << intDip << " to " << this << endl;

  //clear the partons to be checked in this interaction.
  //Only partons involved in this interaction are checked,
  //ther others are assumed to still be ok since last check.
  toCheck.clear();

  //Add interacting dipoles and which partons are interacting.
  interactions.push_back(intDip);
  doesInts.push_back(make_pair(firstInt, secondInt));

  //Set up the real partons.
  //Some settings only add one of partons in the interacting dipole,
  //hence all the "if (firstInt )" etc.
  //Note that addParton recursively calls the parents as well,
  //so all new realparton instances needed are created.
  RealPartonPtr rp1, rp2;
  if ( firstInt ) {
    if ( Debug::level > 5 ) cout << "calling first addparton for interacting parton " << intDip->partons().first->oY() << endl;
    addParton(intDip->partons().first);
    rp1 = getReal(intDip->partons().first);
  }
  if ( secondInt ) {
    if ( Debug::level > 5 ) cout << "calling second addparton for interacting parton " << intDip->partons().second->oY() << endl;
    addParton(intDip->partons().second);
    rp2 = getReal(intDip->partons().second);
  }

  //Find the interaction distance of the interacting partons.
  //This is needed when checking pT-ordering (or almost equivalently dipole-
  //size ordering) later, and the last-before-interaction dipole needs
  //to be checked. The picture here is to compare along the colour chains
  //AFTER the interacting dipole has swinged, which is a more frame-
  //independent approach, as otherwise we would have to compare the
  //last-before-interacting dipole with a dipole that is not present in
  //the final state.
  //Also, since the interacting dipole has (roughly) the weight of a
  //dipole that has emitted, it would mess up the weights if it entered
  //the pT-max reweighting scheme.
  if ( rp1 ) {
    if ( !rp1->interacting ) {
      if ( rec1 == ZERO ) rec1 = rec2;
      if ( rec1 == ZERO ) cerr << "zero recoil in realpartonstate::newinteraction" << endl;
      rp1->intDist = rp1->theParton->pTScale()/rec1;
    }
    else {
      if ( rec1 == ZERO ) rec1 = rec2;
      if ( rec1 == ZERO ) cerr << "zero recoil in realpartonstate::newinteraction" << endl;
      rp1->intDist = min(rp1->intDist, rp1->theParton->pTScale()/rec1);
    }

    //book keeping.
    rp1->interacting = intDip;
    rp1->secondInt = true;
    RealPartonSet a1 = rp1->ancestors();
    ancestors.insert(a1.begin(), a1.end());
  }
  if ( rp2 ) {
    if ( !rp2->interacting ) {
      if ( rec2 == ZERO ) rec2 = rec1;
      if ( rec2 == ZERO ) cerr << "zero recoil in realpartonstate::newinteraction" << endl;
      rp2->intDist = rp2->theParton->pTScale()/rec2;
    }
    else {
      if ( rec2 == ZERO ) rec2 = rec1;
      if ( rec2 == ZERO ) cerr << "zero recoil in realpartonstate::newinteraction" << endl;
      rp2->intDist = min(rp2->intDist, rp2->theParton->pTScale()/rec2);
    }

    //book keeping for future referens
    rp2->interacting = intDip;
    rp2->firstInt = true;
    RealPartonSet a2 = rp2->ancestors();
    ancestors.insert(a2.begin(), a2.end());
  }

  //and make all the ancestors remember that they contributed.
  for ( RealPartonSet::iterator it = ancestors.begin(); it != ancestors.end(); it++ ) {
    (*it)->interactions.insert(intDip);
  }
}

tRealPartonPtr RealPartonState::getReal(tPartonPtr p) {
  //If the parton already exists in the real state, just retur
  //the real parton.
  map<tPartonPtr,RealPartonPtr>::iterator it = partons.lower_bound(p);
  if ( it->first == p ) return it->second;

  //Otherwise create a new realparton for the parton.
  RealPartonPtr rp = new_ptr(RealParton(p));
  rp->realState = this;

  //fluct -1 means that they are not a fluctuation,
  //ie they are not (yet) set up to be merged with an other partons.
  rp->fluct = -1;
  rp->cFluct = -1;

  //insert in the right place in the ordered map of parton -> realparton.
  if ( it != partons.begin() ) partons.insert(--it, make_pair(p, rp));
  else partons.insert(make_pair(p, rp));
  return rp;
}
  

tRealPartonPtr RealPartonState::addParton(tPartonPtr p) {
  if ( !p ) return RealPartonPtr();
  if ( partons.lower_bound(p)->first == p ) return getReal(p);
  RealPartonPtr rp = getReal(p);
  rp->nMothers = 2;
  if ( !(p->valence()) && rp->interactions.size() == 0 ) {
    if ( Current<DipoleEventHandler>()->eventFiller().singleMother() &&
	 !p->swingedEmission() ) {
      rp->nMothers = 1;
      double P1 = p->dist2(*p->parents().second)/
	(p->dist2(*p->parents().second) + p->dist2(*p->parents().first));
      if ( P1 > UseRandom::rnd() ) {//pick feynman diagram
	rp->setOMother(addParton(p->parents().first));
      }
      else
	rp->setOMother(addParton(p->parents().second));
    }
    else {
      rp->setFirstOMother(addParton(p->parents().first));
      rp->setSecondOMother(addParton(p->parents().second));
    }
  }
  if ( p->valence() && rp->interactions.size() == 0 ) {
    valence.insert(rp);
    rp->setYES();
  }
  else toCheck.insert(rp);
  return rp;
}

void RealPartonState::addValence(DipoleState & state) {
  if ( Debug::level > 5 ) cout << "entering addvalence" << endl;
  plus = state.plus();
  minus = state.minus();
  if ( Debug::level > 5 ) cout << "  total plus and minus set to " << plus/GeV << ", " << minus/GeV << endl;
  monitored = false;
  totalRecoil = TransverseMomentum(ZERO, ZERO);
  if ( Debug::level > 5 ) cout << "  total recoil set to 0" << endl;
  vector<DipolePtr> valenceDip = state.initialDipoles();
  for ( int i = 0; i < int(valenceDip.size()); i++) {
    tRealPartonPtr p1 = addParton(valenceDip[i]->partons().first);
    tRealPartonPtr p2 = addParton(valenceDip[i]->partons().second);
    p1->cValence     = true;
    p1->cValencePlus = p1->theParton->valencePlus();
    p1->cValencePT   = p1->theParton->valencePT();
    p2->cValence     = true;
    p2->cValencePlus = p2->theParton->valencePlus();
    p2->cValencePT   = p2->theParton->valencePT();
    oValence.insert(p1);
    oValence.insert(p2);
  }
}

void RealPartonState::merge(RealPartonPtr rp1, RealPartonPtr rp2) {
      if ( Debug::level > 5 )
	cout << "merging partons at y = " << rp1->oY() << ", " << rp2->oY() << " or ("
	     << rp1->theParton->position().x()*GeV << ", " << rp1->theParton->position().y()*GeV << "), ("
	     << rp2->theParton->position().x()*GeV << ", " << rp2->theParton->position().y()*GeV << ")\n";
  if ( rp1->fluct != -1 && rp2->fluct != -1 ) {
    //merge flucts if different (can they already be in same?)
    if ( rp1->fluct != rp2->fluct ) {
      flucts[rp1->fluct].insert(flucts[rp2->fluct].begin(), flucts[rp2->fluct].end());
      int old = rp2->fluct;
      for ( RealPartonSet::iterator it = flucts[old].begin(); it != flucts[old].end(); it++ )
	(*it)->fluct = rp1->fluct;
      flucts[old].clear();
    }
  }
  else if ( rp1->fluct != -1 ) {
    rp2->fluct = rp1->fluct;
    flucts[rp1->fluct].insert(rp2);
  }
  else if ( rp2->fluct != -1 ) {
    rp1->fluct = rp2->fluct;
    flucts[rp2->fluct].insert(rp1);
  }
  else {
    flucts.push_back(RealPartonSet());
    int i = flucts.size() - 1;
    flucts[i].insert(rp1);
    flucts[i].insert(rp2);
    rp1->fluct = i;
    rp2->fluct = i;
  }

  // redistributePlus(flucts[rp1->fluct]);
  //dela upp p+ och pT jamnt mellan de virtuella
  makeCollinear(flucts[rp1->fluct]);
  if ( Debug::level > 5 ) cout << "  done merging" << endl;
}

void RealPartonState::makeCollinear(const RealPartonSet & rps) {
  if ( rps.empty() ) return;
  Energy avPlus = ZERO;
  TransverseMomentum avPT = TransverseMomentum();
  for ( RealPartonSet::const_iterator it = rps.begin(); it != rps.end(); it++ ) {
    avPlus += (*it)->plus;
    avPT += (*it)->pT;
  }
  avPlus /= rps.size();
  avPT /= rps.size();
  for ( RealPartonSet::const_iterator it = rps.begin();
	it != rps.end(); it++ ) {
    (*it)->plus = avPlus;
    (*it)->pT = avPT;
    (*it)->updateYMinus();
  }
}

void RealPartonState::redistributePlus(const RealPartonSet & rps) {
  if ( rps.empty() ) return;
  Energy totalPlus = ZERO;
  double totalWeight = ZERO;
  for ( RealPartonSet::const_iterator it = rps.begin();
	it != rps.end(); it++ ) {
    totalWeight += exp(-(*it)->oY());
    totalPlus += (*it)->plus;
  }
  for ( RealPartonSet::const_iterator it = rps.begin();
	it != rps.end(); it++ ) {
    (*it)->plus = totalPlus*exp(-(*it)->oY())/totalWeight;
    (*it)->updateYMinus();
    // cout << "    redistributed plus at " << (*it)->oY() << " to " << (*it)->plus/GeV
    // 	 << ", relative y is " << -log(exp(-(*it)->oY())/totalWeight) << endl;
  }
}

void RealPartonState::splitFluct(RealPartonPtr rp1, RealPartonPtr rp2) {
  // cout << "entered splitFluct with rps " << rp1->oY() << ", " << rp2->oY() << "!\n";
  // diagnosis(true);
  if ( rp1->fluct == -1 || rp2->fluct == -1 ) {
    // cout << "asked to split a fluct where the partons werent in a fluct" << endl;
    return;
  }
  else if ( rp1->fluct != rp2->fluct ) {
    // cout << "asked to split a fluct that is already split" << endl;
    return;
  }
  int i = rp1->fluct;
  //create two new flucts with rp1 and rp2
  flucts.push_back(RealPartonSet());
  int i1 = flucts.size() - 1;
  flucts[i1].insert(rp1);
  flucts.push_back(RealPartonSet());
  int i2 = flucts.size() - 1;
  flucts[i2].insert(rp2);
  // cout << "    created new flucts " << i1 << " and " << i2 << ", old was " << i << endl;
  //move all other partons to the fluct belonging to the closest parton.
  for ( RealPartonSet::iterator it = flucts[i].begin(); it != flucts[i].end(); it++) {
    if ((*it)->theParton->dist2(*rp1->theParton) < (*it)->theParton->dist2(*rp2->theParton) ) {
      // cout << "   move " << (*it)->oY() << " to " << i1 << endl;
      (*it)->fluct = i1;
      flucts[i1].insert(*it);
    }
    else {
      // cout << "   move " << (*it)->oY() << " to " << i2 << endl;
      (*it)->fluct = i2;
      flucts[i2].insert(*it);
    }
  }
  // cout << "   empty old fluct" << endl;
  flucts[i].clear();
  
  //find highest oY parton.
  RealPartonPtr firstRP = (rp1->oY() < rp2->oY() ) ? rp1:rp2;
  RealPartonPtr secondRP = (rp1->oY() < rp2->oY() ) ? rp2:rp1;
  //let the higher oY recoil to the lower if they are connected. pT only, no plus transfer.
  if ( secondRP->mothers.first == firstRP ) {
    // cout << firstRP->oY() << " is first mother of " << secondRP << endl;
    secondRP->pT += secondRP->theParton->recoil(firstRP->theParton);
    secondRP->updateYMinus();
    InvEnergy range = sqrt(min(secondRP->theParton->dist2(*firstRP->theParton),
			       firstRP->theParton->dist2(*secondRP->mothers.second->theParton)/4.0));
    secondRP->doEffectiveRecoil(firstRP, range, true, ZERO,
				secondRP->theParton->recoil(firstRP->theParton));
  }
  else if ( secondRP->mothers.second == firstRP ) {
    secondRP->pT += secondRP->theParton->recoil(firstRP->theParton);
    secondRP->updateYMinus();
    InvEnergy range = sqrt(min(secondRP->theParton->dist2(*firstRP->theParton),
			       firstRP->theParton->dist2(*secondRP->mothers.first->theParton)/4.0));
    secondRP->doEffectiveRecoil(firstRP, range, false, ZERO,
				secondRP->theParton->recoil(firstRP->theParton));
  }
  makeCollinear(flucts[i1]);
  makeCollinear(flucts[i2]);
  // diagnosis(true);
  //given yminus should be ok I hope?
}

bool RealPartonState::isMotherStepUp(tRealPartonPtr mom, bool firstSide, InvEnergy scale) {
  if ( mom->theParton->valence() ) return false;
  tRealPartonPtr grandMom = (firstSide ? mom->mothers.first:mom->mothers.second);
  if ( mom->nMothers == 1 ) grandMom = mom->mother;

  if ( firstSide && (!mom->children.first.empty() || !mom->firstInt) ) return false;
  if ( !firstSide && (!mom->children.second.empty() || !mom->secondInt) ) return false;
  
  if ( mom->theParton->dist2(*grandMom->theParton) > sqr(scale) ) return false;
  else return true;
}

bool RealPartonState::inFSRRegion(tRealPartonPtr rp) {
  //valence are always ok
  if ( rp->theParton->valence() ) return false;

  //if parton is outside, but not to an interaction
  if ( isOutside(rp, true) ) {
    RealPartonPtr mom = (rp->nMothers == 1 ? rp->mother:rp->mothers.first);

    //merged partons will share the momentum, and thus only have 1/N of the pt.
    int nMerged = 1;
    if ( rp->fluct != -1 ) nMerged = flucts[rp->fluct].size();

    //check ordering with outside partons without effective partons
    //dont chek if already merged, and account for that plus and minus is shared
    //between merged partons
    if ( !(mom->fluct != -1 && mom->fluct == rp->fluct) ) {

      if ( Debug::level > 5 ) cout << "checking FSR PS for " << rp->oY() << ", mom: "
			       << mom->oY() << endl;

      int nMomMerged = 1;
      if ( mom->fluct != -1 ) nMomMerged = flucts[mom->fluct].size();
      if ( mom->minus*double(nMomMerged) > rp->minus*double(nMerged) ) return true;
    }
    //dont bother to check child if rp is interacting since
    //we dont know where the other parton is yet
    if ( !rp->secondInt && !rp->children.second.empty() ) {
      RealPartonPtr child = *rp->children.second.rbegin();
      if ( !(child->fluct != -1 && child->fluct == rp->fluct) ) {

	if ( Debug::level > 5 ) cout << "checking FSR PS for " << rp->oY()
				 << ", child: " << child->oY() << endl;

	int nChildMerged = 1;
	if ( child->fluct != -1 ) nChildMerged = flucts[child->fluct].size();
	if ( child->plus*double(nChildMerged) > rp->plus*double(nMerged) ) return true;
      }
    }
  }

  //and check for outside on the other side in the same way.
  if ( isOutside(rp, false) ) {
    RealPartonPtr mom = (rp->nMothers == 1 ? rp->mother:rp->mothers.second);

    //merged partons will share the momentum, and thus only have 1/N of the pt.
    int nMerged = 1;
    if ( rp->fluct != -1 ) nMerged = flucts[rp->fluct].size();

    //check ordering with outside partons without effective partons
    //dont chek if already merged, and account for that plus and minus is shared
    //between merged partons
    if ( !(mom->fluct != -1 && mom->fluct == rp->fluct) ) {

      if ( Debug::level > 5 ) cout << "checking FSR PS for " << rp->oY() << ", mom: "
			       << mom->oY() << endl;

      int nMomMerged = 1;
      if ( mom->fluct != -1 ) nMomMerged = flucts[mom->fluct].size();
      if ( mom->minus*double(nMomMerged) > rp->minus*double(nMerged) ) return true;
    }
    //dont bother to check child if rp is interacting since
    //we dont know where the other parton is yet
    if ( !rp->firstInt && !rp->children.first.empty() ) {
      RealPartonPtr child = *rp->children.first.rbegin();
      if ( !(child->fluct != -1 && child->fluct == rp->fluct) ) {

	if ( Debug::level > 5 ) cout << "checking FSR PS for " << rp->oY()
				 << ", child: " << child->oY() << endl;

	int nChildMerged = 1;
	if ( child->fluct != -1 ) nChildMerged = flucts[child->fluct].size();
	if ( child->plus*double(nChildMerged) > rp->plus*double(nMerged) ) return true;
      }
    }
  }
  return false;
}

void RealPartonState::fixFSRDoubleCounting(tRealPartonPtr rp) {
  RealPartonPtr mom;
  if ( rp->nMothers == 1 ) mom = rp->mother;
  else if ( isOutside(rp, true) ) mom = rp->mothers.first;
  else mom = rp->mothers.second;
  if ( rp->fluct != -1 && rp->fluct == mom->fluct ) return;
  if ( Debug::level > 5 ) cout << rp->oY() << " is FSR Double Counted with respect to "
			   << mom->oY() << ". Merge." << endl;
  if ( !rp->recoils.empty() ) {
    TransverseMomentum pT;
    if ( rp->nMothers == 1 ) pT = rp->theParton->recoil(mom->theParton);
    else pT = rp->theParton->recoil(rp->mothers.first->theParton) +
	   rp->theParton->recoil(rp->mothers.second->theParton);
    if ( Debug::level > 5 ) cout << "removing pT (" << pT.x()/GeV << ", " << pT.y()/GeV << ")" << endl;

    rp->pT -= pT;
    rp->plus -= pT.pt()*exp(-rp->oY());
    rp->undoRecoils();
    if ( rp->plus < ZERO ) fixNegativePlus(rp);
  }

  merge(rp, mom);
  if ( Debug::level > 5 ) cout << "fixed FSR double counting." << endl;
}

Energy RealPartonState::fixNegativePlus(RealPartonPtr rp) {
  //if already nonnegative, nothing has to be done.
  if ( rp->plus >= ZERO ) return rp->plus;
  //valence partons cant be helped, since there are no parents.
  if ( rp->theParton->valence() ) {
    cerr << "valence parton with negative plus in RealPartonState::fixNegativePlus" << endl;
    return rp->plus;
  }

  if ( rp->nMothers == 1 ) {
    InvEnergy range = sqrt(rp->theParton->dist2(*rp->mother->theParton));

    Energy totalPlus = rp->mother->effectivePlusMinus(range, true).first + rp->plus;
    //check that there is enough plus in parents, otherwise return false.
    if ( totalPlus < ZERO ) {
      Throw<RealPartonKinematicException>()
	<< "not enough plus in single parent to balance negative plus in "
	<< "RealPartonState::fixNegativePlus" << Exception::warning;
      return rp->plus;
    }

    //if enough, divide up so that the parton takes 10% of the effective parents plus.
    rp->doEffectiveRecoil(rp->mother, range, true,
			  -rp->plus + totalPlus/10.0, TransverseMomentum());
    rp->plus = totalPlus/10.0;
    return totalPlus/10.0;
  }
  else if ( rp->nMothers == 2 ) {
    RealPartonPtr mom1 = rp->mothers.first;
    RealPartonPtr mom2 = rp->mothers.second;
    InvEnergy range1 = sqrt(min(rp->theParton->dist2(*mom1->theParton),
				mom2->theParton->dist2(*mom1->theParton)/4.0));
    InvEnergy range2 = sqrt(min(rp->theParton->dist2(*mom2->theParton),
				mom1->theParton->dist2(*mom2->theParton)/4.0));
    Energy plus1 = mom1->effectivePlusMinus(range1, true).first;
    Energy plus2 = mom2->effectivePlusMinus(range2, false).first;
    Energy parentPlus = plus1 + plus2;
    Energy totalPlus = parentPlus + rp->plus;

    //check that there is enough plus in parents, otherwise return false.
    if ( parentPlus + rp->plus < ZERO ) {
      Throw<RealPartonKinematicException>()
	<< "not enough plus in parents to balance negative plus in "
	<< "RealPartonState::fixNegativePlus" << Exception::warning;
      return rp->plus;
    }

    //if enough, give enough to put plus to 0
    rp->doEffectiveRecoil(mom1, range1, true,
			  -rp->plus*plus1/parentPlus + totalPlus*plus1/parentPlus/10.0,
			  TransverseMomentum());
    rp->doEffectiveRecoil(mom2, range2, false,
			  -rp->plus*plus2/parentPlus + totalPlus*plus2/parentPlus/10.0,
			  TransverseMomentum());
    rp->plus = totalPlus/10.0;
    return totalPlus/10.0;
  }

  return rp->plus;
}

bool RealPartonState::isOutside(tRealPartonPtr rp, bool firstSide) {
  //to be outside, there can be neither kids not an interaction on one side
  //the other side has to have either an interaction or a child. also valence dont count.
  if ( (firstSide || rp->nMothers == 1) && rp->children.first.empty() && !rp->firstInt &&
       ( rp->secondInt || !rp->children.second.empty() ) && !rp->theParton->valence() )
    return true;
  if ( (!firstSide || rp->nMothers == 1) && rp->children.second.empty() && !rp->secondInt &&
       ( rp->firstInt || !rp->children.first.empty() ) && !rp->theParton->valence() )
    return true;
  return false;
}

InvEnergy RealPartonState::childScale(tRealPartonPtr rp, bool firstSide) {
  InvEnergy ret = ZERO;
  if ( !rp->interacting ) { //if no interaction, then take last child
    if ( firstSide ? rp->children.second.empty():rp->children.first.empty() ) return ZERO;

    RealPartonPtr child = (firstSide ? *rp->children.second.rbegin():
			   *rp->children.first.rbegin());
    if ( Debug::level > 5 ) cout << "  child is " << child->oY() << endl;
    ret = sqrt(rp->theParton->dist2(*child->theParton));
    if ( Debug::level > 5 ) cout << "    child dist is " << ret*GeV << endl;
    if (rp->fluct != -1 && child->fluct == rp->fluct) {
      ret = childScale(child, firstSide);
    if ( Debug::level > 5 ) cout << "    in same fluct, recur " << endl;
    }
    else {
      InvEnergy childDist = min(child->intDist, Current<DipoleEventHandler>()->coherenceRange());
      if ( child->interacting && child->resolutionScales.lower_bound(rp)->first == rp &&
	   child->resolutionScales[rp] <
	   sqr(childDist)*Current<DipoleEventHandler>()->alphaS(childDist) ) {
      ret = child->intDist;
      if ( Debug::level > 5 ) cout << "    child is merged int, new scale: " << ret*GeV << endl;
      }
    }
  }
  else {
    ret = rp->intDist;
    if ( Debug::level > 5 ) cout << "  intdist is " << ret*GeV << endl;
  }
  return ret;
}

bool RealPartonState::checkFixSingleDGLAP(RealPartonPtr rp, InvEnergy scale) {
  if ( rp->DGLAPchecked )  return true;
  rp->DGLAPchecked = true;

  tRealPartonPtr mom = rp->mother;

  if ( Debug::level > 5 ) cout << rp->oY() << " is outside (single). mother: " << mom->oY() << endl;

  //if the parton is already tagged to merge with the mother, then dont do anything
  if ( rp->fluct != -1 && mom->fluct == rp->fluct ) {
    if ( Debug::level > 5 ) cout << "  already merged, nm." << endl;
    return false;
  }

  if ( isOutside(rp, true) ) {

    //if no scale (distance to mother and child) was supplied, then calculate them
    if ( scale == ZERO ) scale = childScale(rp, true);
    // *** ATTENTION *** not used:   InvEnergy motherScale = sqrt(rp->theParton->dist2(*(mom->theParton)));

    //check if mom is a further step up in pt, then check her first, with this scale
    if ( isMotherStepUp(mom, true, scale) || isMotherStepUp(mom, false, scale) ) {
      //if the higher backwards pt max is safe, then this is also ok.
      if ( checkFixDGLAP(mom, scale) ) {
	if ( Debug::level > 5 ) cout << mom->oY() << " is safe. then so is " << rp->oY() << endl;
	return true;
      }
    }
    //compare to mother and fix if not safe
    if ( !rp->DGLAPSafe(scale) ) {
      fixDGLAP(rp);
      return false;
    }
  }

  if ( isOutside(rp, false) ) {

    //if no scale (distance to mother and child) was supplied, then calculate them
    if ( scale == ZERO ) scale = childScale(rp, false);
    // *** ATTENTION *** not used:   InvEnergy motherScale = sqrt(rp->theParton->dist2(*(mom->theParton)));

    //check if mom is a further step up in pt, then check her first, with this scale
    if ( isMotherStepUp(mom, true, scale) || isMotherStepUp(mom, false, scale) ) {
      //if the higher backwards pt max is safe, then this is also ok.
      if ( checkFixDGLAP(mom, scale) ) {
	if ( Debug::level > 5 ) cout << mom->oY() << " is safe. then so is " << rp->oY() << endl;
	return true;
      }
    }
    //compare to mother and fix if not safe
    if ( !rp->DGLAPSafe(scale) ) {
      fixDGLAP(rp);
      return false;
    }
  }
  return true;
}

void RealPartonState::fixSingleDGLAP(RealPartonPtr rp) {
  if ( rp->fluct != -1 && rp->fluct == rp->mother->fluct )
    return; //already merged with the mother.
  else if ( rp->fluct != -1 && rp->recoils.empty() ) {
    if ( !rp->interacting )
      merge(rp, rp->mother); //no recoils that has to be undone.
    return;
  }

  TransverseMomentum pT = rp->theParton->recoil(rp->mother->theParton);
  rp->pT -= pT; //remove rec from itself
  Energy originalPlus = (pT).pt()*exp(-rp->oY());
  rp->plus -= originalPlus;
  rp->undoRecoils();

  //now readd some recoil limited by the DGLAP scale.
  InvEnergy r = ZERO;
  if ( rp->interacting ) r = rp->intDist;
  else if ( !rp->children.first.empty() )
    r = sqrt((*rp->children.first.rbegin())->theParton->dist2(*rp->theParton));
  else if ( !rp->children.second.empty() )
    r = sqrt((*rp->children.second.rbegin())->theParton->dist2(*rp->theParton));
  else cerr << "fixing single DGLAP, but no children, or interaction scales!" << endl;

  //find from the resolutionscale (alpha_s(r)*r^2) what the minimum r couldve been.
  //since alpha_s(r) has no analytic inverse, do 3 recursive steps from cohernce range.
  if ( r > Current<DipoleEventHandler>()->coherenceRange() )
    r = Current<DipoleEventHandler>()->coherenceRange();
  InvEnergy resScale = sqrt(rp->resolutionScales[rp->mother]/
			    Current<DipoleEventHandler>()->alphaS(r));
  resScale = sqrt(rp->resolutionScales[rp->mother]/
		  Current<DipoleEventHandler>()->alphaS(resScale));
  resScale = sqrt(rp->resolutionScales[rp->mother]/
		  Current<DipoleEventHandler>()->alphaS(resScale));


  double penalty = resScale/r;
  if ( !rp->interacting || true )
    penalty = 0.0;
  else if ( Debug::level > 5 ) {
    cout << "  SINGLE: penalty for " << rp->oY() << " is " << penalty << ", PT1 to " << rp->mother->oY()
  	 << " goes from " << pT.pt()/GeV << " to "
  	 << penalty*pT.pt()/GeV << endl;
  }
  Energy plusRec = (penalty*pT).pt()*exp(-rp->oY());
  rp->pT += penalty*pT;
  rp->plus += plusRec;
  rp->doEffectiveRecoil(rp->mother, sqrt(rp->theParton->dist2(*rp->mother->theParton)),
			true, plusRec, penalty*pT);


  if ( !rp->interacting || true )
    merge(rp, rp->mother);
  else if ( rp->plus < ZERO || rp->mother->plus < ZERO ) {
    merge(rp, rp->mother);
  }

  if ( (rp->plus < ZERO || rp->mother->plus < ZERO) )
      cerr << "negative p+ after fixsingleDGLAP!    :o" << endl;
}

bool RealPartonState::checkFixDGLAP(RealPartonPtr rp, InvEnergy scale) {
  if ( rp->DGLAPchecked )  return true;
  rp->DGLAPchecked = true;

  //DGLAP should be checked only for partons on the "outside" of a chain.
  //if it is
  if ( isOutside(rp, true) ) { //if kids/interaction on only second side
    tRealPartonPtr mom;
    if ( rp->nMothers == 1 ) mom = rp->mother;
    else mom = rp->mothers.first;

    if ( Debug::level > 5 ) cout << rp->oY() << " is outside 1. mother: " << mom->oY() << endl;

    //if the parton is already tagged to merge with the mother, then dont do anything
    if ( rp->fluct != -1 && mom->fluct == rp->fluct ) {
      if ( Debug::level > 5 ) cout << rp->oY() << "  already merged, nm." << endl;
      return false;
    }

    //if no scale (distance to mother and child) was supplied, then calculate them
    if ( scale == ZERO ) scale = childScale(rp, true);
    InvEnergy motherScale = sqrt(rp->theParton->dist2(*(mom->theParton)));

    if ( Debug::level > 5 ) cout << rp->oY() << "  scale: " << scale*GeV
			     << ", motherscale: " << motherScale*GeV << endl;

    //if grandmother is a further step up in pt, then go there and check first
    //then compare to THIS childs scale, to get full suppression
    if ( isMotherStepUp(mom, true, motherScale) || isMotherStepUp(mom, false, motherScale) ) {
      if ( Debug::level > 5 ) cout << "  mother is step up. recur to " << mom->oY() << endl;
      //if the higher backwards pt max is safe, then this is also ok.
      if ( checkFixDGLAP(mom, scale) ) {
	if ( Debug::level > 5 ) cout << mom->oY() << " is safe. then so is " << rp->oY() << endl;
	return true;
      }
      //if the backwards max is removed, then this is the new max, and should be checked
    }

    if ( Debug::level > 5 ) cout << "  checkDGLAPsafe." << endl;

    //check if the recoil should be removed or not
    if ( !rp->checkDGLAPSafe(mom, rp, scale) ) {
      fixDGLAP(rp);
      return false;
    }
  }

  //the same thing, if kids/interaction on only first side
  if ( isOutside(rp, false) ) {
    tRealPartonPtr mom;
    if ( rp->nMothers == 1 ) mom = rp->mother;
    else mom = rp->mothers.second;

    if ( Debug::level > 5 ) cout << rp->oY() << " is outside 2. mother: " << mom->oY() << endl;

    if ( rp->fluct != -1 && mom->fluct == rp->fluct ) {
      if ( Debug::level > 5 ) cout << rp->oY() << "  already merged, nm." << endl;
      return false;
    }

    if ( scale == ZERO ) scale = childScale(rp, false);
    InvEnergy motherScale = sqrt(rp->theParton->dist2(*(mom->theParton)));

    if ( Debug::level > 5 ) cout << rp->oY() << "  scale: " << scale*GeV
			     << ", motherscale: " << motherScale*GeV << endl;

    if ( isMotherStepUp(mom, true, motherScale) || isMotherStepUp(mom, false, motherScale) ) {
      if ( Debug::level > 5 ) cout << "  mother is step up. recur to " << mom->oY() << endl;
      if ( checkFixDGLAP(mom, scale) ) {
	if ( Debug::level > 5 ) cout << mom->oY() << " is safe. then so is " << rp->oY() << endl;
	return true;
      }
    }

    if ( Debug::level > 5 ) cout << "  checkDGLAPsafe." << endl;

    if ( !rp->checkDGLAPSafe(mom, rp, scale) ) {
      fixDGLAP(rp);
      return false;
    }
  }
  if ( Debug::level > 5 ) cout << rp->oY() << " is not outside" << endl;
  return true;
}

void RealPartonState::fixDGLAP(RealPartonPtr rp) {
  if ( Debug::level > 5 && monitored ) cout << rp->oY() << " is unordered. Remove." << endl;
  if ( rp->nMothers == 1 ) {
    fixSingleDGLAP(rp);
    return;
  }
  RealPartonPtr mother = ((rp->children.first.empty() && !rp->firstInt) ? rp->mothers.first:rp->mothers.second);
  if ( rp->fluct != -1 && flucts[rp->fluct].find(mother) != flucts[rp->fluct].end() ) {
    return; //already merged with the mother.
  }
  else if ( rp->fluct != -1 && rp->recoils.empty() ) {
    if ( !rp->interacting )
      merge(rp, mother); //recoils already undone, dont reundo them.
    return;
  }
  if ( (rp->plus < ZERO || mother->plus < ZERO) && Debug::level > 5 )
    cerr << "negative p+ before fixdglap" << endl;

  //remove recoil between itself and other mother, but save recs from future emissions.
  //recoil with merging mother doesnt matter, since averaged anyways.
  RealPartonPtr otherMother = (mother == rp->mothers.first) ?
    rp->mothers.second:rp->mothers.first;
  TransverseMomentum pT1 = rp->theParton->recoil(mother->theParton);
  TransverseMomentum pT2 = rp->theParton->recoil(otherMother->theParton);
  rp->pT -= pT1 + pT2; //remove rec from itself
  rp->plus -= (pT1 + pT2).pt()*exp(-rp->oY());
  rp->undoRecoils();

  //now readd some recoil limited by the DGLAP scale.
  InvEnergy r2;
  if ( rp->interacting ) r2 = rp->intDist;
  else {
    if ( mother == rp->mothers.first )
      r2 = sqrt((*rp->children.second.rbegin())->theParton->dist2(*rp->theParton));
    else
      r2 = sqrt((*rp->children.first.rbegin())->theParton->dist2(*rp->theParton));
  }
  if ( r2 > Current<DipoleEventHandler>()->coherenceRange() )
    r2 = Current<DipoleEventHandler>()->coherenceRange();
  InvEnergy resScale = sqrt(rp->resolutionScales[mother]/
			    Current<DipoleEventHandler>()->alphaS(r2));
  resScale = sqrt(rp->resolutionScales[mother]/
		  Current<DipoleEventHandler>()->alphaS(resScale));
  resScale = sqrt(rp->resolutionScales[mother]/
		  Current<DipoleEventHandler>()->alphaS(resScale));
  double penalty = resScale/r2;
  if ( !rp->interacting )
    penalty = 0.0;
  else if ( Debug::level > 5 ) {
    cout << "  penalty for " << rp->oY() << " is " << penalty << ", PT1 to " << mother->oY()
  	 << " goes from " << pT1.pt()/GeV << " to "
  	 << penalty*pT1.pt()/GeV << endl;
  }
  double P1 = rp->theParton->dist2(*otherMother->theParton)/
    (rp->theParton->dist2(*mother->theParton) + rp->theParton->dist2(*otherMother->theParton));
  double P2 = 1.0 - P1;
  Energy plusRec = (penalty*pT1 + pT2).pt()*exp(-rp->oY());
  rp->pT += penalty*pT1 + pT2;
  rp->plus += plusRec;
  if ( Debug::level > 5 ) {
    cout << "  readding pT of (" << (penalty*pT1 + pT2).x()/GeV << ", "
  	 << (penalty*pT1 + pT2).y()/GeV << ")" << endl;
    cout << "  pT is (" << rp->pT.x()/GeV << ", " << rp->pT.y()/GeV << endl;
    cout << "second plus recoil is " << P2*plusRec/GeV << endl;
  }
  rp->doEffectiveRecoil(mother, sqrt(min(rp->theParton->dist2(*mother->theParton), mother->theParton->dist2(*otherMother->theParton)/4.0)), mother == rp->mothers.first, P1*plusRec, penalty*pT1);
  rp->doEffectiveRecoil(otherMother, sqrt(min(rp->theParton->dist2(*otherMother->theParton), mother->theParton->dist2(*otherMother->theParton)/4.0)), otherMother == rp->mothers.first, P2*plusRec, pT2);

  if ( rp->plus + mother->plus < ZERO ) { //handle this case. should be kindof rare.
    // if ( Debug::level > 5 ) cout << "returned too much plus in fixDGLAP. taking some back again... :/" << endl;
    // if ( Debug::level > 5 ) cout << " plus is " << rp->plus/GeV << endl;
    rp->doEffectiveRecoil(otherMother, sqrt(min(rp->theParton->dist2(*otherMother->theParton), mother->theParton->dist2(*otherMother->theParton)/4.0)), otherMother == rp->mothers.first, -rp->plus, TransverseMomentum());
    rp->plus = ZERO;
  }

  if ( !rp->interacting || true )
    merge(rp, mother);
  else if ( rp->plus < ZERO || mother->plus < ZERO ) {
    merge(rp, mother);
    if ( (rp->plus < ZERO || mother->plus < ZERO) )
      Throw<RealPartonKinematicException>()
	<< "even after refix negative p+ after fixdglap with interacting!" << Exception::warning;
  }
  if ( otherMother->plus < ZERO ) {
    merge(rp, mother);
    merge(rp, otherMother);
  }
  if ( otherMother->plus < ZERO )
      Throw<RealPartonKinematicException>()
	<< "other mother in RealPartonState::fixDGLAP() has negative p+!" << Exception::warning;
}

void RealPartonState::fixUnOrdered(RealPartonPtr rp, bool forced) {
  if ( !rp ) {
    cerr << "fixUnOrdered called without a realparton!" << endl;
  }

  if ( rp->theParton->valence() ) return;
  RealPartonPtr closestMother;
  if ( rp->nMothers == 1 )  closestMother = rp->mother;
  else
    closestMother = ( (rp->theParton->dist2(*rp->mothers.first->theParton)
		       < rp->theParton->dist2(*rp->mothers.second->theParton))
		      ? rp->mothers.first:rp->mothers.second);

  //dont merge interacting unless forced
  if ( rp->interacting && !forced && rp->nMothers == 2 ) {

    //force the unordered parents to merge
    if ( !rp->isFirstOrdered() && !rp->mothers.first->theParton->valence() )
      fixUnOrdered(rp->mothers.first, true);
    if ( !rp->isSecondOrdered() && !rp->mothers.second->theParton->valence() )
      fixUnOrdered(rp->mothers.second, true);
    return;
  }
  else if ( rp->interacting && !forced && rp->nMothers == 1 ) {
    if ( !rp->isSingleOrdered() && !rp->mother->theParton->valence() )
      fixUnOrdered(rp->mother, true);
    return;
  }

  if ( rp->fluct != -1 && rp->fluct == closestMother->fluct ) {
    if ( Debug::level > 5 && closestMother->theParton->valence() )
      cout << "WARNING, recuring back to valence in RPS::fixUnOrdered" << endl;
    fixUnOrdered(closestMother, forced);
    return;
  }

  //remove recoil between itself and other mother, but save recs from future emissions.
  //recoil with merging mother doesnt matter, since averaged anyways.
  RealPartonPtr otherMother = (closestMother == rp->mothers.first) ?
    rp->mothers.second:rp->mothers.first;
  TransverseMomentum pT1 = rp->theParton->recoil(closestMother->theParton);
  TransverseMomentum pT2;
  if ( rp->nMothers == 1 ) pT2 = TransverseMomentum();
  else pT2 = rp->theParton->recoil(otherMother->theParton);
  rp->pT -= pT1 + pT2; //remove rec from itself
  rp->plus -= (pT1 + pT2).pt()*exp(-rp->oY());
  rp->undoRecoils();

  merge(rp, closestMother);
}

bool RealPartonState::checkPlusMinus() {
  Energy totalPlus = ZERO;
  Energy totalMinus = ZERO;
  Energy onShellMinus = ZERO;
  for ( map<tPartonPtr,RealPartonPtr>::iterator it = partons.begin();
	it != partons.end(); it++ ) {
    tRealPartonPtr rp = (*it).second;
    if ( rp->keep == RealParton::NO ) continue;
    totalPlus += rp->plus;
    totalMinus += rp->minus;
    onShellMinus += rp->givenMinus;
    if ( rp->theParton->valence() )
      onShellMinus += sqr(rp->theParton->valencePT().pt())/rp->theParton->valencePlus();
  }
  if ( Debug::level > 5 ) cout << "   total plus: " << totalPlus/GeV << endl
			       << "  total minus: " << totalMinus/GeV << endl
			       << "onshell minus: " << onShellMinus/GeV << endl;
  return true;
}

bool RealPartonState::checkForNegatives() {
  bool ret = false;
  for ( map<tPartonPtr,RealPartonPtr>::iterator it = partons.begin();
	it != partons.end(); it++ ) {
    tRealPartonPtr rp = (*it).second;
    if ( rp->keep == RealParton::NO ) continue;
    if ( rp->minus < ZERO ) {
      // Throw<RealPartonKinematicException>()
      // 	<< "negative minus in real parton at oY " << rp->oY() << Exception::warning;
      ret = true;
    }
    if ( rp->plus < ZERO ) {
      // Throw<RealPartonKinematicException>()
      // 	<< "negative plus in real parton at oY " << rp->oY() << Exception::warning;
      ret = true;
    }
  }
  return ret;
}

bool RealPartonState::checkFlucts() {
  bool ok = true;
  for ( vector< RealPartonSet>::iterator it = flucts.begin(); it != flucts.end(); it++ ) {
    RealPartonSet & fluct = *it;
    if ( fluct.empty() ) continue;
    Energy plus = (*fluct.begin())->plus;
    Energy minus = (*fluct.begin())->minus;
    Energy ptx = (*fluct.begin())->pT.x();
    Energy pty = (*fluct.begin())->pT.y();
    double oy = (*fluct.begin())->oY();

    for( RealPartonSet::iterator jt = ++fluct.begin(); jt != fluct.end(); jt++) {
      RealPartonPtr rp = *jt;
      if ( abs((rp->plus - plus)/plus) > 0.01 ) {
	ok = false;
	if ( Debug::level > 5 )
	  cout << "real parton at " << rp->oY() << " has plus " << rp->plus/GeV << ", while "
	       << oy << " in the same fluct has " << plus/GeV << endl;
      }
      if ( abs((rp->minus - minus)/minus) > 0.01 ) {
	ok = false;
	if ( Debug::level > 5 )
	  cout << "real parton at " << rp->oY() << " has minus " << rp->minus/GeV << ", while "
	       << oy << " in the same fluct has " << minus/GeV << endl;
      }
      if ( abs((rp->pT.x() - ptx)/ptx) > 0.01 ) {
	ok = false;
	if ( Debug::level > 5 )
	  cout << "real parton at " << rp->oY() << " has ptx " << rp->pT.x()/GeV << ", while "
	       << oy << " in the same fluct has " << ptx/GeV << endl;
      }
      if ( abs((rp->pT.y() - pty)/pty) > 0.01 ) {
	ok = false;
	if ( Debug::level > 5 )
	  cout << "real parton at " << rp->oY() << " has pty " << rp->pT.y()/GeV << ", while "
	       << oy << " in the same fluct has " << pty/GeV << endl;
      }
    }
  }
  return ok;
}

void RealPartonState::mergeVirtuals() {
  if ( !checkFlucts() ) {
      Throw<RealPartonKinematicException>()
	<< "the fluctuations are not collinear in mergevirtuals!!" << Exception::warning;
  }
  for ( vector< RealPartonSet>::iterator it = flucts.begin(); it != flucts.end(); it++ ) {
    RealPartonSet & fluct = *it;
    if ( fluct.empty() ) continue;
    RealPartonPtr mother = *fluct.begin();
    if ( Debug::level > 5 ) cout << "merging a fluct from " << mother->oY() << " with "
			     << fluct.size() << " partons" << endl;
    for( RealPartonSet::iterator jt = ++fluct.begin(); jt != fluct.end(); jt++) {
      RealPartonPtr rp = *jt;
      if ( rp->theParton->valence() ) {
	if ( Debug::level > 5 ) cout << "  dont merge valence partons!" << endl;
	continue;
      }
      if ( Debug::level > 5 ) cout << "  removing " << rp->oY() << endl;
      if ( abs(mother->y - rp->y) > 0.00001 ) {
	cout << "  different y!!!! " << endl;
	plotState(true);
      }
      mother->plus += rp->plus;
      mother->pT += rp->pT;
      rp->theParton->onShell(false);
      rp->keep = RealParton::NO;
      for ( RealPartonSet::iterator kt = rp->children.first.begin();
	    kt != rp->children.first.end(); kt++ ) {
	if ( (*kt)->nMothers == 2 ) {
	  (*kt)->mothers.second = mother;
	  mother->children.first.insert(*kt);
	}
	if ( (*kt)->nMothers == 1 ) {
	  (*kt)->mother = mother;
	  mother->children.first.insert(*kt);
	}
      }
      for ( RealPartonSet::iterator kt = rp->children.second.begin();
	    kt != rp->children.second.end(); kt++ ) {
	if ( (*kt)->nMothers == 2 ) {
	  (*kt)->mothers.first = mother;
	  mother->children.second.insert(*kt);
	}
	if ( (*kt)->nMothers == 1 ) {
	  (*kt)->mother = mother;
	  mother->children.second.insert(*kt);
	}
      }
      if ( !rp->theParton->valence() ) {
	if ( rp->nMothers == 2 ) {
	  rp->mothers.first->children.second.erase(rp);
	  rp->mothers.second->children.first.erase(rp);
	  rp->mothers.first = RealPartonPtr();
	  rp->mothers.second = RealPartonPtr();
	}
	if ( rp->nMothers == 1 ) {
	  rp->mother->children.second.erase(rp);
	  rp->mother->children.first.erase(rp);
	  rp->mother = RealPartonPtr();
	}
      }
    }
    mother->updateYMinus();
    // plotState(true);
  }
  // cout << "done merging\n";
}

void RealPartonState::doEvolution() {
  // cout << "entering doEvolution" << endl;
  for ( RealPartonSet::iterator it = toCheck.begin(); it != toCheck.end(); it++ ) {
    RealPartonPtr rp = *it;
    currentY = rp->theParton->oY();
    if ( !rp->quickSetYES() ) {
      rp->checkEmissionProblem();
      break;
    }
    rp->checkEmissionProblem();
  }
}

bool RealPartonState::findConsistentEvolution(RealPartonSet::iterator it) {
  if ( it == toCheck.end() ) {
    return true;}
  partonCalls++;
  if ( partonCalls > 1000000 ) {return false;}
  tRealPartonPtr rp = *it;
  currentY = rp->theParton->oY();
  if ( !(rp->setYES()) ) {
    unOrdered++;
    if ( monitored && interactions.size() > 1 )
      cout << "unordered" << endl;
    // cout << rp->oY() << " unordered" << endl;
    return rp->setNO();
  }
  if ( !rp->theParton->valence() && (!rp->mothers.first || !rp->mothers.second) ) {
    if ( monitored && interactions.size() > 1 )
      cout << "non-valence without mothers, NO" << endl;
    // cout << rp->oY() << "  mother structure error" << endl;
    return rp->setNO();
  }
  if ( rp->mothers.first && rp->mothers.first == rp->mothers.second ) {
    if ( monitored && interactions.size() > 1 )
      cout << "samme mothers" << endl;
    // cout << rp->oY() << "  mother structure error" << endl;
    return rp->setNO();
  }
  // cout << rp->oY() << " passed first tests. check." << endl;
  bool existHistory = false;
  while ( !existHistory ) {
    existHistory = findConsistentEvolution(++it);
    currentY = rp->theParton->oY();
    if ( existHistory && !(rp->DGLAPSafe()) ) {
      if ( monitored && interactions.size() > 1 )
	cout << "unDGLAP" << endl;
      // cout << rp->oY() << ", (" << rp->theParton->position().x()*GeV << ", "
      // 	   << rp->theParton->position().y()*GeV << ") is unDGLAP" << endl;
      unDGLAP++;
      rp->unDGLAP++;
      existHistory = false;
    }
    if ( existHistory && !(rp->interacting) && !(rp->theParton->valence()) && !(rp->hasChild()) )
      existHistory = false;
    // if ( !existHistory && it != toCheck.end() ) {
    //   cout << rp->oY() << " found no history at " << (*it)->oY() << ", set future NO" << endl;
    //   // plotState(true);
    // }
    if ( !existHistory ) setFutureNO(rp); //behovs denna verkligen??
  if ( !existHistory && (it == toCheck.end() || (*it)->interacting) ) return rp->setNO();
  if ( !existHistory && (*it)->theParton->valence() ) return rp->setNO();
  }
  return existHistory;
}

bool RealPartonState::controlEvolution(tDipolePtr intDip, tDipolePtr otherIntDip) {
  if ( failedInts.find(intDip) != failedInts.end() ) return false;
  //think through how to do this without redoing the evolution!
  if ( intDip->interacted() )  {
    rescatter = true;
//     return false; //removes rescattering
  }
  else  rescatter = false;
  unOrdered = 0;
  unDGLAP = 0;	
  partonCalls = 0;			   
  newInteraction(intDip, otherIntDip, true, true);
  // checkOnlyNew(intDip);
  checkHistory(getReal(intDip->partons().first));
  checkHistory(getReal(intDip->partons().second));
  checkFuture();
  cout << "new evo------------------- " << toCheck.size() << " to check ---------------------" << endl;
  cout << "-- interaction partons are " <<  getReal(intDip->partons().first)->oY() << " and "
       << getReal(intDip->partons().second)->oY() << endl;
  // checkInteraction(intDip);
  // cout << "starting real evo" << endl;
  // diagnosis(true);
  if ( (getReal(intDip->partons().first)->interactions.size() > 1 &&
	getReal(intDip->partons().first)->keep == RealParton::NO) ||
       (getReal(intDip->partons().second)->interactions.size() > 1 &&
	getReal(intDip->partons().second)->keep == RealParton::NO) ) {
    cout << "interacting valence NO" << endl;
    revertToPrevious(intDip);
    failedInts.insert(intDip);
    return false;
  }
  if ( !toCheck.empty() && !rescatter ) { //think through rescatter. will intpT still be ok?
    for ( RealPartonSet::iterator it = toCheck.begin(); !findConsistentEvolution(it); it++ ) {
      if ( (*it)->interacting || (*it)->theParton->valence() ) {
	cout << "failed evo, unordered: " << unOrdered << ", unDGLAP: " << unDGLAP << endl;
	revertToPrevious(intDip);
	failedInts.insert(intDip);
	return false;
      }
    }
  }
  if ( !toCheck.empty() )
    currentY = (**(--toCheck.end())).theParton->oY();
  cout << "found evo" << endl;
  // for ( RealPartonSet::iterator it = toCheck.begin(); it != toCheck.end(); it++ )
  //   (*it)->checkInteractionRecoils();
  return true;
}

bool RealPartonState::singleControlEvolution(tDipolePtr intDip, tDipolePtr otherIntDip,
					     bool firstInt, bool secondInt,
					     Energy rec1, Energy rec2) {
  ostream & log = CurrentGenerator::log();
  static DebugItem printstepdetails("DIPSY::PrintStepDetails", 6);
  static DebugItem checkkinematics("DIPSY::CheckKinematics", 6);
  if ( failedInts.find(intDip) != failedInts.end() ) return false;
  if ( intDip->interacted() )  {
    if ( monitored )
      cout << "  rescatter, autopass evolution" << endl;
    rescatter = true;
  }
  else  rescatter = false;
  if ( printstepdetails ) 
    log << "<DIPSY> enter singlecontrolevolution." << endl;
  unOrdered = 0;
  unDGLAP = 0;	
  partonCalls = 0;
  newInteraction(intDip, otherIntDip, firstInt, secondInt, rec1, rec2);
  reset();
  checkAllInteracting();

  for ( RealPartonSet::iterator it = toCheck.begin(); it != toCheck.end(); it++ )
    if ( !(*it)->setYES() )  fixUnOrdered(*it);

  if ( printstepdetails ) log << "<DIPSY> done with forward sweep." << endl;
  if ( monitored ) plotState(true);

  for ( RealPartonSet::reverse_iterator it = toCheck.rbegin(); it != toCheck.rend(); it++ )
    checkFixDGLAP(*it);

  if ( printstepdetails )
    log << "<DIPSY> done with DGLAP corrections." << endl;
  if ( monitored ) plotState(true);


  for ( RealPartonSet::iterator it = toCheck.begin(); it != toCheck.end(); it++ ) {
    if ( inFSRRegion(*it) )  {
      if ( printstepdetails )
	log << "<DIPSY> FSR conflict found (check) at " << (*it)->oY() << endl;
      fixFSRDoubleCounting(*it);
      if ( printstepdetails ) log << "      fixed, check." << endl;
    }
  }

  if ( printstepdetails ) log << "<DIPSY> done with FSR matching" << endl;
  if ( monitored ) plotState(true);

  //make sure all fluctuations are collinear. In general not needed.
  for ( int i = 0; i < int(flucts.size()); i++ ) {
    makeCollinear(flucts[i]);
  }


  if ( checkkinematics && checkForNegatives() )
    Throw<RealPartonKinematicException>()
      << "negatives found at end of cascade evo!" << Exception::warning;
  return true;
}

bool RealPartonState::nonRecursiveControlEvolution(tDipolePtr intDip, tDipolePtr otherIntDip) {
  if ( failedInts.find(intDip) != failedInts.end() ) return false;
  if ( intDip->interacted() )  {
    if ( monitored )
      cout << "  rescatter, autopass evolution" << endl;
    rescatter = true;
  }
  else  rescatter = false;
  // cout << "entering nonrecursivecontrol\n";
  // diagnosis(false);
  unOrdered = 0;
  unDGLAP = 0;	
  partonCalls = 0;
  newInteraction(intDip, otherIntDip, true, true);
  reset();
  checkAllInteracting();
  // cout << "enter loop to find evo------------------------------------------" << endl;
  while ( true ) {
    undoEmissions();
    doEvolution();
    // cout << "did evo" << endl;
    RealPartonPtr worstProblem = findWorstProblem();
    // if ( worstProblem )
    //   cout << "worst problem at " << worstProblem->oY() << ", ("
    // 	   << worstProblem->theParton->position().x()*GeV << ", "
    // 	   << worstProblem->theParton->position().y()*GeV << ")" << endl;
    // else {
    //   cout << "no problem! evo ok for " << this << "" << endl;
    //   // diagnosis(false);
    // }
    if ( !worstProblem ) return true;
    if ( !fix(worstProblem) ) {
      // cout << "failed fixing theproblem, interaction failed" << endl;
      revertToPrevious(intDip);
      failedInts.insert(intDip);
      return false;
    }
  }
}

bool RealPartonState::fullControlEvolution(tDipolePtr intDip, tDipolePtr otherIntDip,
					   bool firstInt, bool secondInt,
					   Energy rec1, Energy rec2) {
  if ( intDip->dipoleState().handler().eventFiller().mode() == 1 )
    return controlEvolution(intDip, otherIntDip);
  if ( intDip->dipoleState().handler().eventFiller().mode() == 2 )
    return singleControlEvolution(intDip, otherIntDip, firstInt, secondInt, rec1, rec2);
  if ( intDip->dipoleState().handler().eventFiller().mode() == 3 )
    return nonRecursiveControlEvolution(intDip, otherIntDip);
  if ( failedInts.find(intDip) != failedInts.end() ) return false;
  //think through how to do this without redoing the evolution!
  if ( intDip->interacted() )  {
    if ( monitored )
      cout << "  rescatter, autopass evolution" << endl;
    rescatter = true;
    cout << "rescatter!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n";
//     return false; //removes rescattering
  }
  else  rescatter = false;
  unOrdered = 0;
  unDGLAP = 0;	
  partonCalls = 0;			   
  newInteraction(intDip, otherIntDip, true, true);
  reset();
  checkAllInteracting();
  cout << "new evo---------------------------------------------------------" << endl;
  for ( RealPartonSet::iterator it = toCheck.begin(); !findConsistentEvolution(it); it++ ) {
    if ( (*it)->interacting || (*it)->theParton->valence() ) {
      // cout << "failed evo, unordered: " << unOrdered << ", unDGLAP: " << unDGLAP << endl;
      // for ( RealPartonSet::iterator jt = toCheck.begin(); jt != toCheck.end(); jt++ ) {
      // 	(*jt)->setYES();
      // }
      // plotState(true);
      revertToPrevious(intDip);
      failedInts.insert(intDip);
      return false;
    }
  }
  cout << "  found consistent evolution!" << endl;
  // diagnosis(false);
  // cout << "checked" << endl;
  currentY = (**(--toCheck.end())).theParton->oY();
  // cout << "  did " << partonCalls << " calls on " << toCheck.size() << " partons.";
  return true;
}

void RealPartonState::reset() {
  for ( RealPartonSet::iterator it = oValence.begin(); it != oValence.end(); it++ ) {
    (*it)->setYES();
  }
  for ( map<tPartonPtr,RealPartonPtr>::iterator it = partons.begin();
	it != partons.end(); it++ ) {
    (*it).second->givenMinus = ZERO;
    (*it).second->fluct = -1;
    (*it).second->DGLAPchecked = false;
    if ( (*it).second->keep == RealParton::YES && !(*it).second->theParton->valence() ) {
      (*it).second->setNO();
    }
    (*it).second->intRecoil = TransverseMomentum();
  }
  // for ( RealPartonSet::iterator it = valence.begin(); it != valence.end(); it++ ) {
  //   (*it)->setYES();
  // }
  totalRecoil = TransverseMomentum();
  backwardsPartons.clear();
  flucts.clear();
  suspects.clear();
  toCheck.clear();
  minusDeficit = (*interactions.begin())->dipoleState().minusDeficit();
}

void RealPartonState::undoEmissions() {
  for ( RealPartonSet::reverse_iterator it = toCheck.rbegin();
	it != toCheck.rend(); it++ ) {
    RealPartonPtr rp = *it;
    if ( !rp->theParton->valence() )
      rp->setNO();
    else {
      if ( rp->mothers.first ) rp->mothers.first->children.second.erase(rp);
      if ( rp->mothers.second ) rp->mothers.second->children.first.erase(rp);
      rp->quickSetYES();
    }
  }
}

void RealPartonState::checkOnlyNew(tDipolePtr intDip) {
  cout << "implement checkonlynew if you really want to use it!!" << endl;
  // toCheck.clear();
  // //interacting partons cannot be left NO, even if earlier evos did so.
  // if ( getReal(intDip->partons().first)->keep == RealParton::NO )
  //   toCheck.insert(getReal(intDip->partons().first));
  // if ( getReal(intDip->partons().second)->keep == RealParton::NO )
  //   toCheck.insert(getReal(intDip->partons().second));
  // for ( map<tPartonPtr,RealPartonPtr>::iterator it = partons.begin();
  // 	it != partons.end(); it++ ) {
  //   if ( ((*it).second->interactions.size() == 1) &&
  // 	 (*it).second->interactions.find(intDip) != (*it).second->interactions.end() ) {
  //     toCheck.insert(it->second);
  //   }
  // }
  // currentY = 0.0;
}

void RealPartonState::checkInteraction(tDipolePtr intDip) {
  toCheck.clear();
  RealPartonPtr rp1 = getReal(intDip->partons().first);
  checkHistory(rp1);
  RealPartonPtr rp2 = getReal(intDip->partons().second);
  checkHistory(rp2);
}

void RealPartonState::checkHistory(tRealPartonPtr rp) {
  toCheck.insert(rp);
  if ( rp->oMothers.first )
    checkHistory(rp->oMothers.first);
  if ( rp->oMothers.second )
    checkHistory(rp->oMothers.second);
}

void RealPartonState::checkFuture() {
  RealPartonSet future;
  for ( RealPartonSet::iterator it = toCheck.begin(); it != toCheck.end(); it++ ) {
    RealPartonPtr rp = *it;
    if ( !rp->future.empty() )
      future.insert(rp->future.begin(), rp->future.end());
    }
  if ( !future.empty() )
    toCheck.insert(future.begin(), future.end());
}

void RealPartonState::checkAllInteracting() { //goes through all partons.
  for ( map<tPartonPtr,RealPartonPtr>::iterator it = partons.begin();
	it != partons.end(); it++ ) {
    if ( !((*it).second->interactions.empty()) ) {
      toCheck.insert(it->second);
    }
    if ( it->first->valence() ) it->second->setValenceMomentum();
  }
  currentY = 0.0;
}

void RealPartonState::setFutureNO(tRealPartonPtr rp) { //remade to go from high y to low/
  for ( RealPartonSet::reverse_iterator it = toCheck.rbegin(); *it != rp; it++) {
    if ( (*it)->keep == RealParton::YES && !(*it)->theParton->valence() )
      (*it)->setNO();
  }
}

void RealPartonState::setNO(const RealPartonSet & rps) {
  for ( RealPartonSet::const_reverse_iterator it = rps.rbegin();
	it != rps.rend(); it++)
    if ( (*it)->keep == RealParton::YES &&
	 !(*it)->theParton->valence() ) (*it)->setNO();
}

RealPartonPtr RealPartonState::findWorstProblem() {
  // cout << "entering worst problem" << endl;
  Energy largestScale = ZERO;
  RealPartonPtr ret;
  for ( RealPartonSet::iterator it = toCheck.begin(); it != toCheck.end(); it++) {
    if ( (*it)->keep == RealParton::NO ) break;
    Energy scale = (*it)->problemScale();
    // cout << "got scale " << scale/GeV << " at " << (*it)->oY() << ", ("
    // 	 << (*it)->theParton->position().x()*GeV << ", "
    // 	 << (*it)->theParton->position().y()*GeV << ")" << endl;
    if ( scale > largestScale ) {
      largestScale = scale;
      ret = *it; 
    }
  }
  return ret;
}

bool RealPartonState::fix(RealPartonPtr problem) {
  // cout << "entered fix" << endl;
  RealPartonPtr cause = problem->findCause();
  if ( !cause ) return false;
  cause->setNO();
  toCheck.erase(cause);
  return true;
}

void RealPartonState::saveState() {
  if ( Debug::level > 5 ) cout << "entering RPS::saveState" << endl;
  for (   map<tPartonPtr, RealPartonPtr>::iterator it = partons.begin();
	  it != partons.end(); it++) {
    (*it).second->saveState();
  }
  for (   RealPartonSet::iterator it = valence.begin(); //can be removed?
	  it != valence.end(); it++) {
    (*it)->saveState();
  }
  failedInts.clear(); //when interaction found, the fails may now be ok.
  if ( Debug::level > 5 ) cout << "  failedInts cleared" << endl;
  cTotalRecoil = totalRecoil;
  cFlucts = flucts;
}

void RealPartonState::revertToPrevious(DipolePtr intDip) {
  if ( Debug::level > 5 ) cout << "revert to previous called with intdip " << intDip << endl;
  // if ( rescatter ) {
  //   //if rescattering, some extra care has to be taken in
  //   //removing the interaction tags from the partons,
  //   //so dont remove any tags the normal way, and call removeLastRescatter instead. 
  //   removeLastRescatter(intDip);
  //   intDip = DipolePtr();
  // }
  for (   map<tPartonPtr,RealPartonPtr>::iterator it = partons.begin();
	  it != partons.end(); it++) {
    it->second->revert(intDip);
  }
  interactions.pop_back();
  doesInts.pop_back();
  totalRecoil = cTotalRecoil;
  flucts = cFlucts;
}

void RealPartonState::removeLastRescatter(tDipolePtr intDip) {
  if ( Debug::level > 5 ) cout << "entering removeLastRescatter for intdip " << intDip << endl;

  //first check which of the two partons has been interacting in previous interactions
  bool firstOldInt = false;
  bool secondOldInt = false;
  //loop over all but the last interaction, to not include this last rescattering
  if ( Debug::level > 5 ) cout << "  number of interactions total: " << interactions.size() << endl;
  list<pair<bool, bool> > ::iterator jt = doesInts.begin();
  for ( list<tDipolePtr>::iterator it = interactions.begin(); it != --interactions.end();
	it++, jt++ ) {
    if ( Debug::level > 5 ) cout << "  old int " << *it << endl;
    if ( *it == intDip ) {
      if ( Debug::level > 5 ) cout << "  found old rescattering" << endl;
      if ( jt->first ) firstOldInt = true;
      if ( jt->second ) secondOldInt = true;
    }
  }

  //if one of the partons was not tagged before, but got tagged by this last rescatter,
  //then those new tags should be removed, otherwise nothing has to be done.
  tRealPartonPtr oldInt;
  if ( firstOldInt == false && doesInts.rbegin()->first ) {
    if ( Debug::level > 5 ) cout << "  reverting tags" << endl;

    //the tags should be only from the second parton, so remove first, and readd second
    //only removing ancestors of first does not work, since they may share ancestors
    removeInt(getReal(intDip->partons().first), intDip);
    addInt(getReal(intDip->partons().second), intDip);
  }
  else if ( secondOldInt == false && doesInts.rbegin()->second ) {
    if ( Debug::level > 5 ) cout << "  reverting tags" << endl;
    removeInt(getReal(intDip->partons().second), intDip);
    addInt(getReal(intDip->partons().first), intDip);
  }

  if ( Debug::level > 5 ) cout << "done in removeLastRescatter" << endl;
}

void RealPartonState::removeInt(tRealPartonPtr rp, DipolePtr intDip) {
  rp->interactions.erase(intDip);
  if ( rp->oMothers.first ) removeInt(rp->oMothers.first, intDip);
  if ( rp->oMothers.second ) removeInt(rp->oMothers.second, intDip);
  if ( rp->oMother ) removeInt(rp->oMother, intDip);
}

void RealPartonState::addInt(tRealPartonPtr rp, DipolePtr intDip) {
  rp->interactions.insert(intDip);
  if ( rp->oMothers.first ) addInt(rp->oMothers.first, intDip);
  if ( rp->oMothers.second ) addInt(rp->oMothers.second, intDip);
  if ( rp->oMother ) addInt(rp->oMother, intDip);
}

Energy RealPartonState::neededValenceMinus() {
  if ( interactions.empty() ) cerr << "neededvalenceminus called for a state "
				" without interaction!" << endl;
  Energy ret = minusDeficit;
  minusDeficit = ZERO;
  for ( RealPartonSet::iterator it = valence.begin();
	it != valence.end(); it++ ) {
    ret += (**it).giveMinus();
  }
  return ret;
}

void RealPartonState::setOnShell(tDipolePtr intDip) {
  getReal(intDip->partons().first)->setOnShell();
  getReal(intDip->partons().second)->setOnShell();
}

void RealPartonState::removeBackwards() {
  for ( RealPartonSet::const_iterator it = backwardsPartons.begin();
	it != backwardsPartons.end(); it++ ) {

  }
}

Energy RealPartonState::totalMinus() {
  Energy ret = ZERO;
  for (   map<tPartonPtr,RealPartonPtr>::const_iterator it = partons.begin();
	  it != partons.end(); it++) {
    tRealPartonPtr rp = (*it).second;
    if ( rp->keep == RealParton::YES ) {
      ret += rp->minus;
    }
  }
  return ret;
}

Energy RealPartonState::totalPlus() {
  Energy ret = ZERO;
  for (   map<tPartonPtr,RealPartonPtr>::const_iterator it = partons.begin();
	  it != partons.end(); it++) {
    tRealPartonPtr rp = (*it).second;
    if ( rp->keep == RealParton::YES ) {
      ret += rp->plus;
    }
  }
  return ret;
}

void RealPartonState::addRecoilers() {
  if ( Debug::level > 5 ) {
    cout << "before recoilers" << endl;
    // plotState(true);
  }
  //loop over all active partons
  for (   map<tPartonPtr,RealPartonPtr>::const_iterator it = partons.begin();
	  it != partons.end(); it++) {
    if ( (*it).second->keep == RealParton::NO ) continue;
    //hit holes in the parton from the recoils
    addRecoiler2((*it).second);
  }
  if ( Debug::level > 5 ) {
    cout << "with recoilers" << endl;
    // plotState(true);
  }
}

void RealPartonState::addRecoiler2(RealPartonPtr rp) {
  TransverseMomentum opt = rp->theParton->valencePT();
  if ( rp->mothers.first ) opt += rp->theParton->recoil(rp->mothers.first->theParton);
  if ( rp->mothers.second ) opt += rp->theParton->recoil(rp->mothers.second->theParton);
  if ( rp->nMothers == 1 && rp->mother ) opt = rp->theParton->valencePT() +
					   rp->theParton->recoil(rp->mother->theParton);
  if ( opt.pt() == ZERO ) opt = rp->theParton->pT();

  TransverseMomentum rec = rp->pT - opt;

  if ( rec.pt() > 2.0*opt.pt() ) {
    double plusRatio = opt.pt()/rec.pt();
    rp->emitRecoiler(rec, plusRatio);
    if ( Debug::level > 5 ) {
      cout << "added recoiler for " << rp->oY() << endl;
      // plotState(true);
    }
  }
}

void RealPartonState::addRecoiler(RealPartonPtr rp) {
  PartonPtr p = rp->theParton;
  //loop over all children.

  //if any first children, create extra gluon in both real and dipole state
  if ( !rp->children.first.empty() || rp->firstInt ) {
    if ( Debug::level > 5 ) {
      cout << "try to create first recoiler to " << rp->oY() << endl;
    }
    //set up
    DipolePtr dip;
    if ( p->dipoles().first ) dip = p->dipoles().first;
    else {
      Throw<RealPartonKinematicException>()
      	<< "Quark wants to create Recoiler on the wrong side!" << Exception::warning;
      return;
    }
    TransverseMomentum pt = TransverseMomentum();
    double ymax = rp->oY();

    //move pt from rp to the new recoiler
    if ( !rp->children.first.empty() ) {
      for ( RealPartonSet::iterator it = rp->children.first.begin();
	    it != rp->children.first.end(); it++ ) {
	RealPartonPtr child = *it;
	//loop through "recoils" and move all pt recoils to the recoiler
	for ( list< RealParton::Recoil >::iterator jt = child->recoils.begin();
	      jt != child->recoils.end(); jt++ ) {
	  cout << child->oY() << " has a recoil of " << jt->second.first.pt()/GeV
	       << " to " << jt->first->oY() << endl;
	  if ( jt->first != rp ) continue;
	  if ( ymax < (*it)->y ) ymax = (*it)->y;
	  pt -= jt->second.first;
	  if ( Debug::level > 5 )
	    cout << child->oY() << " gives a pt kick of " << jt->second.first.pt()/GeV
		 << " to " << rp->oY() << endl;
	  if ( Debug::level > 5 )
	    cout << child->oY() << " recoilers pt is now " << pt.pt()/GeV << endl;
	}
	// if ( child->nMothers == 1 ) child->mother = realRecoiler;
	// if (child->nMothers == 2 ) child->mothers.second = realRecoiler;
      }
    }

    if ( rp->firstInt ) {
      vector<RealPartonPtr> ints;
      //do the same, but use intDist to set the recoil
      for ( list< RealParton::Recoil >::iterator it = rp->recoils.begin();
	    it != rp->recoils.end(); it++ ) {
	if ( it->first->realState != this ) {
	  ints.push_back(it->first);
	}
	for (int i = 0; i < int(ints.size()); i++ ) {
	  //loop through "recoils" and move all pt recoils to the recoiler
	  for ( list< RealParton::Recoil >::iterator jt = ints[i]->recoils.begin();
		jt != ints[i]->recoils.end(); jt++ ) {
	    cout << -ints[i]->oY() << " has an int recoil of " << jt->second.first.pt()/GeV
		 << " to " << jt->first->oY() << endl;
	    if ( jt->first != rp ) continue;
	    if ( ymax < -ints[i]->y ) ymax = -ints[i]->y;
	    pt -= jt->second.first;
	    if ( Debug::level > 5 )
	      cout << ints[i]->oY() << " from other state gives a pt kick of " << jt->second.first.pt()/GeV
		   << " to " << rp->oY() << endl;
	    if ( Debug::level > 5 )
	      cout << ints[i]->oY() << " recoilers pt is now " << pt.pt()/GeV << endl;
	  }
	}
      }


      //maybe take direction from current pt?
      //can "recoils" be used? does it store int recoils?
    }

    //if the hole is smaller than the mother, kick out a recoiler
    if ( pt.pt() > (rp->pT - pt).pt() ) {
      //create gluons
      PartonPtr recoiler = new_ptr(Parton());
      recoiler->onShell(true);
      dip->generatedGluon(recoiler);
      //fix colour flow for dipole state (let it be connected to the original one)
      dip->splitDipole(0.5); //use better approximation for colour choice?
      //initialise real recoiler
      RealPartonPtr realRecoiler = getReal(recoiler);
      realRecoiler->nMothers = 1;
      realRecoiler->setOMother(rp);
      realRecoiler->mother = rp;
      realRecoiler->keep = RealParton::YES;

      //set momentum of recoiler
      realRecoiler->pT = pt;
      rp->pT -= pt;
      //otherwise use ymax to limit how far the recoiler can go in rapidity
      double plusRatio = rp->pT.pt()/(realRecoiler->pT.pt() + rp->pT.pt());
      realRecoiler->plus = plusRatio*rp->plus;
      rp->plus = (1.0 - plusRatio)*rp->plus;
      realRecoiler->updateYMinus();
      rp->updateYMinus();
      recoiler->oY(realRecoiler->y);

      rp->saveState();
      realRecoiler->saveState();

      if ( Debug::level > 5 )
	cout << "recoiler takes plus ratio " << plusRatio << endl;

      //set position of recoiler
      recoiler->position(p->position() + p->pTScale()*pt/sqr(pt.pt()));

      //set mother structure
      if ( !rp->firstInt ) {
	RealPartonPtr lastChild = *rp->children.first.rbegin();
	if ( lastChild->nMothers == 1 ) lastChild->mother = realRecoiler;
	if ( lastChild->nMothers == 2 ) lastChild->mothers.second = realRecoiler;
	realRecoiler->children.first.insert(lastChild);
      }
      else {
	realRecoiler->interacting = rp->interacting;
	realRecoiler->intDist = rp->intDist;
	realRecoiler->firstInt = true;
      }

      if ( Debug::level > 5 ) {
	cout << "recoiler places at transverse distance " << p->pTScale()*rp->pT.pt()/sqr(rp->pT.pt())*GeV << endl;
	cout << "recoiler got pt " << realRecoiler->pT.pt()/GeV << endl;
	cout << "done with first recoiler of " << rp->oY() << endl;
	// dip->dipoleState().plotState(false);
	// plotState(true);
      }
    }
    else {
      if ( Debug::level > 5 ) {
	cout << "not enough pt, dont do recoiler" << endl;
	cout << "recoiler pt was " << pt.pt()/GeV << ", canceled parent pt "
	     << (rp->pT - pt).pt()/GeV << endl;
	// plotState(true);
      }
    }

  }

  //same for other side
  //if any first children, create extra gluon in both real and dipole state
  if ( !rp->children.second.empty() || rp->secondInt ) {
    if ( Debug::level > 5 ) {
      cout << "try to create second recoiler to " << rp->oY() << endl;
    }
    //set up
    DipolePtr dip;
    if ( p->dipoles().second ) dip = p->dipoles().second;
    else {
      Throw<RealPartonKinematicException>()
      	<< "Quark wants to create Recoiler on the wrong side!" << Exception::warning;
      return;
    }
    TransverseMomentum pt = TransverseMomentum();
    double ymax = rp->oY();

    //move pt from rp to the new recoiler
    if ( !rp->children.second.empty() ) {
      for ( RealPartonSet::iterator it = rp->children.second.begin();
	    it != rp->children.second.end(); it++ ) {
	RealPartonPtr child = *it;
	if ( rp->children.first.find(child) != rp->children.first.end() ) {
	  cout << "found second child at " << child->oY()
	       << "  that is also first. skip!" << endl;
	  plotState(true);
	  continue;
	}
	//loop through "recoils" and move all pt recoils to the recoiler
	for ( list< RealParton::Recoil >::iterator jt = child->recoils.begin();
	      jt != child->recoils.end(); jt++ ) {
	  cout << child->oY() << " has a recoil of " << jt->second.first.pt()/GeV
	       << " to " << jt->first->oY() << endl;
	  if ( jt->first != rp ) continue;
	  if ( ymax < (*it)->y ) ymax = (*it)->y;
	  pt -= jt->second.first;
	  if ( Debug::level > 5 )
	    cout << child->oY() << "gives a pt kick of " << jt->second.first.pt()/GeV
		 << " to " << rp->oY() << endl;
	  if ( Debug::level > 5 )
	    cout << child->oY() << "recoilers pt is now " << pt.pt()/GeV << endl;
	}
	// if ( child->nMothers == 1 ) child->mother = realRecoiler;
	// if (child->nMothers == 2 ) child->mothers.second = realRecoiler;
      }
    }

    if ( rp->secondInt ) {
      //do the same, but use intDist to set the recoil
      //maybe take direction from current pt?
      //can "recoils" be used? does it store int recoils?
    }

    //if the hole is smaller than the mother, kick out a recoiler
    if ( pt.pt() > (rp->pT - pt).pt() ) {
      //create gluons
      PartonPtr recoiler = new_ptr(Parton());
      recoiler->onShell(true);
      dip->generatedGluon(recoiler);
      //fix colour flow for dipole state (let it be connected to the original one)
      dip->splitDipole(0.5); //use better approximation for colour choice?
      //initialise real recoiler
      RealPartonPtr realRecoiler = getReal(recoiler);
      realRecoiler->nMothers = 1;
      realRecoiler->setOMother(rp);
      realRecoiler->mother = rp;
      realRecoiler->keep = RealParton::YES;

      //set momentum of recoiler
      realRecoiler->pT = pt;
      rp->pT -= pt;
      //otherwise use ymax to limit how far the recoiler can go in rapidity
      double plusRatio = rp->pT.pt()/(realRecoiler->pT.pt() + rp->pT.pt());
      realRecoiler->plus = plusRatio*rp->plus;
      rp->plus = (1.0 - plusRatio)*rp->plus;
      realRecoiler->updateYMinus();
      rp->updateYMinus();
      recoiler->oY(realRecoiler->y);

      rp->saveState();
      realRecoiler->saveState();

      if ( Debug::level > 5 )
	cout << "recoiler takes plus ratio " << plusRatio << endl;

      //set position of recoiler
      recoiler->position(p->position() + p->pTScale()*pt/sqr(pt.pt()));

      //set mother structure
      if ( !rp->secondInt ) {
	RealPartonPtr lastChild = *rp->children.second.rbegin();
	if ( lastChild->nMothers == 1 ) lastChild->mother = realRecoiler;
	if ( lastChild->nMothers == 2 ) lastChild->mothers.first = realRecoiler;
	realRecoiler->children.second.insert(lastChild);
      }
      else {
	realRecoiler->interacting = rp->interacting;
	realRecoiler->intDist = rp->intDist;
	realRecoiler->secondInt = true;
      }

      if ( Debug::level > 5 ) {
	cout << "recoiler places at transverse distance " << p->pTScale()*rp->pT.pt()/sqr(rp->pT.pt())*GeV << endl;
	cout << "recoiler got pt " << realRecoiler->pT.pt()/GeV << endl;
	cout << "done with first recoiler of " << rp->oY() << endl;
	// dip->dipoleState().plotState(false);
	// plotState(true);
      }
    }
    else {
      if ( Debug::level > 5 ) {
	cout << "not enough pt, dont do recoiler" << endl;
	cout << "recoiler pt was " << pt.pt()/GeV << ", canceled parent pt "
	     << (rp->pT - pt).pt()/GeV << endl;
	// plotState(true);
      }
    }

  }

  cout << "done with " << rp->oY() << endl;
  plotState(true);

}

void RealPartonState::plotState(bool pause) const {
  cout << "print state to realState.dat" << endl;
  ofstream file ("realState.dat");
  plotStateToFile(file, false);
  file.close();
  if ( pause && Debug::level > 5 ) {
    cout << "press any key to continue..." << endl;
    getchar();
  }
}

void RealPartonState::plotBothStates(RealPartonStatePtr rrs, bool pause) const {
  cout << "print state to realState.dat" << endl;
  ofstream file ("realState.dat");
  plotStateToFile(file, false);
  rrs->plotStateToFile(file, true);
  file.close();
  if ( pause && Debug::level > 5 ) {
    cout << "press any key to continue..." << endl;
    getchar();
  }
}

void RealPartonState::plotStateToFile(ostream & file, bool m) const {
  double ptmax = -100;
  double ptmin = 100000;

  double ymax = 0.0;
  PartonPtr p = partons.begin()->second->theParton;
  if ( p->dipoles().first ) ymax = p->dipoles().first->dipoleState().ymax();
  if ( p->dipoles().second ) ymax = p->dipoles().second->dipoleState().ymax();

  for ( map<tPartonPtr,RealPartonPtr>::const_iterator it = partons.begin();
	  it != partons.end(); it++) {
    RealPartonPtr rp = it->second;
    if ( rp->pT.pt()/GeV > ptmax ) ptmax = rp->pT.pt()/GeV;
    if ( rp->pT.pt()/GeV < ptmin ) ptmin = rp->pT.pt()/GeV;
    bool val = rp->theParton->valence();
    TransverseMomentum opt = rp->theParton->valencePT();
    if ( rp->mothers.first ) opt += rp->theParton->recoil(rp->mothers.first->theParton);
    if ( rp->mothers.second ) opt += rp->theParton->recoil(rp->mothers.second->theParton);
    if ( rp->nMothers == 1 && rp->mother ) opt = rp->theParton->valencePT() +
					     rp->theParton->recoil(rp->mother->theParton);
    if ( opt.pt() == ZERO ) opt = rp->theParton->pT();
    //the main real partons
    file << ((!m) ? rp->theParton->oY():-rp->theParton->oY()) << '\t'
	 << rp->theParton->position().y()*GeV << '\t'
	 << rp->keep << '\t'
	 << ((!m) ? rp->y:-rp->y) << '\t'
	 << val << '\t'
	 << opt.pt()/GeV << '\t'
	 << bool(toCheck.find(rp) != toCheck.end()) << '\t'
	 << rp->theParton->position().x()*GeV << endl << endl;

    //info on the consistent state, N/A if state never saved.
    if ( interactions.size() > 1 || cTotalRecoil == totalRecoil ) {
      file << ((!m) ? rp->theParton->oY():-rp->theParton->oY()) << '\t'
	   << opt.pt()/GeV << '\t'
	   << 3 << '\t' <<  endl;
      file << ((!m) ? rp->theParton->y():abs(rp->theParton->y())) << '\t'
	   << rp->theParton->pT().pt()/GeV << '\t'
	   << 4 << '\t' << endl << endl;
    }

    //the interacting marks
    if ( rp->interacting )
      file << ((!m) ? rp->theParton->oY():-rp->theParton->oY()) << '\t'
	   << rp->theParton->position().y()*GeV << '\t'
	   << 2 << '\t'
	   << ((!m) ? rp->y:-rp->y) << '\t'
	   << val << '\t'
	   << opt.pt()/GeV << '\t'
	   << rp->theParton->position().x()*GeV << endl
	   << ymax << '\t'
	   << rp->theParton->position().y()*GeV << '\t'
	   << 2 << '\t'
	   << ((!m) ? rp->y:-rp->y) << '\t'
	   << val << '\t'
	   << rp->theParton->pTScale()/rp->intDist/GeV << '\t'
	   << rp->theParton->position().x()*GeV << endl << endl;

    //line to current state
    file << ((!m) ? rp->theParton->oY():-rp->theParton->oY()) << '\t'
	 << rp->theParton->position().y()*GeV << '\t'
	 << rp->keep << '\t'
	 << 7 << '\t'
	 << opt.pt()/GeV << endl
	 << ((!m) ? rp->y:-rp->y) << '\t'
	 << rp->theParton->position().y()*GeV << '\t'
	 << rp->keep << '\t'
	 << 8 << '\t'
	 << rp->pT.pt()/GeV <<  endl << endl;

    //the pT
    file << ((!m) ? rp->theParton->oY():-rp->theParton->oY()) << '\t'
	 << rp->theParton->position().y()*GeV << '\t'
	 << ((rp->keep == RealParton::YES) ? 9:10) << '\t'
	 << rp->theParton->position().x()*GeV << '\t'
	 << 0.0 << '\t'
	 << 0.0 << endl
	 << ((!m) ? rp->theParton->oY():-rp->theParton->oY()) << '\t'
	 << rp->theParton->position().y()*GeV << '\t'
	 << ((rp->keep == RealParton::YES) ? 9:10) << '\t'
	 << rp->theParton->position().x()*GeV << '\t'
	 << rp->pT.x()/GeV << '\t'
	 << rp->pT.y()/GeV << endl << endl;

    vector<RealPartonPtr> mots;
      mots.push_back(rp->oMothers.first);
      mots.push_back(rp->oMothers.second);
      mots.push_back(rp->mothers.first);
      mots.push_back(rp->mothers.second);
      mots.push_back(rp->oMother);
      mots.push_back(rp->mother);

    //lines to mothers and original mothers
    for ( int i = 0; i < int(mots.size()); i++) {
      RealPartonPtr rp2 = rp;
      if ( mots[i] ) rp2 = mots[i];
    TransverseMomentum opt2 = rp2->theParton->valencePT();
    if ( rp2->mothers.first ) opt2 += rp2->theParton->recoil(rp2->mothers.first->theParton);
    if ( rp2->mothers.second ) opt2 += rp2->theParton->recoil(rp2->mothers.second->theParton);
    if ( rp2->nMothers == 1 && rp2->mother ) opt2 = rp2->theParton->valencePT() +
					     rp2->theParton->recoil(rp2->mother->theParton);
    if ( opt2.pt() == ZERO ) opt2 = rp2->theParton->pT();
      file << ((!m) ? rp->theParton->oY():-rp->theParton->oY()) << '\t'
	   << rp->theParton->position().y()*GeV << '\t'
	   << rp->keep << '\t'
	   << i+3 + 100 << '\t'
	   << opt.pt()/GeV << '\t'
	   << rp->theParton->position().x()*GeV << '\t'
	   << rp->nMothers + 10 << endl
	   << ((!m) ? rp2->theParton->oY():-rp2->theParton->oY()) << '\t'
	   << rp2->theParton->position().y()*GeV << '\t'
	   << rp2->keep << '\t'
	   << i+3 + 100 << '\t'
	   << opt2.pt()/GeV << '\t'
	   << rp2->theParton->position().x()*GeV << '\t'
	   << rp->nMothers + 10 << endl;
    }
  }

  file << ymax << '\t' << ptmax << '\t' << 10 << '\t' << 10 << endl
       << ymax << '\t' << ptmin << '\t' << 10 << '\t' << 10 << endl;

  file << ((!m) ? currentY:-currentY) << '\t' << ptmax << '\t' << 10 << '\t' << 9 << endl
       << ((!m) ? currentY:-currentY) << '\t' << ptmin << '\t' << 10 << '\t' << 9 << endl;
}

pair<int, int> RealPartonState::countYESNO() const {
  int yes = 0;
  int  no = 0;
  for (   map<tPartonPtr,RealPartonPtr>::const_iterator it = partons.begin();
	  it != partons.end(); it++) {
    if ( (*it).second->keep == RealParton::YES ) yes++;
    else no++;
  }
  return make_pair(yes, no);
}

double RealPartonState::avYInEvo() const {
  double number = 0.0;
  double length = 0.0;
  for (   map<tPartonPtr,RealPartonPtr>::const_iterator it = partons.begin();
	  it != partons.end(); it++) {
    RealPartonPtr rp = (*it).second;
    if ( rp->keep == RealParton::YES && !(rp->theParton->valence()) ) {
      number++;
      if ( rp->nMothers == 2 )
	length += (abs(rp->y - rp->mothers.first->y) + abs(rp->y - rp->mothers.second->y))/2.0;
      if ( rp->nMothers == 1 )
	length += abs(rp->y - rp->mother->y);
    }
  }
  if (number == 0.0 ) return 0.0;
  return length/number;
}

bool RealPartonState::diagnosis(bool pause) const {
  calls++;
  bool ok = true;
//   cout << "diagnosis call number " << calls  << " for " << this << endl;
  TransverseMomentum pT = TransverseMomentum();
  Energy dplus = ZERO;
  Energy dminus = ZERO;
  TransverseMomentum cpT = TransverseMomentum();
  Energy cplus = ZERO;
  Energy cminus = ZERO;
  Energy gminus = ZERO;
  int yes = 0;
  int consInt = 0;
  int vals = 0;
  int mothers = 0;
  int kids = 0;
  int exchanges = 0;
  int interactingPartons = 0;
  TransverseMomentum valPT = TransverseMomentum();
  TransverseMomentum intPT = TransverseMomentum();
  Energy valPlus = ZERO;
  Energy valMinus = ZERO;
  for (   map<tPartonPtr,RealPartonPtr>::const_iterator it = partons.begin();
	  it != partons.end(); it++) {
    tRealPartonPtr rp = (*it).second;
    if ( !rp->exchanges.empty() )
      exchanges += rp->exchanges.size();
    if ( rp->cInteracting ) consInt++;
    if ( rp->theParton->valence() ) {
      vals++;
      valPT += rp->theParton->valencePT();
      valPlus += rp->theParton->valencePlus();
      if ( rp->movedValence )
	valMinus += rp->theParton->valencePlus()*exp(2.0*rp->movedValence->oY());
      else
	valMinus += rp->theParton->valencePlus()*exp(2.0*rp->oY());
      if ( rp->keep == RealParton::NO ) {
	pT += rp->theParton->valencePT();
	dplus += rp->theParton->valencePlus();
	dminus += sqr(rp->theParton->valencePT().pt())/rp->theParton->valencePlus();
	cout << "* valence parton that is NO!!" << endl;
	ok = false;
      }
    }
    if ( rp->keep == RealParton::YES ) {
      yes++;
      pT += rp->pT;
      dplus += rp->plus;
      dminus += rp->minus;
      gminus += rp->givenMinus;
      if ( rp->mothers.first ) mothers++;
      if ( rp->mothers.second ) mothers++;
      if ( rp->mother ) mothers++;
      kids += rp->children.first.size();
      kids += rp->children.second.size();
      if ( rp->interacting ) interactingPartons++;
    }
    if ( rp->cKeep == RealParton::YES ) {
      cpT += rp->theParton->pT();
      cplus += rp->theParton->plus();
      cminus += rp->theParton->minus();
    }
    if ( rp->interacting ) {
      intPT += rp->intRecoil;
    }
//     if ( rp->interacting && toCheck.find(rp) == toCheck.end() ) {
//       cout << "interacting parton at " << rp->theParton->oY() << " is not in tocheck" << endl;
//       ok = false;
//     }
  }

  for (   map<tPartonPtr,RealPartonPtr>::const_iterator it = partons.begin();
	  it != partons.end(); it++) {
    tRealPartonPtr rp = (*it).second;
    if ( rp->keep != RealParton::NO && rp->keep != RealParton::YES )
      Throw<RealPartonKinematicException>()
	<< "Real parton at oY " << rp->oY() << " is neither YES nor NO..." << Exception::warning;
    if ( rp->keep == RealParton::YES ) {
//       if ( !rp->checkMomentum() && !(rp->interacting && totalRecoil.pt() != ZERO )
// 	   && rp->theParton->oY() <= currentY )
// 	ok = false;
      if ( rp->plus < ZERO ) {
	Throw<RealPartonKinematicException>()
	  << "negative plus in real parton diagnosis at oY " << rp->oY() << Exception::warning;
	ok = false;
      }
      if ( rp->minus < ZERO ) {
	Throw<RealPartonKinematicException>()
	  << "negative minus in real parton diagnosis at oY " << rp->oY() << Exception::warning;
	ok = false;
      }
      for ( RealPartonSet::const_iterator it = rp->children.first.begin();
	    it != rp->children.first.end(); it++ ) {
	if ( (*it)->mothers.second != rp && (*it)->mother != rp ) {
	  cout << "\n* parton at " << (*it)->theParton->oY() << " doesnt have " << rp->theParton->oY()
	       << " as second mother! child-mother broken!" << endl;
	  ok = false;
	}
	if ( (*it)->keep == RealParton::NO ) {
	  cout << "\n* parton at " << rp->theParton->oY() << "has a NO-child at "
	       << (*it)->theParton->oY() << endl;
	  ok = false;
	}
      }
      for ( RealPartonSet::const_iterator it = rp->children.second.begin();
	    it != rp->children.second.end(); it++ ) {
	if ( (*it)->mothers.first != rp && (*it)->mother != rp ) {
	  cout << "\n* parton at " << (*it)->theParton->oY() << " doesnt have " << rp->theParton->oY()
	       << " as first mother! child-mother broken!" << endl;
	  ok = false;
	}
	if ( (*it)->keep == RealParton::NO ) {
	  cout << "\n* parton at " << rp->theParton->oY() << "has a NO-child at "
	       << (*it)->theParton->oY() << endl;
	  ok = false;
	}
      }
      if ( rp->mothers.first && (rp->mothers.first->children.second.find(rp)
				 == rp->mothers.first->children.second.end()) ) {
	cout << "\n* first mother of " << rp->oY() << " doesnt have this as second child!" << endl;
	ok = false;
      }
      if ( rp->mothers.second && (rp->mothers.second->children.first.find(rp)
				 == rp->mothers.second->children.first.end()) ) {
	cout << "\n* second mother of " << rp->oY() << " doesnt have this as first child!" << endl;
	ok = false;
      }
      if ( rp->mother &&
	   (rp->mother->children.first.find(rp) == rp->mother->children.first.end()) &&
	   (rp->mother->children.second.find(rp) == rp->mother->children.second.end()) ) {
	cout << "\n* single mother of " << rp->oY()
	     << " doesnt have this as first or second child!" << endl;
	ok = false;
      }
      if ( !rp->interacting && !rp->theParton->valence() &&
	   rp->children.first.empty() && rp->children.second.empty() ) {
	cout << "\n* " << rp->oY() << " has no kids, and isnt interacting! where to get p-??\n";
	ok = false;
      }
    }
  }
  if ( (intPT - totalRecoil).pt() > 0.000000001*GeV ) {
    cout << "* sum of intPT != totalRecoil" << endl;
    ok = false;
  }
  if ( (pT - totalRecoil).pt() > 0.000000001*GeV ) {
    cout << "* total pT != totalRecoil" << endl;
    ok = false;
  }
//   if ( valence.size() != 3 ) {
//     cout << "* size of valence is not 3" << endl;
//     ok = false;
//   }
  if ( valPT.pt() > 0.000000001*GeV ) {
    cout << "* valencePT does dot add up to zero!" << endl;
    ok = false;
  }
  if ( vals != int(valence.size()) ) {
    cout << "* there are " << vals << " valence partons, but size of valence is " 
	 << valence.size() << endl;
    ok = false;
  }
  if ( kids != mothers ) {
    cout << "* not same number of kids as mothers!" << endl;
    ok = false;
  }
  if ( abs(valPlus - dplus) > 0.00000001*GeV && totalRecoil.pt() == ZERO ) {
    cout << "* initial state plus not conserved" << endl;
    ok = false;
  }
  if ( dplus == ZERO ) {
    cout << "* no plus!" << endl;
    ok = false;
  }
  if ( isnan(dplus/GeV) || isnan(dminus/GeV) ) {
    cout << "* plus or minus is nan!!" << endl;
    ok = false;
  }
  if ( !ok && Debug::level > 5 ) {
    cout << "*\n*       REAL STATE NOT OK!!!\n*\n";
  }
  if ((!ok || pause) && Debug::level > 5 ) {
    cout << setprecision(10);
    cout << "********* Diagnosis for " << this << " ******* call nr " << calls << " ********" << endl;
    cout << "* number of partons:      " << partons.size() << endl;
    cout << "* size of valence:        " << valence.size() << endl;
    cout << "* number of valence part: " << vals << endl;
    cout << "* number of toCheck:      " << toCheck.size() << endl;
    cout << "* number of interacting:  " << interacting.size() << endl;
    cout << "* number of int rps:      " << interactingPartons << endl;
    cout << "* number of cInteracting: " << consInt << endl;
    cout << "* number of interactions: " << interactions.size() << endl;
    cout << "* number gluon exchanges: " << exchanges << endl;
    cout << "* consistent: " << consistent << endl;
    cout << "* total no of kids  = " << kids << endl;
    cout << "* total no mothers  = " << mothers << endl;
    cout << "* number of YES or INT:   " << yes << endl 
	 << "* " << endl;
    cout << "* total transverse  = (" << pT.x()/GeV << ", " << pT.y()/GeV << ")" << endl;
    cout << "* total recoil      = (" << totalRecoil.x()/GeV << ", " 
	 << totalRecoil.y()/GeV << ")" << endl;
    cout << "* sum of int recoil = (" << intPT.x()/GeV << ", " 
	 << intPT.y()/GeV << ")" << endl;
    cout << "* original plus     = " << plus/GeV << endl;
    cout << "* total plus        = " << dplus/GeV << endl;
    cout << "* 1800 - total plus + " << 1800.0 - dplus/GeV << endl;
    cout << "* original minus    = " << minus/GeV << endl;
    cout << "* total minus       = " << dminus/GeV << endl;
    cout << "* given +val minus  = " << (gminus + valMinus)/GeV << endl;
    cout << "* valence minus     = " << valMinus/GeV << endl;
    cout << "* total ctransverse = (" << cpT.x()/GeV << ", " << cpT.y()/GeV << ")" << endl;
    cout << "* total cplus       = " << cplus/GeV << endl;
    cout << "* total cminus      = " << cminus/GeV << endl;
    cout << "* total valencePT   = (" << valPT.x()/GeV << ", " << valPT.y()/GeV << ")" << endl;
    cout << "* total valencePlus = " << valPlus/GeV << endl;
    cout << "* number of flucts  = " << flucts.size() << endl;

    for (   map<tPartonPtr,RealPartonPtr>::const_iterator it = partons.begin();
	    it != partons.end(); it++) {
      tRealPartonPtr rp = (*it).second;
      if ( rp->fluct != -1 ) {
	cout << "*  " << rp->oY() << " belongs to fluct no " << rp->fluct << endl;
      }
//       cout << "* rp at " << rp->theParton->oY() << " has vetoed p+1: " << rp->plus1veto
// 	   << ", p+2: " << rp->plus2veto << ", p-1: " << rp->minus1veto << ", p-2: " << rp->minus2veto
// 	   << ", DGLAP: " << rp->unDGLAP << endl;
//       if ( rp->theParton->valence() )
// 	cout << "* valencepT: (" << rp->theParton->valencePT().x()/GeV << ", "
// 	     << rp->theParton->valencePT().y()/GeV 
// 	     << "), valenceplus: " << rp->theParton->valencePlus()/GeV <<  endl;
    }
    cout << "*********************************************************" << endl;
    plotState(pause || !ok);
  }
  if ( Debug::level > 5 ) cout << "real state tested. ok: " << ok << endl;
  return ok;
}

