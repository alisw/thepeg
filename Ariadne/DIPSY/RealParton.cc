// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RealParton class.
//

#include "Dipole.h"
#include "RealParton.h"
#include "RealPartonState.h"
#include "DipoleEventHandler.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Utilities/Current.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/Debug.h"
#include "ThePEG/Utilities/DebugItem.h"
#include "ThePEG/Utilities/Throw.h"

#include <iostream>


using namespace DIPSY;

RealParton::RealParton(tPartonPtr p):
  keep(NO), cKeep(NO), interacting(DipolePtr()), cInteracting(DipolePtr()), firstInt(false),
  secondInt(false), cFirstInt(false), cSecondInt(false), theParton(p),
  interactionRecoil(make_pair(TransverseMomentum(), ZERO)), fluct(-1), cFluct(-1), cValence(false),
  givenMinus(ZERO), DGLAPchecked(false), intRecoil(TransverseMomentum()), cIntRecoil(TransverseMomentum()) {}

RealParton::~RealParton() {}

void RealParton::setFirstOMother(tRealPartonPtr rp) {
  if (!(oMothers.first) ) {
    oMothers.first = rp;
    if ( rp )
      rp->oChildren.second.insert(this);
  }
}

void RealParton::setSecondOMother(tRealPartonPtr rp) {
  if (!(oMothers.second) ) {
    oMothers.second = rp;
    if ( rp )
      rp->oChildren.first.insert(this);
  }
}

void RealParton::setOMother(tRealPartonPtr rp) {
  if (!oMother ) {
    oMother = rp;
    if ( rp && theParton->parents().first == rp->theParton ) {
      oMothers.first = rp;
      oMothers.second = rp->oMothers.second;
      rp->oChildren.second.insert(this);
    }
    else {
      oMothers.second = rp;
      oMothers.first = rp->oMothers.first;
      rp->oChildren.first.insert(this);
    }
  }
}

void RealParton::addInteraction(tDipolePtr intDip) {
  interactions.insert(intDip);
  plus1veto = 0;
  plus2veto = 0;
  minus1veto = 0;
  minus2veto = 0;
  unDGLAP = 0;
  if ( oMothers.first ) oMothers.first->addInteraction(intDip);
  if ( oMothers.second ) oMothers.second->addInteraction(intDip);
  if ( oMother ) oMother->addInteraction(intDip);
}

void RealParton::eraseMothers() {
  tPartonPtr p = theParton;
  if ( mothers.second ) {
    mothers.second->eraseFirstChild(this);
  }
  if ( mothers.first ) {
    mothers.first->eraseSecondChild(this);
  }
  if ( mother ) {
    mother->eraseFirstChild(this);
    mother->eraseSecondChild(this);
  }
  y = theParton->oY();
}

void RealParton::eraseChildren() {
  for ( RealPartonSet::const_iterator it = children.first.begin();
	it != children.first.end(); it++ )
    eraseFirstChild(*it);
  for ( RealPartonSet::const_iterator it = children.second.begin();
	it != children.second.end(); it++ )
    eraseSecondChild(*it);
}

RealPartonPtr RealParton::firstColourNeighbor() {
  //if children on the other side, they will have the colour flow
  if ( !children.first.empty() )
    return *(--children.first.end());
  //if only one mother, that will have the color flow
  else if ( nMothers == 1 && mother )
    return mother;
  else
    return mothers.first;
  //this will not find valence neighbors. but maybe that's ok?
}

RealPartonPtr RealParton::secondColourNeighbor() {
  if ( !children.second.empty() )
    return *(--children.second.end());
  else if ( nMothers == 1 && mother )
    return mother;
  else
    return mothers.second;
  //this will not find valence dipoles as effective partons. but maybe that's good?
}

RealParton::RealPartonSet RealParton::effectiveParton(InvEnergy range, bool firstSide) {
  RealPartonSet ret;
  ret.insert(this);
  if ( range > Current<DipoleEventHandler>()->coherenceRange() )
    range = Current<DipoleEventHandler>()->coherenceRange();
  // if ( range == ZERO ) return ret;
  
  if ( firstSide ) {
    RealPartonPtr last = secondColourNeighbor();
    for ( RealPartonPtr rp = firstColourNeighbor();
	  rp && rp->theParton->dist2(*theParton) < sqr(0.99*range);
	  rp = rp->firstColourNeighbor() ) {
      if ( rp != this && ret.find(rp) != ret.end() ) {
	break;
      }
      ret.insert(rp);
    }
  }
  else {
    RealPartonPtr last = firstColourNeighbor();
    for ( RealPartonPtr rp = secondColourNeighbor();
	  rp && rp->theParton->dist2(*theParton) < sqr(0.99*range);
	  rp = rp->secondColourNeighbor() ) {
      if ( rp != this && ret.find(rp) != ret.end() ) {
	break;
      }
      ret.insert(rp);
    }
  }
  //check for split virtuals
  for ( RealPartonSet::iterator it = ret.begin(); it != ret.end(); it++) {
    if ( (*it)->fluct != -1 )
      ret.insert(realState->flucts[(*it)->fluct].begin(), realState->flucts[(*it)->fluct].end());
  }
  return ret;
}

void RealParton::doRecoil( RealPartonPtr p, Energy recPlus, TransverseMomentum recPT ) {
  DipoleStatePtr state;
  //first catch some cases that no recoiler wants.
  if ( !theParton->valence() && (!mothers.first || !mothers.second) &&
       !Current<DipoleEventHandler>()->eventFiller().singleMother() )
    return;
  if ( !theParton->valence() && mothers.first == mothers.second &&
       !Current<DipoleEventHandler>()->eventFiller().singleMother() ) return;
  if ( theParton->dipoles().first ) state = & theParton->dipoles().first->dipoleState();
  else state = & theParton->dipoles().second->dipoleState();
  //do the recoil
  p->future.insert(this);
  switch (state->handler().eventFiller().recoilScheme()) {
  case 0:
    doXmasRecoil(p, recPlus, recPT);
    break;
  case 1:
    doRainfallRecoil(p, recPlus, recPT);
    break;
  case 2:
    sunshineRecoil(p, recPlus, recPT);
    break;
  default:
    doXmasRecoil(p, recPlus, recPT);
  }
}

void RealParton::doXmasRecoil( RealPartonPtr p, Energy recPlus, TransverseMomentum recPT ) {
  p->pT -= recPT;
  p->plus -= recPlus;
  // if ( Debug::level > 5 )
  //   cout << "  do recoil of " << recPlus/GeV << ", (" << -recPT.x()/GeV << ", "
  // 	 << -recPT.y()/GeV << ") at "
  // 	 << p->oY() << ", (" << p->theParton->position().x()*GeV << ", "
  // 	 << p->theParton->position().y()*GeV << ")\n";
  if ( p->realState != realState ) {
    p->interactionRecoil.first -= recPT;
    p->interactionRecoil.second -= recPlus;
    recoils.push_back(make_pair(p, make_pair(recPT, recPlus)));
  }
  else
    recoils.push_back(make_pair(p, make_pair(recPT, recPlus)));
  p->updateYMinus();
  // if ( Debug::level > 5 )
  //   cout << oY() << " is recoiling to " << p->oY() << " with (" << (recPT).x()/GeV << ", "
  // 	 << (recPT).y()/GeV << ")" << endl;
}

void RealParton::doRainfallRecoil( RealPartonPtr p, Energy recPlus, TransverseMomentum recPT ) {
  cout << "Rainfall recoil scheme not implemented yet, using Xmas instead." << endl;
  doRecoil(p, recPlus, recPT);
}

void RealParton::sunshineRecoil(RealPartonPtr p, Energy recPlus, TransverseMomentum recPT) {
  cerr << "note, sunshine recoil not consistent with fast realMode atm" << endl;
  //valencepartons cant borrow extra p+, so do xmas recoil
  if ( p->theParton->valence() ) {
    doXmasRecoil(p, recPlus, recPT);
    return;
  }

  //find CoM energy.
  DipoleStatePtr state;
  if ( theParton->dipoles().first ) state = & theParton->dipoles().first->dipoleState();
  else state = & theParton->dipoles().second->dipoleState();
  Energy CoM = (state->plus() + state->collidingEnergy())/2.0;

  //calculate the right angle (plus^a/minus^b=const) to update y and p- in.
  if ( p->minus < ZERO || p->plus < ZERO ) {
    cerr << "negative minus or plus!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
    cerr << oY() << " is giving recoil to " << p->oY() << " with plus,minus "
	 << p->plus/GeV << ", " << p->minus/GeV << endl;
    realState->diagnosis(true);
  }
  double a = log(p->minus/CoM)/log(p->plus/CoM);
  double C = pow(CoM/GeV,a - 1.0);

  p->pT -= recPT;
  recoils.push_back(make_pair(p, make_pair(recPT, ZERO)));

  //calculate new plus for mother
  Energy newPlus = pow(C*sqr(p->pT.pt()/GeV),1.0/(1.0 + a))*GeV;
  Energy extraPlus = newPlus - p->plus;
  
  //check that grandmothers have enough p+, otherwise leave it as xmas recoil
  double P1 = p->mothers.first->plus/(p->mothers.first->plus + p->mothers.second->plus);
  double P2 = 1.0 - P1;
  if ( p->mothers.first->plus < P1*extraPlus || p->mothers.second->plus < P2*extraPlus ) {
    p->plus -= recPlus;
    recoils.push_back(make_pair(p, make_pair(TransverseMomentum(), recPlus)));
    p->updateYMinus();
    return;
  }
  
  //set the new plus and store the recoils
  p->plus = newPlus - recPlus;
  recoils.push_back(make_pair(p, make_pair(TransverseMomentum(), -extraPlus + recPlus)));
  p->updateYMinus();
  p->mothers.first->plus -= P1*extraPlus;
  recoils.push_back(make_pair(p->mothers.first, make_pair(TransverseMomentum(), P1*extraPlus)));
  p->mothers.first->updateYMinus();
  p->mothers.second->plus -= P2*extraPlus;
  recoils.push_back(make_pair(p->mothers.second, make_pair(TransverseMomentum(), P2*extraPlus)));
  p->mothers.second->updateYMinus();
  // realState->plotState(true);
}

void RealParton::doEffectiveRecoil( RealPartonPtr p, InvEnergy range, bool firstSide,
				    Energy recPlus, TransverseMomentum recPT ) {
  RealPartonSet partons = p->effectiveParton(range, firstSide);

  //find the ratio each individual parton should take of the recoil.
  list<double> ratios = effectiveWeights(partons);

  //do the recoils.
  list<double>::iterator ratio = ratios.begin();
  for ( RealPartonSet::iterator it = partons.begin();
	it != partons.end(); it++, ratio++ ) {
    if ( Current<DipoleEventHandler>()->eventFiller().effectiveWeights() == 0 ) {
      if ( *ratio != 0.0 )
	doRecoil(*it, recPlus*(*ratio), recPT*(*ratio));
    }
    else if ( Current<DipoleEventHandler>()->eventFiller().effectiveWeights() == 1 ) {
      if ( *ratio != 0.0 )
	doRecoil(*it, recPlus*(*ratio), recPT/double(partons.size()));
    }
    else if ( Current<DipoleEventHandler>()->eventFiller().effectiveWeights() == 2 ) {
      if ( *ratio != 0.0 )
	doRecoil(*it, recPlus*(*ratio), ((p == *it) ? recPT:TransverseMomentum()));
    }
    else cerr << "no effective weights specified in RealParton" << endl;
  }
}

void RealParton::doPlusWeightedRecoil( RealPartonPtr p, InvEnergy range, bool firstSide,
				    Energy recPlus, TransverseMomentum recPT ) {
  RealPartonSet partons = p->effectiveParton(range, firstSide);
  Energy totalPlus = ZERO;
  //find the total plus in the effective parton to be recoiled.
  for ( RealPartonSet::const_iterator it = partons.begin(); it != partons.end(); it++ )
    totalPlus += (*it)->plus;

  //find the ratio each individual parton should take of the recoil.
  list<double> ratios;
  for ( RealPartonSet::const_iterator it = partons.begin(); it != partons.end(); it++ )
    ratios.push_back((*it)->plus/totalPlus);

  //do the recoils, using Xmas recoil scheme.
  list<double>::iterator ratio = ratios.begin();
  for ( RealPartonSet::iterator it = partons.begin(); it != partons.end(); it++, ratio++ ) {
    Energy diff = -(*it)->minus;
    doXmasRecoil(*it, recPlus*(*ratio), recPT*(*ratio));
    diff += (*it)->minus;
  }
}

void RealParton::undoRecoils() {
  if ( Debug::level > 5 ) cout << "undoing " << recoils.size() << " recoils." << endl;
  for ( list< Recoil >::iterator it = recoils.begin(); it != recoils.end(); it++ ) {
    if ( it->first->realState != realState ) {
      continue;
    }
    //save interaction recoils in case this parton is swithced on again.
    it->first->pT += it->second.first;
    it->first->plus += it->second.second;
    it->first->updateYMinus();
    it->first->future.erase(this);
    if ( it->first->fluct != -1 ) {
      if ( it->first->fluct >= int(realState->flucts.size()) || it->first->fluct < -1 ) {
	cerr << "fluctuation index out of bounds at " << oY() << "!" << endl;
	if ( Debug::level > 5 ) realState->plotState(true);
      }
      // if ( it->first->fluct < -1 ) realState->diagnosis(true); 
      realState->makeCollinear(realState->flucts[it->first->fluct]);
    }
  }
  recoils.clear();
}

void RealParton::undoInteractionRecoils() {
  if ( Debug::level > 5 ) cout << "undointrecs called at " << oY() << endl;
  for ( list< Recoil >::iterator it = recoils.begin(); it != recoils.end(); it++ ) {
    if ( it->first->realState == realState ) continue;
    if ( Debug::level > 5 ) cout << "removing interaction recoil (" << it->second.first.x()/GeV
			     << ", " << it->second.first.y()/GeV << ") at " << it->first->oY() << endl;
    it->first->pT += it->second.first;
    it->first->plus += it->second.second;
    it->first->updateYMinus();
    it->first->future.erase(this);
    recoils.erase(it);
    it--;
  }
}

// void RealParton::checkInteractionRecoils() {
//   if ( keep == NO || savedInteractionRecoils.empty() ) return;
//   // This parton has been switchen off and then on again and should
//   // therefore reobtain its interaction recoils.
//   for ( list< Recoil >::iterator it = savedInteractionRecoils.begin();
// 	it != savedInteractionRecoils.end(); it++ ) {
//     pT += it->second.first;
//     plus += it->second.second;
//     updateYMinus();
//     recoils.insert(*it);
//   }
//   savedInteractionRecoils.clear();
// }

bool RealParton::setNO() {
  if ( theParton->valence() ) moveValenceStatus();
  future.clear();
  eraseMothers();
  undoRecoils();
  eraseExchanges();
  keep = NO;
  return false;
}

bool RealParton::setYES() {
  if ( movedValence && !theParton->valence() ) reclaimValenceStatus();
  if ( nMothers == 1 ) setMother();
  else  setMothers();
  keep = YES;
  return doRecoil();
}

bool RealParton::quickSetYES() {
  if ( movedValence && !theParton->valence() ) reclaimValenceStatus();
  setMothers();
  keep = YES;
  return quickDoRecoil();
}

void RealParton::reclaimValenceStatus() {
  if ( theParton->valence() ) return;
  if ( !movedValence )
    cerr << "RealParton::reclaimValenceStatus called without movedValence" << endl;
  movedValence->theParton->valence(false);
  realState->valence.erase(movedValence);
  movedValence->movedValence = tRealPartonPtr();
  movedValence = tRealPartonPtr();
  theParton->valence(true);
  realState->valence.insert(this);
}

void RealParton::eraseExchanges() {
  if ( exchanges.empty() )
    return;
  for ( set<pair<pair<tRealPartonPtr, tRealPartonPtr>,
	  pair<tRealPartonPtr, tRealPartonPtr> > >::iterator it = exchanges.begin();
	it != exchanges.end(); it++ )  {
    RealPartonPtr rp11 = it->first.first;
    RealPartonPtr rp12 = it->first.second;
    RealPartonPtr rp21 = it->second.first;
    RealPartonPtr rp22 = it->second.second;
    DipoleXSec::InteractionRecoil rec = Current<DipoleEventHandler>()->xSecFn().
      recoil(make_pair(rp11->theParton, rp12->theParton),
	     make_pair(rp21->theParton, rp22->theParton), ImpactParameters());
    rp11->pT -= rec.first.first;
    rp12->pT -= rec.first.second;
    rp21->pT -= rec.second.first;
    rp22->pT -= rec.second.second;
    rp11->updateYMinus();
    rp12->updateYMinus();
    rp21->updateYMinus();
    rp22->updateYMinus();
  }
  exchanges.clear();
}

pair<Energy, Energy> RealParton::effectivePlusMinus(InvEnergy range, bool firstSide) {
  RealPartonSet ep = effectiveParton(range, firstSide);
  Energy ePlus = ZERO;
  TransverseMomentum ePT = TransverseMomentum();
  for ( RealPartonSet::iterator it = ep.begin();
	it != ep.end(); it++ ) {
    ePlus += (*it)->plus;
    ePT += (*it)->pT;
  }
  Energy eMinus = sqr(ePT.pt())/ePlus;
  return make_pair(ePlus, eMinus);
}

Energy RealParton::inRangeMinus(InvEnergy range, bool firstSide) {
  RealPartonSet ep = effectiveParton(range, firstSide);
  Energy eMinus = ZERO;
  for ( RealPartonSet::iterator it = ep.begin(); it != ep.end(); it++ ) {
    eMinus += (*it)->minus;
  }
  return eMinus;
}

bool RealParton::searchNegative(InvEnergy range, bool firstSide) {
  RealPartonSet ep = effectiveParton(range, firstSide);
  for ( RealPartonSet::iterator it = ep.begin(); it != ep.end(); it++ ) {
    if ( (*it)->plus < ZERO || (*it)->minus < ZERO ) return true;
  }
  return false;
}

TransverseMomentum RealParton::opT() {
  TransverseMomentum ret = theParton->valencePT();
  if ( mothers.first ) ret += theParton->recoil(mothers.first->theParton);
  if ( mothers.second ) ret += theParton->recoil(mothers.second->theParton);
  if ( nMothers == 1 && mother )
    ret = theParton->valencePT() + theParton->recoil(mother->theParton);
  if ( ret.pt() == ZERO ) ret = theParton->pT();
  return ret;
}

bool RealParton::isOrdered() {
  if ( nMothers == 1 ) return isSingleOrdered();

  double PSInf = Current<DipoleEventHandler>()->emitter().PSInflation();
  double PMOrd = Current<DipoleEventHandler>()->emitter().PMinusOrdering();
  if ( theParton->valence() ) {
    return true;
  }
  double P1 = theParton->dist2(*mothers.second->theParton)/
    (theParton->dist2(*mothers.first->theParton)+theParton->dist2(*mothers.second->theParton));
  double P2 = 1 - P1;

  if ( y < max(mothers.first->y, mothers.second->y) ) return false;
  else if ( mothers.first->plus < ZERO || mothers.second->plus < ZERO ) return false;
  // else return true;

  if ( fluct == -1 || fluct != mothers.first->fluct ) {
    InvEnergy range1 = sqrt(min(theParton->dist2(*mothers.first->theParton),
				mothers.second->theParton->dist2(*mothers.first->theParton)/4.0));
    if ( Current<DipoleEventHandler>()->emitter().rangeMode() == 1 )
      range1 = sqrt(mothers.second->theParton->dist2(*mothers.first->theParton)/4.0);

    pair<Energy, Energy> pm1 = mothers.first->effectivePlusMinus(range1, true);

    if ( Current<DipoleEventHandler>()->emitter().minusOrderingMode() == 1 ) {
      pm1.second = sqrt(pm1.second*mothers.first->minus);
    }

    if ( Current<DipoleEventHandler>()->emitter().bothOrderedFS() ) {
      P1 = 1.0;   //leaves the ordering unweighted
      P2 = 1.0;
    }
    
    if ( true ) { //Current<DipoleEventHandler>()->emitter().recoilerVeto() ) {
      TransverseMomentum opT1 = mothers.first->opT();
      if ( (pT - opT1).pt() > 2.0*opT1.pt() )  {
    	P1 *= pT.pt()/(opT1.pt());
      }
    }

    if ( P1*plus/PSInf > pm1.first )  {
      if ( Debug::level > 5 ) cout << oY() << " is plus unordered to " << mothers.first->oY() << endl;
      plus1veto++;
      return false;}
    if ( minus*PSInf/P1 < pm1.second*PMOrd ) {
      if ( Debug::level > 5 ) cout << oY() << " is minus unordered to " << mothers.first->oY() << endl;
      minus1veto++;
      return false;}
  }

  if ( fluct == -1 || fluct != mothers.second->fluct ) {
    InvEnergy range2 = sqrt(min(theParton->dist2(*mothers.second->theParton),
				mothers.second->theParton->dist2(*mothers.first->theParton)/4.0));
    if ( Current<DipoleEventHandler>()->emitter().rangeMode() == 1 )
      range2 = sqrt(mothers.second->theParton->dist2(*mothers.first->theParton)/4.0);

    pair<Energy, Energy> pm2 = mothers.second->effectivePlusMinus(range2, false);

    if ( Current<DipoleEventHandler>()->emitter().minusOrderingMode() == 1 ) {
      pm2.second = sqrt(pm2.second*mothers.second->minus);
    }

    if ( true ) { //Current<DipoleEventHandler>()->emitter().recoilerVeto() ) {
      TransverseMomentum opT2 = mothers.second->opT();
      if ( pT.pt() > 2.0*opT2.pt() )  {
    	P2 *= pT.pt()/(opT2.pt());
      }
    }

    if ( P2*plus/PSInf > pm2.first )  {
      if ( Debug::level > 5 ) cout << oY() << " is plus unordered to " << mothers.second->oY() << endl;
      plus2veto++;
      return false;}
    if ( minus*PSInf/P2 < pm2.second*PMOrd ) {
      if ( Debug::level > 5 ) cout << oY() << " is minus unordered to " << mothers.second->oY() << endl;
      minus2veto++;
      return false;}
  }

  //check that the recoil didnt push any parton into negative p+
  for ( list< pair<tRealPartonPtr, pair<TransverseMomentum, Energy> > >::iterator it =
	  recoils.begin(); it != recoils.end(); it++ ) {
    if ( it->first->plus < ZERO ) return false;
  }
  return true;
}

bool RealParton::isFirstOrdered() {
  double PMOrd = Current<DipoleEventHandler>()->emitter().PMinusOrdering();
  double PSInf = Current<DipoleEventHandler>()->emitter().PSInflation();
  if ( realState->monitored && realState->interactions.size() > 1) realState->plotState(false);
  if ( theParton->valence() ) {
    return true;
  }

  double P1 = theParton->dist2(*mothers.second->theParton)/
    (theParton->dist2(*mothers.first->theParton)+theParton->dist2(*mothers.second->theParton));
  InvEnergy range1 = sqrt(min(theParton->dist2(*mothers.first->theParton),
			      mothers.second->theParton->dist2(*mothers.first->theParton)/4.0));
  pair<Energy, Energy> pm1 = mothers.first->effectivePlusMinus(range1, true);
  if ( Current<DipoleEventHandler>()->emitter().bothOrderedFS() ) {
    P1 = 1.0;   //leaves the ordering unweighted
  }

  if ( P1*plus/PSInf > pm1.first ) return false;
  if ( minus*PSInf/P1 < pm1.second*PMOrd ) return false;

  //check that the recoil didnt push any parton into negative p+
  for ( list< pair<tRealPartonPtr, pair<TransverseMomentum, Energy> > >::iterator it =
	  recoils.begin(); it != recoils.end(); it++ ) {
    if ( it->first->plus < ZERO ) return false;
  }
  return true;
}

bool RealParton::isSecondOrdered() {
  double PMOrd = Current<DipoleEventHandler>()->emitter().PMinusOrdering();
  double PSInf = Current<DipoleEventHandler>()->emitter().PSInflation();
  if ( realState->monitored && realState->interactions.size() > 1) realState->plotState(false);
  if ( theParton->valence() ) {
    return true;
  }
  double P1 = theParton->dist2(*mothers.second->theParton)/
    (theParton->dist2(*mothers.first->theParton)+theParton->dist2(*mothers.second->theParton));
  double P2 = 1 - P1;
  if ( Current<DipoleEventHandler>()->emitter().bothOrderedFS() )
    P2 = 1.0;

  InvEnergy range2 = sqrt(min(theParton->dist2(*mothers.second->theParton),
			      mothers.second->theParton->dist2(*mothers.first->theParton)/4.0));
  pair<Energy, Energy> pm2 = mothers.second->effectivePlusMinus(range2, false);

  if ( P2*plus/PSInf > pm2.first ) return false;
  if ( minus*PSInf/P2 < pm2.second*PMOrd ) return false;

  //check that the recoil didnt push any parton into negative p+
  for ( list< pair<tRealPartonPtr, pair<TransverseMomentum, Energy> > >::iterator it =
	  recoils.begin(); it != recoils.end(); it++ ) {
    if ( it->first->plus < ZERO ) return false;
  }
  return true;
}

bool RealParton::isSingleOrdered() {
  if ( fluct != -1 && fluct == mother->fluct ) return true;
  double PMOrd = Current<DipoleEventHandler>()->emitter().PMinusOrdering();
  double PSInf = Current<DipoleEventHandler>()->emitter().PSInflation();
  if ( theParton->valence() ) return true;
  InvEnergy range = sqrt(theParton->dist2(*mother->theParton));
  pair<Energy, Energy> pm = mother->effectivePlusMinus(range, true);

  if ( y < mother->y ) return false;
  else if ( mother->plus < ZERO ) return false;
  else return true;

  if ( Current<DipoleEventHandler>()->emitter().minusOrderingMode() == 1 ) {
    pm.second = sqrt(pm.second*mother->minus);
  }

  if ( plus/PSInf > pm.first )  {
    if ( Debug::level > 5 ) cout << oY() << " is plus unordered to " << mother->oY() << endl;
    plus1veto++;
    return false;
  }
  if ( minus*PSInf < pm.second*PMOrd ) {
    if ( Debug::level > 5 ) cout << oY() << " is minus unordered to " << mother->oY() << endl;
    minus1veto++;
    return false;
  }
  //check that the recoil didnt push any parton into negative p+
  for ( list< pair<tRealPartonPtr, pair<TransverseMomentum, Energy> > >::iterator it =
	  recoils.begin(); it != recoils.end(); it++ )
    if ( it->first->plus < ZERO ) return false;
  return true;
}

bool RealParton::hasEnergy() {
  for ( list< pair<tRealPartonPtr, pair<TransverseMomentum, Energy> > >::iterator it =
	  recoils.begin(); it != recoils.end(); it++ )
    if ( it->first->plus < ZERO ) {
      return false;
    }
  return true;
}

Energy RealParton::problemScale() {
  return max(max(orderedScale(),DGLAPScale()),motherScale());
}

Energy RealParton::orderedScale() {
  return emissionScale;
}

Energy RealParton::DGLAPScale() {
  if ( DGLAPSafe() ) return ZERO;
  if ( children.first.empty() ) return 1.0/sqrt(theParton->dist2(*mothers.first->theParton));
  if ( children.second.empty() ) return 1.0/sqrt(theParton->dist2(*mothers.second->theParton));
  cerr << "non DGLAP safe with children on both sides in DGLAPScale!!!!!" << endl;
  return ZERO;
}

Energy RealParton::motherScale() {
  Energy ret = ZERO;
  if ( !theParton->valence() ) {
    Energy firstOMother = 1.0/sqrt(theParton->dist2(*oMothers.first->theParton));
    Energy secondOMother = 1.0/sqrt(theParton->dist2(*oMothers.second->theParton));
    if ( !mothers.first )
      ret = max(ret, firstOMother);
    if ( !mothers.second )
      ret = max(ret, secondOMother);
    if ( mothers.first && mothers.first == mothers.second )
      ret = max(ret, 1.0/sqrt(theParton->dist2(*mothers.first->theParton)));
    if ( !interacting && children.first.empty() && children.second.empty() )
      ret = max(ret, pT.pt()); //since this is a virtual fluct, the pT of the parton sets the scale
  }
  return ret;
}

RealPartonPtr RealParton::findCause() {
  Energy oScale = orderedScale();
  Energy DScale = DGLAPScale();
  Energy mScale = motherScale();
  RealPartonPtr cause;
  if ( oScale > max(DScale, mScale) ) {
    // cout << "  unordered" << endl;
    cause = orderedCause();
  }
  else if ( DScale > mScale ) {
    // cout << "  dglap is the cause" << endl;
    cause = DGLAPCause();
  }
  else if ( mScale > ZERO ) {
    // cout << "  bad mother structure" << endl;
    cause = motherCause();
  }
  else cerr << "findCause() called without any errors... >_>" << endl;
  if ( !cause ) {
    cerr << "all error scales ZERO for primary cause" << endl;
    realState->plotState(true);
    return cause;
  }
  // cout << "  cause at " << cause->oY() << ", ("
  // 	   << cause->theParton->position().x()*GeV << ", "
  // 	   << cause->theParton->position().y()*GeV << ")" << endl;

  //now try to handle if an interacting parton is the cause.
  while ( cause->interacting ) {
    // cout << "  interacting cause, move" << endl;
    if ( cause->theParton->valence() ) { //trying to remove interacting valence --> screwed
      if ( !interacting ) return this;
      // cout << "interacting valence parton is the cause, screwed." << endl;
      if ( oScale > max(DScale, mScale) )
	cout << "  failed evo due to unordered\n";
      else if ( DScale > mScale )
	cout << "  failed evo due to DGLAP\n";
      else
	cout << "  failed evo due to mother structure\n";
      realState->plotState(true);
      return RealPartonPtr();
    }
    if ( !cause->mothers.first && !cause->mothers.second ) {
      // cout << "want to remove interacting nonvalence wihtout mothers at " << cause->oY() << endl;
      // realState->plotState(true);
      return RealPartonPtr();
    }
    if ( !cause->mothers.first ) cause = cause->mothers.second; //take only existing mother
    else if ( !cause->mothers.second ) cause = cause->mothers.first;
    else if ( cause->mothers.first->interacting && !cause->mothers.second->interacting )
      cause = cause->mothers.second; //otherwise take non-interacting mother
    else if ( !cause->mothers.first->interacting && cause->mothers.second->interacting )
      cause = cause->mothers.first;
    else if ( cause->theParton->dist2(*cause->mothers.first->theParton) <
	      cause->theParton->dist2(*cause->mothers.second->theParton) )
      cause = cause->mothers.first; //otherwise take closest mother
    else
      cause = cause->mothers.second;
  }
  // cout << "  final cause at " << cause->oY() << ", ("
  // 	   << cause->theParton->position().x()*GeV << ", "
  // 	   << cause->theParton->position().y()*GeV << ")" << endl;
  if ( cause &&  realState->toCheck.find(cause) == realState->toCheck.end() )
    cout << "cause is outside toCheck!!" << endl;
  return cause;
}

void RealParton::checkEmissionProblem() {
  double PSInf = Current<DipoleEventHandler>()->emitter().PSInflation();
  emissionScale = ZERO;
  emissionCause = RealPartonPtr();
  if ( theParton->valence() || !mothers.first || !mothers.second ) { //what to do if 1 mother?
    return;
  }
  double P1 = theParton->dist2(*mothers.second->theParton)/
    (theParton->dist2(*mothers.first->theParton)+theParton->dist2(*mothers.second->theParton));
  double P2 = 1 - P1;
  InvEnergy range1 = sqrt(min(theParton->dist2(*mothers.first->theParton),
			      mothers.second->theParton->dist2(*mothers.first->theParton)/4.0));
  pair<Energy, Energy> pm1 = mothers.first->effectivePlusMinus(range1, true);
  InvEnergy range2 = sqrt(min(theParton->dist2(*mothers.second->theParton),
			      mothers.second->theParton->dist2(*mothers.first->theParton)/4.0));
  pair<Energy, Energy> pm2 = mothers.second->effectivePlusMinus(range2, false);
  if ( P1*plus/PSInf > pm1.first && pT.pt() > emissionScale )  {
    emissionScale = pT.pt();
    emissionCause = this;
  }
  if ( minus*PSInf/P1 < pm1.second && mothers.first->pT.pt() > emissionScale ) {
    emissionScale = mothers.first->pT.pt();
    emissionCause = mothers.first;
  }
  if ( P2*plus/PSInf > pm2.first && pT.pt() > emissionScale )  {
    emissionScale = pT.pt();
    emissionCause = this;
  }
  if ( minus*PSInf/P2 < pm2.second && mothers.second->pT.pt() > emissionScale )  {
    emissionScale = mothers.second->pT.pt();
    emissionCause = mothers.second;
  }
  if ( !hasEnergy() ) {
    emissionScale = pT.pt();
    emissionCause = this;
  }
}

RealPartonPtr RealParton::orderedCause() {
  return emissionCause;
}

RealPartonPtr RealParton::DGLAPCause() {
  if ( DGLAPSafe() ) cerr << "asked for DGLAPCause from dglap safe parton." << endl;
  return this;
  // if ( children.first.empty() ) return mothers.first;
  // if ( children.second.empty() ) return mothers.second;
  cerr << "non DGLAP safe with children on both sides in DGLAPScale!!!!!" << endl;
  return RealPartonPtr();
}

RealPartonPtr RealParton::motherCause() {
  return this;
}

RealPartonPtr RealParton::findSuspect() {
  if ( theParton->valence() ) {
    cout << "valence, couldnt find suspect." << endl;
    return RealPartonPtr();
  }
  if ( !mothers.first || !mothers.second ) {
    cout << "  no mother at " << oY() << endl;
    realState->plotState(true);
  }
  if ( mothers.first && !mothers.second ) return mothers.first;
  if ( !mothers.first && mothers.second ) return mothers.second;
  if ( !mothers.first && !mothers.second ) return RealPartonPtr();
  double PSInf = Current<DipoleEventHandler>()->emitter().PSInflation();
  double P1 = theParton->dist2(*mothers.second->theParton)/
    (theParton->dist2(*mothers.first->theParton)+theParton->dist2(*mothers.second->theParton));
  double P2 = 1 - P1;
  InvEnergy range1 = sqrt(min(theParton->dist2(*mothers.first->theParton),
			      mothers.second->theParton->dist2(*mothers.first->theParton)/4.0));
  pair<Energy, Energy> pm1 = mothers.first->effectivePlusMinus(range1, true);
  InvEnergy range2 = sqrt(min(theParton->dist2(*mothers.second->theParton),
			      mothers.second->theParton->dist2(*mothers.first->theParton)/4.0));
  pair<Energy, Energy> pm2 = mothers.second->effectivePlusMinus(range2, false);
  bool firstSuspect = false;
  bool secondSuspect = false;
  if ( P1*plus/PSInf > pm1.first || minus*PSInf/P1 < pm1.second ) firstSuspect = true;
  if ( P2*plus/PSInf > pm2.first || minus*PSInf/P2 < pm2.second ) secondSuspect = true;

  cout << "first, second suspect: " << firstSuspect << ", " << secondSuspect << endl;

  if ( firstSuspect && !secondSuspect ) return mothers.first;
  else if ( !firstSuspect && secondSuspect ) return mothers.second;
  else if ( mothers.first->oY() > mothers.second->oY() ) return mothers.first;
  else return mothers.second;
  //switch of the youngest of the mothers when in doubt.
}

void RealParton::setValenceMomentum() {
  if ( !theParton->valence() ) return;
  pT = theParton->valencePT();
  plus = theParton->valencePlus();
  y = log(pT.pt()/plus);
  minus = pT.pt()*exp(y);
  // cout << "  set valencemomentum (" << pT.x()/GeV << ", " << pT.y()/GeV << ") at "
  //      << oY() << ", (" << theParton->position().x()*GeV << ", "
  //      << theParton->position().y()*GeV << ")\n";
}

bool RealParton::doRecoil() {
  if ( theParton->valence() ) {
    setValenceMomentum();
    return true;
  }
  recoils.clear();
  if ( nMothers == 1 )
    return doSingleRecoil();
  tRealPartonPtr m1 = mothers.first;
  tRealPartonPtr m2 = mothers.second;
  bool ok = true;
  if ( !m1 && !m2 ) return false;
  tPartonPtr p = theParton;
  pT = TransverseMomentum(ZERO, ZERO);
  y = p->oY();
  if ( m1 && !m2 ) {
    tPartonPtr p1 = m1->theParton;
    pT += p->recoil(p1);
    plus = pT.pt()*exp(-p->oY());
    minus = pT.pt()*exp(p->oY());
    // cout << "  emit1 with pt (" << pT.x()/GeV << ", " << pT.y()/GeV << ") at "
    // 	 << oY() << ", (" << theParton->position().x()*GeV << ", "
    // 	 << theParton->position().y()*GeV << ")\n";
    ok = false;
    InvEnergy range1 = sqrt(p->dist2(*p1));
    doEffectiveRecoil(m1, range1, true, plus, p->recoil(p1));
  }
  else if ( !m1 && m2 ) {
    tPartonPtr p2 = m2->theParton;
    ok = false;
    pT += p->recoil(p2);
    plus = pT.pt()*exp(-p->oY());
    minus = pT.pt()*exp(p->oY());
    // cout << "  emit2 with pt (" << pT.x()/GeV << ", " << pT.y()/GeV << ") at "
    // 	 << oY() << ", (" << theParton->position().x()*GeV << ", "
    // 	 << theParton->position().y()*GeV << ")\n";
    InvEnergy range2 = sqrt(p->dist2(*p2));
    doEffectiveRecoil(m2, range2, false, plus, p->recoil(p2));
  }
  else {
    tPartonPtr p1 = m1->theParton;
    tPartonPtr p2 = m2->theParton;
    pT += p->recoil(p1) + p->recoil(p2);
    plus = pT.pt()*exp(-p->oY());
    minus = pT.pt()*exp(p->oY());

    // cout << "  emit3 with pt (" << pT.x()/GeV << ", " << pT.y()/GeV << ") at "
    // 	 << oY() << ", (" << theParton->position().x()*GeV << ", "
    // 	 << theParton->position().y()*GeV << ")\n";

    //potentially put this one back eventually
    // doSwingRecoil();

//     checkPS();
    double P1 = p->dist2(*m2->theParton)/
      (p->dist2(*m1->theParton) + p->dist2(*m2->theParton));
    double P2 = 1.0 - P1;
    InvEnergy range1 = sqrt(min(p->dist2(*p1), p2->dist2(*p1)/4.0));
    if ( Current<DipoleEventHandler>()->emitter().rangeMode() == 1 )
      range1 = sqrt(p2->dist2(*p1)/4.0);
    doEffectiveRecoil(m1, range1, true, plus*P1, p->recoil(p1));
    InvEnergy range2 = sqrt(min(p->dist2(*p2), p1->dist2(*p2)/4.0));
    if ( Current<DipoleEventHandler>()->emitter().rangeMode() == 1 )
      range2 = sqrt(p1->dist2(*p2)/4.0);
    doEffectiveRecoil(m2, range2, false, plus*P2, p->recoil(p2));
    ok = isOrdered(); //check ordering  recoil
  }
  return ok;
}

bool RealParton::doSingleRecoil() {
  tPartonPtr p = theParton;
  y = p->oY();
  pT = p->recoil(mother->theParton);
  plus = pT.pt()*exp(-p->oY());
  minus = pT.pt()*exp(p->oY());
  InvEnergy range = sqrt(p->dist2(*mother->theParton));
  doEffectiveRecoil(mother, range, true, plus, p->recoil(mother->theParton));
  return isSingleOrdered();
}

bool RealParton::quickDoRecoil() {
  if ( theParton->valence() ) {
    setValenceMomentum();
    return true;
  }
  recoils.clear();
  tRealPartonPtr m1 = mothers.first;
  tRealPartonPtr m2 = mothers.second;
  bool ok = true;
  if ( !m1 && !m2 ) return false;
  tPartonPtr p = theParton;
  pT = TransverseMomentum(ZERO, ZERO);
  y = p->oY();
  if ( m1 && !m2 ) { //can I let these through really? will rest of evo be ok?
    tPartonPtr p1 = m1->theParton;
    pT += p->recoil(p1);
    // cout << "  pt at creation1 (" << pT.x()/GeV << ", " << pT.y()/GeV << ") at "
    // 	 << oY() << ", (" << theParton->position().x()*GeV << ", "
    // 	 << theParton->position().y()*GeV << ")\n";
    plus = pT.pt()*exp(-p->oY());
    minus = pT.pt()*exp(p->oY());
    InvEnergy range1 = sqrt(p->dist2(*p1));
    doEffectiveRecoil(m1, range1, true, plus, p->recoil(p1));
  }
  else if ( !m1 && m2 ) {
    tPartonPtr p2 = m2->theParton;
    pT += p->recoil(p2);
    // cout << "  pt at creation2 (" << pT.x()/GeV << ", " << pT.y()/GeV << ") at "
    // 	 << oY() << ", (" << theParton->position().x()*GeV << ", "
    // 	 << theParton->position().y()*GeV << ")\n";
    plus = pT.pt()*exp(-p->oY());
    minus = pT.pt()*exp(p->oY());
    InvEnergy range2 = sqrt(p->dist2(*p2));
    doEffectiveRecoil(m2, range2, false, plus, p->recoil(p2));
  }
  else {
    tPartonPtr p1 = m1->theParton;
    tPartonPtr p2 = m2->theParton;
    pT += p->recoil(p1) + p->recoil(p2);
    // cout << "  pt at creation3 (" << pT.x()/GeV << ", " << pT.y()/GeV << ") at "
    // 	 << oY() << ", (" << theParton->position().x()*GeV << ", "
    // 	 << theParton->position().y()*GeV << ")\n";
    plus = pT.pt()*exp(-p->oY());
    minus = pT.pt()*exp(p->oY());

    //potentially put this one back eventually
    // doSwingRecoil();

//     checkPS();
    double P1 = p->dist2(*m2->theParton)/
      (p->dist2(*m1->theParton) + p->dist2(*m2->theParton));
    double P2 = 1.0 - P1;
    InvEnergy range1 = sqrt(min(p->dist2(*p1), p2->dist2(*p1)/4.0));
    doEffectiveRecoil(m1, range1, true, plus*P1, p->recoil(p1));
    InvEnergy range2 = sqrt(min(p->dist2(*p2), p1->dist2(*p2)/4.0));
    doEffectiveRecoil(m2, range2, false, plus*P2, p->recoil(p2));
    ok = hasEnergy(); //check ordering  recoil
  }
  return ok;
}

void RealParton::checkPS() const {
  cout << "checking PS" << endl;
  if ( mothers.first == mothers.second ) {
    cout << "same mothers" << endl;
    return;
  }
  DipolePtr dip = new_ptr(Dipole());
  PartonPtr parton1 = new_ptr(Parton());
  PartonPtr parton2 = new_ptr(Parton());
  parton1->position(mothers.first->theParton->position());
  parton1->y(mothers.first->y);
  parton1->pT(mothers.first->pT);
  parton1->plus(mothers.first->plus);
  parton1->minus(mothers.first->minus);
  parton2->position(mothers.second->theParton->position());
  parton2->y(mothers.second->y);
  parton2->pT(mothers.second->pT);
  parton2->plus(mothers.second->plus);
  parton2->minus(mothers.second->minus);
  dip->partons(make_pair(parton1, parton2));
  DipoleStatePtr teststate = new_ptr(DipoleState());
  teststate->handler(&Current<DipoleEventHandler>::current());
  dip->dipoleState(teststate);
  teststate->collidingEnergy(2000.0*GeV);

  Current<DipoleEventHandler>()->emitter().testingPS = true;
  // double inflation =   Current<DipoleEventHandler>()->emitter().PSInflation();
  // Current<DipoleEventHandler>()->emitter().PSInflation = 1.0; REMOVED DUE TO INTERFACED
  Current<DipoleEventHandler>()->emitter().fixY = oY();
  InvEnergy minr1 = ZERO, minr2 = ZERO, maxr1 = ZERO, maxr2 = ZERO;
  maxr2 = ZERO;
  for ( int i = 0; i < 10000; i++ ) {
    Current<DipoleEventHandler>()->emitter().generate(*dip, oY() - 1.0, oY() + 1.0);
    if ( dip->generatedGluon() ) {
      InvEnergy r1 = (dip->generatedGluon()->position() - dip->partons().first->position()).pt();
      InvEnergy r2 = (dip->generatedGluon()->position() - dip->partons().second->position()).pt();
      if ( maxr2 == ZERO || r1 < minr1 )
	minr1 = r1;
      if ( maxr2 == ZERO || r2 < minr2 )
	minr2 = r2;
      if ( maxr2 == ZERO || r1 > maxr1 )
	maxr1 = r1;
      if ( maxr2 == ZERO || r2 > maxr2 )
	maxr2 = r2;
    }
  }
  cout << "done testing!" << endl;
  Current<DipoleEventHandler>()->emitter().testingPS = false;
  // Current<DipoleEventHandler>()->emitter().PSInflation = inflation; REMOVED DUE TO INTERFACED
  cout << "at oY = " << oY() << ", r1 = " 
       << (theParton->position() - mothers.first->theParton->position()).pt()*GeV << ", r2 = "
       << (theParton->position() - mothers.second->theParton->position()).pt()*GeV << endl;
  cout << "orig emission: m1pT.pt() =  " << theParton->m1pT.pt()/GeV
       << ", m1plus =  " << theParton->m1plus/GeV
       << ", m2pT.pt() =  " << theParton->m2pT.pt()/GeV
       << ", m2plus =  " << theParton->m2plus/GeV << endl;
  cout << "orig emission: minr1 = " << theParton->minr1*GeV << ", minr2 = " << theParton->minr2*GeV
       << ", maxr1 = " << theParton->maxr1*GeV << ", maxr2 = " << theParton->maxr2*GeV << endl;
  cout << "real emission: m1pT.pt() =  " << mothers.first->pT.pt()/GeV
       << ", m1plus =  " << mothers.first->plus/GeV
       << ", m2pT.pt() =  " << mothers.second->pT.pt()/GeV
       << ", m2plus =  " << mothers.second->plus/GeV << endl
       << "real  emission: minr1 = " << minr1*GeV << ", minr2 = " << minr2*GeV
       << ", maxr1 = " << maxr1*GeV << ", maxr2 = " << maxr2*GeV << endl << endl;
  cout << "ratios: min1: " << minr1/theParton->minr1
       << ", min2: " << minr2/theParton->minr2
       << ", max1: " << theParton->maxr1/maxr1
       << ", max2: " << theParton->maxr2/maxr2 << endl;
}

bool RealParton::doSwingRecoil() {
  if ( !mothers.first || !mothers.second )
    cerr << "doswingrecoil in realparton called without mother" << endl;
  RealPartonPtr m1 = mothers.first;
  RealPartonPtr m2 = mothers.second;
  RealPartonPtr gm12 = mothers.first->mothers.second;
  RealPartonPtr gm21 = mothers.second->mothers.first;
  //normal non-swing parent structure.
  if ( !gm12 || !gm21 || gm12 == m2 || gm21 == m1 || m1 == m2 )
    return false;
  //standard swing structure
  if ( gm12 != gm21 ) {
    //maybe do recoil m1,gm12 vs gm21,m2
//     if ( 9.0*Current<DipoleEventHandler>()->xSecFn().
// 	 fij(make_pair(m1->theParton, gm12->theParton),
// 	     make_pair(gm21->theParton, m2->theParton), ImpactParameters())/
// 	 Current<DipoleEventHandler>()->swinger().
// 	 swingAmp(make_pair(m1->theParton, gm12->theParton),
// 		  make_pair(gm21->theParton, m2->theParton) ) > UseRandom::rnd() ) {
//       Current<DipoleEventHandler>()->xSecFn().
// 	recoil(make_pair(m1->theParton, gm12->theParton),
// 	       make_pair(gm21->theParton, m2->theParton), ImpactParameters());
//       cout << "did extra recoil in realparton!!" << endl;
//       return true;
//     }
    if ( 1.0/3.0 > UseRandom::rnd() ) {
      // if ( realState->monitored )
      // 	cout << "doing swing recoil!" << endl;
//       realState->diagnosis(true);
      DipoleXSec::InteractionRecoil rec = Current<DipoleEventHandler>()->xSecFn().
	recoil(make_pair(m1->theParton, gm12->theParton),
	       make_pair(gm21->theParton, m2->theParton), ImpactParameters());
      m1->pT += rec.first.first;
      gm12->pT += rec.first.second;
      gm21->pT += rec.second.first;
      m2->pT += rec.second.second;
      m1->updateYMinus();
      m2->updateYMinus();
      gm12->updateYMinus();
      gm21->updateYMinus();
      exchanges.insert(make_pair(make_pair(m1, gm12), make_pair(gm21, m2)));
//       cout << "done swing recoil!" << endl;
//       realState->diagnosis(true);
      //       cout << "did extra recoil in realparton!!" << endl;
      return true;
    }
    else {
//       cout << "tested and failed recoilswing in realparton" << endl;
      return false;
    }
  }
  //M-emission, must have been 2 swings, but do only last
  else {
//     cout << "M-emission in realparton, cant find the recoil. :(" << endl;
    return false;
  }
  return false;
}

void RealParton::updateYMinus() {
  if ( plus <= ZERO || pT.pt2() <= ZERO ) {
    y = 0.0;
  }
  else {
    y = log(pT.pt()/plus);
    minus = pT.pt2()/plus;
  }
}

void RealParton::eraseFirstChild(tRealPartonPtr rp) {
  children.first.erase(rp);
}

void RealParton::eraseSecondChild(tRealPartonPtr rp) {
  children.second.erase(rp);
}

bool RealParton::singleDGLAPSafe(InvEnergy scale) {
  cout << "entering singledglap" << endl;
  if ( secondInt || firstInt ) 
    return checkDGLAPSafe(mother, RealPartonPtr(), scale);
  if ( children.first.empty() ) return true;
  else return checkDGLAPSafe(mother, *(--children.first.end()), scale);
}

bool RealParton::DGLAPSafe(InvEnergy scale) {
  if ( nMothers == 1 )
    return singleDGLAPSafe(scale);
  bool DGLAPSuppresion = true; //should probably be interfaced. switches dglap suppression off.
  if ( !DGLAPSuppresion ) return true;

  //check for the right structure: no childs or int on one side
  //and a child or an int (or both) on the other. pass on if right structure.
  if ( children.first.empty() && !firstInt ) {
    if ( secondInt ) 
      return checkDGLAPSafe(mothers.first, RealPartonPtr(), scale);
    else if ( !children.second.empty() ) 
      return checkDGLAPSafe(mothers.first, *(--children.second.end()), scale);
    else return true;
  }
  if ( children.second.empty() && !secondInt ) {
    if ( firstInt ) 
      return checkDGLAPSafe(mothers.second, RealPartonPtr(), scale);
    else if ( !children.first.empty() ) 
      return checkDGLAPSafe(mothers.second, *(--children.first.end()), scale);
    else return true;
  }
  return true;
}

bool RealParton::checkDGLAPSafe(tRealPartonPtr mother, tRealPartonPtr child, InvEnergy scale) {
  if ( !mother ) return true;

  //if the small dipole is interacting, the dipole is always safe.
  if ( interacting && interacting == mother->interacting ) return true; //not perfect solution

  //the size of the small (backwards) dipole
  InvEnergy r1 = sqrt(mother->theParton->dist2(*theParton));

  //take the large (forward) size from the interaction, or otherwise last child.
  InvEnergy r2;
  if ( !child ) r2 = intDist;
  else r2 = sqrt(child->theParton->dist2(*theParton));

  //set the forward scale to what is supplied, if something is supplied.
  //this is if there are several steps down in pt after each other. compare to bottom.
  if ( scale != ZERO ) r2 = scale;

  // *** ATTENTION *** this is to avoid divide by zero in alphaS.
  if ( r2 == ZERO ) return true;

  //never suppress above coherence range
  if ( r2 > Current<DipoleEventHandler>()->coherenceRange() )
    r2 = Current<DipoleEventHandler>()->coherenceRange();

  // Calculate fudge factor to emulate ME-corrections
  double fudgeME = Current<DipoleEventHandler>()->fudgeME()?
    1.0 - 1.0/(1.0 + cosh(theParton->y() - mother->theParton->y())): 1.0;

  //check if the parton has been checked with this mother before.
  //if so, reuse the randomised resolution scale r1^2*alphas(r1)*R
  map<tRealPartonPtr, InvEnergy2>::iterator it = resolutionScales.lower_bound(mother);
  if ( (*it).first == mother )
    return (*it).second > sqr(r2)*Current<DipoleEventHandler>()->alphaS(r2);
  //otherwise generate a new resolution scale r1^2*alphas(r1)*R
  resolutionScales.insert(pair<tRealPartonPtr, InvEnergy2>
			  (mother, fudgeME*(sqr(r1)*Current<DipoleEventHandler>()->alphaS(r1))
			   /(UseRandom::rnd())));

  if ( Debug::level > 5 ) {
    cout << "mother scale: " << r1*GeV << ", child scale: " << r2*GeV;
    if ( resolutionScales[mother] > (sqr(r2)*Current<DipoleEventHandler>()->alphaS(r2)) )
				      cout << " safe!" << endl;
    else cout << " NOT safe!!" << endl;
  }

  //return safe if the forwards dipole is small enough to see the backward dipole.
  return resolutionScales[mother] > (sqr(r2)*Current<DipoleEventHandler>()->alphaS(r2));
}

tRealPartonPtr RealParton::youngestFirstChild(tRealPartonPtr rp, bool checkMother) {
  if ( oChildren.first.size() > 0 ) {
    RealPartonSet::iterator it = oChildren.first.lower_bound(rp);
    if ( it != oChildren.first.begin() ) {
      if ( (*(--it))->keep == YES ) return *it;
      else {
	tRealPartonPtr cand = (*it)->youngestFirstChild(rp, false);
	if ( cand ) return cand;
      }
    }
  }
  if ( !checkMother ) return RealPartonPtr();
  if ( !(oMothers.first) ) {
    if ( movedValence && movedValence->oY() < rp->oY() ) {
      // cout << "returning " << movedValence->oY() << ", ("
      // 	   << movedValence->theParton->position().x()*GeV << ", "
      // 	   << movedValence->theParton->position().y()*GeV << ") as first mother to "
      // 	   << rp->oY() << ", (" << rp->theParton->position().x()*GeV << ", "
      // 	   << rp->theParton->position().y()*GeV << ")\n";
      return movedValence;
    }
    return RealPartonPtr();
  }
  if ( oMothers.first->keep == YES ) return oMothers.first;
  else return oMothers.first->youngestFirstChild(rp, true);
}

tRealPartonPtr RealParton::youngestSecondChild(tRealPartonPtr rp, bool checkMother) {
  if ( oChildren.second.size() > 0 ) {
    RealParton::RealPartonSet::iterator it = oChildren.second.lower_bound(rp);
    if ( it != oChildren.second.begin() ) {
      if ( (*(--it))->keep == YES ) return *it;
      else {
	tRealPartonPtr cand = (*it)->youngestSecondChild(rp, false);
	if ( cand ) return cand;
      }
    }
  }
  if ( !checkMother ) return RealPartonPtr();
  if ( !(oMothers.second) ) {
    if ( movedValence && movedValence->oY() < rp->oY() ) {
      // cout << "returning " << movedValence->oY() << ", ("
      // 	   << movedValence->theParton->position().x()*GeV << ", "
      // 	   << movedValence->theParton->position().y()*GeV << ") as first mother to "
      // 	   << rp->oY() << ", (" << rp->theParton->position().x()*GeV << ", "
      // 	   << rp->theParton->position().y()*GeV << ")\n";
      return movedValence;
    }
    return RealPartonPtr();
  }
  if ( oMothers.second->keep == YES ) return oMothers.second;
  else return oMothers.second->youngestSecondChild(rp, true);
}

void RealParton::setMothers() {
  if ( nMothers == 1 ) {
    setMother();
    return;
  }
  tRealPartonPtr first;
  tRealPartonPtr second;
  if ( theParton->valence() ) {
    mothers = pair<tRealPartonPtr, tRealPartonPtr>(first, second);
    return;
  }
  if ( oMothers.first->keep == YES ) first = oMothers.first;
  else first = oMothers.first->youngestFirstChild(this, true);
  if ( oMothers.second->keep == YES ) second = oMothers.second;
  else second = oMothers.second->youngestSecondChild(this, true);
  mothers = pair<tRealPartonPtr, tRealPartonPtr>(first, second);
  if ( first ) first->children.second.insert(this);
  if ( second ) second->children.first.insert(this);
}

void RealParton::setMother() {
  if ( theParton->valence() ) {
    mother = RealPartonPtr();
    return;
  }
  if ( !oMother ) cerr << "non-valence parton without mother at y0 = " << oY() << "! :o" << endl;
  if ( oMother->keep == YES )  {
    mother = oMother;
    if ( oMothers.first ) mother->children.second.insert(this);
    else if ( oMothers.second ) mother->children.first.insert(this);
  }
  else {
    if ( oMothers.first ) {
      mother = oMother->youngestFirstChild(this, true);
      mother->children.second.insert(this);
    }
    else if ( oMothers.second ) {
      mother = oMother->youngestSecondChild(this, true);
      mother->children.first.insert(this);
    }
    else
      cerr << "single mother, but neither first or second side!! :(" << endl;
  }
}

void RealParton::setOnShell() {
  if ( keep == NO ) return;
  theParton->onShell(true);
  if ( mothers.first ) mothers.first->setOnShell();
  if ( mothers.second ) mothers.second->setOnShell();
  if ( mother ) mother->setOnShell();
}

bool RealParton::moveValenceStatus() {
  if ( !theParton->valence() )
    cerr << "realparton::moveValenceStatus called for non valence parton." << endl;
  if ( !movedValence ) return moveValenceStatus(oldestHeir(this));
  else return moveValenceStatus(movedValence->oldestHeir(this));
}

tRealPartonPtr RealParton::oldestHeir(tRealPartonPtr rp) {
  tRealPartonPtr ret = tRealPartonPtr();
  for ( RealPartonSet::iterator it = oChildren.first.begin();
	it != oChildren.first.end(); it++ ) {
    tRealPartonPtr cand;
    if ( (*it)->interactions.empty() ) continue;
    if ( (*it)->theParton->oY() <= rp->theParton->oY() || (*it)->theParton->valence() )
      cand = (*it)->oldestHeir(rp);
    else
      cand = *it;
    if ( cand && (!ret || cand->theParton->oY() < ret->theParton->oY()) )
      ret = cand;
  }
  for ( RealPartonSet::iterator it = oChildren.second.begin();
	it != oChildren.second.end(); it++ ) {
    tRealPartonPtr cand;
    if ( (*it)->interactions.empty() ) continue;
    if ( (*it)->theParton->oY() <= rp->theParton->oY() || (*it)->theParton->valence() )
      cand = (*it)->oldestHeir(rp);
    else
      cand = *it;
    if ( cand && (!ret || cand->theParton->oY() < ret->theParton->oY()) )
      ret = cand;
  }
  if ( ret && realState->toCheck.find(ret) == realState->toCheck.end() ) {
    // cout << "found oldest heir not in toCheck, search more..." << endl;
    ret = oldestHeir(ret);
  }
  return ret;
}

RealParton::RealPartonSet RealParton::ancestors() {
  RealPartonSet ret;
  ret.insert(this);
  if ( oMothers.first ) {
    RealPartonSet m1 = oMothers.first->ancestors();
    ret.insert(m1.begin(), m1.end());
  }
  if ( oMothers.second ) {
    RealPartonSet m2 = oMothers.second->ancestors();
    ret.insert(m2.begin(), m2.end());
  }
  if ( oMother ) {
    RealPartonSet m = oMother->ancestors();
    ret.insert(m.begin(), m.end());
  }
  return ret;
}

bool RealParton::moveValenceStatus(RealPartonPtr rp) {
  if ( !rp ) return false;
  if ( realState->toCheck.find(rp) == realState->toCheck.end() ) return false;
  if ( rp->theParton->valence() )
    cerr << "RealParton::moveValenceStatus moving valence to a valence" << endl;
  rp->theParton->valencePT(theParton->valencePT());
  rp->theParton->valencePlus(theParton->valencePlus());
  rp->theParton->valence(true);
  realState->valence.insert(rp);
  theParton->valence(false);
  realState->valence.erase(this);
  if ( movedValence ) {
    rp->movedValence = movedValence;
    movedValence->movedValence = rp;
    movedValence = RealPartonPtr();
  }
  else {
    rp->movedValence = this;
    movedValence = rp;
  }
  // cout << "moved valence status from " << oY() << ", (" << theParton->position().x()*GeV << ","
  //      << theParton->position().y()*GeV <<  ") to " << rp->oY() << ", ("
  //      << rp->theParton->position().x()*GeV << "," << rp->theParton->position().y()*GeV <<  ")\n";
  return true;
}

Energy RealParton::effectiveCheckMinus(InvEnergy range, bool firstSide, RealPartonSet checked) {
  RealPartonSet ep = effectiveParton(range, firstSide);
  Energy ret = ZERO;
  for ( RealPartonSet::iterator it = ep.begin();
	it != ep.end(); it++ ) {
    ret += (*it)->checkMinus(checked);
  }
  return ret;
}

Energy RealParton::checkMinus(RealPartonSet checked) {
  if ( checked.find(this) != checked.end() )
    return ZERO;
  checked.insert(this);
  if ( keep == NO ) {
    cout << "NO parton " << oY() << " encountered in RealParton::checkminus()" << endl;
    realState->plotState(true);
    return ZERO;
  }
  Energy ret = minus - givenMinus;
  if ( theParton->valence() ) {
    Energy valenceMinus = sqr(theParton->valencePT().pt())/theParton->valencePlus();
    ret = (minus - valenceMinus) - givenMinus;
    return ret;
  }
  return ret + mothers.first->giveMinus() + mothers.second->giveMinus();
}

Energy RealParton::effectiveGiveMinus(InvEnergy range, bool firstSide) {
  RealPartonSet ep = effectiveParton(range, firstSide);
  Energy ret = ZERO;
  for ( RealPartonSet::iterator it = ep.begin();
	it != ep.end(); it++ ) {
    ret += (*it)->giveMinus();
  }
  return ret;
}

Energy RealParton::giveMinus() {
  Energy ret = ZERO;
  if ( theParton->valence() ) {
    Energy newminus = minus -
      sqr(theParton->theValencePT.pt())/theParton->theValencePlus;
    ret = newminus - givenMinus;
    givenMinus = newminus;
    return ret;
  }
  if ( keep == YES ) {
    ret = minus - givenMinus;
    givenMinus = minus;
  }
  if ( oMothers.first ) ret += oMothers.first->giveMinus();
  if ( oMothers.second ) ret += oMothers.second->giveMinus();
  if ( oMother ) ret += oMother->giveMinus();
  return ret;


  // if ( keep == NO ) {
  //   cout << "NO parton " << oY() << " encountered in RealParton::giveminus()" << endl;
  //   realState->diagnosis(true);
  //   return ZERO;
  // }
  // Energy ret = minus - givenMinus;
  // if ( theParton->valence() ) {
  //   Energy valenceMinus = sqr(theParton->valencePT().pt())/theParton->valencePlus();
  //   ret = (minus - valenceMinus) - givenMinus;
  //   givenMinus = (minus - valenceMinus);
  //   return ret;
  // }
  // givenMinus = minus;
  // if ( nMothers == 2 ) {
  //   if ( mothers.first ) ret += mothers.first->giveMinus();
  //   if ( mothers.second ) ret += mothers.second->giveMinus();
  // }
  // if ( nMothers == 1 )
  //   if ( mother ) ret += mother->giveMinus();
  // return ret;
}

void RealParton::emitRecoiler(TransverseMomentum rec, double plusRatio) {
  //create gluons
  DipolePtr dip;
  PartonPtr p = theParton;
  bool firstSide = false;
  if ( p->dipoles().first && p->dipoles().second ) {
    if ( p->dipoles().first->partons().first->y() > p->dipoles().second->partons().second->y() )
      firstSide = true;
    else firstSide = false;
  }
  else if ( p->dipoles().first )   firstSide = true;
  else  firstSide = false;
  if ( firstSide ) dip = p->dipoles().first;
  else dip = p->dipoles().second;
  PartonPtr recoiler = new_ptr(Parton());
  recoiler->onShell(true);
  dip->generatedGluon(recoiler);
  //fix colour flow for dipole state (let it be connected to the original one)
  dip->splitDipole(0.5); //use better approximation for colour choice?
  //initialise real recoiler
  RealPartonPtr realRecoiler = realState->getReal(recoiler);
  realRecoiler->nMothers = 1;
  realRecoiler->setOMother(this);
  realRecoiler->mother = this;
  realRecoiler->keep = RealParton::YES;

  //set momentum of recoiler
  realRecoiler->pT = rec;
  pT -= rec;
  //otherwise use ymax to limit how far the recoiler can go in rapidity
  realRecoiler->plus = plusRatio*plus;
  plus *= 1.0 - plusRatio;
  realRecoiler->updateYMinus();
  updateYMinus();
  recoiler->oY(realRecoiler->y);

  while ( realRecoiler->minus > dip->dipoleState().collidingEnergy()/10.0 ) {
    plus -= realRecoiler->plus;
    realRecoiler->plus += realRecoiler->plus;
    realRecoiler->updateYMinus();
    updateYMinus();
    recoiler->oY(realRecoiler->y);
  }

  if ( Debug::level > 5 )
    cout << "recoiler takes plus ratio " << plusRatio << endl;

  //set position of recoiler
  recoiler->position(p->position() + p->pTScale()*rec/sqr(rec.pt()));

  //set mother structure
  for ( RealPartonSet::iterator it = children.first.begin();
	it != children.first.end(); it++ ) {
    RealPartonPtr child = *it;
    realRecoiler->children.first.insert(child);
    if ( child->nMothers == 1 ) child->mother = realRecoiler;
    if ( child->nMothers == 2 ) child->mothers.second = realRecoiler;
  }
  children.first.clear();
  if ( firstSide ) children.first.insert(realRecoiler);

  for ( RealPartonSet::iterator it = children.second.begin();
	it != children.second.end(); it++ ) {
    RealPartonPtr child = *it;
    realRecoiler->children.second.insert(child);
    if ( child->nMothers == 1 ) child->mother = realRecoiler;
    if ( child->nMothers == 2 ) child->mothers.first = realRecoiler;
  }
  children.second.clear();
  if ( !firstSide ) children.second.insert(realRecoiler);

  //fix interacting tags
  if ( firstInt ) {
    realRecoiler->interacting = interacting;
    interacting = DipolePtr();
    firstInt = false;
    realRecoiler->firstInt = true;
    realRecoiler->intDist = intDist;
  }
  if ( secondInt ) {
    realRecoiler->interacting = interacting;
    interacting = DipolePtr();
    secondInt = false;
    realRecoiler->secondInt = true;
    realRecoiler->intDist = intDist;
  }
  saveState();
  realRecoiler->saveState();

  if ( Debug::level > 5 ) cout << "added recoiler at " << realRecoiler->oY() << endl;
}

list<double> RealParton::effectiveWeights(const RealPartonSet & partons) {
  Energy totalPlus = ZERO;
  Energy totalMinus = ZERO;
  double maxy = (*partons.begin())->oY();
  //find the total plus in the effective parton to be recoiled.
  for ( RealPartonSet::const_iterator it = partons.begin();
	it != partons.end(); it++ ) {
    totalPlus += (*it)->plus;
    totalMinus += (*it)->minus;
    if ( (*it)->oY() > maxy )
      maxy = (*it)->oY();
  }

  list<double> ret;
  for ( RealPartonSet::const_iterator it = partons.begin();
	it != partons.end(); it++ ) {
    // ret.push_back( ((*it)->oY()==maxy) ? 1.0:0.0 );
    ret.push_back((*it)->plus/totalPlus);
    // ret.push_back((*it)->minus/totalMinus);
    // ret.push_back(1.0/double(partons.size()));
  }
  return ret;
}

void RealParton::saveState() {
  theParton->plus(plus);
  theParton->minus(minus);
  theParton->pT(pT);
  theParton->y(y);
  cMothers = mothers;
  cChildren = children;
  cKeep = keep;
  cInteracting = interacting;
  // if ( Debug::level > 5 && interacting ) cout << "saving " << oY() << " as interacting" << endl; 
  cFirstInt = firstInt;
  cSecondInt = secondInt;
  cIntRecoil = intRecoil;
  cInteractionRecoil.first = interactionRecoil.first;
  cInteractionRecoil.second = interactionRecoil.second;
  cGivenMinus = givenMinus;
  cValence = theParton->valence();
  cValencePlus = theParton->valencePlus();
  cValencePT = theParton->valencePT();
  cMovedValence = movedValence;
  cIntDist = intDist;
  cExchanges = exchanges;
  cFluct = fluct;
  cInteractions.clear();
  cInteractions.insert(interactions.begin(), interactions.end());
}

void RealParton::revert(tDipolePtr intDip) {
  plus1veto = 0;
  plus2veto = 0;
  minus1veto = 0;
  minus2veto = 0;
  unDGLAP = 0;
  plus = theParton->plus();
  minus = theParton->minus();
  pT = theParton->pT();
  y = theParton->y();
  mothers = cMothers;
  children = cChildren;
  keep = cKeep;
  // if ( Debug::level > 5 && interacting && !cInteracting ) cout << "reverting " << oY()
  // 							   << " to non-interacting" << endl;
  interacting = cInteracting;
  firstInt = cFirstInt;
  secondInt = cSecondInt;
  intRecoil = cIntRecoil;
  interactionRecoil.first = cInteractionRecoil.first;
  interactionRecoil.second = cInteractionRecoil.second;
  givenMinus = cGivenMinus;
  if ( cValence && !theParton->valence() )
    realState->valence.insert(this);
  if ( !cValence && theParton->valence() )
    realState->valence.erase(this);
  theParton->valence(cValence);
  if ( cValence ) {
    theParton->valencePT(cValencePT);
    theParton->valencePlus(cValencePlus);
  }
  movedValence = cMovedValence;
  interactions.clear();
  interactions.insert(cInteractions.begin(), cInteractions.end());
  intDist = cIntDist;
  exchanges = cExchanges;
  fluct = cFluct;
}

bool RealParton::checkMomentum() const {
  if ( keep == NO ) return true; //NO partons are not expected to be anything.
  bool ok = true;
  TransverseMomentum expPT = theParton->valencePT() + intRecoil;
  if ( !theParton->valence() ) {
    if ( mothers.first )
      expPT += theParton->recoil(mothers.first->theParton);
    if ( mothers.second )
      expPT += theParton->recoil(mothers.second->theParton);
  }
  for ( RealPartonSet::const_iterator it = children.first.begin();
	it != children.first.end(); it++ ) {
    if ( !(*it)->theParton->valence() )
      expPT += theParton->recoil((*it)->theParton);
  }
  for ( RealPartonSet::const_iterator it = children.second.begin();
	it != children.second.end(); it++ ) {
    if ( !(*it)->theParton->valence() )
      expPT += theParton->recoil((*it)->theParton);
  }
  if ( (expPT - pT).pt() > 0.000000001*GeV ) {
    cout << "* has unexpected pT = (" << pT.x()/GeV << ", " 
	 << pT.y()/GeV << ") at oY = " << theParton->oY() << endl
	 << "* expected was   pT = (" << expPT.x()/GeV << ", " << expPT.y()/GeV << ")" << endl;
    ok = false;
  }
  if ( abs(pT.pt()*exp(y) - minus) > 0.000000001*GeV ||
       abs(pT.pt()*exp(-y) - plus) > 0.000000001*GeV ) {
    cout << "* not on shell, pT = " << pT.pt()/GeV << ", y = " << y << endl
	 << "* expected plus = " << pT.pt()*exp(-y)/GeV 
	 << ", expected minus = " << pT.pt()*exp(y)/GeV << endl
	 << "* plus = " << plus/GeV << ", minus = " << minus/GeV << endl;
      ok = false;
  }
  return ok;
}
