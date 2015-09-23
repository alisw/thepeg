// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DipoleState class.
//

#include "DipoleState.h"
#include "DipoleEventHandler.h"
#include "EventFiller.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/Current.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Utilities/Debug.h"
#include "ThePEG/Utilities/MaxCmp.h"
#include "ThePEG/Utilities/ObjectIndexer.h"
#include "ThePEG/Utilities/UtilityBase.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ParticleInfo.h"
#include "Sum20Momentum.h"
#include "CPUTimer.h"

#include <stdio.h>

#include <iostream>
#include <fstream>

using namespace DIPSY;

DipoleState::DipoleState(const DipoleState & x)
  : thePlus(x.thePlus), theMinus(x.theMinus), theMinusDeficit(x.theMinusDeficit),
    theHandler(x.theHandler), theInitialDipoles(x.theInitialDipoles),
    theSwingCandidates(x.theSwingCandidates),
    theTouchedSwingCandidates(x.theTouchedSwingCandidates),
    theWFInfo(x.theWFInfo),
    theWeight(x.theWeight), doTakeHistory(x.doTakeHistory), theYmax(x.theYmax),
    theCollidingEnergy(x.theCollidingEnergy), theHistory(x.theHistory),
    allDipoles(x.allDipoles), theIncoming(x.theIncoming),
    theProducedParticles(x.theProducedParticles),
    theIncomingShadows(x.theIncomingShadows),
    theShadowPropagators(x.theShadowPropagators) {
  for ( set<DipolePtr>::iterator it = allDipoles.begin();
	it != allDipoles.end(); ++it )
    (**it).theDipoleState = this;
}

DipoleState::~DipoleState() {}

DipoleStatePtr DipoleState::clone() {
  DipoleStatePtr copy = new_ptr(*this);

  for ( set<DipolePtr>::iterator it = allDipoles.begin();
	it != allDipoles.end(); ++it )
    (**it).theDipoleState = this;

  Dipole::CloneSet toclone;
  for ( set<DipolePtr>::iterator it = allDipoles.begin();
	it != allDipoles.end(); ++it ) {
    toclone.insert(*it);
    (**it).fillReferences(toclone);
  }

  Dipole::TranslationMap trans;
  vector<Ariadne5::ClonePtr> clones;
  for ( Dipole::CloneSet::iterator it = toclone.begin();
	it != toclone.end(); ++it )
    if ( *it ) clones.push_back(trans[*it] = (**it).clone());
  for ( int i = 0, N = clones.size(); i < N; ++i ) {
    clones[i]->rebind(trans);
    if ( tDipolePtr d = dynamic_ptr_cast<tDipolePtr>(clones[i]) ) {
      d->theDipoleState = copy;
    }
  }

  copy->allDipoles.clear();
  trans.translate(inserter(copy->allDipoles),
		  allDipoles.begin(), allDipoles.end());

  copy->theInitialDipoles.clear();
  trans.translate(inserter(copy->theInitialDipoles),
		  theInitialDipoles.begin(), theInitialDipoles.end());

  for ( int i = 0, N = theSwingCandidates.size(); i < N; ++i ) {
    copy->theSwingCandidates =
      vector< vector<tDipolePtr> >(theSwingCandidates.size());
    trans.translate(inserter(copy->theSwingCandidates[i]),
		  theSwingCandidates[i].begin(), theSwingCandidates[i].end());
  }

  for ( int i = 0, N = theTouchedSwingCandidates.size(); i < N; ++i ) {
    copy->theTouchedSwingCandidates =
      vector< vector<tDipolePtr> >(theTouchedSwingCandidates.size());
    trans.translate(inserter(copy->theTouchedSwingCandidates[i]),
		    theTouchedSwingCandidates[i].begin(),
		    theTouchedSwingCandidates[i].end());
  }

  return copy;

  // *** TODO *** fix shadows here
}

void DipoleState::sortDipolesFS() {
  theSwingCandidates.clear();
  theTouchedSwingCandidates.clear();
  theSwingCandidates.resize(Current<DipoleEventHandler>()->nColours());
  theTouchedSwingCandidates.resize(Current<DipoleEventHandler>()->nColours());

  list<DipolePtr> dips = getDipoles();
  for ( list<DipolePtr>::const_iterator it = dips.begin(); it != dips.end(); it++ ) {
    sortDipoleFS(**it);
  }
}

void dummycollect(const Dipole & d, int & ndip, Energy & sump, Energy & summ,
		  TransverseMomentum & sumpt, set<tcDipolePtr> & dips) {
  if ( d.children().first && d.children().second ) {
    dummycollect(*d.children().first, ndip, sump, summ, sumpt, dips);
    dummycollect(*d.children().second, ndip, sump, summ, sumpt, dips);
  } else if ( d.children().first ) {
    dummycollect(*d.children().first, ndip, sump, summ, sumpt, dips);
  } else if ( d.children().second ) { //dead end??
    //    dummycollect(*d.children().second, ndip, sump, summ, sumpt, dips);
  } else {
    ++ndip;
    dips.insert(&d);
    if ( d.partons().first->flavour() == ParticleID::g ) {
      sumpt += d.partons().first->pT()/2.0;
      sump += d.partons().first->plus()/2.0;
      summ += d.partons().first->minus()/2.0;
    } else {
      sumpt += d.partons().first->pT();
      sump += d.partons().first->plus();
      summ += d.partons().first->minus();
    }
    if ( d.partons().second->flavour() == ParticleID::g ) {
      sumpt += d.partons().second->pT()/2.0;
      sump += d.partons().second->plus()/2.0;
      summ += d.partons().second->minus()/2.0;
    } else {
      sumpt += d.partons().second->pT();
      sump += d.partons().second->plus();
      summ += d.partons().second->minus();
    }
  }
}

pair<double,int> DipoleState::lambdaMeasure(Energy2 scale,
				 FactoryBase::tH1DPtr histlength,
				 FactoryBase::tH1DPtr histmass) const {
  pair<double,int> lam (0.0, 0);
  for ( set<DipolePtr>::const_iterator it = allDipoles.begin();
	it != allDipoles.end(); it++ ) {
    Dipole & d = **it;
    if ( !d.children().first && !d.children().second && d.participating() ) {
      Energy2 m2 = (d.partons().first->momentum() +
		    d.partons().second->momentum()).m2();
      if ( m2 > ZERO ) lam.first += log(m2/scale);
      ++lam.second;
      if ( histmass && m2 >= ZERO ) histmass->fill(sqrt(m2/GeV2));
      if ( histlength )
	histlength->fill((d.partons().first->position() -
			  d.partons().second->position()).pt()*GeV);
    }
  }
  return lam;
}


void DipoleState::sortFSDipoles() {
  theSwingCandidates.clear();
  theSwingCandidates.resize(Current<DipoleEventHandler>()->nColours());
  for ( set<DipolePtr>::const_iterator it = allDipoles.begin();
	it != allDipoles.end(); it++ ) {
    Dipole & d = **it;
    if ( !d.children().first && !d.children().second && d.participating() ) {
      if ( ( d.neighbors().first &&
	     d.neighbors().first->colour() == d.colour() ) ||
	   ( d.neighbors().second &&
	     d.neighbors().second->colour() == d.colour() ) )
	generateColourIndex(&d);
      forceAt(theSwingCandidates, (**it).colour()).push_back(*it);
    }
  }
  if ( theSwingCandidates.size() != 0 ) return;
  set<tcDipolePtr> dipset1;
  TransverseMomentum sumpt;
  Energy sump = ZERO;
  Energy summ = ZERO;
  int ndip = 0;
  for ( int ic = 0, Nc = theSwingCandidates.size(); ic < Nc; ++ ic )
    for ( int i = 0, N = theSwingCandidates[ic].size(); i < N; ++i ) {
      const Dipole & d = *theSwingCandidates[ic][i];
      dipset1.insert(&d);
      ++ndip;
      if ( d.partons().first->flavour() == ParticleID::g ) {
	sumpt += d.partons().first->pT()/2.0;
	sump += d.partons().first->plus()/2.0;
	summ += d.partons().first->minus()/2.0;
      } else {
	sumpt += d.partons().first->pT();
	sump += d.partons().first->plus();
	summ += d.partons().first->minus();
      }
      if ( d.partons().second->flavour() == ParticleID::g ) {
	sumpt += d.partons().second->pT()/2.0;
	sump += d.partons().second->plus()/2.0;
	summ += d.partons().second->minus()/2.0;
      } else {
	sumpt += d.partons().second->pT();
	sump += d.partons().second->plus();
	summ += d.partons().second->minus();
      }
    }
  cerr << "<DIPSY> total momentum in sortFSDipoles." << endl
       << "        Plus:  " << sump/GeV << endl
       << "        Minus: " << summ/GeV << endl
       << "        x:     " << sumpt.x()/GeV << endl
       << "        y:     " << sumpt.y()/GeV << endl
       << "        Ndip:  " << ndip << " (" << dipset1.size() << ")" << endl;
  sump = ZERO;
  summ = ZERO;
  sumpt = TransverseMomentum();
  ndip = 0;
  set<tcDipolePtr> dipset2;
  for ( int i = 0, N = initialDipoles().size(); i < N; ++i )
    dummycollect(*initialDipoles()[i], ndip, sump, summ, sumpt, dipset2);
  cerr << "<DIPSY> total momentum in getFSSwinger." << endl
       << "        Plus:  " << sump/GeV << endl
       << "        Minus: " << summ/GeV << endl
       << "        x:     " << sumpt.x()/GeV << endl
       << "        y:     " << sumpt.y()/GeV << endl
       << "        Ndip:  " << ndip << " (" << dipset2.size() << ")" << endl;
  set<tcDipolePtr> diffset;
  set_difference(dipset1.begin(), dipset1.end(), dipset2.begin(), dipset2.end(),
		 inserter(diffset));
  cerr << "        " << diffset.size()
       << " dipoles are in sortFSDipoles but not in getFSSwinger" << endl;
  diffset.clear();
  set_difference(dipset2.begin(), dipset2.end(), dipset1.begin(), dipset1.end(),
		 inserter(diffset));
  cerr << "        " << diffset.size()
       << " dipoles are in getFSSwinger but not in sortFSDipoles" << endl;

}



void DipoleState::sortDipoleFS(Dipole & d){
  if(d.children().second && d.children().first) {
    sortDipole(*(d.children().first));
    sortDipole(*(d.children().second)); }
  else if(d.children().first)
    sortDipole(*(d.children().first));
  else if ( d.children().second );
  else if ( !(d.participating()) );
  else forceAt(theSwingCandidates, d.colour()).push_back(&d);
}

void DipoleState::sortDipoles() {
  theSwingCandidates.clear();
  theSwingCandidates.resize(Current<DipoleEventHandler>()->nColours());
  theTouchedSwingCandidates.clear();
  theTouchedSwingCandidates.resize(Current<DipoleEventHandler>()->nColours());

  for ( int i = 0, N = theInitialDipoles.size(); i < N; ++i ) {
    sortDipole(*(theInitialDipoles[i]));
  }
}

void DipoleState::sortDipole(Dipole & d){
  if(d.children().second && d.children().first) {
    sortDipole(*(d.children().first));
    sortDipole(*(d.children().second)); }
  else if(d.children().first)
    sortDipole(*(d.children().first));
  else if ( d.children().second );
  else if ( d.interacted() );
  else if ( !(d.participating()) );
  else if ( !(d.isOn) );
  else {
    forceAt(theSwingCandidates, d.colour()).push_back(& d);
    if ( !d.hasGen() )
      forceAt(theTouchedSwingCandidates, d.colour()).push_back(& d);
  }
}

void DipoleState::evolve(double ymin, double ymax) {
  if ( Current<DipoleEventHandler>()->effectivePartonMode() < 0 && !hasShadows() )
    setupShadows();
  save();
  for ( set<DipolePtr>::const_iterator it = allDipoles.begin(); it != allDipoles.end(); it++ ) {
    (*it)->reset();
    (*it)->touch();
  }
  while ( true ) {
    sortDipoles();
    tDipolePtr sel = tDipolePtr();
    for ( int i = 0, N = initialDipoles().size(); i < N; ++i ) {
      tDipolePtr si = initialDipoles()[i]->getEmitter(ymin, ymax);
      if ( si && si->generatedY() < ymax && si->hasGen() &&
	   ( !sel || si->generatedY() < sel->generatedY() ))
	sel = si;
    }
    if ( sel ) {
      ymin = sel->generatedY();
      sel->emit();
      save();
      continue;
    }
    break;
  }
  theYmax = ymax;
}



void DipoleState::swingFS(double ymin, double ymax) {
  if ( !handler().swingPtr() ) return;
  const Swinger & swinger = handler().swinger();
  for ( set<DipolePtr>::const_iterator it = allDipoles.begin();
	it != allDipoles.end(); it++ ) {
    (*it)->reset();
    (*it)->touch();
  }
  sortFSDipoles();
  double miny = ymin;
  int lastcol = -1;
  while ( true ) {
    MinCmp<double, Dipole *> sel;
    for ( int ic = 0, Nc = theSwingCandidates.size(); ic < Nc; ++ic ) {
      const vector<tDipolePtr> & candidates = swingCandidates(ic);
      if ( lastcol < 0 || lastcol == ic )
	for ( int i = 0, N = candidates.size(); i < N; ++i ) {
	  Dipole & d = *candidates[i];
	  InvEnergy2 saver = d.swingCache;
	  if ( saver < ZERO )
	    saver = swinger.swingDistanceFS(*d.partons().first,
					    *d.partons().second, ymin/GeV2);
	  if ( d.swingDipole() && d.swingDipole()->touched() ) d.reset();
	  d.swingCache = saver;
	}
      for ( int i = 0, N = candidates.size(); i < N - 1; ++i ) {
	Dipole & d1 = *candidates[i];
	bool safe = d1.hasGen();
	if ( lastcol < 0 || lastcol == ic )
	  for ( int j = i + 1; j < N; ++j ) {
	    Dipole & d2 = *candidates[j];
	    if ( safe && !d2.touched() ) continue;
	    InvEnergy2 c =
	      swinger.swingDistanceFS(*d1.partons().first,
				      *d2.partons().second, ymin/GeV2);
	    InvEnergy2 d =
	      swinger.swingDistanceFS(*d2.partons().first,
				      *d1.partons().second, ymin/GeV2);
	    double R = -log( UseRandom::rnd() );
	    double amp = swinger.swingAmpFS(d1.swingCache, d2.swingCache, c, d);
	    double yi = Constants::MaxRapidity;
	    if ( miny*amp + R < Constants::MaxRapidity*amp ) yi = miny + R/amp;
	    if ( yi < d1.generatedY() || !d1.hasGen() ) {
	      d1.swingDipole(&d2);
	      d1.generatedY(yi);
	      d1.recoilSwing(false);
	    }
	  }
	sel(d1.generatedY(), &d1);
	d1.untouch();
      }
      if ( candidates.size() ) candidates.back()->untouch();
    }
    if ( !sel.index() || sel > ymax ) break;
    miny = sel;
    lastcol = sel.index()->colour();
    swinger.recombineFS(*sel.index());
    static DebugItem tupleswing("DIPSY::SwingTuple", 6);
    static ofstream swingtuple;
    static bool isopen = false;
    static long lastevent = 0;
    if ( tupleswing ) {
      if ( !isopen ) {
	string filename = CurrentGenerator::current().filename() + "-swing.tuple";
	swingtuple.open(filename.c_str());
	isopen = true;
      }
      if ( lastevent != CurrentGenerator::current().currentEventNumber() ) {
	lastevent = CurrentGenerator::current().currentEventNumber();
	swingtuple << "# " << lastevent << endl;
      }
      double lrho = -log10(sqrt(miny));
      pair<double,int> lam = lambdaMeasure(0.36*GeV2);
      swingtuple << lrho << '\t' << lam.first/lam.second << '\t' << lam.first << endl;
    }

  }
  theYmax = ymax;

}

void DipoleState::restoreNonparticipants() {
  for ( int i = 0; i < int(theInitialDipoles.size()); i++ ) {
    tPartonPtr p1 = theInitialDipoles[i]->partons().first;
    tPartonPtr p2 = theInitialDipoles[i]->partons().second;
    tPartonPtr p3;
    if ( i != int(theInitialDipoles.size())-1 &&
	 theInitialDipoles[i+1]->partons().first == p2  )
      p3 = theInitialDipoles[i+1]->partons().second;
    else if ( i != int(theInitialDipoles.size())-1 &&
	      theInitialDipoles[i+1]->partons().second == p1  )
      p3 = theInitialDipoles[i+1]->partons().first;
    else if ( i != 0 && theInitialDipoles[i-1]->partons().first == p2  )
      p3 = theInitialDipoles[i-1]->partons().second;
    else if ( i != 0 && theInitialDipoles[i-1]->partons().second == p1  )
      p3 = theInitialDipoles[i-1]->partons().first;

    if ( !(p1->interacted()) &&
	 !(p2->interacted()) &&
	 !( p3 && p3->interacted() )  ) {
      theInitialDipoles[i]->participating(false);
      restoreDipole(p1, p2);
    }
    else
      theInitialDipoles[i]->participating(true);
  }
}

bool DipoleState::restoreDipole(PartonPtr p1, PartonPtr p2) {
  //there must be dipoles to reconnect.
  if ( !p1->dipoles().second || !p2->dipoles().first ) return false;

  //if already connect, do nothing
  if ( p1->dipoles().second == p2->dipoles().first ) return true;

  //if p1 and p2 are connected to different partons, just swing.
  if ( p1->dipoles().second->partons().second != p2->dipoles().first->partons().first ) {
    p1->dipoles().second->swingDipole(p2->dipoles().first);
    p1->dipoles().second->recombine();
    return true;
  }
  //if p1 and p2 are separated by only one parton, then first swing away one of
  //the dipoles, then repeat.
  else {
    sortDipoles();
    //if no swing can be found, then leave it.
    if ( !forceSwing(p1->dipoles().second->partons().second, 0.0, 0.0) ) return false;
    //otherwise go ahead and reconnect the dipoles.
    p1->dipoles().second->swingDipole(p2->dipoles().first);
    p1->dipoles().second->recombine();
    return true;
  }
}

void DipoleState::normaliseValenceCharge(int mode) {
  if ( mode == 0 ) return;
  if ( double(int(theInitialDipoles.size())/3) != double(theInitialDipoles.size())/3.0 ) {
    return;
  }
  if ( mode == 1 ) {
    set<pair<tPartonPtr, tPartonPtr> > valencePartons;
    for ( int i = 0; i < int(theInitialDipoles.size()); i++) {
      tPartonPtr p1 = theInitialDipoles[i]->partons().first;
      tPartonPtr p2 = theInitialDipoles[i]->partons().second;

      //if there are no dipoles to swing, do nothing
      if ( !p1->dipoles().second || !p2->dipoles().first )
	continue;

      //if already connected, do nothing.
      if ( p1->dipoles().second->partons().second == p2 )
	continue;

      //if connected to the same parton, do nothing.
      if ( p1->dipoles().second->partons().second == p2->dipoles().first->partons().first )
	continue;

      //otherwise swing with 50% prob
      if ( UseRandom::rnd() > 0.5 ) {
	p1->dipoles().second->swingDipole(p2->dipoles().first);
	p1->dipoles().second->recombine();
      }
    }
  }
}

void DipoleState::generateColourIndex(tDipolePtr d) {
  int isys = -1;
  if ( d->neighbors().first )
    isys = d->neighbors().first->colourSystem();
  if (  d->neighbors().second ) {
    int isys2 = d->neighbors().second->colourSystem();
    if ( isys >= 0 && isys2 != isys )
      cerr << "Dipole bewteen different colour systems!" << endl;
    isys = isys2;
  }

  do {
    d->colour(UseRandom::irnd(handler().nColours()), isys);
  } while ( ( d->neighbors().first &&
	      d->colour() == d->neighbors().first->colour() ) ||
	    ( d->neighbors().second &&
	      d->colour() == d->neighbors().second->colour() ) );
}

const list<PartonPtr> DipoleState::virtualPartons() const {
  list<PartonPtr> partons;
  for ( set<DipolePtr>::const_iterator it = allDipoles.begin() ; it != allDipoles.end(); it++ ) {
    if ( (*it)->children().first || (*it)->children().second ) continue;
    if ( !((*it)->neighbors().first) && !((*it)->partons().first->interacted()) )
      if ( !((*it)->partons().first->valence()) )
	partons.push_back((*it)->partons().first);
    if ( !((*it)->partons().second->interacted()) && !((*it)->partons().second->valence()) )
      partons.push_back((*it)->partons().second);
  }
  return partons;
}

const set< list<PartonPtr> > DipoleState::loops() const {
  list<PartonPtr> partonsList = getPartons();
  set<PartonPtr> partons;
  set< list<PartonPtr> > ret;
  for( list<PartonPtr>::const_iterator it = partonsList.begin(); 
       it != partonsList.end(); it++ ) {
    partons.insert(*it);
  }

  //go through all partons and add them to colour chains.
  while ( !partons.empty() ) {
    list<PartonPtr> loop;
    PartonPtr head = *partons.begin();
    bool forward = true;
    PartonPtr p = head;

    //follow the colour flow first in forward direction from head until it comes back to the start
    //if a quark is hit, then continue in other direction.
    do {
      //add at end/start of list depending of if we are moving forward or beckwards in colour
      if ( forward ) {
        loop.push_back(p);
      }
      else {
	loop.push_front(p);
      }

      //remove from partons, so we know its added.
      //could potentially be better to do this without destroy the dipoleState...
      if ( !partons.erase(p) ) {
      Throw<DipoleConnectionException>()
	<< "DIPSY found inconsistent DipoleState when extracting "
	"colour singlets!" << Exception::eventerror;
      }

      //move forward in colour if there are more partons, otherwise turn around
      if ( forward ) {
	if ( p->dipoles().second ) {
	  p = p->dipoles().second->partons().second;
	  if ( Debug::level > 5 ) cout << "  continue forward to " << p->y() << endl;
	}
	else {
	  if ( head->dipoles().first ) {
	    //continue at other side, going in other direction
	    p = head->dipoles().first->partons().first;
	  }
	  else {
	    //want to turn around, but "head" is at the other end, so we are done.
	    p = head;

	  }
	  forward = false;
	}
      }
      else {
	if ( p->dipoles().first ) {
	  p = p->dipoles().first->partons().first;
	}
	else {
	  break;
	}
      }
    } while ( p!= head );
    ret.insert(loop);
  }
  return ret;
}

const double DipoleState::highestY() const {
  MaxCmp<double> ret;
  list< PartonPtr > partons = getPartons();
  for ( list<PartonPtr>::const_iterator it = partons.begin();
	it != partons.end(); it++ ) {
    ret((*it)->y());
    ret((*it)->y());
  }
  return ret;
}

const double DipoleState::lowestY() const {
  MinCmp<double> ret;
  list< PartonPtr > partons = getPartons();
  for ( list<PartonPtr>::const_iterator it = partons.begin();
	it != partons.end(); it++ ) {
    ret((*it)->y());
    ret((*it)->y());
  }
  return ret;
}

const double DipoleState::ymax() const {
  return theYmax;
}

const Energy DipoleState::collidingEnergy() const {
  return theCollidingEnergy;
}

void DipoleState::collidingEnergy(Energy E) {
  theCollidingEnergy = E;
}

bool DipoleState::diagnosis(bool print) const {
  bool ok = true;
  int activeDipoles = 0;
  int partons = 0;
  int interactingDips = 0;
  Energy plus = 0.0*GeV;
  Energy minus = 0.0*GeV;
  Energy plusNeeded = 0.0*GeV;
  Energy minusNeeded = 0.0*GeV;
  Energy plusMissing = 0.0*GeV;
  Energy minusMissing = 0.0*GeV;
  bool hasQuark = false;
  TransverseMomentum totalpT;
  Energy maxpT = 0.0*GeV;
  InvEnergy dipSize = ZERO;

  for(set< DipolePtr >::iterator it = allDipoles.begin();
      it != allDipoles.end(); it++) {
    if( ((*it)->children().first || (*it)->children().second) )
      continue;
    activeDipoles++;
    dipSize += (*it)->size();
    if ( (*it)->DGLAPsafe() ) {
      if ( (*it)->neighbors().first &&
	   (*it)->neighbors().first->size() > (*it)->size() && 
           !((*it)->neighbors().first->DGLAPsafe()) )
	Throw<DipoleDGLAPSafeException>()
	  << "DIPSY found broken DGLAPchain!" << Exception::eventerror;

      if ( (*it)->neighbors().second &&
	   (*it)->neighbors().second->size() > (*it)->size() && 
           !((*it)->neighbors().second->DGLAPsafe()) ) 
	Throw<DipoleDGLAPSafeException>()
	  << "DIPSY found broken DGLAPchain!" << Exception::eventerror;
    }
    if ( (*it)->interacted() ) interactingDips++;

    if ( ((*it)->neighbors().first &&
	  (*it)->neighbors().first->children().first) ||
	 ((*it)->neighbors().first &&
	  (*it)->neighbors().first->children().second) ) 
      Throw<DipoleConnectionException>()
	<< "DIPSY found dipole where negiboring dipole has decayed!"
	<< Exception::eventerror;

    if( (*it)->partons().first->dipoles().second != (*it) )
      Throw<DipoleConnectionException>()
	<< "DIPSY found dipole where first parton "
	<< "doesn't point back at the dipole!" << Exception::eventerror;
    if( (*it)->partons().second->dipoles().first != (*it) )
      Throw<DipoleConnectionException>()
	<< "DIPSY found dipole where second parton "
	<< "doesn't point back at the dipole!" << Exception::eventerror;

    if ( (*it)->neighbors().first &&
	 (*it)->neighbors().first->neighbors().second != (*it) )
      Throw<DipoleConnectionException>()
	<< "DIPSY found dipole where neighbouring dipole "
	<< "doesn't point back at the dipole!" << Exception::eventerror;
    
    if ( (*it)->neighbors().second &&
	 (*it)->neighbors().second->neighbors().first != (*it) )
       Throw<DipoleConnectionException>()
	<< "DIPSY found dipole where neighbouring dipole "
	<< "doesn't point back at the dipole!" << Exception::eventerror;

    if ( (*it)->partons().first == (*it)->partons().second )
      Throw<DipoleConnectionException>()
	<< "DIPSY found dipole looping back to itself!" << Exception::eventerror;

    PartonPtr p = (*it)->partons().first;
    partons++;

    if ( p->dipoles().first && p->dipoles().first->partons().second != p ) 
      Throw<DipoleConnectionException>()
	<< "DIPSY found parton where first dipole "
	<< "doesn't point back at the parton!" << Exception::eventerror;

    if ( p->dipoles().second && p->dipoles().second->partons().first != p ) 
      Throw<DipoleConnectionException>()
	<< "DIPSY found parton where second dipole "
	<< "doesn't point back at the parton!" << Exception::eventerror;

    if ( p->dipoles().first && p->dipoles().second &&
	 p->dipoles().first == p->dipoles().second ) 
      Throw<DipoleConnectionException>()
	<< "DIPSY found parton colour connected to itself!"
	<< Exception::eventerror;

    if ( p->dipoles().first && (p->dipoles().first->children().first ||
                                p->dipoles().first->children().second) ) 
      Throw<DipoleConnectionException>()
	<< "DIPSY found parton where first dipole has decayed!"
	<< Exception::eventerror;

    if ( p->dipoles().second && (p->dipoles().second->children().first ||
				 p->dipoles().second->children().second) ) 
      Throw<DipoleConnectionException>()
	<< "DIPSY found parton where second dipole has decayed!"
	<< Exception::eventerror;

    if ( p->y() < -100.0 || p->y() > 100.0 )
      Throw<DipoleKinematicsException>()
	<< "DIPSY found parton impossible rapidity!" << Exception::eventerror;

    if ( p->plus() < ZERO )
      Throw<DipoleKinematicsException>()
	<< "DIPSY found parton at oY " << p->oY() << " has negative plus "
	<< p->plus()/GeV << "!" << Exception::eventerror;

    if ( p->minus() < ZERO ) 
      Throw<DipoleKinematicsException>()
	<< "DIPSY found parton at oY " << p->oY() << " has negative minus "
	<< p->minus()/GeV << "!" << Exception::eventerror;

    totalpT += p->pT();
    if ( p->pT().pt() > maxpT ) {
      maxpT = p->pT().pt();
      if ( maxpT > 1000.0*GeV )
	Throw<DipoleKinematicsException>()
	  << "DIPSY found parton with very large transverse momentum: "
	  << maxpT/GeV << " (p+ = " << p->plus()/GeV << ", p- = "
	  << p->minus()/GeV << ", y = " << p->y() << ")" << Exception::warning;
    }

    if( abs((p->plus()*p->minus() - sqr(p->pT().pt()))/sqr(GeV)) >
	0.001*sqr(p->pT().pt())/sqr(GeV) ) {
      Throw<DipoleKinematicsException>()
	<< "DIPSY found off-shell parton with invariant mass square "
	<< (p->plus()*p->minus() - sqr(p->pT().pt()))/sqr(GeV)
	<< "!" << Exception::eventerror;
      ok = false;
    }

    if( p->rightMoving() ) {
      plus += p->plus();
      minusMissing += p->minus();
      minusNeeded += p->pT().pt()*exp( p->y() );
    }
    else {
      minus += p->minus();
      plusMissing += p->plus();
      plusNeeded += p->pT().pt()*exp( -p->y() );
    }
    if( !((*it)->neighbors().second) ) {
      hasQuark = true;
      p = (*it)->partons().second;
      partons++;
    totalpT += p->pT();
      if( p->rightMoving() ) {
	plus += p->plus();
	minusMissing += p->minus();
	minusNeeded += p->pT().pt()*exp( p->y() );
      }
      else {
	minus += p->minus();
	plusMissing += p->plus();
	plusNeeded += p->pT().pt()*exp( -p->y() );
      }
    }
  }

  if( totalpT.pt() > 0.00000001*GeV ) {
    Throw<DipoleKinematicsException>()
      << "DIPSY found transverse momentum imbalance: " << totalpT.pt()/GeV
      << Exception::warning;
    ok = false;
  }

  Energy acceptableError = 0.1*GeV;
  if ( handler().eventFiller().mode() == 1 )
    acceptableError = 0.000001*GeV;
  if ( abs(plus - thePlus) > 2.0*acceptableError &&
       abs(plus + minus + minusMissing + plusMissing - thePlus - theMinus) >
       2.0*acceptableError ) {
    Throw<DipoleKinematicsException>()
      << "DIPSY found energy non-conservation: "
      << (plus + minus + minusMissing + plusMissing
	  - thePlus - theMinus)/GeV/2.0 << " GeV" << Exception::warning;
    ok = false;
  }

  if ( isnan(plus/GeV) || isnan(minus/GeV) )
    Throw<DipoleKinematicsException>()
      << "DIPSY found parton with nan kinematics!" << Exception::runerror;

  if ( print ) {
    cout << "------------------ state " << this << " ------------------------------\n";
    if (!ok )
      cout << "| NOT OK!!!" << endl;
    cout << setprecision(10);
    cout << 
    "| original plus = " << thePlus/GeV << endl <<
    "| originl minus = " << theMinus/GeV << endl <<
    "|       plus    = " << plus/GeV << endl <<
    "|       minus   = " << minus/GeV << endl <<
    "| missing minus = " << minusMissing/GeV << endl <<
    "| needed minus  = " << minusNeeded/GeV << endl <<
    "| missing plus  = " << plusMissing/GeV << endl <<
    "| needed plus   = " << plusNeeded/GeV << endl <<
    "| total plus    = " << (plus + plusNeeded - plusMissing)/GeV << endl <<
    "| total minus   = " << (minus + minusNeeded - minusMissing)/GeV << endl <<
    "| total pT      = " << "( "<<totalpT.x()/GeV<<", "<<totalpT.y()/GeV<<" )" << endl <<
    "| originl enrgy = " << (thePlus + theMinus)/GeV/2.0 << endl <<
    "| total energy  = " << (plus + minus + minusMissing + plusMissing)/GeV/2.0 << endl <<
    "| missing enrgy = " << (plus+minus+minusMissing+plusMissing-thePlus-theMinus)/GeV/2.0 << endl <<
    "| allDipoles has " << allDipoles.size() << " dipoles." << endl <<
    "| number of active dipoles: " << activeDipoles << endl <<
    "| number of partons " << partons << endl <<
      "| hasQuark: " << hasQuark << endl <<
      "| number of interacting dipoles: " << interactingDips << endl <<
      "| max pT : " << maxpT/GeV << endl <<
      "| average dipole size: " << GeV*dipSize/double(activeDipoles) << endl;

    cout << "----------------------------------------------------------------------\n";
    if (ok ) cout << "| Found no errors! :)\n";
    else cout << "| STATE NOT OK!! :o\n";
    cout << "----------------------------------------------------------------------\n";

  }
  return ok;
}

bool DipoleState::forceSwing(PartonPtr p, double ymax1, double ymax2) {
  DipolePtr d1 = p->dipoles().first;
  DipolePtr d2 = p->dipoles().second;
  if ( !d1 || !d2 )
    cout << "forcing a swing for a parton that doesnt have 2 dipoles!" << endl;
  //check for swings of same colour up to ymax1
  d1->reset();
  d1->touch();
  bool foundFirstDipole = d1->forceGenerateRec(ymax1);
  d2->reset();
  d2->touch();
  bool foundSecondDipole = d2->forceGenerateRec(ymax1);
  if( foundFirstDipole || foundSecondDipole ) {
    if ( !foundFirstDipole )
      d2->recombine();
    else if ( !foundSecondDipole )
      d1->recombine();
    else {
      if ( d1->generatedY() < d2->generatedY() )
	d1->recombine();
      else
	d2->recombine();
    }
      return true;
  }

  //check for swings of other colours up to ymax2
  int trueColour = d1->colour();
  int trueSys = d1->colourSystem();
  for ( int c=0;c < Current<DipoleEventHandler>()->nColours();c++ ) {
    d1->colour(c, trueSys);
    d1->forceGenerateRec(ymax2);
  }
  d1->colour(trueColour);
  foundFirstDipole = d1->swingDipole();
  trueColour = d2->colour();
  trueSys = d2->colourSystem();
  for ( int c=0;c < Current<DipoleEventHandler>()->nColours();c++ ) {
    d2->colour(c, trueSys);
    d2->forceGenerateRec(ymax2);
  }
  d2->colour(trueColour);
  foundSecondDipole = d2->swingDipole();
  if( foundFirstDipole || foundSecondDipole ) {
    if ( !foundFirstDipole )
      d2->recombine();
    else if ( !foundSecondDipole )
      d1->recombine();
    else {
      if ( d1->generatedY() < d2->generatedY() )
	d1->recombine();
      else
	d2->recombine();
      return true;
    }
  }
  return false;
}

void DipoleState::checkFSMomentum() const {
  static DebugItem checkkinematics("DIPSY::CheckKinematics", 6);
  if ( ! checkkinematics ) return;
  Sum20Momentum sum20(lightCone(thePlus, theMinus));
  list<PartonPtr> pl = getPartons();
  for ( list<PartonPtr>::iterator it = pl.begin(); it != pl.end(); ++it )
    if ( (**it).onShell() ) sum20 -= (**it).momentum();

  if ( !sum20 ) {
    Throw<DipoleKinematicsException>()
      << "DIPSY found energy-momentum non-conservation in final DipoleState."
      << Exception::warning;
    debugShadowTree();
  }
}

void DipoleState::checkFSMomentum(const Step & step) const {
  static DebugItem checkkinematics("DIPSY::CheckKinematics", 6);
  if ( ! checkkinematics ) return;
  Sum20Momentum sum20(step.collision()->incoming().first->momentum() +
		      step.collision()->incoming().second->momentum());

  tPVector fs = step.getFinalState();
  for ( int i = 0, N = fs.size(); i < N; ++i ) sum20 -= fs[i]->momentum();

  if ( !sum20 )
    Throw<DipoleKinematicsException>()
      << "DIPSY found energy-momentum non-conservation in initial Step."
      << Exception::warning;
}

DipoleStatePtr DipoleState::merge(DipoleStatePtr otherState) {
 //reasign all other dipoles to belong to this state.
  for( set< DipolePtr >::iterator it = otherState->allDipoles.begin(); 
      it != otherState->allDipoles.end(); it++ ) {
    (*it)->dipoleState(this);
  }
  //fill up initialDipoles and allDipoles with the other state.
  theInitialDipoles.insert(theInitialDipoles.begin(),
                           otherState->initialDipoles().begin(),
                           otherState->initialDipoles().end());
  theIncoming.insert(otherState->incoming().begin(), otherState->incoming().end());
  allDipoles.insert(otherState->allDipoles.begin(),otherState->allDipoles.end());
  theProducedParticles.insert(otherState->theProducedParticles.begin(),
			      otherState->theProducedParticles.end());
  theIncomingShadows.insert(theIncomingShadows.end(),
			    otherState->theIncomingShadows.begin(),
			    otherState->theIncomingShadows.end());
  //add up the original particles momenta
  thePlus += otherState->thePlus;
  theMinus += otherState->theMinus;

  return this;
}

DipoleStatePtr DipoleState::collide( DipoleStatePtr otherState,
				     const vector<FList::const_iterator> & sel,
                                     const ImpactParameters & b ) {
  //mirror the other state in y
  otherState->mirror(0.0);

  //move the other state according to the impact parameter
  otherState->translate(b);

  //reasign all other dipoles to belong to this state.
  for( set< DipolePtr >::iterator it = otherState->allDipoles.begin(); 
      it != otherState->allDipoles.end(); it++ ) {
    (*it)->dipoleState(this);
  }

  //fill up initialDipoles and allDipoles with the other state.
  theInitialDipoles.insert(theInitialDipoles.begin(),
                           otherState->initialDipoles().begin(),
                           otherState->initialDipoles().end());
  allDipoles.insert(otherState->allDipoles.begin(),otherState->allDipoles.end());

  //decide if each collision in sel is elastic or not, and recouple the non-elastic ones
  for( vector< FList::const_iterator>::const_iterator it = sel.begin(); 
      it != sel.end(); it++ ) {
    if( (*it)->first.second < 2.0*UseRandom::rnd() ) {
      (*it)->second.first->swingDipole( (*it)->second.second );
      (*it)->second.first->recombine();
    }
  }
  return this;
}

void DipoleState::mirror(double y0) {
  for ( set< DipolePtr >::iterator it = allDipoles.begin(); 
	it != allDipoles.end(); it++) {
    if ( !((*it)->children().first || (*it)->children().second) ) {
      (*it)->partons().first->mirror(y0);
      if ( !((*it)->neighbors().second) ) (*it)->partons().second->mirror(y0);
    }
  }

  swap(thePlus, theMinus);

  for ( map<PPtr, vector<PartonPtr> >::iterator it = theIncoming.begin();
	it != theIncoming.end(); ++it )
    it->first->setMomentum(lightCone(it->first->momentum().minus(),
				     it->first->momentum().plus()));

  for ( int i = 0, N = theIncomingShadows.size(); i < N; ++i )
    theIncomingShadows[i]->mirror(y0);

}

void DipoleState::translate(const ImpactParameters & b) {
  for(set< DipolePtr >::iterator it = allDipoles.begin(); 
      it != allDipoles.end(); it++) {
    if( !((*it)->children().first || (*it)->children().second) ) {
      b.translate((*it)->partons().first);
      if( !((*it)->neighbors().second) )
	b.translate((*it)->partons().second);
    }
  }
  for ( int i = 0, N = theIncomingShadows.size(); i < N; ++i )
    theIncomingShadows[i]->translate(b);
}

void DipoleState::fixValence(Step & step) const {
  for ( map<PPtr, vector<PartonPtr> >::const_iterator it = theIncoming.begin();
	it != theIncoming.end(); ++it ) {
    tcWFInfoPtr wfi = WFInfo::getWFInfo(*it->first);
    if ( !wfi ) continue;
    PVector valence;
    for ( int i = 0, N = it->second.size(); i < N; ++i ) {
      valence.push_back(getParticle(it->second[i], true));
      if ( !valence.back() ) valence.pop_back();
    }
    wfi->wf().fixValence(step, it->first, valence);
    checkFSMomentum(step);
  }
  LorentzRotation Rshift = Utilities::boostToCM(step.collision()->incoming());
  Utilities::transform(step.particles(), Rshift);  
}

LorentzMomentum DipoleState::changeMass(LorentzMomentum p, Energy m, LorentzMomentum ref,
  LorentzRotation * Rshift) {
  Energy2 s = (p + ref).m2();
  LorentzRotation R(0.0, 0.0, -(p.z() + ref.z())/(p.e() + ref.e()));
  Energy z = Math::sign(SimplePhaseSpace::getMagnitude(s, ref.mt(), sqrt(sqr(m) + p.perp2())), (R*p).z());
  LorentzMomentum pnew(p.x(), p.y(), z, sqrt(sqr(m) + sqr(z) + p.perp2()));
  LorentzMomentum refnew(ref.x(), ref.y(), -z, sqrt(sqr(z) + ref.mt2()));
  s = (pnew + refnew).m2();
  LorentzRotation shift = LorentzRotation(0.0, 0.0, ref.z()/ref.e())*
    LorentzRotation(0.0, 0.0, -refnew.z()/refnew.e());
  pnew *= shift;
  refnew *= shift;
  if ( Rshift ) *Rshift = R.boostZ((pnew.z() + ref.z())/(pnew.e() + ref.e()));
  return pnew;
}



vector<pair<Parton::Point, InvEnergy> > DipoleState::points() {
  list< PartonPtr > partons = getPartons();
  vector<pair<Parton::Point, InvEnergy> > ret;
  for ( list<PartonPtr>::iterator it = partons.begin(); it != partons.end(); it++) {
    InvEnergy r = max((*it)->dipoles().first->size(), (*it)->dipoles().first->size());
    ret.push_back(make_pair((*it)->position(), r));
  }
  return ret;
}

list< PartonPtr > DipoleState::getPartons() const {
  list<PartonPtr> partons;
  for ( set<DipolePtr>::const_iterator it = allDipoles.begin() ; it != allDipoles.end(); it++ ) {
    if ( (*it)->children().first || (*it)->children().second ) continue;
      partons.push_back((*it)->partons().first);
      if ( !((*it)->neighbors().second) )   //end of chain --> add also right parton
      partons.push_back((*it)->partons().second);
  }
  return partons;
}

list< DipolePtr > DipoleState::getDipoles() const {
  list<DipolePtr> dipoles;
  for ( set<DipolePtr>::const_iterator it = allDipoles.begin() ; it != allDipoles.end(); it++ ) {
    if ( (*it)->children().first || (*it)->children().second ) continue;
    dipoles.push_back(*it);
  }
  return dipoles;
}

void DipoleState::makeOriginal() {
  //clear the original dipoles
  theInitialDipoles.clear(); //owns the dipoles
  theSwingCandidates.clear(); //trancient
  theTouchedSwingCandidates.clear(); //trancient

  //go through allDipoles and set them to on or off.
  for ( set<DipolePtr>::iterator it = allDipoles.begin();
	it != allDipoles.end(); it++ ) {
    DipolePtr d = *it;
    if ( d->children().first || d->children().second ) d->isOn = false;
    else d->isOn = true;
  }

  //go through allDipoles again
  for ( set<DipolePtr>::iterator it = allDipoles.begin();
	it != allDipoles.end(); it++ ) {
    DipolePtr d = *it;

    //remove neighbors that are off, both from dipole and parton
    if ( d->isOn ) {
      if ( d->neighbors().first && !d->neighbors().first->isOn ) {
	d->neighbors(make_pair(DipolePtr(), d->neighbors().second));
	d->partons().first->dipoles(make_pair(DipolePtr(), d->partons().first->dipoles().second) );
      }
      if ( d->neighbors().second && !d->neighbors().second->isOn ) {
	d->neighbors(make_pair(d->neighbors().first, DipolePtr()));
	d->partons().second->dipoles(make_pair(d->partons().second->dipoles().first, DipolePtr()) );
      }
    }
  }

  //go through allDipoles again
  for ( set<DipolePtr>::iterator it = allDipoles.begin();
	it != allDipoles.end(); ) {
    DipolePtr d = *it;
    //remove all pointers to off dipoles (not needed, as transient pointers)

    //insert in original dipoles if the dipole is on.
    if ( d->isOn ) {
      theInitialDipoles.push_back(d);
      d->partons().first->parents() = make_pair(PartonPtr(), PartonPtr());
      d->partons().second->parents() = make_pair(PartonPtr(), PartonPtr());
      d->partons().first->removeAllChildren();
      d->partons().second->removeAllChildren();
      d->reset();
      it++;
    }

    //remove the dipole from allDipoles if it is off.
    if ( !d->isOn )  {
      set<DipolePtr>::iterator jt = it;
      it++;
      allDipoles.erase(jt);
    }
  }

  //set parton data as if they were valence partons
  list<PartonPtr> partons = getPartons();
  int i = 0;
  for( list<PartonPtr>::iterator it = partons.begin();
       it != partons.end(); it++ ) {
    PartonPtr p = *it;
    p->number(i);
    p->ordered(true);
    p->valencePlus(p->plus());
    p->valencePT(p->pT());
    p->oY(p->y());
    p->parents(make_pair(PartonPtr(), PartonPtr()));
    i++;
  }
}

vector<DipoleState::String> DipoleState::strings() {
  LorentzMomentum sum = lightCone(thePlus, theMinus);
  set< list<PartonPtr> > colourLoops = loops();
  vector<DipoleState::String> ret;
  for ( set< list<PartonPtr> >::const_iterator loop = colourLoops.begin();
	loop != colourLoops.end(); loop++ ) {
    DipoleState::String s;
    for ( list<PartonPtr>::const_iterator p = (*loop).begin(); p != (*loop).end(); p++ ) {
      if ( !(**p).onShell() ) cerr << "The FUCK?!? (loops()!)" << endl;
      s.push_back( *p );
      sum -= (**p).momentum();
    }
    if ( Debug::level > 5 ) cout << "new string of size " << s.size() << endl;
    ret.push_back( s );
  }
  static DebugItem checkkinematics("DIPSY::CheckKinematics", 6);
  if ( ! checkkinematics ) return ret;

  if ( abs(sum.plus()) > 0.001*MeV || abs(sum.minus()) > 0.001*MeV ||
       abs(sum.x()) > 0.001*MeV || abs(sum.x()) > 0.001*MeV )
    Throw<DipoleKinematicsException>()
      << "DIPSY found energy-momentum non-conservation in final strings."
      << Exception::warning;

  return ret;
}

//Assume incoming particles massless for now. Valid as long higest pT 
//is much larger than the mass.
void DipoleState::balanceMomenta() {
  //Sum up how much p-/p+ is needed on the left/right side.
  //Assumes that no particles got pushed over to the other side during reabsorbtion.
  Energy PlmMissing = 0.0*GeV;
  Energy PlmNeeded = 0.0*GeV;
  Energy PrpMissing = 0.0*GeV;
  Energy PrpNeeded = 0.0*GeV;
  Energy Plp = 0.0*GeV;
  Energy Prm = 0.0*GeV;
  PartonPtr p;
  for(set< DipolePtr >::iterator it = allDipoles.begin();
      it != allDipoles.end(); it++) {
    if( ((*it)->children().first || (*it)->children().second) )
      continue;
    //Do first parton always
    p = (*it)->partons().first;
    if( p->rightMoving() ) {
      PlmMissing += p->minus();
      PlmNeeded += p->pT().pt()*exp( p->y() );
      Plp += p->plus();
    }
    else {
      Prm += p->minus();
      PrpMissing += p->plus();
      PrpNeeded += p->pT().pt()*exp( -p->y() );
    }
    if( !((*it)->neighbors().second) ) {
      //do second parton as well if at the end of a chain
      p = (*it)->partons().second;
      if( p->rightMoving() ) {
        PlmMissing += p->minus();
        PlmNeeded += p->pT().pt()*exp( p->y() );
        Plp += p->plus();
      }
      else {
        Prm += p->minus();
        PrpMissing += p->plus();
        PrpNeeded += p->pT().pt()*exp( -p->y() );
      }
    }
  }

  Energy PpShuffled = PrpNeeded - PrpMissing;
  Energy PmShuffled = PlmNeeded - PlmMissing;

  double a = -(1 + PmShuffled/Prm);
  double b = -PlmNeeded/Prm;
  double c = -(1 + PpShuffled/Plp);
  double d = -PrpNeeded/Plp;
  double e = (b - d)/c + a;
  //how much p+/p- has to be scaled down on left/right side.
  double q2 = - e/2 + sqrt( sqr(e)/4 + d*a/c );
  double q1 = (d - c*q2)/q2;
  if( max(abs(log(q1)),abs(log(q2))) > 0.1 ) {
  cout << "yl shifted " << log(q2) << endl <<
    "yr shifted " << log(q1) << endl;
  }

  if ( PlmMissing > Prm ) {
      cout << "not enough p- in right state! VETO!" << endl;
  }
  else if ( PrpMissing > Plp ) {
      cout << "not enough p+ in left state! VETO!" << endl;
  }

  //now shift all the partons in rapidity with ln q1 or ln q2 on left/right side.
  for(set< DipolePtr >::iterator it = allDipoles.begin();
      it != allDipoles.end(); it++) {
    if( ((*it)->children().first || (*it)->children().second) )
      continue;
    //Do first parton always
    p = (*it)->partons().first;
    if( p->rightMoving() ) {
      p->plus( p->plus()*q1 );
      p->y( log(p->pT().pt()/p->plus() ) );
      p->minus( p->pT().pt()*exp( p->y() ) );
    }
    else {
      p->minus( p->minus()*q2 );
      p->y( log( p->minus()/p->pT().pt() ) );
      p->plus( p->pT().pt()*exp( -p->y() ) );
    }
    if( !((*it)->neighbors().second) ) {
      //do second parton as well if at the end of a chain
      p = (*it)->partons().second;
      if( p->rightMoving() ) {
	p->plus( p->plus()*q1 );
	p->y( log(p->pT().pt()/p->plus() ) );
	p->minus( p->pT().pt()*exp( p->y() ) );
      }
      else {
	p->minus( p->minus()*q2 );
	p->y( log( p->minus()/p->pT().pt() ) );
	p->plus( p->pT().pt()*exp( -p->y() ) );
      }
    }
  }
}

void DipoleState::saveGluonsToFile(double weight) const {
  EventPtr event = handler().currentEvent();
  string filename = handler().generator()->filename();

    //set up outfile
  ostringstream os;

  PPair incoming = event->incoming();
  LorentzPoint p1 = incoming.first->vertex();
  LorentzPoint p2 = incoming.second->vertex();
  double bx = (p1-p2).x()/femtometer;
  double by = (p1-p2).y()/femtometer;

  //TODO: interface!
  os << "/home/christoffer/DIPSYevents/" << filename << "/event" << event->number() << ".dat";

  // os << "/scratch/parton/christof/DIPSY/events/" << filename << "/event" << weight << ".dat";
  weight = 2.0*int(sqrt(sqr(bx)+sqr(by))) + 1.0;

  ofstream outfile (os.str().c_str());

  //general information
  outfile << "Generated with DIPSY 1.1" << endl;
  outfile << "These are the real gluons only. No FSR or hadronisation has been done." << endl;

  //print weight and impact parameter.

  outfile << endl << "eventweight: " << weight << endl << endl;

  outfile << "impact_parameter(fm):" << " " << bx << " " << by << endl;

  vector<PartonPtr> partonVector;
  int n = 0;
  for(set< DipolePtr >::iterator it = allDipoles.begin();
      it != allDipoles.end(); it++) {
    if( ((*it)->children().first || (*it)->children().second) )
      continue;
    PartonPtr p = (*it)->partons().first;

    if ( p->y() < -3 || p->y() > 3 ) continue;

    partonVector.push_back(p);
    p->number(n);
    n++;
  }

  outfile << "number_of_particles:" << " " << partonVector.size() << endl;
  outfile << "number of spectating nucleons: " << numberOfSpectators() << endl << endl;
  
  //print gluons
  outfile << "particle_number" << " "
	  << "transverse_position_x(fm)" << " "
	  << "transverse_position_y(fm)" << " "
	  << "transverse_momentum_x(GeV)" << " "
	  << "transverse_momentum_y(GeV)" << " "
	  << "rapidity" << " "
	  << "colour_neighbour_number" << " "
	  << "anticolour_neighbour_number" << endl << endl;


  for(int i = 0; i < int(partonVector.size()); i++) {
    PartonPtr p = partonVector[i];
    outfile << p->number() << " "
	    << p->position().x()*hbarc/femtometer << " " 
	    << p->position().y()*hbarc/femtometer << " " 
	    << p->pT().x()/GeV << " " 
	    << p->pT().y()/GeV << " " 
	    << p->y() << " ";
    if ( p->dipoles().first )
      outfile << p->dipoles().first->partons().first->number() << " ";
    else outfile << -1 << " ";
    if ( p->dipoles().second )
      outfile << p->dipoles().second->partons().second->number() << " ";
    else outfile << -1 << " ";
    outfile  << endl;
  }
  outfile.close();
  cout << "printed gluons to file " << os.str().c_str() << endl;
}

int DipoleState::numberOfSpectators() const {
  int nSpectators = 0;
  for ( int i = 0; i < int(theInitialDipoles.size()); i++ ) {
    tPartonPtr p1 = theInitialDipoles[i]->partons().first;
    tPartonPtr p2 = theInitialDipoles[i]->partons().second;
    tPartonPtr p3;
    if ( i != int(theInitialDipoles.size())-1 &&
	 theInitialDipoles[i+1]->partons().first == p2  )
      p3 = theInitialDipoles[i+1]->partons().second;
    else if ( i != int(theInitialDipoles.size())-1 &&
	      theInitialDipoles[i+1]->partons().second == p1  )
      p3 = theInitialDipoles[i+1]->partons().first;
    else if ( i != 0 && theInitialDipoles[i-1]->partons().first == p2  )
      p3 = theInitialDipoles[i-1]->partons().second;
    else if ( i != 0 && theInitialDipoles[i-1]->partons().second == p1  )
      p3 = theInitialDipoles[i-1]->partons().first;

    if ( !(p1->interacted()) &&
	 !(p2->interacted()) &&
	 !( p3 && p3->interacted() )  ) {
      nSpectators++;
    }
  }
  return nSpectators/3;
}

void DipoleState::unifyColourSystems(int sys) {
  for ( set<DipolePtr>::iterator it = allDipoles.begin();
	it != allDipoles.end(); ++it ) (**it).colourSystem(sys);
}

tPPtr DipoleState::getParticle(tcPartonPtr parton, bool nocreate) const {
  map<tcPartonPtr,PPtr>::iterator pit = theProducedParticles.find(parton);
  if ( pit != theProducedParticles.end() ) return pit->second;
  if ( nocreate ) return tPPtr();
  if( parton->plus() <= 0.0*GeV ) cout << "p+ is not positive :(" << endl <<
			    "p+ = " << parton->plus()/GeV << endl <<
			    "p- = " << parton->minus()/GeV << endl <<
			    "pt = " << parton->pT().pt()/GeV << endl <<
			    "y = " << parton->y() << endl;
  PPtr p = CurrentGenerator::current().getParticle(parton->flavour());
  p->setMomentum(parton->momentum());
  p->setVertex(LorentzPoint(parton->position().x()*hbarc
			    , parton->position().y()*hbarc,
			    ZERO, ZERO));
  p->getInfo().push_back(new_ptr(ParticleInfo(parton)));
  theProducedParticles[parton] = p;
  
  return p;
}

void DipoleState::setupShadows() {
  set<tPartonPtr> valence;
  for ( int i = 0, N = initialDipoles().size(); i < N ; ++i ) {
    valence.insert(initialDipoles()[i]->partons().first);
    valence.insert(initialDipoles()[i]->partons().second);
  }
  for ( set<tPartonPtr>::iterator it = valence.begin();
	it != valence.end(); ++it ) {
    theIncomingShadows.push_back(ShadowParton::createValence(**it));
    LorentzMomentum prop = lightCone(plus(), minus(), TransverseMomentum());
    for ( set<tPartonPtr>::iterator ito = valence.begin();
	  ito != valence.end(); ++ito )
      if ( ito != it ) prop -=  (**ito).momentum();
    theShadowPropagators[theIncomingShadows.back()] = prop;
  }
}
  
void DipoleState::resetShadows() {
  for ( int i = 0, N = theIncomingShadows.size(); i < N; ++i )
    theIncomingShadows[i]->reset();
}

void DipoleState::resetInteractedShadows() {
  for ( int i = 0, N = theIncomingShadows.size(); i < N; ++i )
    theIncomingShadows[i]->resetInteracted();
}

LorentzMomentum DipoleState::incomingMomentum(tcSPartonPtr valence, int mode) {
  if ( mode >= 0 ) {
    for ( int i = 0, N = theIncomingShadows.size(); i < N; ++i )
      if ( theIncomingShadows[i] != valence ) {
	theIncomingShadows[i]->resetMomentum0();
	theIncomingShadows[i]->setOnShell(mode);
      }
  }
  return theShadowPropagators[valence];
}

void DipoleState::debugShadowTree() const {
  cerr << "=== Shadow Tree Structure ===" << endl;
    for ( int i = 0, N = theIncomingShadows.size(); i < N; ++i )
      theIncomingShadows[i]->debugTree("");
}

void DipoleState::checkShadowMomentum(Sum20Momentum & sum20,
				      const ImpactParameters * b) const {
  sum20 += b? lightCone(theMinus, thePlus): lightCone(thePlus, theMinus);
    for ( int i = 0, N = theIncomingShadows.size(); i < N; ++i )
      theIncomingShadows[i]->checkMomentum(sum20, b);
}


void DipoleState::persistentOutput(PersistentOStream & os) const {
  os << ounit(thePlus, GeV) << ounit(theMinus, GeV)
     << ounit(theMinusDeficit, GeV) << theHandler << theInitialDipoles
     << theSwingCandidates << theTouchedSwingCandidates
     << theWFInfo << theWeight << doTakeHistory 
     << theYmax << ounit(theCollidingEnergy, GeV) << allDipoles << theIncoming
     << theProducedParticles << theIncomingShadows;
}

void DipoleState::persistentInput(PersistentIStream & is, int) {
  is >> iunit(thePlus, GeV)>> iunit(theMinus, GeV)
     >> iunit(theMinusDeficit, GeV) >> theHandler >> theInitialDipoles
     >> theSwingCandidates >> theTouchedSwingCandidates
     >> theWFInfo >> theWeight >> doTakeHistory 
     >> theYmax >> iunit(theCollidingEnergy, GeV) >> allDipoles >> theIncoming
     >> theProducedParticles >> theIncomingShadows;
}


// Static variable needed for the type description system in ThePEG.
#include "ThePEG/Utilities/DescribeClass.h"
DescribeClass<DipoleState,PersistentBase>
  describeDIPSYDipoleState("DIPSY::DipoleState", "libAriadne5.so libDIPSY.so");


void DipoleState::Init() {}

