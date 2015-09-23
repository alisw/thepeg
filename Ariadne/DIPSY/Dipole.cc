// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Dipole class.
//

#include "Dipole.h"
#include "DipoleState.h"
#include "DipoleEventHandler.h"
#include "ThePEG/Utilities/Current.h"
#include "ThePEG/Utilities/ObjectIndexer.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Utilities/DebugItem.h"

using namespace DIPSY;

Dipole::~Dipole() {
  if ( effmap ) delete effmap;
}

Ariadne5::ClonePtr Dipole::clone() const {
  return new_ptr(*this);
}

void Dipole::fillReferences(CloneSet & cs) const {
  cs.insert(thePartons.first);
  cs.insert(thePartons.second);
  cs.insert(theGeneratedGluon);
}

void Dipole::rebind(const TranslationMap & trans) {
  thePartons.first = trans.translate(thePartons.first);
  thePartons.second = trans.translate(thePartons.second);
  theNeighbors.first = trans.translate(theNeighbors.first);
  theNeighbors.second = trans.translate(theNeighbors.second);
  theChildren.first = trans.translate(theChildren.first);
  theChildren.second = trans.translate(theChildren.second);
  theGeneratedGluon = trans.translate(theGeneratedGluon);
  theSwingDipole = trans.translate(theSwingDipole);
}

void Dipole::interact(Dipole & x, bool onegluon) {
  theInteracted = &x;
  if ( onegluon && partons().first && partons().second ) {
    InvEnergy2 d11 = x.partons().first?
      partons().first->dist2(*x.partons().first): InvEnergy2();
    InvEnergy2 d12 = x.partons().second?
      partons().first->dist2(*x.partons().second): InvEnergy2();
    InvEnergy2 d21 = x.partons().first?
      partons().second->dist2(*x.partons().first): InvEnergy2();
    InvEnergy2 d22 = x.partons().second?
      partons().second->dist2(*x.partons().second): InvEnergy2();
    if ( UseRandom::rndbool((d11 + d12)/(d11 + d12 + d21 + d22)) )
      partons().first->interact(true);
    else
      partons().second->interact(true);
    return;
  }
  if ( partons().first ) partons().first->interact();
  if ( partons().second ) partons().second->interact();
}

pair <InvEnergy, InvEnergy> Dipole::interactionLengths() const {
  if ( !interacted() ) return make_pair( 0.0/GeV, 0.0/GeV );
  return make_pair( sqrt(partons().first->dist2
			 (*(interacted()->partons().second))),
		    sqrt(partons().second->dist2
			 (*(interacted()->partons().first))) );
}

InvEnergy Dipole::rMax() const {
  return theDipoleState->handler().rMax();
}

void Dipole::reset() {
  theGeneratedGluon = PartonPtr();
  theSwingDipole = tDipolePtr();
  theGeneratedY = Constants::MaxRapidity;
  swingCache = -1.0/GeV2;
  if ( effmap ) delete effmap;
  effmap = 0;
}

void Dipole::colExtract(vector< vector<tDipolePtr> > & colvec) {
  if ( children().first ) children().first->colExtract(colvec);
  if ( children().second ) children().second->colExtract(colvec);
  if ( !children().first && !children().second )
    forceAt(colvec, colour()).push_back(this);
}

tDipolePtr Dipole::getFSSwinger(double miny, double maxy) {
  // *** TODO *** Check if dead-end before recursion
  if ( children().first && children().second ) {
    tDipolePtr em1 = children().first->getFSSwinger(miny, maxy);
    tDipolePtr em2 = children().second->getFSSwinger(miny, maxy);
    return em1->generatedY() < em2->generatedY()? em1: em2;
  } else if ( children().first ) {
    return children().first->getFSSwinger(miny, maxy);
  } else if ( children().second ) { //dead end??
    //    return children().second->getFSSwinger(miny, maxy);
  }
  generateFSRec(miny, maxy);
  return this;
}

void Dipole::generateFSRec(double miny, double maxy) {
  if ( dipoleState().handler().swingPtr() && participating() &&
       !children().second )
    dipoleState().handler().swinger().generateFS(*this, miny, maxy);
}

tDipolePtr Dipole::getEmitter(double miny, double maxy) {
  if ( !isOn ) return DipolePtr();
  if ( children().first && children().second ) {
    tDipolePtr em1 = children().first->getEmitter(miny, maxy);
    tDipolePtr em2 = children().second->getEmitter(miny, maxy);
    if ( !em1 ) return em2? em2:em1;
    else if ( !em2 ) return em1;
    return em1->generatedY() < em2->generatedY()? em1: em2;
  } else if ( children().first ) {
    return children().first->getEmitter(miny, maxy);
  } else if ( children().second ) {
    return DipolePtr();
  }
  if ( !hasGen() || ( effectivePartons().first && effectivePartonsChanged() ) ) {
    generate(miny, maxy);
  }
  else {
    generateRec(miny, min(generatedY(), maxy), false);
  }
  return this;
}

void Dipole::generate(double miny, double maxy) {
  reset();
  dipoleState().handler().emitter().generate(*this, miny, maxy);
  generateRec(miny, min(generatedY(), maxy), true);
}

void Dipole::generateRec(double miny, double maxy, bool force) {
  if ( dipoleState().handler().swingPtr() )
    dipoleState().handler().swinger().generate(*this, miny, maxy, force);
}

bool Dipole::forceGenerateRec(double ymax) {
  return dipoleState().handler().swingPtr()?
    dipoleState().handler().swinger().forceGenerate(*this, ymax): false;
}

void Dipole::recombine() {
  if ( dipoleState().handler().swingPtr() )
    dipoleState().handler().swinger().recombine(*this);
}

void Dipole::absorb() {
  PartonPtr p = new_ptr(Parton());
  PartonPtr p1 = partons().first;
  PartonPtr p2 = partons().second;
  DipolePtr d1 = neighbors().first;
  DipolePtr d2 = neighbors().second;

  p->pT( p1->pT() + p2->pT() );
  p->plus( p1->plus() + p2->plus() );
  p->y( log(p->pT().pt()/p->plus()) );
  p->minus( p->pT().pt()*exp(p->y()) );

  d1->neighbors( make_pair( d1->neighbors().first, d2 ) );
  d2->neighbors( make_pair( d1, d2->neighbors().second ) );
  d1->partons( make_pair( d1->partons().first, p ) );
  d2->partons( make_pair( p, d2->partons().second ) );
  p->dipoles( make_pair( d1, d2 ) );

  theChildren.second = d2;
}

void Dipole::emit() {
  static DebugItem trace("DIPSY::Trace", 9);
  if ( tDipolePtr sd = swingDipole() ) {
    if ( trace ) cerr << "Swing " << tag() << " & " << sd->tag();
    recombine();
    if ( trace ) cerr << " -> " << children().first->tag()
		      << " & " << sd->children().first->tag() << endl;
    return;
  }

  dipoleState().handler().emitter().emit(*this);
  //  splitDipole();

  // *** ATTENTION *** implement this.
}

void Dipole::splitDipole(double colsel) {
  if ( !generatedGluon() )
    Throw<NothingGenerated>()
      << "Tried to perform an emission where none was generated."
      << Exception::abortnow;
  theChildren = make_pair(dipoleState().createDipole(),
			  dipoleState().createDipole());
  thePartons.first->theChildren.insert(generatedGluon());
  thePartons.second->theChildren.insert(generatedGluon());

  if ( neighbors().first ) {
    children().first->theNeighbors.first = neighbors().first;
    neighbors().first->theNeighbors.second = children().first;
    neighbors().first->reset();
    neighbors().first->touch();
  }
  children().first->theNeighbors.second = children().second;
  children().second->theNeighbors.first = children().first;
  partons().first->theDipoles.second = children().first;
  partons().second->theDipoles.first = children().second;
  if ( neighbors().second ) {
    children().second->theNeighbors.second = neighbors().second;
    neighbors().second->theNeighbors.first = children().second;
    neighbors().second->reset();
    neighbors().second->touch();
  }
  children().first->partons(make_pair(partons().first,generatedGluon()));
  children().second->partons(make_pair(generatedGluon(),partons().second));
  generatedGluon()->parents(partons());
  generatedGluon()->dipoles(children());
  generatedGluon()->rightMoving( partons().first->rightMoving() );

  if ( UseRandom::rndbool(colsel) ) {
    children().first->theColour = theColour;
    dipoleState().generateColourIndex(children().second);
    generatedGluon()->mainParent(partons().second);
  } else {
    children().second->theColour = theColour;
    dipoleState().generateColourIndex(children().first);
    generatedGluon()->mainParent(partons().first);
  }
}

int Dipole::colourSystem() const {
  return colour()/Current<DipoleEventHandler>()->nColours();
}

void Dipole::colourSystem(int sys) {
  colour(colour()%Current<DipoleEventHandler>()->nColours() +
	 sys*Current<DipoleEventHandler>()->nColours());
}

tEffectivePartonPtr Dipole::getEff(tcPartonPtr p, InvEnergy range) const {
  tEffectivePartonPtr ret;
  if ( Current<DipoleEventHandler>()->effectivePartonMode() < 0 )
    Throw<Exception>()
      << "Tried to use effective partons when we want the shadow strategy!"
      << Exception::abortnow;
  if ( Current<DipoleEventHandler>()->effectivePartonMode() == 3 ) {
    if ( p == partons().first ) {
      effectivePartons().first->setRange(range);
      return effectivePartons().first;
    }
    if ( p == partons().second ) {
      effectivePartons().second->setRange(range);
      return effectivePartons().second;
    }
  }
  if ( !effmap ) buildEffMap();
  if ( p == partons().first ) {
    EffMap::const_iterator it = effmap->first.begin();
    ret = it->second;
    while ( ++it != effmap->first.end() && it->first < range ) ret = it->second;
  }
  if ( p == partons().second ) {
    EffMap::const_iterator it = effmap->second.begin();
    ret = it->second;
    while ( ++it != effmap->second.end() && it->first < range ) ret = it->second;
  }
  return ret;
}

void Dipole::buildEffMap() const {
  if ( effmap ) delete effmap;
  effmap = new pair<EffMap,EffMap>();
  // First take the first partonand create the largest one immaginable.
  EffectivePartonPtr ep =
    EffectiveParton::create(*partons().first, size()/2.0);
  // Then we loop through the included partons in order and create all
  // sub effective partons.
  set<tPartonPtr> intp = ep->internalPartons();
  set<InvEnergy> dists;
  if ( intp.empty() ) {
    const vector<tPartonPtr> & intp = ep->cachedPartons();
    for ( vector<tPartonPtr>::const_iterator it = intp.begin();
	  it != intp.end(); ++it )
      dists.insert(sqrt(ep->originalParton()->dist2(**it)));
  } else {
    for ( set<tPartonPtr>::iterator it = intp.begin(); it != intp.end(); ++it )
      dists.insert(sqrt(ep->originalParton()->dist2(**it)));
  }
  effmap->first.resize(intp.size());
  for ( set<InvEnergy>::iterator it = dists.begin();
	it != dists.end(); ++it ) {
    effmap->first.push_back
      (make_pair(*it, EffectiveParton::create(*partons().first, *it)));
  }
  // Then do the same for the second
  ep = EffectiveParton::create(*partons().second, size()/2.0);
  // Then we loop through the included partons in order and create all
  // sub effective partons.
  intp = ep->internalPartons();
  dists.clear();
  if ( intp.empty() ) {
    const vector<tPartonPtr> & intp = ep->cachedPartons();
    for ( vector<tPartonPtr>::const_iterator it = intp.begin();
	  it != intp.end(); ++it )
      dists.insert(sqrt(ep->originalParton()->dist2(**it)));
  } else {
    for ( set<tPartonPtr>::iterator it = intp.begin(); it != intp.end(); ++it )
      dists.insert(sqrt(ep->originalParton()->dist2(**it)));
  }
  effmap->second.reserve(intp.size());
  for ( set<InvEnergy>::iterator it = dists.begin(); it != dists.end(); ++it ) {
    effmap->second.push_back
      (make_pair(*it, EffectiveParton::create(*partons().second, *it)));
  }
}

void Dipole::persistentOutput(PersistentOStream & os) const {
  os << theDipoleState << thePartons << theNeighbors << theChildren
     << theGeneratedGluon << theSwingDipole << theGeneratedY << theColour
     << theInteracted << isParticipating << isRecoilSwing << isTouched
     << isOn << ounit(swingCache, 1.0/GeV2);
}

void Dipole::persistentInput(PersistentIStream & is, int) {
  if ( effmap ) delete effmap;
  is >> theDipoleState >> thePartons >> theNeighbors >> theChildren
     >> theGeneratedGluon >> theSwingDipole >> theGeneratedY >> theColour
     >> theInteracted >> isParticipating >> isRecoilSwing >> isTouched
     >> isOn >> iunit(swingCache, 1.0/GeV2);
}

string Dipole::tag() const {
  static ObjectIndexer<int, const Dipole> di;
  di(tcDipolePtr());
  static ObjectIndexer<int, const Parton> pi;
  pi(tcPartonPtr());
  static ObjectIndexer<int, const ShadowParton> si;
  si(tcSPartonPtr());
  ostringstream os;
  os << "D" << di(this) << "(P" << pi(partons().first);
  if ( partons().first->shadow() ) os << "[S" << si(partons().first->shadow()) << "]";
  os << ",P" << pi(partons().second);
  if ( partons().second->shadow() ) os << "[S" << si(partons().second->shadow()) << "]";
  os << ")";
  return os.str();
}

// Static variable needed for the type description system in ThePEG.
#include "ThePEG/Utilities/DescribeClass.h"
DescribeClass<Dipole,PersistentBase>
  describeDIPSYDipole("DIPSY::Dipole", "libAriadne5.so libDIPSY.so");

void Dipole::Init() {}

