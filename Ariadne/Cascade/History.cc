// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the History class.
//

#include "History.h"
#include "AriadneHandler.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Utilities/Current.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#define pnew(clss,expr) Ptr<clss>::pointer(new clss expr)


using namespace Ariadne5;

History::History() {}

History::History(int depth, DipoleStatePtr statein)
  : theState(statein), mother(tHistoryPtr()), emission(tcEmPtr()),
    foundOrderedPath(false), foundCompletePath(false),
    theScale(sqrt(statein->totalMomentum().m2())), prob(1.0) {
  recurse(depth, true);
}

History::
History(int depth, Energy scalein, DipoleStatePtr statein,
	bool isOrdered, double probin, tHistoryPtr mothin, tcEmPtr emin)
  : theState(statein), mother(mothin), emission(emin), foundOrderedPath(false),
    foundCompletePath(false), theScale(scalein), prob(probin) {
  recurse(depth, isOrdered);
}

History::~History() {}

void History::recurse(int depth, bool isOrdered) {

  set<EmPtr,ScaleCmp> clusterings;

  // If this is not the fully clustered state, try to find possible
  // clusterings.
  if ( depth > 0 ) clusterings = getClusterings();

  // If no clusterings were found, the recursion is done and we
  // register this node.
  if ( clusterings.empty() ) {
    Energy muF = Current<AriadneHandler>()->checkBornState(*state());
    // If we should have reached the Born-level, but no reasonable
    // state could be identified, we simply skip this history.
    if ( depth == 0 && muF < ZERO ) return;
    registerPath(tHistoryPtr(this), abs(muF), isOrdered, depth == 0);
    return;
  }
    
  for ( set<EmPtr,ScaleCmp>::iterator it = clusterings.begin();
	it != clusterings.end(); ++it ) {
    EmPtr em = *it;
    bool ordered = isOrdered;
    if ( !ordered || ( mother && em->rho < mother->scale() ) ) {
      if ( onlyOrderedPaths()  ) continue;
      ordered = false;
    }

    // Perform the clustering and recurse and construct the next
    // history node.
    DipoleStatePtr newstate = cluster(em);
    if ( newstate ) children.push_back
    		      (HistoryPtr(new History(depth - 1, em->rho, newstate,
    					      ordered, prob*em->prob,
    					      tHistoryPtr(this), em)));
  }
  // Check unlikely case of all reconstructions failing for some
  // reason, n which we register this as an incomplete path.
  if ( children.empty() )
    registerPath(tHistoryPtr(this),
		 abs(Current<AriadneHandler>()->checkBornState(*state())),
		 isOrdered, false);
}

tHistoryPtr History::select(Reordering reor) {
  sel = paths[UseRandom::rnd()];
  if ( !foundOrderedPath ) reorder(reor);
  theScale = scales[sel];
  paths.clear();
  scales.clear();
  sel->clean();
  return sel;
}

void History::clean(HistoryPtr h) {
  if ( h ) {
    children.clear();
    children.push_back(h);
  }
  if ( mother ) mother->clean(HistoryPtr(this));
}

double History::weight(double as0, const AriadneHandler & hdl) const {
  if ( !mother ) return 1.0;
  double w = hdl.alphaS(sqr(scale()))/as0;

  for ( int i = 0, N = state()->remnants().first.size(); i < N; ++i ) {
    if ( mother->state()->remnants().first[i]->hard() ) continue;
    w *= mother->state()->remnants().first[i]->xfx(sqr(scale()));
    w /= state()->remnants().first[i]->xfx(sqr(scale()));
  }
  for ( int i = 0, N = state()->remnants().second.size(); i < N; ++i ) {
    if ( mother->state()->remnants().second[i]->hard() ) continue;
    w *= mother->state()->remnants().second[i]->xfx(sqr(scale()));
    w /= state()->remnants().second[i]->xfx(sqr(scale()));
  }
    
  return w*mother->weight(as0, hdl);
}

double History::weight(double as0, Energy2 muF2) const {
  const AriadneHandler & hdl = Current<AriadneHandler>::current();
  if ( !sel ) return 0.0;
  double w = sel->trialShower(hdl.startingScale(*sel->state()));
  if ( w <= 0.0 ) return 0.0;

  for ( int i = 0, N = state()->remnants().first.size(); i < N; ++i ) {
    tRemParPtr rnew = sel->state()->remnants().first[i];
    if ( rnew->hard() ) continue;
    PDF pdf = hdl.pdf<PDF>(rnew->originalExtracted());
    w *= pdf.xfx(tcPDPtr(&rnew->extractedData()), sqr(scale()), rnew->x());
    tRemParPtr rold = state()->remnants().first[i];
    w /= pdf.xfx(tcPDPtr(&rold->extractedData()), muF2, rold->x());
  }
    
  for ( int i = 0, N = state()->remnants().second.size(); i < N; ++i ) {
    tRemParPtr rnew = sel->state()->remnants().second[i];
    if ( rnew->hard() ) continue;
    PDF pdf = hdl.pdf<PDF>(rnew->originalExtracted());
    w *= pdf.xfx(tcPDPtr(&rnew->extractedData()), sqr(scale()), rnew->x());
    tRemParPtr rold = state()->remnants().second[i];
    w /= pdf.xfx(tcPDPtr(&rold->extractedData()), muF2, rold->x());
  }

  return w*sel->weight(as0, hdl);
}

bool History::onlyOrderedPaths() const {
    if ( !mother || foundOrderedPath ) return foundOrderedPath;
    return  foundOrderedPath = mother->onlyOrderedPaths();
}

bool History::registerPath(tHistoryPtr h, Energy rho0,
			   bool isOrdered, bool isComplete) {
    // We are not interested in improbable paths.
    if ( h->prob <= 0.0) return false;

    // We only register paths in the initial node.
    if ( mother ) return mother->registerPath(h, rho0, isOrdered, isComplete);

    if ( foundOrderedPath && !isOrdered ) return false;
    if ( foundCompletePath && !isComplete ) return false;
    if ( isOrdered && isComplete ) {
      // If this is the first complete, ordered path, discard the
      // old, non-ordered or incomplete ones.
      if ( !foundOrderedPath || !foundCompletePath ) paths.clear();
      foundOrderedPath = true;
      foundCompletePath = true;
    }
    else if ( isComplete ) {
      // If this is the first complete path, discard the old,
      // incomplete ones.
      if ( !foundCompletePath )	paths.clear();
      foundCompletePath = true;
    }
    paths.insert(h->prob, h);
    scales[h] = rho0;

    return true;
}

set<EmPtr,History::ScaleCmp> History::getClusterings() {
  list<EmPtr> clusterings;
  const vector<EmitterPtr> & potential = Current<AriadneHandler>()->emitters();
  for ( int i = potential.size() - 1; i >= 0; --i ) {
    if ( !potential[i]->reconstructor() ) continue;
    list<EmPtr>::iterator it = clusterings.begin();
    while ( it != clusterings.end() ) {
      if ( potential[i]->overrideInverse(**it) )
	it = clusterings.erase(it);
      else
	++it;
    }
    vector<EmPtr> clusts = potential[i]->inverseEmissions(*state());
    clusterings.insert(clusterings.end(), clusts.begin(), clusts.end());
  }

  return set<EmPtr,ScaleCmp>(clusterings.begin(), clusterings.end());

}

DipoleStatePtr History::cluster(tEmPtr emission) {
  DipoleState::TranslationMap trans;
  DipoleStatePtr newstate = state()->preclone(trans);
  emission->rebind(trans);
  newstate->postclone(trans);
  if ( emission->model->performInverse(*emission, *newstate) ) return newstate;
  else return DipoleStatePtr();
}

void History::reorder(Reordering reor) {
  vector<tHistoryPtr> seq;
  tHistoryPtr h = sel;
  while ( h->mother ) {
    seq.push_back(h);
    h = h->mother;
  }
  if ( reor == UseMax ) {
    Energy scalemax = ZERO;
    for ( int i = seq.size() - 1; i >= 0; --i )
      if ( seq[i]->scale() < scalemax ) seq[i]->theScale = scalemax;
      else scalemax = seq[i]->scale();
  }
  else if ( reor == UseMin ) {
    Energy scalemin = sel->scale();
    for ( int i = 1, N = seq.size(); i < N; ++i )
      if ( seq[i]->scale() > scalemin ) seq[i]->theScale = scalemin;
      else scalemin = seq[i]->scale();
  }
}

double History::trialShower(Energy maxscale) const {
  Energy rhomax = maxscale;
  Energy rhomin = scale();
  while ( rhomax > rhomin && 
	  ( rhomax = state()->select(rhomin, rhomax) ) > rhomin ) {
    SaveDipoleState backup(state());
    if ( state()->perform() ) return 0.0;
    theState = backup.revert();
    state()->selected()->dipole->touch();
  }
  if ( mother ) return mother->trialShower(scale());
  return 1.0;
}

void History::persistentOutput(PersistentOStream & os) const {
  os << theState << mother << emission << children << paths
     << ounit(scales, GeV) << foundOrderedPath << foundCompletePath
     << ounit(theScale, GeV) << prob << sel;
}

void History::persistentInput(PersistentIStream & is, int) {
  is >> theState >> mother >> emission >> children >> paths
     >> iunit(scales, GeV) >> foundOrderedPath >> foundCompletePath
     >> iunit(theScale, GeV) >> prob >> sel;
}


// The following static variable is needed for the type description
// system in ThePEG.
DescribeClass<History,Base>
  describeAriadne5History("Ariadne5::History", "libAriadne5.so");

void History::Init() {}

