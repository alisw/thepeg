// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DiffractiveEventFiller class.
//

#include "DipoleEventHandler.h"
#include "EventFiller.h"
#include "DiffractiveEventFiller.h"
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

#include <iostream>
#include <fstream>

using namespace DIPSY;

DiffractiveEventFiller::DiffractiveEventFiller()
  {}

DiffractiveEventFiller::~DiffractiveEventFiller() {}

IBPtr DiffractiveEventFiller::clone() const {
  return new_ptr(*this);
}

IBPtr DiffractiveEventFiller::fullclone() const {
  return new_ptr(*this);
}

/*
 * Fills the step with a diffractive event.
 */
double DiffractiveEventFiller::fill(Step & step, DipoleEventHandler & eh, tPPair inc,
			 DipoleState & dr, DipoleState & de,
			 const ImpactParameters & b) const {
  //Some debugging set up.
  static DebugItem printsteps("DIPSY::PrintSteps", 6);
  static DebugItem printstepdetails("DIPSY::PrintStepDetails", 6);

  //Sets debugging output.
  //TODO: interface
  // bool printstepdetails = true;
  if ( printstepdetails ) cout << setprecision(6);
  if ( printstepdetails ) cout << "entered fill()" << endl;

  //evolve the valence state to maximum allowed rapidity for the excited state.
  //turn swing of for this. then revert for the virtual cascades later.
  double lambda = eh.swingPtr()->lambda();
  //TODO: can this be set to 0?
  eh.swingPtr()->lambda(0.0000000000001);
  dr.evolve(dr.lowestY(), maxY());
  eh.swingPtr()->lambda(lambda);

  if ( printstepdetails ) cout << "evolved real state" << endl;

  //to get the correct weight for the real cascade (specially the absence of
  //a "no more emissions after this" term for the real cascade), we need to
  //return an event for every subcascade. Or equivalently, and avoiding correlations
  //between the final states, returning one random subcascade weighted by the number
  //of subcascades.
  //
  //The selectionProbability is only used when running with nInteractions == -1,
  //and should be set to about 1/(the largest number of subcascades you expect in a
  //real cascade). Too low means that it will run through a lot of real cascades before
  //picking one (not too bad), and too high will give occasional events with very high
  //weight (potentially making the result untrustworthy).
  double selectionProbability = 1;
  //the weight the event receives in association with the choice of subcascade.
  double choiceWeight = 0.0;

  //this choice runs with any number of interactions. While producing an exact result,
  //it converges incredibly slowly, as the cross section will be dominated by a very small
  //fraction of the events (namely the ones with a single interacting large dipole).
  if ( nInteractions() < 0 ) {
    choiceWeight = selectSubcascade(dr, selectionProbability);
  }
  //No interacting dipoles means an elastic event.
  else if ( nInteractions() == 0 ) {
    pickSubcascade(dr, 0);
    increaseNSubCascades(1.0);
    choiceWeight = 1.0;
  }
  //A single interacting dipole.
  else if ( nInteractions() == 1 ) {
    choiceWeight = selectSingleSubcascade(dr);
  }
  //a fix number of interaction larger than 1.
  else if ( nInteractions() > 1 ) {
    choiceWeight = selectNSubcascade(dr, nInteractions());
  }
  //if no subcascade chosen, start over.
  if ( choiceWeight == 0.0 ) {
    if ( printstepdetails ) cout << "choiceweight is 0, return 0 and start over." << endl;
    return 0.0;
  }
  double weight = choiceWeight;

  //identify the ends of the chains. These are treated differently than the rest in the
  //amplitude, so it is easier to identify them one time for all here. These are the circled
  //partons in fig 15 in the diffractive paper.
  vector<PartonPtr> chainCaps = childlessChildren(&dr);

  if ( printstepdetails ) cout << "found chaincaps " << chainCaps.size() << endl;

  //We now clean up the subcascade into a proper real finalstate.
  //Ie, do all the things you have to do to a virtual cascade to make it real,
  //like reweight the pt maxima, double check ordering etc.
  makeReal(dr, chainCaps, de);

  //The cleaning can mess up the interacting partons, although it shouldn't be very frequent.
  //but better to recheck and avoid troubles later on.
  chainCaps = childlessChildren(&dr);

  //check that the number of chaincaps (ie, interactions) in the real state is still the correct one.
  //they may go away if some interacting parton was not ordered etc.
  //if we run all numbers of interaction (-1 or -2), then it doesn't matter and we should always continue.
  if ( int(chainCaps.size()) != nInteractions() && nInteractions() > -1 ) {
    if ( printstepdetails ) cout << "got " << chainCaps.size() << " chainCaps, veto." << endl;
    return 0.0;
  }

  //if the subcascade doesnt reach up to minimum y (too low M_X^2), start over.
  if ( dr.highestY() < minY() ) {
    if ( printstepdetails ) cout << "too low highestY, return 0 and start over. needed "
				 << minY() << ", got " << dr.highestY() << endl;
    return 0.0;
  }

  //estimate the amplitude and veto/reweight accordingly
  //this is to avoid running through all the elastic cascade for all excited states
  //that will get a very low amplitude anyways.
  //
  //I have some doubts on this method, as if the estimate is not accurate for all cascades
  //it may throw away high-amplitude cascade often, and give some events with very high weight.
  //However, tests has shown that the weight distribution looks a bit better with reweighting
  //like this, so I left it in. But keep track on the weight distribution and consider
  //switching this off (-2 interactions) if it gets worse.
  double estimate = reweightCascade(dr, chainCaps);
  if ( nInteractions() != -2 ) {
    if ( estimate == 0.0 ) {
      if ( printstepdetails ) cout << "estimated amplitude 0, probably not good." << endl;
      return 0.0;
    }

    //to estimate weight, we need to square the estimate of the amplitude.
    double passProb = sqr(estimate);

    //to tone down the reweighting a bit (specially to avoid spikes to high weights)
    //the square root of the estimated weight is used.
    passProb = sqrt(passProb);

    //throw away most of the cascades that seems to have low amplitude, but increase the weight
    //of the ones we do decide to run.
    if ( UseRandom::rnd() > passProb ) {
      return 0.0;
    }
    else weight /= passProb;
    //TODO: proper error message
    if ( isnan(weight) ) cout << "weight is not a number after interaction estimate." << endl;
    
  }


  //calculate the amplitude for this real cascade, by summing over virtual cascades.
  double amp = amplitude(dr, chainCaps, b);

  //probability is proportional to the square of the amplitude.
  weight *= sqr(amp);

  //TODO: proper error message
  if ( isnan(weight) ) cout << "weight is not a number after amplitude." << endl;

  if ( printstepdetails ) cout << "calculated amp^2 " << sqr(amp) << endl;
  if ( printstepdetails ) cout << "weight is now " << weight << endl;
  if ( printstepdetails ) cout << "number of gluons: " << dr.getPartons().size() << endl;

  //balance momenta by taking some energy from the elastic state to put
  //the excited state on shell. If this fails (usually due to not enough energy),
  //then throw the event away.
  if ( !balanceMomenta(&dr, &de) ) {
    if ( printstepdetails ) cout << "failed to balance momenta" << endl;
    return 0.0;
  }

  //add finalstate swing
  double yrange = FSSwingTime();
  double ystep = FSSwingTimeStep();
  for ( double y = 0.0; y < yrange; y += ystep ) {
    dr.swingFS(y, y + ystep);
  }

  //extract strings from the excited state.
  vector<String> strings = dr.strings();

  //mirrer elastic state to make it move in the opposite direction of the
  //excited state and turn it into a proton
  //This is will obviously have to be changed if diffracting on something
  //else than a proton.
  de.mirror(0.0);
  DipoleState::String protonString = makeProton(& de);

  //add elastic particle to the strings.
  strings.push_back(protonString);


  //fill step. Scrap event if something doesn't work.
  if ( strings.empty() || ! fillStep(step, inc, strings) ) {
    weight = 0.0;
    currentWeight = 0.0;
    //TODO: proper error message
    if ( printstepdetails ) cout << "empty string or failed fillstep, return 0" << endl;
    return weight;
  }

  if ( printstepdetails ) dr.diagnosis(true);

  //For really long runs, one could imagine that the elastic cascades should be
  //regenerated to reduce the error. Due to limited RAM on the computer
  //one cannot increase the number of elastic events to eb stored over a certain
  //number. Didn't turn out to really be needed yet though, thus commented out.
  // int CEN = generator()->currentEventNumber();
  // int NEC = nElasticCascades();
  // if ( double(CEN/NEC) == double(CEN)/double(NEC) ) {
  //   cout << "generating new elastic cascades" << endl;
  //   DipoleEventHandlerPtr eh = Current<DipoleEventHandler>().ptr();
  //   DipoleStatePtr dummyState = eh->WFR().generate(*eh, 1000.0*GeV);
  //   Energy W = dummyState->handler().lumiFn().maximumCMEnergy();
  //   Energy2 a = eh->WFR().m2() - eh->WFL().m2() + sqr(W);
  //   Energy PR = (a + sqrt(sqr(a) - 4.0*eh->WFR().m2()*sqr(W)))*0.5/W;
  //   initialiseVirtualCascades(*eh, eh->WFR(), PR, -maxY());
  //   cout << "done generating new elastic cascades" << endl;
  // }
  
  if ( printstepdetails )  cout << "done in DiffractiveEventFiller::fill(...), leaving." << endl;
  return weight;
}

/*
 * Some comments on this in the writeup. It is based on the version in the ND collisions.
 */
bool DiffractiveEventFiller::balanceMomenta(DipoleStatePtr excited, DipoleStatePtr elastic) const {
  //First calculate the p+ and p- of the excited state.
  Energy intPlus1 = ZERO;
  Energy intMinus1 = -excited->minus();
  list<PartonPtr> partons = excited->getPartons();
  for ( list<PartonPtr>::iterator it = partons.begin(); it != partons.end(); it++) {
    intPlus1 += (*it)->plus();
    intMinus1 += (*it)->minus();
  }

  //The p+ and p- of the elastic state. Reverse the + and -, due to not being mirrored yet, ie, both
  //the excited and elastic states have large p+ and low p- right now.
  Energy intPlus2 = -elastic->minus();
  //TODO: should it really be "excited" here, not "elastic"? Compare to ND.
  Energy intMinus2 = excited->plus();

  //solve the equation for:
  //x is the fraction of the p+ the excited state has to give to the elastic state.
  //y is the fraction of p- that it has to give to the excited state. Note however, that the elastic
  //state is not yet mirrored, so the roles of p+ and p- should be reversed.

  //Just set up some paramaters that apear in the solution.
  Energy2 A = intPlus1*intMinus2;
  Energy2 B = intPlus1*(intMinus1 - intMinus2) - intPlus2*intMinus2;
  Energy2 C = intPlus2*intMinus2;

  //This is what we have to take the square root of, so check it is positive. If not enough
  //energy, there will be no real solutions to get on shell.
  //TODO: proper error message.
  if ( sqr(B/(2*A)) - C/A < 0.0 ) {
    return false;
  }

  //The solution we are ineterested in. It is a 2:nd degree eq, but I am pretty sure this is the
  //solution we want. The other solution will make the two states bounce out almost the way they came
  //rather than just glancing of, as they should do (low-x, etc). The other has "- sqrt" instead.
  double y = -B/(2*A) + sqrt(sqr(B/(2*A)) - C/A);
  double x = 1.0 - intPlus2/(y*intPlus1);

  //this should not happen, but...
  if ( isnan(x) || isnan(y) || x <= 0.0 || y <= 0.0 ) {
    //TODO: proper error message.
    cout << "found negative or nan boosts in DifractiveEventFiller." << endl;
    return false;
  }

  //With x and y found, change the 4-momenta.
  for ( list<PartonPtr>::iterator it = partons.begin(); it != partons.end(); it++) {
    (*it)->plus(x*(*it)->plus());
    (*it)->updateYMinus();
  }

  //Again, elastic state not mirrored yet.
  elastic->plus(y*elastic->plus());
  elastic->minus(elastic->minus()/y);
  return true;
}

/*
 * Had this in to calculate the inclusive cross section but I guess it should be called from
 * somewhere else.
 * TODO: is this done somewhere else? Comment if needed.
 */
void DiffractiveEventFiller::calculateInclusive(DipoleStatePtr real, CrossSection weight,
						const ImpactParameters & b) const {
  DipoleEventHandlerPtr eh = Current<DipoleEventHandler>().ptr();
  double prob = 0.0;
  for (int i = 0; i < int(virtualCascades.size()); i++  ) {
    double F = eh->xSecFn().sumf(*real, *virtualCascades[i], b);
    prob += eh->xSecFn().unitarize(F);
  }
  prob /= double(virtualCascades.size());
  totalXSec += weight*sqr(prob);
  totalSqr += sqr(weight*sqr(prob));
  neve++;
}

/*
 * 
 */
double DiffractiveEventFiller::selectSubcascade(DipoleState & dr,
					      double selectionProbability) const {
  static DebugItem printstepdetails("DIPSY::PrintStepDetails", 6);
  if ( printstepdetails ) cout << "entering DiffractiveEventFiller::selectSubcascade(...)" << endl;

  //count the elastic subcascade this many times extra
  int elasticBoost = 10;  //interface! tested that results do not depend on this variable
  //which gives that many extra subcascades effectively
  int nSub = nSubcascades(dr) + elasticBoost;
  if ( printstepdetails ) cout << "number of subcascade choices: " << nSub << endl;

  //choose if using this cascade, nsub*selProb
  double weight = nSub*selectionProbability;
  if ( printstepdetails ) cout << "weight is " << weight << endl;
  if ( isnan(weight) ) cout << "weight is not a number after cascade selection." << endl;
  double R = UseRandom::rnd();
  if ( R > weight ) {
    if ( printstepdetails ) cout << "veto, try again with another real cascade" << endl;
    return 0.0;
  }
  //if passed veto with a 1 or less prob, set prob to 1
  //larger probability than 1 has to go into the returned weight
  if ( weight < 1.0 ) weight = 1.0;
  if ( weight > 1.0 ) R *= weight;
  if ( printstepdetails ) cout << "passed selection veto, weight set to " << weight << endl;

  //select which subcascade using an integer between 0 and nSub
  int choice = int(R/selectionProbability);
  if ( printstepdetails ) cout << "pick subcascade no " << choice << endl;

  //identify it and modify dr accordingly
  //first check if it was one of the elastic ones, otherwise choose one normally
  //note that a choice of 0 corresponds to the elastic cascade
  if ( choice < elasticBoost + 1 ) {
    choice = 0;
    weight /= double(elasticBoost) + 1.0;
    if ( printstepdetails ) cout << "pick elastic, weight divided by " << double(elasticBoost) + 1.0 << ", now " << weight << endl;
  }
  else choice -= elasticBoost;
  if ( printstepdetails ) cout << "choice is " << choice << endl;
  if ( isnan(weight) ) cout << "weight is not a number after elastic boost" << endl;
  weight *= pickSubcascade(dr, choice);

  return weight;
}

void DiffractiveEventFiller::addPartonsToSet(DipolePtr d, set<PartonPtr> x) const {
  x.insert(d->partons().first);
  x.insert(d->partons().second);
}

double DiffractiveEventFiller::addRandomChain(DipoleState & dr) const {
  static DebugItem printstepdetails("DIPSY::PrintStepDetails", 6);
  if ( printstepdetails ) cout << "entering DiffractiveEventFiller::selectSubcascade(...)" << endl;


  //first find the old partons and dipoles in a set
  list<PartonPtr> partons = dr.getPartons();
  set<PartonPtr> oldP;
  for ( list<PartonPtr>::iterator it = partons.begin();
	it != partons.end(); it++ )
    oldP.insert(*it);

  //find the set of all partons
  set<PartonPtr> allP;
  set<DipolePtr> allD = dr.getAllDipoles();
  for ( set<DipolePtr>::iterator it = allD.begin();
	it != allD.end(); it++ )
    addPartonsToSet(*it, allP);

  //remove the old partons
  for ( set<PartonPtr>::iterator it = oldP.begin();
	it != oldP.end(); it++ )
    allP.erase(*it);

  set<PartonPtr> newP = allP;

  if ( newP.empty() ) return 0.0;

  //pick one of the leftover partons as the new chain end (ie interaction)
  int choice = UseRandom::rnd()*newP.size();

  set<PartonPtr>::iterator it = newP.begin();
  while (choice > 0 ) {
    it++;
    choice--;
  }
  // addChainFromInteraction(dr,*it);

  return 0.0;
}

double DiffractiveEventFiller::selectSingleSubcascade(DipoleState & dr) const {
  static DebugItem printstepdetails("DIPSY::PrintStepDetails", 6);
  if ( printstepdetails ) cout << "entering DiffractiveEventFiller::selectSubcascade(...)" << endl;

  list<PartonPtr> partons = dr.getPartons();
  int nSub = countSingleStates(dr);
  if ( nSub == 0 ) return 0.0;
  if ( printstepdetails ) cout << "number of subcascade choices: " << nSub << endl;

  //count the number of possible subcascades.
  increaseNSubCascades(double(nSub));


  //first find the valence partons
  vector<DipolePtr> dips = dr.initialDipoles();
  set<PartonPtr> valP;
  for ( int i = 0; i < int(dips.size());i++) {
    valP.insert(dips[i]->partons().first);
    valP.insert(dips[i]->partons().second);
  }

  //remove the valence partons from the set of partons
  for ( list<PartonPtr>::iterator it = partons.begin(); it != partons.end(); it++ ) {
    PartonPtr p = *it;
    if ( valP.find(p) != valP.end() ) {
      it = partons.erase(it);
      it--;
    }
  }

  //pick one of the partons as the last emission (ie interaction)
  int choice = UseRandom::rnd()*partons.size();

  list<PartonPtr>::iterator it = partons.begin();
  while (choice > 0 ) {
    it++;
    choice--;
  }
  pickSubcascadeFromInteraction(dr,*it);

  return nSub;
}

double DiffractiveEventFiller::selectNSubcascade(DipoleState & dr, int N) const {
  static DebugItem printstepdetails("DIPSY::PrintStepDetails", 6);
  if ( printstepdetails ) cout << "entering DiffractiveEventFiller::selectSubcascade(...)" << endl;

  if ( N < 1 ) cerr << "DiffractiveEventFiller::selectNSubcascade called with less than one interaction" << endl;

  double ret = selectSingleSubcascade(dr);

  int n = 1;
  while ( n < N ) {
    ret*= addRandomChain(dr);
  }

  //count the number of possible subcascades.
  increaseNSubCascades(double(ret));

  return ret;
}

int DiffractiveEventFiller::countSingleStates(DipoleState & dr) const {
  list<PartonPtr> partons = dr.getPartons();
  //count the valence partons
  vector<DipolePtr> dips = dr.initialDipoles();
  set<PartonPtr> valP;
  for ( int i = 0; i < int(dips.size());i++) {
    valP.insert(dips[i]->partons().first);
    valP.insert(dips[i]->partons().second);
  }
  //return the number of excited states (partons - valP)
  return partons.size() - valP.size();
}

int DiffractiveEventFiller::nSubcascades(DipoleState & dl) const {
  static DebugItem printstepdetails("DIPSY::PrintStepDetails", 6);
  if ( printstepdetails ) cout << "entering DiffractiveEventFiller::nSubcascades(...)" << endl;
  int ret = 1;
  //iterate through valence dipoles and multiply number of combinations
  vector<DipolePtr> valDips = dl.initialDipoles();

  if ( printstepdetails ) cout << "looking for subcascades of state:" << endl;
  // dl.plotState(true);

  for ( int i = 0; i < int(valDips.size()); i++ ) {
    ret *= nSubcascades(valDips[i]);
  }

  if ( printstepdetails ) cout << "found " << ret << " subcascsades" << endl;

  return ret;
}

int DiffractiveEventFiller::nSubcascades(DipolePtr dip) const {
  static DebugItem printstepdetails("DIPSY::PrintStepDetails", 6);
  //a = 1 + b*c if 2 kids
  //a = b if one kid
  //a = 1 if no kids
  if ( dip->children().first && !dip->children().second ) {
    if ( printstepdetails ) cout << "found swung dipole in nsubcascades..." << endl;
    return nSubcascades(dip->children().first);
  }
  else if ( !dip->children().first && dip->children().second ) {
    return nSubcascades(dip->children().second);
  }
  else if ( dip->children().first && dip->children().second ) {
    return 1 + nSubcascades(dip->children().second)*nSubcascades(dip->children().first);
  }
  else {
    return 1;
  }
}

double DiffractiveEventFiller::pickSubcascade(DipoleState & dr, int choice) const {
  static DebugItem printstepdetails("DIPSY::PrintStepDetails", 6);
  int configuration = choice;
  double weight = 1.0;
  //set all partons off shell
  setOffShell(dr);

  //mark the ones in the subcascade as on shell, by looping over the valence dips
  //this MAY have to be corrected with a swing.
  vector<DipolePtr> valDips = dr.initialDipoles();
  for ( int i = 0; i < int(valDips.size()); i++ ) {
    //valence dips always real
    setPartonsOnShell(valDips[i]);

    //find the configuration number for this valence dipole
    //note that the max configuration number is the product of
    //the number of combinations in each valence dipole
    int N = nSubcascades(valDips[i]);
    int subconfiguration = configuration % N;
    //the remaining valence dips can choose between the config/N combinations.
    configuration /= N;

    if ( printstepdetails ) cout << "for this valence dip, use subconfiguration " << subconfiguration
	 << " out of " << N << endl;
    //now look at this valence dipole.
    weight *= pickSubcascade(valDips[i], subconfiguration);
  }

  if ( printstepdetails ) cout << "done tagging onshells. now remove virtuals" << endl;
  // dr.plotState(true);

  //remove the off-shell ones
  removeVirtuals(& dr);
  if ( printstepdetails ) cout << "removed virtuals" << endl;

  //fix the parent structure that may get broken (REALLY?) when removing the virtuals
  //this is so that the real state will not run into troubles
  // updateParentStructure(&dr);
  return weight;
}

void DiffractiveEventFiller::pickSubcascadeFromInteraction(DipoleState & dr, PartonPtr p) const {
  static DebugItem printstepdetails("DIPSY::PrintStepDetails", 6);
  //set all partons off shell
  setOffShell(dr);

  setValenceOnShell(dr);

  recursiveSetOnShell(p);

  //remove the off-shell ones
  removeVirtuals(& dr);
  if ( printstepdetails ) cout << "removed virtuals" << endl;

  //fix the parent structure that may get broken (REALLY?) when removing the virtuals
  //this is so that the real state will not run into troubles
  // updateParentStructure(&dr);
}


void DiffractiveEventFiller::recursiveSetOnShell(PartonPtr p) const {
  p->onShell(true);
  if ( p->parents().first ) recursiveSetOnShell(p->parents().first);
  if ( p->parents().second ) recursiveSetOnShell(p->parents().second);
}

  double DiffractiveEventFiller::reweightCascade(DipoleState & dr,
						 vector<PartonPtr> chainCaps) const {
  static DebugItem printstepdetails("DIPSY::PrintStepDetails", 6);
    if ( dr.getDipoles().size() == dr.initialDipoles().size() )
      return min(1.0, estimateF(&dr));
    vector<DipoleStatePtr> realI = extractI(chainCaps);
    vector<DipoleStatePtr> realH = extractH(chainCaps);
    double ret = 1.0;
    for ( int i = 0; i < int(realI.size()); i++ ) {
      double Fi = estimateF(realI[i]);
      double Gi = estimateF(realH[i]);
      //cout << ", (fi,Gi): (" << Fi << ", " << Gi << ")";
      ret *= Fi-Gi;
    }
    //cout << endl;
    return min(1.0, ret);
  }

  double DiffractiveEventFiller::estimateF(DipoleStatePtr I) const {
    InvEnergy scale = 3.0/GeV;
    list<DipolePtr> dips = I->getDipoles();
    double ret = 0.0;
    for ( list<DipolePtr>::iterator it = dips.begin(); it != dips.end(); it++ ) {
      DipolePtr dip = *it;
      ret += min(1.0, sqr(dip->size()/scale));
    }
    return ret;
  }

void DiffractiveEventFiller::setOffShell(DipoleState & dr) const {
  list<PartonPtr> partons = dr.getPartons();
  for ( list<PartonPtr>::iterator it = partons.begin(); it != partons.end(); it++) {
    tPartonPtr p = *it;
    p->onShell(false);
  }
}

void DiffractiveEventFiller::setValenceOnShell(DipoleState & dr) const {
  vector<DipolePtr> dips = dr.initialDipoles();
  for ( int i = 0; i < int(dips.size()); i++) {
    tPartonPtr p1 = dips[i]->partons().first;
    tPartonPtr p2 = dips[i]->partons().second;
    p1->onShell(true);
    p2->onShell(true);
  }
}

void DiffractiveEventFiller::setPartonsOnShell(DipolePtr dip) const {
  dip->partons().first->onShell(true);
  dip->partons().second->onShell(true);
}

double DiffractiveEventFiller::pickSubcascade(DipolePtr dip, int choice) const {
  static DebugItem printstepdetails("DIPSY::PrintStepDetails", 6);
  //this dipole is always on shell.
  setPartonsOnShell(dip);
  int configuration = choice;

  //find out how many subcascades come from first child
  int firstN = dip->children().first? nSubcascades(dip->children().first):0;
  //  int secondN = dip->children().second? nSubcascades(dip->children().second):0;

  //how much the no emissions case is boosted. INTERFACE!!
  //1 means default no emission, 2 means double probability.
  //double noEmBo = 1.0 + double(firstN*secondN)/3.0;
  //double noEmBo = 1.0;

  //as the no emission weights are changed, but the average weight has
  //to be maintained, a renormalisation is necessary
  //double weightNorm = (noEmBo + double(firstN*secondN))/
  //(1.0+double(firstN*secondN));
  // cout << "weightnorm: " << weightNorm << " options: "
  //      << firstN*secondN << endl;

  if ( printstepdetails )
    cout << "entering picksubcascade for a dipole with configuration number " << choice << endl;

  //configuration 0 is the one with no more emissions.
  if ( configuration == 0 ) {
    if ( printstepdetails ) cout << "this dipole has no further emissions" << endl;
    //reduce weight, as this option is more probable, see below
    // cout << "no emission 0, returning weight " << weightNorm/noEmBo << endl;
    //return weightNorm/noEmBo;
    return 1.0;
  }

  if ( !dip->children().first || !dip->children().second ) {
    Throw<DiffractiveCombinatorics>()
      << "Found a non-zero configuration number, but no children!! Will pretend configuration number is zero... :(" << Exception::warning;
    //return weightNorm/noEmBo;
    return 1.0;
  }
  
  //sometimes chose no emission anyways.
  //compensate this by reducing weight a factor 2 if that happens.
  //if ( UseRandom::rnd() < noEmBo/(noEmBo + double(firstN*secondN)) ) {
    // cout << "no emission 1, returning weight " << weightNorm/noEmBo << endl;
    //return weightNorm/noEmBo;
    //}
  // cout << "yes emission, emission options " << firstN*secondN << endl;
  //since no emission is excluded, remove 1 so that configuration is
  //between 0 and firstN*secondN
  configuration -= 1;
  if ( printstepdetails ) cout << "this dipole has (reduced) config number " << configuration << endl;

  //split up the configuration number into two independent numbers
  //one between 0 and firstN and one between 0 and secondN
  int firstConfiguration = configuration % firstN;
  int secondConfiguration = configuration/firstN;

  if ( printstepdetails ) cout << "split up into firstconf " << firstConfiguration << " (out of " << firstN << ")"
       << ", and second conf "<< secondConfiguration << " (out of "
       << nSubcascades(dip->children().second) << ")" << endl;

  //recur to each child with their respective configuration numbers
  //double weight = weightNorm;
  double weight = 1.0;
  weight *= pickSubcascade(dip->children().first, firstConfiguration);
  weight *= pickSubcascade(dip->children().second, secondConfiguration);
  // cout << "returning emitted weight " << weight << endl;
  //return weight;
  return 1.0;
}

void DiffractiveEventFiller::updateParentStructure(DipoleStatePtr dr) const {
  list<PartonPtr> partons = dr->getPartons();
  for ( list<PartonPtr>::iterator it = partons.begin(); it != partons.end(); it++) {
    tPartonPtr p = *it;
    //dont bother with the offshell ones
    if ( !p->onShell() ) continue;

    //if the parent is off shell, take the grandparent as new parent
    while ( p->parents().first && !p->parents().first->onShell() )
      p->parents(make_pair(p->parents().first->parents().first,
			   p->parents().second));
    while ( p->parents().second && !p->parents().second->onShell() )
      p->parents(make_pair(p->parents().first,
			   p->parents().second->parents().second));
  }
}



//checks through all partons, and returns the ones without children.
vector<PartonPtr> DiffractiveEventFiller::childlessChildren(DipoleStatePtr dr) const {
  static DebugItem printstepdetails("DIPSY::PrintStepDetails", 6);
  list<PartonPtr> partons = dr->getPartons();
  list<PartonPtr> childless = partons;
  //start with in the list, and remove the parents.
  //this is necessary as the partons do not remember their kids.
  while ( !partons.empty() ) {
    list<PartonPtr> offShells;
    for ( list<PartonPtr>::iterator it = partons.begin(); it != partons.end(); it++) {
      tPartonPtr p = *it;
      //remove parents from childless
      if ( p->parents().first ) {
	childless.remove(p->parents().first);
      }
      if ( p->parents().second ) {
	childless.remove(p->parents().second);
      }
      //save the off shell parents, as their parents in turn should not be counted as childless
      if ( p->parents().first && !p->parents().first->onShell() ) {
	offShells.push_back(p->parents().first);
	if ( printstepdetails )
	  cout << "found off shell at oy " << p->parents().first->oY() << endl;
      }
      if ( p->parents().second && !p->parents().second->onShell() ) {
	offShells.push_back(p->parents().second);
	if ( printstepdetails ) cout << "found off shell at oy " << p->parents().second->oY() << endl;
      }
    }
    partons = offShells;
    if ( printstepdetails ) cout << "found " << partons.size() << " off shell parents." << endl;
  }
  vector<PartonPtr> ret;
  for ( list<PartonPtr>::iterator it = childless.begin(); it != childless.end(); it++) {
    //first check that we don't add things without parents (ie valence)
    if ( !((*it)->parents().first) && !((*it)->parents().second) ) continue;
    ret.push_back(*it);
  }
  return ret;
}

//removes all valence partons in the vector
void DiffractiveEventFiller::removeValence(vector<PartonPtr> partons) const {
  for ( vector<PartonPtr>::iterator it = partons.begin(); it != partons.end(); it++ ) {
    if ( (*it)->valence() ) {
      it = partons.erase(it);
      it--;
    }
  }
}

void DiffractiveEventFiller::makeReal(DipoleState & dr, vector<PartonPtr> chainCaps,
				      DipoleState & de) const {
  static DebugItem printstepdetails("DIPSY::PrintStepDetails", 6);
  //prepare by setting all partons off shell, they will then be set on shell
  //if they make it through the machinery
  setOffShell(dr);

  //create the real state
  RealPartonStatePtr rs = new_ptr(RealPartonState());
  rs->addValence(dr);

  //run dr through a real state reweighting, with no interaction recoils
  //set interaction scale to that of the "interacting" dipoles
  for ( int i = 0; i < int(chainCaps.size()); i++ ) {
    if ( printstepdetails ) cout << "adding chainCap no " << i << endl;
    PartonPtr p = chainCaps[i];
    //if already on-shell, nothing has to be done
    if ( p->onShell() ) continue;
    DipolePtr d1 = p->dipoles().first;
    DipolePtr d2 = p->dipoles().second;
    //first find the scale of the cap no i from the smallest connected dipole
    //this will set the pt scale the chain starts the non-DGLAP suppression from.
    Energy iScale = ZERO;
    if ( d1 && d2 )
      iScale = p->pTScale()/min(d1->size(), d2->size());
    else if ( d1 && !d2 )
      iScale = p->pTScale()/d1->size();
    else if ( !d1 && d2 )
      iScale = p->pTScale()/d2->size();
    else   Throw<DiffractiveCombinatorics>()
      << "Found (childless) parton without dipoles in DiffractiveEventFiller::makeReal"
      << Exception::warning;
      
    //then run the machinery for reweighting and ordering
    if ( d1 )
      rs->fullControlEvolution(d1, DipolePtr(), true, true, iScale, iScale);
    if ( d2 )
      rs->fullControlEvolution(d2, DipolePtr(), true, true, iScale, iScale);

    //set the partons on shell
    if ( d1 )
      rs->setOnShell(d1);
    if ( d2 )
      rs->setOnShell(d2);
  }

  for ( RealParton::RealPartonSet::iterator it = rs->valence.begin();
	it != rs->valence.end(); it++ ) {
    (*it)->setOnShell();
  }

  //balance p- (and pt if overall recoil done) using dr
  //this then has to be inserted into the returned elastic particle
  RealPartonStatePtr ers = new_ptr(RealPartonState());
  ers->addValence(de);
  fixBoost(rs, ers);

  if ( printstepdetails ) cout << "balanced p+-" << endl;

  //remove the virtuals in the dipole states
  rs->mergeVirtuals();
  rs->saveState();
  removeVirtuals(&dr);
  ers->saveState();
}





double DiffractiveEventFiller::amplitude(DipoleState & dr, vector<PartonPtr> chainCaps,
					 const ImpactParameters & b) const {
  static DebugItem printstepdetails("DIPSY::PrintStepDetails", 6);
  if ( printstepdetails ) cout << setprecision(10);

  //make sure dr is in a state to be evolved and interact.
  //the state gets a bit messed up by reabsorbtion
  list<DipolePtr> dips = dr.getDipoles();
  for ( list<DipolePtr>::const_iterator it = dips.begin(); it != dips.end(); it++ ) {
    DipolePtr d = *it;
    d->reset();
    if ( d->neighbors().first && d->neighbors().first->colour() == d->colour() )
      dr.generateColourIndex(d);
    if ( d->neighbors().second && d->neighbors().second->colour() == d->colour() )
      dr.generateColourIndex(d);
  }

  DipoleEventHandlerPtr eh = Current<DipoleEventHandler>().ptr();

  // cout << "entering amplitude calculation for real state with " << chainCaps.size()
  //      << " interactions" << endl;

  //first check if the state is the valence one, if so return the elastic amplitude. otherwise do the diffractive excitation amplitude
  if ( dr.getDipoles().size() == dr.initialDipoles().size() )
    return elasticAmplitude(dr, b);

  //create one non-cap dipolestate (C) (these are the parents in the paper)
  DipoleStatePtr realC = steriliseCaps(&dr);

  //and for each cap:
  //one state of the two last dipoles (I_i) and one state of the second last dipole (X_i).
  vector<DipoleStatePtr> realI = extractI(chainCaps);
  vector<DipoleStatePtr> realH = extractH(chainCaps);
  //and the rapidities at which the caps were emitted.
  vector<double> intY = extractIntY(chainCaps);

  //loop over virtual cascades from these real subcascades, pair each with a stored
  //virtual cascades from the elastic side, and calculate the averages of
  //cascadeNonInt (e^{F_U})
  //hiddenInt: 1 - exp(-F_X_i) (for each i)
  //interactingInt: 1 - exp(-F_W_i) (for each i)
  double cascadeNonInt = 0.0;
  vector<double> hiddenInt(realI.size(), 0.0);
  vector<double> interactingInt(realI.size(), 0.0);
  vector<double> hiddenInthalf(realI.size(), 0.0);
  vector<double> interactingInthalf(realI.size(), 0.0);

  // cout << "vectors set up, starting loop" << endl;

  for (int i = 0; i < int(virtualCascades.size()); i++  ) {
    DipoleStatePtr virtU = realC->clone();
    virtU->makeOriginal();
    virtU->evolve(virtU->lowestY(), maxY());
    double FU = eh->xSecFn().sumf(*virtU, *virtualCascades[i], b);

    cascadeNonInt += exp(-FU);

    for ( int j = 0; j < int(realI.size()); j++ ) {
      DipoleStatePtr virtW = realI[j]->clone();
      virtW->makeOriginal();
      virtW->evolve(intY[j], maxY());
      double FWi = eh->xSecFn().sumf(*virtW, *virtualCascades[i], b);

      DipoleStatePtr virtX = realH[j]->clone();
      virtX->makeOriginal();
      virtX->evolve(intY[j], maxY());
      double FXi = eh->xSecFn().sumf(*virtX, *virtualCascades[i], b);

      hiddenInt[j] += eh->xSecFn().unitarize(FXi);
      interactingInt[j] += eh->xSecFn().unitarize(FWi);
      if ( double(i) < double(virtualCascades.size())/2.0 ) {
      hiddenInthalf[j] += eh->xSecFn().unitarize(FXi);
      interactingInthalf[j] += eh->xSecFn().unitarize(FWi);
      }
    }
  }

  //normalise
  cascadeNonInt /= virtualCascades.size();
  for ( int j = 0; j < int(realI.size()); j++ ) {
    hiddenInt[j] /= virtualCascades.size();
    interactingInt[j] /= virtualCascades.size();
    // cout << setprecision(10);
    // cout << "hidden: " << hiddenInt[j] << ", last: " << interactingInt[j]
    // 	 << ", diff: " << interactingInt[j] - hiddenInt[j] << ", relative error: "
    // 	 << (2.0*(interactingInthalf[j] - hiddenInthalf[j])/virtualCascades.size() - (interactingInt[j] - hiddenInt[j]))/(interactingInt[j] - hiddenInt[j])
    // 	 << ", parents: " << cascadeNonInt << endl;
  }

  //the final amplitude then is CascadeNonInt*hiddenNonInt*prod_i(interactingInti)
  double amp = cascadeNonInt;
  // cout << "cascade sudakov = " << cascadeNonInt;
  for ( int j = 0; j < int(realI.size()); j++ ) {
    // cout << ", interacting amplitude = " << interactingInt[j]
    //  	 << ", hidden amplitude = " << hiddenInt[j] << endl;
    amp *= interactingInt[j] - hiddenInt[j];
  }

  // cout << "returning amp " << amp << endl;

  return amp;
}

double DiffractiveEventFiller::elasticAmplitude(DipoleState & dr, const ImpactParameters & b) const {
  static DebugItem printstepdetails("DIPSY::PrintStepDetails", 6);
  if ( printstepdetails ) cout << setprecision(10);
  cout << "elastic state. calculate elastic amplitude. dips: "
       << dr.getDipoles().size() << ", partons: " << dr.getPartons().size() << endl;

  DipoleEventHandlerPtr eh = Current<DipoleEventHandler>().ptr();

  // //testing
  // DipoleStatePtr newState = eh->WFR().generate(*eh, 200.0*GeV);
  // cout << "A completely new state:" << endl;
  // newState->coutData();

  double amp = 0.0;
  double amphalf = 0.0;

  for (int i = 0; i < int(virtualCascades.size()); i++  ) {
    DipoleStatePtr virt = dr.clone();
    virt->makeOriginal();
    //DipoleStatePtr virt = copyActiveState(&dr);
    // DipoleStatePtr virtNew = eh->WFR().generate(*eh, 200.0*GeV);

    // set<DipolePtr> dipset = dr.allDipoles;
    //   for ( set<DipolePtr>::const_iterator it = allDipoles.begin();
    // 	    it != allDipoles.end(); it++ ) {

    //   }

    list<DipolePtr> newdips = virt->getDipoles();
    // cout << "copied unevolved dips:" << endl;
    // for ( list<DipolePtr>::const_iterator it = newdips.begin(); it != newdips.end(); it++ ) {
    //   DipolePtr d = *it;
    //   cout << "address: " << d
    // 	   << ", y1: " << d->partons().first->y()
    // 	   << ", y2: " << d->partons().first->y()
    // 	   << ", effective Parton 1: " << d->effectivePartons().first
    // 	   << ", effective Parton 2: " << d->effectivePartons().second
    // 	   << ", colour: " << d->colour() << endl;
    // }
    // newdips = virtNew->getDipoles();
    // cout << "new dips:" << endl;
    // for ( list<DipolePtr>::const_iterator it = newdips.begin(); it != newdips.end(); it++ ) {
    //   DipolePtr d = *it;
    //   cout << "address: " << d
    // 	   << ", y1: " << d->partons().first->y()
    // 	   << ", child1: " << d->children().first
    // 	   << ", child2: " << d->children().second
    // 	   << ", generated gluon: " << d->generatedGluon()
    // 	   << ", generatedy: " << d->generatedY()
    // 	   << ", swingdip: " << d->swingDipole()
    // 	   << ", colour: " << d->colour() << endl;
    // }

    if ( printstepdetails )  cout << "evolve from " << virt->lowestY() << " to " << maxY() << endl;
    virt->evolve(virt->lowestY(), maxY());
    // cout << "copied state has " << virt->getPartons().size() << " partons." << endl;
    // virtNew->evolve(virtNew->lowestY(), maxY());
    // cout << "new state has " << virtNew->getPartons().size() << " partons." << endl;

    double F = eh->xSecFn().sumf(*virt, *virtualCascades[i], b);
    F = eh->xSecFn().unitarize(F);

    if ( printstepdetails )
      cout << "for elastic cascade " << i << ", elastic amp is = " << F << endl;

    amp += F;
    if ( double(i) < double(virtualCascades.size())/2.0 ) {
      amphalf += F;  
    }
  }
  amp /= double(virtualCascades.size());
  double relError = abs(2.0*amphalf/virtualCascades.size() - amp)/amp;
  double correctionFactor = 1.0/sqrt(1.0+sqr(relError));


  //cout << "elastic amp: " << amp << ", relative error: "
  //      << relError << ", correction factor: " << correctionFactor << endl;

  return amp*correctionFactor;
}

// DipoleStatePtr DiffractiveEventFiller::copyActiveState(DipoleStatePtr dr) const {
//   DipoleStatePtr copy = DipoleState(*dr);
//   list<DipolePtr> dips = dr->getDipoles();
//   list<PartonPtr> ps = dr->getPartons();
//   vector<PartonPtr> newPs;
//   for ( list<DipolePtr>::const_iterator it = dips.begin();
// 	it != dips.end(); it++ ) {
    
//   }
// }

//removes the chain caps from the state dr.
DipoleStatePtr DiffractiveEventFiller::steriliseCaps(DipoleStatePtr dr) const {
  DipoleStatePtr state = dr->clone();
  vector<PartonPtr> caps = childlessChildren(state);
  removeValence(caps);

  list<DipolePtr> parents = state->getDipoles();
  // cout << "entering sterilisecaps, dipoles before: " << parents.size()
  //      << ", caps: " << caps.size() << endl;
  for ( int i = 0; i < int(caps.size()); i++ ) {
    parents.remove(caps[i]->dipoles().first);
    parents.remove(caps[i]->dipoles().second);
  }
  // cout << "dipoles after: " << parents.size() << endl;
  if ( parents.size() == 0 ) {
    DipoleStatePtr ret = new_ptr(DipoleState());
    ret->handler(Current<DipoleEventHandler>().ptr());
    return ret;
  }

  DipolePtr parent = *parents.begin();
  // cout << "parent dipole is " << parent << endl;
  for ( int i = 0; i < int(caps.size()); i++ ) {
    DipolePtr d1 = caps[i]->dipoles().first;
    DipolePtr d2 = caps[i]->dipoles().second;
    d1->children().second = d2;
    d2->children().second = d1;
    PartonPtr p1 = d1->partons().first;
    PartonPtr p2 = d2->partons().second;
    // cout << "cap: " << caps[i]->oY();
    // cout << ", p1: " << p1->oY();
    // cout << ", p2: " << p2->oY() << endl;
    p1->removeChild(caps[i]);
    p2->removeChild(caps[i]);
    d1->turnOff();
    d2->turnOff();
  }
  return state;
}

vector<DipoleStatePtr> DiffractiveEventFiller::extractI(vector<PartonPtr> chainCaps) const {
  static DebugItem printstepdetails("DIPSY::PrintStepDetails", 6);
  if ( printstepdetails ) cout << "entering extract I" << endl;
  vector<DipoleStatePtr> ret;
  for ( int i = 0; i < int(chainCaps.size()); i++ ) {
    PartonPtr c = chainCaps[i];
    if ( printstepdetails ) cout << "set up cap" << endl;
    if ( !c->dipoles().first || !c->dipoles().second )
      cerr << "quark chaincap in DEF::extractI!!!! :o" << endl; //make proper error message
    else {
      PartonPtr c1 = c->dipoles().first->partons().first;
      PartonPtr c2 = c->dipoles().second->partons().second;

      if ( printstepdetails ) cout << "create state, partons and dipoles" << endl;
      DipoleStatePtr state = new_ptr(DipoleState());
      state->handler(Current<DipoleEventHandler>().ptr());
      PartonPtr p = new_ptr(Parton());
      PartonPtr p1 = new_ptr(Parton());
      PartonPtr p2 = new_ptr(Parton());
      DipolePtr d1 = state->createDipole();
      DipolePtr d2 = state->createDipole();
      d1->partons(make_pair(p1, p));
      d2->partons(make_pair(p, p2));
      d1->neighbors(make_pair(DipolePtr(), d2));
      d2->neighbors(make_pair(d1, DipolePtr()));
      d1->colour(c->dipoles().first->colour());
      d2->colour(c->dipoles().second->colour());

      p->plus(c->plus());
      p->pT(c->pT());
      p->updateYMinus();
      p->oY(p->y());
      p->position(c->position());
      p->dipoles(make_pair(d1, d2));

      p1->plus(c1->plus());
      p1->pT(c1->pT());
      p1->updateYMinus();
      p1->oY(p1->y());
      p1->position(c1->position());
      p1->dipoles(make_pair(DipolePtr(), d1));

      p2->plus(c2->plus());
      p2->pT(c2->pT());
      p2->updateYMinus();
      p2->oY(p2->y());
      p2->position(c2->position());
      p2->dipoles(make_pair(d2, DipolePtr()));

      state->addDipole(*d1);
      state->addDipole(*d2);

      // cout << "L dip sizes: " << d1->size()*GeV << ", "
      // 	   << d2->size()*GeV;

      if ( printstepdetails ) cout << "extracted I" << endl;
      // state->plotState(true);

      ret.push_back(state);
    }
  }
  return ret;
}

vector<DipoleStatePtr> DiffractiveEventFiller::extractH(vector<PartonPtr> chainCaps) const {
  static DebugItem printstepdetails("DIPSY::PrintStepDetails", 6);
  vector<DipoleStatePtr> ret;
  for ( int i = 0; i < int(chainCaps.size()); i++ ) {
    PartonPtr c = chainCaps[i];
    if ( !c->dipoles().first || !c->dipoles().second )
      cerr << "quark chaincap in DEF::extractH!!!! :o" << endl; //make proper error message
    else {
    PartonPtr c1 = c->dipoles().first->partons().first;
    PartonPtr c2 = c->dipoles().second->partons().second;

    DipoleStatePtr state = new_ptr(DipoleState());
    state->handler(Current<DipoleEventHandler>().ptr());
    PartonPtr p1 = new_ptr(Parton());
    PartonPtr p2 = new_ptr(Parton());
    DipolePtr d = state->createDipole();
    d->partons(make_pair(p1, p2));
    d->neighbors(make_pair(DipolePtr(), DipolePtr()));

    InvEnergy r1 = (c->position() - c1->position()).pt();
    InvEnergy r2 = (c->position() - c2->position()).pt();
    double P1 = sqr(r2)/(sqr(r1) + sqr(r2));
    double P2 = 1.0 - P1;

    if ( P1 > 0.5 )
      d->colour(c->dipoles().first->colour());
    else
      d->colour(c->dipoles().second->colour());

    p1->plus(c1->plus() + P1*c->plus());
    p1->pT(c1->pT() + P1*c->pT());
    p1->updateYMinus();
    p1->oY(p1->y());
    p1->position(c1->position());
    p1->dipoles(make_pair(DipolePtr(), d));

    p2->plus(c2->plus() + P2*c->plus());
    p2->pT(c2->pT() + P2*c->pT());
    p2->updateYMinus();
    p2->oY(p2->y());
    p2->position(c2->position());
    p2->dipoles(make_pair(d, DipolePtr()));

    state->addDipole(*d);

    // cout << ", H dip sizes: " << d->size()*GeV << endl;

    if ( printstepdetails ) cout << "extracted H" << endl;
    // state->plotState(true);

    ret.push_back(state);
    }
  }
  return ret;
}

vector<double> DiffractiveEventFiller::extractIntY(vector<PartonPtr> chainCaps) const {
  vector<double> ret;
  for ( int i = 0; i < int(chainCaps.size()); i++ ) {
    ret.push_back(chainCaps[i]->oY());
    // cout << "intY is " << chainCaps[i]->oY() << ", ";
  }
  return ret;
}

DipoleState::String DiffractiveEventFiller::makeProton(DipoleStatePtr state) const {
  DipoleState::String ret;
  PartonPtr proton = new_ptr(Parton());
  proton->plus(state->plus());
  proton->minus(state->minus());
  proton->pT(TransverseMomentum());
  ret.push_back(proton);
  proton->flavour(ParticleID::pplus);
  return ret;
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void DiffractiveEventFiller::doinit() throw(InitException) {
  EventFiller::doinit();
  totalXSec = ZERO;
  totalSqr = ZERO;
  neve = 0;
}

void DiffractiveEventFiller::dofinish() {
  EventFiller::dofinish();
  totalXSec /= double(neve);
  totalSqr /= double(neve);
  CrossSection err = sqrt(totalSqr - sqr(totalXSec))/double(neve);
  generator()->log()
    << "inclusive elastic plus SD (up to this y) cross section is "
    << ouniterr(totalXSec, err, nanobarn) << " nb." << endl;

  generator()->log()
    << "total number of subcascades encountered: " << nSubCascades() << endl;
}

void DiffractiveEventFiller::doinitrun() {
  EventFiller::doinitrun();
  static DebugItem printstepdetails("DIPSY::PrintStepDetails", 6);
  if ( printstepdetails ) cout << "entering DiffractiveEventHandler doinitrun" << endl;
  DipoleEventHandlerPtr eh = Current<DipoleEventHandler>().ptr();
  DipoleStatePtr dummyState = eh->WFR().generate(*eh, 1000.0*GeV);
  Energy W = dummyState->handler().lumiFn().maximumCMEnergy();
  Energy2 a = eh->WFR().m2() - eh->WFL().m2() + sqr(W);
  Energy PR = (a + sqrt(sqr(a) - 4.0*eh->WFR().m2()*sqr(W)))*0.5/W;
  double elasticY = log(PR/sqrt(eh->WFR().m2()));
  maxY(elasticY - minRapGap());
  minY(elasticY - maxRapGap());
  initialiseVirtualCascades(*eh, eh->WFR(), PR, -maxY());
}

void DiffractiveEventFiller::persistentOutput(PersistentOStream & os) const {
  os << theNElasticCascades << theMinRapGap << theMaxRapGap << theMinY << theMaxY << theNInteractions;
}

void DiffractiveEventFiller::persistentInput(PersistentIStream & is, int) {
  is >> theNElasticCascades >> theMinRapGap >> theMaxRapGap >> theMinY >> theMaxY >> theNInteractions;
}

DescribeClass<DiffractiveEventFiller,EventFiller>
describeDIPSYDiffractiveEventFiller("DIPSY::DiffractiveEventFiller", "libAriadne5.so libDIPSY.so");

void DiffractiveEventFiller::Init() {

  static ClassDocumentation<DiffractiveEventFiller> documentation
    ("The DiffractiveEventFiller class is able to produce an initial ThePEG::Step "
     "from two colliding DipoleStates.");

  static Parameter<DiffractiveEventFiller,int> interfaceNElasticCascades
    ("NElasticCascades",
     "How many pregenerated elastic cascades should be used in the average over virtual cascades.",
     &DiffractiveEventFiller::theNElasticCascades, 0, 1, 0,
     true, false, Interface::lowerlim);

  static Parameter<DiffractiveEventFiller,double> interfaceMinRapGap
    ("MinRapGap",
     "The minimum rapidity gap in the generated (gluonic) events. There may be some leak due to "
     "FSR and hadronisation",
     &DiffractiveEventFiller::theMinRapGap, 3.0, 0.0, 0,
     true, false, Interface::lowerlim);

  static Parameter<DiffractiveEventFiller,double> interfaceMaxRapGap
    ("MaxRapGap",
     "The maximum rapidity gap in the generated (gluonic) events. There may be some leak due to "
     "FSR and hadronisation",
     &DiffractiveEventFiller::theMaxRapGap, 5.0, 0.0, 0,
     true, false, Interface::lowerlim);

  static Parameter<DiffractiveEventFiller, int> interfaceNInteractions
    ("NInteractions",
     "The number of interactions in the events being sampled."
     "0 is elastic"
     "-1 is the old algorithm with reweighting (seems to converge a bit faster)"
     "-2 is the old algorithm without reweighting (slower for the tested runs, but "
     "may be better for other runs.",
     &DiffractiveEventFiller::theNInteractions, -2, 0, 0,
     true, false, Interface::lowerlim);

}







//ymax should be the M_x limit, and the real state should be evolved to meet the virtuals
void  DiffractiveEventFiller::initialiseVirtualCascades
(DipoleEventHandler & eh, WaveFunction & WFR, Energy E, double ymax) {
  static DebugItem printstepdetails("DIPSY::PrintStepDetails", 6);
  if ( printstepdetails ) cout << "entering initialise virtual cascades" << endl;
  virtualCascades.clear();

  int virtualSamples = nElasticCascades();
  generator()->log() << "evolving elastic from "
		     << -eh.WFR().generate(eh, E)->lowestY()
		     << " to " << -ymax << endl;
  for ( int i = 0; i < virtualSamples; i ++ ) {
    DipoleStatePtr virt = eh.WFR().generate(eh, E);
    virt->evolve(virt->lowestY(), ymax);
    virtualCascades.push_back(virt);
    // if ( printstepdetails ) cout << "done virtual state no " << i << endl;
    // virt->plotState(true);
  }

}
