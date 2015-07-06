// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DipoleEventHandler class.
//

#include "DipoleEventHandler.h"
#include "SimpleProtonState.h"
#include "PhotonDipoleState.h"
#include "Parton.h"
#include "OldStyleEmitter.h"
#include "ThePEG/Interface/Switch.h"
#include "PT1DEmitter.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Handlers/LuminosityFunction.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Config/algorithm.h"
#include "gsl/gsl_sf_bessel.h"
#include "ThePEG/Utilities/Current.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/PDF/PartonExtractor.h"
#include "ThePEG/Utilities/Debug.h"
#include "ThePEG/Utilities/DebugItem.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "CPUTimer.h"

#include <iostream>
#include <fstream>

using namespace DIPSY;

DipoleEventHandler::DipoleEventHandler()
  : EventHandler(false), theNColours(9), theRMax(3.5*InvGeV),
    theBaryonSize(0.0*InvGeV),
    theCoherenceRange(0.5*InvGeV), theEffectivePartonMode(0), theCollisionType(0),
    doShowHistory(false), theLambdaQCD(0.22*GeV), theNF(0),
    theFixedAlphaS(0.0), theExternalAlphaS(ASPtr()),
    theWFR(WaveFunctionPtr()), theWFL(WaveFunctionPtr()),
    theBGen(ImpactParameterGeneratorPtr()), theYFrame(0.5), theFudgeME(false),
    thePreSamples(1000), thePreSampleL(1), thePreSampleR(1), thePreSampleB(1),
    theXSecFn(DipoleXSecPtr()), theEmitter(EmitterPtr()),
theSwinger(SwingerPtr()) {}

DipoleEventHandler::~DipoleEventHandler() {}

IBPtr DipoleEventHandler::clone() const {
  return new_ptr(*this);
}

IBPtr DipoleEventHandler::fullclone() const {
  return new_ptr(*this);
}

double DipoleEventHandler::elapsed() {
  static long last = 0;
  long el = clock() - last;
  last += el;
  return double(el)*0.000001;
}

void DipoleEventHandler::presample() {
  //redirect diffractive runs, as the cross section calculation
  //then requires more than a single loop over collisions.
  if ( collisionType() == 1 ) return diffractivePresample();

  //Set up before the loop
  Current<DipoleEventHandler> current(this);
  for_each(analyses, mem_fun(&DipoleAnalysisHandler::initialize));
  CrossSection xmax = 0.0*picobarn;
  CrossSection xnmax = 0.0*picobarn;
  elapsed();
  generator()->log()
    << endl << "Starting DIPSY run at CoM energy "
    << lumiFn().maximumCMEnergy()/GeV << " GeV" << endl;

  //Loop over collisions
  for ( int i = 0, N = preSamples(); i < N; ++i ) {
    //Setup before collision. Energy, jacobian, luminosity.
    Energy W = lumiFn().maximumCMEnergy();
    double jac = 1.0;
    pair<double,double> ll(1.0, 1.0);
    if ( int ndim = lumiFn().nDim(make_pair(WFL().particle(), WFR().particle())) ) {
      vector<double> r = UseRandom::rndvec(ndim);
      ll = lumiFn().generateLL(&(r[0]), jac);
      W *= exp(-0.5*(ll.first + ll.second));
    }
    
    //Find the p+ (or p-) of each state.
    Energy2 a = WFL().m2() - WFR().m2() + sqr(W);
    Energy PL = (a + sqrt(sqr(a) - WFL().m2()*sqr(W)))*0.5/W;
    a = WFR().m2() - WFL().m2() + sqr(W);
    Energy PR = (a + sqrt(sqr(a) - WFR().m2()*sqr(W)))*0.5/W;

    //Create the dipole states and impact parameters.
    vector<DipoleStatePtr> vl(preSampleL());
    vector<DipoleStatePtr> vr(preSampleR());
    vector<ImpactParameters> vb(preSampleB());
    vector< vector< vector<double> > >
      probs(preSampleL(),
	    vector< vector<double> >(preSampleR(), vector<double>(preSampleB(), 0.0)));

    // Setup a number of left- and right-moving states
    double highYL = 0.0;
    for ( int il = 0; il < preSampleL(); ++il ) {
      DipoleStatePtr dl = WFL().generate(*this, PL);
      vl[il] = dl;
      dl->collidingEnergy(PR);
      highYL += dl->highestY();
    }
    double highYR = 0.0;
    for ( int ir = 0; ir < preSampleR(); ++ir ) {
      DipoleStatePtr dr = WFR().generate(*this, PR);
      vr[ir] = dr;
      dr->collidingEnergy(PL);
      highYR += dr->highestY();
    }

    // Evolve the states
    double y0 = interactionFrame(highYL/double(preSampleL()), -highYR/double(preSampleR()));
    for ( int il = 0; il < preSampleL(); ++il ) {
      vl[il]->evolve(vl[il]->lowestY(), y0);
      vl[il]->unifyColourSystems();
    }
    for ( int ir = 0; ir < preSampleR(); ++ir ) {
      vr[ir]->evolve(vr[ir]->lowestY(), -y0);
      vr[ir]->unifyColourSystems();
    }

    for ( int ib = 0; ib < preSampleB(); ++ib ) {
      ImpactParameters b = bGen().generate();
      vb[ib] = b;
    }

    for ( int il = 0; il < preSampleL(); ++il )
      for ( int ir = 0; ir < preSampleR(); ++ir )
	for ( int ib = 0; ib < preSampleB(); ++ib ) {
	  DipoleStatePtr dl =vl[il];
	  DipoleStatePtr dr =vr[ir];
	  ImpactParameters b = vb[ib];
	  //Calculate interaction probability, and apply all the weights
	  double prob = xSecFn().sumf(*dr, *dl, b);
	  probs[il][ir][ib] = prob;
	  CrossSection weight = sqr(hbarc)*dr->weight()*dl->weight()*b.weight()*jac;
	  //This is the contribution to the non-diffractive
	  //cross section 1 - exp(-2(amp)^2)
	  CrossSection x = weight*xSecFn().unitarize(2.0*prob);
	  if ( x > xmax ) {
	    xnmax = xmax;
	    xmax = x;
	  } else
	    xnmax = max(x, xnmax);

	  //Call the analysis (which calculates cross sections)
	  for ( int i = 0, N = analyses.size(); i < N; ++i )
	    analyses[i]->analyze(*dr, *dl, b, xSecFn(), prob, weight);
	}
    // Call the analysis (which calculates cross sections) for all
    // systems and impact parameters.
    for ( int i = 0, N = analyses.size(); i < N; ++i )
      analyses[i]->analyze(vr, vl, vb, xSecFn(), probs, jac);    
  }

  //Write output (also from the called analyses) to log file.
  if ( preSamples() > 0 ) {
    generator()->log()
      << "Presampled " << preSamples() << " collisions ("
      << elapsed()/max(preSamples(),0) << " seconds per collisions)" << endl
      << "Maximum cross section: "
      << ouniterr(xmax, xmax-xnmax, nanobarn) << " nb." <<endl;
    for_each(analyses,
	     bind2nd(mem_fun(&DipoleAnalysisHandler::finalize),
		     preSamples()*preSampleL()*preSampleR()*preSampleB()));
    generator()->log() << endl;
  }
  stats = XSecStat(xmax > 0.0*nanobarn? xmax: 1.0*millibarn);
}

void DipoleEventHandler::diffractivePresample() {

  //Setup before 
  Current<DipoleEventHandler> current(this);
  for_each(analyses, mem_fun(&DipoleAnalysisHandler::initialize));
  CrossSection xmax = 0.0*picobarn;
  CrossSection xnmax = 0.0*picobarn;
  elapsed();
  generator()->log()
    << endl << "Starting Diffractive DIPSY run at CoM energy "
    << lumiFn().maximumCMEnergy()/GeV << " GeV" << endl;
  Energy W = lumiFn().maximumCMEnergy();
  double jac = 1.0;
  pair<double,double> ll(1.0, 1.0);

  //Create the starting states from their wavefunction.
  if ( int ndim = lumiFn().nDim(make_pair(WFL().particle(), WFR().particle())) ) {
    vector<double> r = UseRandom::rndvec(ndim);
    ll = lumiFn().generateLL(&(r[0]), jac);
    W *= exp(-0.5*(ll.first + ll.second));
  }

  //Calculate the energy of each incoming state in this frame.
  Energy2 a = WFL().m2() - WFR().m2() + sqr(W);
  Energy PL = (a + sqrt(sqr(a) - WFL().m2()*sqr(W)))*0.5/W;
  a = WFR().m2() - WFL().m2() + sqr(W);
  Energy PR = (a + sqrt(sqr(a) - WFR().m2()*sqr(W)))*0.5/W;

  //Set up rapidity intervals to evolve over.
  double y0 = interactionFrame(-log(PL/sqrt(abs(WFL().m2()))), log(PR/sqrt(abs(WFR().m2()))));

  //Pre-generate the virtual cascade, ie the ones from the
  //elastic side (right side by default).
  vector<DipoleStatePtr> virtualCascades;
  for ( int i = 0, N = preSamples(); i < N; i ++ ) {
    DipoleStatePtr virt = WFR().generate(*this, PR);
    virt->collidingEnergy(PL);
    virt->evolve(virt->lowestY(), -y0);
    virtualCascades.push_back(virt);
  }

  CrossSection sigmaSD = ZERO;
  QTY<4,0,0>::Type sigmaSD2 = ZERO;

  //Now run main loop over the real (excited in final state) cascades.
  for ( int i = 0, N = preSamples(); i < N; ++i ) {

    DipoleStatePtr dl = WFL().generate(*this, PL);
    dl->collidingEnergy(PR);
    dl->evolve(dl->lowestY(), y0);
    ImpactParameters b = bGen().generate();

    double sum = 0.0;
    double sumWeights = 0.0;

    //Loop over the pregenerated virtual  elastic
    for ( int j = 0; j < N; j ++ ) {
      double amp = xSecFn().unitarize(xSecFn().sumf(*virtualCascades[j], *dl, b));
      sum += amp*virtualCascades[j]->weight();
      sumWeights += virtualCascades[j]->weight();
    }
    double average = sum/sumWeights;

    //Take statistics on inclusive single diffractive cross section.
    CrossSection weight = sqr(hbarc)*dl->weight()*b.weight()*jac;
    CrossSection x = weight*average;
    sigmaSD += weight*sqr(average);
    sigmaSD2 += sqr(weight*sqr(average));
    if ( x > xmax ) {
      xnmax = xmax;
      xmax = x;
    } else
      xnmax = max(x, xnmax);
    //Pass on to analysers.
    for ( int j = 0, M = analyses.size(); j < M; ++j )
      analyses[j]->analyze(*WFR().generate(*this, PR), *dl, b, xSecFn(), average, weight);
  }

  //Output results to log file.
  if ( preSamples() > 0 ) {
    sigmaSD /= preSamples();
    sigmaSD2 /= preSamples();
    //Estimate error from the fluctuations.
    CrossSection err = sqrt((sigmaSD2 - sqr(sigmaSD))/preSamples());
    generator()->log()
      << "Presampled " << preSamples() << " real cascades, collided with "
      << preSamples() << " virtual cascades ("
      << elapsed()/max(preSamples(),0) << " seconds per real cascade)" << endl
      << "Maximum cross section: "
      << ouniterr(xmax, xmax-xnmax, nanobarn) << " nb." <<endl
      << "Integrated single diffractive cross section up to y = " << y0
      << " is " << ouniterr(sigmaSD, err, nanobarn) << " nb." <<endl;
    for_each(analyses,
	     bind2nd(mem_fun(&DipoleAnalysisHandler::finalize), preSamples()));
    generator()->log() << endl;
  }
  stats = XSecStat(xmax > 0.0*nanobarn? xmax: 1.0*millibarn);
}


void DipoleEventHandler::initialize() {
  Current<DipoleEventHandler> current(this);
  theWFL->initialize(*this);
  theWFR->initialize(*this);
}

EventPtr DipoleEventHandler::generateEvent() {
  //Redirect diffractive events.
  if ( collisionType() == 1 ) return generateDiffractiveEvent();

  while ( true ) {
    Current<DipoleEventHandler> current(this);

    //Set up energy, jacobians and wavefunctions.
    Energy W = lumiFn().maximumCMEnergy();
    double jac = 1.0;
    pair<double,double> ll(1.0, 1.0);
    if ( int ndim = lumiFn().nDim(make_pair(WFL().particle(), WFR().particle())) ) {
      vector<double> r = UseRandom::rndvec(ndim);
      ll = lumiFn().generateLL(&(r[0]), jac);
    }

    //Find lightcone momenta.
    Energy2 a = WFL().m2() - WFR().m2() + sqr(W);
    Energy PL = (a + sqrt(sqr(a) - 4.0*WFL().m2()*sqr(W)))*0.5/W;
    a = WFR().m2() - WFL().m2() + sqr(W);
    Energy PR = (a + sqrt(sqr(a) - 4.0*WFR().m2()*sqr(W)))*0.5/W;

    PPair inc(WFL().particle()->produceParticle(lightCone(PL, WFL().m2()/PL)),
	      WFR().particle()->produceParticle(lightCone(WFR().m2()/PR, PR)));
    LorentzRotation cmboost =
      lumiFn().getBoost()*LorentzRotation(0.0, 0.0, tanh(0.5*(ll.second - ll.first)));
    inc.first->transform(cmboost);
    inc.second->transform(cmboost);

    //Create dipole states.
    DipoleStatePtr dr = WFR().generate(*this, PR);
    DipoleStatePtr dl = WFL().generate(*this, PL);

    //in some settings, the states need to know about the
    //energy of the other state.
    dr->collidingEnergy(PL);
    dl->collidingEnergy(PR);

    double y0 = interactionFrame(dl->highestY(), -dr->highestY());

    dr->evolve(dr->lowestY(), -y0);
    dr->unifyColourSystems();
    dl->evolve(dl->lowestY(), y0);
    dl->unifyColourSystems();

    // This generates the impact parameter depending on the position
    //of the partons. Intention was to generate high-pt events more
    //often by selecting a b such that two partons would end up
    //very close to each other more often. Should work in theory, but
    //didn't make as big difference as hoped, maybe due to too many
    //of the high pTs coming from the cascade.
    // vector<pair<Parton::Point, InvEnergy> > points1 = dl->points();
    // vector<pair<Parton::Point, InvEnergy> > points2 = dr->points();
    // ImpactParameters b = bGen().generateDynamic(points1, points2);

    //Standard impact parameter generation.
    ImpactParameters b = bGen().generate();
    inc.first->setVertex
      (LorentzPoint(hbarc*b.bVec().x()/2, hbarc*b.bVec().y()/2, ZERO, ZERO));
    inc.second->setVertex
      (LorentzPoint(-hbarc*b.bVec().x()/2, -hbarc*b.bVec().y()/2, ZERO, ZERO));
 
    currentEvent(new_ptr(Event(inc, this, generator()->runName(),
			       generator()->currentEventNumber(), 1.0)));
    currentCollision(new_ptr(Collision(inc, currentEvent(), this)));
    if ( currentEvent() ) currentEvent()->addCollision(currentCollision());
    currentStep(new_ptr(Step(currentCollision())));
    currentCollision()->addStep(currentStep());

    double genweight = sqr(hbarc)*dr->weight()*dl->weight()*b.weight()*jac/
      stats.maxXSec();

    //Pass to filler. This is where the final state is decided.
    double prob = eventFiller().fill(*currentStep(), *this, inc, *dl, *dr, b);

    double weight = prob*genweight;
    stats.select(weight);
    currentEvent()->weight(weight);
    if ( weight <= 0.0 ) continue;
    stats.accept();

    try {
      initGroups();
      continueCollision();
      return currentEvent();
    }
    catch (Veto) {
      stats.reject(weight);
    }
    catch (Stop) {
      break;
    }
    catch (Exception &) {
      stats.reject(weight);
      throw;
    }
  }
  return currentEvent();
}

EventPtr DipoleEventHandler::generateDiffractiveEvent() {

  //Very similar to generateEvent.
  while ( true ) {
    //Same setup as in the non-diffractive case.
    Current<DipoleEventHandler> current(this);
    Energy W = lumiFn().maximumCMEnergy();
    double jac = 1.0;
    pair<double,double> ll(1.0, 1.0);
    if ( int ndim = lumiFn().nDim(make_pair(WFL().particle(), WFR().particle())) ) {
      vector<double> r = UseRandom::rndvec(ndim);
      ll = lumiFn().generateLL(&(r[0]), jac);
    }
    Energy2 a = WFL().m2() - WFR().m2() + sqr(W);
    Energy PL = (a + sqrt(sqr(a) - 4.0*WFL().m2()*sqr(W)))*0.5/W;
    a = WFR().m2() - WFL().m2() + sqr(W);
    Energy PR = (a + sqrt(sqr(a) - 4.0*WFR().m2()*sqr(W)))*0.5/W;
    PPair inc(WFL().particle()->produceParticle(lightCone(PL, WFL().m2()/PL)),
	      WFR().particle()->produceParticle(lightCone(WFL().m2()/PR, PR)));
    LorentzRotation cmboost =
      lumiFn().getBoost()*LorentzRotation(0.0, 0.0, tanh(0.5*(ll.second - ll.first)));
    inc.first->transform(cmboost);
    inc.second->transform(cmboost);

    //create real excited valence. Evolution is done in diffFill.
    DipoleStatePtr dr = WFL().generate(*this, PL);
    //create elastic
    DipoleStatePtr de = WFR().generate(*this, PR);

    ImpactParameters b = bGen().generate();

    inc.first->setVertex
      (LorentzPoint(hbarc*b.bVec().x()/2, hbarc*b.bVec().y()/2, ZERO, ZERO));
    inc.second->setVertex
      (LorentzPoint(-hbarc*b.bVec().x()/2, -hbarc*b.bVec().y()/2, ZERO, ZERO));
 
    currentEvent(new_ptr(Event(inc, this, generator()->runName(),
    			       generator()->currentEventNumber(), 1.0)));
    currentCollision(new_ptr(Collision(inc, currentEvent(), this)));
    if ( currentEvent() ) currentEvent()->addCollision(currentCollision());
    currentStep(new_ptr(Step(currentCollision())));
    currentCollision()->addStep(currentStep());

    double genweight = sqr(hbarc)*dr->weight()*de->weight()*b.weight()*jac/
      stats.maxXSec();

    //Call filler, which will be redirected to diffractive version.
    double prob = eventFiller().fill(*currentStep(), *this, inc, *dr, *de, b);
    double weight = prob*genweight;

    if ( weight <= 0.0 ) continue;
    if ( isnan(weight) ) {
      continue;
    }

    stats.select(weight);
    currentEvent()->weight(weight);
    stats.accept();

    try {
      initGroups();
      continueCollision();
      return currentEvent();
    }
    catch (Veto) {
      stats.reject(weight);
    }
    catch (Stop) {
      break;
    }
    catch (Exception &) {
      stats.reject(weight);
      throw;
    }
  }
  return currentEvent();
}

CrossSection DipoleEventHandler::integratedXSec() const {
  return stats.xSec();
}

CrossSection DipoleEventHandler::integratedXSecErr() const {
  return stats.xSecErr();
}

CrossSection DipoleEventHandler::maxXSec() const {
  return stats.maxXSec();
}

CrossSection DipoleEventHandler::histogramScale() const {
  return stats.xSec()/stats.sumWeights();
}

void DipoleEventHandler::statistics(ostream & os) const {
  string line = "======================================="
    "=======================================\n";

  if ( stats.accepted() <= 0 ) {
    os << line << "No events generated by event handler '" << name() << "'."
       << endl;
      return;
  }

  os << line << "Statistics for event handler \'" << name() << "\':\n"
     << "                                       "
     << "generated    number of    Cross-section\n"
     << "                                       "
     << "   events     attempts             (nb)\n";

  os << line << "Total:" << setw(42) << stats.accepted() << setw(13)
     << stats.attempts() << setw(17)
     << ouniterr(stats.xSec(), stats.xSecErr(), nanobarn)
     << endl << line;


}

void DipoleEventHandler::doinit() throw(InitException) {

  Current<DipoleEventHandler> current(this);
	
  EventHandler::doinit();
}

void DipoleEventHandler::dofinish() {
  generator()->log() << endl << "Ending DIPSY run at CoM energy "
       << lumiFn().maximumCMEnergy()/GeV << " GeV" << endl;
  theEventFiller->finish();
  EventHandler::dofinish();
  CPUClock::dump(generator()->log());
}

void DipoleEventHandler::doinitrun() {
  Current<DipoleEventHandler> current(this);

  EventHandler::doinitrun();
  theNErr = 0;
  for_each(analyses, mem_fun(&DipoleAnalysisHandler::initrun));
  if ( theSwinger ) theSwinger->initrun();
  theEmitter->initrun();
  theBGen->initrun();
  theXSecFn->initrun();
  theEventFiller->initrun();
  theWFL->initrun();
  theWFR->initrun();
  if ( externalAlphaS() ) externalAlphaS()->initrun();
  presample();

  // Fix up a dummy XComb object
  Energy emax = lumiFn().maximumCMEnergy();
  cuts()->initialize(sqr(emax), 0.0);
  cPDPair incoming = make_pair(WFL().particle(), WFR().particle());
  PartonPairVec bins = partonExtractor()->getPartons(emax, incoming, *cuts());
  theLastXComb =  new_ptr(XComb(emax, incoming, this, partonExtractor(),
				tCascHdlPtr(), bins[0], cuts()));
  theLastXComb->lastSHat(sqr(emax));
  theLastXComb->lastY(0.0);
  theLastXComb->lastP1P2(make_pair(0.0, 0.0));
  theLastXComb->lastL1L2(make_pair(0.0, 0.0));
  theLastXComb->lastX1X2(make_pair(1.0, 1.0));
  theLastXComb->lastScale(sqr(emax));
  theLastXComb->lastAlphaS(1.0);
  theLastXComb->lastAlphaEM(1.0/137.0);
  partonExtractor()->select(theLastXComb);
}

void DipoleEventHandler::rebind(const TranslationMap & trans) throw(RebindException) {
  // dummy = trans.translate(dummy);
  EventHandler::rebind(trans);
}

IVector DipoleEventHandler::getReferences() {
  IVector ret = EventHandler::getReferences();
  // ret.push_back(dummy);
  return ret;
}


double DipoleEventHandler::alphaS(InvEnergy r) const {
  if ( fixedAlphaS() > 0.0 )
    return fixedAlphaS();
  else if ( LambdaQCD() > 0.0*GeV && LambdaQCD() < 2.0/rMax() )
    return 6.0*Constants::pi/((33.0 - 2.0*nF())*
			      log(2.0/(min(r,rMax())*LambdaQCD())));
  else if ( externalAlphaS() )
    return externalAlphaS()->value(sqr(2.0/min(r,rMax())), SM());
  else
    return SM().alphaS(sqr(2.0/min(r,rMax())));
}

double DipoleEventHandler::interactionFrame(double ymin, double ymax) const {
  double yframe = yFrame();
  if ( yFrame() < -1.0 && yFrame() >= -1.5 ) {
    double dum = -(yFrame() + 1.0);
    yframe = UseRandom::rnd(dum, 1.0 - dum);
  }
  double y0 = ymin + yframe*(ymax - ymin);
  if ( y0 > ymax ) y0 = max(ymin, ymax + (1.0 - yframe));
  else if ( y0 < ymin ) y0 = min(ymax, ymin + yframe);
  return y0;
}

void DipoleEventHandler::persistentOutput(PersistentOStream & os) const {
  os << theNColours << ounit(theRMax, InvGeV) << ounit(theBaryonSize, InvGeV)
     << ounit(theCoherenceRange, InvGeV)
     << theEffectivePartonMode << theCollisionType << doShowHistory
     << ounit(theLambdaQCD, GeV) << theNF << theFixedAlphaS << theExternalAlphaS
     << theWFR << theWFL << theBGen << theYFrame << theFudgeME << thePreSamples << thePreSampleL
     << thePreSampleR << thePreSampleB << theXSecFn << theEmitter << theSwinger
     << theEventFiller << analyses << stats;
}

void DipoleEventHandler::persistentInput(PersistentIStream & is, int) {
  is >> theNColours >> iunit(theRMax, InvGeV) >> iunit(theBaryonSize, InvGeV)
     >> iunit(theCoherenceRange, InvGeV)
     >> theEffectivePartonMode >> theCollisionType >> doShowHistory
     >> iunit(theLambdaQCD, GeV) >> theNF >> theFixedAlphaS >> theExternalAlphaS
     >> theWFR >> theWFL >> theBGen >> theYFrame >> theFudgeME >> thePreSamples >> thePreSampleL
     >> thePreSampleR >> thePreSampleB >> theXSecFn >> theEmitter >> theSwinger
     >> theEventFiller >> analyses >> stats;
}


// Static variable needed for the type description system in ThePEG.
#include "ThePEG/Utilities/DescribeClass.h"
DescribeClass<DipoleEventHandler,EventHandler>
  describeDIPSYDipoleEventHandler("DIPSY::DipoleEventHandler",
				  "libAriadne5.so libDIPSY.so");


void DipoleEventHandler::Init() {

  static ClassDocumentation<DipoleEventHandler> documentation
    ("The DipoleEventHandler is a special EventHandler capable of generating "
     "minimum-bias events using the Lund version of the Mueller dipole model.");

  static Parameter<DipoleEventHandler,int> interfaceNColours
    ("NColours",
     "The number of different colour indices in the swing mechanism. Should be "
     "\\f$N_c^2\\f$ but can be varied to investigate swing effects.",
     &DipoleEventHandler::theNColours, 9, 3, 0,
     true, false, Interface::lowerlim);

  static Parameter<DipoleEventHandler,InvEnergy> interfaceRMax
    ("RMax",
     "The general hadronic size in units of inverse GeV.",
     &DipoleEventHandler::theRMax, InvGeV, 3.5*InvGeV, 0.0*InvGeV, 0*InvGeV,
     true, false, Interface::lowerlim);

  static Parameter<DipoleEventHandler,InvEnergy> interfaceCoherenceRange
    ("CoherenceRange",
     "The maximum range at which partons are allowed to emit coherently."
     "This is also used as maximum range in the DGLAP suppression",
     &DipoleEventHandler::theCoherenceRange, InvGeV, 0.5*InvGeV, 0.0*InvGeV, 0*InvGeV,
     true, false, Interface::lowerlim);

  static Switch<DipoleEventHandler,int> interfaceEffectivePartonMode
    ("EffectivePartonMode",
     "How the partons are grouped inteo effective partons.",
     &DipoleEventHandler::theEffectivePartonMode, 0, true, false);
  static SwitchOption interfaceEffectivePartonModeColours
    (interfaceEffectivePartonMode,
     "Colours",
     "Groups with colour neighbours. Makes more sence since the "
     "colour flow determines how the emissions are made, but the "
     "swing makes this mode not always group up recoiling partons, "
     "which can mess up availible phase space and ordering. "
     "This is default.",
     0);
  static SwitchOption interfaceEffectivePartonModeFastColours
    (interfaceEffectivePartonMode,
     "FastColours",
     "Groups with colour neighbours as for the \"Colour\" option, "
     "but speed up the generation by caching different ranges for "
     "the same parton.",
     2);
  static SwitchOption interfaceEffectivePartonModeFastColours2
    (interfaceEffectivePartonMode,
     "FastColours2",
     "Groups with colour neighbours as for the \"FastColour\" option, "
     "but speed up the generation by caching different ranges for "
     "the same parton.",
     3);
  static SwitchOption interfaceEffectivePartonModeRelatives
    (interfaceEffectivePartonMode,
     "Relatives",
     "Groups with parents and childs. Always pairs up with the "
     "recoiling pt, and should get the phase space and ordering "
     "correct, but makes the emissions depend on the history"
     ", not only current state.",
     1);

  static Switch<DipoleEventHandler,int> interfaceCollisionType
    ("CollisionType",
     "What type of collison.",
     &DipoleEventHandler::theCollisionType, 0, true, false);
  static SwitchOption interfaceCollisionTypeNonDiffractive
    (interfaceCollisionType,
     "NonDiffractive",
     "Should be self explanatory. :P",
     0);
  static SwitchOption interfaceTypeSingleDiffractive
    (interfaceCollisionType,
     "SingleDiffractive",
     "Should be self explanatory. :P",
     1);

  static Parameter<DipoleEventHandler,bool> interfaceShowHistory
    ("ShowHistory",
     "If the history of every event should be studied manually."
     "Note that only final state events get studied, not the presamples.",
     &DipoleEventHandler::doShowHistory, false, false, false,
     true, false, Interface::lowerlim);

  static Parameter<DipoleEventHandler,Energy> interfaceLambdaQCD
    ("LambdaQCD",
     "The value of \\f$\\Lambda_{QCD}\\f$ to be used in the running coupling. "
     "If zero, the <interface>AlphaS</interface> object will be used for the "
     "coupling instead.",
     &DipoleEventHandler::theLambdaQCD, GeV, 0.22*GeV, 0.0*GeV, 0*GeV,
     true, false, Interface::lowerlim);

  static Parameter<DipoleEventHandler,double> interfaceFixedAlphaS
    ("FixedAlphaS",
     "The value of the constant coupling. If zero, a running coupling is "
     "assumed.",
     &DipoleEventHandler::theFixedAlphaS, 0.0, 0.0, 0,
     true, false, Interface::lowerlim);

  static Reference<DipoleEventHandler,AlphaSBase> interfaceExternalAlphaS
    ("ExternalAlphaS",
     "An external \\f$\\alpha_S\\f$ object to be used if "
     "<interface>LambdaQCD</interface> is zero. If null the object "
     "specified in the overall StandardModel object will be used insted.",
     &DipoleEventHandler::theExternalAlphaS, true, false, true, true, false);

  static Parameter<DipoleEventHandler,int> interfaceNF
    ("NF",
     "The number of flavours to be used in the running coupling. Not active "
     "if an external running coupling (<interface>ExternalAlphaS</interface> "
     "is used.",
     &DipoleEventHandler::theNF, 3, 1, 0,
     true, false, Interface::lowerlim);

  static Reference<DipoleEventHandler,WaveFunction> interfaceWFR
    ("WFR",
     "The wave function of the incoming particle along the positive z-axis.",
     &DipoleEventHandler::theWFR, true, false, true, false, false);

  static Reference<DipoleEventHandler,WaveFunction> interfaceWFL
    ("WFL",
     "The wave function of the incoming particle along the negative z-axis.",
     &DipoleEventHandler::theWFL, true, false, true, false, false);

  static Parameter<DipoleEventHandler,double> interfaceYFrame
    ("YFrametest",
     "Indicate in which frame the dipole systems should collide. A value of "
     "0.5 means that both systems will be evolved an equal rapidity distance. "
     "0.0 (1.0) means that the right(left)-moving system will not be evolved "
     "at all while the left(right)-moving system will be evolved as much as "
     "possible. A value larger than 1.0 (less than 0.0) means the right-(left-)moving "
     "system will be evolved a fixed rapidity interval YFrametest - 1 (-YFrametest)",
     &DipoleEventHandler::theYFrame, 0.5, 0.0, 1.0,
     true, false, Interface::nolimits);


  static Switch<DipoleEventHandler,bool> interfaceFudgeME
    ("FudgeME",
     "Indicate whether a fudge factor should be included to tame the high-pt tail "
     "according to an approximate matrix element correction.",
     &DipoleEventHandler::theFudgeME, false, true, false);
  static SwitchOption interfaceFudgeMEFudge
    (interfaceFudgeME,
     "Fudge",
     "Include fudge factor",
     true);
  static SwitchOption interfaceFudgeMENoFudge
    (interfaceFudgeME,
     "NoFudge",
     "Do not include fudge factor",
     false);

  static Parameter<DipoleEventHandler,int> interfacePreSamples
    ("PreSamples",
     "The number of collisions to analyze in the presampling.",
     &DipoleEventHandler::thePreSamples, 1000, 0, 0,
     true, false, Interface::lowerlim);

  static Parameter<DipoleEventHandler,int> interfacePreSampleL
    ("PreSampleL",
     "The number of left-moving systems to generate for each presample.",
     &DipoleEventHandler::thePreSampleL, 1, 1, 0,
     true, false, Interface::lowerlim);

  static Parameter<DipoleEventHandler,int> interfacePreSampleR
    ("PreSampleR",
     "The number of right-moving systems to generate for each presample.",
     &DipoleEventHandler::thePreSampleR, 1, 1, 0,
     true, false, Interface::lowerlim);

  static Parameter<DipoleEventHandler,int> interfacePreSampleB
    ("PreSampleB",
     "The number ofimpact parameters to generate for each presample.",
     &DipoleEventHandler::thePreSampleB, 1, 1, 0,
     true, false, Interface::lowerlim);

  static Reference<DipoleEventHandler,DipoleXSec> interfaceXSecFn
    ("XSecFn",
     "The object responsible for calculating the cross section for two "
     "colliding dipole systems.",
     &DipoleEventHandler::theXSecFn, true, false, true, false, false);

  static Reference<DipoleEventHandler,ImpactParameterGenerator> interfaceBGen
    ("BGen",
     "The object responsible for generating the impact parameters.",
     &DipoleEventHandler::theBGen, true, false, true, false, false);

  static Reference<DipoleEventHandler,Emitter> interfaceEmitter
    ("Emitter",
     "The object responsible for generating and performing dipole emissions "
     "of gluons.",
     &DipoleEventHandler::theEmitter, true, false, true, false, false);

  static Reference<DipoleEventHandler,Swinger> interfaceSwinger
    ("Swinger",
     "The object responsible for generating and performing dipole swings.",
     &DipoleEventHandler::theSwinger, true, false, true, true, false);

  static Reference<DipoleEventHandler,EventFiller> interfaceEventFiller
    ("EventFiller",
     "The object responsible for filling an event with final state gluons.",
     &DipoleEventHandler::theEventFiller, true, false, true, false, false);

  static RefVector<DipoleEventHandler,DipoleAnalysisHandler>
    interfaceAnalysisHandlers
    ("AnalysisHandlers",
     "A list of analysis to be performed in the presample phase.",
     &DipoleEventHandler::analyses, -1, true, false, true, false, false);


  static Parameter<DipoleEventHandler,InvEnergy> interfaceBaryonSize
    ("BaryonSize",
     "The typical size of a baryon to be used by wave functions not "
     "defining their own. If zero, <interface>RMax</interface> will "
     "be used instead.",
     &DipoleEventHandler::theBaryonSize, InvGeV, 0.0*InvGeV, 0.0*InvGeV, 0*InvGeV,
     true, false, Interface::lowerlim);

}

