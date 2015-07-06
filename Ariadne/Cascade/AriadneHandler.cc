// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AriadneHandler class.
//

#include "AriadneHandler.h"
#include "ConsistencyChecker.h"
#include "QCDDipoleFinder.h"
#include "EMDipoleFinder.h"
#include "ReweightBase.h"
#include "DipoleState.h"
#include "EmitterBase.h"
#include "BornCheckerBase.h"
#include "ScaleSetter.h"
#include "DISFinder.h"
#include "ResonanceFinder.h"
#include "Ariadne/Cascade/Models/RemnantModel.h"
#include "Ariadne/Cascade/Models/ColourResonanceModel.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/RefVector.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/EventRecord/Collision.h"
#include "ThePEG/EventRecord/SubProcess.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "ThePEG/Handlers/Hint.h"
#include "ThePEG/Handlers/XComb.h"
#include "ThePEG/Utilities/UtilityBase.h"
#include "ThePEG/Utilities/DescriptionList.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Cuts/Cuts.h"
#include "ThePEG/Utilities/EnumIO.h"
#include "ThePEG/Utilities/Current.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Utilities/DebugItem.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/StandardModel/StandardModelBase.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include "Models/FSGluonEmission.h"
#include "Models/FSQQEmission.h"
#include "Models/DipoleSwing.h"
#include "ThePEG/Analysis/FactoryBase.h"
#ifndef LWH_AIAnalysisFactory_H
#ifndef LWH 
#define LWH ThePEGLWH
#endif
#include "ThePEG/Analysis/LWH/AnalysisFactory.h"
#endif


using namespace Ariadne5;

AriadneHandler::AriadneHandler()
  : theRunningCoupling(simpleRunning), theAlpha0(0.2), theLambdaQCD(0.22*GeV),
    scaleFactor(1.0), thePTCut(0.6*GeV), theAlphaEM0(1.0/137.0),
    thePTCutEM(0.6*GeV), theNCol(8), theNFlav(5), thePhotonEmissions(false),
    theSoftMu(0.6*GeV), theSoftAlpha(1.0), theHardAlpha(-1.0), theBeta(2.0),
    theMaxEmissions(0), suspendConsistencyChecks(false), thePurgeStrategy(neverpurge),
    thePurgeFactor(0.0) {}

AriadneHandler::~AriadneHandler() {}

IBPtr AriadneHandler::clone() const {
  return new_ptr(*this);
}

IBPtr AriadneHandler::fullclone() const {
  return new_ptr(*this);
}

bool AriadneHandler::preInitialize() const {
  if ( CascadeHandler::preInitialize() ) return true;
  return runningCoupling() == internalRunning && !internalAlphaS();
}

void AriadneHandler::doinit() throw(InitException) {
  CascadeHandler::doinit();
  if ( runningCoupling() != internalRunning || internalAlphaS() ) return;
  theInternalAlphaS = dynamic_ptr_cast<Ptr<AlphaSBase>::pointer>
    (generator()->preinitCreate("ThePEG::O1AlphaS", fullName() + "/AlphaS",
				"O1AlphaS.so"));
  ostringstream os;
  os << ounit(lambdaQCD(), GeV);
  generator()->preinitInterface(theInternalAlphaS,
				"LambdaQCD", "set", os.str());
  generator()->preinitInterface(theInternalAlphaS, "LambdaFlav", "set", "5");
  theInternalAlphaS->update();

  // Remove duplicate emitters.
  vector<EmitterPtr> cleaned;
  set<const ClassDescriptionBase *> cleanset;
  for ( int i = 0, N = emitters().size(); i < N; ++i ) {
    const ClassDescriptionBase * cd =
      DescriptionList::find(typeid(*emitters()[i]));
    if ( cleanset.find(cd) != cleanset.end() ) continue;
    cleaned.push_back(emitters()[i]);
    cleanset.insert(cd);
  }
  theEmitters.swap(cleaned);

}

void AriadneHandler::doinitrun() {
  static DebugItem histem("Ariadne5::HistEm", 6);
  CascadeHandler::doinitrun();
  theFlavourThresholds.clear();
  vector<Energy2> thrsh;
  switch ( runningCoupling() ) {
  case noRunning:
    return;
  case simpleRunning:
  case externalRunning:
    thrsh = SM().alphaSPtr()->flavourThresholds();
    break;
  case internalRunning:
    thrsh = internalAlphaS()->flavourThresholds();
  }
  for ( int i = 0, N = thrsh.size(); i < N; ++i ) thrsh[i] /= scaleFactor;
  theFlavourThresholds.insert(thrsh.begin(), thrsh.end());

  if ( histem ) {
    generator()->histogramFactory()->initrun();
    generator()->histogramFactory()->registerClient(this);
    generator()->histogramFactory()->mkdir("/Ariadne5Debug");
    histall = generator()->histogramFactory()->createHistogram1D
      ("/Ariadne5Debug/All", 100, -1.0, 9.0);
    histswing = generator()->histogramFactory()->createHistogram1D
      ("/Ariadne5Debug/Swings", 100, -1.0, 9.0);
    histglue = generator()->histogramFactory()->createHistogram1D
      ("/Ariadne5Debug/Glue", 100, -1.0, 9.0);
    histqq = generator()->histogramFactory()->createHistogram1D
      ("/Ariadne5Debug/QQ", 100, -1.0, 9.0);
    histlam = generator()->histogramFactory()->createHistogram1D
      ("/Ariadne5Debug/Lambda", 100, -1.0, 9.0);
  }
}

void AriadneHandler::cascade() {
  static DebugItem histem("Ariadne5::HistEm", 6);
  static DebugItem tupleswing("Ariadne5::SwingTuple", 6);
  static ofstream swingtuple;
  if ( tupleswing ) {
    static bool isopen = false;
    if ( !isopen ) {
      string filename = CurrentGenerator::current().filename() + "-swing.tuple";
      swingtuple.open(filename.c_str());
      isopen = true;
    }
  }

  Current<AriadneHandler> current(this);

  DipoleStatePtr state;

  tSubProPtr sub;
  tPVector final;

  Energy rhomax = ZERO;
  //  Energy rhomin = pTCut();
  Energy rhomin = ZERO;
  bool perf = false;
  int emnbr = 0;

  // Check if the state has already been generated in reweightCKKW
  if ( theCKKWMap.count(lastXCombPtr()) ){
    CKKWState ckkw = theCKKWMap[lastXCombPtr()];
    sub = subProcess();
    state = ckkw.history->state();
    checkState(*state);
    rhomax = startingScale(*state);
    while ( emnbr == 0 &&
	    ( rhomax = state->select(rhomin, rhomax) ) > rhomin ) {
      SaveDipoleState backup(state);
      if ( state->perform() ) {
	if ( !ckkw.maxMult && checkTreeState(*state) && passCuts(state) ) {
	  theCKKWMap.clear();
	  throw Veto();
	}
	perf = true;
	emnbr++;
      } else {
	state = backup.revert();
	state->selected()->dipole->touch();
      }
    }
  }
  else {

    tCollPtr coll = eventHandler()->currentCollision();
    state = new_ptr(DipoleState(coll->incoming()));
    sub = coll->primarySubProcess();
    if ( sub->decayed() ) sub = tSubProPtr();

    // Decide whether we are cascadeing a whole sub-process or if we
    // only need to deal with final-state shower.
    if ( ( hint().tagged() || !sub ) && !tagged().empty() ) {
      final = tagged();
    }
    else if ( sub ) {
      state->setup(*sub);
    }
    else
      return;
    if ( !final.empty() ) {
      state->setup(set<tPPtr>(final.begin(), final.end()));
    }
    checkState(*state);
    if ( !state->checkIntegrity() )
      Throw<Exception>()
	<< "Ariadne failed to setup dipole system. This is a serious error,"
	<< "Please inform the author." << Exception::runerror;

    rhomax = startingScale(*state);
  }
  theCKKWMap.clear();

  if ( tupleswing ) {
    swingtuple << "# " << setw(4) << CurrentGenerator::current().currentEventNumber()
	       << setw(8) << "rho"
	       << setw(14) << "lam/dip"
	       << setw(14) << "lambda"
	       << setw(6) << "cross"
	       << setw(6) << "below"
	       << setw(14) << "folding"
	       << setw(6) << "geno"
	       << setw(6) << "emno" << endl;
    pair<double,int> lam = state->lambdaMeasure(sqr(pTCut()));
    int cross = state->crossings();
    swingtuple << setw(14) << log10(rhomax/GeV)
	       << setw(14) << lam.first/lam.second
	       << setw(14) << lam.first
	       << setw(6) << cross
	       << setw(6) << state->gluonsBelow()
	       << setw(14) << state->folding() << endl;
  }

  if ( purgeStrategy() == onlybefore || purgeStrategy() > neverpurge )
    state->purgeGluons(pTCut()*purgeFactor());

  if ( tupleswing ) {
    pair<double,int> lam = state->lambdaMeasure(sqr(pTCut()));
    int cross = state->crossings();
    swingtuple << setw(14) << log10(rhomax/GeV)
	       << setw(14) << lam.first/lam.second
	       << setw(14) << lam.first
	       << setw(6) << cross
	       << setw(6) << state->gluonsBelow()
	       << setw(14) << state->folding() << endl;
  }

  while ( ( rhomax = state->select(rhomin, rhomax) ) > rhomin  && 
	  ( maxEmissions() == 0 || emnbr < maxEmissions() ) ) {
    SaveDipoleState backup(state);
    if ( state->perform() ) {
      perf = true;
      emnbr++;
      if ( tupleswing ) {
	double lrho = log10(state->selected()->rho/GeV);
	pair<double,int> lam = state->lambdaMeasure(sqr(pTCut()));
	int cross = state->crossings();
	swingtuple << setw(14) << lrho
		   << setw(14) << lam.first/lam.second
		   << setw(14) << lam.first
		   << setw(6) << cross
		   << setw(6) << state->gluonsBelow()
		   << setw(14) << state->folding()
		   << setw(6) << state->selected()->geno
		   << setw(6) << state->selected()->emno;
	if ( DipoleSwing * em = dynamic_cast<DipoleSwing*>((Emission*)(state->selected())) ) {
	  swingtuple << "s" << setw(6) << state->index(em->dipoles.first) << setw(6)
		     << state->index(em->dipoles.second);
	}
	else if ( dynamic_cast<FSGluonEmission*>((Emission*)(state->selected())) )
	  swingtuple << "g" << setw(6)<< state->index(state->selected()->dipole);
	else if ( dynamic_cast<FSQQEmission*>((Emission*)(state->selected())) )
	  swingtuple << "q" << setw(6) << state->index(state->selected()->dipole);

	swingtuple << endl;
      }
      if ( histem ) {
	double w = 1.0/double(state->activeDipoles().size());
	double lrho = log10(state->selected()->rho/GeV);
	histall->fill(lrho, w);
	histlam->fill(lrho, state->lambdaMeasure().first);
	if ( dynamic_cast<DipoleSwing*>((Emission*)(state->selected())) )
	  histswing->fill(lrho, w);
	else if ( dynamic_cast<FSGluonEmission*>((Emission*)(state->selected())) )
	  histglue->fill(lrho, w);
	else if ( dynamic_cast<FSQQEmission*>((Emission*)(state->selected())) )
	  histqq->fill(lrho, w);
	else
	  state->selected()->debug();
      }

      if ( purgeStrategy() == everystep ) state->purgeGluons(pTCut()*purgeFactor());

    } else {
      state = backup.revert();
      state->selected()->dipole->touch();
    }
    if ( !state->checkIntegrity() ) {
      throw IntegretyException() << "AriadneHandler::cascade "
        << "The dipole state is not self consistent."
        << Exception::eventerror;
    }
  }

  if ( perf && ( purgeStrategy() == onlyafter || purgeStrategy() == beforeandafter ) )
     state->purgeGluons(pTCut()*purgeFactor());
  if ( tupleswing ) {
    pair<double,int> lam = state->lambdaMeasure(sqr(pTCut()));
    int cross = state->crossings();
    swingtuple << setw(14) << log10(pTCut()/GeV)
	       << setw(14) << lam.first/lam.second
	       << setw(14) << lam.first
	       << setw(6) << cross
	       << setw(6) << state->gluonsBelow()
	       << setw(14) << state->folding() << endl;
  }


  if ( final.empty() ) {
    if ( perf ) state->fill(*newStep());
    sub->decayed(true);
  } else {
    if ( perf ) state->fill(*newStep());
  }

}

double AriadneHandler::reweightCKKW(int minMult, int maxMult) {
  if(minMult == maxMult){
    return 1.0;
  }

  CKKWState ckkw;
  DipoleStatePtr state = new_ptr(DipoleState());
  tXCPtr lastXC = lastXCombPtr();

  int outgoing = lastXC->subProcess()->outgoing().size();
  if ( outgoing < minMult || outgoing > maxMult ) {
        throw CKKWMultiplicityException() 
          << "Ariadne5::CascadeHandler::reweightCKKW "
          << "Number of outgoing particles out of range."
          << Exception::eventerror;
  }
  int steps = outgoing - minMult;
  ckkw.maxMult = (outgoing == maxMult);
  state->setup(*(lastXC->subProcess()));
  checkState(*state);
  if ( !state->checkIntegrity() )
    Throw<Exception>()
      << "Ariadne failed to setup dipole system. This is a serious error,"
      << "Please inform the author." << Exception::runerror;

  ckkw.history = HistoryPtr(new History(steps, state));
  if ( !ckkw.history->select() ) return 0.0;

  double w = ckkw.history->weight(lastXC->lastAlphaS(), lastXC->lastScale());

  if ( w > 0.0 ) theCKKWMap[lastXC] = ckkw;
  return w;
}

bool AriadneHandler::passCuts(tcDipoleStatePtr state) {
  /* *** ATTENTION ***
  tCutsPtr cuts = lastCutsPtr();
  Lorentz5Momentum ph = state->hardFS().momentum();
  cuts->initSubProcess(ph.mass2(), ph.rapidity());

  tcPDVector pdata;
  vector< LorentzMomentum > p;
  typedef HardSubSys::PartonSet PartonSet;
  const PartonSet & active = state->hardSubSys().coloured();
  const PartonSet & produced = state->hardSubSys().produced();

  for(PartonSet::const_iterator it = active.begin(); it != active.end();
      it++){
    pdata.push_back((*it)->dataPtr());
    p.push_back((*it)->momentum());
  }
  for(PartonSet::const_iterator it = produced.begin(); it != produced.end();
      it++){
    pdata.push_back((*it)->dataPtr());
    p.push_back((*it)->momentum());
  }
  LorentzRotation R(0.0, 0.0, - ph.z()/ph.e());
  Utilities::transform(p, R);

  return cuts->passCuts(pdata, p, state->particles().first->dataPtr(),
      state->particles().second->dataPtr());
  */ return false;
}

double AriadneHandler::alphaS(Energy2 scale) const {
  scale *= scaleFactor;
  int Nf = SM().Nf(scale);
  switch ( runningCoupling() ) {
  case noRunning:
    return alpha0();
  case simpleRunning:
    return 12.0*Constants::pi/
      ((33.0 - 2.0*min(int(SM().Nf(scale)),5))*log(scale/sqr(lambdaQCD())));
  case internalRunning:
    return internalAlphaS()->value(scale, SM());
  case externalRunning:
    return SM().alphaS(scale);
  }
  return ZERO*Nf;
}

Energy AriadneHandler::checkBornState(const DipoleState & ds) const {
  Energy scale = ZERO;
  for ( int i = 0, N = theBornCheckers.size(); i < N; ++i ) {
    Energy mu = theBornCheckers[i]->check(ds);
    if ( mu > ZERO ) return mu;
    if ( scale == ZERO ) scale = mu;
  }
  if ( scale == ZERO )
    Throw<CKKWBornException>()
      << "Could not reconstruct the given event in the CKKW-L algorithm. "
      << *(eventHandler()->currentEvent()) << Exception::runerror;
  return scale;
}

bool AriadneHandler::checkTreeState(const DipoleState & ds) const {
  for ( int i = 0, N = theBornCheckers.size(); i < N; ++i )
    if ( theBornCheckers[i]->checkTree(ds) ) return true;
  return false;
}

vector<tQCDPtr> AriadneHandler::findQCDDipoles(DipoleState & state) const {
  return theQCDFinder->findDipoles(state);
}

vector<tEMDipPtr> AriadneHandler::findEMDipoles(DipoleState & state) const {
  return theEMFinder? theEMFinder->findDipoles(state): vector<tEMDipPtr>();
}

Energy AriadneHandler::startingScale(const DipoleState & state) const {
    return theScaleSetter->scale(state);
}



bool  AriadneHandler::checkState(DipoleState & state, tcEmPtr e) {
  if ( !consistency ) return true;
  if ( !e ) suspendConsistencyChecks = false;
  else if ( suspendConsistencyChecks ) return true;
  if ( consistency->check(state, e) ) return true;
  if ( !e ) suspendConsistencyChecks = true;
  return false;
}

pair<tRemParPtr,tRemParPtr>
AriadneHandler::findDISLeptons(SubProcess & sub, DipoleState & state) const {
  return theDISFinder? theDISFinder->findDISLeptons(sub, state):
    pair<tRemParPtr,tRemParPtr>();
}

pair<tRemParPtr,tRemParPtr>
AriadneHandler::findDISQuarks(pair<tRemParPtr,tRemParPtr> leptons,
			      SubProcess & sub, DipoleState & state) const {
  return theDISFinder? theDISFinder->findDISQuarks(leptons, sub, state):
    pair<tRemParPtr,tRemParPtr>();
}

tPVector AriadneHandler::resonances(SubProcess & sub) const {
  return theResonanceFinder->resonances(sub);
}

void AriadneHandler::persistentOutput(PersistentOStream & os) const {
  os << oenum(theRunningCoupling) << theAlpha0 << ounit(theLambdaQCD, GeV)
     << theInternalAlphaS << scaleFactor << theFlavourThresholds.size()
     << ounit(thePTCut, GeV) << theAlphaEM0 << ounit(thePTCutEM, GeV)
     << theNCol << theNFlav << thePhotonEmissions << ounit(theSoftMu, GeV)
     << theSoftAlpha << theHardAlpha << theBeta << theReweighters
     << theMaxEmissions << theEmitters << theBornCheckers << theScaleSetter
     << theDISFinder << theResonanceFinder << theQCDFinder << theEMFinder
     << theColourResonanceModel << theRemnantModel
     << consistency << suspendConsistencyChecks << oenum(thePurgeStrategy) << thePurgeFactor;
  for ( set<Energy2>::const_iterator it = theFlavourThresholds.begin();
	it != theFlavourThresholds.end(); ++it ) os << ounit(*it, GeV2);
}

void AriadneHandler::persistentInput(PersistentIStream & is, int) {
  int size = 0;
  is >> ienum(theRunningCoupling) >> theAlpha0 >> iunit(theLambdaQCD, GeV)
     >> theInternalAlphaS >> scaleFactor >> size
     >> iunit(thePTCut, GeV) >> theAlphaEM0 >> iunit(thePTCutEM, GeV)
     >> theNCol >> theNFlav >> thePhotonEmissions >> iunit(theSoftMu, GeV)
     >> theSoftAlpha >> theHardAlpha >> theBeta >> theReweighters
     >> theMaxEmissions >> theEmitters >> theBornCheckers >> theScaleSetter
     >> theDISFinder >> theResonanceFinder >> theQCDFinder >> theEMFinder
     >> theColourResonanceModel >> theRemnantModel
     >> consistency >> suspendConsistencyChecks >> ienum(thePurgeStrategy) >> thePurgeFactor;
  theFlavourThresholds.clear();
  Energy2 t = ZERO;
  for ( int i = 0; i < size; ++i ) {
    is >> iunit(t, GeV2);
    theFlavourThresholds.insert(t);
  }
}

DescribeClass<AriadneHandler,ThePEG::CascadeHandler> 
describeAriadne5("Ariadne5::AriadneHandler", "libAriadne5.so");

Energy AriadneHandler::minPTCut() const {
  return ( runningCoupling() == simpleRunning ||
	  ( runningCoupling() == internalRunning && !theInternalAlphaS ) )?
    lambdaQCD()/sqrt(scaleFactor): 0.0*GeV;
}

Energy AriadneHandler::maxLambdaQCD() const {
  return ( runningCoupling() == simpleRunning ||
	  ( runningCoupling() == internalRunning && !theInternalAlphaS ) )?
    pTCut()*sqrt(scaleFactor): Constants::MaxEnergy;
}

void AriadneHandler::Init() {

  static ClassDocumentation<AriadneHandler> documentation
    ("The AriadneHandler class administers the Ariadne dipole "
     "cascade.",
     "Parton cascades performed by Ariadne\\cite{Lav05} according to the "
     "Dipole Cascade Model\\cite{Gustafson:1986db,Gustafson:1988rq,"
     "Andersson:1989gp,Andersson:1990ki}.",
     "\\bibitem{Lav05}"
     "Nils Lavesson and Leif L\\\"onnblad, Preprint in preparation.\n"
     "\\bibitem{Gustafson:1986db}"
     "G\\\"osta Gustafson, Phys.~Lett.~{\\bf B175} (1986) 453.\n"
     "\\bibitem{Gustafson:1988rq}"
     "G\\\"osta Gustafson and Ulf Pettersson, "
     "Nucl.~Phys.~{\\bf B306} (1988) 746.\n"
     "\\bibitem{Andersson:1989gp}"
     "Bo Andersson, et al., Z.~Phys.~{\\bf C43} (1989) 625.\n"
     "\\bibitem{Andersson:1990ki}"
     "Bo Andersson, et al., Nucl.~Phys.~{\\bf B339} (1990) 393.");

  static Switch<AriadneHandler,RunningOption> interfaceRunningCoupling
    ("RunningCoupling",
     "Strategy for handling \\f$\\alpha_S\\f$.",
     &AriadneHandler::theRunningCoupling, simpleRunning, true, false);
  static SwitchOption interfaceRunningCouplingSimpleRunning
    (interfaceRunningCoupling,
     "SimpleRunning",
     "Use a one loop running coupling with \\f$\\Lambda_{QCD}\\f$ given "
     "by <interface>LambdaQCD</interface>.",
     simpleRunning);
  static SwitchOption interfaceRunningCouplingConstant
    (interfaceRunningCoupling,
     "Constant",
     "Use a constant coupling given by <interface>Alpha0</interface>.",
     noRunning);
  static SwitchOption interfaceRunningCouplingExternalRunning
    (interfaceRunningCoupling,
     "ExternalRunning",
     "Use whatever coupling is specified by the current StandardModelBase "
     "object.",
     externalRunning);
  static SwitchOption interfaceRunningCouplingInternalRUnning
    (interfaceRunningCoupling,
     "InternalRunning",
     "Use the coupling specified by <interface>InternalAlphaS</interface>.",
     internalRunning);

  static Parameter<AriadneHandler,double> interfaceAlpha0
    ("Alpha0",
     "The constant \\f$\\alpha_S\\f$ to use if "
     "<interface>RunningCoupling</interface> is set to Constant.",
     &AriadneHandler::theAlpha0, 0.2, 0.0, 0,
     true, false, Interface::lowerlim);

  static Parameter<AriadneHandler,double> interfaceAlphaEM0
    ("AlphaEM0",
     "The constant \\f$\\alpha_{EM}\\f$ to use. If zero, use whatever "
    "is specified in the current StandardModelBase object.",
     &AriadneHandler::theAlphaEM0, 1.0/137.0, 0.0, 0,
     true, false, Interface::lowerlim);

  static Parameter<AriadneHandler,Energy> interfaceLambdaQCD
    ("LambdaQCD",
     "The \\f$\\Lambda_{QCD}\\f$ to use in the one loop running "
     "\\f$\\alpha_S\\f$ if <interface>RunningCoupling</interface> "
     "is set to Running.",
     &AriadneHandler::theLambdaQCD, GeV, 0.22*GeV, 0.0*GeV,
     Constants::MaxEnergy,
     true, false, Interface::limited,
     0, 0, 0, &AriadneHandler::maxLambdaQCD, 0);


  static Parameter<AriadneHandler,double> interfaceScaleFactor
    ("ScaleFactor",
     "Scale factor used to multiply the emission scales in the argument of "
     "\\f$\alpha_S\\f$.",
     &AriadneHandler::scaleFactor, 1.0, 0.0, 0,
     true, false, Interface::lowerlim);


  static Parameter<AriadneHandler,Energy> interfacePTCut
    ("PTCut",
     "The cutoff in invariant transverse momentum for QCD emissions.",
     &AriadneHandler::thePTCut, GeV, 0.6*GeV, 0.0*GeV, 0*GeV,
     true, false, Interface::lowerlim,
     0, 0, &AriadneHandler::minPTCut, 0, 0);

  static Parameter<AriadneHandler,Energy> interfacePTCutEM
    ("PTCutEM",
     "The cutoff in invariant transverse momentum for QED emissions.",
     &AriadneHandler::thePTCutEM, GeV, 0.6*GeV, 0.0*GeV, 0*GeV,
     true, false, Interface::lowerlim);

  static Parameter<AriadneHandler,int> interfaceDipoleColours
    ("DipoleColours",
     "The number of differently coloured dipoles possible. This should "
     "normally be 8, but may be varied to check the effects of colour "
     "reconnections.",
     &AriadneHandler::theNCol, 8, 3, 0,
     true, false, Interface::lowerlim);

  static Parameter<AriadneHandler,int> interfaceNFlav
    ("NFlav",
     "The number of possible flavours in a \f$g\to q\bar{q}\f$ splitting.",
     &AriadneHandler::theNFlav, 5, 0, 8,
     true, false, Interface::lowerlim);

  static Switch<AriadneHandler,bool> interfacePhotonEmissions
    ("PhotonEmissions",
     "Switches photon emission in the cascade on and off.",
     &AriadneHandler::thePhotonEmissions, false, true, false);
  static SwitchOption interfacePhotonEmissionsOn
    (interfacePhotonEmissions,
     "On",
     "Switch photon emission on",
     true);
  static SwitchOption interfacePhotonEmissionsOff
    (interfacePhotonEmissions,
     "Off",
     "Switch photon emission off.",
     false);

  static Parameter<AriadneHandler,Energy> interfaceSoftMu
    ("SoftMu",
     "The inverse extension of a hadron remnant used in the soft "
     "suppression mechanism. See also <interface>SoftAlphs</interface> and "
     "<interface>Beta</interface>",
     &AriadneHandler::theSoftMu, GeV, 0.6*GeV, 0.0*GeV, 0*GeV,
     true, false, Interface::lowerlim);

  static Parameter<AriadneHandler,double> interfaceSoftAlpha
    ("SoftAlpha",
     "The dimension of the extension of a hadron remnant used in the "
     "soft-suppression mechanism. See also <interface>SoftMu</interface> "
     "and <interface>Beta</interface>.",
     &AriadneHandler::theSoftAlpha, 1.0, 0.0, 0,
     true, false, Interface::lowerlim);

  static Parameter<AriadneHandler,double> interfaceHardAlpha
    ("HardAlpha",
     "The dimension of the extension of a hard remnant used in the "
     "soft-suppression mechanism for radiation off a scattered quark in "
     "DIS at small \\f$Q^2\\f$. If set negative, the value of "
     "<interface>SoftAlpha</interface> will be used instead. See also "
     "<interface>Beta</interface>.",
     &AriadneHandler::theHardAlpha, -1.0, -1.0, 0,
     true, false, Interface::lowerlim);

  static Parameter<AriadneHandler,double> interfaceBeta
    ("Beta",
     "The power in the suppression of radiation from extended dipoles. "
     "The original soft suppression model used a sharp cutoff in the "
     "transverse-momentum--rapidity space of an emitted gluon. This parameter "
     "is used to allow emissions with larger transverse momentum according to "
     "\\f$P(p_\\perp^2>p_{\\perp cut}^2="
     "\\left(\\frac{p_{\\perp cut}^2}{p_\\perp^2}\\right)^\\beta\\f$. if "
     "negative, the sharp cutoff is retained.",
     &AriadneHandler::theBeta, 2.0, -1.0, 0,
     true, false, Interface::lowerlim);

  static RefVector<AriadneHandler,Ariadne5::ReweightBase> interfaceReweighters
    ("Reweighters",
     "A vector of objects implementing reweightings of basic dipole "
     "emissions. Each dipole emission will be reweighted.",
     &AriadneHandler::theReweighters, -1, true, false, true, false, false);

  static Parameter<AriadneHandler,int> interfaceMaxEmissions
    ("MaxEmissions",
     "This number specifies the maximum number of emissions from the "
     "cascade. If it is set to zero an unlimited number is allowed. "
     "This parameter should only be used for debugging purposes.",
     &AriadneHandler::theMaxEmissions, 0, 0, 0,
     true, false, Interface::lowerlim);

  static Reference<AriadneHandler,AlphaSBase> interfaceInternalAlphaS
    ("InternalAlphaS",
     "An internal AlphaSBase object to be used if "
     "<interface>RunningCoupling</interface> is set to InternalRunning. "
     "If no such object is given, A O1AlphaS object will be created and "
     "assigned in the initialization with <interface>LambdaQCD</interface> "
     "used as lambda for five flavours.",
     &AriadneHandler::theInternalAlphaS, true, false, true, true, false);

  static RefVector<AriadneHandler,EmitterBase> interfaceEmitters
    ("Emitters",
     "The vector of EmittorBase objects responsible for generatong and "
     "performing emissions according to the Dipole Cascade Model and its "
     "extentions. For each dipole, all Emitters are tried in turn if to see "
     "if they are able to model emissions. Only one EmittorBase object of "
     "each SubClass will be used (earlier ones will override later ones "
     "in the vector).",
     &AriadneHandler::theEmitters, -1, true, false, true, false, false);


  static RefVector<AriadneHandler,BornCheckerBase> interfaceBornCheckers
    ("BornCheckers",
     " A vector of BornCheckerBase objects which are used in the CKKW-L "
     "algorithm to check if a reclustered dipole state corresponds to a "
     "reasonable Born-level state (the lowest jet-multiplicity state in a "
     "CKKW-L merging). At least one object must be assigned in order for "
     "the CKKW-L algorithm to work.",
     &AriadneHandler::theBornCheckers, -1, true, false, true, false, false);

  static Reference<AriadneHandler,ScaleSetter> interfaceScaleSetter
    ("ScaleSetter",
     "The object responsible for choosing the starting scale of the "
     "Ariadne dipole shower.",
     &AriadneHandler::theScaleSetter, true, false, true, false, true);

  static Reference<AriadneHandler,DISFinder> interfaceDISFinder
    ("DISFinder",
     "The object responsible for identifying DIS-like event.",
     &AriadneHandler::theDISFinder, true, false, true, true, false);

  static Reference<AriadneHandler,ResonanceFinder> interfaceResonanceFinder
    ("ResonanceFinder",
     "The object responsible for identifying DIS-like event.",
     &AriadneHandler::theResonanceFinder, true, false, true, false, false);

  static Reference<AriadneHandler,QCDDipoleFinder> interfaceQCDFinder
    ("QCDFinder",
     "The Object responsible for identifying QCD Diploles.",
     &AriadneHandler::theQCDFinder, true, false, true, false, false);

  static Reference<AriadneHandler,EMDipoleFinder> interfaceEMFinder
    ("EMFinder",
     "The Object responsible for identifying electro-magnetic Diploles. "
     "If null, all electro-magnetic radiation will be switched off.",
     &AriadneHandler::theEMFinder, true, false, true, true, false);

  static Reference<AriadneHandler,ColourResonanceModel>
    interfaceColourResonanceModel
    ("ColourResonanceModel",
     "The object responsible for radiating from decay products from "
     "coloured resonances.",
     &AriadneHandler::theColourResonanceModel, true, false, true, true, false);

  static Reference<AriadneHandler,RemnantModel> interfaceRemnantModel
    ("RemnantModel",
     "The object responsible for radiating from remnants.",
     &AriadneHandler::theRemnantModel, true, false, true, false, false);

  static Reference<AriadneHandler,ConsistencyChecker> interfaceConsistencyChecker
    ("ConsistencyChecker",
     "The object responsible for checking the consistency of a DipoleState.",
     &AriadneHandler::consistency, true, false, true, true, false);

  static Switch<AriadneHandler,PurgeStrategy> interfacePurgeStrategy
    ("PurgeStrategy",
     "The strategy for purging gluons with transverse momentum less than a factor <interface>PurgeFactor</interface> times the cutoff.",
     &AriadneHandler::thePurgeStrategy, neverpurge, true, false);
  static SwitchOption interfacePurgeStrategyOnlyBefore
    (interfacePurgeStrategy,
     "OnlyBefore",
     "Gluons are only purged before the cascade.",
     onlybefore);
  static SwitchOption interfacePurgeStrategyOnlyAfter
    (interfacePurgeStrategy,
     "OnlyAfter",
     "Gluons are only purged after the cascade.",
     onlyafter);
  static SwitchOption interfacePurgeStrategyNever
    (interfacePurgeStrategy,
     "Never",
     "Gluon purging is switched off.",
     neverpurge);
  static SwitchOption interfacePurgeStrategyBeforeAndAfter
    (interfacePurgeStrategy,
     "BeforeAndAfter",
     "Gluons are purged before and after the cascade.",
     beforeandafter);
  static SwitchOption interfacePurgeStrategyEveryStep
    (interfacePurgeStrategy,
     "EveryStep",
     "Gluons are purged before (and after) each step of the cascade.",
     everystep);

  static Parameter<AriadneHandler,double> interfacePurgeFactor
    ("PurgeFactor",
     "The factor used to determine how far below the cutoff a gluon may be before "
     "it is purged. If the transverse momentum of a gluon is less than this factor "
     "times the cutoff it will be purged. Cf. <interface>PurgeStrategy</interface>.",
     &AriadneHandler::thePurgeFactor, 0.0, 0.0, 0,
     true, false, Interface::lowerlim);

  interfacePTCut.rank(10);
  interfaceLambdaQCD.rank(9);
  interfaceNFlav.rank(8);
  interfaceSoftMu.rank(7);
  interfaceSoftAlpha.rank(6);
  interfaceHardAlpha.rank(5);
  interfaceBeta.rank(4);

}


