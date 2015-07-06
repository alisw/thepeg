// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the CascadeHandler class.
//

#include "CascadeHandler.h"
#include "MECorrBase.h"
#include "ReweightBase.h"
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
#include "ThePEG/Cuts/Cuts.h"
#include "Ariadne/DipoleCascade/DipoleState.h"
#include "ThePEG/Utilities/EnumIO.h"
#include "ThePEG/Utilities/Current.h"
#include "ThePEG/Repository/EventGenerator.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "CascadeHandler.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace Ariadne {

CascadeHandler::CascadeHandler()
  : theRunningCoupling(simpleRunning), theAlpha0(0.2), theLambdaQCD(0.22*GeV),
    thePTCut(0.6*GeV), theAlphaEM0(1.0/137.0), thePTCutEM(0.6*GeV),
    theNCol(8), theNFlav(5), thePhotonEmissions(false), theSoftMu(0.6*GeV),
    theSoftAlpha(1.0), theHardAlpha(-1.0), theBeta(2.0),
    theMaxEmissions(0) {}

CascadeHandler::CascadeHandler(const CascadeHandler & x)
  : ThePEG::CascadeHandler(x), theRunningCoupling(x.theRunningCoupling),
    theAlpha0(x.theAlpha0), theLambdaQCD(x.theLambdaQCD),
    theInternalAlphaS(x.theInternalAlphaS),
    thePTCut(x.thePTCut), theAlphaEM0(x.theAlphaEM0),
    thePTCutEM(x.thePTCutEM), theNCol(x.theNCol), theNFlav(x.theNFlav),
    thePhotonEmissions(x.thePhotonEmissions), theSoftMu(x.theSoftMu),
    theSoftAlpha(x.theSoftAlpha), theHardAlpha(x.theHardAlpha),
    theBeta(x.theBeta), theReweighters(x.theReweighters),
    theMECorrectors(x.theMECorrectors), theMaxEmissions(x.theMaxEmissions) {}

CascadeHandler::~CascadeHandler() {}

bool CascadeHandler::preInitialize() const {
  if ( ThePEG::CascadeHandler::preInitialize() ) return true;
  return runningCoupling() == internalRunning && !internalAlphaS();
}

void CascadeHandler::doinit() throw(InitException) {
  ThePEG::CascadeHandler::doinit();
  if ( runningCoupling() != internalRunning || internalAlphaS() ) return;
  theInternalAlphaS = dynamic_ptr_cast<Ptr<AlphaSBase>::pointer>
    (generator()->preinitCreate("ThePEG::O1AlphaS", fullName() + "/AlphaS",
				"O1AlphaS.so"));
  ostringstream os;
  os << lambdaQCD()/GeV;
  generator()->preinitInterface(theInternalAlphaS,
				"LambdaQCD", "set", os.str());
  generator()->preinitInterface(theInternalAlphaS, "LambdaFlav", "set", "5");
  theInternalAlphaS->update();
}

void CascadeHandler::cascade() {

  Current<Ariadne::CascadeHandler> current(this);

  LorentzRotation rot;

  DipoleStatePtr state;

  tSubProPtr sub;
  tPVector final;

  Energy2 pt2max;
  Energy2 pt2min = sqr(pTCut());
  bool perf = false;
  int emnbr = 0;

  // Check if the state has already been generated in reweightCKKW
  if ( theCKKWMap.count(lastXCombPtr()) ){
    CKKWState ckkw = theCKKWMap[lastXCombPtr()];
    sub = subProcess();
    state = ckkw.state;
    rot = ckkw.rotation;
    pt2max = ckkw.pt2;
    while ( emnbr == 0 && ( pt2max = state->select(pt2min, pt2max) ) > pt2min ) {
      DipoleState::TranslationMap trans;
      DipoleStatePtr newstate = state->preclone(trans);
      if ( state->perform() ) {
	if(! ckkw.maxMult && passCuts(state)){
	  theCKKWMap.clear();
	  throw Veto();
	}
	perf = true;
	emnbr++;
      } else {
	newstate->postclone(trans);
	state = newstate;
	state->selected()->touch();
      }
    }
  }
  else {
    state = new_ptr(DipoleState(this));
    Energy2 stot = 0.0*GeV2;

    tCollPtr coll = eventHandler()->currentCollision();
    sub = coll->primarySubProcess();
    if ( sub->decayed() ) sub = tSubProPtr();

    // Decide whether we are cascadeing a whole sub-process or if we
    // only need to deal with final-state shower.
    if ( ( hint().tagged() || !sub ) && !tagged().empty() ) {
      final = tagged();
    }
    else if ( sub ) {
      rot = state->init(sub);
      stot = state->sTot();
    }
    else
      return;
    if ( !final.empty() ) {
      rot = Utilities::boostToCM(final.begin(), final.end());
      state->init(final);
      stot = Utilities::sumMomentum(final.begin(), final.end()).m2();
    }

    pt2max = stot/4.0;
  }
  theCKKWMap.clear();

  while ( ( pt2max = state->select(pt2min, pt2max) ) > pt2min  && 
      ( maxEmissions() == 0 || emnbr < maxEmissions() ) ) {
    DipoleState::TranslationMap trans;
    DipoleStatePtr newstate = state->preclone(trans);
    if ( state->perform() ) {
      perf = true;
      emnbr++;
    } else {
      newstate->postclone(trans);
      state = newstate;
      state->selected()->touch();
    }
    if(! state->checkIntegrety()){
      throw IntegretyException() << "Ariadne::CascadeHandler::cascade "
        << "The dipole state is not self consistant."
        << Exception::eventerror;
    }
  }

  rot.invert();
  if ( final.empty() ) {
    if ( perf ) state->fill(sub, rot, newStep());
    else sub->transform(rot);
    sub->decayed(true);
  } else {
    Utilities::transform(final.begin(), final.end(), rot);
    if ( perf ) state->fill(final, rot, newStep());
  }
}

double CascadeHandler::reweightCKKW(int minMult, int maxMult){
  if(minMult == maxMult){
    return 1.0;
  }

  CKKWState ckkw;
  ckkw.state = new_ptr(DipoleState(this));
  tXCPtr lastXC = lastXCombPtr();

  int outgoing = lastXC->subProcess()->outgoing().size();
  if(outgoing < minMult || outgoing > maxMult){
        throw CKKWMultiplicityException() 
          << "Ariadne::CascadeHandler::reweightCKKW "
          << "Number of outgoing particles out of range."
          << Exception::eventerror;
  }
  int steps = outgoing - minMult;
  ckkw.maxMult = (outgoing == maxMult);
  ckkw.rotation = ckkw.state->init(lastXC->subProcess());

  if(! ckkw.state->constructHistory(steps)){
    return 0.0;
  }
  double weight = ckkw.state->couplingProduct(steps) /
    pow(lastXC->lastAlphaS(), steps) * ckkw.state->PDFRatioProduct(steps);
  ckkw.pt2 = ckkw.state->constructedPT2();
  if(ckkw.state->sudakovVeto(steps)){
    return 0.0;
  }

  theCKKWMap[lastXC] = ckkw;
  return weight;
}

bool CascadeHandler::passCuts(tcDipoleStatePtr state){
  tCutsPtr cuts = lastCutsPtr();
  Lorentz5Momentum ph = state->hardSubSys().momentum();
  cuts->initSubProcess(ph.mass2(), ph.rapidity());

  tcPDVector pdata;
  vector< LorentzMomentum > p;
  typedef HardSubSys::PartonSet PartonSet;
  const PartonSet & active = state->hardSubSys().active();
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
}

void CascadeHandler::persistentOutput(PersistentOStream & os) const {
  os << oenum(theRunningCoupling) << theAlpha0 << ounit(theLambdaQCD, GeV)
     << theInternalAlphaS
     << ounit(thePTCut, GeV) << theAlphaEM0 << ounit(thePTCutEM, GeV)
     << theNCol << theNFlav << thePhotonEmissions << ounit(theSoftMu, GeV)
     << theSoftAlpha << theHardAlpha << theBeta << theMECorrectors << theReweighters
     << theMaxEmissions;
}

void CascadeHandler::persistentInput(PersistentIStream & is, int) {
  is >> ienum(theRunningCoupling) >> theAlpha0 >> iunit(theLambdaQCD, GeV)
     >> theInternalAlphaS
     >> iunit(thePTCut, GeV) >> theAlphaEM0 >> iunit(thePTCutEM, GeV)
     >> theNCol >> theNFlav >> thePhotonEmissions >> iunit(theSoftMu, GeV)
     >> theSoftAlpha >> theHardAlpha >> theBeta >> theMECorrectors >> theReweighters
     >> theMaxEmissions;
}

ClassDescription<CascadeHandler> CascadeHandler::initCascadeHandler;
// Definition of the static class description member.

Energy CascadeHandler::minPTCut() const {
  return ( runningCoupling() == simpleRunning ||
	  ( runningCoupling() == internalRunning && !theInternalAlphaS ) )?
    lambdaQCD(): 0.0*GeV;
}

Energy CascadeHandler::maxLambdaQCD() const {
  return ( runningCoupling() == simpleRunning ||
	  ( runningCoupling() == internalRunning && !theInternalAlphaS ) )?
    pTCut(): Constants::MaxEnergy;
}

void CascadeHandler::Init() {

  static ClassDocumentation<CascadeHandler> documentation
    ("The Ariadne::CascadeHandler class administers the Ariadne dipole "
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

  static Switch<CascadeHandler,RunningOption> interfaceRunningCoupling
    ("RunningCoupling",
     "Strategy for handling \\f$\\alpha_S\\f$.",
     &CascadeHandler::theRunningCoupling, simpleRunning, true, false);
  static SwitchOption interfaceRunningCouplingRunning
    (interfaceRunningCoupling,
     "Running",
     "Use a one loop running coupling with \\f$\\Lambda_{QCD}\\f$ given "
     "by <interface>LambdaQCD</interface>.",
     simpleRunning);
  static SwitchOption interfaceRunningCouplingConstant
    (interfaceRunningCoupling,
     "Constant",
     "Use a constant coupling given by <interface>Alpha0</interface>.",
     noRunning);
  static SwitchOption interfaceRunningCouplingExternal
    (interfaceRunningCoupling,
     "External",
     "Use whatever coupling is specified by the current StandardModelBase "
     "object.",
     externalRunning);

  static Parameter<CascadeHandler,double> interfaceAlpha0
    ("Alpha0",
     "The constant \\f$\\alpha_S\\f$ to use if "
     "<interface>RunningCoupling</interface> is set to Constant.",
     &CascadeHandler::theAlpha0, 0.2, 0.0, 0,
     true, false, Interface::lowerlim);

  static Parameter<CascadeHandler,double> interfaceAlphaEM0
    ("AlphaEM0",
     "The constant \\f$\\alpha_{EM}\\f$ to use. If zero, use whatever "
    "is specified in the current StandardModelBase object.",
     &CascadeHandler::theAlphaEM0, 1.0/137.0, 0.0, 0,
     true, false, Interface::lowerlim);

  static Parameter<CascadeHandler,Energy> interfaceLambdaQCD
    ("LambdaQCD",
     "The \\f$\\Lambda_{QCD}\\f$ to use in the one loop running "
     "\\f$\\alpha_S\\f$ if <interface>RunningCoupling</interface> "
     "is set to Running.",
     &CascadeHandler::theLambdaQCD, GeV, 0.22*GeV, 0.0*GeV,
     Constants::MaxEnergy,
     true, false, Interface::limited,
     0, 0, 0, &CascadeHandler::maxLambdaQCD, 0);

  static Parameter<CascadeHandler,Energy> interfacePTCut
    ("PTCut",
     "The cutoff in invariant transverse momentum for QCD emissions.",
     &CascadeHandler::thePTCut, GeV, 0.6*GeV, 0.0*GeV, 0*GeV,
     true, false, Interface::lowerlim,
     0, 0, &CascadeHandler::minPTCut, 0, 0);

  static Parameter<CascadeHandler,Energy> interfacePTCutEM
    ("PTCutEM",
     "The cutoff in invariant transverse momentum for QED emissions.",
     &CascadeHandler::thePTCutEM, GeV, 0.6*GeV, 0.0*GeV, 0*GeV,
     true, false, Interface::lowerlim);

  static Parameter<CascadeHandler,int> interfaceDipoleColours
    ("DipoleColours",
     "The number of differently coloured dipoles possible. This should "
     "normally be 8, but may be varied to check the effects of colour "
     "reconnections.",
     &CascadeHandler::theNCol, 8, 3, 0,
     true, false, Interface::lowerlim);

  static Parameter<CascadeHandler,int> interfaceNFlav
    ("NFlav",
     "The number of possible flavours in a \f$g\to q\bar{q}\f$ splitting.",
     &CascadeHandler::theNFlav, 5, 0, 8,
     true, false, Interface::lowerlim);

  static Switch<CascadeHandler,bool> interfacePhotonEmissions
    ("PhotonEmissions",
     "Switches photon emission in the cascade on and off.",
     &CascadeHandler::thePhotonEmissions, false, true, false);
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

  static Parameter<CascadeHandler,Energy> interfaceSoftMu
    ("SoftMu",
     "The inverse extension of a hadron remnant used in the soft "
     "suppression mechanism. See also <interface>SoftAlphs</interface> and "
     "<interface>Beta</interface>",
     &CascadeHandler::theSoftMu, GeV, 0.6*GeV, 0.0*GeV, 0*GeV,
     true, false, Interface::lowerlim);

  static Parameter<CascadeHandler,double> interfaceSoftAlpha
    ("SoftAlpha",
     "The dimension of the extension of a hadron remnant used in the "
     "soft-suppression mechanism. See also <interface>SoftMu</interface> "
     "and <interface>Beta</interface>.",
     &CascadeHandler::theSoftAlpha, 1.0, 0.0, 0,
     true, false, Interface::lowerlim);

  static Parameter<CascadeHandler,double> interfaceHardAlpha
    ("HardAlpha",
     "The dimension of the extension of a hard remnant used in the "
     "soft-suppression mechanism for radiation off a scattered quark in "
     "DIS at small \\f$Q^2\\f$. If set negative, the value of "
     "<interface>SoftAlpha</interface> will be used instead. See also "
     "<interface>Beta</interface>.",
     &CascadeHandler::theHardAlpha, -1.0, -1.0, 0,
     true, false, Interface::lowerlim);

  static Parameter<CascadeHandler,double> interfaceBeta
    ("Beta",
     "The power in the suppression of radiation from extended dipoles. "
     "The original soft suppression model used a sharp cutoff in the "
     "transverse-momentum--rapidity space of an emitted gluon. This parameter "
     "is used to allow emissions with larger transverse momentum according to "
     "\\f$P(p_\\perp^2>p_{\\perp cut}^2="
     "\\left(\\frac{p_{\\perp cut}^2}{p_\\perp^2}\\right)^\\beta\\f$. if "
     "negative, the sharp cutoff is retained.",
     &CascadeHandler::theBeta, 2.0, -1.0, 0,
     true, false, Interface::lowerlim);

  static RefVector<CascadeHandler,MECorrBase> interfaceMECorrectors
    ("MECorrectors",
     "A vector of objects implementing leading-order matrix-element "
     "corrections of basic dipole emissions. For each initial dipole "
     "in a cascade, each of the object in this vector will be tried in "
     "turn, and the first object to claim to be able to handle the "
     "correction is chosen.",
     &CascadeHandler::theMECorrectors, -1, true, false, true, false, false);


  static RefVector<CascadeHandler,Ariadne::ReweightBase> interfaceReweighters
    ("Reweighters",
     "A vector of objects implementing reweightings of basic dipole "
     "emissions. Each dipole emission will be reweighted.",
     &CascadeHandler::theReweighters, -1, true, false, true, false, false);


  static Parameter<CascadeHandler,int> interfaceMaxEmissions
    ("MaxEmissions",
     "This number specifies the maximum number of emissions from the "
     "cascade. If it is set to zero an unlimited number is allowed. "
     "This parameter should only be used for debugging purposes.",
     &CascadeHandler::theMaxEmissions, 0, 0, 0,
     true, false, Interface::lowerlim);

  static Reference<CascadeHandler,AlphaSBase> interfaceInternalAlphaS
    ("InternalAlphaS",
     "An internal AlphaSBase object to be used if "
     "<interface>RunningCoupling</interface> is set to InternalRunning. "
     "If no such object is given, A O1AlphaS object will be created and "
     "assigned in the initialization with <interface>LambdaQCD</interface> "
     "used as lambda for five flavours.",
     &CascadeHandler::theInternalAlphaS, true, false, true, true, false);

  interfacePTCut.rank(10);
  interfaceLambdaQCD.rank(9);
  interfaceNFlav.rank(8);
  interfaceSoftMu.rank(7);
  interfaceSoftAlpha.rank(6);
  interfaceHardAlpha.rank(5);
  interfaceBeta.rank(4);
  interfaceMECorrectors.rank(2);
}

}
