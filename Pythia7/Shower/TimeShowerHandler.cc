// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TimeShowerHandler class.
//

#include "TimeShowerHandler.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "TimeShowerHandler.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Pythia7;

TimeShowerHandler::~TimeShowerHandler() {
  if ( theShowerModel ) delete theShowerModel;
}

void TimeShowerHandler::persistentOutput(PersistentOStream & os) const {
  os << theAngularOrdering << theNQuark << theAlphaSMode << theMEMode
     << theQEDShower << theInitialCone << thePhiPolAsym << thePhiCoherAsym
     << theRespectScale << ounit(theQ0,GeV) << ounit(theQ0ChgQ,GeV)
     << ounit(theQ0ChgL,GeV) << theAlphaSFix << ounit(theLambda5,GeV)
     << theAlphaEMFix << theQ0FracPS;
}

void TimeShowerHandler::persistentInput(PersistentIStream & is, int) {
  is >> theAngularOrdering >> theNQuark >> theAlphaSMode >> theMEMode
     >> theQEDShower >> theInitialCone >> thePhiPolAsym >> thePhiCoherAsym
     >> theRespectScale >> iunit(theQ0,GeV) >> iunit(theQ0ChgQ,GeV)
     >> iunit(theQ0ChgL,GeV) >> theAlphaSFix >> iunit(theLambda5,GeV)
     >> theAlphaEMFix >> theQ0FracPS;
}

ClassDescription<TimeShowerHandler> TimeShowerHandler::initTimeShowerHandler;
// Definition of the static class description member.

Shower::TimeShower * TimeShowerHandler::getModel() {
  if ( !theShowerModel ) theShowerModel = new Shower::TimeShower;
  theShowerModel->ANGULARORDER = theAngularOrdering;
  theShowerModel->NQUARK = theNQuark;
  theShowerModel->ALPHASMODE = theAlphaSMode;
  theShowerModel->MEMODE = theMEMode;
  theShowerModel->QEDSHOWER = theQEDShower;
  theShowerModel->INITIALCONE = theInitialCone;
  theShowerModel->PHIPOLASYM = thePhiPolAsym;
  theShowerModel->PHICOHERASYM = thePhiCoherAsym;
  theShowerModel->RESPECTSCALE = theRespectScale;
  theShowerModel->Q0 = theQ0/GeV;
  theShowerModel->Q0CHGQ = theQ0ChgQ/GeV;
  theShowerModel->Q0CHGL = theQ0ChgL/GeV;
  theShowerModel->ALPHASFIX = theAlphaSFix;
  theShowerModel->LAMBDA5 = theLambda5/GeV;
  theShowerModel->ALPHAEMFIX = theAlphaEMFix;
  theShowerModel->Q0FRACPS = theQ0FracPS;
  return theShowerModel;
}

void TimeShowerHandler::Init() {

  static ClassDocumentation<TimeShowerHandler> documentation
    ("Administers the time-like shower.");

  static Switch<TimeShowerHandler,int> interfaceAngularOrdering
    ("AngularOrdering",
     "Angular ordering in shower",
     &TimeShowerHandler::theAngularOrdering, 2, true, false);
  static SwitchOption interfaceAngularOrderingOff
    (interfaceAngularOrdering,
     "Off",
     "Don't use angular ordering",
     0);
  static SwitchOption interfaceAngularOrderingGluon
    (interfaceAngularOrdering,
     "Gluon",
     "Use angular ordering for gluon emissions.",
     1);
  static SwitchOption interfaceAngularOrderingAllQCD
    (interfaceAngularOrdering,
     "AllQCD",
     "Use angular ordering for all QCD splittings.",
     2);

  static Parameter<TimeShowerHandler,int> interfaceNQuark
    ("NQuark",
     "Number of allowed quark flavours in \\f$g\\rightarrow q\\bar{q}\\f$ "
     "branching.",
     &TimeShowerHandler::theNQuark, 5, 0, 10,
     true, false, true);

  static Switch<TimeShowerHandler,int> interfaceAlphaSMode
    ("AlphaSMode",
     "Running of \\f$\\alpha_S\\f$ in evolution.",
     &TimeShowerHandler::theAlphaSMode, 2, true, false);
  static SwitchOption interfaceAlphaSModeFix
    (interfaceAlphaSMode,
     "Fix",
     "Use fixed \\f$\\alpha_S\\f$",
     0);
  static SwitchOption interfaceAlphaSModeQ2
    (interfaceAlphaSMode,
     "Q2",
     "Running \\f$\\alpha_S\\f$ with scale \\f$Q^2/4\\f$.",
     1);
  static SwitchOption interfaceAlphaSModePT2
    (interfaceAlphaSMode,
     "PT2",
     "Running \\f$\\alpha_S\\f$ with scale \\f$p_\\perp^2\\f$.",
     2);
  static SwitchOption interfaceAlphaSModeMixed
    (interfaceAlphaSMode,
     "Mixed",
     "Running \\f$\\alpha_S\\f$ with scale \\f$p_\\perp^2\\f$ except for "
     "\\f$g\\rightarrow q\\bar{q}\\f$ where \\f$Q^2/4\\f$ is used.",
     3);

  static Switch<TimeShowerHandler,int> interfaceMEMode
    ("MEMode",
     "Use of matrix element corrections.",
     &TimeShowerHandler::theMEMode, 1, true, false);
  static SwitchOption interfaceMEModeNo
    (interfaceMEMode,
     "No",
     "Don't use matrix element corrections.",
     0);
  static SwitchOption interfaceMEModeYes
    (interfaceMEMode,
     "Yes",
     "Use matrix element corrections.",
     1);

  static Switch<TimeShowerHandler,int> interfaceQEDShower
    ("QEDShower",
     "Allow a QED shower together with QCD ones.",
     &TimeShowerHandler::theQEDShower, 2, true, false);
  static SwitchOption interfaceQEDShowerNo
    (interfaceQEDShower,
     "No",
     "Don't use QED showers",
     0);
  static SwitchOption interfaceQEDShowerYes
    (interfaceQEDShower,
     "Yes",
     "Allow on-shell photon emissions.",
     1);
  static SwitchOption interfaceQEDShowerPhoton
    (interfaceQEDShower,
     "Photon",
     "Allow both photon emissions and photon branchings.",
     2);

  static Switch<TimeShowerHandler,int> interfaceInitialCone
    ("InitialCone",
     "Restrict first emission within cone given by colour flow in "
     "hard process.",
     &TimeShowerHandler::theInitialCone, 2, true, false);
  static SwitchOption interfaceInitialConeNo
    (interfaceInitialCone,
     "No",
     "Don't restric emissions.",
     0);
  static SwitchOption interfaceInitialConeYes
    (interfaceInitialCone,
     "Yes",
     "Restrict emissions with isotropic phi angle inside cone.",
     1);
  static SwitchOption interfaceInitialConeAsym
    (interfaceInitialCone,
     "Asym",
     "Restrict emissions with anisotropic phi angle inside cone.",
     2);

  static Switch<TimeShowerHandler,int> interfacePhiPolAsym
    ("PhiPolAsym",
     "Azimuthal asymmetry induced by gluon polarization.",
     &TimeShowerHandler::thePhiPolAsym, 1, true, false,
     0,
     0, 0);
  static SwitchOption interfacePhiPolAsymNo
    (interfacePhiPolAsym,
     "No",
     "No asymmetry.",
     0);
  static SwitchOption interfacePhiPolAsymYes
    (interfacePhiPolAsym,
     "Yes",
     "Asymmetry",
     1);


  static Switch<TimeShowerHandler,int> interfacePhiCoherAsym
    ("PhiCoherAsym",
     "Azimuthal asymmetry induced by colour coherence.",
     &TimeShowerHandler::thePhiCoherAsym, 1, true, false);
  static SwitchOption interfacePhiCoherAsymNo
    (interfacePhiCoherAsym,
     "No",
     "No asymmetry.",
     0);
  static SwitchOption interfacePhiCoherAsymYes
    (interfacePhiCoherAsym,
     "Yes",
     "Asymmetry",
     1);

  static Switch<TimeShowerHandler,int> interfaceRespectScale
    ("RespectScale",
     "Use the scale variable of original partons to restrict branchings.",
     &TimeShowerHandler::theRespectScale, 0, true, false);
  static SwitchOption interfaceRespectScaleNo
    (interfaceRespectScale,
     "No",
     "No restriction.",
     0);
  static SwitchOption interfaceRespectScaleQ2
    (interfaceRespectScale,
     "Q2",
     "Restrict with \\f$Q^2 <\\f$ scale.",
     1);
  static SwitchOption interfaceRespectScalePT2
    (interfaceRespectScale,
     "PT2",
     "Restrict with \\f$p_\\perp^2 <\\f$ scale.",
     2);
  static SwitchOption interfaceRespectScaleETheta2
    (interfaceRespectScale,
     "ETheta2",
     "Restrict with \\f$(E*\\theta)^2 <\\f$ scale.",
     3);
  static SwitchOption interfaceRespectScaleTheta2
    (interfaceRespectScale,
     "Theta2",
     "Restrict with \\f$\\theta^2 <\\f$ scale.",
     4);

  static Parameter<TimeShowerHandler,Energy> interfaceQ0
    ("Q0",
     "Parton shower cut-off mass for QCD emissions.",
     &TimeShowerHandler::theQ0, GeV, 1.0*GeV, 0.0*GeV, 1000.0*GeV,
     true, false, true);

  static Parameter<TimeShowerHandler,Energy> interfaceQ0ChgQ
    ("Q0ChgQ",
     "Parton shower cut-off mass for photon coupling to coloured particle.",
     &TimeShowerHandler::theQ0ChgQ, GeV, 1.0*GeV, 0.0*GeV, 1000.0*GeV,
     true, false, true);

  static Parameter<TimeShowerHandler,Energy> interfaceQ0ChgL
    ("Q0ChgL",
     "Parton shower cut-off mass for pure QED branchings. Assumed <= Q0ChgQ.",
     &TimeShowerHandler::theQ0ChgL, GeV, 0.001*GeV, 0.0*GeV, 1000.0*GeV,
     true, false, true, 0, 0, 0, &TimeShowerHandler::Q0ChgQ, 0);

  static Parameter<TimeShowerHandler,double> interfaceAlphaSFix
    ("AlphaSFix",
     "Fixed \\f$\\alpha_S\\f$ value for AlphaSMode == 0.",
     &TimeShowerHandler::theAlphaSFix, 0.2, 0.0, 10.0,
     true, false, true);

  static Parameter<TimeShowerHandler,Energy> interfaceLambda5
    ("Lambda5",
     "\\f$\\Lambda_{\\mbox{QCD}}\\f$ (five flavours) in alpha_strong for "
     "AlphaSMode >= 1.",
     &TimeShowerHandler::theLambda5, GeV, 0.25*GeV, 0.0*GeV, 10.0*GeV,
     true, false, true);

  static Parameter<TimeShowerHandler,double> interfaceAlphaEMFix
    ("AlphaEMFix",
     "Fixed \\f$alpha_{\\mbox{EM}}\\f$ value.",
     &TimeShowerHandler::theAlphaEMFix, 0.0073, 0.0, 10.0,
     true, false, true);

  static Parameter<TimeShowerHandler,double> interfaceQ0FracPS
    ("Q0FracPS",
     " Fraction of Q0 cut-off mass used as safety margin in daughter mass "
     "sum. Relevant for total parton multiplicity.",
     &TimeShowerHandler::theQ0FracPS, 0.25, 0.0, 1.0,
     true, false, true);

  interfaceQ0.rank(11);
  interfaceLambda5.rank(10);
  interfaceNQuark.rank(9);
  interfaceAlphaSMode.rank(8);
  interfaceMEMode.rank(7);
  interfaceQEDShower.rank(6);

}

