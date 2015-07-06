// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SpaceShowerHandler class.
//

#include "SpaceShowerHandler.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "SpaceShowerHandler.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Pythia7;

SpaceShowerHandler::~SpaceShowerHandler() {
  if ( theShowerModel ) delete theShowerModel;
}

Shower::SpaceShower * SpaceShowerHandler::getModel() {
  if ( !theShowerModel ) theShowerModel = new Shower::SpaceShower;
  theShowerModel->HADRONSHOWER = theHadronShower;
  theShowerModel->LEPTONSHOWER = theLeptonShower;
  theShowerModel->NQUARK = theNQuarks;
  theShowerModel->ALPHASMODE = theAlphaSMode;
  theShowerModel->MAXVIRTUALITY = theMaxVirtuality;
  theShowerModel->MEMODE = theMEMode;
  theShowerModel->SOFTGLUONRESUM = theSoftGluonResum;
  theShowerModel->FINALCONE = theFinalCone;
  theShowerModel->Q2ORDER = theQ2Order;
  theShowerModel->ANGULARORDER = useAngularOrdering;
  theShowerModel->PHIPOLASYM = thePhiPolAsym;
  theShowerModel->PHICOHERASYM = thePhiCoherAsym;
  theShowerModel->RESPECTSCALE = theRespectScale;
  theShowerModel->Q0 = theQ0/GeV;
  theShowerModel->Q0CHGQ = theQ0ChgQ/GeV;
  theShowerModel->Q0CHGL = theQ0ChgL/GeV;
  theShowerModel->ALPHASFIX = theAlphaSFix;
  theShowerModel->LAMBDA5 = theLambda5/GeV;
  theShowerModel->ALPHAEMFIX = theAlphaEMFix;
  theShowerModel->EMINEMITTED = theEMinEmitted/GeV;
  theShowerModel->ZMINEMITTED = theZMinEmitted;
  theShowerModel->XMINEMITTEDCHG = theXMinEmittedChg;
  theShowerModel->TINYQCHG = theTinyQChg/GeV;
  theShowerModel->TINYPDF = theTinyPDF;
  theShowerModel->TINYKERNELPDF = theTinyKernelPDF;
  theShowerModel->TINYKINPREC = theTinyKinPrec;
  theShowerModel->HEAVYEVOL = theHeavyEvol;
  theShowerModel->EXTRAQEDPREWT = theExtraPreweight;
  theShowerModel->HEAVYXMAX = theHeavyMax;
  theShowerModel->Q2STARTFRAC = theQ2StartFrac;
  return theShowerModel;
}

void SpaceShowerHandler::persistentOutput(PersistentOStream & os) const {
  os << theHadronShower << theLeptonShower << theNQuarks << theAlphaSMode
     << theMaxVirtuality << theMEMode << theSoftGluonResum << theFinalCone
     << theQ2Order << useAngularOrdering << thePhiPolAsym << thePhiCoherAsym
     << theRespectScale << ounit(theQ0, GeV) << ounit(theQ0ChgQ, GeV)
     << ounit(theQ0ChgL, GeV) << theAlphaSFix << ounit(theLambda5, GeV)
     << theAlphaEMFix << ounit(theEMinEmitted, GeV) << theZMinEmitted
     << theXMinEmittedChg << ounit(theTinyQChg, GeV) << theTinyPDF
     << theTinyKernelPDF << theTinyKinPrec << theHeavyEvol
     << theExtraPreweight << theHeavyMax << theQ2StartFrac;
}

void SpaceShowerHandler::persistentInput(PersistentIStream & is, int) {
  is >> theHadronShower >> theLeptonShower >> theNQuarks >> theAlphaSMode
     >> theMaxVirtuality >> theMEMode >> theSoftGluonResum >> theFinalCone
     >> theQ2Order >> useAngularOrdering >> thePhiPolAsym >> thePhiCoherAsym
     >> theRespectScale >> iunit(theQ0, GeV) >> iunit(theQ0ChgQ, GeV)
     >> iunit(theQ0ChgL, GeV) >> theAlphaSFix >> iunit(theLambda5, GeV)
     >> theAlphaEMFix >> iunit(theEMinEmitted, GeV) >> theZMinEmitted
     >> theXMinEmittedChg >> iunit(theTinyQChg, GeV) >> theTinyPDF
     >> theTinyKernelPDF >> theTinyKinPrec >> theHeavyEvol
     >> theExtraPreweight >> theHeavyMax >> theQ2StartFrac;
}

ClassDescription<SpaceShowerHandler> SpaceShowerHandler::initSpaceShowerHandler;
// Definition of the static class description member.

void SpaceShowerHandler::Init() {

  static ClassDocumentation<SpaceShowerHandler> documentation
    ("Administers the space-like shower.");


  static Switch<SpaceShowerHandler,int> interfaceHadronShower
    ("HadronShower",
     "Which kind of ISR shower is allowed for incoming hadrons.",
     &SpaceShowerHandler::theHadronShower, 2, true, false);
  static SwitchOption interfaceHadronShowerNone
    (interfaceHadronShower,
     "None",
     "No ISR allowed.",
     0);
  static SwitchOption interfaceHadronShowerQCD
    (interfaceHadronShower,
     "QCD",
     "Only QCD emissions.",
     1);
  static SwitchOption interfaceHadronShowerPhoton
    (interfaceHadronShower,
     "Photon",
     "QCD and photon emissions.",
     2);
  static SwitchOption interfaceHadronShowerSpacePhoton
    (interfaceHadronShower,
     "SpacePhoton",
     "QCD, photon emissions and space-like photon branchings "
     "(if supported by the PDF).",
     3);

  static Switch<SpaceShowerHandler,int> interfaceLeptonShower
    ("LeptonShower",
     "Which kind of ISR shower is allowed for incoming leptons.",
     &SpaceShowerHandler::theLeptonShower, 1, true, false);
  static SwitchOption interfaceLeptonShowerNone
    (interfaceLeptonShower,
     "None",
     "No ISR allowed.",
     0);
  static SwitchOption interfaceLeptonShowerPhoton
    (interfaceLeptonShower,
     "Photon",
     "Only photon emission.",
     1);
  static SwitchOption interfaceLeptonShowerSpacePhoton
    (interfaceLeptonShower,
     "SpacePhoton",
     "Photon emission and space-like photon branchings "
     "(if supported by the PDF).",
     2);

  static Parameter<SpaceShowerHandler,int> interfaceNQuarks
    ("NQuarks",
     "Number of allowed quark flavours in \\f$g\\rightarrow q\\bar{q}\\f$ "
     "branching.",
     &SpaceShowerHandler::theNQuarks, 5, 0, 8,
     true, false, true);

  static Switch<SpaceShowerHandler,int> interfaceAlphaSMode
    ("AlphaSMode",
     "Running of alpha_S in evolution.",
     &SpaceShowerHandler::theAlphaSMode, 2, true, false);
  static SwitchOption interfaceAlphaSModeFixed
    (interfaceAlphaSMode,
     "Fixed",
     "Use fixed \\f$\\alpha_S\\f$.",
     0);
  static SwitchOption interfaceAlphaSModeQ2
    (interfaceAlphaSMode,
     "Q2",
     "Running \\f$\\alpha_S\\f$ with scale \\f$Q^2\\f$.",
     1);
  static SwitchOption interfaceAlphaSModePT2
    (interfaceAlphaSMode,
     "PT2",
     "Running \\f$\\alpha_S\\f$ with scale \\f$p_\\perp^2\\f$.",
     2);

  static Switch<SpaceShowerHandler,int> interfaceMaxVirtuality
    ("MaxVirtuality",
     "Maximum virtuality setting the starting point of the evolution.",
     &SpaceShowerHandler::theMaxVirtuality, 0, true, false);
  static SwitchOption interfaceMaxVirtualityS
    (interfaceMaxVirtuality,
     "S",
     "Total invariant mass squared.",
     0);
  static SwitchOption interfaceMaxVirtualitySHat
    (interfaceMaxVirtuality,
     "SHat",
     "Invariant mass squared of hard sub-process.",
     1);
  static SwitchOption interfaceMaxVirtualityMT2
    (interfaceMaxVirtuality,
     "MT2",
     "Average \\f$m_\\perp^2\\f$ of hard sub-process.",
     2);
  static SwitchOption interfaceMaxVirtualityMinMT2
    (interfaceMaxVirtuality,
     "MinMT2",
     "Smallest \\f$m_\\perp^2\\f$ in hard sub-process.",
     3);

  static Switch<SpaceShowerHandler,int> interfaceMEMode
    ("MEMode",
     "Use of matrix element corrections.",
     &SpaceShowerHandler::theMEMode, 1, true, false);
  static SwitchOption interfaceMEModeYes
    (interfaceMEMode,
     "Yes",
     "Use matrix element corrections.",
     1);
  static SwitchOption interfaceMEModeNo
    (interfaceMEMode,
     "No",
     "Don't use matrix element corrections.",
     0);

  static Switch<SpaceShowerHandler,int> interfaceSoftGluonResum
    ("SoftGluonResum",
     "Resum the effect of multiple soft gluon emissions.",
     &SpaceShowerHandler::theSoftGluonResum, 1, true, false);
  static SwitchOption interfaceSoftGluonResumYes
    (interfaceSoftGluonResum,
     "Yes",
     "Use resummation.",
     1);
  static SwitchOption interfaceSoftGluonResumNo
    (interfaceSoftGluonResum,
     "No",
     "Don't use resummation.",
     0);

  static Switch<SpaceShowerHandler,int> interfaceFinalCone
    ("FinalCone",
     "Restrict first emission within cone given by colour flow in hard process.",
     &SpaceShowerHandler::theFinalCone, 2, true, false);
  static SwitchOption interfaceFinalConeNo
    (interfaceFinalCone,
     "No",
     "No restriction.",
     0);
  static SwitchOption interfaceFinalConeYes
    (interfaceFinalCone,
     "Yes",
     "Use cone with isotropic phi-angle inside.",
     1);
  static SwitchOption interfaceFinalConeAnisotropic
    (interfaceFinalCone,
     "Anisotropic",
     "Use cone with anisotropic phi-angle inside.",
     2);

  static Switch<SpaceShowerHandler,int> interfaceQ2Order
    ("Q2Order",
     "Usage of strict Q2 ordering.",
     &SpaceShowerHandler::theQ2Order, 0, true, false);
  static SwitchOption interfaceQ2OrderYes
    (interfaceQ2Order,
     "Yes",
     "Use \\f$Q^2\\f$ ordering",
     0);
  static SwitchOption interfaceQ2OrderNo
    (interfaceQ2Order,
     "No",
     "No \\f$Q^2\\f$ ordering, use maximum kinematically allowed  \\f$Q^2\\f$ "
     "in each branching (toy scenario to be developed).",
     1);

  static Switch<SpaceShowerHandler,int> interfaceAngularOrdering
    ("AngularOrdering",
     "Determines whether the space-like shower should use angular ordering "
     "or not.",
     &SpaceShowerHandler::useAngularOrdering, 1, true, false);
  static SwitchOption interfaceAngularOrderingYes
    (interfaceAngularOrdering,
     "Yes",
     "Use angular ordering.",
     1);
  static SwitchOption interfaceAngularOrderingNo
    (interfaceAngularOrdering,
     "No",
     "Don't use angular ordering.",
     0);

  static Switch<SpaceShowerHandler,int> interfacePhiPolAsym
    ("PhiPolAsym",
     "Azimuthal asymmetry induced by gluon polarization.",
     &SpaceShowerHandler::thePhiPolAsym, 1, true, false);
  static SwitchOption interfacePhiPolAsymYes
    (interfacePhiPolAsym,
     "Yes",
     "Use Azimuthal asymmetry.",
     1);
  static SwitchOption interfacePhiPolAsymNo
    (interfacePhiPolAsym,
     "No",
     "Don't use Azimuthal asymmetry.",
     0);

  static Switch<SpaceShowerHandler,int> interfacePhiCoherAsym
    ("PhiCoherAsym",
     "Azimuthal asymmetry induced by colour coherence.",
     &SpaceShowerHandler::thePhiCoherAsym, 1, true, false);
  static SwitchOption interfacePhiCoherAsymYes
    (interfacePhiCoherAsym,
     "Yes",
     "Use Azimuthal asymmetry.",
     1);
  static SwitchOption interfacePhiCoherAsymNo
    (interfacePhiCoherAsym,
     "No",
     "Don't use Azimuthal asymmetry.",
     0);

  static Switch<SpaceShowerHandler,int> interfaceRespectScale
    ("RespectScale",
     "Use the scale variable of original partons to restrict branchings.",
     &SpaceShowerHandler::theRespectScale, 0, true, false);
  static SwitchOption interfaceRespectScaleNo
    (interfaceRespectScale,
     "No",
     "Don't restrict.",
     0);
  static SwitchOption interfaceRespectScaleQ2
    (interfaceRespectScale,
     "Q2",
     "Restrict \\f$Q^2 <\\f$ scale of original parton.",
     1);
  static SwitchOption interfaceRespectScalePT2
    (interfaceRespectScale,
     "PT2",
     "Restrict \\f$p_\\perp^2 <\\f$ scale of original parton.",
     2);
  static SwitchOption interfaceRespectScaleETheta2
    (interfaceRespectScale,
     "ETheta2",
     "Restrict \\f$(E*\\theta)^2 <\\f$ scale of original parton.",
     3);
  static SwitchOption interfaceRespectScaleTheta2
    (interfaceRespectScale,
     "Theta2",
     "Restrict \\f$\\theta^2 <\\f$ scale of original parton.",
     4);

  static Parameter<SpaceShowerHandler,Energy> interfaceQ0
    ("Q0",
     "Parton shower cut-off mass for QCD emissions.",
     &SpaceShowerHandler::theQ0, GeV, 1.0*GeV, 0.0*GeV, 1.0e10*GeV,
     true, false, true);

  static Parameter<SpaceShowerHandler,Energy> interfaceQ0ChgQ
    ("Q0ChgQ",
     "Parton shower cut-off mass for photon coupling to coloured particle.",
     &SpaceShowerHandler::theQ0ChgQ, GeV, 1.0*GeV, 0.0*GeV, 1.0e10*GeV,
     true, false, true);

  static Parameter<SpaceShowerHandler,Energy> interfaceQ0ChgL
    ("Q0ChgL",
     "Parton shower cut-off mass for pure QED branchings. Assumed <= Q0ChgQ.",
     &SpaceShowerHandler::theQ0ChgL, GeV, 0.001*GeV, 0.0*GeV, 1.0e10*GeV,
     true, false, true, 0, 0, 0, 0, &SpaceShowerHandler::Q0ChgQ);

  static Parameter<SpaceShowerHandler,double> interfaceAlphaSFix
    ("AlphaSFix",
     "Fixed alpha_strong value for AlphaSMode == 0.",
     &SpaceShowerHandler::theAlphaSFix, 0.2, 0.0, 10.0,
     true, false, true);

  static Parameter<SpaceShowerHandler,Energy> interfaceLambda5
    ("Lambda5",
     "\\f$\\Lambda_{\\mbox{QCD}}\\f$ (five flavours) in alpha_strong for "
     "AlphaSMode >= 1.",
     &SpaceShowerHandler::theLambda5, GeV, 0.18*GeV, 0.0*GeV, 10.0*GeV,
     true, false, true);

  static Parameter<SpaceShowerHandler,double> interfaceAlphaEMFix
    ("AlphaEMFix",
     "Fixed \\f$\\alpha_{\\mbox{EM}}\\f$ value. ",
     &SpaceShowerHandler::theAlphaEMFix, 0.0073, 0.0, 10.0,
     true, false, true);

  static Parameter<SpaceShowerHandler,Energy> interfaceEMinEmitted
    ("EMinEmitted",
     "Minimum energy of emitted QCD parton in rest frame of subprocess.",
     &SpaceShowerHandler::theEMinEmitted, GeV, 2.0*GeV, 0.0*GeV, 1.0e10*GeV,
     true, false, true);

  static Parameter<SpaceShowerHandler,double> interfaceZMinEmitted
    ("ZMinEmitted",
     " Minimum fraction 1 - z of emitted QCD parton, in addition to other limits.",
     &SpaceShowerHandler::theZMinEmitted, 0.001, 0.0, 1.0,
     true, false, true);

  static Parameter<SpaceShowerHandler,double> interfaceXMinEmittedChg
    ("XMinEmittedChg",
     "Minimum x fraction of emitted photon - matched to treatment of photon PDF.",
     &SpaceShowerHandler::theXMinEmittedChg, 1.0e-10, 0.0, 1.0,
     true, false, true);

  static Parameter<SpaceShowerHandler,Energy> interfaceTinyQChg
    ("TinyQChg",
     "Smallest particle mass for QED evolution (= electron mass).",
     &SpaceShowerHandler::theTinyQChg, GeV, 0.00051*GeV, 0.0*GeV, 1.0*GeV,
     true, true, true);

  static Parameter<SpaceShowerHandler,double> interfaceTinyPDF
    ("TinyPDF",
     "Vanishingly small parton density.",
     &SpaceShowerHandler::theTinyPDF, 1.0e-10, 0.0, 1.0,
     true, true, true);

  static Parameter<SpaceShowerHandler,double> interfaceTinyKernelPDF
    ("TinyKernelPDF",
     "Vanishingly small product of splitting kernels and "
     "parton density ratios.",
     &SpaceShowerHandler::theTinyKernelPDF, 1.0e-5, 0.0, 1.0,
     true, true, true);

  static Parameter<SpaceShowerHandler,double> interfaceTinyKinPrec
    ("TinyKinPrec",
     "Vanishingly small recoil mass in branching kinematics reconstruction.",
     &SpaceShowerHandler::theTinyKinPrec, 1.0e-10, 0.0, 1.0,
     true, true, true);

  static Parameter<SpaceShowerHandler,double> interfaceHeavyEvol
    ("HeavyEvol",
     "Safety margin in x that heavy flavour evolution is at all possible.",
     &SpaceShowerHandler::theHeavyEvol, 0.8, 0.0, 1.0,
     true, true, true);

  static Parameter<SpaceShowerHandler,double> interfaceExtraPreweight
    ("ExtraPreweight",
     "Extra preweight in QED shower evolution, to avoid maximum violation.",
     &SpaceShowerHandler::theExtraPreweight, 0.1, 0.0, 1.0,
     true, true, true);

  static Parameter<SpaceShowerHandler,double> interfaceHeavyMax
    ("HeavyMax",
     "Maximum allowed x when reconstructing back to heavy flavour "
     "from gluon or photon.",
     &SpaceShowerHandler::theHeavyMax, 0.5, 0.0, 1.0,
     true, true, true);

  static Parameter<SpaceShowerHandler,double> interfaceQ2StartFrac
    ("Q2StartFrac",
     "Mimimum gap in \\f$Q^2\\f$ values to allow iteration when parton "
     "density vanishes.",
     &SpaceShowerHandler::theQ2StartFrac, 0.9, 0.0, 1.0,
     true, true, true);

  interfaceHadronShower.rank(10);
  interfaceLeptonShower.rank(9);
  interfaceQ0.rank(8);
  interfaceLambda5.rank(7);
  interfaceNQuarks.rank(6);
  interfaceAlphaSMode.rank(5);
  interfaceMaxVirtuality.rank(4);
  interfaceMEMode.rank(3);

}

