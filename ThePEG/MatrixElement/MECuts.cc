// -*- C++ -*-
//
// MECuts.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MECuts class.
//

#include "MECuts.h"
#include "MECuts.xh"
#include "ThePEG/EventRecord/SubProcess.h"
#include "ThePEG/EventRecord/Collision.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;

MECuts::MECuts()
  : theMHatMin(2.0*GeV), theMHatMax(-1.0*GeV), thePTHatMin(ZERO),
    thePTHatMax(-1.0*GeV), thePTHatSingularMin(1.0*GeV),
    theSingularMassMax(1.0*GeV), theCTHMin(-1.0),
    theCTHMax(1.0), theTHatMin(ZERO), theTHatMax(-1.0*GeV2),
    theUHatMin(ZERO), theUHatMax(-1.0*GeV2),
    theScaleMin(ZERO), theScaleMax(-1.0*GeV2) {}

void MECuts::newcut(const SubProcess &) const {}

void MECuts::cut(const SubProcess & sp) const
{
  newcut(sp);
  const ParticleVector & out = sp.outgoing();
  int N = out.size();
  if ( N < 1 ) throw MECutSetup();
  if ( !sp.incoming().first || !sp.incoming().second ) throw MECutSetup();

  // General cuts:
  if ( sp.shat() < sHatMin() || sp.shat() >= sHatMax() ) throw Veto();

  // Now check 2->2 proceses:
  if ( N != 2 ) return;

  if ( -sp.that() < tHatMin() || -sp.that() >= tHatMax() ) throw Veto();

  if ( -sp.uhat() < uHatMin() || -sp.uhat() >= uHatMax() ) throw Veto();

  double cth = out[0]->momentum().cosTheta();
  if ( cth < cTHMin() || cth >= cTHMax() ) throw Veto();

  Energy pt = out[0]->momentum().perp();
  if ( pt < pTHatMin() || pt >= pTHatMax() ) throw Veto();

}

IBPtr MECuts::clone() const {
  return new_ptr(*this);
}

void MECuts::doupdate() {
  Interfaced::doupdate();
  if ( cTHMax() <= cTHMin() )
    throw MECutZeroInterval(*this, "CTHMax <= CTHMin");
}

void MECuts::persistentOutput(PersistentOStream & os) const {
  os << ounit(theMHatMin, GeV) << ounit(theMHatMax, GeV)
     << ounit(thePTHatMin, GeV) << ounit(thePTHatMax, GeV)
     << ounit(thePTHatSingularMin, GeV) << ounit(theSingularMassMax, GeV)
     << theCTHMin << theCTHMax
     << ounit(theTHatMin, GeV2) << ounit(theTHatMax, GeV2)
     << ounit(theUHatMin, GeV2) << ounit(theUHatMax, GeV2)
     << ounit(theScaleMin, GeV2) << ounit(theScaleMax, GeV2);
}

void MECuts::persistentInput(PersistentIStream & is, int) {
  is >> iunit(theMHatMin, GeV) >> iunit(theMHatMax, GeV)
     >> iunit(thePTHatMin, GeV) >> iunit(thePTHatMax, GeV)
     >> iunit(thePTHatSingularMin, GeV) >> iunit(theSingularMassMax, GeV)
     >> theCTHMin >> theCTHMax >>
    iunit(theTHatMin, GeV2) >> iunit(theTHatMax, GeV2)
     >> iunit(theUHatMin, GeV2) >> iunit(theUHatMax, GeV2)
     >> iunit(theScaleMin, GeV2) >> iunit(theScaleMax, GeV2);
}

ClassDescription<MECuts> MECuts::initMECuts;

void MECuts::Init() {

  static ClassDocumentation<MECuts> documentation
    ("There is no documentation for the ThePEG::MECuts class");

  static Parameter<MECuts,Energy> interfaceMHatMin
    ("SubProcess/MHatMin",
     "The minimum value allowed for \\f$\\hat{m}=\\sqrt{\\hat{s}}\\f$ in GeV "
     "in the hard subprocess. If the upper limit in "
     "<interface>SubProcess/MHatMax</interface> "
     "is less than this, the upper limit is inactive. "
     "This limit is automatically checked by the event handler.",
     &MECuts::theMHatMin, GeV,
     2.0*GeV, 1.0*MeV, Constants::MaxEnergy, false, false, true);

  static Parameter<MECuts,Energy> interfaceMHatMax
    ("SubProcess/MHatMax",
     "The minimum value allowed for \\f$\\hat{m}=\\sqrt{\\hat{s}}\\f$ in GeV "
     "in the hard subprocess If the lower limit in "
     "<interface>SubProcess/MHatMin</interface> "
     "is larger than this, the upper limit is inactive. "
     "This limit is automatically checked by the event handler.",
     &MECuts::theMHatMax, GeV,
     -1.0*GeV, -1.0*GeV, Constants::MaxEnergy, false, false, true);

  static Parameter<MECuts,Energy> interfacePTHatMin
    ("SubProcess/PTHatMin",
     "The minimum value allowed for \\f$\\hat{p_\\perp}\\f$ in GeV "
     "in the hard subprocess. If the upper limit in "
     "<interface>SubProcess/PTHatMax</interface> "
     "is less than this, the upper limit is inactive. "
     "This limit is automatically checked by the event handler",
     &MECuts::thePTHatMin, GeV,
     ZERO, ZERO, Constants::MaxEnergy, false, false, true);

  static Parameter<MECuts,Energy> interfacePTHatMax
    ("SubProcess/PTHatMax",
     "The minimum value allowed for \\f$\\hat{p_\\perp}\\f$ in GeV "
     "in the hard subprocess If the lower limit in "
     "<interface>SubProcess/PTHatMin</interface> "
     "is larger than this, the upper limit is inactive. "
     "This limit is automatically checked by the event handler.",
     &MECuts::thePTHatMax, GeV,
     -1.0*GeV, -1.0*GeV, Constants::MaxEnergy, false, false, true);

  static Parameter<MECuts,Energy> interfacePTHatSingularMin
    ("SubProcess/PTHatSingularMin",
     "The minimum value allowed for \\f$\\hat{p_\\perp}\\f$ in GeV for "
     "processes which are singular in the limit "
     "\\f$\\hat{p_\\perp}\\rightarrow 0\\f$. This "
     "cut is in addition to PTHatMin. Hard \\f$2\\rightarrow 2\\f$ "
     "processes which do not proceed via intermediate resonances are "
     "considered singular if either or both final-state products have a mass "
     "less than SingularMassMax. This limit is not checked "
     "automatically, but is assumed to be checked by the relevant "
     "ThePEG::PartonXSecFn objects.",
     &MECuts::thePTHatSingularMin, GeV,
     1.0*GeV, ZERO, Constants::MaxEnergy,
     false, false, true);


  static Parameter<MECuts,Energy> interfaceSingularMassMax
    ("SubProcess/SingularMassMax",
     "Hard \\f$2\\rightarrow 2\\f$ processes which do not proceed via "
     "intermediate resonances are considered singular if either or both "
     "final-state products have a mass less than this limit (int GeV). "
     "For singular processes the aditional \\f$\\hat{p_\\perp}\\f$ cut in "
     "PTHatSingularMin should be applied. This limit is not "
     "checked automatically, but is assumed to be checked by the relevant "
     "ThePEG::PartonXSecFn objects.",
     &MECuts::theSingularMassMax, GeV,
     1.0*GeV, ZERO, Constants::MaxEnergy,
     false, false, true);

  static Parameter<MECuts,double> interfaceCTHMin
    ("SubProcess/CosThetaHatMin",
     "The minimum allowed value of \\f$\\cos{\\hat{\\theta}}\\f$, where "
     "\\f$\\hat{\\theta}\\f$ is the scattering angle in the restframe of a "
     "hard \\f$2\\rightarrow 2\\f$ scattering. "
     "This limit is automatically checked by the event handler.",
     &MECuts::theCTHMin,
     -1.0, -1.0, 1.0, false, false, true,
     0, 0, 0, &MECuts::cTHMax);

  static Parameter<MECuts,double> interfaceCTHMax
    ("SubProcess/CosThetaHatMax",
     "The maximum allowed value of \\f$\\cos{\\hat{\\theta}}\\f$, where "
     "\\f$\\hat{\\theta}\\f$ is the scattering angle in the restframe of a "
     "hard \\f$2\\rightarrow 2\\f$ scattering. "
     "This limit is automatically checked by the event handler.",
     &MECuts::theCTHMax,
     1.0, -1.0, 1.0, false, false, true,
     0, 0, &MECuts::cTHMin, 0);

  static Parameter<MECuts,Energy2> interfaceTHatMin
    ("SubProcess/THatMin",
     "The minimum allowed value of \\f$|\\hat{t}|=-\\hat{t}\\f$ in "
     "GeV<sup>2</sup> in a hard \\f$2\\rightarrow 2\\f$ scattering. If the "
     "upper limit in <interface>SubProcess/THatMax</interface> "
     "is less than this, the upper limit is inactive. "
     "This limit is automatically checked by the event handler.",
     &MECuts::theTHatMin, GeV2,
     ZERO, ZERO, Constants::MaxEnergy2, false, false, true);

  static Parameter<MECuts,Energy2> interfaceTHatMax
    ("SubProcess/THatMax",
     "The maximum allowed value of \\f$|\\hat{t}|=-\\hat{t}\\f$ in "
     "GeV<sup>2</sup> in a hard \\f$2\\rightarrow 2\\f$ scattering. If the "
     "lower limit in <interface>SubProcess/THatMin</interface> "
     "is larger than this, the upper limit is inactive. "
     "This limit is automatically checked by the event handler.",
     &MECuts::theTHatMax, GeV2,
     -1.0*GeV, -1.0*GeV, Constants::MaxEnergy2, false, false, true);

  static Parameter<MECuts,Energy2> interfaceUHatMin
    ("SubProcess/UHatMin",
     "The minimum allowed value of \\f$|\\hat{u}|=-\\hat{u}\\f$ in "
     "GeV<sup>2</sup> in a hard \\f$2\\rightarrow 2\\f$ scattering. If the "
     "upper limit in <interface>SubProcess/UHatMax</interface> "
     "is less than this, the upper limit is inactive. "
     "This limit is automatically checked by the event handler.",
     &MECuts::theUHatMin, GeV2,
     ZERO, ZERO, Constants::MaxEnergy2, false, false, true);

  static Parameter<MECuts,Energy2> interfaceUHatMax
    ("SubProcess/UHatMax",
     "The maximum allowed value of \\f$|\\hat{u}|=-\\hat{u}\\f$ in "
     "GeV<sup>2</sup> in a hard \\f$2\\rightarrow 2\\f$ scattering. If the "
     "lower limit in <interface>SubProcess/UHatMin</interface> "
     "is larger than this, the upper limit is inactive. "
     "This limit is automatically checked by the event handler.",
     &MECuts::theUHatMax, GeV2,
     -1.0*GeV, -1.0*GeV, Constants::MaxEnergy2, false, false, true);

  static Parameter<MECuts,Energy2> interfaceScaleMin
    ("SubProcess/ScaleMin",
     "The minimum allowed value of the user-defined scale in GeV<sup>2</sup> "
     "in a hard scattering. If the upper limit in "
     "<interface>SubProcess/ScaleMax</interface> "
     "is less than this, the upper limit is "
     "inactive. This limit is automatically checked by the event handler.",
     &MECuts::theScaleMin, GeV2,
     ZERO, ZERO, Constants::MaxEnergy2, false, false, true);

  static Parameter<MECuts,Energy2> interfaceScaleMax
    ("SubProcess/ScaleMax",
     "The maximum allowed value of user defined scale in GeV<sup>2</sup> in a "
     "hard scattering. If the lower limit in "
     "<interface>SubProcess/ScaleMin</interface> "
     "is larger than this, the upper limit is inactive. "
     "This limit is automatically checked by the event handler.",
     &MECuts::theScaleMax, GeV2,
     -1.0*GeV, -1.0*GeV, Constants::MaxEnergy2, false, false, true);

}

MECutSetup::MECutSetup() {
  theMessage << "MECuts was asked to check an icomplete "
	     << "SubProcess/Collision/Event. This should never happen."
	     << " the run will be aborted immediately.";
  severity(abortnow);
}

MECutZeroInterval::MECutZeroInterval(const MECuts & i, string s) {
  theMessage << "The MECuts object '" << i.name() << "' is not "
	     << "Properly set up since an interval was non-existent :"
	     << s << ".";
  severity(warning);
}

