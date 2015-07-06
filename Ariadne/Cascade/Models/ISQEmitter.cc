// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ISQEmitter class.
//

#include "ISQEmitter.h"
#include "ISQEmission.h"
#include "RemnantModel.h"
#include "Ariadne/Cascade/QCDDipole.h"
#include "Ariadne/Cascade/DipoleState.h"
#include "Ariadne/Cascade/AriadneHandler.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Utilities/UtilityBase.h"
#include "ThePEG/Utilities/Current.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Utilities/DebugItem.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Repository/CurrentGenerator.h"


#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Ariadne5;

ISQEmitter::ISQEmitter(): maxPDFX(0.7) {}

ISQEmitter::~ISQEmitter() {}

bool ISQEmitter::overrides(const EmitterBase &, DipoleBase &) const {
  return false;
}

ISQEmitter::PDFLimits
ISQEmitter::
maxPDFRatios(tRemParPtr rem, Energy2 mt2max, Energy2 mt2min, tcPDPtr g,
	     Energy2 s, Energy2 mh2,  Energy2 mq2, Energy2 Q2) const {
  static DebugItem debuglimits("Ariadne5::ISQEmitter::PDFLimits");
  static DebugItem debuglimits0("Ariadne5::ISQEmitter::PDFLimits0");
  static DebugItem debuglimits2("Ariadne5::ISQEmitter::PDFLimits2");
  ISQEmitter::PDFLimits limits;
  double x = rem->x();

  // First we devide in scale intevals so that the denominator does
  // not vary too much in each interval. The minimum interval is given
  // by the cutoff.
  map<Energy2,double> denoms;
  Energy2 step = sqr(handler().pTCut());
  denoms[mt2max] = rem->xfx(mt2max);
  denoms[mt2min] = rem->xfx(mt2min);
  if ( rem->x() > 0.5 )
    denoms[sqrt(mt2max*mt2min)] = rem->xfx(sqrt(mt2max*mt2min));
  bool done = false;
  const double rmax = 4.0;
  while ( !done ) {
    done = true;
    for ( map<Energy2,double>::iterator it1 = denoms.begin();
	  it1 != denoms.end(); ) {
      map<Energy2,double>::iterator it2 = it1++;
      if ( it1 == denoms.end() ) break;
      if ( ( it1->second > it2->second*rmax || it2->second > it1->second*rmax )
	   && it1->first > it2->first + step ) {
	Energy2 mt2 = sqrt(it1->first*it2->first);
	denoms[mt2] = rem->xfx(mt2);
	done = false;
	break;
      }
    }
  }

  // Now for each interval we get the estimate of the largest ratio.
  for ( map<Energy2,double>::reverse_iterator it1 = denoms.rbegin();
	it1 != denoms.rend(); ) {
    map<Energy2,double>::reverse_iterator it2 = it1++;
    if ( it1 == denoms.rend() || it1->second <= 0.0 ) break;
    double z1 = RemnantModel::zmax(it1->first, mq2, mh2, Q2);
    double r1 = rem->xfx(g, it1->first, x/z1)/it1->second;
    double z2 = RemnantModel::zmax(it2->first, mq2, mh2, Q2);
    double r2 = rem->xfx(g, it2->first, x/z2)/it2->second;
    limits.push_back(make_pair(it1->first, 2.0*max(r1, r2)));
  }

  if ( debuglimits ||
       ( debuglimits2 && limits.size() > 1 ) ||
       ( debuglimits0 && limits.size() < 1 ) ) {
    cerr << "Ariadne5::ISQEmitter::PDFLimits: "
	 << g->id() << "/" << rem->data().id() << " x = " << x
	 << " minmax = " << sqrt(mt2max)/GeV << "/" << sqrt(mt2min)/GeV << endl;
    for ( int i = 0, N = limits.size(); i < N; ++i )
      cerr << "    " << sqrt(limits[i].first)/GeV
	   << ": " << limits[i].second << endl;
  }

  reverse(limits.begin(), limits.end());
  return limits;
}

void ISQEmitter::
setMom(QCDDipole & d, const ISQEmission & e, tRemParPtr rem, tcPDPtr qex) const {
  // First set the new extracted parton and other kinematics.
  e.Rh = rem->getHardTransform(e.pold.second, e.genmom.third);
  d.state()->transformHadronicState(e.Rh);
  rem->setMomentum(e.genmom.first, qex);
}

void ISQEmitter::
revertMom(QCDDipole & d, const ISQEmission & e, tRemParPtr rem) const {
  rem->setMomentum(e.pold.first, e.exorig);
  d.state()->transformHadronicState(e.Rh.inverse());
  rem->x(e.x);

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void ISQEmitter::persistentOutput(PersistentOStream & os) const {
  os << maxPDFX;
}

void ISQEmitter::persistentInput(PersistentIStream & is, int) {
  is >> maxPDFX;
}

DescribeAbstractClass<ISQEmitter,EmitterBase>
describeAriadne5ISQEmitter("Ariadne5::ISQEmitter", "libAriadne5.so");

void ISQEmitter::Init() {

  static ClassDocumentation<ISQEmitter> documentation
    ("The ISQEmitter class implements functionality for initial-state "
     "quark emissions.");

  static Parameter<ISQEmitter,double> interfaceMaxPDFXRatio
    ("MaxPDFXRatio",
     "For large x-values, some PDFs suffer from lack of precision and "
     "sometimes negative values. Above the x-value given here, if the "
     "PDF ratio is larger than the estimated maximum, it is simply set "
     "to the maximum ratio.",
     &ISQEmitter::maxPDFX, 0.7, 0.0, 1.0,
     true, false, Interface::limited);

}

