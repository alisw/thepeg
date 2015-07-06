// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FixedImpactGenerator class.
//

#include "FixedImpactGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Repository/UseRandom.h"

#include "DipoleEventHandler.h" //needed for the selector, where does it point? :o

#include "ThePEG/Utilities/Debug.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "FixedImpactGenerator.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace DIPSY;

FixedImpactGenerator::FixedImpactGenerator()
  : theMinB(0.0*femtometer), theMaxB(0.0*femtometer), theBAngle(-1) {}

FixedImpactGenerator::
FixedImpactGenerator(const FixedImpactGenerator & x)
  : ImpactParameterGenerator(x), theMinB(x.theMinB), theMaxB(x.theMaxB), theBAngle(x.theBAngle) {}

IBPtr FixedImpactGenerator::clone() const {
  return new_ptr(*this);
}

IBPtr FixedImpactGenerator::fullclone() const {
  return new_ptr(*this);
}

FixedImpactGenerator::~FixedImpactGenerator() {}

ImpactParameters FixedImpactGenerator::generate(double seed) const {

  InvEnergy b = sqrt(sqr(minB()) + rnd()*(sqr(maxB()) - sqr(minB())));
  b = max(b, minB());
  double phi = bAngle();
  if ( phi == -1 ) phi = 2.0*Constants::pi*rnd();
  return ImpactParameters(Point(b*cos(phi), b*sin(phi)),
  			  2.0*Constants::pi*rnd(), 1.0/GeV2);
}

void FixedImpactGenerator::persistentOutput(PersistentOStream & os) const {
  os << ounit(theMinB, femtometer) << ounit(theMaxB, femtometer) << theBAngle;
}

void FixedImpactGenerator::persistentInput(PersistentIStream & is, int) {
  is >> iunit(theMinB, femtometer) >> iunit(theMaxB, femtometer) >> theBAngle;
}



// Static variable needed for the type description system in ThePEG.
#include "ThePEG/Utilities/DescribeClass.h"
DescribeClass<FixedImpactGenerator,ImpactParameterGenerator>
  describeDIPSYFixedImpactGenerator("DIPSY::FixedImpactGenerator",
					"FixedImpactGenerator.so");


void FixedImpactGenerator::Init() {

  static ClassDocumentation<FixedImpactGenerator> documentation
    ("The objective of the FixedImpactGenerator class is to generate "
     "ImpactParameters objects to be used when two "
     "<code>DipoleState</code>s have been generated w.r.t.  origo in the "
     "transverse to displace and rotate one of them before "
     "collision. This class will generate the absolute value of the impact "
     "parameter to be within a given interval. Note that the cross sections "
     "produced with this objects will not be trustworthy, but the "
     "distibution of events will.");

  static Parameter<FixedImpactGenerator,Length> interfaceMinB
    ("MinB",
     "The Minimum value of the impact parameter (in fm). ",
     &FixedImpactGenerator::theMinB, femtometer, 0.0*femtometer,
     0.0*femtometer, 0*femtometer,
     true, false, Interface::lowerlim);

  static Parameter<FixedImpactGenerator,Length> interfaceMaxB
    ("MaxB",
     "The Maximum value of the impact parameter (in fm). If set less than "
     "<interface>MinB</interface>, a fixed value of the impact parameter "
     "will be generated according to <interface>MinB</interface>",
     &FixedImpactGenerator::theMaxB, femtometer, 0.0*femtometer,
     0.0*femtometer, 0*femtometer,
     true, false, Interface::lowerlim);

  static Parameter<FixedImpactGenerator,double> interfaceBAngle
    ("BAngle",
     "The angle of the distance B in impact parameter space."
     "Set to -1 for random angle.",
     &FixedImpactGenerator::theBAngle, 0.0, -1.,
     -1., -1.,
     true, false, Interface::lowerlim);

}

