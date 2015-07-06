// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GaussianImpactGenerator class.
//

#include "GaussianImpactGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Repository/UseRandom.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "GaussianImpactGenerator.tcc"
#endif

using namespace DIPSY;

GaussianImpactGenerator::GaussianImpactGenerator() {}

IBPtr GaussianImpactGenerator::clone() const {
  return new_ptr(*this);
}

IBPtr GaussianImpactGenerator::fullclone() const {
  return new_ptr(*this);
}

GaussianImpactGenerator::~GaussianImpactGenerator() {}

ImpactParameters GaussianImpactGenerator::generate(double rnd0) const {

  InvEnergy b = sqrt(-2.0*log(rnd0))*width();
  double phi = 2.0*Constants::pi*rnd();
  InvEnergy2 w = 2.0*Constants::pi*sqr(width())*exp(sqr(b/width())/2.0);
  return ImpactParameters(ImpactParameters::Point(b*cos(phi), b*sin(phi)),
   			  2.0*Constants::pi*rnd(), w);

}

// Static variable needed for the type description system in ThePEG.
#include "ThePEG/Utilities/DescribeClass.h"
DescribeNoPIOClass<GaussianImpactGenerator,ImpactParameterGenerator>
  describeDIPSYGaussianImpactGenerator("DIPSY::GaussianImpactGenerator",
					"GaussianImpactGenerator.so");


void GaussianImpactGenerator::Init() {

  static ClassDocumentation<GaussianImpactGenerator> documentation
    ("The objective of the GaussianImpactGenerator class is to generate "
     "ImpactParameters objects to be used when two "
     "<code>DipoleState</code>s have been generated w.r.t.  origo in the "
     "transverse to displace and rotate one of them before "
     "collision. This base class will generate impact parameters "
     "according to a Gaussian distribution, and the weigth in the "
     "produced ImpactParameters objects is set accordingly. Sub-classes "
     "may override the generate() function to use any distribution as "
     "long as the weight is set accordingly.");

}

