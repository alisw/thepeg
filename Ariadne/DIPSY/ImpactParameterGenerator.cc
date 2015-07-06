// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ImpactParameterGenerator class.
//

#include "ImpactParameterGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Repository/UseRandom.h"

#include "DipoleEventHandler.h" //needed for the selector, where does it point? :o

#include "ThePEG/Utilities/Debug.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "ImpactParameterGenerator.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace DIPSY;

ImpactParameterGenerator::ImpactParameterGenerator()
  : theWidth(5.0*InvGeV) {}

ImpactParameterGenerator::
ImpactParameterGenerator(const ImpactParameterGenerator & x)
  : HandlerBase(x), theWidth(x.theWidth) {}

IBPtr ImpactParameterGenerator::clone() const {
  return new_ptr(*this);
}

IBPtr ImpactParameterGenerator::fullclone() const {
  return new_ptr(*this);
}

ImpactParameterGenerator::~ImpactParameterGenerator() {}

ImpactParameters ImpactParameterGenerator::generate(double seed) const {

  // InvEnergy b = sqrt(-2.0*log(rnd()))*width();
  // double phi = 2.0*Constants::pi*rnd();
  // InvEnergy2 w = 2.0*Constants::pi*sqr(width())*exp(sqr(b/width())/2.0);
  // return ImpactParameters(ImpactParameters::Point(b*cos(phi), b*sin(phi)),
  // 			  2.0*Constants::pi*rnd(), w);

  InvEnergy b = -log(seed)*width();
  double phi = 2.0*Constants::pi*rnd();
  InvEnergy2 w = 2.0*Constants::pi*width()*b*exp(b/width());
  return ImpactParameters(Point(b*cos(phi), b*sin(phi)),
  			  2.0*Constants::pi*rnd(), w);
}

ImpactParameters ImpactParameterGenerator::
generateDynamic(vector<pair<Point, InvEnergy> > points1,
		 vector<pair<Point, InvEnergy> > points2) const {
  //set up the selector with the product of the widths as weight
  vector<double> widths;
  Selector<pair<int, int>, InvEnergy2> sel;
  for ( int m = 0; m < int(points1.size()); m++ ) {
    for ( int n = 0; n < int(points2.size()); n++ ) {
      sel.insert(points1[m].second*points2[n].second, make_pair(m, n) );
    }
  }
  pair<int, int> ij = sel[UseRandom::rnd()];
  Point p1 = points1[ij.first].first;
  Point p2 = points2[ij.second].first;

  //set width the geometric mean, let the width set the scale, using the general
  //scale 3/Gev, which is roufly the size of the proton, and the confimenet range
  InvEnergy dynWidth = width()/(3.0/GeV)*sqrt(points1[ij.first].second*points2[ij.second].second);

  // cout << "picked point (" << p1.x()*GeV << ", " << p1.y()*GeV << "), with n,m "
  //      << ij.first << ", " << ij.second << endl;
  // cout << "widths are " << points1[ij.first].second*GeV << " and "
  //      << points2[ij.second].second*GeV << ", mean " << dynWidth*GeV
  //      << ", fij weight: " << sqr(dynWidth)*GeV2 << endl;

  //first generate the rotation of the second state.
  double theta = 2.0*Constants::pi*rnd();
  ImpactParameters ip = ImpactParameters(Point(), theta, 1./GeV2);
  p2 = ip.rotate(p2);

  //generate the distance between the selected points
  //overestimate b*exp(-b/w) with 2/w*exp(-b/(2w)-1)
  InvEnergy b;
  do b = -2.0*dynWidth*log(UseRandom::rnd());
  while ( b/(2.0*dynWidth)*exp(1-b/(2.0*dynWidth)) < UseRandom::rnd() );

  //generate direction between states
  double phi = 2.0*Constants::pi*rnd();

  // cout << "picked a p-p distance " << b*GeV << endl;

  //find the required translation for the centers of each state.
  Point bVec = Point(b*cos(phi), b*sin(phi)) + p1 - p2;

  //find what the probability density in d2b is for that translation,
  //summed over all pairs of particles.
  Energy2 f = ZERO;
  double sumWeights = 0.0;
  for ( int m = 0; m < int(points1.size()); m++ ) {
    for ( int n = 0; n < int(points2.size()); n++ ) {
      Point bmn = bVec - points1[m].first + ip.rotate(points2[n].first);
      //find the width for this pair
      InvEnergy thisWidth = sqrt(points1[m].second*points2[n].second);
      //add the prob density of the dist from this pair, with the extra weight
      f += exp(-bmn.pt()/thisWidth)/(2.0*Constants::pi)*GeV2;
      //keep track of the extra weights of the prob dens
      sumWeights += sqr(thisWidth)*GeV2;
      // if ( exp(-bmn.pt()/thisWidth)/(2.0*Constants::pi*sqr(thisWidth))/GeV2 > 0.1 )
      // 	cout << "m,n: " << m << ", " << n << ", fnm(b): " << exp(-bmn.pt()/thisWidth)/(2.0*Constants::pi) << endl;
    }
  }
  //normalise by dividing by the sum off all the weights for the individual p-p prob dists
  // cout << "sum of all weights: " << sumWeights << endl;
  f /= sumWeights;

  // cout << "total prob dens " << f/GeV2 << endl;

  //weight is 1/f(bVec)
  InvEnergy2 w = 1.0/f;

  // cout << "weight " << w*GeV2 << endl;

  // cout << points1.size() << ", " << points2.size() << ", dynwidth: "
  //      << GeV*width()*sqrt(1+points1.size()*points2.size())/40.0
  //      << ", bprime: " << b*GeV << ", 100xf: " << 100.0*f/GeV2 << endl;

 if (Debug::level > 5 )
    cout << "setting up intercept at (" << p1.x()*GeV << ", " << p1.y()*GeV << ")" << endl;
  if (Debug::level > 5 )
    cout << "other particle should be at (" << (p1+bVec).x()*GeV << ", "
	 << (p1+bVec).y()*GeV << ")" << endl;

  // return ImpactParameters(Point(ZERO, 3.0/GeV), theta, 100./GeV2); //for fix b runs
 
  return ImpactParameters(bVec, theta, w);
}

void ImpactParameterGenerator::persistentOutput(PersistentOStream & os) const {
  os << ounit(theWidth, InvGeV);
}

void ImpactParameterGenerator::persistentInput(PersistentIStream & is, int) {
  is >> iunit(theWidth, InvGeV);
}



// Static variable needed for the type description system in ThePEG.
#include "ThePEG/Utilities/DescribeClass.h"
DescribeClass<ImpactParameterGenerator,HandlerBase>
  describeDIPSYImpactParameterGenerator("DIPSY::ImpactParameterGenerator",
					"libAriadne5.so libDIPSY.so");


void ImpactParameterGenerator::Init() {

  static ClassDocumentation<ImpactParameterGenerator> documentation
    ("The objective of the ImpactParameterGenerator class is to generate "
     "ImpactParameters objects to be used when two "
     "<code>DipoleState</code>s have been generated w.r.t.  origo in the "
     "transverse to displace and rotate one of them before "
     "collision. This base class will generate impact parameters "
     "according to a Gaussian distribution, and the weigth in the "
     "produced ImpactParameters objects is set accordingly. Sub-classes "
     "may override the generate() function to use any distribution as "
     "long as the weight is set accordingly.");

  static Parameter<ImpactParameterGenerator,InvEnergy> interfaceWidth
    ("Width",
     "The width of the generated distribution.",
     &ImpactParameterGenerator::theWidth, InvGeV, 10.0*InvGeV, 0.0*InvGeV,
     0*InvGeV, true, false, Interface::lowerlim);

}

