// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SimpleNucleus class.
//

#include "SimpleNucleus.h"
#include "SimpleNucleusState.h"
#include "DipoleEventHandler.h"
#include "NucleusData.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Repository/Repository.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Utilities/EnumIO.h"
#include "ThePEG/PDT/StandardMatchers.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include "gsl/gsl_sf_erf.h"

using namespace DIPSY;

SimpleNucleus::~SimpleNucleus() {}

void SimpleNucleus::initialize(const DipoleEventHandler & eh) {
  if ( proton() ) proton()->initialize(eh);
  if ( neutron() ) neutron()->initialize(eh);
  if ( nucleon() ) nucleon()->initialize(eh);
}

Energy2 SimpleNucleus::m2() const {
  if ( abs(particle()->id()) > 1000000000 ) return sqr(particle()->mass());
  return sqr(A*particle()->mass());
}

DipoleStatePtr SimpleNucleus::
generate(const DipoleEventHandler & eh, Energy plus) {

  vector<Point> positions;

  while ( int(positions.size()) < A ) {

    Point pos(2.0*UseRandom::rnd(-R, R),
	      2.0*UseRandom::rnd(-R, R),
	      2.0*UseRandom::rnd(-R, R));
    Length r = pos.mag();
    if ( !UseRandom::rndbool(ws(r)) ) continue;
    int ntry = 0;
    do {
      bool overlap = false;
      for ( int i = 0, N = positions.size(); i < N && !overlap; ++i )
        if ( (positions[i] - pos).mag() < abs(Rn) ) overlap = true;

      if ( !overlap ) {
	positions.push_back(pos);
	break;
      }

      if ( nTry == 0 ) positions.clear();
      else if ( ++ntry < nTry ) {
	double theta = acos(UseRandom::rnd(-1.0,1.0));
	double phi   = UseRandom::rnd( 0.,2.*M_PI);
	pos = Point(r*sin(theta)*cos(phi), r*sin(theta)*sin(phi), r*cos(theta));
      }
    } while ( ++ntry < nTry );

  }

  if ( rdist || recenter() ) {
    Point cms;
    for ( int i = 0, N = positions.size(); i < N; ++i )
      cms += positions[i];
    cms /= double(positions.size());
    for ( int i = 0, N = positions.size(); i < N; ++i ) {
      if ( recenter() ) positions[i] -= cms;
      double r = positions[i].mag()/femtometer;
      rdist->fill(r, sqr(1.0/r));
    }
  }


  Energy minus = m2()/plus;
  return new_ptr(SimpleNucleusState(eh, *this, plus, minus, positions,
				    new_ptr(WFInfo(this, R/Constants::hbarc))));

}

void SimpleNucleus::persistentOutput(PersistentOStream & os) const {
  os << A << Z << ounit(R, femtometer) << ounit(a, femtometer)
     << w << wf << wfp << wfn << ounit(Rn, femtometer)
     << oenum(useInterNucleonSwing) << nTry << doRecenter;
}

void SimpleNucleus::persistentInput(PersistentIStream & is, int) {
  is >> A >> Z >> iunit(R, femtometer) >> iunit(a, femtometer)
     >> w >> wf >> wfp >> wfn >> iunit(Rn, femtometer)
     >> ienum(useInterNucleonSwing) >> nTry >> doRecenter;
}

string SimpleNucleus::setParameters(string cmd) {
  string name = StringUtils::car(cmd);
  cmd = StringUtils::cdr(cmd);

  if ( Ptr<NucleusData>::pointer nucl =
       Repository::GetPtr<Ptr<NucleusData>::pointer>(name) ) {
    A = nucl->A();
    Z = nucl->Z();
    setParticle(nucl);
  } else {
    if ( name == "He" ) {
      A = 4;
      Z = 2;
    } else if ( name == "Li" ) {
      A = 6;
      Z = 3;
    } else if ( name == "C" ) {
      A = 12;
      Z = 6;
    } else if ( name == "O" ) {
      A = 16;
      Z = 8;
    } else if ( name == "Al" ) {
      A = 27;
      Z = 13;
    } else if ( name == "S" ) {
      A = 32;
      Z = 16;
    } else if ( name == "Ca" ) {
      A = 40;
      Z = 20;
    } else if ( name == "Ni" ) {
      A = 58;
      Z = 28;
    } else if ( name == "Cu" ) {
      A = 63;
      Z = 29;
    } else if ( name == "W" ) {
      A = 186;
      Z = 74;
    } else if ( name == "Au" ) {
      A = 197;
      Z = 79;
    } else if ( name == "Pb" ) {
      A = 208;
      Z = 82;
    } else if ( name == "U" ) {
      A = 238;
      Z = 92;
    } else {
      return string("Nucleus \"") + name +
	string("\" not available Please set individual parameters explicitly.");
    }
  }

  if (  StringUtils::car(cmd) == "GLISSANDO" ) {
    R = (1.1*pow(double(A),1.0/3.0) - 0.656*pow(double(A),-1.0/3.0))*femtometer;
    a = 0.459*femtometer;
    Rn = -0.9*femtometer;
    w = 0.0;
    return "";
  }

  switch ( abs(Z)*1000+abs(A) ) {
  case 2004:
    R = 1.71*femtometer;
    a = 0.5*femtometer;
    w = 0.0;
    break;
  case 3006:
    R = 1.96*femtometer;
    a = 0.5*femtometer;
    w = 0.0;
    break;
  case 6012:
    R = 2.47*femtometer;
    a = 0.00001*femtometer;
    w = 0.0;
    break;
  case 8016:
    R = 2.608*femtometer;
    a = 0.513*femtometer;
    w = -0.051;
    break;
  case 13027:
    R = 3.07*femtometer;
    a = 0.519*femtometer;
    w = 0.0;
    break;
  case 16032:
    R = 3.458*femtometer;
    a = 0.61*femtometer;
    w = 0.0;
    break;
  case 20040:
    R = 3.76*femtometer;
    a = 0.586*femtometer;
    w = -0.161;
    break;
  case 28058:
    R = 4.309*femtometer;
    a = 0.516*femtometer;
    w = -0.1308;
    break;
  case 29063:
    R = 4.2*femtometer;
    a = 0.596*femtometer;
    w = 0.0;
    break;
  case 74186:
    R = 6.51*femtometer;
    a = 0.535*femtometer;
    w = 0.0;
    break;
  case 79197:
    R = 6.38*femtometer;
    a = 0.535*femtometer;
    w = 0.0;
    break;
  case 82208:
    R = 6.68*femtometer;
    a = 0.546*femtometer;
    w = 0.0;
    break;
  case 92238:
    R = 6.68*femtometer;
    a = 0.6*femtometer;
    w = 0.0;
    break;
  default:
    return string("Nucleus \"") + name +
      string("\" not available. Please set individual parameters explicitly.");
  }
  return "";
}


void SimpleNucleus::dofinish() {
  WaveFunction::dofinish();
}

void SimpleNucleus::doinitrun() {
  WaveFunction::doinitrun();
  if ( HistFacPtr fac = generator()->histogramFactory() ) {
    fac->initrun();
    fac->registerClient(this);
    fac->mkdirs("/DIPSY/SimpleNucleus");
    rdist = fac->createHistogram1D("/DIPSY/SimpleNucleus/RDist",100,0.0,10.0);
  }
}


// Static variable needed for the type description system in ThePEG.
#include "ThePEG/Utilities/DescribeClass.h"
DescribeClass<SimpleNucleus,DIPSY::WaveFunction>
  describeDIPSYSimpleNucleus("DIPSY::SimpleNucleus", "SimpleNucleus.so");

void SimpleNucleus::Init() {

  static ClassDocumentation<SimpleNucleus> documentation
    ("The SimpleNucleus class represents the unevolved necleus wave function "
     "described in terms of particles distributed according to a Wood-Saxon "
     "potential.");

  static Parameter<SimpleNucleus,int> interfaceA
    ("A",
     "The mass number of the nucleus.",
     &SimpleNucleus::A, 208, 0, 0,
     true, false, Interface::lowerlim);

  static Parameter<SimpleNucleus,int> interfaceZ
    ("Z",
     "The charge of the nucleus.",
     &SimpleNucleus::Z, 82, 0, 0,
     true, false, Interface::lowerlim);

  static Parameter<SimpleNucleus,Length> interfaceR
    ("R",
     "The size of the nucleus in units of fermi.",
     &SimpleNucleus::R, femtometer, 6.68*femtometer, 0.0*femtometer,
     0.0*femtometer, true, false, Interface::lowerlim);

  static Parameter<SimpleNucleus,Length> interfacea
    ("a",
     "The suppression at the edge of the nucleus in units of fermi.",
     &SimpleNucleus::a, femtometer, 0.546*femtometer, 0.0*femtometer,
     0.0*femtometer, true, false, Interface::lowerlim);

  static Parameter<SimpleNucleus,double> interfacew
    ("w",
     "Optional quadratic term in Wood-Saxon potential",
     &SimpleNucleus::w, 0.0, 0, 0,
     true, false, Interface::nolimits);

  static Reference<SimpleNucleus,WaveFunction> interfaceWF
    ("WF",
     "The wave function describing individual nucleons.",
     &SimpleNucleus::wf, true, false, true, true, false);

  static Reference<SimpleNucleus,WaveFunction> interfaceWFP
    ("WFP",
     "The wave function describing individual protons.",
     &SimpleNucleus::wfp, true, false, true, true, false);

  static Reference<SimpleNucleus,WaveFunction> interfaceWFN
    ("WFN",
     "The wave function describing individual neutrons.",
     &SimpleNucleus::wfn, true, false, true, true, false);

  static Parameter<SimpleNucleus,Length> interfaceRn
    ("Rn",
     "Size in units of fermi of an individual nucleon for the exclusion "
     "mechanism.",
     &SimpleNucleus::Rn, femtometer, 1.3*femtometer, 0.0*femtometer,
     0.0*femtometer, true, false, Interface::nolimits);

  static Command<SimpleNucleus> interfaceSetNucleus
    ("SetNucleus",
     "Set parameters of the given nucleus according to measured values "
     "if available.",
     &SimpleNucleus::setParameters, true);


  static Switch<SimpleNucleus,bool> interfaceInterNucleonSwing
    ("InterNucleonSwing",
     "Allow for swings between the nucleons in the evolution.",
     &SimpleNucleus::useInterNucleonSwing, true, true, false);
  static SwitchOption interfaceInterNucleonSwingYes
    (interfaceInterNucleonSwing,
     "On",
     "Allow for inter-nucleon swings.",
     true);
  static SwitchOption interfaceInterNucleonSwingNo
    (interfaceInterNucleonSwing,
     "Off",
     "Disallow inter-nucleon swings.",
     false);
  
  static Parameter<SimpleNucleus,int> interfaceNTry
    ("NTry",
     "Number of times to try new angles for the same r-value if "
     "overlapping nucleons are found. If zero, the whole set of "
     "generated positions is discarded if anoverlap is found.",
     &SimpleNucleus::nTry, 1, 0, 0,
     true, false, Interface::lowerlim);

  static Switch<SimpleNucleus,bool> interfaceRecenter
    ("Recenter",
     "Set true if the nucleons should be recentered to put the center of mass at zero impact parameter after distribution of nucleons.",
     &SimpleNucleus::doRecenter, false, true, false);
  static SwitchOption interfaceRecenterNoRecenter
    (interfaceRecenter,
     "NoRecenter",
     "Do not recenter the nucleons.",
     false);
  static SwitchOption interfaceRecenterRecenter
    (interfaceRecenter,
     "Recenter",
     "Recenter the nucleons.",
     true);

  interfaceSetNucleus.rank(10);
  interfaceA.rank(9);
  interfaceZ.rank(8);
  interfaceR.rank(7);
  interfacea.rank(6);
  interfaceWF.rank(5);

}
