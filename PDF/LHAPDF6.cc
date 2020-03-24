// -*- C++ -*-
//
// LHAPDF6.cc is a part of ThePEG - Toolkit for HEP Event Generation
// // Copyright (C) 2014-2019 Leif Lonnblad, David Grellscheid
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHAPDF class.
//

#include "LHAPDF6.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Interface/Deleted.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "LHAPDF/LHAPDF.h"

using std::vector;
using std::string;
using std::pair;

using ThePEG::GeV2;
using ThePEG::cPDVector;

ThePEG::LHAPDF::LHAPDF()
  : thePDF(), thePDFName("THEPEG_NO_PDFSET_CHOSEN"), 
    theMember(0), theMaxFlav(5),
    xMin(0.0), xMax(1.0), Q2Min(ZERO), Q2Max(Constants::MaxEnergy2) {}

ThePEG::LHAPDF::LHAPDF(const LHAPDF & x)
  : PDFBase(x), 
    thePDF(), thePDFName(x.thePDFName), 
    theMember(x.theMember),
    theMaxFlav(x.theMaxFlav),
    xMin(x.xMin), xMax(x.xMax), Q2Min(x.Q2Min), Q2Max(x.Q2Max) {}

ThePEG::IBPtr ThePEG::LHAPDF::clone() const {
  return new_ptr(*this);
}

ThePEG::IBPtr ThePEG::LHAPDF::fullclone() const {
  return new_ptr(*this);
}

void ThePEG::LHAPDF::initPDFptr() {
  ::LHAPDF::setVerbosity(std::max(0, 
                         ThePEG::Debug::level - 1));
  if (    thePDF 
       && thePDF->set().name() == thePDFName 
       && thePDF->memberID() == theMember ) 
    return;
  delete thePDF;
  thePDF = ::LHAPDF::mkPDF(thePDFName, theMember);
  xMin = thePDF->xMin();
  xMax = thePDF->xMax();
  Q2Min = thePDF->q2Min() * GeV2;
  Q2Max = thePDF->q2Max() * GeV2;
}

void ThePEG::LHAPDF::doinit() {
  PDFBase::doinit();
  initPDFptr();
}

void ThePEG::LHAPDF::dofinish() {
  PDFBase::dofinish();
  delete thePDF;
  thePDF = 0;
}

void ThePEG::LHAPDF::doinitrun() {
  PDFBase::doinitrun();
  initPDFptr();
}

void ThePEG::LHAPDF::setPDFName(string name) {
  if ( ::LHAPDF::endswith(name, ".LHgrid") ) {
  	name = name.substr(0, name.size() - 7);
  }
  else if ( ::LHAPDF::endswith(name, ".LHpdf") ) {
  	name = name.substr(0, name.size() - 6);
  }

  // fix the eternal typo
  if ( name == "cteq6ll" ) name = "cteq6l1";

  if ( ::LHAPDF::contains(::LHAPDF::availablePDFSets(), name) ) {
    thePDFName = name;
    theMember = 0;
  }
  else {
    Throw<ThePEG::LHAPDF::NotInstalled>()
    	<< "'set " << fullName() << ":PDFName "
    	<< name << "': PDF not installed. Try 'lhapdf install'.\n"
    	<< Exception::setuperror;
  }
}

void ThePEG::LHAPDF::setPDFMember(int member) {
  try {
    ::LHAPDF::PDFInfo * test = 
  	::LHAPDF::mkPDFInfo(thePDFName, member);
    if ( test )
      theMember = member;
    delete test;
  }
  catch (::LHAPDF::ReadError & e) {
   Throw<ThePEG::LHAPDF::NotInstalled>()
    	<< e.what() << Exception::setuperror;
  }
}

string ThePEG::LHAPDF::doTest(string input) {
  double x = 0;
  Energy2 Q2 = ZERO;
  Energy2 P2 = ZERO;
  istringstream is(input);
  is >> x >> iunit(Q2, GeV2) >> iunit(P2, GeV2);
  initPDFptr();
  ostringstream os;
  for ( int i = 0; i < 13; ++i ) os << " " << thePDF->xfxQ2(i,x,Q2/GeV2);
  return os.str();
}  

bool ThePEG::LHAPDF::canHandleParticle(tcPDPtr particle) const {
  using namespace ParticleID;
  return abs(particle->id()) == pplus || abs(particle->id()) == n0;
}

cPDVector ThePEG::LHAPDF::partons(tcPDPtr particle) const {
  // We assume that all standard partons can be extracted.
  const ::LHAPDF::PDFSet & pdfset = ::LHAPDF::getPDFSet(thePDFName);
  const vector<int> & flavs = 
    pdfset.get_entry_as< vector<int> >("Flavors");

  cPDVector ret;
  ret.reserve( flavs.size() );
  if ( canHandleParticle(particle) ) {
    for ( size_t i=0; i < flavs.size(); ++i ) {
      ret.push_back(getParticleData(flavs[i]));
    }
  }
  assert( ret.size() == flavs.size() );
  return ret;
}

double ThePEG::LHAPDF::xfx(tcPDPtr particle, tcPDPtr parton,
			   Energy2 partonScale,
			   double x, double, Energy2) const {
  // Here we should return the actual density.
  using namespace ThePEG::ParticleID;

  double Q2 = partonScale/GeV2;

   if ( ! thePDF->inRangeXQ2(x, Q2) ) {
       switch ( rangeException ) {
       case rangeThrow: Throw<Exception>()
 	<< "Momentum fraction (x=" << x << ") or scale (Q2=" << Q2
 	<< " GeV^2) was outside of limits in PDF " << name() << "."
 	<< Exception::eventerror;
	 break;
       case rangeZero:
 	 return 0.0;
       case rangeFreeze:
 	x = min(max(x, xMin), xMax);
        Q2 = min(max(Q2, Q2Min/GeV2), Q2Max/GeV2);
       }
   } 

   int pid = parton->id();
   int abspid = abs(pid);
   
  switch ( pid ) {
  case t:
  case tbar:
  case b:
  case bbar:
  case c:
  case cbar:
    return maxFlav() < abspid ? 0.0 : thePDF->xfxQ2(pid,x,Q2);
  case s:
  case sbar:
    return thePDF->xfxQ2(pid,x,Q2);
  case u:
    switch ( particle->id() ) {
    case n0:        return thePDF->xfxQ2(d   ,x,Q2);
    case pbarminus: return thePDF->xfxQ2(ubar,x,Q2);
    case nbar0:     return thePDF->xfxQ2(dbar,x,Q2);
    case pplus:
    default:        return thePDF->xfxQ2(u   ,x,Q2);
    }
  case ubar:
    switch ( particle->id() ) {
    case n0:        return thePDF->xfxQ2(dbar,x,Q2);
    case pbarminus: return thePDF->xfxQ2(u   ,x,Q2);
    case nbar0:     return thePDF->xfxQ2(d   ,x,Q2);
    case pplus:
    default:        return thePDF->xfxQ2(ubar,x,Q2);
    }
  case d:
    switch ( particle->id() ) {
    case n0:        return thePDF->xfxQ2(u   ,x,Q2);
    case pbarminus: return thePDF->xfxQ2(dbar,x,Q2);
    case nbar0:     return thePDF->xfxQ2(ubar,x,Q2);
    case pplus:
    default:        return thePDF->xfxQ2(d   ,x,Q2);
    }
  case dbar:
    switch ( particle->id() ) {
    case n0:        return thePDF->xfxQ2(ubar,x,Q2);
    case pbarminus: return thePDF->xfxQ2(d   ,x,Q2);
    case nbar0:     return thePDF->xfxQ2(u   ,x,Q2);
    case pplus:
    default:        return thePDF->xfxQ2(dbar,x,Q2);
    }
  case g:
    return thePDF->xfxQ2(g,x,Q2);
  case ParticleID::gamma:
    return thePDF->xfxQ2(ParticleID::gamma,x,Q2);
  }
  return 0.0;
}

double ThePEG::LHAPDF::xfvl(tcPDPtr particle, tcPDPtr parton,
			    Energy2 partonScale,
			    double l, Energy2 particleScale) const {
  using Math::exp1m;
  return xfvx(particle, parton, partonScale,
	      exp(-l), exp1m(-l), particleScale);
}

double ThePEG::LHAPDF::xfvx(tcPDPtr particle, tcPDPtr parton,
			    Energy2 partonScale,
			    double x, double, Energy2) const {
  // Here we should return the actual valence density. This will only
  // work properly for nucleons
  using namespace ThePEG::ParticleID;

  double Q2 = partonScale / GeV2;

   if ( ! thePDF->inRangeXQ2(x, Q2) ) {
       switch ( rangeException ) {
       case rangeThrow: Throw<Exception>()
 	<< "Momentum fraction (x=" << x << ") or scale (Q2=" << Q2
 	<< " GeV^2) was outside of limits in PDF " << name() << "."
 	<< Exception::eventerror;
	 break;
       case rangeZero:
 	 return 0.0;
       case rangeFreeze:
 	x = min(max(x, xMin), xMax);
        Q2 = min(max(Q2, Q2Min/GeV2), Q2Max/GeV2);
       }
   } 

  switch ( parton->id() ) {
  case t:
  case tbar:
  case b:
  case bbar:
  case c:
  case cbar:
  case s:
  case sbar:
  case ParticleID::gamma:
    return 0.0;
  case u:
    switch ( particle->id() ) {
    case n0:    return thePDF->xfxQ2(d,x,Q2) - thePDF->xfxQ2(dbar,x,Q2);
    case pbarminus: return 0.0;
    case nbar0: return 0.0;
    case pplus: return thePDF->xfxQ2(u,x,Q2) - thePDF->xfxQ2(ubar,x,Q2);
    default: return 0.0;
    }
  case ubar:
    switch ( particle->id() ) {
    case n0: return 0.0;
    case pbarminus: return thePDF->xfxQ2(u,x,Q2) - thePDF->xfxQ2(ubar,x,Q2);
    case nbar0:     return thePDF->xfxQ2(d,x,Q2) - thePDF->xfxQ2(dbar,x,Q2);
    case pplus:
    default: return 0.0;
    }
  case d:
    switch ( particle->id() ) {
    case n0:    return thePDF->xfxQ2(u,x,Q2) - thePDF->xfxQ2(ubar,x,Q2);
    case pbarminus: return 0.0;
    case nbar0: return 0.0;
    case pplus: return thePDF->xfxQ2(d,x,Q2) - thePDF->xfxQ2(dbar,x,Q2);
    default: return 0.0;
    }
  case dbar:
    switch ( particle->id() ) {
    case n0: return 0.0;
    case pbarminus: return thePDF->xfxQ2(d,x,Q2) - thePDF->xfxQ2(dbar,x,Q2);
    case nbar0:     return thePDF->xfxQ2(u,x,Q2) - thePDF->xfxQ2(ubar,x,Q2);
    case pplus:
    default: return 0.0;
    }
  case ParticleID::g:
    return 0.0;
  }
  return 0.0;
}

double ThePEG::LHAPDF::xfsx(tcPDPtr particle, tcPDPtr parton,
			    Energy2 partonScale, double x,
			    double, Energy2) const {
  // Here we should return the actual density.
  using namespace ThePEG::ParticleID;

  double Q2 = partonScale / GeV2;

   if ( ! thePDF->inRangeXQ2(x, Q2) ) {
       switch ( rangeException ) {
       case rangeThrow: Throw<Exception>()
 	<< "Momentum fraction (x=" << x << ") or scale (Q2=" << Q2
 	<< " GeV^2) was outside of limits in PDF " << name() << "."
 	<< Exception::eventerror;
	 break;
       case rangeZero:
 	 return 0.0;
       case rangeFreeze:
 	x = min(max(x, xMin), xMax);
        Q2 = min(max(Q2, Q2Min/GeV2), Q2Max/GeV2);
       }
   } 
   int pid = parton->id();
   int abspid = abs(pid);
   
  switch ( pid ) {
  case t:
  case tbar:
  case b:
  case bbar:
  case c:
  case cbar:
    return maxFlav() < abspid ? 0.0 : thePDF->xfxQ2(pid,x,Q2);
  case s:
  case sbar:
    return thePDF->xfxQ2(pid,x,Q2);
  case u:
    switch ( particle->id() ) {
    case n0:        return thePDF->xfxQ2(dbar,x,Q2);
    case pbarminus: return thePDF->xfxQ2(ubar,x,Q2);
    case nbar0:     return thePDF->xfxQ2(dbar,x,Q2);
    case pplus:     return thePDF->xfxQ2(ubar,x,Q2);
    default:        return thePDF->xfxQ2(u   ,x,Q2);
    }
  case ubar:
    switch ( particle->id() ) {
    case n0:        return thePDF->xfxQ2(dbar,x,Q2);
    case pbarminus: return thePDF->xfxQ2(ubar,x,Q2);
    case nbar0:     return thePDF->xfxQ2(dbar,x,Q2);
    case pplus:
    default:        return thePDF->xfxQ2(ubar,x,Q2);
    }
  case d:
    switch ( particle->id() ) {
    case n0:        return thePDF->xfxQ2(ubar,x,Q2);
    case pbarminus: return thePDF->xfxQ2(dbar,x,Q2);
    case nbar0:     return thePDF->xfxQ2(ubar,x,Q2);
    case pplus:     return thePDF->xfxQ2(dbar,x,Q2);
    default:        return thePDF->xfxQ2(d   ,x,Q2);
    }
  case dbar:
    switch ( particle->id() ) {
    case n0:        return thePDF->xfxQ2(ubar,x,Q2);
    case pbarminus: return thePDF->xfxQ2(dbar,x,Q2);
    case nbar0:     return thePDF->xfxQ2(ubar,x,Q2);
    case pplus:
    default:        return thePDF->xfxQ2(dbar,x,Q2);
    }
  case g:
    return thePDF->xfxQ2(g,x,Q2);
  case ParticleID::gamma:
    return thePDF->xfxQ2(ParticleID::gamma,x,Q2);
  }
  return 0.0;
}

void ThePEG::LHAPDF::persistentOutput(PersistentOStream & os) const {
  os << thePDFName << theMember << theMaxFlav
     << xMin << xMax << ounit(Q2Min, GeV2) << ounit(Q2Max, GeV2);
}

void ThePEG::LHAPDF::persistentInput(PersistentIStream & is, int) {
  is >> thePDFName >> theMember >> theMaxFlav
     >> xMin >> xMax >> iunit(Q2Min, GeV2) >> iunit(Q2Max, GeV2);
  initPDFptr();
}

ThePEG::DescribeClass<ThePEG::LHAPDF, ThePEG::PDFBase>
describeLHAPDF("ThePEG::LHAPDF", "ThePEGLHAPDF.so");

void ThePEG::LHAPDF::Init() {

  static ClassDocumentation<LHAPDF> documentation
    ("The LHAPDF class inherits from PDFBase and implements an interface "
     "to the LHAPDF library of parton density function parameterizations. "
     "This class is available even if LHAPDF was not properly installed "
     "when ThePEG was installed, but will then produce an error in the "
     "initialization. Note that the valence densities from the xfvx() and "
     "xfvl() function will only work properly for nucleons. All other "
     "particles will have zero valence densities.");

  static Deleted<LHAPDF> interfacePType
    ("PType",
     "The LHAPDFv6 interface currently does not support pi.");

  static Parameter<LHAPDF,string> interfacePDFName
    ("PDFName",
     "The name of the PDF set to be used. Should correspond to "
     "the LHAPDF v6 name.",
     &ThePEG::LHAPDF::thePDFName, "THEPEG_NO_PDFSET_CHOSEN", true, false,
     &ThePEG::LHAPDF::setPDFName);

  static Deleted<LHAPDF> interfacePDFNumber
    ("PDFNumber",
     "Not implemented in the LHAPDFv6 interface. "
     "Use :PDFName and :Member instead.");

  static Parameter<LHAPDF,int> interfaceMember
    ("Member",
     "The chosen member of the selected PDF set.",
     &ThePEG::LHAPDF::theMember, 0, 0, 0, 
     true, false, Interface::lowerlim,
     &ThePEG::LHAPDF::setPDFMember);

  static Deleted<LHAPDF> interfacePDFLIBNumbers
    ("PDFLIBNumbers",
     "Not implemented in the LHAPDFv6 interface. "
     "Use :PDFName and :Member instead.");

  static Deleted<LHAPDF> interfaceEnablePartonicGamma
    ("EnablePartonicGamma",
     "Not required in LHAPDFv6.");

  static Deleted<LHAPDF> interfacePhotonOption
    ("PhotonOption",
     "Not required in LHAPDFv6.");

  static Deleted<LHAPDF> interfaceMaxNSet
    ("MaxNSet",
     "Not required in LHAPDFv6.");

  static Command<LHAPDF> interfaceTest
    ("Test",
     "Write out the values of the chosen PDF set using the x, Q2 and P2 "
     "parameters supplied.",
     &ThePEG::LHAPDF::doTest, true);

  static Deleted<LHAPDF> interfaceVerboseLevel
    ("VerboseLevel",
     "LHAPDFv6 uses general debug level instead.");

  static Parameter<LHAPDF,int> interfaceMaxFlav
    ("MaxFlav",
     "The maximum number of flavours for which non-zero densities are "
     "reported. The actual number of flavours may be less depending on "
     "the chosen PDF set.",
     &ThePEG::LHAPDF::theMaxFlav, 5, 3, 0,
     true, false, Interface::lowerlim);

  interfacePDFName.rank(9);
  interfaceMember.rank(8);

}

