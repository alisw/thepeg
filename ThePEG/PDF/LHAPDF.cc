// -*- C++ -*-
//
// LHAPDF.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LHAPDF class.
//

// macros are passed in from -D compile flag
#ifndef THEPEG_PKGDATADIR
#error Makefile.am needs to define THEPEG_PKGDATADIR
#endif
#ifndef LHAPDF_PKGDATADIR
#error Makefile.am needs to define LHAPDF_PKGDATADIR
#endif

#include "LHAPDF.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Command.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Utilities/EnumIO.h"
#include "ThePEG/Utilities/Debug.h"
#include "ThePEG/Utilities/Throw.h"
#include "config.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#ifdef ThePEG_HAS_FPU_CONTROL
#include <fpu_control.h>
#endif
#ifdef ThePEG_HAS_FENV
#include <fenv.h>
#endif

using namespace ThePEG;

typedef double F77ThePEGDouble;
typedef int F77ThePEGInteger;

extern "C" {
  void setlhaparm_(const char *, F77ThePEGInteger);
  void initpdfsetbynamem_(F77ThePEGInteger &, const char *, F77ThePEGInteger);
  void initpdfm_(F77ThePEGInteger &,F77ThePEGInteger  &);
  void evolvepdfm_(F77ThePEGInteger &, F77ThePEGDouble &,
		   F77ThePEGDouble &, F77ThePEGDouble *);
  void evolvepdfphotonm_(F77ThePEGInteger &, F77ThePEGDouble &,
			 F77ThePEGDouble &, F77ThePEGDouble *, F77ThePEGDouble &);
  void evolvepdfpm_(F77ThePEGInteger &, F77ThePEGDouble &,
		    F77ThePEGDouble &, F77ThePEGDouble &,
		    F77ThePEGInteger &, F77ThePEGDouble *);
  void numberpdfm_(F77ThePEGInteger &, F77ThePEGInteger &);
  void getnfm_(F77ThePEGInteger &, F77ThePEGInteger &);
  void lhaprint_(F77ThePEGInteger &);
}

struct TmpMaskFpuDenorm {
#ifdef ThePEG_HAS_FPU_CONTROL_NEVER
  fpu_control_t oldcw;
  TmpMaskFpuDenorm() {
    volatile fpu_control_t cw;
    _FPU_GETCW(cw);
    oldcw = cw;
    cw |= (_FPU_MASK_DM);
    if ( cw != oldcw ) _FPU_SETCW(cw);
  }
  ~TmpMaskFpuDenorm() {
    volatile fpu_control_t cw;
    _FPU_GETCW(cw);
    if ( cw != oldcw ) {
      cw = oldcw;
      _FPU_SETCW(cw);
    }
  }
#else
#ifdef ThePEG_HAS_FENV
  int oldexcept;
  TmpMaskFpuDenorm() {
    oldexcept = fegetexcept();
    fedisableexcept(FE_INEXACT);
  }
  ~TmpMaskFpuDenorm() {
    feenableexcept(oldexcept);
  }
#else  
  TmpMaskFpuDenorm() {
    Debug::maskFpuDenorm();
  }
#endif
#endif
};

#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/EnumParticles.h"

LHAPDF::LHAPDF()
  : thePType(nucleonType), thePDFName("cteq6ll.LHpdf"), theMember(0),
    thePhotonOption(7), enablePartonicGamma(false),
    theVerboseLevel(0), theMaxFlav(5), nset(-1),
    lastQ2(-1.0*GeV2), lastX(-1.0), lastP2(-1.0*GeV2),
    xMin(0.0), xMax(1.0), Q2Min(ZERO), Q2Max(Constants::MaxEnergy2) {}

LHAPDF::LHAPDF(const LHAPDF & x)
  : PDFBase(x), thePType(x.thePType), thePDFName(x.thePDFName),
    theMember(x.theMember),
    thePhotonOption(x.thePhotonOption),
    enablePartonicGamma(x.enablePartonicGamma),
    theVerboseLevel(x.theVerboseLevel), theMaxFlav(x.theMaxFlav), nset(x.nset),
    lastQ2(-1.0*GeV2), lastX(-1.0), lastP2(-1.0*GeV2),
    xMin(x.xMin), xMax(x.xMax), Q2Min(x.Q2Min), Q2Max(x.Q2Max) {}

IBPtr LHAPDF::clone() const {
  return new_ptr(*this);
}

IBPtr LHAPDF::fullclone() const {
  return new_ptr(*this);
}

void LHAPDF::doinit() {
  PDFBase::doinit();
  setMinMax();
  checkInit();
}

void LHAPDF::dofinish() {
  PDFBase::dofinish();
}

void LHAPDF::doinitrun() {
  PDFBase::doinitrun();
  checkInit();
}

void LHAPDF::throwNotInstalled() {
  throw LHAPDF::NotInstalled()
    << "Tried to initialize a LHAPDF object, "
    << "but the LHAPDF library was not installed" << Exception::runerror;
}

void LHAPDF::initpdfsetm() const {
  TmpMaskFpuDenorm fpumask;
  F77ThePEGInteger iset = nset + 1;
  initpdfsetbynamem_(iset, PDFName().c_str(), PDFName().length());
  lastNames[nset] = PDFName();
}

void LHAPDF::initpdfm() const {
  TmpMaskFpuDenorm fpumask;
  F77ThePEGInteger iset = nset + 1;
  F77ThePEGInteger mem = member();
  initpdfm_(iset, mem);
  lastMem[nset] = member();
  lastReset();
}

void LHAPDF::lastReset() const {
  lastQ2 = -1.0*GeV2;
  lastX = -1.0;
  lastP2 = -1.0*GeV2;
}

void LHAPDF::setnset() const {
  TmpMaskFpuDenorm fpumask;
  F77ThePEGInteger i = !theVerboseLevel;
  lhaprint_(i);
  if ( nset < 0 || nset >= MaxNSet) {
    // First check if any other nset-value is using 'our' pdf set.
    for ( nset = 0; nset < min(lastNSet, MaxNSet); ++nset )
      if ( lastMem[nset] == member() && lastNames[nset] == PDFName() )
	return;
    // Otherwise book a new nset.
    nset = (lastNSet++)%MaxNSet;
  }
}

int LHAPDF::getMaxMember() const {
  TmpMaskFpuDenorm fpumask;
  checkInit();
  F77ThePEGInteger iset = nset + 1;
  F77ThePEGInteger maxmem = 1;
  numberpdfm_(iset, maxmem);
  return maxmem;
}

int LHAPDF::getMaxFlav() const {
  TmpMaskFpuDenorm fpumask;
  checkInit();
  F77ThePEGInteger iset = nset + 1;
  F77ThePEGInteger maxflav = 1;
  getnfm_(iset, maxflav);
  return maxflav>=0 ? min(maxflav, maxFlav()) : maxFlav();
}

void LHAPDF::setMaxNSet(int n) {
  MaxNSet = n;
  lastNames.resize(MaxNSet);
  lastMem.resize(MaxNSet);
}

int LHAPDF::getMaxNSet() const {
  return MaxNSet;
}

void LHAPDF::checkInit() const {
  setlhaparm_("SILENT", 6);
  if ( nset < 0 || nset >= MaxNSet) {
    setnset();
    initpdfsetm();
    initpdfm();
  }
  else if ( PDFName() != lastNames[nset] ) {
    initpdfsetm();
    initpdfm();
  }
  else if ( member() != lastMem[nset] ) {
    initpdfm();
  }
}

std::string LHAPDF::getIndexPath() {
  // macro is passed in from -D compile flag
  return std::string(LHAPDF_PKGDATADIR) + "/PDFsets.index";
}

bool LHAPDF::openLHAIndex(ifstream & is) {
  if ( is.is_open() ) is.close();
  is.open( getIndexPath().c_str() );
  if ( is ) return true;
  is.clear();
  string instpath = std::string(THEPEG_PKGDATADIR) + "/PDFsets.index";
  is.open( instpath.c_str() );
  if ( is ) return true;
  is.clear();
  is.open("../PDF/PDFsets.index");
  if ( is ) return true;
  is.clear();
  is.open("../../ThePEG/PDF/PDFsets.index");
  if ( is ) return true;
  is.clear();
  is.open("./PDFsets.index");
  if ( is ) return true;
  is.clear(); 
  return false;
}

bool LHAPDF::indexLine(istream & is, int & set, int & mem, string & file,
		       int & pdftyp, int & pdfgup, int & pdfsup,
		       double & xmin, double & xmax,
		       double & q2min, double & q2max) const {
  string dummy;
  is >> set >> pdftyp >> pdfgup >> pdfsup >> file >> mem
     >> q2min >> q2max >> xmin >> xmax;
  // workaround for some C++11 implementations
  return bool(getline(is,dummy));
}

void LHAPDF::setMinMax() {
  ifstream is;
  if ( !openLHAIndex(is) ) Throw<InitException>()
    << "Could not open the LHAPDF index file so min/max values of "
    << "x and Q^2 could not be found." << Exception::warning;
  int set = 0;
  int mem = 0;
  string file;
  double xmin = 0.0;
  double xmax = 0.0;
  double q2min = 0.0;
  double q2max = 0.0;
  int pdftyp = 0;
  int pdfgup = 0;
  int pdfsup = 0;

  while ( indexLine(is, set, mem, file, pdftyp, pdfgup, pdfsup,
		    xmin, xmax, q2min, q2max) ) {
    if ( file == thePDFName && mem >= theMember ) {
      xMin = xmin;
      xMax = xmax;
      Q2Min = q2min*GeV2;
      Q2Max = q2max*GeV2;
      return;
    }
  }
}

void LHAPDF::setPDFNumber(int n) {
  ifstream is;
  if ( !openLHAIndex(is) ) Throw<InterfaceException>()
    << "Could not open the LHAPDF index file. The PDF set and member is "
    << "left unchanged." << Exception::warning;
  int set = 0;
  int mem = 0;
  string file;
  double xmin = 0.0;
  double xmax = 0.0;
  double q2min = 0.0;
  double q2max = 0.0;
  int pdftyp = 0;
  int pdfgup = 0;
  int pdfsup = 0;

  while ( indexLine(is, set, mem, file, pdftyp, pdfgup, pdfsup,
		    xmin, xmax, q2min, q2max) ) {
    if ( n == set ) {
      thePDFName = file;
      theMember = mem;
      return;
    }
  }
}

int LHAPDF::getPDFNumber() const {
  ifstream is;
  if ( !openLHAIndex(is) ) Throw<InterfaceException>()
    << "Could not open the LHAPDF index file. The PDF set and member is "
    << "left unchanged." << Exception::warning;
  int set = 0;
  int mem = 0;
  string file;
  double xmin = 0.0;
  double xmax = 0.0;
  double q2min = 0.0;
  double q2max = 0.0;
  int pdftyp = 0;
  int pdfgup = 0;
  int pdfsup = 0;

  while ( indexLine(is, set, mem, file, pdftyp, pdfgup, pdfsup,
		    xmin, xmax, q2min, q2max) )
    if ( thePDFName == file && theMember >= mem ) return set;
  return 0;
}

void LHAPDF::setPDFLIBNumbers(int group, int num) {
  ifstream is;
  if ( !openLHAIndex(is) ) Throw<InterfaceException>()
    << "Could not open the LHAPDF index file. The PDF set and member is "
    << "left unchanged." << Exception::warning;
  int set = 0;
  int mem = 0;
  string file;
  double xmin = 0.0;
  double xmax = 0.0;
  double q2min = 0.0;
  double q2max = 0.0;
  int pdftyp = 0;
  int pdfgup = 0;
  int pdfsup = 0;

  while ( indexLine(is, set, mem, file, pdftyp, pdfgup, pdfsup,
		    xmin, xmax, q2min, q2max) ) {
    if ( pdftyp == thePType && pdfgup == group && pdfsup == num ) {
      thePDFName = file;
      theMember = mem;
      return;
    }
  }
}

string LHAPDF::setPDFLIBNumbers(string cmd) {
  istringstream is(cmd);
  int pdfgup = 0;
  int pdfsup = 0;
  is >> pdfgup >> pdfsup;
  setPDFLIBNumbers(pdfgup, pdfsup);
  return "";
}

pair<int,int> LHAPDF::getPDFLIBNumbers() const {
  ifstream is;
  if ( !openLHAIndex(is) ) Throw<InterfaceException>()
    << "Could not open the LHAPDF index file. The PDF set and member is "
    << "left unchanged." << Exception::warning;
  int set = 0;
  int mem = 0;
  string file;
  double xmin = 0.0;
  double xmax = 0.0;
  double q2min = 0.0;
  double q2max = 0.0;
  int pdftyp = 0;
  int pdfgup = 0;
  int pdfsup = 0;

  while ( indexLine(is, set, mem, file, pdftyp, pdfgup, pdfsup,
		    xmin, xmax, q2min, q2max) )
    if ( thePDFName == file && theMember >= mem )
      return make_pair(pdfgup, pdfsup);
  return make_pair(0, 0);
}

void LHAPDF::setPDFMember(int n) {
  theMember = n;
}

void LHAPDF::setPDFName(string name) {
  thePDFName = name;
}

string LHAPDF::doTest(string input) {
  double x = 0;
  Energy2 Q2 = ZERO;
  Energy2 P2 = ZERO;
  istringstream is(input);
  is >> x >> iunit(Q2, GeV2) >> iunit(P2, GeV2);
  checkUpdate(x, Q2, P2);
  ostringstream os;
  for ( int i = 0; i < 13; ++i ) os << " " << lastXF[i];
  return os.str();
}  

void LHAPDF::checkUpdate(double x, Energy2 Q2, Energy2 P2) const {
  TmpMaskFpuDenorm fpumask;
  checkInit();
  if ( x == lastX && Q2 == lastQ2 && P2 == lastP2 ) return;
  lastX = x;
  lastQ2 = Q2;
  lastP2 = P2;
  vector<F77ThePEGDouble> res(13);

  if ( x < xMin || x > xMax || Q2 < Q2Min || Q2 > Q2Max ) {
      switch ( rangeException ) {
      case rangeThrow: Throw<Exception>()
	<< "Momentum fraction (x=" << x << ") or scale (Q2=" << double(Q2/GeV2)
	<< " GeV^2) was outside of limits in PDF " << name() << "."
	<< Exception::eventerror;
      case rangeZero:
	lastXF = vector<double>(res.begin(), res.end());
	return;
      case rangeFreeze:
	lastX = x = min(max(x, xMin), xMax);
	lastQ2 = Q2 = min(max(Q2, Q2Min), Q2Max);
      }
  } 

  F77ThePEGInteger iset = nset + 1;
  F77ThePEGDouble Q = sqrt(Q2/GeV2);
  F77ThePEGDouble xx = x;
  if ( ptype() == photonType ) {
    F77ThePEGDouble PP2 = P2/GeV2;
    F77ThePEGInteger IP2 = thePhotonOption;
    evolvepdfpm_(iset, xx, Q, PP2, IP2, &res[0]);
  } else {
    if(!enablePartonicGamma) 
      evolvepdfm_(iset, xx, Q, &res[0]);
    else {
      double gamma(0.);
      evolvepdfphotonm_(iset, xx, Q, &res[0],gamma);
      res.push_back(gamma);
    }
  }
  lastXF = vector<double>(res.begin(), res.end());
}

bool LHAPDF::canHandleParticle(tcPDPtr particle) const {
  using namespace ParticleID;
  switch ( ptype() ) {
  case nucleonType:
    return abs(particle->id()) == pplus || abs(particle->id()) == n0;
  case pionType:
    return particle->id() == pi0 || particle->id() == ParticleID::gamma;
  case photonType:
    return particle->id() == ParticleID::gamma;
  default:
    return false;
  }
}

cPDVector LHAPDF::partons(tcPDPtr particle) const {
  // We assume that all standard partons can be extracted.
  cPDVector ret;
  if ( canHandleParticle(particle) ) {
    ret.push_back(getParticleData(ParticleID::g));
    for ( int i = 1, N = getMaxFlav(); i <= N; ++i ) {
      ret.push_back(getParticleData(i));
      ret.push_back(getParticleData(-i));
    }
    // special if needed add photon
    if(enablePartonicGamma) 
      ret.push_back(getParticleData(ParticleID::gamma));
  }
  return ret;
}

namespace LHAPDFIndex {
enum VectorIndices {
  topb = 0, botb = 1, chab = 2, strb = 3, upb = 4, dowb = 5, glu = 6, dow = 7,
  up = 8, str = 9, cha = 10, bot = 11, top = 12, photon = 13 };
}

double LHAPDF::xfx(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
                      double x, double, Energy2 particleScale) const {
  // Here we should return the actual density.
  using namespace ThePEG::ParticleID;
  using namespace LHAPDFIndex;
  checkUpdate(x, partonScale, particleScale);
  switch ( parton->id() ) {
  case t:
    return maxFlav() < 6? 0.0: lastXF[top];
  case tbar:
    return maxFlav() < 6? 0.0: lastXF[topb];
  case b:
    return maxFlav() < 5? 0.0: lastXF[bot];
  case bbar:
    return maxFlav() < 5? 0.0: lastXF[botb];
  case c:
    return maxFlav() < 4? 0.0: lastXF[cha];
  case cbar:
    return maxFlav() < 4? 0.0: lastXF[chab];
  case ParticleID::s:
    return lastXF[str];
  case sbar:
    return lastXF[strb];
  case u:
    switch ( particle->id() ) {
    case n0: return lastXF[dow];
    case pbarminus: return lastXF[upb];
    case nbar0: return lastXF[dowb];
    case pplus:
    default: return lastXF[up];
    }
  case ubar:
    switch ( particle->id() ) {
    case n0: return lastXF[dowb];
    case pbarminus: return lastXF[up];
    case nbar0: return lastXF[dow];
    case pplus:
    default: return lastXF[upb];
    }
  case d:
    switch ( particle->id() ) {
    case n0: return lastXF[up];
    case pbarminus: return lastXF[dowb];
    case nbar0: return lastXF[upb];
    case pplus:
    default: return lastXF[dow];
    }
  case dbar:
    switch ( particle->id() ) {
    case n0: return lastXF[upb];
    case pbarminus: return lastXF[dow];
    case nbar0: return lastXF[up];
    case pplus:
    default: return lastXF[dowb];
    }
  case ParticleID::g:
    return lastXF[glu];
  case ParticleID::gamma:
    return enablePartonicGamma ? lastXF[photon] : 0.;
  }
  return 0.0;
}

double LHAPDF::xfvl(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
		     double l, Energy2 particleScale) const {
  using Math::exp1m;
  return xfvx(particle, parton, partonScale, exp(-l), exp1m(-l), particleScale);
}

double LHAPDF::xfvx(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
		     double x, double, Energy2 particleScale) const {
  // Here we should return the actual valence density. This will only
  // work properly for nucleons
  using namespace ThePEG::ParticleID;
  using namespace LHAPDFIndex;
  checkUpdate(x, partonScale, particleScale);

  switch ( parton->id() ) {
  case t:
  case tbar:
  case b:
  case bbar:
  case c:
  case cbar:
  case ParticleID::s:
  case sbar:
  case ParticleID::gamma:
    return 0.0;
  case u:
    switch ( particle->id() ) {
    case n0: return lastXF[dow] - lastXF[dowb];
    case pbarminus: return 0.0;
    case nbar0: return 0.0;
    case pplus: return lastXF[up] - lastXF[upb];
    default: return 0.0;
    }
  case ubar:
    switch ( particle->id() ) {
    case n0: return 0.0;
    case pbarminus: return lastXF[up] - lastXF[upb];
    case nbar0: return lastXF[dow] - lastXF[dowb];
    case pplus:
    default: return 0.0;
    }
  case d:
    switch ( particle->id() ) {
    case n0: return lastXF[up] - lastXF[upb];
    case pbarminus: return 0.0;
    case nbar0: return 0.0;
    case pplus: return lastXF[dow] - lastXF[dowb];
    default: return 0.0;
    }
  case dbar:
    switch ( particle->id() ) {
    case n0: return 0.0;
    case pbarminus: return lastXF[dow] - lastXF[dowb];
    case nbar0: return lastXF[up] - lastXF[upb];
    case pplus:
    default: return 0.0;
    }
  case ParticleID::g:
    return 0.0;
  }
  return 0.0;
}

double LHAPDF::xfsx(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
		    double x, double, Energy2 particleScale) const {
  // Here we should return the actual density.
  using namespace ThePEG::ParticleID;
  using namespace LHAPDFIndex;
  checkUpdate(x, partonScale, particleScale);

  switch ( parton->id() ) {
  case t:
    return maxFlav() < 6? 0.0: lastXF[top];
  case tbar:
    return maxFlav() < 6? 0.0: lastXF[topb];
  case b:
    return maxFlav() < 5? 0.0: lastXF[bot];
  case bbar:
    return maxFlav() < 5? 0.0: lastXF[botb];
  case c:
    return maxFlav() < 4? 0.0: lastXF[cha];
  case cbar:
    return maxFlav() < 4? 0.0: lastXF[chab];
  case ParticleID::s:
    return lastXF[str];
  case sbar:
    return lastXF[strb];
  case u:
    switch ( particle->id() ) {
    case n0: return lastXF[dowb];
    case pbarminus: return lastXF[upb];
    case nbar0: return lastXF[dowb];
    case pplus: return lastXF[upb];
    default: return lastXF[up];
    }
  case ubar:
    switch ( particle->id() ) {
    case n0: return lastXF[dowb];
    case pbarminus: return lastXF[upb];
    case nbar0: return lastXF[dowb];
    case pplus: return lastXF[upb];
    default: return lastXF[upb];
    }
  case d:
    switch ( particle->id() ) {
    case n0: return lastXF[upb];
    case pbarminus: return lastXF[dowb];
    case nbar0: return lastXF[upb];
    case pplus: return lastXF[dowb];
    default: return lastXF[dow];
    }
  case dbar:
    switch ( particle->id() ) {
    case n0: return lastXF[upb];
    case pbarminus: return lastXF[dowb];
    case nbar0: return lastXF[upb];
    case pplus: return lastXF[dowb];
    default: return lastXF[dowb];
    }
  case ParticleID::g:
    return lastXF[glu];
  case ParticleID::gamma:
    return enablePartonicGamma ? lastXF[photon] : 0.;
  }
  return 0.0;
}

void LHAPDF::persistentOutput(PersistentOStream & os) const {
  os << oenum(thePType) << thePDFName << theMember << thePhotonOption 
     << enablePartonicGamma << theVerboseLevel << theMaxFlav
     << xMin << xMax << ounit(Q2Min, GeV2) << ounit(Q2Max, GeV2);
}

void LHAPDF::persistentInput(PersistentIStream & is, int) {
  is >> ienum(thePType) >> thePDFName >> theMember >> thePhotonOption 
     >> enablePartonicGamma >> theVerboseLevel >> theMaxFlav
     >> xMin >> xMax >> iunit(Q2Min, GeV2) >> iunit(Q2Max, GeV2);
  nset = -1;
  lastReset();
}

ClassDescription<LHAPDF> LHAPDF::initLHAPDF;
// Definition of the static class description member.

int LHAPDF::MaxNSet = 3;

int LHAPDF::lastNSet = 0;

vector<string> LHAPDF::lastNames = vector<string>(3);

vector<int> LHAPDF::lastMem = vector<int>(3);

void LHAPDF::Init() {

  static ClassDocumentation<LHAPDF> documentation
    ("The LHAPDF class inherits from PDFBase and implements an interface "
     "to the LHAPDF library of parton density function parameterizations. "
     "This class is available even if LHAPDF was not properly installed "
     "when ThePEG was installed, but will then produce an error in the "
     "initialization. Note that the valence densities from the xfvx() and "
     "xfvl() function will only work properly for nucleons. All other "
     "particles will have zero valence densities.");

  static Switch<LHAPDF,PType> interfacePType
    ("PType",
     "The type of incoming particles which can be handled by this PDF.",
     &LHAPDF::thePType, nucleonType, true, false);
  static SwitchOption interfacePTypeNucleon
    (interfacePType,
     "Nucleon",
     "Nucleon densities.",
     nucleonType);
  static SwitchOption interfacePTypePionOrVMD
    (interfacePType,
     "PionOrVMD",
     "Pion densities (can also be used for VMD photons).",
     pionType);
  static SwitchOption interfacePTypePhoton
    (interfacePType,
     "Photon",
     "Photon densities.",
     photonType);
  interfacePType.setHasDefault(false);

  static Parameter<LHAPDF,string> interfacePDFName
    ("PDFName",
     "The name if the PDF set to be used. Should be the full name including "
     "the <code>.LHpdf</code> or <code>.LHgrid</code> suffix.",
     &LHAPDF::thePDFName, "cteq6ll.LHpdf", true, false,
     &LHAPDF::setPDFName);
  interfacePDFName.setHasDefault(false);

  static Parameter<LHAPDF,int> interfacePDFNumber
    ("PDFNumber",
     "The number of the PDF set and member to be used.",
     0, 10042, 1, 0,
     true, false, Interface::lowerlim,
     &LHAPDF::setPDFNumber, &LHAPDF::getPDFNumber,
     (int(LHAPDF::*)()const)(0), (int(LHAPDF::*)()const)(0),
     (int(LHAPDF::*)()const)(0));
  interfacePDFNumber.setHasDefault(false);

  static Parameter<LHAPDF,int> interfaceMember
    ("Member",
     "The chosen member of the selected PDF set.",
     &LHAPDF::theMember, 0, 0, Constants::MaxInt,
     true, false, Interface::limited,
     &LHAPDF::setPDFMember, (int(LHAPDF::*)()const)(0),
     (int(LHAPDF::*)()const)(0), &LHAPDF::getMaxMember,
     (int(LHAPDF::*)()const)(0));
  interfaceMember.setHasDefault(false);

  static Command<LHAPDF> interfacePDFLIBNumbers
    ("PDFLIBNumbers",
     "Set the PDF set and member to be used by specifying the old PDFLIB "
     "group and set number.",
     &LHAPDF::setPDFLIBNumbers, true);

  static Switch<LHAPDF,bool> interfaceEnablePartonicGamma
    ("EnablePartonicGamma",
     "Enable the option of having photons as partons inside a hadron",
     &LHAPDF::enablePartonicGamma, false, false, false);
  static SwitchOption interfaceEnablePartonicGammaYes
    (interfaceEnablePartonicGamma,
     "Yes",
     "Include partonic photons",
     true);
  static SwitchOption interfaceEnablePartonicGammaNo
    (interfaceEnablePartonicGamma,
     "No",
     "Don't include them",
     false);

  static Switch<LHAPDF,int> interfacePhotonOption
    ("PhotonOption",
     "Different options for handling off-shell photon distributions.",
     &LHAPDF::thePhotonOption, 7, true, false);
  static SwitchOption interfacePhotonOptionDipoleDampening
    (interfacePhotonOption,
     "DipoleDampening",
     "Dipole dampening by integration (very time consuming).",
     1);
  static SwitchOption interfacePhotonOptionMaxScales
    (interfacePhotonOption,
     "MaxScales",
     "\\f$P_0^2=\\max(Q_0^2,P^2)\\f$",
     2);
  static SwitchOption interfacePhotonOptionAddScales
    (interfacePhotonOption,
     "AddScales",
     "\\f$P_0^{'2}=Q_0^2+p^2\\f$",
     3);
  static SwitchOption interfacePhotonOptionPeffPreserve
    (interfacePhotonOption,
     "PeffPreserve",
     "\\f$P_{eff}\\f$ preserving momentum sum.",
     4);
  static SwitchOption interfacePhotonOptionPintPreserve
    (interfacePhotonOption,
     "PintPreserve",
     "\\f$P_{int}\\f$ preserving momentum sum and average evolution range.",
     5);
  static SwitchOption interfacePhotonOptionPeffMatch
    (interfacePhotonOption,
     "PeffMatch",
     "\\f$P_{eff}\\f$ matched to \\f$P_0\\f$ in \\f$P^2\\rightarrow Q^2\\f$ "
     "limit.",
     6);
  static SwitchOption interfacePhotonOptionPintMatch
    (interfacePhotonOption,
     "PintMatch",
     "\\f$P_{int}\\f$ matched to \\f$P_0\\f$ in \\f$P^2\\rightarrow Q^2\\f$ "
     "limit.",
     7);

  static Parameter<LHAPDF,int> interfaceMaxNSet
    ("MaxNSet",
     "The maximum number of simultaneous pdfs that can be used in LHAPDF. "
     "Should be set to the parameter <code>nmxset</code> in the "
     "<code>parmsetup.inc</code> file compiled into the installed LHAPDF "
     "library you are using (by default this is set to 3)",
     0, 3, 1, 0,
     true, false, Interface::lowerlim,
     &LHAPDF::setMaxNSet, &LHAPDF::getMaxNSet,
     (int(LHAPDF::*)()const)(0), (int(LHAPDF::*)()const)(0),
     (int(LHAPDF::*)()const)(0));


  static Command<LHAPDF> interfaceTest
    ("Test",
     "Write out the values of the chosen PDF set using the x, Q2 and P2 "
     "parameters supplied.",
     &LHAPDF::doTest, true);


  static Switch<LHAPDF,int> interfaceVerboseLevel
    ("VerboseLevel",
     "The verbosity of the output from the LHAPDF library.",
     &LHAPDF::theVerboseLevel, 0, true, false);
  static SwitchOption interfaceVerboseLevelSilent
    (interfaceVerboseLevel,
     "Silent",
     "Trying to inhibit all output from the LHAPDF library "
     "(unfortunately not always possible).",
     0);
  static SwitchOption interfaceVerboseLevelNormal
    (interfaceVerboseLevel,
     "Normal",
     "Normal output from the LHAPDF library "
     "(unfortunately to the standard output).",
     1);
  interfaceVerboseLevel.setHasDefault(false);

  static Parameter<LHAPDF,int> interfaceMaxFlav
    ("MaxFlav",
     "The maximum number of flavours for which non-zero densities are "
     "reported. The actual number of flavours may be less depending on "
     "the chosen PDF set.",
     &LHAPDF::theMaxFlav, 5, 3, 0,
     true, false, Interface::lowerlim);

  interfacePType.rank(10);
  interfacePDFName.rank(9);
  interfaceMember.rank(8);

}

