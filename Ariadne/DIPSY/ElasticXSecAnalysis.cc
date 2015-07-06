// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ElasticXSecAnalysis class.
//

#include "ElasticXSecAnalysis.h"
#include "DipoleXSec.h"
#include "DipoleState.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"

#include "gsl/gsl_sf_bessel.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace DIPSY;

ElasticXSecAnalysis::ElasticXSecAnalysis(): nb(100), db(0.1*InvGeV), nq(200), dq(0.05*GeV) {}

ElasticXSecAnalysis::~ElasticXSecAnalysis() {}

void ElasticXSecAnalysis::initialize() {
  generator()->histogramFactory()->initrun();
  generator()->histogramFactory()->registerClient(this);
  sumTw = vector<CrossSection>(nb + 1, 0.0*picobarn);
  sumT2w = vector<CrossSection>(nb + 1, 0.0*picobarn);
  sumT4w = vector<CrossSection>(nb + 1, 0.0*picobarn);
  sumElAw = vector<CrossSection>(nb + 1, 0.0*picobarn);
  sumElA2w = vector<CrossSection>(nb + 1, 0.0*picobarn);
  sumw = vector<CrossSection>(nb + 1, 0.0*picobarn);
  sumElw = vector<CrossSection>(nb + 1, 0.0*picobarn);
  sumn = vector<long>(nb + 1, 0);
  bmap.clear();

  sumDEw = vector<CrossSection>(nq + 1, 0.0*picobarn);
  sumDE2w = vector<CrossSection>(nq + 1, 0.0*picobarn);

}

void ElasticXSecAnalysis::
analyze(const DipoleState & dl, const DipoleState & dr,
	const ImpactParameters & b, const DipoleXSec & xsec,
	double fsum, CrossSection weight) {
  int ib = min(nb, int(b.bVec().pt()/db));
  sumTw[ib] += xsec.unitarize(fsum)*weight;
  sumT2w[ib] += sqr(xsec.unitarize(fsum))*weight;
  sumT4w[ib] += sqr(sqr(xsec.unitarize(fsum)))*weight;
  sumElAw[ib] += xsec.unitarize(fsum)*weight;
  sumElA2w[ib] += sqr(xsec.unitarize(fsum))*dl.weight()*dr.weight()*weight;
  if (dl.weight()*dr.weight() != 0.0  )
    sumElw[ib] += weight/(dl.weight()*dr.weight());
  sumw[ib] += weight;

  for (int iq = 0; iq < nq; iq++) {
    Energy q = iq*dq;
    double ampb = xsec.unitarize(fsum);
    double ampq = ampb*gsl_sf_bessel_J0(b.bVec().pt()*q);
    sumDEw[iq] += ampq*weight;
    sumDE2w[iq] += sqr(ampq)*weight;
  }

  CrossSection bw = sqr(hbarc)*b.weight();
  double x = xsec.unitarize(fsum)*weight/bw;
  ++sumn[ib];
  bmap[b.bVec().pt()] = make_pair(x, bw);
}

void ElasticXSecAnalysis::finalize(long neve) {
  if ( neve <= 0 ) return;

  CrossSection elXsec = 0.0*picobarn;
  CrossSection diffXsec = 0.0*picobarn;
  CrossSection totXsec = 0.0*picobarn;
  QTY<4,0,0>::Type totErrSqr = 0.0*sqr(picobarn);
  QTY<4,0,0>::Type diffErrSqr = 0.0*sqr(picobarn);
  QTY<4,0,0>::Type elErrSqr = 0.0*sqr(picobarn);

  CrossSection totalWeight = 0.0*picobarn;

  if ( ThePEG_DEBUG_LEVEL )
    generator()->log() << setw(20) << name() + " -----" << " Debug output:" << endl;
  for ( int ib = 0; ib <= nb; ++ib ) {
    if ( sumw[ib] == ZERO ) continue;
    double avT = sumTw[ib]/sumw[ib];
    double avT2 = sumT2w[ib]/sumw[ib];
    double avT4 = sumT4w[ib]/sumw[ib];
    double errT = sqrt(abs(avT2 - sqr(avT))/sumn[ib]);
    if ( sumn[ib] == 1 ) errT = avT;
    double errT2 = sqrt(abs(avT4 - sqr(avT2))/sumn[ib]);
    if ( sumn[ib] == 1 ) errT2 = avT2;

    if ( ThePEG_DEBUG_LEVEL )
      generator()->log() << setw(25) << "b: " << ib*db*GeV << ", avT: " << avT << endl;

    double avElA = sumElAw[ib]/sumElw[ib];
    double avElA2 = sumElA2w[ib]/sumElw[ib];
    double errElA = sqrt(abs(avElA2 - sqr(avElA))/sumn[ib]);
    double errElA2 = 2.0*errElA*avElA;
    if ( sumn[ib] == 1 ) errElA = sqr(avElA);

    elXsec += sumElw[ib]*sqr(avElA);
    elErrSqr += sqr(sumElw[ib]*errElA2);
    diffXsec += sumw[ib]*avT2;
    diffErrSqr += sqr(sumw[ib]*errT2);
    totXsec += sumw[ib]*2.0*avT;
    totErrSqr += sqr(2.0*sumw[ib]*errT);

    totalWeight += sumw[ib];
  }
  if ( ThePEG_DEBUG_LEVEL )
    generator()->log()  << setw(20) <<  "----- " + name() << endl;;

  CrossSection toterr = sqrt(totErrSqr)/double(neve);
  CrossSection differr = sqrt(diffErrSqr)/double(neve);
  CrossSection elerr = sqrt(elErrSqr)/double(neve);
  CrossSection diffexerr = sqrt(sqr(differr) + sqr(elerr));
  CrossSection inelexerr =  sqrt(sqr(differr) + sqr(toterr));
  elXsec /= double(neve);
  diffXsec /= double(neve);
  totXsec /= double(neve);

  generator()->log().setf(ios::left, ios::adjustfield);
  generator()->log()
    << setw(50) << name() + ": Elastic cross section:"
    << ouniterr(elXsec, elerr, nanobarn) << " nb." << endl
    << setw(50) << name() + ": Crosscheck totxsec:"
    << ouniterr(totXsec, toterr, nanobarn) << " nb." << endl
    << setw(50) << name() + ": Diffractive cross section:"
    << ouniterr(diffXsec, differr, nanobarn) << " nb." << endl
    << setw(50) << name() + ": Diffractive excitation cross section:"
    << ouniterr(diffXsec - elXsec, diffexerr, nanobarn) << " nb." << endl
    << setw(50) << name() + ": Non-diffractive inelastic cross section:"
    << ouniterr(totXsec - diffXsec, inelexerr, nanobarn) << " nb." << endl;

  dSigmadq = generator()->histogramFactory()->createHistogram1D
    ("dSigmadq",nq,-dq/GeV,(nq-0.5)*dq/GeV);

  for (int iq = 0; iq < nq; iq++) {
    Energy q = iq*dq;

    CrossSection Aq = sumDEw[iq]/double(neve);

    //this exta factor 1/2pi comes from fourier transform 
    //and how the MC weights are for \int d^2b rather than \int db*b that
    //should go with the bessel funciton
    QTY<3,0,0>::Type sigmaq = q/hbarc*sqr(Aq)/(2.0*Constants::pi);

    dSigmadq->fill(q/GeV, sigmaq/nanobarn*GeV/hbarc);
  }
}

IBPtr ElasticXSecAnalysis::clone() const {
  return new_ptr(*this);
}

IBPtr ElasticXSecAnalysis::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void ElasticXSecAnalysis::persistentOutput(PersistentOStream & os) const {
  os << nb << ounit(db, InvGeV) << nq << ounit(dq, GeV) << sumn;
}

void ElasticXSecAnalysis::persistentInput(PersistentIStream & is, int) {
  is >> nb >> iunit(db, InvGeV) >> nq >> iunit(dq, GeV) >> sumn;
}


// Static variable needed for the type description system in ThePEG.
#include "ThePEG/Utilities/DescribeClass.h"
DescribeClass<ElasticXSecAnalysis,DIPSY::DipoleAnalysisHandler>
  describeDIPSYElasticXSecAnalysis("DIPSY::ElasticXSecAnalysis",
				   "ElasticXSecAnalysis.so");

void ElasticXSecAnalysis::Init() {

  static Parameter<ElasticXSecAnalysis,int> interfaceNBins
    ("NBins",
     "Number of bins in impact parameter.",
     &ElasticXSecAnalysis::nb, 100, 1, 0,
     true, false, Interface::lowerlim);

  static Parameter<ElasticXSecAnalysis,InvEnergy> interfaceDeltaB
    ("DeltaB",
     "Size of intervals in impact parameter.",
     &ElasticXSecAnalysis::db, InvGeV, 0.1*InvGeV, 0.0*InvGeV, 0.0*InvGeV,
     true, false, Interface::lowerlim);

  static Parameter<ElasticXSecAnalysis,int> interfaceNqBins
    ("NqBins",
     "Number of bins in exchanged momentum q.",
     &ElasticXSecAnalysis::nq, 100, 1, 0,
     true, false, Interface::lowerlim);

  static Parameter<ElasticXSecAnalysis,Energy> interfaceDeltaq
    ("Deltaq",
     "Size of intervals in exchanged momentum q.",
     &ElasticXSecAnalysis::dq, GeV, 0.1*GeV, 0.0*GeV, 0.0*GeV,
     true, false, Interface::lowerlim);

  static ClassDocumentation<ElasticXSecAnalysis> documentation
    ("There is no documentation for the ElasticXSecAnalysis class");

}

