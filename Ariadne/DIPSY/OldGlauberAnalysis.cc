// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the OldGlauberAnalysis class.
//

#include "OldGlauberAnalysis.h"
#include "DipoleXSec.h"
#include "DipoleState.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"

#include "gsl/gsl_sf_bessel.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace DIPSY;

OldGlauberAnalysis::OldGlauberAnalysis(): nnSigTot(96.0*millibarn), nnSigEl(24.0*millibarn),
				    nnSigInND(30.0*millibarn) {}

OldGlauberAnalysis::~OldGlauberAnalysis() {}

void OldGlauberAnalysis::initialize() {
  generator()->histogramFactory()->initrun();
  generator()->histogramFactory()->registerClient(this);

  sumNTot = 0.0;
  sumNColl = vector<double>(nstr, 0.0);
  sumWeights = sumsig = sumsigel = sumsignd = sumsigdt = sumcoll = sumcoll2 = sumpart = sumpart2
    = vector<CrossSection>(nstr, ZERO);
  sumsig2 = sumsigel2 = sumsignd2 =  sumsigdt2 = vector<CrossSection2>(nstr, ZERO);
  sum = ZERO;
  sum2 = ZERO;
  sumND = ZERO;
  sumND2 = ZERO;
}

void OldGlauberAnalysis::
analyze(const DipoleState & dr, const DipoleState & dl,
	const ImpactParameters & b, const DipoleXSec & xsec,
	double fsum, CrossSection weight) {

  if ( dl.initialDipoles().size()%3 ) return;
  if ( dr.initialDipoles().size()%3 ) return;
  int Nl = dl.initialDipoles().size()/3;
  int Nr = dr.initialDipoles().size()/3;

  sum += 2.0*xsec.unitarize(fsum)*weight;
  sum2 += sqr(2.0*xsec.unitarize(fsum)*weight);
  sumND += xsec.unitarize(2.0*fsum)*weight;
  sumND2 += sqr(xsec.unitarize(2.0*fsum)*weight);

  int Nl = dl.initialDipoles().size()/3;
  int Nr = dr.initialDipoles().size()/3;
  vector< vector<int> > lpart(2*nstr, vector<int>(Nl));
  vector< vector<int> > rpart(2*nstr, vector<int>(Nr));
  vector<int> ncoll(2*nstr);
  vector<int> coll(2*nstr);
  double alpha = 2.0*(1.0 - nnSigInND/nnSigTot);

  for ( int il = 0; il < Nl; ++il ) {
    Parton::Point pl = (dl.initialDipoles()[il*3]->partons().first->position() +
		       dl.initialDipoles()[il*3 + 1]->partons().first->position() +
		       dl.initialDipoles()[il*3 + 2]->partons().first->position())/3.0;
    for ( int ir = 0; ir < Nr; ++ir ) {
      Parton::Point pr = (dr.initialDipoles()[ir*3]->partons().first->position() +
		  dr.initialDipoles()[ir*3 + 1]->partons().first->position() +
		  dr.initialDipoles()[ir*3 + 2]->partons().first->position())/3.0;

      Length r = b.difference(pl, pr).pt()*hbarc;
      

      // Black disk
      Length R = sqrt((nnSigTot - nnSigEl)/Constants::pi);
      if ( r <= R ) {
	++lpart[0][il];
	++rpart[0][ir];
	++ncoll[0];
	++coll[0];
	++lpart[nstr][il];
	++rpart[nstr][ir];
	++ncoll[nstr];
	++coll[nstr];
      }

      // Black disk new style
      R = sqrt(nnSigTot/(2.0*Constants::pi));
      if ( r <= R ) {
	++lpart[1][il];
	++rpart[1][ir];
	++ncoll[1];
	++coll[1];
	++lpart[nstr + 1][il];
	++rpart[nstr + 1][ir];
	++ncoll[nstr + 1];
	++coll[nstr + 1];
      }

      // Gray disk Andras style
      R = sqrt(nnSigTot/Constants::pi);
      double a = (nnSigTot - nnSigEl)/nnSigTot;
      if ( r <= R && UseRandom::rndbool(a) ) {
	++lpart[2][il];
	++rpart[2][ir];
	++ncoll[2];
	++coll[2];
      }
      if ( r <= R && UseRandom::rndbool(a) ) {
	++lpart[nstr + 2][il];
	++rpart[nstr + 2][ir];
	++ncoll[nstr + 2];
	++coll[nstr + 2];
      }

      // Gray disk New style
      R = nnSigTot/sqrt(4.0*Constants::pi*nnSigEl);
      a = 2.0*nnSigEl/nnSigTot;
      if ( r <= R && UseRandom::rndbool(a) ) {
	++lpart[3][il];
	++rpart[3][ir];
	++ncoll[3];
	++coll[3];
      }
      if ( r <= R && UseRandom::rndbool(a) ) {
	++lpart[nstr + 3][il];
	++rpart[nstr + 3][ir];
	++ncoll[nstr + 3];
	++coll[nstr + 3];
      }

      // Gussian new style
      R = nnSigTot/sqrt(8.0*Constants::pi*nnSigEl);
      a = 4.0*nnSigEl/nnSigTot;
      if ( UseRandom::rndbool(a*exp(-sqr(r/R))) ) {
	++lpart[4][il];
	++rpart[4][ir];
	++ncoll[4];
	++coll[4];
      }
      if ( UseRandom::rndbool(a*exp(-sqr(r/R))) ) {
	++lpart[nstr + 4][il];
	++rpart[nstr + 4][ir];
	++ncoll[nstr + 4];
	++coll[nstr + 4];
      }

      // Gray disk NEW style
      R = nnSigTot/sqrt(4.0*Constants::pi*nnSigEl);
      a = nnSigEl/(nnSigTot - nnSigInND);
      if ( r <= R && UseRandom::rndbool(a) ) {
	++coll[5];
	if ( UseRandom::rndbool(alpha) ) {
	  ++lpart[5][il];
	  ++rpart[5][ir];
	  ++ncoll[5];
	}
      }
      if ( r <= R && UseRandom::rndbool(a) ) {
	++coll[nstr + 5];
	if ( UseRandom::rndbool(alpha) ) {
	  ++lpart[nstr + 5][il];
	  ++rpart[nstr + 5][ir];
	  ++ncoll[nstr + 5];
	}
      }
      // Gray disk NEW style (inel)
      if ( r <= R && UseRandom::rndbool(a) ) {
	++coll[6];
	if ( UseRandom::rndbool(2.0*alpha - sqr(alpha)) ) {
	  ++lpart[6][il];
	  ++rpart[6][ir];
	  ++ncoll[6];
	}
      }
      if ( r <= R && UseRandom::rndbool(a) ) {
	++coll[nstr + 6];
	if ( UseRandom::rndbool(2.0*alpha - sqr(alpha)) ) {
	  ++lpart[nstr + 6][il];
	  ++rpart[nstr + 6][ir];
	  ++ncoll[nstr + 6];
	}
      }
    }
  }

  CrossSection w = sqr(hbarc)*b.weight();
  sumNTot += 1.0;
  for ( int istr = 0; istr < nstr; ++istr ) {
    if ( !coll[istr] ) continue;
    double T = 1.0;
    if ( istr >= 5 ) T = 1.0 - pow(1.0 - alpha, coll[istr]);
    if ( ncoll[istr] ) sumWeights[istr] += w;
    sumNColl[istr] += 1.0;
    sumsig[istr] += T*w;
    sumsig2[istr] += sqr(T*w);
    double T2 = 1.0;

    sumsigdt[istr] += sqr(T)*w;
    sumsigdt2[istr] += sqr(sqr(T)*w);
    if ( istr >= 5 ) T2 = 1.0 - pow(1.0 - alpha, coll[istr + nstr]);
    if ( coll[istr + nstr] ) {
      sumsigel[istr] += T*T2*w;
      sumsigel2[istr] += sqr(T*T2*w);
    }

    sumsignd[istr] += (2.0*T-sqr(T))*w;
    sumsignd2[istr] += sqr((2.0*T-sqr(T))*w);
    double nlpart = 0;
    double nrpart = 0;
    for ( int il = 0; il < Nl; ++il ) if ( lpart[istr][il] ) nlpart += 1.0;
    for ( int ir = 0; ir < Nr; ++ir ) if ( rpart[istr][ir] ) nrpart += 1.0;
    sumcoll[istr] += ncoll[istr]*w;
    sumcoll2[istr] += sqr(ncoll[istr])*w;
    sumpart[istr] += (nlpart + nrpart)*w;
    sumpart2[istr] += sqr(nlpart + nrpart)*w;
  } 
}

template <typename Unit>
double errat(Unit a, Unit b, Unit da, Unit db ) {
  return sqrt(sqr(da/a) + sqr(db/b))*a/b;
}

void OldGlauberAnalysis::finalize(long neve) {
  if ( sumNTot <= 0 ) return;

  sum /= neve;
  sum2 /= neve;
  CrossSection err = sqrt((sum2 - sqr(sum))/neve);
  stub(": Total cross section (DIPSY): ") << ouniterr(sum, err, nanobarn) << " nb." << endl;

  sumND /= neve;
  sumND2 /= neve;
  CrossSection errND = sqrt((sumND2 - sqr(sumND))/neve);
  stub(": Inelastic ND xsec (DIPSY): ") << ouniterr(sumND, errND, nanobarn) << " nb" << endl;

  printstrat(1);
  printstrat(3);
  printstrat(5);
  printstrat(4);
  printstrat(0);
  printstrat(2);
  printstrat(6);

}

void OldGlauberAnalysis::printstrat(int istr) const {
  string strat = " (old black disk): ";
  if ( istr == 1 ) strat = " (black disc): ";
  if ( istr == 2 ) strat = " (old grey disc): ";
  if ( istr == 3 ) strat = " (grey disc): ";
  if ( istr == 4 ) strat = " (Gaussian): ";
  if ( istr == 5 ) strat = " (grey3 disc): ";
  if ( istr == 6 ) strat = " (Inel grey disc): ";
  CrossSection xsectot = 2.0*sumsig[istr]/sumNTot;
  CrossSection xsecerr = 2.0*sqrt((sumsig2[istr]/sumNTot - sqr(sumsig[istr]/sumNTot))/sumNTot);
  stub(": Total cross section" + strat)
    << ouniterr(xsectot, xsecerr, nanobarn) << " nb" << endl;
  CrossSection xsecel = sumsigel[istr]/sumNTot;
  CrossSection xsecelerr = sqrt((sumsigel2[istr]/sumNTot - sqr(sumsigel[istr]/sumNTot))/sumNTot);
  stub(": Elastic xsec" + strat) << ouniterr(xsecel, xsecelerr, nanobarn) << " nb" << endl;
  CrossSection xsecdt = sumsigdt[istr]/sumNTot - xsecel;
  CrossSection xsecdterr =
    sqrt((sumsigdt2[istr]/sumNTot - sqr(sumsigdt[istr]/sumNTot))/sumNTot + sqr(xsecelerr));
  stub(": Diffractive xsec" + strat) << ouniterr(xsecdt, xsecdterr, nanobarn) << " nb" << endl;
  CrossSection xsecnd = sumsignd[istr]/sumNTot;
  CrossSection xsecnderr = sqrt((sumsignd2[istr]/sumNTot - sqr(sumsignd[istr]/sumNTot))/sumNTot);
  stub(": Inelastic ND xsec" + strat)
    << ouniterr(xsecnd, xsecnderr, nanobarn) << " nb" << endl;
  double avpart = sumpart[istr]/sumWeights[istr];
  double av2part = sumpart2[istr]/sumWeights[istr];
  double errpart = sqrt((av2part - sqr(avpart))/sumNTot);
  stub(": #participants" + strat) << ouniterr(avpart, errpart, 1.0) << endl;
  double avcoll = sumcoll[istr]/sumWeights[istr];
  double av2coll = sumcoll2[istr]/sumWeights[istr];
  double errcoll = sqrt((av2coll - sqr(avcoll))/sumNTot);
  stub(": #collisions" + strat) << ouniterr(avcoll, errcoll, 1.0) << endl;
}

IBPtr OldGlauberAnalysis::clone() const {
  return new_ptr(*this);
}

IBPtr OldGlauberAnalysis::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void OldGlauberAnalysis::persistentOutput(PersistentOStream & os) const {
  const CrossSection nb = nanobarn;
  const CrossSection2 nb2= sqr(nb);
  os << ounit(nnSigTot, nb) << ounit(nnSigEl, nb) << ounit(nnSigInND, nb)
     << sumNTot << sumNColl << ounit(sumWeights, nb) << ounit(sumsig, nb)
     << ounit(sumsigel, nb) << ounit(sumsignd, nb) << ounit(sumsigdt, nb)
     << ounit(sumsig2, nb2) << ounit(sumsigel2, nb2) << ounit(sumsignd2, nb2)
     << ounit(sumsigdt2, nb2) << ounit(sumcoll, nb) << ounit(sumcoll2, nb)
     << ounit(sumpart, nb) << ounit(sumpart2, nb) << ounit(sum, nb)
     << ounit(sum2, nb2) << ounit(sumND, nb) << ounit(sumND2, nb2);
}

void OldGlauberAnalysis::persistentInput(PersistentIStream & is, int) {
  const CrossSection nb = nanobarn;
  const CrossSection2 nb2= sqr(nb);
  is >> iunit(nnSigTot, nb) >> iunit(nnSigEl, nb) >> iunit(nnSigInND, nb)
     >> sumNTot >> sumNColl >> iunit(sumWeights, nb) >> iunit(sumsig, nb)
     >> iunit(sumsigel, nb) >> iunit(sumsignd, nb) >> iunit(sumsigdt, nb)
     >> iunit(sumsig2, nb2) >> iunit(sumsigel2, nb2) >> iunit(sumsignd2, nb2)
     >> iunit(sumsigdt2, nb2) >> iunit(sumcoll, nb) >> iunit(sumcoll2, nb)
     >> iunit(sumpart, nb) >> iunit(sumpart2, nb) >> iunit(sum, nb)
     >> iunit(sum2, nb2) >> iunit(sumND, nb) >> iunit(sumND2, nb2);
}




// Static variable needed for the type description system in ThePEG.
#include "ThePEG/Utilities/DescribeClass.h"
DescribeClass<OldGlauberAnalysis,DIPSY::DipoleAnalysisHandler>
  describeDIPSYOldGlauberAnalysis("DIPSY::OldGlauberAnalysis",
				   "OldGlauberAnalysis.so");

void OldGlauberAnalysis::Init() {

  static ClassDocumentation<OldGlauberAnalysis> documentation
    ("There is no documentation for the OldGlauberAnalysis class");

  static Parameter<OldGlauberAnalysis,CrossSection> interfaceTotalnnXSec
    ("TotalnnXSec",
     "The assumed total nucleon-nucleon cross section assumed in the calculation (in mb).",
     &OldGlauberAnalysis::nnSigTot, millibarn, 96.0*millibarn, 0.0*millibarn, 0*millibarn,
     true, false, Interface::lowerlim);

  static Parameter<OldGlauberAnalysis,CrossSection> interfaceElasticnnXSec
    ("ElasticnnXSec",
     "The assumed elastic nucleon-nucleon cross section assumed in the calculation (in mb).",
     &OldGlauberAnalysis::nnSigEl, millibarn, 24.0*millibarn, 0.0*millibarn, 0*millibarn,
     true, false, Interface::lowerlim);

  static Parameter<OldGlauberAnalysis,CrossSection> interfaceElasticnnXSecInND
    ("InElasticnnXSec",
     "The assumed inelastic, non-diffractive nucleon-nucleon cross section "
     "assumed in the calculation (in mb).",
     &OldGlauberAnalysis::nnSigInND, millibarn, 30.0*millibarn, 0.0*millibarn, 0*millibarn,
     true, false, Interface::lowerlim);

}

