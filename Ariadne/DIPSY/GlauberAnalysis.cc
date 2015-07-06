// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the GlauberAnalysis class.
//

#include "GlauberAnalysis.h"
#include "DipoleXSec.h"
#include "DipoleState.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"

#include "gsl/gsl_sf_bessel.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace DIPSY;

template <typename T>
valarray<T> sqr(valarray<T> v) { return v*v; }

GlauberAnalysis::GlauberAnalysis()
  : nnSigTot(50.0*millibarn), nnSigEl(9.2*millibarn), nnSigInND(35.9*millibarn), Nt(1000), tMax(10.0*GeV2) {}

GlauberAnalysis::~GlauberAnalysis() {}

void GlauberAnalysis::initialize() {
  generator()->histogramFactory()->initrun();
  generator()->histogramFactory()->registerClient(this);

  ntot = 0;
  sigtot = signd = sigel = sigdt = sigdr = sigdl = sigdd
    = valarray<double>(ZERO, nstr);
  sigtot2 = signd2 = sigel2 = sigdt2 = sigdr2 = sigdl2 = sigdd2
    = valarray<double>(ZERO, nstr);
  sTwLR = sTw2LR = sTwLR2 = sTwL2R = sTwR2L = valarray<double>(ZERO, nstr);
  sT2wLR = sT2w2LR = sT2wLR2 = sT2wL2R = sT2wR2L = valarray<double>(ZERO, nstr);
  swLR = sw2LR = swLR2 = swL2R = swR2L = 0.0;
  hists = vector<FactoryBase::tH1DPtr>(3*nstr);
  bookHistos();

}

valarray<double> GlauberAnalysis::getT(const DipoleState & dr, const DipoleState & dl,
				      const ImpactParameters & b, const DipoleXSec & xsec,
				      double fsum) const {
  valarray<double> T(0.0, nstr);
  T[0] = xsec.unitarize(fsum);
  int Ndl = dl.initialDipoles().size();
  if ( Ndl%3 && Ndl != 1 ) return T;
  int Ndr = dr.initialDipoles().size();
  if ( Ndr%3 && Ndr != 1 ) return T;
  int Nl = dl.initialDipoles().size()/3;
  if ( Ndl == 1 ) Nl = 1;
  int Nr = dr.initialDipoles().size()/3;
  if ( Ndr == 1 ) Nr = 1;

  double alpha = 2.0*(1.0 - nnSigInND/nnSigTot);
  int nalpha = 0;
    
  for ( int il = 0; il < Nl; ++il ) {
    Parton::Point pl =
      (dl.initialDipoles()[0]->partons().first->position() +
       dl.initialDipoles()[0]->partons().second->position())/2.0;
    if ( Ndl > 1 )
      pl = (dl.initialDipoles()[il*3]->partons().first->position() +
	    dl.initialDipoles()[il*3 + 1]->partons().first->position() +
	    dl.initialDipoles()[il*3 + 2]->partons().first->position())/3.0;

    for ( int ir = 0; ir < Nr; ++ir ) {
      Parton::Point pr =
	(dr.initialDipoles()[0]->partons().first->position() +
	 dr.initialDipoles()[0]->partons().second->position())/2.0;
      if ( Ndr > 1 )
	pr = (dr.initialDipoles()[ir*3]->partons().first->position() +
	      dr.initialDipoles()[ir*3 + 1]->partons().first->position() +
	      dr.initialDipoles()[ir*3 + 2]->partons().first->position())/3.0;

      Length r = b.difference(pl, pr).pt()*hbarc;

      // Black disc
      Length R = sqrt(nnSigTot/(2.0*Constants::pi));
      if ( r <= R ) T[1] = 1.0;
      
      // Gray disc
      R = nnSigTot/sqrt(4.0*Constants::pi*nnSigEl);
      double a = 2.0*nnSigEl/nnSigTot;
      if ( r <= R && UseRandom::rndbool(a) ) T[2] = 1.0;

      // Gray3 disc
      R = nnSigTot/sqrt(4.0*Constants::pi*nnSigEl);
      a = nnSigEl/(nnSigTot - nnSigInND);
      if ( r <= R && UseRandom::rndbool(a) ) ++nalpha;

       // Gussian
      R = nnSigTot/sqrt(8.0*Constants::pi*nnSigEl);
      a = 4.0*nnSigEl/nnSigTot;
      if ( UseRandom::rndbool(a*exp(-sqr(r/R))) ) T[4] = 1.0;

      // Old black Disk
      R = sqrt((nnSigTot - nnSigEl)/Constants::pi);
      if ( r <= R ) T[5] = 1.0;

      // Old gray disk
      R = sqrt(nnSigTot/Constants::pi);
      a = (nnSigTot - nnSigEl)/nnSigTot;
      if ( r <= R && UseRandom::rndbool(a) ) T[6] = 1.0;

      // ND black Disk
      R = sqrt(nnSigInND/Constants::pi);
      if ( r <= R ) T[7] = 1.0;

    }
  }

  T[3] = 1.0 - pow(1.0 - alpha, nalpha);

  return T;
}

void GlauberAnalysis::
analyze(const vector<DipoleStatePtr> & vr, const vector<DipoleStatePtr> & vl,
	const vector<ImpactParameters> & vb, const DipoleXSec & xsec,
	const Vec3D & probs, double jac) {

  int Nr = vr.size();
  int Nl = vl.size();
  if ( Nr*Nl*vb.size() == 0 ) return;

  for ( int ib = 0, Nb = vb.size(); ib < Nb; ++ib ) {
    CrossSection bweight = sqr(hbarc)*vb[ib].weight()*jac;
    double bw = bweight/nanobarn;
    double bw2 = sqr(bw);

    valarray<double> sumLR(0.0, nstr); // T averaged over L and R
    int nLR = 0;
    valarray<double> sum2LR(0.0, nstr); // T^2 average over L and R
    int n2LR = 0;
    valarray<double> sumLR2(0.0, nstr); // (T averaged over L and R )^2
    int nLR2 = 0;
    valarray<double> sumL2R(0.0, nstr); // (T averaged over L)^2 averaged over R
    int nL2R = 0;
    valarray<double> sumR2L(0.0, nstr); // (T averaged over R)^2 averaged over L
    int nR2L = 0;

    vector< vector < valarray<double> > > UT(Nl, vector < valarray<double> >(Nr));
    for ( int ir1 = 0; ir1 < Nr; ++ir1 )
      for ( int il1 = 0; il1 < Nl; ++il1 )
	UT[il1][ir1] = getT(*vl[il1], *vr[ir1], vb[ib], xsec, probs[il1][ir1][ib]);


    for ( int ir1 = 0; ir1 < Nr; ++ir1 ) for ( int il1 = 0; il1 < Nl; ++il1 ) {
	valarray<double> UT1 = UT[il1][ir1];
	double w1 = vr[ir1]->weight()*vl[il1]->weight();
	sumLR += UT1*w1;
	++nLR;
	sTwLR += UT1*w1*bw;
	sT2wLR += sqr(UT1)*sqr(w1*bw);
	swLR += 1.0;
	//	sum2LR += sqr(UT1)*sqr(w1);
	sum2LR += sqr(UT1)*w1;
	++n2LR;
	sTw2LR += sqr(UT1)*w1*bw;
	sT2w2LR += sqr(sqr(UT1))*sqr(w1*bw);
	sw2LR += 1.0;
	for ( int ir2 = 0; ir2 < Nr; ++ir2 ) for ( int il2 = 0; il2 < Nl; ++il2 ) {
	    valarray<double> UT2 = UT[il2][ir2];
	    double w2 = vr[ir2]->weight()*vl[il2]->weight();
	    //	    valarray<double> UT12 = UT1*w1*UT2*w2;
	    valarray<double> UT12 = UT1*UT2;
	    if ( il1 != il2 && ir1 != ir2 ){
	      sumLR2 += UT12*w1*w2;
	      ++nLR2;
	      sTwLR2 += UT12*w1*w2*bw;
	      sT2wLR2 += sqr(UT12)*sqr(w1*w2*bw);
	      swLR2 += 1.0;
	    }
	    if ( il1 != il2 && ir1 == ir2 ) {
	      sumL2R += UT12*w1*w2;
	      ++nL2R;
	      sTwL2R += UT12*w1*w2*bw;
	      sTwL2R += sqr(UT12)*sqr(w1*w2*bw);
	      swL2R += 1.0;
	    }
	    if ( il1 == il2 && ir1 != ir2 ) {
	      sumR2L += UT12*w1*w2;
	      ++nR2L;
	      sTwR2L += UT12*w1*w2*bw;
	      sTwR2L += sqr(UT12)*sqr(w1*w2*bw);
	      swR2L += 1.0;
	    }
	  }
      }

    sumLR /=double(nLR);
    sum2LR /=double(n2LR);
    sumL2R /=double(nL2R);
    sumR2L /=double(nR2L);
    sumLR2 /=double(nLR2);

    sigtot += 2.0*sumLR*bw;
    sigtot2 += 4.0*sqr(sumLR)*bw2;
    signd += (2.0*sumLR - sum2LR)*bw;
    signd2 += (2.0*sumLR - sum2LR)*(2.0*sumLR - sum2LR)*bw2;
    sigel += sumLR2*bw;
    sigel2 += sumLR2*sumLR2*bw2;
    sigdl += (sumL2R - sumLR2)*bw;
    sigdl2 += (sumL2R - sumLR2)*(sumL2R - sumLR2)*bw2;
    sigdr += (sumR2L - sumLR2)*bw;
    sigdr2 += (sumR2L - sumLR2)*(sumR2L - sumLR2)*bw2;
    sigdd += (sum2LR - sumR2L - sumL2R + sumLR2)*bw;
    sigdd2 += (sum2LR - sumR2L - sumL2R + sumLR2)*(sum2LR - sumR2L - sumL2R + sumLR2)*bw2;
    sigdt += (sum2LR - sumLR2)*bw;
    sigdt2 += (sum2LR - sumLR2)*(sum2LR - sumLR2)*bw2;

    InvEnergy b = vb[ib].bVec().pt();
    double tw = vb[ib].weight()*jac*GeV2;
    for ( int it = 0; it < Nt; ++it ) {
      Energy q = sqrt((double(it) + 0.5)*tMax/double(Nt));
      double J0 = gsl_sf_bessel_J0(b*q);
      double sq = sqrt(q/GeV);
      for ( int i = 0; i < nstr; ++i ) {
	hists[i]->fill(sqr(q)/GeV2, sq*J0*sumLR[i]*bw);
	if ( J0 > 0.0 )
	  hists[i + nstr]->fill(sqr(q)/GeV2, J0*sumLR2[i]*tw);
	else
	  hists[i + 2*nstr]->fill(sqr(q)/GeV2, J0*sumLR2[i]*tw);
      }
    }

    ++ntot;
  }
}

string GlauberAnalysis::getStrat(int i) const {
  string strat = "DIPSY";
  if ( i == 1 ) strat = "black disc";
  if ( i == 2 ) strat = "grey disc";
  if ( i == 3 ) strat = "grey3 disc";
  if ( i == 4 ) strat = "Gaussian";
  if ( i == 5 ) strat = "old black disc";
  if ( i == 6 ) strat = "old grey disc";
  if ( i == 7 ) strat = "old blac disc ND";
  return strat;
}

void GlauberAnalysis::bookHistos() {
  generator()->histogramFactory()->mkdirs("/Glauber");
  generator()->histogramFactory()->mkdirs("/tmp");
  for ( int i = 0; i < nstr; ++i ) {
    string I(1, '0' + i);
    hists[i] = generator()->histogramFactory()->createHistogram1D
      ("/tmp/dSigmadt-a-" + I, Nt, 0.0, tMax/GeV2);
    hists[nstr + i] = generator()->histogramFactory()->createHistogram1D
      ("/tmp/dSigmadt-b-" + I, Nt, 0.0, tMax/GeV2);
    hists[2*nstr + i] = generator()->histogramFactory()->createHistogram1D
      ("/tmp/-dSigmadt-b-" + I, Nt, 0.0, tMax/GeV2);
  }
}

void GlauberAnalysis::print(valarray<double> sig, valarray<double> sig2,
			    int ntot, string xstype) const {
  for ( int i = 0; i < nstr; ++i ) {
    string strat = getStrat(i);
    sig[i] /= double(ntot);
    sig2[i] /= double(ntot);
    double err = sqrt((sig2[i] - sqr(sig[i]))/double(ntot));
    stub(": " + xstype + ": (" + strat + "):")
      << ouniterr(sig[i], err, 1.0) << " nb." << endl;
  }
}

void GlauberAnalysis::finalize(long neve) {
  if ( neve <= 0 ) return;
  print(sigtot, sigtot2, ntot, "Total");
  //  print(2.0*sTwLR, 4.0*sT2wLR, swLR, "Total");
  print(signd, signd2, ntot, "Inelastic(ND)");
  print(sigtot - sigel, sigtot2 - sigel2, ntot, "Inelastic(tot)");
  print(sigel, sigel2, ntot, "Elastic");
  //  print(sTwLR2, sT2wLR2, swLR2, "Elastic");
  print(sigdt, sigdt2, ntot, "Diff. total");
  print(sigdr, sigdr2, ntot, "Diff. exc. (R)");
  print(sigdl, sigdl2, ntot, "Diff. exc. (L)");
  print(sigdd, sigdd2, ntot, "Double diff. exc.");

  for (int i = 0; i < nstr; ++i ) {
    string I(1, '0' + i);
    FactoryBase::tH1DPtr tmp =
      generator()->histogramFactory()->histogramFactory().multiply
      ("/Glauber/dSigmadt-a-" + I, *hists[i], *hists[i]);
    tmp->setTitle("dSigma/dt a (" + getStrat(i) + ")"); 
    tmp->scale(nanobarn*GeV2/(sqr(hbarc)*double(ntot*ntot)*2.0*Constants::pi));
    generator()->histogramFactory()->histogramFactory().destroy(hists[i]);
    tmp = generator()->histogramFactory()->histogramFactory().subtract
      ("/Glauber/dSigmadt-b-" + I, *hists[i + nstr], *hists[i + 2*nstr]);
    tmp->setTitle("dSigma/dt b (" + getStrat(i) + ")"); 
    tmp->scale(1.0/(double(ntot)*2.0*Constants::pi));

  }

}

IBPtr GlauberAnalysis::clone() const {
  return new_ptr(*this);
}

IBPtr GlauberAnalysis::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void GlauberAnalysis::persistentOutput(PersistentOStream & os) const {
  const CrossSection nb = nanobarn;
  os << ounit(nnSigTot, nb) << ounit(nnSigEl, nb) << ounit(nnSigInND, nb)
     << ntot << sigtot << signd << sigel
     << sigdt << sigdl << sigdr << sigdd
     << sigtot2 << signd2 << sigel2
     << sigdt2 << sigdl2 << sigdr2 << sigdd2 << Nt << ounit(tMax,GeV2);
}

void GlauberAnalysis::persistentInput(PersistentIStream & is, int) {
  const CrossSection nb = nanobarn;
  is >> iunit(nnSigTot, nb) >> iunit(nnSigEl, nb) >> iunit(nnSigInND, nb)
     >> ntot >> sigtot >> signd >> sigel
     >> sigdt >> sigdl >> sigdr >> sigdd
     >> sigtot2 >> signd2 >> sigel2
     >> sigdt2 >> sigdl2 >> sigdr2 >> sigdd2 >> Nt >> iunit(tMax, GeV2);
  hists.resize(3*nstr);
  //  bookHistos();
}


// Static variable needed for the type description system in ThePEG.
#include "ThePEG/Utilities/DescribeClass.h"
DescribeClass<GlauberAnalysis,DIPSY::DipoleAnalysisHandler>
  describeDIPSYGlauberAnalysis("DIPSY::GlauberAnalysis",
				   "GlauberAnalysis.so");

void GlauberAnalysis::Init() {

  static ClassDocumentation<GlauberAnalysis> documentation
    ("There is no documentation for the GlauberAnalysis class");

  static Parameter<GlauberAnalysis,CrossSection> interfaceTotalnnXSec
    ("TotalnnXSec",
     "The assumed total nucleon-nucleon cross section assumed in the calculation (in mb).",
     &GlauberAnalysis::nnSigTot, millibarn, 96.0*millibarn, 0.0*millibarn, 0*millibarn,
     true, false, Interface::lowerlim);

  static Parameter<GlauberAnalysis,CrossSection> interfaceElasticnnXSec
    ("ElasticnnXSec",
     "The assumed elastic nucleon-nucleon cross section assumed in the calculation (in mb).",
     &GlauberAnalysis::nnSigEl, millibarn, 24.0*millibarn, 0.0*millibarn, 0*millibarn,
     true, false, Interface::lowerlim);

  static Parameter<GlauberAnalysis,CrossSection> interfaceElasticnnXSecInND
    ("InElasticnnXSec",
     "The assumed inelastic, non-diffractive nucleon-nucleon cross section "
     "assumed in the calculation (in mb).",
     &GlauberAnalysis::nnSigInND, millibarn, 30.0*millibarn, 0.0*millibarn, 0*millibarn,
     true, false, Interface::lowerlim);

  static Parameter<GlauberAnalysis,int> interfaceNt
    ("Nt",
     "Number of bins in t-plots.",
     &GlauberAnalysis::Nt, 1000, 0, 0,
     true, false, Interface::lowerlim);

  static Parameter<GlauberAnalysis,Energy2> interfacetMax
    ("tMax",
     "Maximum t-value in t-plot.",
     &GlauberAnalysis::tMax, GeV2, 10.0*GeV2, 0.0*GeV2, 0*GeV2,
     true, false, Interface::lowerlim);

}

