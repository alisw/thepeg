// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SemiInclusiveXSecAnalysis class.
//

#include "SemiInclusiveXSecAnalysis.h"
#include "DipoleXSec.h"
#include "DipoleState.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"

#include "gsl/gsl_sf_bessel.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace DIPSY;

SemiInclusiveXSecAnalysis::SemiInclusiveXSecAnalysis() {}

SemiInclusiveXSecAnalysis::~SemiInclusiveXSecAnalysis() {}

void SemiInclusiveXSecAnalysis::initialize() {
  ntot = 0;
  sigtot = signd = sigel = sigdt = sigdr = sigdl = sigdd = ZERO;
  sigtot2 = signd2 = sigel2 = sigdt2 = sigdr2 = sigdl2 = sigdd2 = ZERO;
}

void SemiInclusiveXSecAnalysis::
analyze(const vector<DipoleStatePtr> & vr, const vector<DipoleStatePtr> & vl,
	const vector<ImpactParameters> & vb, const DipoleXSec & xsec,
	const Vec3D & probs, double jac) {

  int Nr = vr.size();
  int Nl = vl.size();
  if ( Nr*Nl*vb.size() == 0 ) return;

  for ( int ib = 0, Nb = vb.size(); ib < Nb; ++ib ) {
    CrossSection bweight = sqr(hbarc)*vb[ib].weight()*jac;

    double sumLR = 0.0; // T averaged over L and R
    int nLR = 0;
    double sum2LR = 0.0; // T^2 average over L and R
    int n2LR = 0;
    double sumLR2 = 0.0; // (T averaged over L and R )^2
    int nLR2 = 0;
    double sumL2R = 0.0; // (T averaged over L)^2 averaged over R
    int nL2R = 0;
    double sumR2L = 0.0; // (T averaged over R)^2 averaged over L
    int nR2L = 0;

    for ( int ir1 = 0; ir1 < Nr; ++ir1 ) for ( int il1 = 0; il1 < Nl; ++il1 ) {
	double T1 = probs[il1][ir1][ib];
	double w1 = vr[ir1]->weight()*vl[il1]->weight();
	double UT1 = xsec.unitarize(T1);
	sumLR += UT1*w1;
	++nLR;
	sum2LR += sqr(UT1*w1);
	++n2LR;
	for ( int ir2 = 0; ir2 < Nr; ++ir2 ) for ( int il2 = 0; il2 < Nl; ++il2 ) {
	    double T2 = probs[il2][ir2][ib];
	    double w2 = vr[ir2]->weight()*vl[il2]->weight();
	    double UT2 = xsec.unitarize(T2);
	    double UT = UT1*w1*UT2*w2;
	    if ( il1 != il2 && ir1 != ir2 ){
	      sumLR2 += UT;
	      ++nLR2;
	    }
	    if ( il1 != il2 && ir1 == ir2 ) {
	      sumL2R += UT;
	      ++nL2R;
	    }
	    if ( il1 == il2 && ir1 != ir2 ) {
	      sumR2L += UT;
	      ++nR2L;
	    }
	  }
      }

    sumLR /=double(nLR);
    sum2LR /=double(n2LR);
    sumL2R /=double(nL2R);
    sumR2L /=double(nR2L);
    sumLR2 /=double(nLR2);

    CrossSection bw = bweight;
    CrossSection2 bw2 = sqr(bw);

    sigtot += 2.0*sumLR*bw;
    sigtot2 += sqr(2.0*sumLR)*bw2;
    signd += (2.0*sumLR - sum2LR)*bw;
    signd2 += sqr(2.0*sumLR - sum2LR)*bw2;
    sigel += sumLR2*bw;
    sigel2 += sqr(sumLR2)*bw2;
    sigdl += (sumL2R - sumLR2)*bw;
    sigdl2 += sqr(sumL2R - sumLR2)*bw2;
    sigdr += (sumR2L - sumLR2)*bw;
    sigdr2 += sqr(sumR2L - sumLR2)*bw2;
    sigdd += (sum2LR - sumR2L - sumL2R + sumLR2)*bw;
    sigdd2 += sqr(sum2LR - sumR2L - sumL2R + sumLR2)*bw2;
    sigdt += (sum2LR - sumLR2)*bw;
    sigdt2 += sqr(sum2LR - sumLR2)*bw2;

    ++ntot;
  }
}

void SemiInclusiveXSecAnalysis::finalize(long neve) {
  if ( neve <= 0 ) return;

  sigtot /= double(ntot);
  sigtot2 /= double(ntot);
  CrossSection err = sqrt((sigtot2 - sqr(sigtot))/double(ntot));
  stub(": Total:") << ouniterr(sigtot, err, nanobarn) << " nb." << endl;

  signd /= double(ntot);
  signd2 /= double(ntot);
  err = sqrt((signd2 - sqr(signd))/double(ntot));
  stub(": Inelastic (ND):") << ouniterr(signd, err, nanobarn) << " nb." << endl;

  if ( sigel2 > ZERO ) {
    sigel /= double(ntot);
    sigel2 /= double(ntot);
    err = sqrt((sigel2 - sqr(sigel))/double(ntot));
    stub(": Elastic:") << ouniterr(sigel, err, nanobarn) << " nb." << endl;
  }
  
  if ( sigdt2 > ZERO ) {
    sigdt /= double(ntot);
    sigdt2 /= double(ntot);
    err = sqrt((sigdt2 - sqr(sigdt))/double(ntot));
    stub(": Diffractive total:") << ouniterr(sigdt, err, nanobarn) << " nb." << endl;
  }

  if ( sigdr2 > ZERO ) {
    sigdr /= double(ntot);
    sigdr2 /= double(ntot);
    err = sqrt((sigdr2 - sqr(sigdr))/double(ntot));
    stub(": Diffractive excitation (R):")
      << ouniterr(sigdr, err, nanobarn) << " nb." << endl;
  }

  if ( sigdl2 > ZERO ) {
    sigdl /= double(ntot);
    sigdl2 /= double(ntot);
    err = sqrt((sigdl2 - sqr(sigdl))/double(ntot));
    stub(": Diffractive excitation (L):")
      << ouniterr(sigdl, err, nanobarn) << " nb." << endl;
  }

  if ( sigdd2 > ZERO ) {
    sigdd /= double(ntot);
    sigdd2 /= double(ntot);
    err = sqrt((sigdd2 - sqr(sigdd))/double(ntot));
    stub(": Double diffractive excitation:")
      << ouniterr(sigdd, err, nanobarn) << " nb." << endl;
  }
}

IBPtr SemiInclusiveXSecAnalysis::clone() const {
  return new_ptr(*this);
}

IBPtr SemiInclusiveXSecAnalysis::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void SemiInclusiveXSecAnalysis::persistentOutput(PersistentOStream & os) const {
  const CrossSection nb = nanobarn;
  const CrossSection2 nb2= sqr(nb);
  os << ntot << ounit(sigtot, nb) << ounit(signd, nb) << ounit(sigel, nb)
     << ounit(sigdt, nb) << ounit(sigdl, nb) << ounit(sigdr, nb) << ounit(sigdd, nb)
     << ounit(sigtot2, nb2) << ounit(signd2, nb2) << ounit(sigel2, nb2)
     << ounit(sigdt2, nb2) << ounit(sigdl2, nb2) << ounit(sigdr2, nb2) << ounit(sigdd2, nb2);
}

void SemiInclusiveXSecAnalysis::persistentInput(PersistentIStream & is, int) {
  const CrossSection nb = nanobarn;
  const CrossSection2 nb2= sqr(nb);
  is >> ntot >> iunit(sigtot, nb) >> iunit(signd, nb) >> iunit(sigel, nb)
     >> iunit(sigdt, nb) >> iunit(sigdl, nb) >> iunit(sigdr, nb) >> iunit(sigdd, nb)
     >> iunit(sigtot2, nb2) >> iunit(signd2, nb2) >> iunit(sigel2, nb2)
     >> iunit(sigdt2, nb2) >> iunit(sigdl2, nb2) >> iunit(sigdr2, nb2) >> iunit(sigdd2, nb2);
}


// Static variable needed for the type description system in ThePEG.
#include "ThePEG/Utilities/DescribeClass.h"
DescribeClass<SemiInclusiveXSecAnalysis,DIPSY::DipoleAnalysisHandler>
  describeDIPSYSemiInclusiveXSecAnalysis("DIPSY::SemiInclusiveXSecAnalysis",
				   "SemiInclusiveXSecAnalysis.so");

void SemiInclusiveXSecAnalysis::Init() {

  static ClassDocumentation<SemiInclusiveXSecAnalysis> documentation
    ("There is no documentation for the SemiInclusiveXSecAnalysis class");

}

