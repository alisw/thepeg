// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the NoFluctXSecAnalysis class.
//

#include "NoFluctXSecAnalysis.h"
#include "DipoleXSec.h"
#include "DipoleState.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"

#include "gsl/gsl_sf_bessel.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace DIPSY;

NoFluctXSecAnalysis::NoFluctXSecAnalysis() {}

NoFluctXSecAnalysis::~NoFluctXSecAnalysis() {}

void NoFluctXSecAnalysis::initialize() {
  generator()->histogramFactory()->initrun();
  generator()->histogramFactory()->registerClient(this);

  n00 = nL0 = n2L0 = n0R = nLR = 0;
  sum00 = sumL0 = sum2L0 = sum0R = sumLR = ZERO;
  sum002 = sumL02 = sum2L02 = sum0R2 = sumLR2 = ZERO;

}

void NoFluctXSecAnalysis::
analyze(const vector<DipoleStatePtr> & vr, const vector<DipoleStatePtr> & vl,
	const vector<ImpactParameters> & vb, const DipoleXSec & xsec,
	const Vec3D & probs, double jac) {

  int Nr = vr.size();
  int Nl = vl.size();
  if ( Nr*Nl*vb.size() == 0 ) return;

  for ( int ib = 0, Nb = vb.size(); ib < Nb; ++ib ) {
    CrossSection bweight = sqr(hbarc)*vb[ib].weight()*jac;
    vector<double> UTL(vr.size(), 0.0);
    vector<double> avTL(vr.size(), 0.0);
    vector<double> avTL2(vr.size(), 0.0);
    vector<double> UTR(vl.size(), 0.0);
    vector<double> avTR(vl.size(), 0.0);
    double UTLR = 0.0;
    double avTLR = 0.0;
    for ( int ir = 0; ir < Nr; ++ir ) for ( int il = 0; il < Nl; ++il ) {
	double T = probs[il][ir][ib];
	double weight = vr[ir]->weight()*vl[il]->weight();
	double UT = 2.0*xsec.unitarize(T);
	UTL[ir] += UT*weight;
	avTL[ir] += T*weight;
	if ( il > 0 ) avTL2[ir] += T*weight;
	UTR[il] += UT*weight;
	avTR[il] += T*weight;
	UTLR += UT*weight;
	avTLR += T*weight;
      }
 
    UTLR /= double(Nr*Nl);
    sum00 += UTLR*bweight;
    sum002 += sqr(UTLR*bweight);
    n00 += 1;
    if ( Nr > 1 && Nl > 1) {
      avTLR = (2.0*xsec.unitarize(avTLR/double(Nr*Nl))*double(Nr*Nl) - UTLR)/double(Nr*Nl - 1);
      sumLR += avTLR*bweight;
      sumLR2 += sqr(avTLR*bweight);
      nLR += 1;
    }

    if ( Nl > 2 ) for ( int ir = 0; ir < Nr; ++ir ) {
	UTL[ir] /= double(Nl);
	avTL[ir] = (2.0*xsec.unitarize(avTL[ir]/double(Nl))*double(Nl) - UTL[ir])/double(Nl - 1);
	sumL0 += avTL[ir]*bweight;
	sumL02 += sqr(avTL[ir]*bweight);
	nL0 += 1;
	if ( Nl > 2 ) {
	  avTL2[ir] = (2.0*xsec.unitarize(avTL2[ir]/double(Nl - 1))*double(Nl - 1)
		       - UTL[ir])/double(Nl - 2);
	  avTL2[ir] = double(Nl)*avTL[ir] - double(Nl - 1)*avTL2[ir];
	  sum2L0 += avTL2[ir]*bweight;
	  sum2L02 += sqr(avTL2[ir]*bweight);
	  n2L0 += 1;
	}
      }


    if ( Nr > 1 ) for ( int il = 0; il < Nl; ++il ) {
	UTR[il] /= double(Nr);
	avTR[il] = (2.0*xsec.unitarize(avTR[il]/double(Nr))*double(Nr) - UTR[il])/double(Nr - 1);
	sum0R += avTR[il]*bweight;
	sum0R2 += sqr(avTR[il]*bweight);
	n0R += 1;
      }
  }
}

void NoFluctXSecAnalysis::finalize(long neve) {
  if ( neve <= 0 || n00 <= 0 ) return;

  sum00 /= double(n00);
  CrossSection err00 = sqrt((sum002/double(n00) - sqr(sum00))/double(n00));
  generator()->log().setf(ios::left, ios::adjustfield);
  generator()->log()
    << setw(50) << name() + ": Crosscheck totxsec:"
    << ouniterr(sum00, err00, nanobarn) << " nb." << endl;
  if ( nLR > 0 ) {
    sumLR /= double(nLR);
    CrossSection errLR = sqrt((sumLR2/double(nLR) - sqr(sumLR))/double(nLR));
    generator()->log()
      << setw(50) << name() + ": xsec average over L+R:"
      << ouniterr(sumLR, errLR, nanobarn) << " nb." << endl;
  }
  if ( nL0 > 0 ) {
    sumL0 /= double(nL0);
    CrossSection errL0 = sqrt((sumL02/double(nL0) - sqr(sumL0))/double(nL0));
    generator()->log()
      << setw(50) << name() + ": xsec average over L:"
      << ouniterr(sumL0, errL0, nanobarn) << " nb." << endl;
  }
  if ( n2L0 > 0 ) {
    sum2L0 /= double(n2L0);
    CrossSection err2L0 = sqrt((sum2L02/double(n2L0) - sqr(sum2L0))/double(n2L0));
    generator()->log()
      << setw(50) << name() + ": xsec average over L(2):"
      << ouniterr(sum2L0, err2L0, nanobarn) << " nb." << endl;
  }
  if ( n0R > 0 ) {
    sum0R /= double(n0R);
    CrossSection err0R = sqrt((sum0R2/double(n0R) - sqr(sum0R))/double(n0R));
    generator()->log()
      << setw(50) << name() + ": xsec average over R:"
      << ouniterr(sum0R, err0R, nanobarn) << " nb." << endl;
  }
    
}

IBPtr NoFluctXSecAnalysis::clone() const {
  return new_ptr(*this);
}

IBPtr NoFluctXSecAnalysis::fullclone() const {
  return new_ptr(*this);
}


// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void NoFluctXSecAnalysis::persistentOutput(PersistentOStream & os) const {
  const CrossSection2 nb2= sqr(nanobarn);
  os << ounit(sum00, nanobarn) << ounit(sumL0, nanobarn) << ounit(sum2L0, nanobarn) << ounit(sum0R, nanobarn)
     << ounit(sumLR, nanobarn) << ounit(sum002, nb2) << ounit(sumL02, nb2) << ounit(sum2L02, nb2)
     << ounit(sum0R2, nb2) << ounit(sumLR2, nb2) << n00 << nL0 << n2L0 << n0R << nLR;
}

void NoFluctXSecAnalysis::persistentInput(PersistentIStream & is, int) {
  const CrossSection2 nb2= sqr(nanobarn);
  is >> iunit(sum00, nanobarn) >> iunit(sumL0, nanobarn) >> iunit(sum2L0, nanobarn) >> iunit(sum0R, nanobarn)
     >> iunit(sumLR, nanobarn) >> iunit(sum002, nb2) >> iunit(sumL02, nb2) >> iunit(sum2L02, nb2)
     >> iunit(sum0R2, nb2) >> iunit(sumLR2, nb2) >> n00 >> nL0 >> n2L0 >> n0R >> nLR;
}


// Static variable needed for the type description system in ThePEG.
#include "ThePEG/Utilities/DescribeClass.h"
DescribeClass<NoFluctXSecAnalysis,DIPSY::DipoleAnalysisHandler>
  describeDIPSYNoFluctXSecAnalysis("DIPSY::NoFluctXSecAnalysis",
				   "NoFluctXSecAnalysis.so");

void NoFluctXSecAnalysis::Init() {

  static ClassDocumentation<NoFluctXSecAnalysis> documentation
    ("There is no documentation for the NoFluctXSecAnalysis class");

}

