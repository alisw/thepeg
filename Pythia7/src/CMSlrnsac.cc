// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the CMSlrnsac class.
//

#include "CMSlrnsac.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Utilities/DescribeClass.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#ifndef LWH_AIAnalysisFactory_H
#ifndef LWH 
#define LWH ThePEGLWH
#endif
#include "ThePEG/Analysis/LWH/AnalysisFactory.h"
#endif

using namespace ThePEG;

CMSlrnsac::CMSlrnsac(): doRweight(2) {}

CMSlrnsac::CMSlrnsac(const CMSlrnsac & x): AnalysisHandler(x),
  doRweight(x.doRweight), externalFile(x.externalFile) {}

CMSlrnsac::~CMSlrnsac() {}

void CMSlrnsac::analyze(const tPVector & inparticles, double weight) {
  
  PVector fileparticles;
  tPVector particles;

  if ( externalFile.empty() ) {
    particles = inparticles;
  } else {
    while ( file.readline() && !file.find("*** BEGIN EVENT ***") ) {
      if ( !file.find("*** END EVENT ***") )
	generator()->misc() << file.getline();
    }
    int n = 0;
    if ( !file.readline() ) return;
    file >> n >> weight;
    if ( !file ) return;
    for ( int i = 0; i < n; ++i ) {
      int id = 0;
      double px = 0.0;
      double py = 0.0;
      double pz = 0.0;
      double pe = 0.0;
      double pm = 0.0;
      if ( !file.readline() ) return;
      file >> id >> px >> py >> pz >> pe >> pm;
      PPtr p = generator()->getParticle(id);
      p->setMomentum(Lorentz5Momentum(px*GeV, py*GeV, pz*GeV, pe*GeV, pm*GeV));
      fileparticles.push_back(p);
      particles.push_back(p);
    }
  }

  vector<int> nch(5);
  vector< vector<double> > thisphi(5);
  vector< vector<double> > thiseta(5);
  vector<double> allphi;
  vector<double> alleta;
  vector<double> phi13;
  vector<double> eta13;
  int Nch04 = 0;
  for ( int i = 0, N = particles.size(); i < N; ++i ) {
    if ( particles[i]->data().iCharge() == 0 ) continue;
    double eta =particles[i]->eta();
    if ( abs(eta) >= etamax ) continue;
    Energy pt = particles[i]->momentum().perp();
    if ( pt < 0.1*GeV ) continue;
    if ( pt >= 0.4*GeV ) ++Nch04;
    int ipt = ptbin(pt);
    double phi = particles[i]->momentum().phi();
    thisphi[ipt].push_back(phi);
    thiseta[ipt].push_back(eta);
    allphi.push_back(phi);
    alleta.push_back(eta);
    if ( ipt == 1 || ipt == 2 ) {
      phi13.push_back(phi);
      eta13.push_back(eta);
    }
  }

  if ( Nch04 < 2 ) return;

  int inch = nchbin(Nch04);

  if ( externalFile.empty() ) {
    LorentzPoint b1 =
      particles[0]->birthStep()->collision()->incoming().first->vertex();
    LorentzPoint b2 =
      particles[0]->birthStep()->collision()->incoming().second->vertex();
    double b = sqrt(sqr(b1.x() - b2.x()) + sqr(b1.y() - b2.y()))/femtometer;
    bdist[inch]->fill(b, weight);
  }

  if ( refNch04[inch] == 0 ) {
    refNch04[inch] = Nch04;
    refphi[inch] = thisphi;
    refeta[inch] = thiseta;
    refweight[inch] = weight;
    return;
  }
  
  mult->fill(Nch04, weight);

  double rweight = weight;
  switch ( doRweight ) {
  case 1:
    rweight = sqrt(weight*refweight[inch]);
    break;
  case 2:
    rweight = weight*refweight[inch];
    break;
  }
  sumweight[inch] += weight;
  sumrefweight[inch] += rweight;

  // First do phi-projections.
  for ( int ipt = 0; ipt < 4; ++ipt ) {
    multbin[inch][ipt]->fill(thisphi[ipt].size(), weight);
    multbin2[inch][ipt]->fill(thisphi[ipt].size(), weight);
    double wt =
      2.0*weight/max(double(thisphi[ipt].size()*(thisphi[ipt].size() - 1)), 1.0);
    double rwt =
      rweight/max(double(thisphi[ipt].size()*refphi[inch][ipt].size()), 1.0);

    for ( int i = 0, N = thisphi[ipt].size(); i < N; ++i ) {

      // Start with signal sample.
      for ( int j = i + 1; j < N; ++ j ) {
	double deta = abs(thiseta[ipt][i] - thiseta[ipt][j]);
	if ( deta < 2.0 ) continue;
	double dphi = abs(thisphi[ipt][i] - thisphi[ipt][j]);
	if ( dphi > Constants::pi ) dphi = 2.0*Constants::pi - dphi;
	Sphi[inch][ipt]->fill(dphi, wt);
	S2phi[inch][ipt]->fill(deta, dphi, wt);
	S3phi[inch][ipt]->fill(dphi, weight/double(thisphi[ipt].size()));
	//	S3phi[inch][ipt]->fill(dphi, wt);
      }
      // Then background sample.
      for ( int j = 0, M = refphi[inch][ipt].size(); j < M; ++j ) {
	double deta = abs(thiseta[ipt][i] - refeta[inch][ipt][j]);
	if ( deta < 2.0 ) continue;
	double dphi = abs(thisphi[ipt][i] - refphi[inch][ipt][j]);
	if ( dphi > Constants::pi ) dphi = 2.0*Constants::pi - dphi;
	Bphi[inch][ipt]->fill(dphi, rwt);
	B2phi[inch][ipt]->fill(deta, dphi, rwt);
	B3phi[inch][ipt]->fill(dphi, rwt);
      }
	
    }
  }

  // Now do 2d histos. Begin with 1<pt<3 bins
  vector<double> refallphi = refphi[inch][1];
  refallphi.insert(refallphi.end(),
		   refphi[inch][2].begin(), refphi[inch][2].end());
  vector<double> refalleta = refeta[inch][1];
  refalleta.insert(refalleta.end(),
		   refeta[inch][2].begin(), refeta[inch][2].end());

  double wt = 2.0*weight/max(double(phi13.size()*(phi13.size() - 1)), 1.0);
  double rwt = rweight/max(double(phi13.size()*refallphi.size()), 1.0);

  double wt2 = 2.0*weight/max(double(phi13.size()), 1.0);

  multMB13->fill(phi13.size(), weight);
  if ( inch == 3 ) multHN13->fill(phi13.size(), weight);

  for ( int i = 0, N = phi13.size(); i < N; ++i ) {

    // First Signal.
    for ( int j = i + 1; j < N; ++j ) {
      double deta = abs(eta13[i] - eta13[j]);
      double dphi = abs(phi13[i] - phi13[j]);
      if ( dphi > Constants::pi ) dphi = 2.0*Constants::pi - dphi;
      SMB13->fill(deta, dphi, wt);
      SMB13->fill(-deta, dphi, wt);
      S2MB13->fill(deta, dphi, wt2);
      S2MB13->fill(-deta, dphi, wt2);
      if ( inch == 3 ) {
	SHN13->fill(deta, dphi, wt);
	SHN13->fill(-deta, dphi, wt);
	S2HN13->fill(deta, dphi, wt2);
	S2HN13->fill(-deta, dphi, wt2);
      }
      dphi = dphi <= 0.5*Constants::pi? -dphi: 2.0*Constants::pi - dphi;
      SMB13->fill(deta, dphi, wt);
      SMB13->fill(-deta, dphi, wt);
      S2MB13->fill(deta, dphi, wt2);
      S2MB13->fill(-deta, dphi, wt2);
      if ( inch == 3 ) {
	SHN13->fill(deta, dphi, wt);
	SHN13->fill(-deta, dphi, wt);
	S2HN13->fill(deta, dphi, wt2);
	S2HN13->fill(-deta, dphi, wt2);
      }
    }

    // Then background.
    for ( int j = 0, M = refallphi.size(); j < M; ++j ) {
      double deta = abs(eta13[i] - refalleta[j]);
      double dphi = abs(phi13[i] - refallphi[j]);
      if ( dphi > Constants::pi ) dphi = 2.0*Constants::pi - dphi;
      BMB13->fill(deta, dphi, rwt);
      BMB13->fill(-deta, dphi, rwt);
      B2MB13->fill(deta, dphi, rwt);
      B2MB13->fill(-deta, dphi, rwt);
      if ( inch == 3 ) {
	BHN13->fill(deta, dphi, rwt);
	BHN13->fill(-deta, dphi, rwt);
	B2HN13->fill(deta, dphi, rwt);
	B2HN13->fill(-deta, dphi, rwt);
      }
      dphi = dphi <= 0.5*Constants::pi? -dphi: 2.0*Constants::pi - dphi;
      BMB13->fill(deta, dphi, rwt);
      BMB13->fill(-deta, dphi, rwt);
      B2MB13->fill(deta, dphi, rwt);
      B2MB13->fill(-deta, dphi, rwt);
      if ( inch == 3 ) {
	BHN13->fill(deta, dphi, rwt);
	BHN13->fill(-deta, dphi, rwt);
	B2HN13->fill(deta, dphi, rwt);
	B2HN13->fill(-deta, dphi, rwt);
      }
    }
  }

  // Continue 2D with 0.1<pt bins
  refallphi.insert(refallphi.end(),
		   refphi[inch][0].begin(), refphi[inch][0].end());
  refallphi.insert(refallphi.end(),
		   refphi[inch][3].begin(), refphi[inch][3].end());
  refallphi.insert(refallphi.end(),
		   refphi[inch][4].begin(), refphi[inch][4].end());
  refalleta.insert(refalleta.end(),
		   refeta[inch][0].begin(), refeta[inch][0].end());
  refalleta.insert(refalleta.end(),
		   refeta[inch][3].begin(), refeta[inch][3].end());
  refalleta.insert(refalleta.end(),
		   refeta[inch][4].begin(), refeta[inch][4].end());
  
  wt = 2.0*weight/max(double(allphi.size()*(allphi.size() - 1)), 1.0);
  rwt = rweight/max(double(allphi.size()*refallphi.size()), 1.0);
  wt2 = 2.0*weight/max(double(allphi.size()), 1.0);

  multMB01->fill(allphi.size(), weight);
  if ( inch == 3 ) multHN01->fill(allphi.size(), weight);

  for ( int i = 0, N = allphi.size(); i < N; ++i ) {

    // First Signal.
    for ( int j = i + 1; j < N; ++j ) {
      double deta = abs(alleta[i] - alleta[j]);
      double dphi = abs(allphi[i] - allphi[j]);
      if ( dphi > Constants::pi ) dphi = 2.0*Constants::pi - dphi;
      SMB01->fill(deta, dphi, wt);
      SMB01->fill(-deta, dphi, wt);
      S2MB01->fill(deta, dphi, wt2);
      S2MB01->fill(-deta, dphi, wt2);
      if ( inch == 3 ) {
	SHN01->fill(deta, dphi, wt);
	SHN01->fill(-deta, dphi, wt);
	S2HN01->fill(deta, dphi, wt2);
	S2HN01->fill(-deta, dphi, wt2);
      }
      dphi = ( dphi <= 0.5*Constants::pi? -dphi: 2.0*Constants::pi - dphi );
      SMB01->fill(deta, dphi, wt);
      SMB01->fill(-deta, dphi, wt);
      S2MB01->fill(deta, dphi, wt2);
      S2MB01->fill(-deta, dphi, wt2);
     if ( inch == 3 ) {
	SHN01->fill(deta, dphi, wt);
	SHN01->fill(-deta, dphi, wt);
	S2HN01->fill(deta, dphi, wt2);
	S2HN01->fill(-deta, dphi, wt2);
      }
    }

    // Then background.
    for ( int j = 0, M = refallphi.size(); j < M; ++j ) {
      double deta = abs(alleta[i] - refalleta[j]);
      double dphi = abs(allphi[i] - refallphi[j]);
      if ( dphi > Constants::pi ) dphi = 2.0*Constants::pi - dphi;
      BMB01->fill(deta, dphi, rwt);
      BMB01->fill(-deta, dphi, rwt);
      B2MB01->fill(deta, dphi, rwt);
      B2MB01->fill(-deta, dphi, rwt);
      if ( inch == 3 ) {
	BHN01->fill(deta, dphi, rwt);
	BHN01->fill(-deta, dphi, rwt);
	B2HN01->fill(deta, dphi, rwt);
	B2HN01->fill(-deta, dphi, rwt);
      }
      dphi = dphi <= 0.5*Constants::pi? -dphi: 2.0*Constants::pi - dphi;
      BMB01->fill(deta, dphi, rwt);
      BMB01->fill(-deta, dphi, rwt);
      B2MB01->fill(deta, dphi, rwt);
      B2MB01->fill(-deta, dphi, rwt);
      if ( inch == 3 ) {
	BHN01->fill(deta, dphi, rwt);
	BHN01->fill(-deta, dphi, rwt);
	B2HN01->fill(deta, dphi, rwt);
	B2HN01->fill(-deta, dphi, rwt);
      }
    }
  }

  refNch04[inch] = Nch04;
  refphi[inch] = thisphi;
  refeta[inch] = thiseta;
  refweight[inch] = weight;

}

IBPtr CMSlrnsac::clone() const {
  return new_ptr(*this);
}

IBPtr CMSlrnsac::fullclone() const {
  return new_ptr(*this);
}

void CMSlrnsac::dofinish() {
  AnalysisHandler::dofinish();

  if ( !externalFile.empty() ) file.close();

  unitNormalize(mult);
  double sumw = 0.0;
  double sumrefw = 0.0;
  for ( int inch = 0; inch < 4; ++inch ) {
    unitNormalize(bdist[inch]);
    sumw += sumweight[inch];
    sumrefw += sumrefweight[inch];
    for ( int ipt = 0; ipt < 4; ++ipt ) {
      if ( sumweight[inch] > 0.0 && sumrefweight[inch] > 0.0 ) {
	Sphi[inch][ipt]->scale(1.0/sumweight[inch]);
	Bphi[inch][ipt]->scale(1.0/sumrefweight[inch]);
	S2phi[inch][ipt]->scale(1.0/sumweight[inch]);
	B2phi[inch][ipt]->scale(1.0/sumrefweight[inch]);
      }
      double weve = nozero(multbin2[inch][ipt]->sumBinHeights());
      S3phi[inch][ipt]->scale(1.0/weve);
      double avN = multbin2[inch][ipt]->mean();
      double intg = nozero(S3phi[inch][ipt]->sumBinHeights());
      S3phi[inch][ipt]->scale((avN - 1.0)/intg);
      B3phi[inch][ipt]->scale(S3phi[inch][ipt]->sumBinHeights()/
			      nozero(B3phi[inch][ipt]->sumBinHeights()));

      string indx;
      indx += '0' + inch + 1;
      indx += '0' + ipt + 1;
      tH1DPtr tmp = iHistogramFactory().subtract
	(dir + "/SmBphi" + indx, *Sphi[inch][ipt], *Bphi[inch][ipt]);
      iHistogramFactory().divide
	(dir + "/Rphi" + indx, *tmp, *Bphi[inch][ipt]);
      tH2DPtr tmp2 = iHistogramFactory().subtract
	(dir + "/S2mBphi" + indx, *S2phi[inch][ipt], *B2phi[inch][ipt]);
      tmp2 = iHistogramFactory().divide
	(dir + "/R2phi" + indx, *tmp2, *B2phi[inch][ipt]);
      tmp = iHistogramFactory().projectionY
	(dir + "/Rphiproj" + indx, *tmp2);
      tmp->scale(0.4); // The bin width.
      tmp = iHistogramFactory().subtract
	(dir + "/S3mBphi" + indx, *S3phi[inch][ipt], *B3phi[inch][ipt]);
      tmp = iHistogramFactory().divide
	(dir + "/R3phi" + indx, *tmp, *B3phi[inch][ipt]);
      tmp->scale(S3phi[inch][ipt]->sumBinHeights());
      unitNormalize(multbin[inch][ipt]);
    }
  }
  if ( sumw > 0.0 && sumrefw > 0.0 ) {
    SMB01->scale(1.0/sumw);
    BMB01->scale(1.0/sumrefw);
    SMB13->scale(1.0/sumw);
    BMB13->scale(1.0/sumrefw);
  }
  if ( sumweight[3] > 0.0 && sumrefweight[3] > 0.0 ) {
    SHN01->scale(1.0/sumweight[3]);
    BHN01->scale(1.0/sumrefweight[3]);
    SHN13->scale(1.0/sumweight[3]);
    BHN13->scale(1.0/sumrefweight[3]);
  }

  double weve = nozero(multMB01->sumBinHeights());
  S2MB01->scale(1.0/weve);
  B2MB01->scale(S2MB01->sumBinHeights()/nozero(B2MB01->sumBinHeights()));
  weve = nozero(multMB13->sumBinHeights());
  S2MB13->scale(1.0/weve);
  B2MB13->scale(S2MB13->sumBinHeights()/nozero(B2MB13->sumBinHeights()));
  weve = nozero(multHN01->sumBinHeights());
  S2HN01->scale(1.0/weve);
  B2HN01->scale(S2HN01->sumBinHeights()/nozero(B2HN01->sumBinHeights()));
  weve = nozero(multHN13->sumBinHeights());
  S2HN13->scale(1.0/weve);
  B2HN13->scale(S2HN13->sumBinHeights()/nozero(B2HN13->sumBinHeights()));

  tH2DPtr
  tmp = iHistogramFactory().subtract(dir + "/SmBMB01", *SMB01, *BMB01);
  iHistogramFactory().divide(dir + "/RMB01", *tmp, *BMB01);
  tmp = iHistogramFactory().subtract(dir + "/SmBHN01", *SHN01, *BHN01);
  iHistogramFactory().divide(dir + "/RHN01", *tmp, *BHN01);
  tmp = iHistogramFactory().subtract(dir + "/SmBMB13", *SMB13, *BMB13);
  iHistogramFactory().divide(dir + "/RMB13", *tmp, *BMB13);
  tmp = iHistogramFactory().subtract(dir + "/SmBHN13", *SHN13, *BHN13);
  iHistogramFactory().divide(dir + "/RHN13", *tmp, *BHN13);

  tmp = iHistogramFactory().subtract(dir + "/S2mBMB01", *S2MB01, *B2MB01);
  tmp = iHistogramFactory().divide(dir + "/R2MB01", *tmp, *B2MB01);
  tmp->scale(S2MB01->sumBinHeights());
  tmp = iHistogramFactory().subtract(dir + "/S2mBHN01", *S2HN01, *B2HN01);
  tmp =iHistogramFactory().divide(dir + "/R2HN01", *tmp, *B2HN01);
  tmp->scale(S2MB13->sumBinHeights());
  tmp = iHistogramFactory().subtract(dir + "/S2mBMB13", *S2MB13, *B2MB13);
  tmp = iHistogramFactory().divide(dir + "/R2MB13", *tmp, *B2MB13);
  tmp->scale(S2HN01->sumBinHeights());
  tmp = iHistogramFactory().subtract(dir + "/S2mBHN13", *S2HN13, *B2HN13);
  tmp = iHistogramFactory().divide(dir + "/R2HN13", *tmp, *B2HN13);
  tmp->scale(S2HN13->sumBinHeights());

  tmp = iHistogramFactory().createHistogram2D
    (dir + "/S0",
     33, -2.0*etamax, 2.0*etamax,
     30, -0.5*Constants::pi,1.5*Constants::pi);
  for ( int ieta = 0; ieta < 33; ++ieta ) {
    double eta = (ieta + 0.5)*4.0*etamax/33.0 -2.0*etamax;
    double w = (2.0*etamax - abs(eta))/(32.0*Constants::pi*sqr(etamax));
    for ( int iphi = 0; iphi < 30; ++iphi ) {
      double phi = (iphi + 0.5)*2.0*Constants::pi/30.0 - 0.5*Constants::pi;
      tmp->fill(eta, phi, w);
    }
  }
  iHistogramFactory().divide(dir + "/S0MB01", *BMB01, *tmp);
  iHistogramFactory().divide(dir + "/S0HN01", *BHN01, *tmp);
  iHistogramFactory().divide(dir + "/S0MB13", *BMB13, *tmp);
  iHistogramFactory().divide(dir + "/S0HN13", *BHN13, *tmp);


}

void CMSlrnsac::doinitrun() {
  AnalysisHandler::doinitrun();
  etamax = 2.4;
  histogramFactory().registerClient(this);
  dir = "/CMSlrnsac";
  if ( doRweight == 1 ) dir = "/CMSlrnsac1"; 
  if ( doRweight == 2 ) dir = "/CMSlrnsac2"; 
  histogramFactory().mkdirs(dir);
  Bphi = Sphi = B3phi = S3phi = multbin = multbin2  =
    vector< vector<tH1DPtr> >(4, vector<tH1DPtr>(4));
  B2phi = S2phi = vector< vector<tH2DPtr> >(4, vector<tH2DPtr>(4));
  bdist = vector<tH1DPtr>(4);
  for ( int inch = 0; inch < 4; ++inch ) {
    string indx;
    indx += '0' + inch + 1;
    bdist[inch] = histogramFactory().createHistogram1D(dir + "/bdist" + indx,
						       50, 0.0, 5.0);
    for ( int ipt = 0; ipt < 4; ++ ipt ) {
      string indxx = indx;
      indxx += '0' + ipt + 1;
      Sphi[inch][ipt] =
	histogramFactory().createHistogram1D(dir + "/Sphi" + indxx,
					     17, 0.0, Constants::pi);
      Bphi[inch][ipt] =
	histogramFactory().createHistogram1D(dir + "/Bphi" + indxx,
					     17, 0.0, Constants::pi);
      S2phi[inch][ipt] =
	histogramFactory().createHistogram2D(dir + "/S2phi" + indxx,
					     7, 2.0, 2.0*etamax,
					     17, 0.0, Constants::pi);
      B2phi[inch][ipt] =
	histogramFactory().createHistogram2D(dir + "/B2phi" + indxx,
					     7, 2.0, 2.0*etamax,
					     17, 0.0, Constants::pi);
      S3phi[inch][ipt] =
	histogramFactory().createHistogram1D(dir + "/S3phi" + indxx,
					     17, 0.0, Constants::pi);
      B3phi[inch][ipt] =
	histogramFactory().createHistogram1D(dir + "/B3phi" + indxx,
					     17, 0.0, Constants::pi);
      multbin[inch][ipt] =
	histogramFactory().createHistogram1D(dir + "/multbin" + indxx,
					     250, -0.5, 249.5);
      multbin2[inch][ipt] =
	histogramFactory().createHistogram1D(dir + "/multbin2" + indxx,
					     248, 1.5, 249.5);
   }
  }
  SMB01 =
    histogramFactory().createHistogram2D
    (dir + "/SMB01",
     33, -2.0*etamax, 2.0*etamax,
     30, -0.5*Constants::pi,1.5*Constants::pi);
  BMB01 =
    histogramFactory().createHistogram2D
    (dir + "/BMB01",
     33, -2.0*etamax, 2.0*etamax,
     30, -0.5*Constants::pi,1.5*Constants::pi);
  SMB13 =
    histogramFactory().createHistogram2D
    (dir + "/SMB13",
     33, -2.0*etamax, 2.0*etamax,
     30, -0.5*Constants::pi,1.5*Constants::pi);
  BMB13 =
    histogramFactory().createHistogram2D
    (dir + "/BMB13",
     33, -2.0*etamax, 2.0*etamax,
     30, -0.5*Constants::pi,1.5*Constants::pi);
  SHN01 =
    histogramFactory().createHistogram2D
    (dir + "/SHN01",
     33, -2.0*etamax, 2.0*etamax,
     30, -0.5*Constants::pi,1.5*Constants::pi);
  BHN01 =
    histogramFactory().createHistogram2D
    (dir + "/BHN01",
     33, -2.0*etamax, 2.0*etamax,
     30, -0.5*Constants::pi,1.5*Constants::pi);
  SHN13 =
    histogramFactory().createHistogram2D
    (dir + "/SHN13",
     33, -2.0*etamax, 2.0*etamax,
     30, -0.5*Constants::pi,1.5*Constants::pi);
  BHN13 =
    histogramFactory().createHistogram2D
    (dir + "/BHN13",
     33, -2.0*etamax, 2.0*etamax,
     30, -0.5*Constants::pi,1.5*Constants::pi);

  S2MB01 =
    histogramFactory().createHistogram2D
    (dir + "/S2MB01",
     33, -2.0*etamax, 2.0*etamax,
     30, -0.5*Constants::pi,1.5*Constants::pi);
  B2MB01 =
    histogramFactory().createHistogram2D
    (dir + "/B2MB01",
     33, -2.0*etamax, 2.0*etamax,
     30, -0.5*Constants::pi,1.5*Constants::pi);
  multMB01 = 
    histogramFactory().createHistogram1D
    (dir + "/multMB01", 248, 1.5, 249.5);
  S2MB13 =
    histogramFactory().createHistogram2D
    (dir + "/S2MB13",
     33, -2.0*etamax, 2.0*etamax,
     30, -0.5*Constants::pi,1.5*Constants::pi);
  B2MB13 =
    histogramFactory().createHistogram2D
    (dir + "/B2MB13",
     33, -2.0*etamax, 2.0*etamax,
     30, -0.5*Constants::pi,1.5*Constants::pi);
  multMB13 = 
    histogramFactory().createHistogram1D
    (dir + "/multMB13", 248, 1.5, 249.5);
  S2HN01 =
    histogramFactory().createHistogram2D
    (dir + "/S2HN01",
     33, -2.0*etamax, 2.0*etamax,
     30, -0.5*Constants::pi,1.5*Constants::pi);
  B2HN01 =
    histogramFactory().createHistogram2D
    (dir + "/B2HN01",
     33, -2.0*etamax, 2.0*etamax,
     30, -0.5*Constants::pi,1.5*Constants::pi);
  multHN01 = 
    histogramFactory().createHistogram1D
    (dir + "/multHN01", 248, 1.5, 249.5);
  S2HN13 =
    histogramFactory().createHistogram2D
    (dir + "/S2HN13",
     33, -2.0*etamax, 2.0*etamax,
     30, -0.5*Constants::pi,1.5*Constants::pi);
  B2HN13 =
    histogramFactory().createHistogram2D
    (dir + "/B2HN13",
     33, -2.0*etamax, 2.0*etamax,
     30, -0.5*Constants::pi,1.5*Constants::pi);
  multHN13 = 
    histogramFactory().createHistogram1D
    (dir + "/multHN13", 248, 1.5, 249.5);

  mult = histogramFactory().createHistogram1D(dir + "/Nch",
					     200, -0.5, 199.5);
  refphi = vector< vector< vector<double> > >(4);
  refeta = vector< vector< vector<double> > >(4);
  refNch04 = vector<int>(4);
  refweight = vector<double>(4);
  sumrefweight = vector<double>(4);
  sumweight = vector<double>(4);

  if ( !externalFile.empty() ) {
    string::size_type n = externalFile.find("%N");
    if ( n == string::npos )
      file.open(externalFile);
    else {
      ostringstream os;
      os << generator()->N();
      file.open(externalFile.replace(n, 2, os.str()));
    }
      
  }

}

void CMSlrnsac::persistentOutput(PersistentOStream & os) const {
  os << doRweight << externalFile;
}

void CMSlrnsac::persistentInput(PersistentIStream & is, int) {
  is >> doRweight >> externalFile;
}

// ClassDescription<CMSlrnsac> CMSlrnsac::initCMSlrnsac;
// Definition of the static class description member.

DescribeClass<CMSlrnsac,AnalysisHandler>
describeCMSlrnsac("ThePEG::CMSlrnsac", "CMSlrnsac.so");

void CMSlrnsac::Init() {

  static ClassDocumentation<CMSlrnsac> documentation
    ("There is no documentation for the CMSlrnsac class");


  static Switch<CMSlrnsac,int> interfaceRweight
    ("Rweight",
     "How to weight the background samples",
     &CMSlrnsac::doRweight, 2, true, false);
  static SwitchOption interfaceRweightSameAsSignal
    (interfaceRweight,
     "SameAsSignal",
     "Fill with the same weight as the signal irrespectively of the weight "
     "of the reference event.",
     0);
  static SwitchOption interfaceRweightSqrtOfWeights
    (interfaceRweight,
     "SqrtOfWeights",
     "Fill with the square root of the weight of the signal event times the "
     "weight of the reference event.",
     1);
  static SwitchOption interfaceRweightProdOfWeights
    (interfaceRweight,
     "ProdOfWeights",
     "Fill with the weight of the signal event times the weight "
     "of the reference event.",
     2);

  static Parameter<CMSlrnsac,string> interfaceFileName
    ("FileName",
     "If not empty, all events generated with ThePEG are discarded "
     "and events are read from the given file instead (free format, "
     "but all events must start with a line containing "
     "'*** BEGIN EVENT ***' followed by a line containing the number "
     "of particles to follow and the event weight, then followed by the "
     "given number of lines with final state particles represented by "
     "their PDG id and 5-momentum (space separated x, y,z, e and m "
     "components). If the string ends with a <code>|</code> the "
     "preceeding string is interpreted as a command, the output of which "
     "will be read through a pipe. Any occurence of '%N' in the name will be "
     "replaced by the number of events to be read.",
     &CMSlrnsac::externalFile, "", false, false);

}

