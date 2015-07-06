// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the FSAnalysis class.
//

#include "FSAnalysis.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "DipoleEventHandler.h"
#include "EventFiller.h"
#include "DiffractiveEventFiller.h"
#include "ThePEG/Utilities/Current.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#ifndef LWH_AIAnalysisFactory_H
#ifndef LWH 
#define LWH ThePEGLWH
#endif
#include "ThePEG/Analysis/LWH/AnalysisFactory.h"
#endif

using namespace DIPSY;

FSAnalysis::FSAnalysis() {}

FSAnalysis::~FSAnalysis() {}

pair<double, double> FSAnalysis::findCoG(vector < vector<double> > gluons) {
  if ( gluons.empty() ) return make_pair(0.0, 0.0);
  pair<double, double> sum = make_pair(0.0, 0.0);
  for ( int i = 0; i < int(gluons.size()); i++) {
    sum.first += gluons[i][1];
    sum.second += gluons[i][2];
  }
  sum.first /= double(gluons.size());
  sum.second /= double(gluons.size());
  return sum;
}

double FSAnalysis::eccentricity(int n, vector < vector<double> > gluons, pair<double, double> CoG) {
  return eccentricity(n, gluons, CoG, true, false);
}

double FSAnalysis::eccentricityN(int n, vector < vector<double> > gluons, pair<double, double> CoG) {
  return eccentricity(n, gluons, CoG, false, false);
}

double FSAnalysis::eccentricityW(int n, vector < vector<double> > gluons, pair<double, double> CoG) {
  return eccentricity(n, gluons, CoG, true, true);
}

double FSAnalysis::eccentricity(int n, vector < vector<double> > gluons, pair<double, double> CoG, bool always2, bool weighted) {
  if ( n < 1 ) return 0.0;
  if ( int(gluons.size()) < 2 ) return 0.0;
  double avRnCos = 0.0; // <r^n cos(n phi)>
  double avRnSin = 0.0; // <r^n sin(n phi)>
  double avRn = 0.0; // <r^n>
  int power = always2 ? 2:n;
  if ( n == 1 ) power = 3;
  for ( int i = 0; i < int(gluons.size()); i++) {
    pair<double, double> p = make_pair(gluons[i][1] - CoG.first, gluons[i][2] - CoG.second);
    double r = sqrt(sqr(p.first) + sqr(p.second));
    double rn = 1.0;
    for ( int j = 0; j < power; j++ ) rn *= r;
    if ( weighted ) rn *= max(1.0,log(sqrt(sqr(gluons[i][3]) + sqr(gluons[i][4])))+1.0);
    double phi = atan2(p.second, p.first);
    avRnCos += rn*cos(n*phi);
    avRnSin += rn*sin(n*phi);
    avRn += rn;
  }
  avRnCos /= double(gluons.size());
  avRnSin /= double(gluons.size());
  avRn /= double(gluons.size());
  return sqrt(sqr(avRnCos) + sqr(avRnSin))/avRn;
}

double FSAnalysis::participantAngle(int n, vector < vector<double> > gluons, pair<double, double> CoG) {
  return participantAngle(n, gluons, CoG, true, false);
}

double FSAnalysis::participantAngleN(int n, vector < vector<double> > gluons, pair<double, double> CoG) {
  return participantAngle(n, gluons, CoG, false, false);
}

double FSAnalysis::participantAngleW(int n, vector < vector<double> > gluons, pair<double, double> CoG) {
  return participantAngle(n, gluons, CoG, true, true);
}

double FSAnalysis::participantAngle(int n, vector < vector<double> > gluons, pair<double, double> CoG, bool always2, bool weighted) {
  if ( n < 1 ) return 0.0;
  if ( int(gluons.size()) < 2 ) return 0.0;
  double avRnCos = 0.0; // <r^n cos(n phi)>
  double avRnSin = 0.0; // <r^n sin(n phi)>
  int power = always2 ? 2:n;
  if ( n == 1 ) power = 3;
  for ( int i = 0; i < int(gluons.size()); i++) {
    pair<double, double> p = make_pair(gluons[i][1] - CoG.first, gluons[i][2] - CoG.second);
    double r = sqrt(sqr(p.first) + sqr(p.second));
    double rn = 1.0;
    for ( int j = 0; j < power; j++ ) rn *= r;
    if ( weighted ) rn *= max(1.0,log(sqrt(sqr(gluons[i][3]) + sqr(gluons[i][4])))+1.0);
    double phi = atan2(p.second, p.first);
    avRnCos += rn*cos(n*phi);
    avRnSin += rn*sin(n*phi);
  }
  avRnCos /= double(gluons.size());
  avRnSin /= double(gluons.size());
  return atan2(avRnSin, avRnCos)/double(n);
}

vector<double> FSAnalysis::getValues(char line[256]) {
  vector<double> ret;
  for( int i = 0; line[i] != '\0'; i++) {
    char value[256] = "";
    bool found = false;
    for(int j = 0; line[i]!=' '; i++, j++) {
      found = true;
      value[j] = line[i];
    }
    if ( found ) {
      ret.push_back(atof(value));
    }
  }
  return ret;
}

pair<double, double> FSAnalysis::getB(char line[256]) {
  pair<double, double> ret;
  int words = 0;
  for( int i = 0; line[i] != '\0'; i++) {
    char value[256] = "";
    bool found = false;
    for(int j = 0; line[i]!=' '; i++, j++) {
      found = true;
      value[j] = line[i];
    }
    if ( found ) {
      words++;
      if ( words == 2 ) ret.first = atof(value);
      if ( words == 3 ) ret.second = atof(value);
      if ( words > 2 ) return ret;
    }
  }
  return ret;
}

double FSAnalysis::getWeight(char line[256]) {
  double ret = 0.0;
  int words = 0;
  for( int i = 0; line[i] != '\0'; i++) {
    char value[256] = "";
    bool found = false;
    for(int j = 0; line[i]!=' '; i++, j++) {
      found = true;
      value[j] = line[i];
    }
    if ( found ) {
      words++;
      if ( words == 2 ) ret = atof(value);
      if ( words > 1 ) return ret;
    }
  }
  return ret;
}

int FSAnalysis::getNSpectators(char line[256]) {
  double ret = 0;
  int words = 0;
  for( int i = 0; line[i] != '\0'; i++) {
    char value[256] = "";
    bool found = false;
    for(int j = 0; line[i]!=' '; i++, j++) {
      found = true;
      value[j] = line[i];
    }
    if ( found ) {
      words++;
      if ( words == 5 ) ret = atoi(value);
      if ( words > 4 ) return ret;
    }
  }
  return ret;
}

vector < vector<double> > FSAnalysis::extractEvent(string filename) {
  vector < vector<double> > ret;
  ifstream infile(filename.c_str());
  if ( !infile ) {
    cout << "file " << filename << " didnt exist!! :o" << endl;
    getchar();
    return ret;
  }
  char line[256] = "";
  for ( int i = 0; i < 11; i++ ) infile.getline(line,256);
  while( infile.good() ) {
    infile.getline(line,256);
    if ( infile.gcount() < 5 ) continue;
    vector<double> values = getValues(line);
    ret.push_back(values);
  }
  infile.close();
  return ret;
}

pair<double, double> FSAnalysis::extractB(string filename) {
  pair<double, double> ret;
  ifstream infile(filename.c_str());
  if ( !infile ) {
    cout << "file " << filename << " didnt exist!! :o" << endl;
    getchar();
    return ret;
  }
  char line[256] = "";
  for ( int i = 0; i < 5; i++ ) infile.getline(line,256);
  if( infile.good() ) {
    infile.getline(line,256);
    if ( infile.gcount() < 5 ) return ret;
    ret = getB(line);
  }
  infile.close();
  return ret;
}

double FSAnalysis::extractWeight(string filename) {
  double ret = 0.0;
  ifstream infile(filename.c_str());
  if ( !infile ) {
    cout << "file " << filename << " didnt exist!! :o" << endl;
    getchar();
    return ret;
  }
  char line[256] = "";
  for ( int i = 0; i < 3; i++ ) infile.getline(line,256);
  if( infile.good() ) {
    infile.getline(line,256);
    if ( infile.gcount() < 5 ) return ret;
    ret = getWeight(line);
  }
  infile.close();
  return ret;
}

int FSAnalysis::extractNSpectators(string filename) {
  int ret = 0;
  ifstream infile(filename.c_str());
  if ( !infile ) {
    cout << "file " << filename << " didnt exist!! :o" << endl;
    getchar();
    return ret;
  }
  char line[256] = "";
  for ( int i = 0; i < 7; i++ ) infile.getline(line,256);
  if( infile.good() ) {
    infile.getline(line,256);
    if ( infile.gcount() < 5 ) return ret;
    ret = getNSpectators(line);
  }
  infile.close();
  return ret;
}

void FSAnalysis::changeWeight(string filename, double factor) {
  // double ret = 0.0;
  // fstream infile(filename.c_str());
  // cout << "changeing weight in file " << filename << endl;
  // if ( !infile ) {
  //   cout << "file " << filename << " didnt exist!! :o" << endl;
  //   getchar();
  //   return;
  // }
  // char line[256] = "";
  // for ( int i = 0; i < 3; i++ ) infile.getline(line,256);
  // if( infile.good() ) {
  //   infile.getline(line,256);
  //   if ( infile.gcount() < 5 ) return;
  //   ret = getWeight(line);
  //   cout << "found weight " << ret;
  //   infile << "this is inserted text" << endl;
  //   infile.write("this is inserted text with put", 20);
  //   infile.flush();
  // }
  // infile.close();
  // getchar();
  // return;
}

void FSAnalysis::eccentricities(tEventPtr event) {

  //extract the information from file
  ostringstream filename;
  //this determined which directories are analysed.
  //for example "./LHCPbPb/LHCPbPbB" can be used to
  //run the same analysis on LHC lead-lead
  string name = generator()->filename();
  filename << "/scratch/parton/christof/DIPSY/events/" << name << "B";
  int Bbin = int ((event->number()-1)/1000);
  int eventNumber = ((event->number()-1) % 1000) + 1;
    
  filename << Bbin << Bbin+1 << "/event" << eventNumber << ".dat";

  if ( int(event->number()/100) == double(event->number())/100.0 )
    cout << "Analysing event " << filename.str() << endl;
  vector < vector<double> > gluons = extractEvent(filename.str());

  pair<double, double> Bvec = extractB(filename.str());
  double B = sqrt(sqr(Bvec.first) + sqr(Bvec.second));
  double W = extractWeight(filename.str());
  int nSpectators = extractNSpectators(filename.str());
  int nGlue = gluons.size();

  //find the gluons in rapidity [-1,1]
  vector<vector<double> > centralGluons;
  int centralNGlue = 0;
  double centralWeightedNGlue = 0.0;
  for ( int i = 0; i < nGlue; i++ ) {
    if ( abs(gluons[i][5]) < 1.0 ) {
      centralNGlue++;
      double pt2 = sqr(gluons[i][3]) + sqr(gluons[i][4]);
      centralWeightedNGlue += log(max(1.0, pt2));
      centralGluons.push_back(gluons[i]);
    }
  }

  //find all gluons that do not have the valence pt only.
  vector<vector<double> > touched;
  for ( int i = 0; i < nGlue; i++ ) {
    double pt2 = sqr(gluons[i][3]) + sqr(gluons[i][4]);
    if ( pt2 > sqr(0.315) ||
	 pt2 < sqr(0.3) ) {
      touched.push_back(gluons[i]);
    }
  }

  //find the gluons in [-3,-1]
  vector < vector<double> > participantsm3tom1;
  for ( int i = 0; i < nGlue; i++ ) {
    if ( abs(gluons[i][5] + 2.0) < 1.0 ) {
      participantsm3tom1.push_back(gluons[i]);
    }
  }

  //find the gluons in [1,3]
  vector < vector<double> > participantsp1top3;
  for ( int i = 0; i < nGlue; i++ ) {
    if ( abs(gluons[i][5] - 2.0) < 1.0 ) {
      participantsp1top3.push_back(gluons[i]);
    }
  }


  histNGlue->fill(B, centralNGlue*W);
  histNGlueSqr->fill(B, sqr(centralNGlue)*W);
  histNGlueWeights->fill(B, W);

  histNGlueSpect->fill(nSpectators, centralNGlue*W);
  histNGlueSpectSqr->fill(nSpectators, sqr(centralNGlue)*W);
  histNGlueSpectWeights->fill(nSpectators, W);

  histNSpectator->fill(B, nSpectators*W);
  histNSpectatorSqr->fill(B, sqr(nSpectators)*W);
  histNSpectatorWeights->fill(B, W);

  centralGlueMult += nGlue*W;

  //analysis with usual central participants
  pair<double, double> CoG = findCoG(centralGluons);

  double Psi = atan2(Bvec.second, Bvec.first);

  totalWeight += W;
  if ( B < 3 ) centralWeight += W;
  if ( B > 5 && B < 9 ) midWeight += W;

  double e1N = eccentricityN(1, centralGluons, CoG);
  double e1W = eccentricityW(1, centralGluons, CoG);
  double e1 = eccentricity(1, centralGluons, CoG);
  histe1->fill(nSpectators, e1*W);
  histe1Weights->fill(nSpectators, W);
  histe1N->fill(nSpectators, e1N*W);
  histe1NWeights->fill(nSpectators, W);
  histe1W->fill(nSpectators, e1W*W);
  histe1WWeights->fill(nSpectators, W);
  double phi1N = participantAngleN(1, centralGluons, CoG);
  double phi1W = participantAngleW(1, centralGluons, CoG);
  double phi1 = participantAngle(1, centralGluons, CoG);
  double phiPsi1N = M_PI/2.0 - abs(fmod(M_PI - abs(M_PI - abs(phi1N - Psi)),M_PI/1.0)-M_PI/2.0);
  double phiPsi1W = M_PI/2.0 - abs(fmod(M_PI - abs(M_PI - abs(phi1W - Psi)),M_PI/1.0)-M_PI/2.0);
  double phiPsi1 = M_PI/2.0 - abs(fmod(M_PI - abs(M_PI - abs(phi1 - Psi)),M_PI/1.0)-M_PI/2.0);
  if ( nSpectators < 294 ) {
    histphi1->fill(phiPsi1, W);
    if ( nSpectators < 44 ) histphi1central->fill(phiPsi1, W);
    if ( nSpectators > 143 ) histphi1mid->fill(phiPsi1, W);
    histphi1N->fill(phiPsi1N, W);
    histphi1W->fill(phiPsi1W, W);
  }

  double e2N = eccentricityN(2, centralGluons, CoG);
  double e2W = eccentricityW(2, centralGluons, CoG);
  double e2 = eccentricity(2, centralGluons, CoG);
  histe2->fill(nSpectators, e2*W);
  histe2Weights->fill(nSpectators, W);
  histe2N->fill(nSpectators, e2N*W);
  histe2NWeights->fill(nSpectators, W);
  histe2W->fill(nSpectators, e2W*W);
  histe2WWeights->fill(nSpectators, W);
  double phi2N = participantAngleN(2, centralGluons, CoG);
  double phi2W = participantAngleW(2, centralGluons, CoG);
  double phi2 = participantAngle(2, centralGluons, CoG);
  double phiPsi2N = M_PI/2.0 - abs(fmod(M_PI - abs(M_PI - abs(phi2N - Psi)),M_PI/1.0)-M_PI/2.0);
  double phiPsi2W = M_PI/2.0 - abs(fmod(M_PI - abs(M_PI - abs(phi2W - Psi)),M_PI/1.0)-M_PI/2.0);
  double phiPsi2 = M_PI/2.0 - abs(fmod(M_PI - abs(M_PI - abs(phi2 - Psi)),M_PI/1.0)-M_PI/2.0);
  if ( nSpectators < 294 ) {
    histphi2->fill(phiPsi2, W);
    if ( nSpectators < 44 ) histphi2central->fill(phiPsi2, W);
    if ( nSpectators > 143 ) histphi2mid->fill(phiPsi2, W);
    histphi2N->fill(phiPsi2N, W);
    histphi2W->fill(phiPsi2W, W);
  }

  double e3N = eccentricityN(3, centralGluons, CoG);
  double e3W = eccentricityW(3, centralGluons, CoG);
  double e3 = eccentricity(3, centralGluons, CoG);
  histe3->fill(nSpectators, e3*W);
  histe3Weights->fill(nSpectators, W);
  histe3N->fill(nSpectators, e3N*W);
  histe3NWeights->fill(nSpectators, W);
  histe3W->fill(nSpectators, e3W*W);
  histe3WWeights->fill(nSpectators, W);
  double phi3N = participantAngleN(3, centralGluons, CoG);
  double phi3W = participantAngleW(3, centralGluons, CoG);
  double phi3 = participantAngle(3, centralGluons, CoG);
  double phiPsi3N = M_PI/3.0 - abs(fmod(M_PI - abs(M_PI - abs(phi3N - Psi)),M_PI/1.5)-M_PI/3.0);
  double phiPsi3W = M_PI/3.0 - abs(fmod(M_PI - abs(M_PI - abs(phi3W - Psi)),M_PI/1.5)-M_PI/3.0);
  double phiPsi3 = M_PI/3.0 - abs(fmod(M_PI - abs(M_PI - abs(phi3 - Psi)),M_PI/1.5)-M_PI/3.0);
  if ( nSpectators < 294 ) {
    histphi3->fill(phiPsi3, W);
    if ( nSpectators < 44 ) histphi3central->fill(phiPsi3, W);
    if ( nSpectators > 143 ) histphi3mid->fill(phiPsi3, W);
    histphi3N->fill(phiPsi3N, W);
    histphi3W->fill(phiPsi3W, W);
  }

  double e4N = eccentricityN(4, centralGluons, CoG);
  double e4W = eccentricityW(4, centralGluons, CoG);
  double e4 = eccentricity(4, centralGluons, CoG);
  histe4->fill(nSpectators, e4*W);
  histe4Weights->fill(nSpectators, W);
  histe4N->fill(nSpectators, e4N*W);
  histe4NWeights->fill(nSpectators, W);
  histe4W->fill(nSpectators, e4W*W);
  histe4WWeights->fill(nSpectators, W);
  double phi4N = participantAngleN(4, centralGluons, CoG);
  double phi4W = participantAngleW(4, centralGluons, CoG);
  double phi4 = participantAngle(4, centralGluons, CoG);
  double phiPsi4N = M_PI/4.0 - abs(fmod(M_PI - abs(M_PI - abs(phi4N - Psi)),M_PI/2.)-M_PI/4.0);
  double phiPsi4W = M_PI/4.0 - abs(fmod(M_PI - abs(M_PI - abs(phi4W - Psi)),M_PI/2.)-M_PI/4.0);
  double phiPsi4 = M_PI/4.0 - abs(fmod(M_PI - abs(M_PI - abs(phi4 - Psi)),M_PI/2.)-M_PI/4.0);
  if ( nSpectators < 294 ) {
    histphi4->fill(phiPsi4, W);
    if ( nSpectators < 44 ) histphi4central->fill(phiPsi4, W);
    if ( nSpectators > 143 ) histphi4mid->fill(phiPsi4, W);
    histphi4N->fill(phiPsi4N, W);
    histphi4W->fill(phiPsi4W, W);
  }

  double e5N = eccentricityN(5, centralGluons, CoG);
  double e5W = eccentricityW(5, centralGluons, CoG);
  double e5 = eccentricity(5, centralGluons, CoG);
  histe5->fill(nSpectators, e5*W);
  histe5Weights->fill(nSpectators, W);
  histe5N->fill(nSpectators, e5N*W);
  histe5NWeights->fill(nSpectators, W);
  histe5W->fill(nSpectators, e5W*W);
  histe5WWeights->fill(nSpectators, W);
  double phi5N = participantAngleN(5, centralGluons, CoG);
  double phi5W = participantAngleW(5, centralGluons, CoG);
  double phi5 = participantAngle(5, centralGluons, CoG);
  double phiPsi5N = M_PI/5.0 - abs(fmod(M_PI - abs(M_PI - abs(phi5N - Psi)),M_PI/2.5)-M_PI/5.0);
  double phiPsi5W = M_PI/5.0 - abs(fmod(M_PI - abs(M_PI - abs(phi5W - Psi)),M_PI/2.5)-M_PI/5.0);
  double phiPsi5 = M_PI/5.0 - abs(fmod(M_PI - abs(M_PI - abs(phi5 - Psi)),M_PI/2.5)-M_PI/5.0);
  if ( nSpectators < 294 ) {
    histphi5->fill(phiPsi5, W);
    if ( nSpectators < 44 ) histphi5central->fill(phiPsi5, W);
    if ( nSpectators > 143 ) histphi5mid->fill(phiPsi5, W);
    histphi5N->fill(phiPsi5N, W);
    histphi5W->fill(phiPsi5W, W);
  }

  if ( centralGluons.size() > 2 ) {
    Area area = averageOverlapArea(centralGluons);
    histArea->fill(nSpectators, area/sqr(femtometer)*W);
    histAreaWeights->fill(nSpectators, W);
  }

  histe1Sqrd->fill(nSpectators, sqr(e1)*W);
  histe1SqrdWeights->fill(nSpectators, W);
  histe2Sqrd->fill(nSpectators, sqr(e2)*W);
  histe2SqrdWeights->fill(nSpectators, W);
  histe3Sqrd->fill(nSpectators, sqr(e3)*W);
  histe3SqrdWeights->fill(nSpectators, W);
  histe4Sqrd->fill(nSpectators, sqr(e4)*W);
  histe4SqrdWeights->fill(nSpectators, W);
  histe5Sqrd->fill(nSpectators, sqr(e5)*W);
  histe5SqrdWeights->fill(nSpectators, W);

  histe1Quad->fill(nSpectators, sqr(sqr(e1))*W);
  histe1QuadWeights->fill(nSpectators, W);
  histe2Quad->fill(nSpectators, sqr(sqr(e2))*W);
  histe2QuadWeights->fill(nSpectators, W);
  histe3Quad->fill(nSpectators, sqr(sqr(e3))*W);
  histe3QuadWeights->fill(nSpectators, W);
  histe4Quad->fill(nSpectators, sqr(sqr(e4))*W);
  histe4QuadWeights->fill(nSpectators, W);
  histe5Quad->fill(nSpectators, sqr(sqr(e5))*W);
  histe5QuadWeights->fill(nSpectators, W);

  //all touched partons, used for the non-evo case
  pair<double, double> CoGTouched = findCoG(touched);

  double e2Touched = eccentricity(2, touched, CoGTouched);
  histe2Touched->fill(nSpectators, e2Touched*W);
  histe2TouchedWeights->fill(nSpectators, W);
  double e3Touched = eccentricity(3, touched, CoGTouched);
  histe3Touched->fill(nSpectators, e3Touched*W);
  histe3TouchedWeights->fill(nSpectators, W);
  double e4Touched = eccentricity(4, touched, CoGTouched);
  histe4Touched->fill(nSpectators, e4Touched*W);
  histe4TouchedWeights->fill(nSpectators, W);
  double e5Touched = eccentricity(5, touched, CoGTouched);
  histe5Touched->fill(nSpectators, e5Touched*W);
  histe5TouchedWeights->fill(nSpectators, W);

  //all (-3,-1) partons
  pair<double, double> CoGm3tom1 = findCoG(participantsm3tom1);

  double e1m3tom1 = eccentricity(1, participantsm3tom1, CoGm3tom1);
  histe1m3tom1->fill(nSpectators, e1m3tom1*W);
  histe1m3tom1Weights->fill(nSpectators, W);
  double e2m3tom1 = eccentricity(2, participantsm3tom1, CoGm3tom1);
  histe2m3tom1->fill(nSpectators, e2m3tom1*W);
  histe2m3tom1Weights->fill(nSpectators, W);
  double e3m3tom1 = eccentricity(3, participantsm3tom1, CoGm3tom1);
  histe3m3tom1->fill(nSpectators, e3m3tom1*W);
  histe3m3tom1Weights->fill(nSpectators, W);
  double e4m3tom1 = eccentricity(4, participantsm3tom1, CoGm3tom1);
  histe4m3tom1->fill(nSpectators, e4m3tom1*W);
  histe4m3tom1Weights->fill(nSpectators, W);
  double e5m3tom1 = eccentricity(5, participantsm3tom1, CoGm3tom1);
  histe5m3tom1->fill(nSpectators, e5m3tom1*W);
  histe5m3tom1Weights->fill(nSpectators, W);

  //all (1,3) partons
  pair<double, double> CoGp1top3 = findCoG(participantsp1top3);

  double e1p1top3 = eccentricity(1, participantsp1top3, CoGp1top3);
  histe1p1top3->fill(nSpectators, e1p1top3*W);
  histe1p1top3Weights->fill(nSpectators, W);
  double e2p1top3 = eccentricity(2, participantsp1top3, CoGp1top3);
  histe2p1top3->fill(nSpectators, e2p1top3*W);
  histe2p1top3Weights->fill(nSpectators, W);
  double e3p1top3 = eccentricity(3, participantsp1top3, CoGp1top3);
  histe3p1top3->fill(nSpectators, e3p1top3*W);
  histe3p1top3Weights->fill(nSpectators, W);
  double e4p1top3 = eccentricity(4, participantsp1top3, CoGp1top3);
  histe4p1top3->fill(nSpectators, e4p1top3*W);
  histe4p1top3Weights->fill(nSpectators, W);
  double e5p1top3 = eccentricity(5, participantsp1top3, CoGp1top3);
  histe5p1top3->fill(nSpectators, e5p1top3*W);
  histe5p1top3Weights->fill(nSpectators, W);

  //CoG distance
  double CoGdist = sqrt(sqr(CoGm3tom1.first - CoGp1top3.first) +
			sqr(CoGm3tom1.second - CoGp1top3.second));
  histCoGdistance->fill(CoGdist, W);


  //correaltions
  if ( nSpectators < 294 ) {
    over100Weight += W;
    histeCorr1->fill(e1m3tom1, e1p1top3, W);
    histeCorr2->fill(e2m3tom1, e2p1top3, W);
    histeCorr3->fill(e3m3tom1, e3p1top3, W);
    histeCorr4->fill(e4m3tom1, e4p1top3, W);
    histeCorr5->fill(e5m3tom1, e5p1top3, W);
  }

  if ( nSpectators < 294 ) {
    epFepB1 += e1m3tom1*e1p1top3*W;
    epF1 += e1p1top3*W;
    epB1 += e1m3tom1*W;
    epSqrdF1 += sqr(e1p1top3)*W;
    epSqrdB1 += sqr(e1m3tom1)*W;
    epFepB2 += e2m3tom1*e2p1top3*W;
    epF2 += e2p1top3*W;
    epB2 += e2m3tom1*W;
    epSqrdF2 += sqr(e2p1top3)*W;
    epSqrdB2 += sqr(e2m3tom1)*W;
    if ( nSpectators < 44 ) {
      over350Weight += W;
      epFepB2central += e2m3tom1*e2p1top3*W;
      epF2central += e2p1top3*W;
      epB2central += e2m3tom1*W;
      epSqrdF2central += sqr(e2p1top3)*W;
      epSqrdB2central += sqr(e2m3tom1)*W;
    }
    if ( nSpectators > 143 ) {
      over100under250Weight += W;
      epFepB2mid += e2m3tom1*e2p1top3*W;
      epF2mid += e2p1top3*W;
      epB2mid += e2m3tom1*W;
      epSqrdF2mid += sqr(e2p1top3)*W;
      epSqrdB2mid += sqr(e2m3tom1)*W;
    }
    epFepB3 += e3m3tom1*e3p1top3*W;
    epF3 += e3p1top3*W;
    epB3 += e3m3tom1*W;
    epSqrdF3 += sqr(e3p1top3)*W;
    epSqrdB3 += sqr(e3m3tom1)*W;
    epFepB4 += e4m3tom1*e4p1top3*W;
    epF4 += e4p1top3*W;
    epB4 += e4m3tom1*W;
    epSqrdF4 += sqr(e4p1top3)*W;
    epSqrdB4 += sqr(e4m3tom1)*W;
  }

  histe1e1p1top3->fill(nSpectators, e1m3tom1*e1p1top3*W);
  histe1e1p1top3Weights->fill(nSpectators, W);
  histe2e2p1top3->fill(nSpectators, e2m3tom1*e2p1top3*W);
  histe2e2p1top3Weights->fill(nSpectators, W);
  histe3e3p1top3->fill(nSpectators, e3m3tom1*e3p1top3*W);
  histe3e3p1top3Weights->fill(nSpectators, W);
  histe4e4p1top3->fill(nSpectators, e4m3tom1*e4p1top3*W);
  histe4e4p1top3Weights->fill(nSpectators, W);
  histe5e5p1top3->fill(nSpectators, e5m3tom1*e5p1top3*W);
  histe5e5p1top3Weights->fill(nSpectators, W);

  if ( nSpectators < 294 ) {
    double phi1m = participantAngle(1, participantsm3tom1, CoGm3tom1);
    double phi1p = participantAngle(1, participantsp1top3, CoGp1top3);
    double phiDiff1 = M_PI/1.0 - abs(fmod(abs(phi1m - phi1p), 2.0*M_PI/1.0) - M_PI/1.0);
    histphiCorr1->fill(phiDiff1, W);
    double phi2m = participantAngle(2, participantsm3tom1, CoGm3tom1);
    double phi2p = participantAngle(2, participantsp1top3, CoGp1top3);
    double phiDiff2 = M_PI/2.0 - abs(fmod(abs(phi2m - phi2p), 2.0*M_PI/2.0) - M_PI/2.0);
    histphiCorr2->fill(phiDiff2, W);
    if ( nSpectators < 44 )   histphiCorr2central->fill(phiDiff2, W);
    if ( nSpectators > 143 )   histphiCorr2mid->fill(phiDiff2, W);
    double phi3m = participantAngle(3, participantsm3tom1, CoGm3tom1);
    double phi3p = participantAngle(3, participantsp1top3, CoGp1top3);
    double phiDiff3 = M_PI/3.0 - abs(fmod(abs(phi3m - phi3p), 2.0*M_PI/3.0) - M_PI/3.0);
    histphiCorr3->fill(phiDiff3, W);
    double phi4m = participantAngle(4, participantsm3tom1, CoGm3tom1);
    double phi4p = participantAngle(4, participantsp1top3, CoGp1top3);
    double phiDiff4 = M_PI/4.0 - abs(fmod(abs(phi4m - phi4p), 2.0*M_PI/4.0) - M_PI/4.0);
    histphiCorr4->fill(phiDiff4, W);
    double phi5m = participantAngle(5, participantsm3tom1, CoGm3tom1);
    double phi5p = participantAngle(5, participantsp1top3, CoGp1top3);
    double phiDiff5 = M_PI/5.0 - abs(fmod(abs(phi5m - phi5p), 2.0*M_PI/5.0) - M_PI/5.0);
    histphiCorr5->fill(phiDiff5, W);
  }

  //dN/dy in centrality classes
  for ( int i = 0; i < nGlue; i++ ) {
    histdNdnSpectdeta->fill(nSpectators, gluons[i][5], W);
    if (nSpectators < 100)  histdNdeta0to100spect->fill(gluons[i][5], W);
    else if (nSpectators < 200)  histdNdeta100to200spect->fill(gluons[i][5], W);
    else if (nSpectators < 300)  histdNdeta200to300spect->fill(gluons[i][5], W);
    else if (nSpectators < 400)  histdNdeta300to400spect->fill(gluons[i][5], W);
  }
  if (nSpectators < 100)  spect0to100Weight += W;
  else if (nSpectators < 200)  spect100to200Weight += W;
  else if (nSpectators < 300)  spect200to300Weight += W;
  else if (nSpectators < 400)  spect300to400Weight += W;

  //flux tube analysis
  //sort the particles in [-3,-2] and [2,3] in bins of r (using CoGTouched).
  vector< vector < vector<double> > > participantsFluxMinus, participantsFluxPlus;
  for ( int i = 0; i < 15; i++) {
    vector < vector<double> > dummy1, dummy2;
    participantsFluxMinus.push_back(dummy1);
    participantsFluxPlus.push_back(dummy2);
  }
  for ( int i = 0; i < int(gluons.size()); i++ ) {
    if ( abs(gluons[i][5] + 2.5) < 0.5 ) {
      double r = sqrt(sqr(gluons[i][1] - CoGTouched.first) +
		      sqr(gluons[i][2] - CoGTouched.second));
      if ( r < 15 )
	participantsFluxMinus[int(r)].push_back(gluons[i]);
    }
    else if ( abs(gluons[i][5] - 2.5) < 0.5 ) {
      double r = sqrt(sqr(gluons[i][1] - CoGTouched.first) +
		      sqr(gluons[i][2] - CoGTouched.second));
      if ( r < 15 )
	participantsFluxPlus[int(r)].push_back(gluons[i]);
    }
  }

  for ( int i = 0; i < 15; i++ ) {
    for ( int j = 0; j < int(participantsFluxMinus[i].size()); j++ ) {
      double phi1 = atan2(participantsFluxMinus[i][j][3] - CoGTouched.second,
			  participantsFluxMinus[i][j][2] - CoGTouched.first);
      for ( int k = 0; k < int(participantsFluxPlus[i].size()); k++ ) {
	double phi2 = atan2(participantsFluxPlus[i][k][3] - CoGTouched.second,
			  participantsFluxPlus[i][k][2] - CoGTouched.first);
	double diff = min(abs(phi1 - phi2), 2*M_PI - abs(phi1 - phi2));
	histFluxTube->fill(double(i)+0.5, diff, W);
      }
    }
    for ( double j = M_PI/100.0; j < 2.0*M_PI; j += 2.0*M_PI/100.0 ) {
      histFluxTubeWeights->fill(double(i)+0.5, j, W*2.0*M_PI/100.0);
    }
  }
}

double FSAnalysis::reactionPlaneEccentricity(const vector<cPPtr> partons) {
  const int n = partons.size();
  assert(n > 0);
  double xsum  = 0.0;
  double ysum  = 0.0;
  double x2sum = 0.0;
  double y2sum = 0.0;
  for (vector<cPPtr>::const_iterator it = partons.begin(); it !=
	 partons.end(); ++ it) {
    cPPtr p = *it;
    const double x = p->vertex().x()/femtometer;
    const double y = p->vertex().y()/femtometer;
    xsum  += x;
    ysum  += y;
    x2sum += sqr(x);
    y2sum += sqr(y);
  }
  const double xave  = xsum  / (double)n;
  const double yave  = ysum  / (double)n;
  const double x2ave = x2sum / (double)n;
  const double y2ave = y2sum / (double)n;
  const double xx = x2ave - sqr(xave);
  const double yy = y2ave - sqr(yave);
  return (yy - xx) / (yy + xx);
}

/**
 * The definition of the participant eccentricity parameter.
 * PLB641 (2006) 260, Eq. 2.
 */
double FSAnalysis::participantEccentricity(const vector<cPPtr> partons) {
  const int n = partons.size();
  assert(n > 0);
  double xsum  = 0.0;
  double ysum  = 0.0;
  double x2sum = 0.0;
  double y2sum = 0.0;
  double xysum = 0.0;
  for (vector<cPPtr>::const_iterator it = partons.begin(); it !=
	 partons.end(); ++ it) {
    cPPtr p = *it;
    const double x = p->vertex().x()/femtometer;
    const double y = p->vertex().y()/femtometer;
    xsum  += x;
    ysum  += y;
    x2sum += sqr(x);
    y2sum += sqr(y);
    xysum += x * y;
  }
  const double xave  = xsum  / (double)n;
  const double yave  = ysum  / (double)n;
  const double x2ave = x2sum / (double)n;
  const double y2ave = y2sum / (double)n;
  const double xyave = xysum / (double)n;
  const double xx = x2ave - sqr(xave);
  const double yy = y2ave - sqr(yave);
  const double xy = xyave - xave * yave;
  return sqrt(sqr(yy - xx) + 4 * sqr(xy)) / (yy + xx);
}

/**
 * The definition of the average overlap area.
 * PRC81, 024901 (2010), Eq. 11.
 */
Area FSAnalysis::averageOverlapArea(const vector<cPPtr> partons) {
  const int n = partons.size();
  assert(n > 0);
  double xsum  = 0.0;
  double ysum  = 0.0;
  double x2sum = 0.0;
  double y2sum = 0.0;
  double xysum = 0.0;
  for (vector<cPPtr>::const_iterator it = partons.begin(); it !=
	 partons.end(); ++ it) {
    cPPtr p = *it;
    const double x = p->vertex().x()/femtometer;
    const double y = p->vertex().y()/femtometer;
    xsum  += x;
    ysum  += y;
    x2sum += sqr(x);
    y2sum += sqr(y);
    xysum += x * y;
  }
  const double xave  = xsum  / (double)n;
  const double yave  = ysum  / (double)n;
  const double x2ave = x2sum / (double)n;
  const double y2ave = y2sum / (double)n;
  const double xyave = xysum / (double)n;
  const double xx = x2ave - sqr(xave);
  const double yy = y2ave - sqr(yave);
  const double xy = xyave - xave * yave;
  if ( xx * yy - sqr(xy) < 0.0 ) return ZERO;
  return 4*M_PI * sqrt(xx * yy - sqr(xy))*sqr(femtometer);
}

Area FSAnalysis::averageOverlapArea(vector < vector<double> > gluons) {
  const int n = gluons.size();
  assert(n > 0);
  double xsum  = 0.0;
  double ysum  = 0.0;
  double x2sum = 0.0;
  double y2sum = 0.0;
  double xysum = 0.0;
  for (int i = 0; i < n; i++) {
    const double x = gluons[i][1];
    const double y = gluons[i][2];
    xsum  += x;
    ysum  += y;
    x2sum += sqr(x);
    y2sum += sqr(y);
    xysum += x * y;
  }
  const double xave  = xsum  / (double)n;
  const double yave  = ysum  / (double)n;
  const double x2ave = x2sum / (double)n;
  const double y2ave = y2sum / (double)n;
  const double xyave = xysum / (double)n;
  const double xx = x2ave - sqr(xave);
  const double yy = y2ave - sqr(yave);
  const double xy = xyave - xave * yave;
  if ( xx * yy - sqr(xy) < 0.0 ) return ZERO;
  return 4*M_PI * sqrt(xx * yy - sqr(xy))*sqr(femtometer);
}

void FSAnalysis::v2GluonAnalyze(tEventPtr event) {
  SubProPtr sub = event->primarySubProcess();
  double W = event->weight();

  tPVector particles = event->getFinalState();
  vector<cPPtr> centrals;
  int multF = 0;
  for ( int i = 0; i < int(particles.size()); i++ ) {
    if ( particles[i]->data().charge() != ZERO && abs(particles[i]->eta()) < 0.9 ) {
      multF++;
      centrals.push_back(particles[i]);
    }
  }

  double pairs = 0.0;
  double cos2 = 0.0;
  for ( int i = 0; i < int(centrals.size()); i++ ) {
    for ( int j = 0; j < int(centrals.size()); j++ ) {
      if ( i == j ) continue;
      pairs++;
      cos2 += cos(2.0*(centrals[i]->momentum().phi() - centrals[j]->momentum().phi()));
    }
  }
  if ( pairs > 0.0 ) {
    cos2 /= pairs;
    v2final->fill(multF, cos2*W);
    v2finalWeights->fill(multF, W);
    dNchFinal->fill(multF, W);
  }

  int pairsgap = 0;
  double eliptic12 = 0.0;
  for ( int i = 0; i < int(particles.size()); i++ ) {
    if ( particles[i]->data().charge() != ZERO
	 && particles[i]->eta() < -0.5 && particles[i]->eta() > -2.0 ) {
      for ( int j = 0; j < i; j++ ) {
	if ( particles[j]->data().charge() != ZERO
	     && particles[j]->eta() > 0.5 && particles[j]->eta() < 2.0  ) {
	  pairsgap++;
	  Energy pix = particles[i]->momentum().x();
	  Energy piy = particles[i]->momentum().y();
	  Energy pjx = particles[j]->momentum().x();
	  Energy pjy = particles[j]->momentum().y();
	  Energy pti = sqrt(sqr(pix) + sqr(piy));
	  Energy ptj = sqrt(sqr(pjx) + sqr(pjy));
	  eliptic12 += 1.0-2.0*sqr((pix*pjy - pjx*piy)/(pti*ptj));
	}
      }
    }
  }
  if ( pairsgap > 0 ) {
    eliptic12 /= double(pairsgap);
    v2gap->fill(multF, event->weight()*eliptic12);
  }

  //find gluons in [-5,5]. Should cut away the spectators.
  vector<cPPtr> participants;
  int multG = 0.0;
  for ( int i = 0; i < int(sub->outgoing().size()); i++ ) {
    if ( abs(sub->outgoing()[i]->eta()) < 5.0 ) {
      multG++;
      participants.push_back(sub->outgoing()[i]);
    }
  }

  nParticipants->fill(multG, W);
  averageNParticipants += multG*event->weight();

  if ( multG > 2 ) {
    //calculate the eccentricity and area of the participants
    double RPE = reactionPlaneEccentricity(participants);
    double PE = participantEccentricity(participants);
    Area AOA = averageOverlapArea(participants);
    
    RPEcc->fill(RPE, W);
    PartEcc->fill(PE, W);
    OverlapArea->fill(AOA/sqr(femtometer), W);

    EccArea->fill(AOA/sqr(femtometer), PE*W);
    AreaPart->fill(multG, AOA/sqr(femtometer)*W);
    PE2->fill(multG,sqr(PE)*W);
    PE4->fill(multG,sqr(sqr(PE))*W);

    //calculate the square of eq 5 in the v2 paper event by event
    double averageMultG = 7.6;
    double averagedN = 1.5*6.2;
    double dN = multG*averagedN/averageMultG;
    double dNch = dN/1.5;
    double eq5Sqr = sqr(PE*0.2/(1.0 + (5.8/sqr(femtometer))*AOA/dN));
    v2flow->fill(dNch, eq5Sqr*W);
    dNchFlow->fill(dNch, W);
    density->fill(dNch, dN/(AOA/sqr(femtometer))*W);
    if ( dNch > 10. && dN/1.5 < 20. ) density1020->fill(dN/(AOA/sqr(femtometer)), W);
    if ( dNch > 40. && dN/1.5 < 60. ) density4060->fill(dN/(AOA/sqr(femtometer)), W);
    if ( dNch > 60. ) density60plus->fill(dN/(AOA/sqr(femtometer)), W);

    //and the same for v_2(4)^4
    double flowquad = sqr(sqr(PE*0.2/(1.0 + (5.8/sqr(femtometer))*AOA/dN)));
    v4flow->fill(dNch, flowquad*W);

    if ( dN/AOA > 5.8/sqr(femtometer) ) {
      v2mix->fill(dNch*1.8, eq5Sqr*W);
      dNchMix->fill(dNch*1.8, W);
      denseRatio->fill(dNch*1.8, W);
    }
    else {
      v2mix->fill(multF, cos2*W);
      dNchMix->fill(multF, W);
    }
    denseRatioWeights->fill(dNch*1.8, W);
  }
  else {
    v2mix->fill(multF, cos2*W);
    dNchMix->fill(multF, W);
  }

  tPVector chargedMid;
  for ( int i = 0; i < int(particles.size()); i++ )
    if ( particles[i]->data().charge() != ZERO && abs(particles[i]->eta()) < 0.9 )
      chargedMid.push_back(particles[i]);

  double quad = 0.0;
  int quads = 0;
  for ( int i = 0; i < int(chargedMid.size()); i++ ) {
    for ( int j = 0; j < int(chargedMid.size()); j++ ) {
      if ( j == i ) continue;
      for ( int k = 0; k < int(chargedMid.size()); k++ ) {
	if ( k == i || k == j ) continue;
	for ( int l = 0; l < int(chargedMid.size()); l++ ) {
	  if ( l == i || l == j || l == k ) continue;
	  quads++;
	  quad += cos(2.0*(chargedMid[i]->momentum().phi() + chargedMid[j]->momentum().phi() -
			   chargedMid[k]->momentum().phi() - chargedMid[l]->momentum().phi()));
	}
      }
    }
  }
  quad /= double(quads);
  if ( quads > 0 ) {
    v2quad->fill(multF, event->weight()*quad);
    v2quadWeights->fill(multF, event->weight());
  }

}

void FSAnalysis::analyze(tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);

  // exportToSPheRIO(event);
  // getchar();

  PPair incoming = event->incoming();

  // totalWeight += event->weight();
  weights->fill(log(event->weight()), 1.0);

  // First go through the gluons coming directly from DIPSY
  SubProPtr sub = event->primarySubProcess();
  int nglue = 0;
  Energy maxPT = ZERO;
  Energy sumPT = ZERO;
  for ( int i = 0, N = sub->outgoing().size(); i < N; ++i ) {
    PPtr p = sub->outgoing()[i];
    if ( p->id() == ParticleID::g ) ++nglue;
    Energy pt = sqrt(sqr(p->momentum().x()) + sqr(p->momentum().y()));
    DIPSYPT->fill(p->rapidity(), pt/GeV*event->weight());
    if ( pt > maxPT ) maxPT = pt;
    if ( abs(p->eta()) < 1.0 ) sumPT += pt;
  }
  DIPSYSumPT->fill(sumPT/GeV, event->weight());
  ngbefore->fill(nglue, event->weight());
  maxPTBefore->fill(maxPT/GeV, event->weight());
  int nGlueDIPSY = nglue;

  // Now find the last step before string fragmentation
  int colstep = 0;
  for ( int is = 0, Ns = event->primaryCollision()->steps().size(); is < Ns; ++is ) {
    StepPtr step = event->primaryCollision()->steps()[is];
    for ( ParticleSet::iterator it = step->particles().begin();
	  it != step->particles().end(); ++it ) {
      PPtr p = *it;
      if ( p->hasColour() ) colstep = is;
    }
  }

  // Go through the gluons after the final state shower, before hadronisation.
  StepPtr step = event->primaryCollision()->steps()[colstep];
  nglue = 0;
  for ( ParticleSet::iterator it = step->particles().begin();
	it != step->particles().end(); ++it ) {
    PPtr p = *it;
    Energy pt = sqrt(sqr(p->momentum().x()) + sqr(p->momentum().y()));
    if ( p->id() == ParticleID::g ) ++nglue;
    secondLastPT->fill(p->rapidity(), pt/GeV*event->weight());

    if ( p->id() == ParticleID::g ) {
      etaDiffAriadne += event->weight()*abs(step->colourNeighbour(p)->eta() - p->eta());
      NGlueAriadne += event->weight();
    }
  }
  ngafter->fill(nglue, event->weight());

  //now look at final state
  bool isINEL = false;
  int Ncharge = 0;
  int nGamma = 0;
  tPVector particles = event->getFinalState();
  Energy sumET = ZERO;
  for ( int i = 0; i < int(particles.size()); i++ ) {
    PPtr p = particles[i];
    Energy pt = sqrt(sqr(p->momentum().x()) + sqr(p->momentum().y()));
    Energy Et = sqrt(sqr(p->momentum().x()) + sqr(p->momentum().y()) +sqr(p->mass()));
    finalPT->fill(p->rapidity(), pt/GeV*event->weight());
    finalPTeta->fill(p->eta(), pt/GeV*event->weight());
    finalETeta->fill(p->eta(), Et/GeV*event->weight());
    if ( p->data().charge() != ZERO ) {
      chargedPT->fill(p->rapidity(), pt/GeV*event->weight());
      finalChargedPTeta->fill(p->eta(), pt/GeV*event->weight());
      finalNcheta->fill(p->eta(), event->weight());
      chargedN->fill(p->rapidity(), event->weight());
      if ( abs(particles[i]->eta()) < 1 ) {
	sumET += Et;
	isINEL = true; 
	chargedNALICE->fill(p->eta(), event->weight());
	Ncharge++;
      }
      if ( abs(particles[i]->eta()) < 1.0 )
	EtDensPHENIX += event->weight()*Et;
      if ( p->eta() > 2.0 ) nGamma++;
    }
  }

  //for the SD analysis. first find largest rapidity gap
  double gapY = 0.0;
  double maxGap = 0.0;
  PSet ordered;
  for ( int i = 0; i < int(particles.size()); i++ )
    ordered.insert(particles[i]);
  for ( PSet::iterator it = ordered.begin(); it != ordered.end(); it++ ) {
    if ( it == ordered.begin() ) continue;
    PPtr p = *it;
    it--;
    PPtr prev = *it;
    it++;
    double eta1 = p->eta();
    if ( abs(p->eta()) > 1000.0 ) {
      if ( p->momentum().z() > ZERO ) eta1 = 10.0;
      else eta1 = -10.0;
    }
    double eta2 = prev->eta();
    if ( abs(prev->eta()) > 1000.0 ) {
      if ( prev->momentum().z() > ZERO ) eta2 = 10.0;
      else eta2 = -10.0;
    }
    if ( eta1 - eta2 > maxGap ) {
      gapY = (eta1 + eta2)/2.0;
      maxGap = eta1 - eta2;
    }
  }
  //UA4 trigger: one charged particle in 2.5 - 5.6 eta.
  bool AU4SD = false;
  for ( int i = 0; i < int(particles.size()); i++ ) {
    PPtr p = particles[i];
    if ( p->data().charge() != ZERO && p->eta() > 2.5 && p->eta() < 5.6 )
      AU4SD = true;
  }
  if ( AU4SD ) {
    AU4SDweights->fill(log(event->weight()), 1.0);
    PSet X;
    for ( PSet::iterator it = ordered.begin(); it != ordered.end(); it++ ) {
      PPtr p = *it;
      double eta = p->eta();
      if ( abs(p->eta()) > 1000.0 ) {
	if ( p->momentum().z() > ZERO ) eta = 10.0;
	else eta = -10.0;
      }
      // cout << ", " << eta;
      if ( eta > gapY )  X.insert(p);
    }
    // cout << endl;
    int NchDiff = 0;
    Lorentz5Momentum diffMom;
    for ( PSet::iterator it = X.begin(); it != X.end(); it++ ) {
      PPtr p = *it;
      diffMom += p->momentum();
      if ( p->data().charge() != ZERO ) NchDiff++;
    }
    diffMom.rescaleMass();
    // cout << "weight is " << event->weight()
    //  	 << ", M_X is " << diffMom.mass()/GeV << ", <n> is " << NchDiff << endl;
    // cout << "-----------------------------------------------------" << endl;
    double MX = diffMom.mass()/GeV;

    for ( PSet::iterator it = X.begin(); it != X.end(); it++ ) {
      PPtr p = *it;
      if ( p->eta() != ZERO && p->data().charge() != ZERO ) {
	if ( MX < 50 ) UA4dndetaMX20->fill(p->eta(), event->weight());
	else if ( MX < 70 ) UA4dndetaMX60->fill(p->eta(), event->weight());
	else if ( MX < 90 ) UA4dndetaMX80->fill(p->eta(), event->weight());
	else if ( MX < 110 ) UA4dndetaMX100->fill(p->eta(), event->weight());
	else if ( MX < 120 ) UA4dndetaMX115->fill(p->eta(), event->weight());
	else if ( MX < 160 ) UA4dndetaMX140->fill(p->eta(), event->weight());

      }
    }
    if ( MX < 50 ) {
      for ( double eta = -10+0.5*0.1; eta < 10; eta += 0.1 ) {
	UA4dndetaMX20Weights->fill(eta, event->weight()*0.1);
      }
    } else if ( MX < 70 ) {
      for ( double eta = -10+0.5*0.1; eta < 10; eta += 0.1 ) {
	UA4dndetaMX60Weights->fill(eta, event->weight()*0.1);
      }
    } else if ( MX < 90 ) {
      for ( double eta = -10+0.5*0.1; eta < 10; eta += 0.1 ) {
	UA4dndetaMX80Weights->fill(eta, event->weight()*0.1);
      }
    } else if ( MX < 110 ) {
      for ( double eta = -10+0.5*0.1; eta < 10; eta += 0.1 ) {
	UA4dndetaMX100Weights->fill(eta, event->weight()*0.1);
      }
    } else if ( MX < 120 ) {
      for ( double eta = -10+0.5*0.1; eta < 10; eta += 0.1 ) {
	UA4dndetaMX115Weights->fill(eta, event->weight()*0.1);
      }
    } else if ( MX < 160 ) {
      for ( double eta = -10+0.5*0.1; eta < 10; eta += 0.1 ) {
	UA4dndetaMX140Weights->fill(eta, event->weight()*0.1);
      }
    }

    double lnMX2 = -1000.0;
    if ( MX > 0.001 ) lnMX2 = 2.0*log(MX);
    SDavN->fill(lnMX2, NchDiff*event->weight());
    SDavNWeights->fill(lnMX2, event->weight());
    SDavN2->fill(lnMX2, sqr(NchDiff)*event->weight());
    SDavN2Weights->fill(lnMX2, event->weight());

  }
  if ( maxGap > 3 ) { //require rapidity gap of 3 units to do SD analysis

    //DipoleEventHandlerPtr eh = Current<DipoleEventHandler>().ptr();
    //EventFiller & def = eh->eventFiller();
    //cout << "ncascades: " << eh->eventFiller().nSubcascades() << endl;

    PSet X;
    for ( PSet::iterator it = ordered.begin(); it != ordered.end(); it++ ) {
      PPtr p = *it;
      double eta = p->eta();
      if ( abs(p->eta()) > 1000.0 ) {
	if ( p->momentum().z() > ZERO ) eta = 10.0;
	else eta = -10.0;
      }
      // cout << ", " << eta;
      if ( eta > gapY )  X.insert(p);
    }
    // cout << endl;
    int NchDiff = 0;
    Lorentz5Momentum diffMom;
    for ( PSet::iterator it = X.begin(); it != X.end(); it++ ) {
      PPtr p = *it;
      diffMom += p->momentum();
      if ( p->data().charge() != ZERO ) NchDiff++;
    }
    diffMom.rescaleMass();
    // cout << "weight is " << event->weight()
    // 	 << ", M_X is " << diffMom.mass()/GeV << ", <n> is " << NchDiff << endl;
    double MX = diffMom.mass()/GeV;
    // if ( MX > 120 ) getchar();
    double lnMX2 = -1000.0;
    if ( MX > 0.001 )
      lnMX2 = 2.0*log(MX);
    lnMX2dist->fill(lnMX2, event->weight());
    lnMX2distFrac->fill(lnMX2, event->weight());
    if ( nGlueDIPSY == 0 ) lnMX2dist0->fill(lnMX2, event->weight());
    if ( nGlueDIPSY == 1 ) lnMX2dist1->fill(lnMX2, event->weight());
    if ( nGlueDIPSY > 1 ) lnMX2dist2->fill(lnMX2, event->weight());
    for ( double y = -1.0+0.5*0.1; y < 19; y += 0.1 ) {
      lnMX2distWeights->fill(y, event->weight());
      lnMX2dist0Weights->fill(y, event->weight());
      lnMX2dist1Weights->fill(y, event->weight());
      lnMX2dist2Weights->fill(y, event->weight());
    }
    // cout << "gap: " << gapY-maxGap/2.0 << " to " << gapY+maxGap/2.0 << ", MX: " << MX
    // 	 << ", ln(MX^2): " << lnMX2 << ", N: " << NchDiff << endl;
    //the rapidity of the X system
    MXmaxY->fill(gapY + maxGap/2.0, event->weight()*MX);
    MXmaxYWeights->fill(gapY + maxGap/2.0, event->weight());
    MX2maxY->fill(gapY + maxGap/2.0, event->weight()*sqr(MX));
    MX2maxYWeights->fill(gapY + maxGap/2.0, event->weight());
    double CMS = 0.0;
    if ( MX > 0.001 )
      CMS = 0.5*log(diffMom.plus()/diffMom.minus());
    double x = sqrt(diffMom.plus()/diffMom.minus());
    SDweights->fill(log(event->weight()), 1.0);
    // cout << "CMS y: " << CMS << endl;
    for ( PSet::iterator it = X.begin(); it != X.end(); it++ ) {
      PPtr p = *it;
      SDNchdeta->fill(p->eta(), event->weight());
      double y = 0.5*log(p->momentum().plus()/p->momentum().minus()) - CMS;
      //boost to rest system of X
      Energy plusBoost = p->momentum().plus()/x;
      Energy minusBoost = p->momentum().minus()*x;
      Energy eBoost = (plusBoost + minusBoost)/2.0;
      Energy pzBoost = (plusBoost - minusBoost)/2.0;
      Energy pT = sqrt(sqr(p->momentum().x()) + sqr(p->momentum().y()));
      double etaStar = -log(tan(atan2(pT,pzBoost)/2.0));
      // cout << "this particle, rel. y: " << y << ", rel. eta: " << etaStar << endl;
      // cout << "pT: " << pT/GeV << ", pzBoost: " << pzBoost/GeV
      // 	   << ", mass: " << p->mass()/GeV << endl;
      //dndy
      if ( MX > 3.0 && MX < 8.0 && p->data().charge() != ZERO )
	SDNchdyMX38->fill(y, event->weight());
      if ( MX > 8.0 && MX < 15.0 && p->data().charge() != ZERO )
	SDNchdyMX815->fill(y, event->weight());
      if ( MX > 15.0 && MX < 30.0 && p->data().charge() != ZERO )
	SDNchdyMX1530->fill(y, event->weight());
      //dE/deta
      if ( MX > 3.0 && MX < 8.0 )
	SDdedetaMX38->fill(etaStar, eBoost/GeV*event->weight());
      if ( MX > 8.0 && MX < 18.0 )
	SDdedetaMX818->fill(etaStar, eBoost/GeV*event->weight());
      if ( MX > 18.0 && MX < 30.0 )
	SDdedetaMX1830->fill(etaStar, eBoost/GeV*event->weight());
    }
    if ( MX > 3.0 && MX < 8.0 ) {
      for ( double y = -10+0.5*0.1; y < 10; y += 0.1 ) {
	SDNchdyMX38Weights->fill(y, event->weight());
      }
      for ( double eta = -10+0.5*0.1; eta < 10; eta += 0.1 ) {
	SDdedetaMX38Weights->fill(eta, event->weight());
      }
    }
    if ( MX > 8.0 && MX < 15.0 ) {
      for ( double y = -10+0.5*0.1; y < 10; y += 0.1 ) {
	SDNchdyMX815Weights->fill(y, event->weight());
      }
    }
    if ( MX > 8.0 && MX < 18.0 ) {
      for ( double eta = -10+0.5*0.1; eta < 10; eta += 0.1 ) {
	SDdedetaMX818Weights->fill(eta, event->weight());
      }
    }
    if ( MX > 15.0 && MX < 30.0 ) {
      for ( double y = -10+0.5*0.1; y < 10; y += 0.1 ) {
	SDNchdyMX1530Weights->fill(y, event->weight());
      }
    }
    if ( MX > 18.0 && MX < 30.0 ) {
      for ( double eta = -10+0.5*0.1; eta < 10; eta += 0.1 ) {
	SDdedetaMX1830Weights->fill(eta, event->weight());
      }
    }

    // SDavN->fill(lnMX2, NchDiff*event->weight());
    // SDavNWeights->fill(lnMX2, event->weight());
    // SDavN2->fill(lnMX2, sqr(NchDiff)*event->weight());
    // SDavN2Weights->fill(lnMX2, event->weight());
  }

  finalSumET->fill(sumET/GeV, event->weight());
  if ( isINEL ) {
    INEL += event->weight();
    chargeMultALICE->fill(Ncharge, event->weight());

  }

  //track the 2D energy density in a central rapidity slice
  sub = event->primarySubProcess();
  for ( int i = 0, N = sub->outgoing().size(); i < N; ++i ) {
    PPtr p = sub->outgoing()[i];
    if ( p->rapidity() < -1.0 || p->rapidity() > 1.0 ) continue;
    Energy E = sqrt(sqr(p->momentum().x()) + sqr(p->momentum().y()) + sqr(p->momentum().z()) +sqr(p->mass()));
    double x1 = p->vertex().x()/femtometer;
    double x2 = p->vertex().y()/femtometer;
    test2D->fill(x1, x2, E/GeV*event->weight());
    test2DWeights->fill(x1, x2, event->weight());
  }

  //   v2GluonAnalyze(event);

   // eccentricities(event);

   // saveGluonsToFile(event);
}

void FSAnalysis::saveGluonsToFile(tEventPtr event) {
  string filename = generator()->filename();

    //set up outfile
  ostringstream os;
  os << "events/" << filename << "/event" << event->number() << ".dat";
  ofstream outfile (os.str().c_str());

  //general information
  outfile << "Generated with DIPSY 1.1" << endl;
  outfile << "These are the real gluons only. No FSR or hadronisation has been done." << endl;

  //print impact parameter.
  PPair incoming = event->incoming();
  LorentzPoint p1 = incoming.first->vertex();
  LorentzPoint p2 = incoming.second->vertex();
  double bx = (p1-p2).x()/femtometer;
  double by = (p1-p2).y()/femtometer;
  outfile << "impact_parameter(fm):" << " " << bx << " " << by << endl;

  SubProPtr sub = event->primarySubProcess();
  outfile << "number_of_particles:" << " " <<sub->outgoing().size() << endl << endl;
  outfile << "number of spectating nucleons: " << numberOfSpectators(event) << endl;
  
  //print gluons
  outfile << "particle_number" << " "
	  << "transverse_position_x(fm)" << " "
	  << "transverse_position_y(fm)" << " "
	  << "transverse_momentum_x(GeV)" << " "
	  << "transverse_momentum_y(GeV)" << " "
	  << "rapidity" << " "
	  << "colour_neighbour_number" << " "
	  << "anticolour_neighbour_number" << endl << endl;
  for ( int i = 0, N = sub->outgoing().size(); i < N; ++i ) {
    PPtr p = sub->outgoing()[i];
    outfile << p->number() << " "
	    << p->vertex().x()/femtometer << " " 
	    << p->vertex().y()/femtometer << " " 
	    << p->momentum().x()/GeV << " " 
	    << p->momentum().y()/GeV << " " 
	    << p->eta() << " ";
    if ( event->primaryCollision()->steps()[0]->colourNeighbour(p) )
      outfile << event->primaryCollision()->steps()[0]->colourNeighbour(p)->number() << " ";
    else outfile << -1 << " ";
    if ( event->primaryCollision()->steps()[0]->antiColourNeighbour(p) )
      outfile << event->primaryCollision()->steps()[0]->antiColourNeighbour(p)->number() << " ";
    else outfile << -1 << " ";
    outfile  << endl;
  }

  outfile.close();
  cout << "printed gluons to file " << "events/" << filename << "/event" << event->number()
       << ".dat" << endl;
}

void FSAnalysis::exportToSPheRIO(tEventPtr event) {
  string filename = generator()->filename();

    //set up outfile
  ostringstream os;
  os << "SPheRIO/" << filename << "/event" << event->number() << ".dat";
  ofstream outfile (os.str().c_str());

  PPair incoming = event->incoming();

  //Energy
  outfile << (incoming.first->momentum().t() + incoming.first->momentum().t())/GeV << endl;

  //nuclei info CALL WF SOMEHOW!! FIX!!
  int targetN = 79;
  int targetA = 197;
  int projectileN = 79;
  int projectileA = 197;
  outfile << targetN << " " << targetA << endl;
  outfile << projectileN << " " << projectileA << endl;

  //print impact parameter.
  LorentzPoint p1 = incoming.first->vertex();
  LorentzPoint p2 = incoming.second->vertex();
  double bx = (p1-p2).x()/femtometer;
  double by = (p1-p2).y()/femtometer;
  outfile << bx << " " << by << endl;

  //set up output matrix
  int xBins = 25;
  double xMin = -10.0;
  double xMax = 10.0;
  double xBinSize = (xMax-xMin)/double(xBins);
  int yBins = 25;
  double yMin = -10.0;
  double yMax = 10.0;
  double yBinSize = (yMax-yMin)/double(yBins);
  int zBins = 25;
  double zMin = -10.0;
  double zMax = 10.0;
  double zBinSize = (zMax-zMin)/double(zBins);
  vector<double> bin(23, 0.0);
  vector<vector<double> > xline(xBins, bin);
  vector<vector<vector<double> > > xyPlane(yBins, xline);
  vector<vector<vector<vector<double> > > > theMatrix(zBins, xyPlane);

  //the transverse length the particles are streamed
  Length freeStream = 1.0*femtometer;

  //fill the matrix
  tPVector particles = event->getFinalState();
  for ( int i = 0; i < int(particles.size()); i++ ) {
    PPtr p = particles[i];
    Lorentz5Momentum q = p->momentum();
    //free stream
    Length deltaX = freeStream*p->momentum().x()/
      sqrt(sqr(p->momentum().x())+sqr(p->momentum().y()));
    Length deltaY = freeStream*p->momentum().y()/
      sqrt(sqr(p->momentum().x())+sqr(p->momentum().y()));
    double x = (p->vertex().x()+deltaX)/femtometer;
    double y = (p->vertex().y()+deltaY)/femtometer;
    double z = p->eta();
    int xBin = int((x-xMin)/xBinSize);
    if ( xBin < 0 || xBin > xBins-1 ) continue;
    int yBin = int((y-yMin)/yBinSize);
    if ( yBin < 0 || yBin > yBins-1 ) continue;
    int zBin = int((z-zMin)/zBinSize);
    if ( zBin < 0 || zBin > zBins-1 ) continue;
    //T00, T01, ..., T33
    theMatrix[xBin][yBin][zBin][3] += q.t()/GeV; //T00
    theMatrix[xBin][yBin][zBin][4] += q.x()/GeV;
    theMatrix[xBin][yBin][zBin][5] += q.y()/GeV;
    theMatrix[xBin][yBin][zBin][6] += q.z()/GeV;
    theMatrix[xBin][yBin][zBin][7] += q.x()/GeV; //T10
    theMatrix[xBin][yBin][zBin][8] += q.x()*q.x()/q.t()/GeV;
    theMatrix[xBin][yBin][zBin][9] += q.x()*q.y()/q.t()/GeV;
    theMatrix[xBin][yBin][zBin][10] += q.x()*q.z()/q.t()/GeV;
    theMatrix[xBin][yBin][zBin][11] += q.y()/GeV; //T20
    theMatrix[xBin][yBin][zBin][12] += q.y()*q.x()/q.t()/GeV;
    theMatrix[xBin][yBin][zBin][13] += q.y()*q.y()/q.t()/GeV;
    theMatrix[xBin][yBin][zBin][14] += q.y()*q.z()/q.t()/GeV;
    theMatrix[xBin][yBin][zBin][15] += q.z()/GeV;   //T30
    theMatrix[xBin][yBin][zBin][16] += q.z()*q.x()/q.t()/GeV;
    theMatrix[xBin][yBin][zBin][17] += q.z()*q.y()/q.t()/GeV;
    theMatrix[xBin][yBin][zBin][18] += q.z()*q.z()/q.t()/GeV;
    //energy density
    theMatrix[xBin][yBin][zBin][19] += q.t()/(xBinSize*yBinSize*zBinSize)/GeV;
  }
  
  outfile << "1 1 1" << endl;
  outfile << xBins << " " << yBins << " " << zBins << endl;
  outfile << xMin << " " << xMax << " "
	  << yMin << " " << yMax << " "
	  << zMin << " " << zMax << " " << endl;
  outfile << endl;
  for ( int i = 0; i < int(xBins); i++ ) {
    for ( int j = 0; j < int(yBins); j++ ) {
      for ( int k = 0; k < int(zBins); k++ ) {
	theMatrix[i][j][k][0] = xMin + xBinSize*(double(i)+0.5);
	theMatrix[i][j][k][1] = yMin + yBinSize*(double(j)+0.5);
	theMatrix[i][j][k][2] = zMin + zBinSize*(double(k)+0.5);
	if ( theMatrix[i][j][k][3] != 0.0 ) {
	  double E = theMatrix[i][j][k][3];
	  double px = theMatrix[i][j][k][4];
	  double py = theMatrix[i][j][k][5];
	  double pz = theMatrix[i][j][k][6];
	  theMatrix[i][j][k][20] = px/E;
	  theMatrix[i][j][k][21] = py/E;
	  theMatrix[i][j][k][22] = pz/E;
	}
	for (int l = 0; l < int(theMatrix[i][j][k].size()); l++) {
	  outfile << theMatrix[i][j][k][l] << " ";
	}
	outfile << endl;
      }
    }
  }

  outfile.close();
  cout << "printed gluons to file " << "events/" << filename << "/event"
       << event->number() << ".dat" << endl;
}

int FSAnalysis::numberOfSpectators(tEventPtr event) {
  SubProPtr sub = event->primarySubProcess();
  vector<cPPtr> unTouched;
  for ( int i = 0; i < int(sub->outgoing().size()); i++ ) {
    if ( sub->outgoing()[i]->momentum().perp() < 0.47*GeV &&
	 sub->outgoing()[i]->momentum().perp() > 0.46*GeV ) {
      unTouched.push_back(sub->outgoing()[i]);
    }
  }
  int nSpectators = 0;
  for ( int i = 0, N = unTouched.size(); i < N; ++i ) {
    //check if third colour neigbour is itself, then spectator.
    cPPtr p = unTouched[i];
    cPPtr p1 = event->primaryCollision()->steps()[0]->colourNeighbour(p);
    cPPtr p2 = event->primaryCollision()->steps()[0]->colourNeighbour(p1);
    cPPtr p3 = event->primaryCollision()->steps()[0]->colourNeighbour(p2);
    if ( p3 == p ) nSpectators++;
  }
  return nSpectators/3;
}

LorentzRotation FSAnalysis::transform(tEventPtr event) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void FSAnalysis::analyze(const tPVector & particles) {
  AnalysisHandler::analyze(particles);
  // Calls analyze() for each particle.
}

void FSAnalysis::analyze(tPPtr) {}

IBPtr FSAnalysis::clone() const {
  return new_ptr(*this);
}

IBPtr FSAnalysis::fullclone() const {
  return new_ptr(*this);
}


void FSAnalysis::doinitrun() {
  AnalysisHandler::doinitrun();
  histogramFactory().initrun();
  histogramFactory().registerClient(this);
  epsRP = histogramFactory().createHistogram1D
    ("epsRP",200,-0.0,1.0);
  epsPart = histogramFactory().createHistogram1D
    ("epsPart",200,-0.0,1.0);
  epsPartOver20 = histogramFactory().createHistogram1D
    ("epsPartOver20",200,-0.0,1.0);
  epsPartOver40 = histogramFactory().createHistogram1D
    ("epsPartOver40",200,-0.0,1.0);
  mult12 = histogramFactory().createHistogram1D
    ("mult12",200,-0.5,199.5);
  weights = histogramFactory().createHistogram1D
    ("weights",200,-30.0,5.0);
  SDweights = histogramFactory().createHistogram1D
    ("SDweights",200,-20.0,20.0);
  AU4SDweights = histogramFactory().createHistogram1D
    ("AU4SDweights",200,-20.0,20.0);
  SDavN = histogramFactory().createHistogram1D
    ("SDavN",30, 0.0, 10.0);
  SDavNWeights = histogramFactory().createHistogram1D
    ("SDavNWeights",30, 0.0, 10.0);
  SDavN2 = histogramFactory().createHistogram1D
    ("SDavN2",30, 0.0, 10.0);
  SDavN2Weights = histogramFactory().createHistogram1D
    ("SDavN2Weights",30, 0.0, 10.0);
  ngbefore = histogramFactory().createHistogram1D
    ("ngbefore",100,-0.5,99.5);
  ngafter = histogramFactory().createHistogram1D
    ("ngafter",100,-0.5,99.5);
  DIPSYPT = histogramFactory().createHistogram1D
    ("DIPSYPT",200,-10.0,10.0);
  finalPT = histogramFactory().createHistogram1D
    ("finalPT",200,-10.0,10.0);
  finalPTeta = histogramFactory().createHistogram1D
    ("finalPTeta",200,-10.0,10.0);
  finalChargedPTeta = histogramFactory().createHistogram1D
    ("finalChargedPTeta",200,-10.0,10.0);
  finalETeta = histogramFactory().createHistogram1D
    ("finalETeta",200,-10.0,10.0);
  finalNcheta = histogramFactory().createHistogram1D
    ("finalNcheta",200,-10.0,10.0);
  lnMX2dist = histogramFactory().createHistogram1D
    ("lnMX2dist",200,-1.0,19.0);
  lnMX2distFrac = histogramFactory().createHistogram1D
    ("lnMX2distFrac",200,-1.0,19.0);
  lnMX2distWeights = histogramFactory().createHistogram1D
    ("lnMX2distWeights",200,-1.0,19.0);
  lnMX2dist0 = histogramFactory().createHistogram1D
    ("lnMX2dist0",200,-1.0,19.0);
  lnMX2dist0Weights = histogramFactory().createHistogram1D
    ("lnMX2dist0Weights",200,-1.0,19.0);
  lnMX2dist1 = histogramFactory().createHistogram1D
    ("lnMX2dist1",200,-1.0,19.0);
  lnMX2dist1Weights = histogramFactory().createHistogram1D
    ("lnMX2dist1Weights",200,-1.0,19.0);
  lnMX2dist2 = histogramFactory().createHistogram1D
    ("lnMX2dist2",200,-1.0,19.0);
  lnMX2dist2Weights = histogramFactory().createHistogram1D
    ("lnMX2dist2Weights",200,-1.0,19.0);
  SDNchdeta = histogramFactory().createHistogram1D
    ("SDNchdeta",200,-10.0,10.0);
  MXmaxY = histogramFactory().createHistogram1D
    ("MXmaxY",200,-10.0,10.0);
  MXmaxYWeights = histogramFactory().createHistogram1D
    ("MXmaxYWeights",200,-10.0,10.0);
  MX2maxY = histogramFactory().createHistogram1D
    ("MX2maxY",200,-10.0,10.0);
  MX2maxYWeights = histogramFactory().createHistogram1D
    ("MX2maxYWeights",200,-10.0,10.0);

  //UA4
  UA4dndetaMX20 = histogramFactory().createHistogram1D
    ("UA4dndetaMX20",200,-10.0,10.0);
  UA4dndetaMX20Weights = histogramFactory().createHistogram1D
    ("UA4dndetaMX20Weights",200,-10.0,10.0);
  UA4dndetaMX60 = histogramFactory().createHistogram1D
    ("UA4dndetaMX60",200,-10.0,10.0);
  UA4dndetaMX60Weights = histogramFactory().createHistogram1D
    ("UA4dndetaMX60Weights",200,-10.0,10.0);
  UA4dndetaMX80 = histogramFactory().createHistogram1D
    ("UA4dndetaMX80",200,-10.0,10.0);
  UA4dndetaMX80Weights = histogramFactory().createHistogram1D
    ("UA4dndetaMX80Weights",200,-10.0,10.0);
  UA4dndetaMX100 = histogramFactory().createHistogram1D
    ("UA4dndetaMX100",200,-10.0,10.0);
  UA4dndetaMX100Weights = histogramFactory().createHistogram1D
    ("UA4dndetaMX100Weights",200,-10.0,10.0);
  UA4dndetaMX115 = histogramFactory().createHistogram1D
    ("UA4dndetaMX115",200,-10.0,10.0);
  UA4dndetaMX115Weights = histogramFactory().createHistogram1D
    ("UA4dndetaMX115Weights",200,-10.0,10.0);
  UA4dndetaMX140 = histogramFactory().createHistogram1D
    ("UA4dndetaMX140",200,-10.0,10.0);
  UA4dndetaMX140Weights = histogramFactory().createHistogram1D
    ("UA4dndetaMX140Weights",200,-10.0,10.0);

  SDNchdyMX38 = histogramFactory().createHistogram1D
    ("SDNchdyMX38",200,-10.0,10.0);
  SDNchdyMX815 = histogramFactory().createHistogram1D
    ("SDNchdyMX815",200,-10.0,10.0);
  SDNchdyMX1530 = histogramFactory().createHistogram1D
    ("SDNchdyMX1530",200,-10.0,10.0);
  SDNchdyMX38Weights = histogramFactory().createHistogram1D
    ("SDNchdyMX38Weights",200,-10.0,10.0);
  SDNchdyMX815Weights = histogramFactory().createHistogram1D
    ("SDNchdyMX815Weights",200,-10.0,10.0);
  SDNchdyMX1530Weights = histogramFactory().createHistogram1D
    ("SDNchdyMX1530Weights",200,-10.0,10.0);

  SDdedetaMX38 = histogramFactory().createHistogram1D
    ("SDdedetaMX38",200,-10.0,10.0);
  SDdedetaMX818 = histogramFactory().createHistogram1D
    ("SDdedetaMX818",200,-10.0,10.0);
  SDdedetaMX1830 = histogramFactory().createHistogram1D
    ("SDdedetaMX1830",200,-10.0,10.0);
  SDdedetaMX38Weights = histogramFactory().createHistogram1D
    ("SDdedetaMX38Weights",200,-10.0,10.0);
  SDdedetaMX818Weights = histogramFactory().createHistogram1D
    ("SDdedetaMX818Weights",200,-10.0,10.0);
  SDdedetaMX1830Weights = histogramFactory().createHistogram1D
    ("SDdedetaMX1830Weights",200,-10.0,10.0);

  secondLastPT = histogramFactory().createHistogram1D
    ("secondLastPT",200,-10.0,10.0);
  etaDiffAriadne = 0.0;
  NGlueAriadne = 0.0;
  chargedPT = histogramFactory().createHistogram1D
    ("chargedPT",200,-10.0,10.0);
  chargedN = histogramFactory().createHistogram1D
    ("chargedN",200,-10.0,10.0);
  chargedNALICE = histogramFactory().createHistogram1D
    ("chargedNALICE",10,-1.0,1.0);
  INEL = 0.0;
  EtDensPHENIX = ZERO;
  totalWeight = 0.0;
  centralWeight = 0.0;
  midWeight = 0.0;
  over100Weight = 0.0;
  over100under250Weight = 0.0;
  over350Weight = 0.0;
  chargeMultALICE = histogramFactory().createHistogram1D
    ("chargeMultALICE",101,-0.5,100.5);
  maxPTBefore = histogramFactory().createHistogram1D
    ("maxPTBefore",200,0.0,100.0);
  DIPSYSumPT = histogramFactory().createHistogram1D
    ("DIPSYSumPT",200,0.0,100.0);
  finalSumET = histogramFactory().createHistogram1D
    ("finalSumET",200,0.0,100.0);

  //v2
  nParticipants = histogramFactory().createHistogram1D
    ("nParticipants",200,-1.0,199.5);

  RPEcc = histogramFactory().createHistogram1D
    ("RPEcc",200,-1.0,1.0);
  PartEcc = histogramFactory().createHistogram1D
    ("PartEcc",200,-1.0,1.0);
  OverlapArea = histogramFactory().createHistogram1D
    ("OverlapArea",200,0.0,10.0);

  EccArea = histogramFactory().createHistogram1D
    ("EccArea",200,0.0,10.0);
  AreaPart = histogramFactory().createHistogram1D
    ("AreaPart",200,-1.0,199.5);
  PE2 = histogramFactory().createHistogram1D
    ("PE2",200,-1.0,199.5);
  PE4 = histogramFactory().createHistogram1D
    ("PE4",200,-1.0,199.5);

  v2flow = histogramFactory().createHistogram1D
    ("v2flow",200, -1.0, 199.5);
  v4flow = histogramFactory().createHistogram1D
    ("v4flow",200, -1.0, 199.5);
  dNchFlow = histogramFactory().createHistogram1D
    ("dNchFlow",200, -1.0, 199.5);

  v2mix = histogramFactory().createHistogram1D
    ("v2mix",200, -1.0, 199.5);
  dNchMix = histogramFactory().createHistogram1D
    ("dNchMix",200, -1.0, 199.5);

  v2final = histogramFactory().createHistogram1D
    ("v2final",200, -1.0, 199.5);
  v2finalWeights = histogramFactory().createHistogram1D
    ("v2finalWeights",200, -1.0, 199.5);
  dNchFinal = histogramFactory().createHistogram1D
    ("dNchFinal",200, -1.0, 199.5);

  denseRatio = histogramFactory().createHistogram1D
    ("denseRatio",200, -1.0, 199.5);
  denseRatioWeights = histogramFactory().createHistogram1D
    ("denseRatioWeights",200, -1.0, 199.5);

  density = histogramFactory().createHistogram1D
    ("density",200, -1.0, 199.5);
  density1020 = histogramFactory().createHistogram1D
    ("density1020",200, -1.0, 199.5);
  density4060 = histogramFactory().createHistogram1D
    ("density4060",200, -1.0, 199.5);
  density60plus = histogramFactory().createHistogram1D
    ("density60plus",200, -1.0, 199.5);

  v2quad = histogramFactory().createHistogram1D
    ("v2quad",200, -1.0, 199.5);
  v2quadWeights = histogramFactory().createHistogram1D
    ("v2quadWeights",200, -1.0, 199.5);
  v2gap = histogramFactory().createHistogram1D
    ("v2gap",200, -1.0, 199.5);

  test2D = histogramFactory().createHistogram2D
    ("test2D",150, -10.0, 20.0, 100, -10.0, 10.0);
  test2DWeights = histogramFactory().createHistogram2D
    ("test2DWeights",150, -20.0, 10.0, 100, -10.0, 10.0);

  averageNParticipants = 0.0;
  //end v2

  //eccentricities
  centralGlueMult = 0.0;

  histNGlue = histogramFactory().createHistogram1D
    ("nGlue",100, 0.0, 20.0);
  histNGlueSqr = histogramFactory().createHistogram1D
    ("nGlueSqr",100, 0.0, 20.0);
  histNGlueWeights = histogramFactory().createHistogram1D
    ("nGlueWeights",100, 0.0, 20.0);

  histNGlueSpect = histogramFactory().createHistogram1D
    ("nGlueSpect",100, 0.0, 500.0);
  histNGlueSpectSqr = histogramFactory().createHistogram1D
    ("nGlueSpectSqr",100, 0.0, 500.0);
  histNGlueSpectWeights = histogramFactory().createHistogram1D
    ("nGlueSpectWeights",100, 0.0, 500.0);

  histNSpectator = histogramFactory().createHistogram1D
    ("nSpectator",100, 0.0, 20.0);
  histNSpectatorSqr = histogramFactory().createHistogram1D
    ("nSpectatorSqr",100, 0.0, 20.0);
  histNSpectatorWeights = histogramFactory().createHistogram1D
    ("nSpectatorWeights",100, 0.0, 20.0);

  histArea = histogramFactory().createHistogram1D
    ("area",200, 0.0, 600.0);
  histAreaWeights = histogramFactory().createHistogram1D
    ("areaWeights",200, 0.0, 600.0);

  //eps1
  histe1 = histogramFactory().createHistogram1D
    ("e1",200, 0.0, 600.0);
  histe1Weights = histogramFactory().createHistogram1D
    ("e1Weights",200, 0.0, 600.0);
  histe1N = histogramFactory().createHistogram1D
    ("e1N",200, 0.0, 600.0);
  histe1NWeights = histogramFactory().createHistogram1D
    ("e1NWeights",200, 0.0, 600.0);
  histe1W = histogramFactory().createHistogram1D
    ("e1W",200, 0.0, 600.0);
  histe1WWeights = histogramFactory().createHistogram1D
    ("e1WWeights",200, 0.0, 600.0);
  histphi1 = histogramFactory().createHistogram1D
    ("phi1",200, 0.0, M_PI);
  histphi1central = histogramFactory().createHistogram1D
    ("phi1central",200, 0.0, M_PI);
  histphi1mid = histogramFactory().createHistogram1D
    ("phi1mid",200, 0.0, M_PI);
  histphi1N = histogramFactory().createHistogram1D
    ("phi1N",200, 0.0, M_PI);
  histphi1W = histogramFactory().createHistogram1D
    ("phi1W",200, 0.0, M_PI);
  //eps2
  histe2 = histogramFactory().createHistogram1D
    ("e2",200, 0.0, 600.0);
  histe2Weights = histogramFactory().createHistogram1D
    ("e2Weights",200, 0.0, 600.0);
  histe2N = histogramFactory().createHistogram1D
    ("e2N",200, 0.0, 600.0);
  histe2NWeights = histogramFactory().createHistogram1D
    ("e2NWeights",200, 0.0, 600.0);
  histe2W = histogramFactory().createHistogram1D
    ("e2W",200, 0.0, 600.0);
  histe2WWeights = histogramFactory().createHistogram1D
    ("e2WWeights",200, 0.0, 600.0);
  histphi2 = histogramFactory().createHistogram1D
    ("phi2",200, 0.0, M_PI);
  histphi2central = histogramFactory().createHistogram1D
    ("phi2central",200, 0.0, M_PI);
  histphi2mid = histogramFactory().createHistogram1D
    ("phi2mid",200, 0.0, M_PI);
  histphi2N = histogramFactory().createHistogram1D
    ("phi2N",200, 0.0, M_PI);
  histphi2W = histogramFactory().createHistogram1D
    ("phi2W",200, 0.0, M_PI);
  //eps3
  histe3 = histogramFactory().createHistogram1D
    ("e3",200, 0.0, 600.0);
  histe3Weights = histogramFactory().createHistogram1D
    ("e3Weights",200, 0.0, 600.0);
  histe3N = histogramFactory().createHistogram1D
    ("e3N",200, 0.0, 600.0);
  histe3NWeights = histogramFactory().createHistogram1D
    ("e3NWeights",200, 0.0, 600.0);
  histe3W = histogramFactory().createHistogram1D
    ("e3W",200, 0.0, 600.0);
  histe3WWeights = histogramFactory().createHistogram1D
    ("e3WWeights",200, 0.0, 600.0);
  histphi3 = histogramFactory().createHistogram1D
     ("phi3",200, 0.0, M_PI);
  histphi3central = histogramFactory().createHistogram1D
     ("phi3central",200, 0.0, M_PI);
  histphi3mid = histogramFactory().createHistogram1D
     ("phi3mid",200, 0.0, M_PI);
  histphi3N = histogramFactory().createHistogram1D
    ("phi3N",200, 0.0, M_PI);
  histphi3W = histogramFactory().createHistogram1D
    ("phi3W",200, 0.0, M_PI);
  //eps4
  histe4 = histogramFactory().createHistogram1D
    ("e4",200, 0.0, 600.0);
  histe4Weights = histogramFactory().createHistogram1D
    ("e4Weights",200, 0.0, 600.0);
  histe4N = histogramFactory().createHistogram1D
    ("e4N",200, 0.0, 600.0);
  histe4NWeights = histogramFactory().createHistogram1D
    ("e4NWeights",200, 0.0, 600.0);
  histe4W = histogramFactory().createHistogram1D
    ("e4W",200, 0.0, 600.0);
  histe4WWeights = histogramFactory().createHistogram1D
    ("e4WWeights",200, 0.0, 600.0);
  histphi4 = histogramFactory().createHistogram1D
    ("phi4",200, 0.0, M_PI);
  histphi4central = histogramFactory().createHistogram1D
    ("phi4central",200, 0.0, M_PI);
  histphi4mid = histogramFactory().createHistogram1D
    ("phi4mid",200, 0.0, M_PI);
  histphi4N = histogramFactory().createHistogram1D
    ("phi4N",200, 0.0, M_PI);
  histphi4W = histogramFactory().createHistogram1D
    ("phi4W",200, 0.0, M_PI);
  //eps5
  histe5 = histogramFactory().createHistogram1D
    ("e5",200, 0.0, 600.0);
  histe5Weights = histogramFactory().createHistogram1D
    ("e5Weights",200, 0.0, 600.0);
  histe5N = histogramFactory().createHistogram1D
    ("e5N",200, 0.0, 600.0);
  histe5NWeights = histogramFactory().createHistogram1D
    ("e5NWeights",200, 0.0, 600.0);
  histe5W = histogramFactory().createHistogram1D
    ("e5W",200, 0.0, 600.0);
  histe5WWeights = histogramFactory().createHistogram1D
    ("e5WWeights",200, 0.0, 600.0);
  histphi5 = histogramFactory().createHistogram1D
    ("phi5",200, 0.0, M_PI);
  histphi5central = histogramFactory().createHistogram1D
    ("phi5central",200, 0.0, M_PI);
  histphi5mid = histogramFactory().createHistogram1D
    ("phi5mid",200, 0.0, M_PI);
  histphi5N = histogramFactory().createHistogram1D
    ("phi5N",200, 0.0, M_PI);
  histphi5W = histogramFactory().createHistogram1D
    ("phi5W",200, 0.0, M_PI);

  //eps^2
  histe1Sqrd = histogramFactory().createHistogram1D
    ("e1Sqrd",200, 0.0, 600.0);
  histe1SqrdWeights = histogramFactory().createHistogram1D
    ("e1SqrdWeights",200, 0.0, 600.0);
  histe2Sqrd = histogramFactory().createHistogram1D
    ("e2Sqrd",200, 0.0, 600.0);
  histe2SqrdWeights = histogramFactory().createHistogram1D
    ("e2SqrdWeights",200, 0.0, 600.0);
  histe3Sqrd = histogramFactory().createHistogram1D
    ("e3Sqrd",200, 0.0, 600.0);
  histe3SqrdWeights = histogramFactory().createHistogram1D
    ("e3SqrdWeights",200, 0.0, 600.0);
  histe4Sqrd = histogramFactory().createHistogram1D
    ("e4Sqrd",200, 0.0, 600.0);
  histe4SqrdWeights = histogramFactory().createHistogram1D
    ("e4SqrdWeights",200, 0.0, 600.0);
  histe5Sqrd = histogramFactory().createHistogram1D
    ("e5Sqrd",200, 0.0, 600.0);
  histe5SqrdWeights = histogramFactory().createHistogram1D
    ("e5SqrdWeights",200, 0.0, 600.0);

  //eps^4
  histe1Quad = histogramFactory().createHistogram1D
    ("e1Quad",200, 0.0, 600.0);
  histe1QuadWeights = histogramFactory().createHistogram1D
    ("e1QuadWeights",200, 0.0, 600.0);
  histe2Quad = histogramFactory().createHistogram1D
    ("e2Quad",200, 0.0, 600.0);
  histe2QuadWeights = histogramFactory().createHistogram1D
    ("e2QuadWeights",200, 0.0, 600.0);
  histe3Quad = histogramFactory().createHistogram1D
    ("e3Quad",200, 0.0, 600.0);
  histe3QuadWeights = histogramFactory().createHistogram1D
    ("e3QuadWeights",200, 0.0, 600.0);
  histe4Quad = histogramFactory().createHistogram1D
    ("e4Quad",200, 0.0, 600.0);
  histe4QuadWeights = histogramFactory().createHistogram1D
    ("e4QuadWeights",200, 0.0, 600.0);
  histe5Quad = histogramFactory().createHistogram1D
    ("e5Quad",200, 0.0, 600.0);
  histe5QuadWeights = histogramFactory().createHistogram1D
    ("e5QuadWeights",200, 0.0, 600.0);

  //for touched analysis
  //eps2
  histe2Touched = histogramFactory().createHistogram1D
    ("e2Touched",200, 0.0, 600.0);
  histe2TouchedWeights = histogramFactory().createHistogram1D
    ("e2TouchedWeights",200, 0.0, 600.0);
  //eps3
  histe3Touched = histogramFactory().createHistogram1D
    ("e3Touched",200, 0.0, 600.0);
  histe3TouchedWeights = histogramFactory().createHistogram1D
    ("e3TouchedWeights",200, 0.0, 600.0);
  //eps4
  histe4Touched = histogramFactory().createHistogram1D
    ("e4Touched",200, 0.0, 600.0);
  histe4TouchedWeights = histogramFactory().createHistogram1D
    ("e4TouchedWeights",200, 0.0, 600.0);
  //eps5
  histe5Touched = histogramFactory().createHistogram1D
    ("e5Touched",200, 0.0, 600.0);
  histe5TouchedWeights = histogramFactory().createHistogram1D
    ("e5TouchedWeights",200, 0.0, 600.0);

  //for (-3,-1) analysis
  //eps1
  histe1m3tom1 = histogramFactory().createHistogram1D
    ("e1m3tom1",200, 0.0, 600.0);
  histe1m3tom1Weights = histogramFactory().createHistogram1D
    ("e1m3tom1Weights",200, 0.0, 600.0);
  //eps2
  histe2m3tom1 = histogramFactory().createHistogram1D
    ("e2m3tom1",200, 0.0, 600.0);
  histe2m3tom1Weights = histogramFactory().createHistogram1D
    ("e2m3tom1Weights",200, 0.0, 600.0);
  //eps3
  histe3m3tom1 = histogramFactory().createHistogram1D
    ("e3m3tom1",200, 0.0, 600.0);
  histe3m3tom1Weights = histogramFactory().createHistogram1D
    ("e3m3tom1Weights",200, 0.0, 600.0);
  //eps4
  histe4m3tom1 = histogramFactory().createHistogram1D
    ("e4m3tom1",200, 0.0, 600.0);
  histe4m3tom1Weights = histogramFactory().createHistogram1D
    ("e4m3tom1Weights",200, 0.0, 600.0);
  //eps5
  histe5m3tom1 = histogramFactory().createHistogram1D
    ("e5m3tom1",200, 0.0, 600.0);
  histe5m3tom1Weights = histogramFactory().createHistogram1D
    ("e5m3tom1Weights",200, 0.0, 600.0);

  //for (1,3) analysis
  //eps1
  histe1p1top3 = histogramFactory().createHistogram1D
    ("e1p1top3",200, 0.0, 600.0);
  histe1p1top3Weights = histogramFactory().createHistogram1D
    ("e1p1top3Weights",200, 0.0, 600.0);
  //eps2
  histe2p1top3 = histogramFactory().createHistogram1D
    ("e2p1top3",200, 0.0, 600.0);
  histe2p1top3Weights = histogramFactory().createHistogram1D
    ("e2p1top3Weights",200, 0.0, 600.0);
  //eps3
  histe3p1top3 = histogramFactory().createHistogram1D
    ("e3p1top3",200, 0.0, 600.0);
  histe3p1top3Weights = histogramFactory().createHistogram1D
    ("e3p1top3Weights",200, 0.0, 600.0);
  //eps4
  histe4p1top3 = histogramFactory().createHistogram1D
    ("e4p1top3",200, 0.0, 600.0);
  histe4p1top3Weights = histogramFactory().createHistogram1D
    ("e4p1top3Weights",200, 0.0, 600.0);
  //eps5
  histe5p1top3 = histogramFactory().createHistogram1D
    ("e5p1top3",200, 0.0, 600.0);
  histe5p1top3Weights = histogramFactory().createHistogram1D
    ("e5p1top3Weights",200, 0.0, 600.0);

  //CoG distance
  histCoGdistance = histogramFactory().createHistogram1D
    ("CoGdistance",200, 0.0, 20.0);

  //correlations
  histeCorr1 = histogramFactory().createHistogram2D
    ("eCorr1",50, 0.0, 1.0, 50, 0.0, 1.0);
  histeCorr2 = histogramFactory().createHistogram2D
    ("eCorr2",50, 0.0, 1.0, 50, 0.0, 1.0);
  histeCorr3 = histogramFactory().createHistogram2D
    ("eCorr3",50, 0.0, 1.0, 50, 0.0, 1.0);
  histeCorr4 = histogramFactory().createHistogram2D
    ("eCorr4",50, 0.0, 1.0, 50, 0.0, 1.0);
  histeCorr5 = histogramFactory().createHistogram2D
    ("eCorr5",50, 0.0, 1.0, 50, 0.0, 1.0);
  histphiCorr1 = histogramFactory().createHistogram1D
    ("phiCorr1",200, 0.0, M_PI);
  histphiCorr2 = histogramFactory().createHistogram1D
    ("phiCorr2",200, 0.0, M_PI);
  histphiCorr2mid = histogramFactory().createHistogram1D
    ("phiCorr2mid",200, 0.0, M_PI);
  histphiCorr2central = histogramFactory().createHistogram1D
    ("phiCorr2central",200, 0.0, M_PI);
  histphiCorr3 = histogramFactory().createHistogram1D
    ("phiCorr3",200, 0.0, M_PI);
  histphiCorr4 = histogramFactory().createHistogram1D
    ("phiCorr4",200, 0.0, M_PI);
  histphiCorr5 = histogramFactory().createHistogram1D
    ("phiCorr5",200, 0.0, M_PI);

  //<eps(1,3)*eps(-3,-1)>
  epFepB1 = epF1 = epB1 = epSqrdF1 = epSqrdB1 = 0.0;
  epFepB2 = epF2 = epB2 = epSqrdF2 = epSqrdB2 = 0.0;
  epFepB2mid = epF2mid= epB2mid = epSqrdF2mid = epSqrdB2mid = 0.0;
  epFepB2central = epF2central= epB2central = epSqrdF2central = epSqrdB2central = 0.0;
  epFepB3 = epF3 = epB3 = epSqrdF3 = epSqrdB3 = 0.0;
  epFepB4 = epF4 = epB4 = epSqrdF4 = epSqrdB4 = 0.0;

  histe1e1p1top3 = histogramFactory().createHistogram1D
    ("e1e1p1top3",200, 0.0, 600.0);
  histe1e1p1top3Weights = histogramFactory().createHistogram1D
    ("e1e1p1top3Weights",200, 0.0, 600.0);
  histe2e2p1top3 = histogramFactory().createHistogram1D
    ("e2e2p1top3",200, 0.0, 600.0);
  histe2e2p1top3Weights = histogramFactory().createHistogram1D
    ("e2e2p1top3Weights",200, 0.0, 600.0);
  histe3e3p1top3 = histogramFactory().createHistogram1D
    ("e3e3p1top3",200, 0.0, 600.0);
  histe3e3p1top3Weights = histogramFactory().createHistogram1D
    ("e3e3p1top3Weights",200, 0.0, 600.0);
  histe4e4p1top3 = histogramFactory().createHistogram1D
    ("e4e4p1top3",200, 0.0, 600.0);
  histe4e4p1top3Weights = histogramFactory().createHistogram1D
    ("e4e4p1top3Weights",200, 0.0, 600.0);
  histe5e5p1top3 = histogramFactory().createHistogram1D
    ("e5e5p1top3",200, 0.0, 600.0);
  histe5e5p1top3Weights = histogramFactory().createHistogram1D
    ("e5e5p1top3Weights",200, 0.0, 600.0);

  //dN/deta
  histdNdnSpectdeta = histogramFactory().createHistogram2D
    ("dNdnSpectdeta",25, 0.0, 500.0, 100, -10.0, 10.0);

  histdNdeta0to100spect = histogramFactory().createHistogram1D
    ("dNdeta0to100spect",200, -10, 10);
  histdNdeta100to200spect = histogramFactory().createHistogram1D
    ("dNdeta100to200spect",200, -10, 10);
  histdNdeta200to300spect = histogramFactory().createHistogram1D
    ("dNdeta200to300spect",200, -10, 10);
  histdNdeta300to400spect = histogramFactory().createHistogram1D
    ("dNdeta300to400spect",200, -10, 10);

  spect0to100Weight = 0.0;
  spect100to200Weight = 0.0;
  spect200to300Weight = 0.0;
  spect300to400Weight = 0.0;

  //flux tubes
  histFluxTube = histogramFactory().createHistogram2D
    ("fluxTube",15, 0, 15.0, 100, 0.0, 2.0*M_PI);
  histFluxTubeWeights = histogramFactory().createHistogram2D
    ("fluxTubeWeights",15, 0, 15.0, 100, 0.0, 2.0*M_PI);
}


void FSAnalysis::dofinish() {
  AnalysisHandler::dofinish();
  double sumw = weights->sumBinHeights();
  double SDsumw = SDweights->sumBinHeights();
  double AU4SDsumw = AU4SDweights->sumBinHeights();
  double xsecw = generator()->integratedXSec()/millibarn;
  double totalWeightScale = 0.0;
  double over100WeightScale = 0.0;
  double over100under250WeightScale = 0.0;
  double over350WeightScale = 0.0;
  double spect0to100WeightScale = 0.0;
  double spect100to200WeightScale = 0.0;
  double spect200to300WeightScale = 0.0;
  double spect300to400WeightScale = 0.0;
  if ( totalWeight != 0.0 ) totalWeightScale = 1.0/totalWeight;
  if ( over100Weight != 0.0 ) over100WeightScale = 1.0/over100Weight;
  if ( over100under250Weight != 0.0 ) over100under250WeightScale = 1.0/over100under250Weight;
  if ( over350Weight != 0.0 ) over350WeightScale = 1.0/over350Weight;
  if ( spect0to100Weight != 0.0 ) spect0to100WeightScale = 1.0/spect0to100Weight;
  if ( spect100to200Weight != 0.0 ) spect100to200WeightScale = 1.0/spect100to200Weight;
  if ( spect200to300Weight != 0.0 ) spect200to300WeightScale = 1.0/spect200to300Weight;
  if ( spect300to400Weight != 0.0 ) spect300to400WeightScale = 1.0/spect300to400Weight;
  generator()->log()
    << setprecision(10)
    << "integrated xsec in microbarn: " << 1000.0*xsecw << endl;
  if ( xsecw <= 0.0 || sumw <= 0.0 ) return;
  normalize(ngbefore); 
  normalize(ngafter);
  if ( sumw > 0.0 ) weights->scale(1.0/sumw);
  if ( SDsumw > 0.0 ) SDweights->scale(1.0/SDsumw);
  if ( AU4SDsumw > 0.0 ) AU4SDweights->scale(1.0/SDsumw);
  normalize(epsRP, millibarn);
  normalize(epsPart, millibarn);
  normalize(epsPartOver20, millibarn);
  normalize(epsPartOver40, millibarn);
  normalize(mult12, millibarn);
  normalize(DIPSYPT, millibarn);
  DIPSYPT->scale( 1.0/generator()->integratedXSec()*(millibarn) );
  normalize(finalPT, millibarn);
  finalPT->scale( 1.0/generator()->integratedXSec()*(millibarn) );
  normalize(finalPTeta, millibarn);
  finalPTeta->scale( 1.0/generator()->integratedXSec()*(millibarn) );
  normalize(finalChargedPTeta, millibarn);
  finalChargedPTeta->scale( 1.0/generator()->integratedXSec()*(millibarn) );
  normalize(finalETeta, millibarn);
  finalETeta->scale( 1.0/generator()->integratedXSec()*(millibarn) );
  normalize(finalNcheta, millibarn);
  finalNcheta->scale( 1.0/generator()->integratedXSec()*(millibarn) );
  generator()->histogramFactory()->normalizeToXSecFraction(lnMX2distFrac);
  generator()->histogramFactory()->histogramFactory().
    divide("/lnMX2distNorm", *lnMX2dist, *lnMX2distWeights);
  generator()->histogramFactory()->histogramFactory().
    divide("/lnMX2dist0Norm", *lnMX2dist0, *lnMX2dist0Weights);
  generator()->histogramFactory()->histogramFactory().
    divide("/lnMX2dist1Norm", *lnMX2dist1, *lnMX2dist1Weights);
  generator()->histogramFactory()->histogramFactory().
    divide("/lnMX2dist2Norm", *lnMX2dist2, *lnMX2dist2Weights);
  normalize(SDNchdeta, millibarn);
  generator()->histogramFactory()->histogramFactory().
    divide("/MXmaxYNorm", *MXmaxY, *MXmaxYWeights);
  generator()->histogramFactory()->histogramFactory().
    divide("/MX2maxYNorm", *MX2maxY, *MX2maxYWeights);
  SDNchdeta->scale( 1.0/generator()->integratedXSec()*(millibarn) );

  generator()->histogramFactory()->histogramFactory().
    divide("/UA4dndetaMX20Norm", *UA4dndetaMX20, *UA4dndetaMX20Weights);
  generator()->histogramFactory()->histogramFactory().
    divide("/UA4dndetaMX60Norm", *UA4dndetaMX60, *UA4dndetaMX60Weights);
  generator()->histogramFactory()->histogramFactory().
    divide("/UA4dndetaMX80Norm", *UA4dndetaMX80, *UA4dndetaMX80Weights);
  generator()->histogramFactory()->histogramFactory().
    divide("/UA4dndetaMX100Norm", *UA4dndetaMX100, *UA4dndetaMX100Weights);
  generator()->histogramFactory()->histogramFactory().
    divide("/UA4dndetaMX115Norm", *UA4dndetaMX115, *UA4dndetaMX115Weights);
  generator()->histogramFactory()->histogramFactory().
    divide("/UA4dndetaMX140Norm", *UA4dndetaMX140, *UA4dndetaMX140Weights);

  generator()->histogramFactory()->histogramFactory().
    divide("/SDNchdyMX38Norm", *SDNchdyMX38, *SDNchdyMX38Weights);
  generator()->histogramFactory()->histogramFactory().
    divide("/SDNchdyMX815Norm", *SDNchdyMX815, *SDNchdyMX815Weights);
  generator()->histogramFactory()->histogramFactory().
    divide("/SDNchdyMX1530Norm", *SDNchdyMX1530, *SDNchdyMX1530Weights);

  generator()->histogramFactory()->histogramFactory().
    divide("/SDdedetaMX38Norm", *SDdedetaMX38, *SDdedetaMX38Weights);
  generator()->histogramFactory()->histogramFactory().
    divide("/SDdedetaMX818Norm", *SDdedetaMX818, *SDdedetaMX818Weights);
  generator()->histogramFactory()->histogramFactory().
    divide("/SDdedetaMX1830Norm", *SDdedetaMX1830, *SDdedetaMX1830Weights);


  normalize(secondLastPT, millibarn);
  secondLastPT->scale( 1.0/generator()->integratedXSec()*(millibarn) );
  generator()->log()
    << "average eta distance from FSR gluon to colour neighbour: "
    << etaDiffAriadne/NGlueAriadne << endl;
  normalize(chargedPT, millibarn);
  chargedPT->scale( 1.0/generator()->integratedXSec()*(millibarn) );
  normalize(chargedN, millibarn);
  chargedN->scale( 1.0/generator()->integratedXSec()*(millibarn) );
  chargedNALICE->scale(5.0/INEL);
  chargeMultALICE->scale(1.0/INEL);
  EtDensPHENIX *= totalWeightScale;
  generator()->log() << "average et density is " << EtDensPHENIX/GeV << endl;
  normalize(maxPTBefore, millibarn);
  normalize(DIPSYSumPT, millibarn);
  DIPSYSumPT->scale( 1.0/generator()->integratedXSec()*(millibarn) );
  normalize(finalSumET, millibarn);
  finalSumET->scale( 1.0/generator()->integratedXSec()*(millibarn) );

  //v2
  generator()->histogramFactory()->histogramFactory().
    divide("/AreaPartNorm", *AreaPart, *nParticipants);
  tH1DPtr PE2Norm = generator()->histogramFactory()->histogramFactory().
    divide("/PE2Norm", *PE2, *nParticipants);
  tH1DPtr PE4Norm = generator()->histogramFactory()->histogramFactory().
    divide("/PE4Norm", *PE4, *nParticipants);

  tH1DPtr PE22 = generator()->histogramFactory()->histogramFactory().
    multiply("/PE22", *PE2Norm, *PE2Norm);
  tH1DPtr PE22d = generator()->histogramFactory()->histogramFactory().
    add("/PE22d", *PE22, *PE22);
  generator()->histogramFactory()->histogramFactory().
    subtract("/eps4", *PE4Norm, *PE22d);

  generator()->histogramFactory()->histogramFactory().
    divide("/EccAreaNorm", *EccArea, *OverlapArea);

  normalize(RPEcc, millibarn);
  normalize(PartEcc, millibarn);
  normalize(OverlapArea, millibarn);

  tH1DPtr v22 = generator()->histogramFactory()->histogramFactory().
    divide("/v22Flow", *v2flow, *dNchFlow);
  generator()->histogramFactory()->histogramFactory().
    divide("/v22Mix", *v2mix, *dNchMix);
  generator()->histogramFactory()->histogramFactory().
    divide("/v2finalNorm", *v2final, *v2finalWeights);

  tH1DPtr v4flowNorm = generator()->histogramFactory()->histogramFactory().
    divide("/v42Flow", *v4flow, *dNchFlow);

  tH1DPtr v24 = generator()->histogramFactory()->histogramFactory().
    multiply("/v24", *v22, *v22);

  tH1DPtr temp = generator()->histogramFactory()->histogramFactory().
    subtract("/temp", *v4flowNorm, *v24);
  generator()->histogramFactory()->histogramFactory().
    subtract("/v4", *temp, *v24);

  generator()->histogramFactory()->histogramFactory().
    divide("/v2quadNorm", *v2quad, *v2quadWeights);

  generator()->histogramFactory()->histogramFactory().
    divide("/v2gapMean", *v2gap, *dNchFinal);

  generator()->histogramFactory()->histogramFactory().
    divide("/denseRatioNorm", *denseRatio, *denseRatioWeights);
  generator()->histogramFactory()->histogramFactory().
    divide("/densityNorm", *density, *dNchFlow);

  //SD
  generator()->histogramFactory()->histogramFactory().
    divide("/SDavNNorm", *SDavN, *SDavNWeights);
  generator()->histogramFactory()->histogramFactory().
    divide("/SDavN2Norm", *SDavN2, *SDavN2Weights);

  normalize(dNchFlow, millibarn);
  normalize(dNchMix, millibarn);
  normalize(dNchFinal, millibarn);
  normalize(nParticipants, millibarn);
  normalize(density1020, millibarn);
  normalize(density4060, millibarn);
  normalize(density60plus, millibarn);

  test2D->scale( totalWeightScale );

  generator()->log()
    << "average number of participants: "
    << averageNParticipants*totalWeightScale << endl
    << "average number of gluons in rapidity [-1,1]: "
    << centralGlueMult*totalWeightScale << endl;

  generator()->log()
    << "/totalWeight "
    << totalWeight << endl;
  //end v2

  //eccentricities
  generator()->histogramFactory()->histogramFactory().
    divide("/nGlueNorm", *histNGlue, *histNGlueWeights);
  generator()->histogramFactory()->histogramFactory().
    divide("/nGlueSqrNorm", *histNGlueSqr, *histNGlueWeights);

  generator()->histogramFactory()->histogramFactory().
    divide("/nGlueSpectNorm", *histNGlueSpect, *histNGlueSpectWeights);
  generator()->histogramFactory()->histogramFactory().
    divide("/nGlueSpectSqrNorm", *histNGlueSpectSqr, *histNGlueSpectWeights);

  generator()->histogramFactory()->histogramFactory().
    divide("/nSpectatorNorm", *histNSpectator, *histNSpectatorWeights);
  generator()->histogramFactory()->histogramFactory().
    divide("/nSpectatorSqrNorm", *histNSpectatorSqr, *histNSpectatorWeights);

  generator()->histogramFactory()->histogramFactory().
    divide("/areaNorm", *histArea, *histAreaWeights);

  //eps1
  generator()->histogramFactory()->histogramFactory().
    divide("/e1Norm", *histe1, *histe1Weights);
  generator()->histogramFactory()->histogramFactory().
    divide("/e1NNorm", *histe1N, *histe1NWeights);
  generator()->histogramFactory()->histogramFactory().
    divide("/e1WNorm", *histe1W, *histe1WWeights);
  histphi1->scale( over100WeightScale );
  histphi1central->scale( over350WeightScale );
  histphi1mid->scale( over100under250WeightScale );
  histphi1W->scale( over100WeightScale );
  histphi1N->scale( over100WeightScale );
  //eps2
  generator()->histogramFactory()->histogramFactory().
    divide("/e2Norm", *histe2, *histe2Weights);
  generator()->histogramFactory()->histogramFactory().
    divide("/e2NNorm", *histe2N, *histe2NWeights);
  generator()->histogramFactory()->histogramFactory().
    divide("/e2WNorm", *histe2W, *histe2WWeights);
  histphi2->scale( over100WeightScale );
  histphi2central->scale( over350WeightScale );
  histphi2mid->scale( over100under250WeightScale );
  histphi2W->scale( over100WeightScale );
  histphi2N->scale( over100WeightScale );
  //eps3
  generator()->histogramFactory()->histogramFactory().
    divide("/e3Norm", *histe3, *histe3Weights);
  generator()->histogramFactory()->histogramFactory().
    divide("/e3NNorm", *histe3N, *histe3NWeights);
  generator()->histogramFactory()->histogramFactory().
    divide("/e3WNorm", *histe3W, *histe3WWeights);
  histphi3->scale( over100WeightScale );
  histphi3central->scale( over350WeightScale );
  histphi3mid->scale( over100under250WeightScale );
  histphi3W->scale( over100WeightScale );
  histphi3N->scale( over100WeightScale );
  //eps4
  generator()->histogramFactory()->histogramFactory().
    divide("/e4Norm", *histe4, *histe4Weights);
  generator()->histogramFactory()->histogramFactory().
    divide("/e4NNorm", *histe4N, *histe4NWeights);
  generator()->histogramFactory()->histogramFactory().
    divide("/e4WNorm", *histe4W, *histe4WWeights);
  histphi4->scale( over100WeightScale );
  histphi4central->scale( over350WeightScale );
  histphi4mid->scale( over100under250WeightScale );
  histphi4W->scale( over100WeightScale );
  histphi4N->scale( over100WeightScale );
  //eps5
  generator()->histogramFactory()->histogramFactory().
    divide("/e5Norm", *histe5, *histe5Weights);
  generator()->histogramFactory()->histogramFactory().
    divide("/e5NNorm", *histe5N, *histe5NWeights);
  generator()->histogramFactory()->histogramFactory().
    divide("/e5WNorm", *histe5W, *histe5WWeights);
  histphi5->scale( over100WeightScale );
  histphi5central->scale( over350WeightScale );
  histphi5mid->scale( over100under250WeightScale );
  histphi5W->scale( over100WeightScale );
  histphi5N->scale( over100WeightScale );

  //eps^2
  generator()->histogramFactory()->histogramFactory().
    divide("/e1SqrdNorm", *histe1Sqrd, *histe1SqrdWeights);
  generator()->histogramFactory()->histogramFactory().
    divide("/e2SqrdNorm", *histe2Sqrd, *histe2SqrdWeights);
  generator()->histogramFactory()->histogramFactory().
    divide("/e3SqrdNorm", *histe3Sqrd, *histe3SqrdWeights);
  generator()->histogramFactory()->histogramFactory().
    divide("/e4SqrdNorm", *histe4Sqrd, *histe4SqrdWeights);
  generator()->histogramFactory()->histogramFactory().
    divide("/e5SqrdNorm", *histe5Sqrd, *histe5SqrdWeights);

  //eps^4
  generator()->histogramFactory()->histogramFactory().
    divide("/e1QuadNorm", *histe1Quad, *histe1QuadWeights);
  generator()->histogramFactory()->histogramFactory().
    divide("/e2QuadNorm", *histe2Quad, *histe2QuadWeights);
  generator()->histogramFactory()->histogramFactory().
    divide("/e3QuadNorm", *histe3Quad, *histe3QuadWeights);
  generator()->histogramFactory()->histogramFactory().
    divide("/e4QuadNorm", *histe4Quad, *histe4QuadWeights);
  generator()->histogramFactory()->histogramFactory().
    divide("/e5QuadNorm", *histe5Quad, *histe5QuadWeights);

  //eps for touched
  generator()->histogramFactory()->histogramFactory().
    divide("/e2TouchedNorm", *histe2Touched, *histe2TouchedWeights);
  generator()->histogramFactory()->histogramFactory().
    divide("/e3TouchedNorm", *histe3Touched, *histe3TouchedWeights);
  generator()->histogramFactory()->histogramFactory().
    divide("/e4TouchedNorm", *histe4Touched, *histe4TouchedWeights);
  generator()->histogramFactory()->histogramFactory().
    divide("/e5TouchedNorm", *histe5Touched, *histe5TouchedWeights);

  //eps in (-3,-1)
  generator()->histogramFactory()->histogramFactory().
    divide("/e1m3tom1Norm", *histe1m3tom1, *histe1m3tom1Weights);
  generator()->histogramFactory()->histogramFactory().
    divide("/e2m3tom1Norm", *histe2m3tom1, *histe2m3tom1Weights);
  generator()->histogramFactory()->histogramFactory().
    divide("/e3m3tom1Norm", *histe3m3tom1, *histe3m3tom1Weights);
  generator()->histogramFactory()->histogramFactory().
    divide("/e4m3tom1Norm", *histe4m3tom1, *histe4m3tom1Weights);
  generator()->histogramFactory()->histogramFactory().
    divide("/e5m3tom1Norm", *histe5m3tom1, *histe5m3tom1Weights);

  //eps in (1,3)
  generator()->histogramFactory()->histogramFactory().
    divide("/e1p1top3Norm", *histe1p1top3, *histe1p1top3Weights);
  generator()->histogramFactory()->histogramFactory().
    divide("/e2p1top3Norm", *histe2p1top3, *histe2p1top3Weights);
  generator()->histogramFactory()->histogramFactory().
    divide("/e3p1top3Norm", *histe3p1top3, *histe3p1top3Weights);
  generator()->histogramFactory()->histogramFactory().
    divide("/e4p1top3Norm", *histe4p1top3, *histe4p1top3Weights);
  generator()->histogramFactory()->histogramFactory().
    divide("/e5p1top3Norm", *histe5p1top3, *histe5p1top3Weights);

  //CoG distance
  histCoGdistance->scale( totalWeightScale );

  //correlations
  histeCorr1->scale( over100WeightScale );
  histeCorr2->scale( over100WeightScale );
  histeCorr3->scale( over100WeightScale );
  histeCorr4->scale( over100WeightScale );
  histeCorr5->scale( over100WeightScale );
  histphiCorr1->scale( over100WeightScale );
  histphiCorr2->scale( over100WeightScale );
  histphiCorr2mid->scale( over100under250WeightScale );
  histphiCorr2central->scale( over350WeightScale );
  histphiCorr3->scale( over100WeightScale );
  histphiCorr4->scale( over100WeightScale );
  histphiCorr5->scale( over100WeightScale );

  //<eps(1,3)*eps(-3,-1)>
  epFepB1 *= over100WeightScale;
  epF1 *= over100WeightScale;
  epB1 *= over100WeightScale;
  epSqrdF1 *= over100WeightScale;
  epSqrdB1 *= over100WeightScale;
  epFepB2 *= over100WeightScale;
  epF2 *= over100WeightScale;
  epB2 *= over100WeightScale;
  epSqrdF2 *= over100WeightScale;
  epSqrdB2 *= over100WeightScale;
  epFepB2mid *= over100under250WeightScale;
  epF2mid *= over100under250WeightScale;
  epB2mid *= over100under250WeightScale;
  epSqrdF2mid *= over100under250WeightScale;
  epSqrdB2mid *= over100under250WeightScale;
  epFepB2central *= over350WeightScale;
  epF2central *= over350WeightScale;
  epB2central *= over350WeightScale;
  epSqrdF2central *= over350WeightScale;
  epSqrdB2central *= over350WeightScale;
  epFepB3 *= over100WeightScale;
  epF3 *= over100WeightScale;
  epB3 *= over100WeightScale;
  epSqrdF3 *= over100WeightScale;
  epSqrdB3 *= over100WeightScale;
  epFepB4 *= over100WeightScale;
  epF4 *= over100WeightScale;
  epB4 *= over100WeightScale;
  epSqrdF4 *= over100WeightScale;
  epSqrdB4 *= over100WeightScale;
  if ( over100WeightScale != 0.0 ) {
    rho1 = (epFepB1-epF1*epB1)/(sqrt(epSqrdF1 - sqr(epF1))*sqrt(epSqrdB1 - sqr(epB1)));
    rho2 = (epFepB2-epF2*epB2)/(sqrt(epSqrdF2 - sqr(epF2))*sqrt(epSqrdB2 - sqr(epB2)));
    rho2mid = (epFepB2mid-epF2mid*epB2mid)/
      (sqrt(epSqrdF2mid - sqr(epF2mid))*sqrt(epSqrdB2mid - sqr(epB2mid)));
    rho2central = (epFepB2central-epF2central*epB2central)/
      (sqrt(epSqrdF2central - sqr(epF2central))*sqrt(epSqrdB2central - sqr(epB2central)));
    rho3 = (epFepB3-epF3*epB3)/(sqrt(epSqrdF3 - sqr(epF3))*sqrt(epSqrdB3 - sqr(epB3)));
    rho4 = (epFepB4-epF4*epB4)/(sqrt(epSqrdF4 - sqr(epF4))*sqrt(epSqrdB4 - sqr(epB4)));
    // cout << "rho1 = " << rho1 << endl;
    // cout << "rho2 = " << rho2 << endl;
    // cout << "rho2mid = " << rho2mid << endl;
    // cout << "rho2central = " << rho2central << endl;
    // cout << "rho3 = " << rho3 << endl;
    // cout << "rho4 = " << rho4 << endl;
  }

  generator()->histogramFactory()->histogramFactory().
    divide("/e1e1p1top3Norm", *histe1e1p1top3, *histe1e1p1top3Weights);
  generator()->histogramFactory()->histogramFactory().
    divide("/e2e2p1top3Norm", *histe2e2p1top3, *histe2e2p1top3Weights);
  generator()->histogramFactory()->histogramFactory().
    divide("/e3e3p1top3Norm", *histe3e3p1top3, *histe3e3p1top3Weights);
  generator()->histogramFactory()->histogramFactory().
    divide("/e4e4p1top3Norm", *histe4e4p1top3, *histe4e4p1top3Weights);
  generator()->histogramFactory()->histogramFactory().
    divide("/e5e5p1top3Norm", *histe5e5p1top3, *histe5e5p1top3Weights);

  //dN/deta/dnSpect
  histdNdnSpectdeta->scale( totalWeightScale );

  histdNdeta0to100spect->scale( spect0to100WeightScale );
  histdNdeta100to200spect->scale( spect100to200WeightScale );
  histdNdeta200to300spect->scale( spect200to300WeightScale );
  histdNdeta300to400spect->scale( spect300to400WeightScale );

  //flux tube analysis
  generator()->histogramFactory()->histogramFactory().
    divide("/fluxTubeNorm", *histFluxTube, *histFluxTubeWeights);
}

void FSAnalysis::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void FSAnalysis::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}


// Static variable needed for the type description system in ThePEG.
#include "ThePEG/Utilities/DescribeClass.h"
DescribeClass<FSAnalysis,AnalysisHandler>
  describeDIPSYFSAnalysis("DIPSY::FSAnalysis", "FSAnalysis.so");

void FSAnalysis::Init() {

  static ClassDocumentation<FSAnalysis> documentation
    ("There is no documentation for the FSAnalysis class");

}

