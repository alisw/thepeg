// Header file for the TimeParticle and TimeShower classes.

#ifndef TimeShower_H
#define TimeShower_H

#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>
#include "Pythia7/Config/Pythia7.h"
#include "Basics.h"
#include "Beam.h"

namespace Pythia7 {
namespace Shower {

using namespace std;
 
//**************************************************************************


/**
 * This class holds info on a timelike particle in the shower evolution.
 * It derives from the Particle class.
 * Note that specific shower info is initialized in an off state.
 * Used by the internal Pythia7 Shower classes.
 */
class TimeParticle : public Particle {

public:
  /**
   * Constructors.
   */
  TimeParticle()
    : canBranch(false), hasBranched(false),
      shouldEvolveMore(false), idDaughter(0), copyOfPrimary(-1), coneSide(0),
      sister(-1), Q2Now(0.0), eMax(0.0), eNow(0.0), pAbsNow(0.0),
      zDaughter(0.0), phiDaughter(0.0), hasColVec4(false),
      hasAntiColVec4(false) {}  
  /** NOT DOCUMENTED */
  TimeParticle(long idin, long statusin = 0, long mother1in = -1, 
    long mother2in = -1, long colin = 0, long anticolin = 0, 
    Vec4 pin = Vec4(0.,0.,0.,0.), double min = 0., double scalein = 0.)
    : Particle(idin, statusin, mother1in, mother2in, colin, anticolin, pin, 
	       min, scalein),
      canBranch(false), hasBranched(false), shouldEvolveMore(false),
      idDaughter(0), copyOfPrimary(-1), coneSide(0), sister(-1), Q2Now(0.0),
      eMax(0.0), eNow(0.0), pAbsNow(0.0), zDaughter(0.0), phiDaughter(0.0),
      hasColVec4(false), hasAntiColVec4(false)  {}
  /** NOT DOCUMENTED */
  TimeParticle(const TimeParticle& pt)
    : Particle(pt), canBranch(pt.canBranch), hasBranched(pt.hasBranched),
      shouldEvolveMore(pt.shouldEvolveMore), idDaughter(pt.idDaughter),
      copyOfPrimary(pt.copyOfPrimary), coneSide(pt.coneSide), sister(pt.sister),
      Q2Now(pt.Q2Now), eMax(pt.eMax), eNow(pt.eNow), pAbsNow(pt.pAbsNow),
      zDaughter(pt.zDaughter), phiDaughter(pt.phiDaughter),
      hasColVec4(pt.hasColVec4), hasAntiColVec4(pt.hasAntiColVec4),
      colVec4(pt.colVec4), antiColVec4(pt.antiColVec4) {}  
  /** NOT DOCUMENTED */
  TimeParticle& operator=(const TimeParticle& pt) { 
    if (this != &pt) { Particle::operator=(pt); 
    canBranch = pt.canBranch; hasBranched = pt.hasBranched; 
    shouldEvolveMore = pt.shouldEvolveMore; idDaughter = pt.idDaughter; 
    copyOfPrimary = pt.copyOfPrimary; coneSide = pt.coneSide; 
    sister = pt.sister; Q2Now = pt.Q2Now; eMax = pt.eMax; eNow = pt.eNow; 
    pAbsNow = pt.pAbsNow; zDaughter = pt.zDaughter; 
    phiDaughter = pt.phiDaughter; hasColVec4 = pt.hasColVec4; 
    hasAntiColVec4 = pt.hasAntiColVec4; colVec4 = pt.colVec4; 
    antiColVec4 = pt.antiColVec4;}  return *this; } 
  /** NOT DOCUMENTED */
  TimeParticle(const Particle& pt) : Particle(pt) {
    canBranch = false;  hasBranched = false; sister = -1; 
    hasColVec4 = false; hasAntiColVec4 = false;}  
  /** NOT DOCUMENTED */
  TimeParticle& operator=(const Particle& pt) {Particle::operator=(pt);   
    canBranch = false;  hasBranched = false; sister = -1; 
    hasColVec4 = false; hasAntiColVec4 = false; return *this; }  


  /**
   * Member functions.
   */
  double mNow() const {if (hasBranched) return sqrt(m()*m() + Q2Now); 
    else return m(); }
  /** NOT DOCUMENTED */
  double m2Now() const {if (hasBranched) return m()*m() + Q2Now; 
    else return m()*m(); }
  /** NOT DOCUMENTED */
  double mOff() const {return sqrt(m()*m() + Q2Now); }
  /** NOT DOCUMENTED */
  double m2Off() const {return m()*m() + Q2Now; }
  /** NOT DOCUMENTED */
  void setColVec4() {hasColVec4 = false;}
  /** NOT DOCUMENTED */
  void setColVec4(Vec4 colVec4in) 
    {colVec4 = colVec4in; hasColVec4 = true;}
  /** NOT DOCUMENTED */
  void setAntiColVec4() {hasAntiColVec4 = false;}
  /** NOT DOCUMENTED */
  void setAntiColVec4(Vec4 antiColVec4in) 
    {antiColVec4 = antiColVec4in; hasAntiColVec4 = true;}
  /** NOT DOCUMENTED */
  void rot(double theta, double phi) {
    Particle::rot(theta, phi); 
    if (hasColVec4) colVec4.rot(theta, phi);
    if (hasAntiColVec4) antiColVec4.rot(theta, phi); }
  /** NOT DOCUMENTED */
  void bst(double betaX, double betaY, double betaZ) {
    Particle::bst(betaX, betaY, betaZ);
    if (hasColVec4) colVec4.bst(betaX, betaY, betaZ);  
    if (hasAntiColVec4) antiColVec4.bst(betaX, betaY, betaZ); }
  /** NOT DOCUMENTED */
  void bst(const Vec4& vec) {
    Particle::bst(vec);
    if (hasColVec4) colVec4.bst(vec);  
    if (hasAntiColVec4) antiColVec4.bst(vec); }
  /** NOT DOCUMENTED */
  void rotbst(const RotBstMatrix& M) {
    Particle::rotbst(M); 
    if (hasColVec4) colVec4.rotbst(M);
    if (hasAntiColVec4) antiColVec4.rotbst(M); }


  /**
   * Private members to be accessible from TimeShower.
   */
  friend class TimeShower; 


private:
  /** @cond NEVERTOBEDOCUMENTED */
  /** NOT DOCUMENTED */
  bool canBranch, hasBranched, shouldEvolveMore;
  /** NOT DOCUMENTED */
  long idDaughter, copyOfPrimary, coneSide, sister;
  /** NOT DOCUMENTED */
  double Q2Now, eMax, eNow, pAbsNow, zDaughter, phiDaughter; 
  /** NOT DOCUMENTED */
  bool hasColVec4, hasAntiColVec4;
  /** NOT DOCUMENTED */
  Vec4 colVec4, antiColVec4;
  /** @endcond */
};
 
//**************************************************************************


/**
 * The TimeShower class does timelike showers.
 * Used by the internal Pythia7 Shower classes.
 */
class TimeShower {

public:
  /** @cond NEVERTOBEDOCUMENTED */
  /**
   * Constants (possibly to be changed externally).
   */
  static long ANGULARORDER, NQUARK, ALPHASMODE, MEMODE, QEDSHOWER, 
    INITIALCONE, PHIPOLASYM, PHICOHERASYM, RESPECTSCALE;
  /**
   * Constants (possibly to be changed externally).
   */
  static double Q0, Q0CHGQ,Q0CHGL, ALPHASFIX, LAMBDA5, ALPHAEMFIX, 
    Q0FRACPS;
  /** @endcond */

  /**
   * Constructor.
   */
  TimeShower(long capacity = 100) {entry.reserve(capacity);}


  /**
   * Top-level driver routine to do a single time-like shower.
   */
  void shower(Event&, vector<long> = vector<long>(1,long(-1)), 
    double = -1., long = 0, long = 0, double = 0.5);

private:
  /** @cond NEVERTOBEDOCUMENTED */
  /** NOT DOCUMENTED */
  vector<TimeParticle> entry;
  /** NOT DOCUMENTED */
  long nPrimary, maxColIndx, inFlavour1, inFlavour2;
  /** NOT DOCUMENTED */
  RotBstMatrix bstMother;
  /** NOT DOCUMENTED */
  bool hasME, MEorder, MEsplit, MEgluinoDau;
  /** NOT DOCUMENTED */
  long MEkind, MEcombi;
  /** NOT DOCUMENTED */
  double MEmix;
  /** @endcond */
 
  /**
   * Read in info on system to be treated.
   */
  void read(Event&, vector<long>);
  /**
   * Write back shower after treatment.
   */
  void write(Event&);
  /**
   * Set up primary partons for evolution.
   */
  bool setUpPrimary(double, long, long, double);
  /**
   * Check and set kinematics for primary partons after evolution.
   */
  bool kinemPrimary();
  /**
   * Set up daughters of parton branching for subsequent evolution.
   */
  void setUpBranching(long);
  /**
   * Check and set kinematics for branching after daughter evolution.
   */
  bool kinemBranching(long, long, long);
  /**
   * Evolve a single parton; main physics routine.
   */
  void evolveParton(long);
  /**
   * Pick one of the partons for further evolution when required.
   */
  long pickParton(long, long);
  /**
   * Check that z and Q2 choices of branching are consistent.
   */
  bool zQcheck(long);
  /**
   * Construct branching kinematics for setUpBranching & kinemBranching.
   */
  void kinemConstruct(long, long, long);
  /**
   * Find class of QCD ME correction.
   */
  void findMEkind(long, long, double);
  /**
   * Find type of particle; used by findMEkind.
   */
  long findMEparticle(long);
  /**
   * Find mixture of V and A in gamma/Z: energy- and flavour-dependent. 
   */
  double gammaZmix();
  /**
   * Set up to calculate QCD ME correction with calcMEcorr.
   */
  double findMEcorr(long);
  /**
   * Calculate value of QCD ME correction.
   */
  static double calcMEcorr(long, long, double, double, double, double, double);

};

//**************************************************************************

}
}

#endif // TimeShower_H
