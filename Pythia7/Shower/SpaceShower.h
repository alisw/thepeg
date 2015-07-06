// Header file for the SpaceParticle and SpaceShower classes.

#ifndef SpaceShower_H
#define SpaceShower_H

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
 * This class holds info on a spacelike particle in the shower evolution.
 * It derives from the Particle class.
 * Note that specific shower info is initialized in an off state.
 * Used by the internal Pythia7 Shower classes.
 */
class SpaceParticle : public Particle {

public:
  /**
   * Constructors.
   */
  SpaceParticle()
    : Particle(), canBranch(false), hasBranched(false),
      shouldEvolveMore(false), idMother(-1), idSister(-1),
      coneSide(0), copyOfPrimary(-1), Q2Now(0.0), x(0.0),
      zMother(0.0), phiDaughter(0.0), sin2thetaMax(0.0),
      sin2thetaNow(0.0), hasColVec4(false), hasAntiColVec4(false) {}  
  /** NOT DOCUMENTED */
  SpaceParticle(long idin, long statusin = 0, long mother1in = -1, 
		long mother2in = -1, long colin = 0, long anticolin = 0, 
		Vec4 pin = Vec4(0.,0.,0.,0.),
		double min = 0., double scalein = 0.)
    : Particle(idin, statusin, mother1in, mother2in,
	       colin, anticolin, pin, min), canBranch(false),
      hasBranched(false), shouldEvolveMore(false), idMother(-1),
      idSister(-1), coneSide( 0), copyOfPrimary(-1), Q2Now(0.0),
      x(0.0), zMother(0.0), phiDaughter(0.0), sin2thetaMax(0.0),
      sin2thetaNow(0.0), hasColVec4(false), hasAntiColVec4(false) {}
  /** NOT DOCUMENTED */
  SpaceParticle(const SpaceParticle& pt)
    : Particle(pt), canBranch(pt.canBranch), hasBranched(pt.hasBranched),
      shouldEvolveMore(pt.shouldEvolveMore), idMother(pt.idMother),
      idSister(pt.idSister), coneSide(pt.coneSide),
      copyOfPrimary(pt.copyOfPrimary), Q2Now(pt.Q2Now), x(pt.x),
      zMother(pt.zMother), phiDaughter(pt.phiDaughter),
      sin2thetaMax(pt.sin2thetaMax), sin2thetaNow(pt.sin2thetaNow),
    hasColVec4(pt.hasColVec4), hasAntiColVec4(pt.hasAntiColVec4),
    colVec4(pt.colVec4), antiColVec4(pt.antiColVec4) {}  
  /** NOT DOCUMENTED */
  SpaceParticle& operator=(const SpaceParticle& pt) {
    if (this != &pt) { Particle::operator=(pt); 
    canBranch = pt.canBranch; hasBranched = pt.hasBranched; 
    shouldEvolveMore = pt.shouldEvolveMore; idMother = pt.idMother;
    idSister = pt.idSister; coneSide = pt.coneSide; 
    copyOfPrimary = pt.copyOfPrimary; Q2Now = pt.Q2Now;
    x = pt.x; zMother = pt.zMother; phiDaughter = pt.phiDaughter; 
    sin2thetaMax = pt.sin2thetaMax; sin2thetaNow = pt.sin2thetaNow; 
    hasColVec4 = pt.hasColVec4; hasAntiColVec4 = pt.hasAntiColVec4; 
    colVec4 = pt.colVec4; antiColVec4 = pt.antiColVec4;} return *this; } 
  /**
   * Construct a SpaceParticle from a Particle.
   */
  SpaceParticle(const Particle& pt) : Particle(pt) {
    canBranch = false; hasColVec4 = false; hasAntiColVec4 = false;}  
  /** NOT DOCUMENTED */
  SpaceParticle& operator=(const Particle& pt) {Particle::operator=(pt);   
    canBranch = false; hasColVec4 = false; hasAntiColVec4 = false;
    return *this; }  


  /**
   * Member functions.
   */
  double mNow() const {if (!hasBranched) return 0.; return -sqrt(Q2Now);}
  /** NOT DOCUMENTED */
  double m2Now() const {if (!hasBranched) return 0.; return -Q2Now;}
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
   * Private members to be accessible from SpaceShower.
   */
  friend class SpaceShower; 


private:
  /** @cond NEVERTOBEDOCUMENTED */
  /** NOT DOCUMENTED */
  bool canBranch, hasBranched, shouldEvolveMore;
  /** NOT DOCUMENTED */
  long idMother, idSister, coneSide, copyOfPrimary;
  /** NOT DOCUMENTED */
  double Q2Now, x, zMother, phiDaughter, sin2thetaMax, sin2thetaNow;  
  /** NOT DOCUMENTED */
  bool hasColVec4, hasAntiColVec4;
  /** NOT DOCUMENTED */
  Vec4 colVec4, antiColVec4;
  /** @endcond */
};
 
//**************************************************************************


/**
 * The SpaceShower class does spacelike showers.
 * Used by the internal Pythia7 Shower classes.
 */
class SpaceShower {

public:
  /** @cond NEVERTOBEDOCUMENTED */
  /**
   * Constants (possibly to be changed externally).
   */
  static long HADRONSHOWER, LEPTONSHOWER, ANGULARORDER, NQUARK, ALPHASMODE, 
    Q2ORDER, MAXVIRTUALITY, MEMODE, SOFTGLUONRESUM, FINALCONE, PHIPOLASYM, 
    PHICOHERASYM, RESPECTSCALE;
  /**
   * Constants (possibly to be changed externally).
   */
  static double Q0, Q0CHGQ,Q0CHGL, ALPHASFIX, LAMBDA5, ALPHAEMFIX, 
    EMINEMITTED, ZMINEMITTED, XMINEMITTEDCHG, TINYQCHG, TINYPDF, 
    TINYKERNELPDF, TINYKINPREC, HEAVYEVOL, EXTRAQEDPREWT, HEAVYXMAX, 
    Q2STARTFRAC;
  /** @endcond */

  /**
   * Constructor.
   */
  SpaceShower(long capacity = 20) {entry1.reserve(capacity); 
    entry2.reserve(capacity);}


  /**
   * Top-level driver routine to do a single time-like shower.
   */
  void shower(Event&, BeamParticle&, BeamParticle&, long = -1, long = -1,
    long = -1);

private: 
  /** @cond NEVERTOBEDOCUMENTED */
  /** NOT DOCUMENTED */
  BeamParticle* beam1;
  /** NOT DOCUMENTED */
  BeamParticle* beam2;
  /** NOT DOCUMENTED */
  vector<SpaceParticle> entry1;
  /** NOT DOCUMENTED */
  vector<SpaceParticle> entry2;
  /** NOT DOCUMENTED */
  bool hasME;
  /** NOT DOCUMENTED */
  long in1, in2, maxColIndx, MEkind, nFinal, finalId1, finalId2;
  /** NOT DOCUMENTED */
  double s, eCM, sHat, mHat, sHatFinalState, avgMT2, minMT2;
  /** @endcond */

  /**
   * Read in info on system to be treated.
   */
  void read(Event&, BeamParticle&, BeamParticle&, long, long);
  /**
   * Write back shower after treatment.
   */
  void write(Event&);
  /**
   * Set up primary partons for evolution.
   */
  bool setUpPrimary(long);
  /**
   * Check and set kinematics for primary partons after evolution.
   */
  bool kinemPrimary();
  /**
   * Pick one of the sides for further evolution when required.
   */
  long pickSide();
  /**
   * Set up daughters of parton branching for subsequent evolution.
   */
  void setUpBranching(long);
  /**
   * Check and set kinematics for branching after daughter evolution.
   */
  bool kinemBranching(long);
  /**
   * Evolve a single parton; main physics routine.
   */
  void evolveParton(long);
  /**
   * Find class of ME correction.
   */
  void findMEkind(long);
  /**
   * Provide maximum of expected ME weight; for preweighting of evolution.
   */
  double calcMEmax(long, long);
  /**
   * Provide actual ME weight for current branching.
   */
  double calcMEcorr(long, long, double, double); 

};
 
//**************************************************************************

}
}

#endif // SpaceShower_H
