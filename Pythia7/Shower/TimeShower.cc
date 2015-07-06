// Function definitions (not found in the header) for the
// TimeShower class, which does timelike showers. 

#include "Basics.h"
#include "Beam.h"
#include "Shower.h"
#include "TimeShower.h"

using namespace Pythia7;
using namespace Shower;

//**************************************************************************

// Angular ordering; = 0: off, = 1: on for g emission, 
// = 2: on also for g -> q qbar splitting.
long TimeShower::ANGULARORDER = 2;

// Number of allowed quark flavours in g -> q qbar branching.
long TimeShower::NQUARK = 5;

// Running of alpha_strong in evolution:
// = 0: fixed; = 1: scale Q^2/4; = 2: scale pT^2; 
// = 3: scale pT2, except for g -> q qbar, where it is Q^2/4 . 
long TimeShower::ALPHASMODE = 2;

// Use of matrix element corrections:
// = 0: no; = 1: yes.
long TimeShower::MEMODE = 1;

// Also allow a QED shower together with QCD ones: = 0: no; 
// = 1: radiate on-shell photons; = 2 : also allow photons to branch.
long TimeShower::QEDSHOWER = 2;

// Restrict first emission within cone given by colour flow in hard process.
// = 0: no; = 1: yes, isotropic phi angle inside cone; 
// = 2: yes, also with anisotropic phi angle inside cone.
long TimeShower::INITIALCONE = 2;

// Azimuthal asymmetry induced by gluon polarization.
// = 0: no; = 1: yes.
long TimeShower::PHIPOLASYM = 1;

// Azimuthal asymmetry induced by colour coherence.
// = 0: no; = 1: yes.
long TimeShower::PHICOHERASYM = 1;

// Use the scale variable of original partons to restrict branchings.
// = 0: no; = 1: yes, the Q2 < scale; = 2: yes, the pT2 < scale,
// = 3: yes, the (E*theta)^2 < scale; = 4: yes, the theta^2 < scale.
// (In all cases relations are approximate.) 
long TimeShower::RESPECTSCALE = 0;

// Parton shower cut-off mass for QCD emissions.
double TimeShower::Q0 = 1.;

// Parton shower cut-off mass for photon coupling to coloured particle.
double TimeShower::Q0CHGQ = 1.;

// Parton shower cut-off mass for pure QED branchings. Assumed <= Q0CHGQ.
double TimeShower::Q0CHGL = 0.001; 

// Fixed alpha_strong value for ALPHASMODE == 0. 
double TimeShower::ALPHASFIX = 0.2;

// Lambda_QCD(five flavours) in alpha_strong for ALPHASMODE >= 1 .
double TimeShower::LAMBDA5 = 0.25;

// Fixed alpha_em value. 
double TimeShower::ALPHAEMFIX = 0.0073;

// Fraction of Q0 cut-off mass used as safety margin in 
// daughter mass sum. Relevant for total parton multiplicity.
double TimeShower::Q0FRACPS = 0.25;

//*********

// Top-level driver routine to do a single time-like shower.

void TimeShower::shower(Event& event, vector<long> primary, 
  double Q2maxin, long MEkindin, long MEcombiin, double MEmixin) {

  // Read in info on system to be treated.
  read(event, primary);
  
  // Prepare to do shower; skip if no particles may shower.
  if ( setUpPrimary(Q2maxin, MEkindin, MEcombiin, MEmixin) ) {

    // Evolve primary partons (>= 2) independently of each other.
    for (long i = 1; i <= nPrimary; ++i) {
      if (entry[i].canBranch) evolveParton(i); 
    }

    // Check correlated set of primary parton branchings; 
    // construct primary kinematics if consistent.
    while ( !kinemPrimary() ) {

      // Else pick one parton and evolve it further; iterate as required.
      long i = pickParton(1, nPrimary);
      evolveParton(i);
    }   

    // Now loop over all partons that branched, until end of shower.
    long iMax = nPrimary;
    long iNow = 0;
    while (iNow < iMax) {
      ++iNow;
      if (entry[iNow].hasBranched) {

        // Construct kinematics of previously selected branching,
        // initially assuming the daughters are on mass shell. 
        setUpBranching(iNow);
        iMax += 2; 

        // This gives a pair of sister partons that is to be evolved.
        evolveParton(iMax-1); 
        evolveParton(iMax); 

        // If kinematics is inconsistent, evolve one parton further.
        while ( !kinemBranching(iNow, iMax-1, iMax) ) {
          long i = pickParton(iMax-1, iMax);
          evolveParton(i);
        }
      }
    }  

  }
  // Write back the generated shower.
  write(event);
  
}

//*********

// Read in info on system to be treated.

void TimeShower::read(Event& event, vector<long> primary) {

  // Primary particles either given as input or are all with status = 1.
  if (primary.size() > 0 && primary[0] < 0) primary.resize(0);
  if (primary.size() == 0) {
    for (long i = 0; i < event.size(); ++i) {
      if (event[i].status() == 1) primary.push_back(i); 
    }
  } 
  nPrimary = primary.size();
  entry.resize(nPrimary+1);
  TimeParticle& origin = entry[0];

  // Read in primary particles from external event.
  Vec4 pSum;
  for (long i = 1; i <= nPrimary; ++i) {
    TimeParticle& current = entry[i];
    long inev = primary[i-1];
    current = event[inev];
    current.mother2(inev);
    pSum += current.p();

    // Attach vectors of colour and anticolour mother or partner.
    if (current.col() > 0) {
      for (long j = 0; j < inev; ++j) {
        if (event[j].status() < 0 && event[j].col() == current.col() )
          current.setColVec4(event[j].p());
      }
      for (long j = 0; j < nPrimary; ++j) {
        if (event[primary[j]].anticol() == current.col() )
          current.setColVec4(event[primary[j]].p());
      }
    }
    if (current.anticol() > 0) {
      for (long j = 0; j < inev; ++j) {
        if (event[j].status() < 0 && event[j].anticol() == current.anticol() )
          current.setAntiColVec4(event[j].p());
      }
      for (long j = 0; j < nPrimary; ++j) {
        if (event[primary[j]].col() == current.anticol() )
          current.setAntiColVec4(event[primary[j]].p());
      }
    }
  }

  // Fill line 0 with one mother or sum of all primary.
  double eCM = pSum.mcalc(); 
  long mother1 = entry[1].mother1();
  long mother2 = entry[1].mother2();
  origin = (mother1 >= 0 && (mother2 == mother1 || mother2 < 0)) 
    ? event[mother1] : Particle( 0, 2, -1, -1, 0, 0, pSum, eCM );
  
  // Boost to rest frame of showering system.
  bstMother.reset(); 
  bstMother.bstback(origin.p());
  for (long i = 0; i < long(entry.size()); ++i) {entry[i].rotbst(bstMother);}
  bstMother.invert();

  // Read current largest colour index and incoming flavours.
  maxColIndx = event.colIndx();
  inFlavour1 = inFlavour2 = 0;
  if (mother2 != mother1 && mother1 >= 0 && mother2 >= 0) {
    inFlavour1 = event[mother1].id();
    inFlavour2 = event[mother2].id();
  } else if (origin.id() != 0) {
    long grandmother1 = event[mother1].mother1();
    if (grandmother1 >= 0) inFlavour1 = event[grandmother1].id(); 
    long grandmother2 = event[mother1].mother2();
    if (grandmother2 >= 0) inFlavour2 = event[grandmother2].id(); 
  }
  
}     

//*********

// Write back shower after treatment.

void TimeShower::write(Event& event) {

  // Mark primary partons as treated.
  for (long i = 1; i <= nPrimary; ++i) {
    long iPrim = entry[i].mother2();
    if (event[iPrim].status() > 0) event[iPrim].addstatus(3);
  }

  // Boost shower particles, shift mother flags, and add to event record.
   long shift = event.size() - 1;
  for (long i = 1; i < long(entry.size()); ++i) {
    TimeParticle& current = entry[i]; 
    Particle temp = current;
    temp.rotbst(bstMother); 
    temp.m(current.mNow());
    long shifted = (i <= nPrimary) ? current.mother2() 
      : current.mother1() + shift;
    if (current.copyOfPrimary > 0 && current.copyOfPrimary <= nPrimary) 
      temp.mothers(shifted, shifted);
    else temp.mothers(shifted, -1);
    if (current.hasBranched) temp.scale(current.Q2Now);
    else temp.scale(0.);
    event.append(temp);
  }

  // Also write current largest colour index in use.
  event.colIndx(maxColIndx);

}
     
//*********
  
// Set up primary partons for evolution.

bool TimeShower::setUpPrimary(double Q2maxin, long MEkindin, 
  long MEcombiin, double MEmixin) {

  // Check that number of partons not too high or energy not too low.
  if (nPrimary < 2 || nPrimary > 10) return false;
  double eKin = entry[0].m();
  double Q0NowSum = 0.;
  for (long i = 1; i < long(entry.size()); ++i) {
    TimeParticle& current = entry[i]; 
    eKin -= current.m();
    Q0NowSum += (iColour(current.id()) > 1) ? min(Q0, Q0CHGQ) : Q0CHGL;
  } 
  if (eKin < Q0NowSum) return false;  

  // If two initial partons, find ME correction kind.
  hasME = false;
  if (nPrimary == 2 && MEMODE >= 1) {
    findMEkind(MEkindin, MEcombiin, MEmixin); 
    if (MEkind > 0) hasME = true;
  }

  // Define sisters for pair.
  if (nPrimary == 2) {entry[1].sister = 2; entry[2].sister = 1;}
    
  // Check nature of initial partons, to see if they can shower.
  for (long i = 1; i < long(entry.size()); ++i) {
    TimeParticle& current = entry[i];     
    current.copyOfPrimary = i;
    bool isColoured = (iColour(current.id()) > 1) ? true : false;
    bool isCharged = (QEDSHOWER >=1 && iCharge(current.id()) != 0) ? 
      true : false;
    bool isPhoton = (QEDSHOWER >= 2 && current.id() == 22) ? true : false;

    // If yes, set maximum virtualities.
    if (current.status() > 0 && (isColoured || isCharged || isPhoton) ) {
      current.canBranch = true;
      double m = current.m();
      current.Q2Now = eKin * (eKin + 2.*m);
      if (Q2maxin > 0. && Q2maxin < current.Q2Now) current.Q2Now = Q2maxin; 
      if (RESPECTSCALE == 1 && current.scale() < current.Q2Now) 
        current.Q2Now = current.scale(); 
      current.eNow = current.e();
      current.eMax = eKin + m;

      // Decide which of colour and anticolour should set maximum emission cone.
      current.coneSide = 0;
      if (INITIALCONE >= 1 && (current.hasColVec4 || current.hasAntiColVec4)) {
        current.coneSide = 1;
        if (current.hasAntiColVec4 && (!current.hasColVec4 
          || Rndm::flat() > 0.5)) current.coneSide = 2;
      }

    // Mark partons that cannot branch. Finish loop.
    } else current.canBranch = false;
    current.hasBranched = false;
  }
  return true;
}

//*********

// Check and set kinematics for primary partons after evolution.

bool TimeShower::kinemPrimary() {

  // Check that sum of daughter masses not above mother mass.
  double eCM = entry[0].m();
  double eCMs = eCM * eCM;
  double eKin = eCM;
  double Q0NowSum = 0.;
  for (long i = 1; i <= nPrimary; ++i) {
    TimeParticle& current = entry[i]; 
    current.shouldEvolveMore = (current.hasBranched) ? true : false;
    eKin -= current.mNow();
    Q0NowSum += (iColour(current.id()) > 1) ? min(Q0, Q0CHGQ) : Q0CHGL;
  }
  if (eKin < Q0NowSum) return false;

  // Find particle energies (in cm frame) for two-particle configurations.
  if (nPrimary == 2) {
    double m1s = entry[1].m2Now(); 
    double m2s = entry[2].m2Now(); 
    entry[1].eNow = 0.5 * (eCMs + m1s - m2s) / eCM; 
    entry[2].eNow = 0.5 * (eCMs + m2s - m1s) / eCM; 
    entry[1].pAbsNow = entry[2].pAbsNow = pAbsLambdaK(eCMs, m1s, m2s); 

  // Find particle energies (in cm frame) for multi-particle configurations.
  } else {
    double eSum = 0.;
    double rSum = 0.;
    for (long i = 1; i <= nPrimary; ++i) {
      TimeParticle& current = entry[i]; 
      current.pAbsNow = current.pAbs();
      current.eNow = sqrt( current.p2() + current.m2Now() );
      eSum = eSum + current.eNow;
      rSum = rSum + current.m2Now() / current.eNow; 
    }
    long loop = 0;
    do {
      ++loop;
      double fac = (eCM - rSum) / (eSum - rSum);
      eSum = 0.;
      rSum = 0.;
      for (long i = 1; i <= nPrimary; ++i) {
        TimeParticle& current = entry[i]; 
        current.pAbsNow = fac * current.pAbsNow;
        current.eNow = sqrt( pow(current.pAbsNow,2) + current.m2Now() );
        eSum = eSum + current.eNow;
        rSum = rSum + current.m2Now() / current.eNow; 
      }
    } while ( loop < 10 && abs(eSum - eCM) > 1e-12 * eCM ); 
  }    

  // Check whether each branching corresponds to acceptable (z, Q2) pair.
  for (long i = 1; i <= nPrimary; ++i) entry[i].shouldEvolveMore = zQcheck(i);
  for (long i = 1; i <= nPrimary; ++i) {
    if (entry[i].shouldEvolveMore) return false;
  } 

  // When kinematics is fine, rescale momenta to correct energy.
  for (long i = 1; i <= nPrimary; ++i) {
    TimeParticle& current = entry[i]; 
    double fac = current.pAbsNow / current.pAbs();
    current.rescalep(fac);
    current.e(current.eNow);     
  }

  return true;
}

//*********

// Set up daughters of parton branching for subsequent evolution.
  
void TimeShower::setUpBranching(long i0) {

  // Add daughter entries to bottom of shower.
  entry.push_back(TimeParticle());
  long i1 = entry.size() - 1; 
  entry.push_back(TimeParticle());
  long i2 = entry.size() - 1; 
  TimeParticle& mother = entry[i0]; 
  TimeParticle& dau1 = entry[i1];
  TimeParticle& dau2 = entry[i2];

  // Define daughter status and properties.
  dau1.status(1);
  dau1.mother1(i0); 
  dau1.sister = i2;
  dau1.canBranch = true;
  dau1.hasBranched = false;
  dau1.copyOfPrimary = mother.copyOfPrimary;
  if (mother.id() == 21 || mother.id() == 22) dau1.copyOfPrimary = 0;
  dau1.coneSide = 0;
  dau2.status(1);
  dau2.mother1(i0); 
  dau2.sister = i1;
  dau2.canBranch = true;
  dau2.hasBranched = false;
  dau2.copyOfPrimary = 0;
  dau2.coneSide = 0;

  // Define daughter flavours, kind, masses, and colours. 
  // Gluon emission case.
  if (mother.idDaughter == 21) {
    dau1.status(mother.status());
    dau1.id(mother.id());
    dau2.id(21);
    dau1.m(mother.m());
    ++maxColIndx;
    if (mother.col() > 0 && mother.anticol() > 0) {
      if( Rndm::flat() < 0.5) {
        dau1.cols( mother.col(), maxColIndx );
        dau2.cols( maxColIndx, mother.anticol() );
      } else {
        dau1.cols( maxColIndx, mother.anticol() );
        dau2.cols( mother.col(), maxColIndx );
      }
    } else if (mother.col() > 0) {
      dau1.cols( maxColIndx, 0);
      dau2.cols( mother.col(), maxColIndx );
    } else { 
      dau1.cols( 0, maxColIndx);
      dau2.cols( maxColIndx, mother.anticol() );
    }

  // Photon emission case.
  } else if (mother.idDaughter == 22) {
    dau1.status(mother.status());
    dau1.id(mother.id());
    dau2.id(22);
    dau1.m(mother.m());
    dau1.cols( mother.col(), mother.anticol() );
    dau2.cols( 0, 0 );

  // Gluon or photon splitting case.
  } else {
    long idDau =  mother.idDaughter;
    dau1.id( idDau);
    dau2.id(-idDau);
    double mDau = Mass(idDau);
    dau1.m(mDau);
    dau2.m(mDau);
    if (mother.id() == 21) {
      dau1.cols( mother.col(), 0 );
      dau2.cols( 0, mother.anticol() );
    } else if (idDau < 10) {
      ++maxColIndx;
      dau1.cols( maxColIndx, 0 );
      dau2.cols( 0, maxColIndx );
    } else {
      dau1.cols( 0, 0 );
      dau2.cols( 0, 0 );
    }
    dau1.copyOfPrimary = i1;
    dau2.copyOfPrimary = i2;
  }

  // Construct kinematics assuming daughters on mass shell.
  kinemConstruct( i0, i1, i2);

  // Set maximal virtualities and energies for daughter evolution.
  double Q0NowSum = (iColour(dau1.id()) > 1) ? min(Q0, Q0CHGQ) : Q0CHGL;
  Q0NowSum += (iColour(dau2.id()) > 1) ? min(Q0, Q0CHGQ) : Q0CHGL;
  double m0s = pow(mother.mNow() - Q0FRACPS * Q0NowSum, 2);
  dau1.Q2Now = min(m0s, pow(dau1.e(),2)) - dau1.m2();
  dau2.Q2Now = min(m0s, pow(dau2.e(),2)) - dau2.m2();
  double z = mother.zDaughter;
  dau1.eMax = z * mother.e();  
  dau2.eMax = (1.-z) * mother.e();
  dau1.scale(mother.scale());
  dau2.scale(mother.scale());
 
  // Update status code of mother that decayed.
  mother.addstatus(3);

}

//*********

// Check and set kinematics for branching after daughter evolution;
// especially pick non-isotropic azimuthal angle.

bool TimeShower::kinemBranching(long i0, long i1, long i2) {

  // Check that sum of daughter masses not above mother mass.
  TimeParticle& mother = entry[i0];
  TimeParticle& dau1 = entry[i1];
  TimeParticle& dau2 = entry[i2];
  dau1.shouldEvolveMore = (dau1.hasBranched) ? true : false;
  dau2.shouldEvolveMore = (dau2.hasBranched) ? true : false;
  double Q0NowSum = (iColour(dau1.id()) > 1) ? min(Q0, Q0CHGQ) : Q0CHGL;
  Q0NowSum += (iColour(dau2.id() > 1)) ? min(Q0, Q0CHGQ) : Q0CHGL;
  if ( (dau1.hasBranched || dau2.hasBranched) 
    && mother.mNow() - dau1.mNow() - dau2.mNow() < Q0FRACPS * Q0NowSum) 
    return false;

  // Construct kinematics of branching.
  kinemConstruct(i0, i1, i2);
  long loop = 0;

  // Check whether each branching corresponds to acceptable (z, Q2) pair.
  dau1.shouldEvolveMore = zQcheck(i1);
  dau2.shouldEvolveMore = zQcheck(i2);
  if (dau1.shouldEvolveMore || dau2.shouldEvolveMore) return false; 

  // Grandmother; softer of two gluons in g -> g g, else only gluon.
  TimeParticle& grandmother = entry[mother.mother1()];    
  long iSoftG = i2;
  if (mother.id() == 21 && dau1.e() < dau2.e()) iSoftG = i1; 
  TimeParticle& softG = entry[iSoftG];  

  // Azimuthal asymmetry coefficient from interference with initial colours.
  bool haveConeAsym = false;
  double asymCone = 0.;
  Vec4 coneP;
  if (INITIALCONE >= 2 && mother.coneSide >=1 && mother.idDaughter == 21) {
    coneP = (mother.coneSide == 1) ? mother.colVec4 : mother.antiColVec4; 
    double thetaGM = theta(softG.p(), mother.p());
    double thetaCM = theta(coneP, mother.p());
    asymCone = min( 0.95, (thetaGM/thetaCM) * (1. - thetaCM/M_PI));
    if (asymCone > 0.0001) haveConeAsym = true;
  }

  // Azimuthal asymmetry coefficient from gluon/photon polarization.
  bool havePolAsym = false;
  double asymPol = 0.;
  if (PHIPOLASYM == 1 && (mother.id() == 21 || mother.id() == 22) &&
    i0 > nPrimary) {
    double zAunt = grandmother.zDaughter;     
    if (mother.sister > i0) zAunt = 1. - zAunt; 
    asymPol = (grandmother.id() != 21 && grandmother.id() != 22) 
      ? 2. * zAunt / (1. + zAunt*zAunt) 
      : pow( zAunt / (1. - zAunt*(1.-zAunt)), 2);
    double z = mother.zDaughter;
    asymPol *= (dau1.id() != 21) ? - z * (1.-z) / (1. - 2.*z*(1.-z)) 
      : pow( z * (1.-z) / (1. - z*(1.-z)), 2);  
    if (asymPol > 0.0001) havePolAsym = true;
  }

  // Azimuthal asymmetry coefficient from soft gluon coherence.
  bool haveCoherAsym = false;
  double asymCoher = 0.;
  if (PHICOHERASYM == 1 && mother.idDaughter == 21 && i0 > nPrimary) { 
    double zSoft = (iSoftG == i1) ? mother.zDaughter : 1. - mother.zDaughter;
    double zMother = grandmother.zDaughter;
    if (mother.sister < i0) zMother = 1. - zMother;
    asymCoher = min(0.95,  (mother.mNow()/grandmother.mNow()) 
      * sqrt( (1.-zSoft) * (1.-zMother) / (zSoft * zMother) ));     
    if (asymCoher > 0.0001) haveCoherAsym = true;
    if (entry[mother.sister].id() == 22) haveCoherAsym = false;
  }

  // Reweight relative phi angle by rejection and renewed kinematics.
  if (haveConeAsym || havePolAsym || haveCoherAsym) {
    double wtPhi;
    do {
      if (loop > 0) kinemConstruct(i0, i1, i2);
      ++loop;
      wtPhi = 1.;

      // Azimuthal asymmetry weight from interference with initial colours.
      if (haveConeAsym) {
        double cosPhiRel = cosphi(softG.p(), coneP, mother.p());
        wtPhi *= (1. - asymCone) * (1. - asymCone * cosPhiRel) 
          / (1. + asymCone*asymCone - 2. * asymCone * cosPhiRel);      
      }

      // Azimuthal asymmetry weight from gluon/photon polarization.
      if (havePolAsym) {
        double cosPhiRel = cosphi(dau1.p(), grandmother.p(), mother.p());     
        wtPhi *= (1. + asymPol * (2. * cosPhiRel*cosPhiRel - 1.))
	  / (1. + abs(asymPol)); 
      }

      // Azimuthal asymmetry weight from soft gluon coherence.
      if (haveCoherAsym) {
        double cosPhiRel = cosphi(softG.p(), grandmother.p(), mother.p());
        wtPhi *= (1. - asymCoher) * (1. - asymCoher * cosPhiRel) 
          / (1. + asymCoher*asymCoher - 2. * asymCoher * cosPhiRel);      
      }
  
    // Appropriate rejection of phi angles.
    } while (wtPhi < Rndm::flat());   
  }

  // Kinematics completed and accepted.  
  return true;
}

//*********

// Evolve a single parton; main physics routine.

void TimeShower::evolveParton(long i) {

  // Reset flag. 
  TimeParticle& evolving = entry[i];
  evolving.hasBranched = false;

  // Check whether QCD and QED emissions allowed, separately.
  long intColour = iColour(evolving.id());
  bool isColoured = (intColour > 1) ? true : false;
  bool isGluon = (evolving.id() == 21) ? true : false;
  long intCharge = iCharge(evolving.id());
  double charge = intCharge/3.;  
  bool isCharged = (QEDSHOWER >=1 && intCharge != 0) ? true : false;
  bool isPhoton = (QEDSHOWER >= 2 && evolving.id() == 22) ? true : false;

  // Set Q0 scale based on the allowed branchings
  double Q0Now;
  if (isColoured && isCharged) Q0Now = min(Q0, Q0CHGQ);
  else if (isColoured) Q0Now = Q0;
  else if (isCharged || isPhoton) Q0Now = Q0CHGL;
  else return; 
  
  // Find upper estimate of allowed z range.
  double m = evolving.m();
  double m2 = m*m;
  double Q2 = evolving.Q2Now; 
  if (Q2 < Q0Now*(Q0Now+m)) return; 
  double eMax = evolving.eMax;
  if (eMax < 1.01 * Q0Now) return; 
  double zCut = 0.5 * (1. - sqrt(1. - pow(Q0Now/eMax, 2)));
  if (zCut < 1e-7) zCut = pow(0.5*Q0Now/eMax, 2);
  
  // Integral of QCD AP kernel: (s)quark, gluino, gluon.
  // (Rate increased also for (s)quark if sister is gluino!) 
  double evolCol = 0.;
  double LAMBDA2 = LAMBDA5 * LAMBDA5; 
  double alphaLogMin = log(0.25 * Q0Now*Q0Now / LAMBDA2);
  double gqqFrac = 0.;
  double nQuark5 = min(5., double(NQUARK));
  if (isColoured) {
    double kernelInt = 0.;
    if (!isGluon) {
      kernelInt = (8./3.) * log((1.-zCut)/zCut);
      long iCopy = evolving.copyOfPrimary;
      if (intColour == 8 || (hasME && MEgluinoDau && iCopy > 0 
        && iCopy < nPrimary) ) kernelInt *= 9./4.;
    } else {
      kernelInt = 6. * log((1.-zCut)/zCut) + 0.5 * nQuark5;
      gqqFrac = 0.5 * nQuark5 / kernelInt;
    } 
    // Evolution coefficient depends on alpha_s running scheme.
    if (ALPHASMODE <= 0) evolCol = 2. * M_PI / (ALPHASFIX * kernelInt);
    else if (ALPHASMODE == 1) evolCol = 23. / (6. * kernelInt); 
    else evolCol = 23. * alphaLogMin / (6. * kernelInt);
  }

  // Integral of QED AP kernel.
  double evolChg = 0.; 
  if (isCharged) {
    double kernelInt = charge*charge * 2. * log((1.-zCut)/zCut);
    evolChg = 2. * M_PI / (ALPHAEMFIX * kernelInt); 
  } else if (isPhoton) {
    double kernelInt = (20./3.); 
    evolChg = 2. * M_PI / (ALPHAEMFIX * kernelInt); 
  }

  // Optional: set up maximal cone of emission for primaries.
  double thetaInitMax = 100.;
  if (evolving.coneSide == 1) thetaInitMax = 
    theta(evolving.p(), evolving.colVec4);
  else if (evolving.coneSide == 2) thetaInitMax = 
    theta(evolving.p(), evolving.antiColVec4);

  // Begin loop over evolution in Q2.
  double z, wtZQ=0.0, zRange; 
  long idDau;
  bool isOKbranch;
  do {

    // First pick a Q2 for QCD shower, if relevant.
    isOKbranch = false;
    double Q2col = Q2;
    if (isColoured) {
      // Evolution with fixed alpha_s (or upper estimate for pT^2 scale).
      if (ALPHASMODE != 1) Q2col *= pow(Rndm::flat(), evolCol); 
      // Evolution with alpha_s(m^2/4).
      else Q2col = 4. * LAMBDA2 * pow( 0.25 * Q2col / LAMBDA2, 
        pow(Rndm::flat(), evolCol));
    }

    // Then pick a Q2 for QED shower, if relevant.
    double Q2chg = Q2;
    if (isCharged || isPhoton) {
      // Evolution with fixed alpha_em.
      //      Q2chg *= pow(Rndm::flat(), evolChg); 
      Q2chg *= exp(max(-100.0, evolChg*log(Rndm::flat())));
    } 

    // Pick between QCD and QED branching where required.
    // May need to reassign if one evolution below cut while other above.
    bool pickedQCD = (isColoured) ? true : false;
    if (isColoured && isCharged) {
      if (Q2chg > Q2col) pickedQCD = false;
      if (pickedQCD && Q0 > Q0CHGQ && Q2col < Q0*(Q0+m)) pickedQCD = false;
      if (!pickedQCD && Q0CHGQ > Q0 && Q2chg < Q0CHGQ*(Q0CHGQ+m)) 
        pickedQCD = true;
    } 
    Q2 = (pickedQCD) ? Q2col : Q2chg;  

    // Evolution finished if below cut-off. 
    if (pickedQCD && Q2 < Q0*(Q0+m)) return;
    if (!pickedQCD && isColoured && Q2 < Q0CHGQ*(Q0CHGQ+m)) return;  
    if (!pickedQCD && !isColoured &&  Q2 < Q0CHGL*(Q0CHGL+m)) return;

    // Select z value of branching quark/squark/gluino -> ditto + g.
    if (pickedQCD && !isGluon) {
      idDau = 21; 
      z = 1. - (1.-zCut) * pow(zCut/(1.-zCut), Rndm::flat()) ;
      // Only use splitting kernel when no matrix element later on.
      if (hasME && evolving.copyOfPrimary > 0) wtZQ = 1.;
      else wtZQ = 0.5 * (1. + z*z); 
    
    // Select z value of branching g -> g + g.
    } else if (pickedQCD && gqqFrac < Rndm::flat()) {
      idDau = 21; 
      z = (1.-zCut) * pow(zCut/(1.-zCut), Rndm::flat()) ;
      if (Rndm::flat() > 0.5) z = 1. - z;
      wtZQ = pow(1. - z * (1.-z), 2);
     
    // Select flavour and z value of branching g -> q + qbar.
    // Use definition z = (1 + cos(theta*))/2 here, for ME match.
    } else if (pickedQCD) {
      idDau = min(5L, 1 + long(floor(nQuark5 *Rndm::flat()))); 
      double mDau2 = pow( Mass(idDau), 2);
      double mRat2 = 4. * mDau2 / (m2 + Q2);
      if (mRat2 > 0.999999) continue; 
      z = Rndm::flat();       
      wtZQ = ( z*z + (1.-z)*(1.-z) + 2. * mRat2 * z * (1.-z) )
        * sqrt(1. - mRat2); 

    // Select z value of charged -> charged + gamma. 
    } else if (!isPhoton) {
      idDau = 22;
      z = 1. - (1.-zCut) * pow(zCut/(1.-zCut), Rndm::flat()) ;
      // Not any matrix element correction here (yet).
      wtZQ = 0.5 * (1. + z*z); 

    } else {
    // Select flavour and z value of branching gamma -> q + qbar, l + lbar.
    // Use definition z = (1 + cos(theta*))/2 here, for ME match.
      long idTmp = long(floor(20. *Rndm::flat()));
      if (idTmp < 1) idDau = 1;
      else if (idTmp < 5) idDau = 2;
      else if (idTmp < 6) idDau = 3;
      else if (idTmp < 10) idDau = 4;
      else if (idTmp < 11) idDau = 5;
      else if (idTmp < 14) idDau = 11;
      else if (idTmp < 17) idDau = 13;
      else idDau = 15;
      if (idDau < 10 && idDau > NQUARK) continue; 
      if (idDau < 10 && Q2 < Q0CHGQ*Q0CHGQ) continue; 
      if (idDau > 10 && Q2 < Q0CHGL*Q0CHGL) continue; 
      double mDau2 = pow( Mass(idDau), 2);
      double mRat2 = 4. * mDau2 / (m2 + Q2);
      if (mRat2 > 0.999999) continue; 
      z = Rndm::flat();       
      wtZQ = ( z*z + (1.-z)*(1.-z) + 2. * mRat2 * z * (1.-z) )
        * sqrt(1. - mRat2); 
    }

    // Save temporary values (used by ME correction). 
    evolving.Q2Now = Q2;
    evolving.idDaughter = idDau;
    evolving.zDaughter = z;
    
    // Correct to alpha_strong of pT^2, alternatively of Q^2/4 for g -> q qbar.
    if (pickedQCD && ALPHASMODE >= 2) {
      if (ALPHASMODE >= 3 && evolving.idDaughter <= 5) {
        if (alphaLogMin / log(0.25 * Q2 / LAMBDA2) < Rndm::flat()) continue;
      } else {
        double pT2 = z * (1. - z) * Q2;
        if (pT2 < 0.25 * Q0*Q0) continue;
        if (alphaLogMin / log(pT2 / LAMBDA2) < Rndm::flat()) continue;
      }
    }

    // Check that z consistent with chosen mass (for preliminary energy).
    // No check for g/gamma -> q + qbar, l + lbar,  where angular z works.
    if (idDau == 21 || idDau == 22) {
      double eHere = eMax;
      if (i <= nPrimary && nPrimary == 2) eHere = 0.5 * (entry[0].m2() 
        + m2 + Q2 - entry[3-i].m2()) / entry[0].m();
      zRange = sqrtpos(1. - (m2 + Q2) / (eHere * eHere));
      if ( z <= 0.5 * (1. - zRange) || z >= 0.5 * (1. + zRange) ) continue; 
    }

    // Perform QCD matrix element correction where allowed. 
    if (pickedQCD && hasME && evolving.copyOfPrimary > 0 && !isGluon) 
      wtZQ *= findMEcorr(i);

    // Impose angular ordering constraint on QCD emissions (but not QED ones).
    if (pickedQCD && i > nPrimary && ( (ANGULARORDER == 1 && idDau == 21) 
      || ANGULARORDER >= 2)) {
      long mo = evolving.mother1();
      while (entry[mo].idDaughter == 22 && mo > nPrimary) 
        mo = entry[mo].mother1(); 
      double zProd = entry[mo].zDaughter;
      double theta2Prod = entry[mo].Q2Now / (zProd * (1.-zProd) * 
        pow(entry[mo].e(), 2));
      double theta2Dec = Q2 / (z * (1.-z) * pow(evolving.eNow, 2));
      if (entry[mo].idDaughter != 22 && theta2Dec > theta2Prod ) continue; 
    } 

    // Restrict maximal cone of emission for primaries (optional).
    if (evolving.coneSide >= 1) {
      double theta2Max = max( z/(1.-z), (1.-z)/z) * Q2 / pow(evolving.eNow, 2); 
      if (theta2Max > thetaInitMax*thetaInitMax) continue;  
    }

    // Respect user-set maximum scale for emissions - different alternatives.
    if (RESPECTSCALE == 2 && z * (1. - z) * Q2 > evolving.scale()) wtZQ = 0.;
    if (RESPECTSCALE == 3 && Q2 / (z * (1.-z)) > evolving.scale()) wtZQ = 0.;
    if (RESPECTSCALE == 4 && Q2 / (z * (1.-z) * pow(evolving.eNow, 2)) 
      > evolving.scale()) wtZQ = 0.;

    // Monte Carlo rejection; continue evolution in Q2.
    isOKbranch = true;
  } while (!isOKbranch || wtZQ < Rndm::flat() );
  evolving.hasBranched = true;
}

//*********

// Pick one of the partons for further evolution when required.
  
long TimeShower::pickParton(long iMin, long iMax) {

  // The simple choices when only one parton can be picked.
  long picked = 0;
  long nPossible = 0;
  for (long i = iMin; i <= iMax; ++i) {
    if (entry[i].shouldEvolveMore) {     
      picked = i;
      ++nPossible;
    }
  } 
  if (nPossible == 1) return picked;

  // Else pick one at random from possible ones.
  double randPicked = nPossible * Rndm::flat();
  for (long i = iMin; i <= iMax; ++i) {
    if (entry[i].shouldEvolveMore) { 
      --randPicked;
      if (randPicked <= 0.) return i;
    }
  }    
  return -1;
}

//*********

// Check whether z and Q2 choices of branching are inconsistent.
// Gives value to shouldEvolveMore, i.e. true = inconsistent.

bool TimeShower::zQcheck(long i) {

  // Check (z, Q2) pair for gluon or photon emission.
  TimeParticle& current = entry[i];
  if (current.hasBranched && (current.idDaughter == 21 
    || current.idDaughter == 22)) {
    double z = current.zDaughter; 
    double zRange = sqrtpos(1. - current.m2Now() / pow(current.eNow, 2));
    if (z < 0.5 * (1. - zRange) || z > 0.5 * (1. + zRange)) return true;
  }
  // But gluon of photon splitting is always OK (at this stage).
  return false;

}

//*********

// Construct branching kinematics for setUpBranching & kinemBranching.

void TimeShower::kinemConstruct(long i0, long i1, long i2) {

  // Info on decaying mother and decay daughters.
  TimeParticle& mother = entry[i0];
  TimeParticle& dau1 = entry[i1];
  TimeParticle& dau2 = entry[i2];
  double m0s = mother.m2Now();
  double e0 = mother.e();
  double p0 = mother.pAbs();
  double m1s = dau1.m2Now(); 
  double m2s = dau2.m2Now();

  // Begin kinematics construction.
  double lambda = lambdaKRoot(m0s, m1s, m2s);
  double pT, e1, pz1;
  double z = mother.zDaughter;

  // For gluon or photon emission z is energy fraction.
  if (mother.idDaughter == 21 || mother.idDaughter == 22) {
    pT = lambda * sqrtpos(z * (1.-z) * e0*e0 / m0s - 0.25) / p0;
    e1 = e0 * (0.5 * (m0s - lambda + m1s - m2s) + lambda * z) / m0s;
    pz1 = (e0 * e1 - 0.5 * (m0s + m1s - m2s)) / p0; 

  // For gluon or photon splitting z = (1 + cos(theta*))/2.
  } else {
    pT = lambda * sqrt( z * (1.-z) / m0s);
    double epz1 = (e0 + p0) * (0.5 * (m0s - lambda + m1s - m2s) 
      + lambda * z) / m0s;
    e1 = 0.5 * ( epz1 + (m1s + pT*pT) / epz1); 
    pz1 = 0.5 * ( epz1 - (m1s + pT*pT) / epz1);
  }

  // Pick azimuthal phi angle evenly - modified by rejection in kinemBranching.
  double phi = 2. * M_PI * Rndm::flat();
  mother.phiDaughter = phi;

  // Fill info on momenta, including rotation to proper direction.
  dau1.p( pT * cos(phi), pT * sin(phi), pz1 , e1) ;
  dau1.eNow = e1; 
  dau2.p( -pT * cos(phi), -pT * sin(phi), p0 - pz1, e0 - e1) ;
  dau2.eNow = e0 - e1;
  double theta0 = mother.theta();
  double phi0 = mother.phi();
  dau1.rot(theta0, phi0);
  dau2.rot(theta0, phi0);

}

//*********

// Find class of QCD ME correction.

void TimeShower::findMEkind(long MEkindin, long MEcombiin, 
  double MEmixin) {
  
  // Find type of mother and daughters.
  long idMother = entry[0].id();
  long motherType = findMEparticle(idMother);
  long dau1Type = findMEparticle(entry[1].id());
  long dau2Type = findMEparticle(entry[2].id());
  long minDauType = min(dau1Type, dau2Type);
  long maxDauType = max(dau1Type, dau2Type);
  MEorder = (dau2Type > dau1Type) ? true : false;
  MEsplit = (maxDauType <= 2 || maxDauType == 6) ? true : false; 
  MEgluinoDau = (maxDauType == 6) ? true : false;
 
  // If nonvanishing input then respect that.
  if (MEkindin > 0) {
    MEkind = MEkindin;
    MEcombi = MEcombiin;
    MEmix = MEmixin;
    return;
  }
        
  // Else start from default, which is no ME corrections, 
  // and try to find matching ME cases below.
  MEkind = 0;
  MEcombi = 0;
  MEmix = 0.5;
  if (minDauType == 0) return;

  // Vector/axial vector -> q + qbar; q -> q + V.
  if (minDauType == 1 && maxDauType == 1 && 
    (motherType == 3 || motherType == 0) ) {
    MEkind = 2;
    if (idMother == 21 || idMother == 22) MEcombi = 1;
    else if (idMother == 23) {MEcombi = 3; MEmix = gammaZmix();}
    else if (idMother == 24 || idMother == 0) MEcombi = 4;
  }
  // For chi -> chi q qbar, use V/A -> q qbar as first approximation.
  else if (minDauType == 1 && maxDauType == 1 && motherType == 5)
    MEkind =2;
  else if (minDauType == 1 && maxDauType == 3 && (motherType == 0
    || motherType == 1)) MEkind = 3;
 
  // Scalar/pseudoscalar -> q + qbar; q -> q + S.
  else if (minDauType == 1 && maxDauType == 1 && motherType == 4) {
    MEkind =4;
    if (idMother == 25 || idMother == 35 || idMother == 37) MEcombi = 1;
    else if (idMother == 36) MEcombi = 2;
  } 
  else if (minDauType == 1 && maxDauType == 4 && 
    (motherType == 0 || motherType == 1) ) MEkind = 5;
 
  // V -> ~q + ~qbar; ~q -> ~q + V; S -> ~q + ~qbar; ~q -> ~q + S.
  else if (minDauType == 2 && maxDauType == 2 && 
    (motherType == 0 || motherType == 3) ) MEkind = 6;
  else if (minDauType == 2 && maxDauType == 3 && 
    (motherType == 0 || motherType == 2) ) MEkind = 7;
  else if (minDauType == 2 && maxDauType == 2 && motherType == 4)
    MEkind = 8;
  else if (minDauType == 2 && maxDauType == 4 && 
    (motherType == 0 || motherType == 2) ) MEkind = 9;
 
  // chi -> q + ~qbar; ~q -> q + chi; q -> ~q + chi.
  else if (minDauType == 1 && maxDauType == 2 && 
    (motherType == 0 || motherType == 5) ) MEkind = 10;
  else if (minDauType == 1 && maxDauType == 5 && 
    (motherType == 0 || motherType == 2) ) MEkind = 11;
  else if (minDauType == 2 && maxDauType == 5 && 
    (motherType == 0 || motherType == 1) ) MEkind = 12;
 
  // ~g -> q + ~qbar; ~q -> q + ~g; q -> ~q + ~g.
  else if (minDauType == 1 && maxDauType == 2 && motherType == 6)
    MEkind = 13;
  else if (minDauType == 1 && maxDauType == 6 && 
    (motherType == 0 || motherType == 2) ) MEkind = 14;
  else if (minDauType == 2 && maxDauType == 6 && 
    (motherType == 0 || motherType == 1) ) MEkind = 15;

  // g -> ~g + ~g (eikonal approximation).
  else if (minDauType == 6 && maxDauType == 6 && motherType == 0)
    MEkind = 16;

}

//*********
 
// Find type of particle for ME type: 0 = unknown, 1 = quark,
// 2 = squark, 3 = vector boson (also g), 4 = colourless scalar,
// 5 = colourless neutralino/chargino, 6 = gluino.

long TimeShower::findMEparticle(long id) {

  // find colour and spin of particle.
  long type = 0;
  long col = iColour(id); 
  long spin = iSpin(id);

  // Find particle type from colour and spin.
       if (col == 3 && spin == 2) type = 1;
  else if (col == 3 && spin == 1) type = 2;
  else if (col == 1 && spin == 3) type = 3;
  else if (col == 8 && spin == 3) type = 3;
  else if (col == 1 && spin == 1) type = 4;
  else if (col == 1 && spin == 2) type = 5;
  else if (col == 8 && spin == 2) type = 6;
  return type;
}  

//*********

// Find mixture of V and A in gamma/Z: energy- and flavour-dependent. 

double TimeShower::gammaZmix() {

  // Some parameters that eventually will be fed in externally.
  double SIN2W = 0.232;
  double MZ = 91.188;
  double GAMMAZ = 2.495;
  double mix = 0.5;

  // Initial flavours and couplings; return if don't make sense.
  if (inFlavour1 + inFlavour2 != 0) return mix;
  long inFlav = abs(inFlavour1);
  double eIn;
  if (inFlav < 10 && inFlav%2 == 1) eIn = -1./3.;
  else if (inFlav < 10) eIn = 2./3.;
  else if (inFlav < 20 && inFlav%2 == 1) eIn = -1.;
  else if (inFlav < 20) eIn = 0.;
  else return mix;
  double aIn = (eIn < -0.1) ? -1. : 1.;
  double vIn = aIn - 4. * SIN2W * eIn; 

  // Final flavours and couplings; return if don't make sense.
  if (entry[1].id() + entry[2].id() != 0) return mix;
  long outFlav = abs(entry[1].id());
  double eOut;
  if (outFlav < 10 && outFlav%2 == 1) eOut = -1./3.;
  else if (outFlav < 10) eOut = 2./3.;
  else if (outFlav < 20 && outFlav%2 == 1) eOut = -1.;
  else if (outFlav < 20) eOut = 0.;
  else return mix;
  double aOut = (eOut < -0.1) ? -1. : 1.;
  double vOut = aOut - 4. * SIN2W * eOut; 

  // Calculate vector and axial expressions and find mix.
  if (inFlav > 0) {
    double xwc = 1. / (16. * SIN2W * (1.-SIN2W));
    double sHat = entry[0].m2();
    double prop = 1. / ( pow(sHat - MZ*MZ, 2) + sHat * GAMMAZ*GAMMAZ); 
    double vect = eIn*eIn * eOut*eOut 
    + 2. * eIn*vIn * eOut*vOut * xwc * sHat * (sHat - MZ*MZ) * prop 
    + (vIn*vIn + aIn*aIn) * vOut*vOut * xwc*xwc * sHat*sHat * prop;
    double axiv = (vIn*vIn + aIn*aIn) * aOut*aOut * xwc*xwc * sHat*sHat * prop;
    mix = vect / (vect + axiv);
  } else { 
    mix = vOut*vOut / ( vOut*vOut + aOut*aOut);
  } 

  return mix;
}

//*********

// Set up to calculate QCD ME correction with calcMEcorr.
// Normally for primary particles, but also from g/gamma -> f fbar.
  
double TimeShower::findMEcorr(long i) {
  
  // Some initial values.
  double wtME = 1.;
  double wtPS = 1.; 
  double m1 = entry[i].m();
  double mNow1 = entry[i].mOff();
  double z = entry[i].zDaughter;
  long iCopy = entry[i].copyOfPrimary;
  double m2, eCM;
  
  // For primary partons the cm energy is well defined.
  if (i <= nPrimary) {
    m2 = entry[3-i].m();
    eCM = entry[0].m();
  // Else define effective reduced energy of remaining system.
  } else {
    m2 = (iCopy <= nPrimary) ? entry[iCopy].m() : m1;
    long iMother = entry[i].mother1();
    double mMother = entry[iMother].mNow();
    double eMother = entry[iMother].e();
    double zMother = entry[iMother].zDaughter;
    double eDau1 = max( mNow1, eMother * ( zMother + (1.-zMother) 
      * pow(mNow1/mMother,2) ) );
    eCM = eDau1 + sqrtpos(eDau1 * eDau1 - mNow1 * mNow1 + m2 * m2);
  }

  // Construct ME variables.
  double r1 = m1 / eCM;
  double r2 = m2 / eCM; 
  double x1 = (1. + pow(mNow1/eCM, 2) - r2*r2) 
    * (z + (1.-z) * pow(m1/mNow1, 2) );   
  double x2 = 1. + r2*r2 - pow(mNow1/eCM, 2);

  // Evaluate normal ME, for proper order of particles.
  if ( (i == 1 && MEorder) || (i == 2 && !MEorder) ) {
    wtME = calcMEcorr(MEkind, MEcombi, MEmix, x1, x2, r1, r2);
  } else if (i <= nPrimary) {
    wtME = calcMEcorr(MEkind, MEcombi, MEmix, x2, x1, r2, r1);
  // Evaluate ME for secondary g or gamma -> f fbar (f = fermion).
  // Note temporary choice for gluon !!
  } else if (entry[entry[iCopy].mother1()].id() == 21) {
    wtME = calcMEcorr(13, 1, 0., x2, x1, r2, r1);
  } else {
    wtME = calcMEcorr(2, 1, 0., x2, x1, r2, r1);
  }

  // Split up total ME when two radiating particles.
  if (MEsplit) wtME = wtME * max(1e-10, 1. + r1*r1 - r2*r2 - x1) 
    / max(1e-10, 2. - x1 - x2);

  // Evaluate shower rate to be compared with.
  wtPS = 2. / ( max(1e-10, 2.-x1-x2) * max(1e-10, 1.+r2*r2-r1*r1-x2) );
  if (MEgluinoDau) wtPS *= 9./4.;
  
  // Return ratio of actual ME to assumed PS rate of emission.
  return wtME/wtPS; 
}

//*********

// Matrix elements for gluon (or photon) emission from
// a two-body state; to be used by the parton shower routine.
// Here x_i = 2 E_i/E_cm, r_i = m_i/E_cm and
// 1/sigma_0 d(sigma)/d(x_1)d(x_2) = (alpha-strong/2 pi) * CF * PYMAEL,
// i.e. normalization is such that one recovers the familiar
// (X1**2+X2**2)/((1-X1)*(1-X2)) for the massless case.
// Coupling structure:
// kind =  1 : eikonal soft-gluon expression (spin-independent)
//      =  2 : V -> q qbar (V = vector/axial vector colour singlet)
//      =  3 : q -> q V
//      =  4 : S -> q qbar (S = scalar/pseudoscalar colour singlet)
//      =  5 : q -> q S
//      =  6 : V -> ~q ~qbar (~q = squark)
//      =  7 : ~q -> ~q V
//      =  8 : S -> ~q ~qbar
//      =  9 : ~q -> ~q S
//      = 10 : chi -> q ~qbar (chi = neutralino/chargino)
//      = 11 : ~q -> q chi
//      = 12 : q -> ~q chi
//      = 13 : ~g -> q ~qbar
//      = 14 : ~q -> q ~g
//      = 15 : q -> ~q ~g
//      = 16 : (9/4)*(eikonal) for gg -> ~g ~g
// Note that the order of the decay products is important.
// combi = 1 : pure non-gamma5, i.e. vector/scalar/...
//       = 2 : pure gamma5, i.e. axial vector/pseudoscalar/....
//       = 3 : mixture mix*(combi=1) + (1-mix)*(combi=2)
//       = 4 : mixture (combi=1) +- (combi=2)

double TimeShower::calcMEcorr(long kind, long combiin, double mixin, 
  double x1, double x2, double r1, double r2) {

  // Frequent variable combinations.
  double x3 = 2. - x1 - x2;
  double x1s = x1 * x1;
  double x2s = x2 * x2;
  double x3s = x3 * x3;
  double x1c = x1 * x1s;
  double x2c = x2 * x2s;
  double x3c = x3 * x3s;
  double r1s = r1 * r1;
  double r2s = r2 * r2;
  double r1c = r1 * r1s;
  double r2c = r2 * r2s;
  double r1q = r1s * r1s;
  double r2q = r2s * r2s;
  double prop1 = 1. + r1s - r2s - x1; 
  double prop2 = 1. + r2s - r1s - x2;
  double prop1s = prop1 * prop1;
  double prop2s = prop2 * prop2;
  double prop12 = prop1 * prop2;
  double prop13 = prop1 * x3;
  double prop23 = prop2 * x3;

  // Check input values. Return zero outside allowed phase space.
  if (x1 - 2.*r1 < 1e-10 || prop1 < 1e-10) return 0.;
  if (x2 - 2.*r2 < 1e-10 || prop2 < 1e-10) return 0.;
  if (x1 + x2 - 1. - pow(r1+r2,2) < 1e-10) return 0.;
  if ((x1s - 4.*r1s) * (x2s - 4.*r2s) - pow(2.*(1. - x1 - x2 + r1s + r2s) 
    + x1*x2, 2) < 1e-10) return 0.;

  // Initial values; phase space.
  long combi = combiin; 
  if (combiin < 1) combi = 1;
  if (combiin > 4) combi = 4;   
  double mix = max(0., min(1., mixin));
  bool isSet1 = false;
  bool isSet2 = false;
  bool isSet4 = false;
  double ps = lambdaKRoot(1., r1*r1, r2*r2);
  double rLO=0.0, rFO=0.0, rLO1=0.0, rFO1=0.0,
    rLO2=0.0, rFO2=0.0, rLO4=0.0, rFO4=0.0;
  double offset = 0;
 
  // Select which kind of ME to use.
  switch (kind) {

    // case 1 is equal to default, i.e. eikonal expression.

    // V -> q qbar (V = gamma*/Z0/W+-/...).
    case 2:
      if (combi == 1 || combi == 3) {
        rLO1 = ps*(2.-r1s-r1q+6.*r1*r2-r2s+2.*r1s*r2s-r2q)/2.;
        rFO1 = -(3.+6.*r1s+r1q-6.*r1*r2+6.*r1c*r2-2.*r2s-6.*r1s*r2s
        +6.*r1*r2c+r2q-3.*x1+6.*r1*r2*x1+2.*r2s*x1+x1s-2.*r1s*x1s
        +3.*r1s*x3+6.*r1*r2*x3-r2s*x3-2.*x1*x3-5.*r1s*x1*x3
        +r2s*x1*x3+x1s*x3-3.*x3s-3.*r1s*x3s+r2s*x3s
        +2.*x1*x3s+x3c-x2)
        /prop2s
        -2.*(-3.+r1s-6.*r1*r2+6.*r1c*r2+3.*r2s-4.*r1s*r2s
        +6.*r1*r2c+2.*x1+3.*r1s*x1+r2s*x1-x1s-r1s*x1s
        -r2s*x1s+4.*x3+2.*r1s*x3+3.*r1*r2*x3-r2s*x3-3.*x1*x3
        -2.*r1s*x1*x3+x1s*x3-x3s-r1s*x3s+r1*r2*x3s+x1*x3s)
        /prop12
        -(-1.+2.*r1s+r1q+6.*r1*r2+6.*r1c*r2-2.*r2s-6.*r1s*r2s
        +6.*r1*r2c+r2q-x1-2.*r1s*x1-6.*r1*r2*x1+8.*r2s*x1+x1s
        -2.*r2s*x1s-r1s*x3+r2s*x3-r1s*x1*x3+r2s*x1*x3+x1s*x3+x2)
        /prop1s;
        rFO1 = rFO1/2.;
        isSet1 = true;
      }
      if (combi == 2 || combi == 3) {
        rLO2 = ps*(2.-r1s-r1q-6.*r1*r2-r2s+2.*r1s*r2s-r2q)/2.;
        rFO2 = -(3.+6.*r1s+r1q+6.*r1*r2-6.*r1c*r2-2.*r2s-6.*r1s*r2s    
        -6.*r1*r2c+r2q-3.*x1-6.*r1*r2*x1+2.*r2s*x1+x1s-2.*r1s*x1s
        +3.*r1s*x3-6.*r1*r2*x3-r2s*x3-2.*x1*x3-5.*r1s*x1*x3
        +r2s*x1*x3+x1s*x3-3.*x3s-3.*r1s*x3s+r2s*x3s+2.*x1*x3s+x3c-x2)
        /prop2s
        -2.*(-3+r1s+6.*r1*r2-6.*r1c*r2+3.*r2s-4.*r1s*r2s-6.*r1*r2c
        +2.*x1+3.*r1s*x1+r2s*x1-x1s-r1s*x1s-r2s*x1s+4.*x3+2.*r1s*x3
        -3.*r1*r2*x3-r2s*x3-3.*x1*x3-2.*r1s*x1*x3+x1s*x3-x3s-r1s*x3s
        -r1*r2*x3s+x1*x3s)
        /prop12
        -(-1.+2.*r1s+r1q-6.*r1*r2-6.*r1c*r2-2.*r2s-6.*r1s*r2s
        -6.*r1*r2c+r2q-x1-2.*r1s*x1+6.*r1*r2*x1+8.*r2s*x1+x1s
        -2.*r2s*x1s-r1s*x3+r2s*x3-r1s*x1*x3+r2s*x1*x3+x1s*x3+x2)
        /prop1s;
        rFO2 = rFO2/2.;
        isSet2 = true;
      }
      if (combi == 4) {
        rLO4 = ps*(2.-r1s-r1q-r2s+2.*r1s*r2s-r2q)/2.;
        rFO4 = (1.-r1q+6.*r1s*r2s-r2q+x1+3.*r1s*x1-9.*r2s*x1-3.*x1s
        -r1s*x1s+3.*r2s*x1s+x1c-x2-r1s*x2+r2s*x2-r1s*x1*x2+r2s*x1*x2
        +x1s*x2)
        /prop1s 
        -2.*(1.+r1s+r2s-4.*r1s*r2s+r1s*x1+2.*r2s*x1-x1s-r2s*x1s
        +2.*r1s*x2+r2s*x2-3.*x1*x2+x1s*x2-x2s-r1s*x2s+x1*x2s)
        /prop12
        +(1.-r1q+6.*r1s*r2s-r2q-x1+r1s*x1-r2s*x1+x2-9.*r1s*x2
        +3.*r2s*x2+r1s*x1*x2-r2s*x1*x2-3.*x2s+3.*r1s*x2s-r2s*x2s
        +x1*x2s+x2c)
        /prop2s;
        rFO4 = rFO4/2.;
        isSet4 = true;
      }
      break; 
 
    // q -> q V.
    case 3:
      if (combi == 1 || combi == 3) {
        rLO1 = ps*(1.-2.*r1s+r1q+r2s-6.*r1*r2s+r1s*r2s-2.*r2q);
        rFO1 = -2.*(-1.+r1-2.*r1s+2.*r1c-r1q+pow(r1,5)-r2s+r1*r2s
        -5.*r1s*r2s+r1c*r2s-2.*r1*r2q+2.*x1-2.*r1*x1+2.*r1s*x1
        -2.*r1c*x1+2.*r2s*x1+5.*r1*r2s*x1+r1s*r2s*x1+2.*r2q*x1
        -x1s+r1*x1s-r2s*x1s+3.*x2+4.*r1s*x2+r1q*x2+2.*r2s*x2
        +2.*r1s*r2s*x2-4.*x1*x2-2.*r1s*x1*x2-r2s*x1*x2+x1s*x2
        -2.*x2s-2.*r1s*x2s+x1*x2s)
        /prop23
        +(2.*r2s+6.*r1*r2s-6.*r1s*r2s+6.*r1c*r2s+2.*r2q+6.*r1*r2q
        -r2s*x1+r1s*r2s*x1-r2q*x1+x2-r1q*x2-3.*r2s*x2-6.*r1*r2s*x2
        +9.*r1s*r2s*x2-2.*r2q*x2-x1*x2+r1s*x1*x2-x2s-3.*r1s*x2s
        +2.*r2s*x2s+x1*x2s)
        /prop2s
        +(-4.-8.*r1s-4.*r1q+4.*r2s-4.*r1s*r2s+8.*r2q+9.*x1+10.*r1s*x1
        +r1q*x1-3.*r2s*x1+6.*r1*r2s*x1+r1s*r2s*x1-2.*r2q*x1-6.*x1s-
        2.*r1s*x1s+x1c+7.*x2+8.*r1s*x2+r1q*x2-7.*r2s*x2+6.*r1*r2s*x2
        +r1s*r2s*x2-2.*r2q*x2-9.*x1*x2-3.*r1s*x1*x2+2.*r2s*x1*x2
        +2.*x1s*x2-3.*x2s-r1s*x2s+2.*r2s*x2s+x1*x2s)
	/x3s;
        isSet1 = true;
      }
      if (combi == 2 || combi == 3) {
        rLO2 = ps*(1.-2.*r1s+r1q+r2s+6.*r1*r2s+r1s*r2s-2.*r2q);
        rFO2 = 2*(1.+r1+2.*r1s+2.*r1c+r1q+pow(r1,5)+r2s+r1*r2s
        +5.*r1s*r2s+r1c*r2s-2.*r1*r2q-2.*x1-2.*r1*x1-2.*r1s*x1
        -2.*r1c*x1-2.*r2s*x1+5.*r1*r2s*x1-r1s*r2s*x1-2.*r2q*x1+x1s
        +r1*x1s+r2s*x1s-3.*x2-4.*r1s*x2-r1q*x2-2.*r2s*x2
        -2.*r1s*r2s*x2+4.*x1*x2+2.*r1s*x1*x2+r2s*x1*x2-x1s*x2
        +2.*x2s+2.*r1s*x2s-x1*x2s)
        /prop23
        +(2.*r2s-6.*r1*r2s-6.*r1s*r2s-6.*r1c*r2s+2.*r2q-6.*r1*r2q
        -r2s*x1+r1s*r2s*x1-r2q*x1+x2-r1q*x2-3.*r2s*x2+6.*r1*r2s*x2
        +9.*r1s*r2s*x2-2.*r2q*x2-x1*x2+r1s*x1*x2-x2s-3.*r1s*x2s
        +2.*r2s*x2s+x1*x2s)
        /prop2s
        +(-4.-8.*r1s-4.*r1q+4.*r2s-4.*r1s*r2s+8.*r2q+9.*x1+10.*r1s*x1
        +r1q*x1-3.*r2s*x1-6.*r1*r2s*x1+r1s*r2s*x1-2.*r2q*x1-6.*x1s
        -2.*r1s*x1s+x1c+7.*x2+8.*r1s*x2+r1q*x2-7.*r2s*x2-6.*r1*r2s*x2
        +r1s*r2s*x2-2.*r2q*x2-9.*x1*x2-3.*r1s*x1*x2+2.*r2s*x1*x2
        +2.*x1s*x2-3.*x2s-r1s*x2s+2.*r2s*x2s+x1*x2s)
	/x3s;
        isSet2 = true;
      }
      if (combi == 4) {
        rLO4 = ps*(1.-2.*r1s+r1q+r2s+r1s*r2s-2.*r2q);
        rFO4 = 2*(1.+2.*r1s+r1q+r2s+5.*r1s*r2s-2.*x1-2.*r1s*x1
        -2.*r2s*x1-r1s*r2s*x1-2.*r2q*x1+x1s+r2s*x1s-3.*x2-4.*r1s*x2
        -r1q*x2-2.*r2s*x2-2.*r1s*r2s*x2+4.*x1*x2+2.*r1s*x1*x2+r2s*x1*x2
        -x1s*x2+2.*x2s+2.*r1s*x2s-x1*x2s)
        /prop23
        +(2.*r2s-6.*r1s*r2s+2.*r2q-r2s*x1+r1s*r2s*x1-r2q*x1+x2-r1q*x2
        -3.*r2s*x2+9.*r1s*r2s*x2-2.*r2q*x2-x1*x2+r1s*x1*x2-x2s-3.*r1s*x2s
        +2.*r2s*x2s+x1*x2s)
        /prop2s
        +(-4.-8.*r1s-4.*r1q+4.*r2s-4.*r1s*r2s+8.*r2q+9.*x1+10.*r1s*x1
        +r1q*x1-3.*r2s*x1+r1s*r2s*x1-2.*r2q*x1-6.*x1s-2.*r1s*x1s+x1c
        +7.*x2+8.*r1s*x2+r1q*x2-7.*r2s*x2+r1s*r2s*x2-2.*r2q*x2-9.*x1*x2
        -3.*r1s*x1*x2+2.*r2s*x1*x2+2.*x1s*x2-3.*x2s-r1s*x2s+2.*r2s*x2s
        +x1*x2s)
        /x3s;
        isSet4 = true;
      }
      break; 
 
    // S -> q qbar    (S = h0/H0/A0/H+-/...).
    case 4:
      if (combi == 1 || combi == 3) {
        rLO1 = ps*(1.-r1s-r2s-2.*r1*r2);
        rFO1 = -(-1.+r1q-2.*r1*r2-2.*r1c*r2-6.*r1s*r2s-2.*r1*r2c+r2q+x1
        -r1s*x1+2.*r1*r2*x1+3.*r2s*x1+x2+r1s*x2-r2s*x2-x1*x2)
        /prop1s
        -2.*(r1s+r1q-2.*r1c*r2+r2s-6.*r1s*r2s-2.*r1*r2c+r2q-r1s*x1
        +r1*r2*x1+2.*r2s*x1+2.*r1s*x2+r1*r2*x2-r2s*x2-x1*x2)
        /prop12
        -(-1.+r1q-2.*r1*r2-2.*r1c*r2-6.*r1s*r2s-2.*r1*r2c+r2q+x1-r1s*x1
        +r2s*x1+x2+3.*r1s*x2+2.*r1*r2*x2-r2s*x2-x1*x2)
        /prop2s;
        isSet1 = true;
      }
      if (combi == 2 || combi == 3) {
        rLO2 = ps*(1.-r1s-r2s+2.*r1*r2);
        rFO2 = -(-1.+r1q+2.*r1*r2+2.*r1c*r2-6.*r1s*r2s+2.*r1*r2c+r2q+x1
        -r1s*x1-2.*r1*r2*x1+3.*r2s*x1+x2+r1s*x2-r2s*x2-x1*x2)
        /prop1s
        -(-1.+r1q+2.*r1*r2+2.*r1c*r2-6.*r1s*r2s+2.*r1*r2c+r2q+x1
        -r1s*x1+r2s*x1+x2+3.*r1s*x2-2.*r1*r2*x2-r2s*x2-x1*x2)
        /prop2s
        +2.*(-r1s-r1q-2.*r1c*r2-r2s+6.*r1s*r2s-2.*r1*r2c-r2q+r1s*x1
        +r1*r2*x1-2.*r2s*x1-2.*r1s*x2+r1*r2*x2+r2s*x2+x1*x2)
        /prop12;
        isSet2 = true;
      }
      if (combi == 4) {
        rLO4 = ps*(1.-r1s-r2s);
        rFO4 = -(-1.+r1q-6.*r1s*r2s+r2q+x1-r1s*x1+3.*r2s*x1+x2
        +r1s*x2-r2s*x2-x1*x2)
        /prop1s
        -2.*(r1s+r1q+r2s-6.*r1s*r2s+r2q-r1s*x1
        +2.*r2s*x1+2.*r1s*x2-r2s*x2-x1*x2)
        /prop12
        -(-1.+r1q-6.*r1s*r2s+r2q+x1-r1s*x1+r2s*x1
        +x2+3.*r1s*x2-r2s*x2-x1*x2)
        /prop2s;
        isSet4 = true;
      }
      break; 
 
    // q -> q S.
    case 5:
      if (combi == 1 || combi == 3) {
        rLO1 = ps*(1.+r1s-r2s+2.*r1);
        rFO1 = (4.-4.*r1s+4.*r2s-3.*x1-2.*r1*x1+r1s*x1-r2s*x1-5.*x2
        -2.*r1*x2+r1s*x2-r2s*x2+x1*x2+x2s)
        /x3s
        -2.*(3.-r1-5.*r1s-r1c+3.*r2s+r1*r2s-2.*x1-r1*x1
        +r1s*x1-4.*x2+2.*r1s*x2-r2s*x2+x1*x2+x2s)
        /prop23
        +(2.-2.*r1-6.*r1s-2.*r1c+2.*r2s-2.*r1*r2s-x1+r1s*x1
        -r2s*x1-3.*x2+2.*r1*x2+3.*r1s*x2-r2s*x2+x1*x2+x2s)
        /prop2s;
        isSet1 = true;
      }
      if (combi == 2 || combi == 3) {
        rLO2 = ps*(1.+r1s-r2s-2.*r1);
        rFO2 = (4.-4.*r1s+4.*r2s-3.*x1+2.*r1*x1+r1s*x1-r2s*x1-5.*x2
        +2.*r1*x2+r1s*x2-r2s*x2+x1*x2+x2s)
        /x3s
        -2.*(3.+r1-5.*r1s+r1c+3.*r2s-r1*r2s-2.*x1+r1*x1
        +r1s*x1-4.*x2+2.*r1s*x2-r2s*x2+x1*x2+x2s)
        /prop23
        +(2.+2.*r1-6.*r1s+2.*r1c+2.*r2s+2.*r1*r2s-x1+r1s*x1
        -r2s*x1-3.*x2-2.*r1*x2+3.*r1s*x2-r2s*x2+x1*x2+x2s)
        /prop2s;
        isSet2 = true;
      }
      if (combi == 4) {
        rLO4 = ps*(1.+r1s-r2s);
        rFO4 = (4.-4.*r1s+4.*r2s-3.*x1+r1s*x1-r2s*x1-5.*x2+r1s*x2
        -r2s*x2+x1*x2+x2s)
        /x3s
        -2.*(3.-5.*r1s+3.*r2s-2.*x1+r1s*x1-4.*x2+2.*r1s*x2
        -r2s*x2+x1*x2+x2s)
        /prop23
        +(2.-6.*r1s+2.*r2s-x1+r1s*x1-r2s*x1-3.*x2+3.*r1s*x2
        -r2s*x2+x1*x2+x2s)
        /prop2s;
        isSet4 = true;
      }
      break; 
 
    // V -> ~q ~qbar  (~q = squark).
    case 6:
      rLO1 = ps*(1.-2.*r1s+r1q-2.*r2s-2.*r1s*r2s+r2q);
      rFO1 = 2.*3.+(1.+r1s+r2s-x1)*(4.*r1s-x1s)
      /prop1s
      +2.*(-1.-3.*r1s-r2s+x1+x1s*0.5+x2-x1*x2*0.5)
      /prop1
      +(1.+r1s+r2s-x2)*(4.*r2s-x2s)
      /prop2s
      +2.*(-1.-r1s-3.*r2s+x1+x2-x1*x2*0.5+x2s*0.5)
      /prop2
      -(-4.*r1s-4.*r1q-4.*r2s-8.*r1s*r2s-4.*r2q+2.*x1+6.*r1s*x1
      +6.*r2s*x1-2.*x1s+2.*x2+6.*r1s*x2+6.*r2s*x2-4.*x1*x2
      -2.*r1s*x1*x2-2.*r2s*x1*x2+x1s*x2-2.*x2s+x1*x2s)
      /prop12;
      isSet1 = true;
      break; 
 
    // ~q -> ~q V.
    case 7:
      rLO1 = ps*(1.-2.*r1s+r1q-2.*r2s-2.*r1s*r2s+r2q);
      rFO1 = 16.*r2s-8.*(4.*r2s+2.*r2s*x1+x2+r1s*x2+r2s*x2-x1*x2
      -2.*x2s)
      /(3.*prop2)
      +8.*(1.+r1s+r2s-x2)*(4.*r2s-x2s)
      /(3.*prop2s)
      +8.*(x1+x2)*(-1.-2.*r1s-r1q-2.*r2s+2.*r1s*r2s-r2q+2.*x1
      +2.*r1s*x1+2.*r2s*x1-x1s+2.*x2+2.*r1s*x2+2.*r2s*x2-2.*x1*x2-x2s)
      /(3.*x3s)
      +8.*(-1.-r1s+r2s-x1)*(2.*r2s*x1+x2+r1s*x2+r2s*x2-x1*x2-x2s)
      /(3.*prop2*x3)
      -8.*(1.+2.*r1s+r1q+2.*r2s-2.*r1s*r2s+r2q-2.*x1-2.*r1s*x1
      -4.*r2s*x1+x1s-3.*x2-3.*r1s*x2-3.*r2s*x2+3.*x1*x2+2.*x2s)
      /(3.*x3);
      rFO1 = 3.*rFO1/8.;
      isSet1 = true;
      break; 
 
    // S -> ~q ~qbar.
    case 8:
      rLO1 = ps;
      rFO1 = (-1.-2.*r1s-r1q-2.*r2s+2.*r1s*r2s-r2q+2.*x1+2.*r1s*x1
      +2.*r2s*x1-x1s-r2s*x1s+2.*x2+2.*r1s*x2+2.*r2s*x2-3.*x1*x2
      -r1s*x1*x2-r2s*x1*x2+x1s*x2-x2s-r1s*x2s+x1*x2s)
      /(prop1s*prop2s);
      rFO1 = 2.*rFO1;
      isSet1 = true;
      break; 
 
    // ~q -> ~q S.
    case 9:
      rLO1 = ps;
      rFO1 = (-1.-r1s-r2s+x2)
      /prop2s
      +(1.+r1s-r2s+x1)
      /prop23
      -(x1+x2)
      /x3s;
      isSet1 = true;
      break; 
 
    // chi -> q ~qbar   (chi = neutralino/chargino).
    case 10:
      if (combi == 1 || combi == 3) {
        rLO1 = ps*(1.+r1s-r2s+2.*r1);
        rFO1 = (2.*r1+x1)*(-1.-r1s-r2s+x1)
        /prop1s
        +2.*(-1.-r1s-2.*r1c-r2s-2.*r1*r2s+3.*x1*0.5+r1*x1
        -r1s*x1*0.5-r2s*x1*0.5+x2+r1*x2+r1s*x2-x1*x2*0.5)
        /prop12
        +(2.-2.*r1-6.*r1s-2.*r1c+2.*r2s-2.*r1*r2s-x1+r1s*x1
        -r2s*x1-3.*x2+2.*r1*x2+3.*r1s*x2-r2s*x2+x1*x2+x2s)
        /prop2s;
        isSet1 = true;
      }
      if (combi == 2 || combi == 3) {
        rLO2 = ps*(1.-2.*r1+r1s-r2s);
        rFO2 = (2.*r1-x1)*(1.+r1s+r2s-x1)
        /prop1s
        +2.*(-1.-r1s+2.*r1c-r2s+2.*r1*r2s+3.*x1*0.5-r1*x1
        -r1s*x1*0.5-r2s*x1*0.5+x2-r1*x2+r1s*x2-x1*x2*0.5)
        /prop12
        +(2.+2.*r1-6.*r1s+2.*r1c+2.*r2s+2.*r1*r2s-x1+r1s*x1
        -r2s*x1-3.*x2-2.*r1*x2+3.*r1s*x2-r2s*x2+x1*x2+x2s)/
        prop2s;
        isSet2 = true;
      }
      if (combi == 4) {
        rLO4 = ps*(1.+r1s-r2s);
        rFO4 = x1*(-1.-r1s-r2s+x1)
        /prop1s
        +2.*(-1.-r1s-r2s+3.*x1*0.5-r1s*x1*0.5-r2s*x1*0.5
        +x2+r1s*x2-x1*x2*0.5)
        /prop12
        +(2.-6.*r1s+2.*r2s-x1+r1s*x1-r2s*x1-3.*x2+3.*r1s*x2
        -r2s*x2+x1*x2+x2s)
        /prop2s;
        isSet4 = true;
      }
      break; 
 
    // ~q -> q chi.
    case 11:
      if (combi == 1 || combi == 3) {
        rLO1 = ps*(1.-pow(r1+r2,2));
        rFO1 = (1.+r1s+2.*r1*r2+r2s-x1-x2)*(x1+x2)
        /x3s
        -(-1.+r1q-2.*r1*r2-2.*r1c*r2-6.*r1s*r2s-2.*r1*r2c+r2q+x1
        -r1s*x1+r2s*x1+x2+3.*r1s*x2+2.*r1*r2*x2-r2s*x2-x1*x2)
        /prop2s
        +(-1.-2.*r1s-r1q-2.*r1*r2-2.*r1c*r2+2.*r1*r2c+r2q+x1+r1s*x1
        -2.*r1*r2*x1-3.*r2s*x1+2.*r1s*x2-2.*r2s*x2+x1*x2+x2s)
        /prop23;
        isSet1 = true;
      }
      if (combi == 2 || combi == 3) {
        rLO2 = ps*(1.-pow(r1-r2,2));
        rFO2 = (1.+r1s-2.*r1*r2+r2s-x1-x2)*(x1+x2)
        /x3s
        -(-1.+r1q+2.*r1*r2+2.*r1c*r2-6.*r1s*r2s+2.*r1*r2c+r2q+x1
        -r1s*x1+r2s*x1+x2+3.*r1s*x2-2.*r1*r2*x2-r2s*x2-x1*x2)
        /prop2s
        +(-1.-2.*r1s-r1q+2.*r1*r2+2.*r1c*r2-2.*r1*r2c+r2q+x1+r1s*x1
        +2.*r1*r2*x1-3.*r2s*x1+2.*r1s*x2-2.*r2s*x2+x1*x2+x2s)
        /prop23;
        isSet2 = true;
      }
      if (combi == 4) {
        rLO4 = ps*(1.-r1s-r2s);
        rFO4 = (1.+r1s+r2s-x1-x2)*(x1+x2)
        /x3s
        -(-1.+r1q-6.*r1s*r2s+r2q+x1-r1s*x1+r2s*x1+x2
        +3.*r1s*x2-r2s*x2-x1*x2)
        /prop2s
        +(-1.-2.*r1s-r1q+r2q+x1+r1s*x1-3.*r2s*x1
        +2.*r1s*x2-2.*r2s*x2+x1*x2+x2s)
        /prop23;
        isSet4 = true;
      }
      break; 
 
    // q -> ~q chi.
    case 12:
      if (combi == 1 || combi == 3) {
        rLO1 = ps*(1.-r1s+r2s+2.*r2);
        rFO1 = (2.*r2+x2)*(-1.-r1s-r2s+x2)
        /prop2s
        +(4.+4.*r1s-4.*r2s-5.*x1-r1s*x1-2.*r2*x1+r2s*x1+x1s
        -3.*x2-r1s*x2-2.*r2*x2+r2s*x2+x1*x2)
        /x3s
        +2.*(-1.-r1s+r2+r1s*r2-r2s-r2c+x1+r2*x1+r2s*x1+2.*x2
        +r1s*x2-x1*x2*0.5-x2s*0.5)
        /prop23;
        isSet1 = true;
      }
      if (combi == 2 || combi == 3) {
        rLO2 = ps*(1.-r1s+r2s-2.*r2);
        rFO2 = (2.*r2-x2)*(1.+r1s+r2s-x2)
        /prop2s
        +(4.+4.*r1s-4.*r2s-5.*x1-r1s*x1+2.*r2*x1+r2s*x1+x1s
        -3.*x2-r1s*x2+2.*r2*x2+r2s*x2+x1*x2)
        /x3s
        +2.*(-1.-r1s-r2-r1s*r2-r2s+r2c+x1-r2*x1+r2s*x1+2.*x2
        +r1s*x2-x1*x2*0.5-x2s*0.5)
        /prop23;
        isSet2 = true;
      }
      if (combi == 4) {
        rLO4 = ps*(1.-r1s+r2s);
        rFO4 = x2*(-1.-r1s-r2s+x2)
        /prop2s
        +(4.+4.*r1s-4.*r2s-5.*x1-r1s*x1+r2s*x1+x1s
        -3.*x2-r1s*x2+r2s*x2+x1*x2)
        /x3s
        +2.*(-1.-r1s-r2s+x1+r2s*x1+2.*x2
        +r1s*x2-x1*x2*0.5-x2s*0.5)
        /prop23;
        isSet4 = true;
      }
      break; 
 
    // ~g -> q ~qbar.
    case 13:
      if (combi == 1 || combi == 3) {
        rLO1 = ps*(1.+r1s-r2s+2.*r1);
        rFO1 = 4.*(2.*r1+x1)*(-1.-r1s-r2s+x1)
        /(3.*prop1s)
        -(-1.-r1s-2.*r1c-r2s-2.*r1*r2s+3.*x1*0.5+r1*x1-r1s*x1*0.5
        -r2s*x1*0.5+x2+r1*x2+r1s*x2-x1*x2*0.5)
        /(3.*prop12)
        +3.*(-1.+r1-r1s-r1c-r2s+r1*r2s+2.*x1+r2s*x1-x1s*0.5+x2+r1*x2
        +r1s*x2-x1*x2*0.5)
        /prop13
        +3.*(4.-4.*r1s+4.*r2s-3.*x1-2.*r1*x1+r1s*x1-r2s*x1-5.*x2
        -2.*r1*x2+r1s*x2-r2s*x2+x1*x2+x2s)
        /x3s
        -3.*(3.-r1-5.*r1s-r1c+3.*r2s+r1*r2s-2.*x1-r1*x1+r1s*x1
        -4.*x2+2.*r1s*x2-r2s*x2+x1*x2+x2s)
        /prop23
        +4.*(2.-2.*r1-6.*r1s-2.*r1c+2.*r2s-2.*r1*r2s-x1+r1s*x1-r2s*x1
        -3.*x2+2.*r1*x2+3.*r1s*x2-r2s*x2+x1*x2+x2s)
        /(3.*prop2s);
        rFO1 = 3.*rFO1/4.;
        isSet1 = true;
      }
      if (combi == 2 || combi == 3) {
        rLO2 = ps*(1.+r1s-r2s-2.*r1);
        rFO2 = 4.*(2.*r1-x1)*(1.+r1s+r2s-x1)
        /(3.*prop1s)
        +3.*(-1.-r1-r1s+r1c-r2s-r1*r2s+2.*x1+r2s*x1-x1s*0.5
        +x2-r1*x2+r1s*x2-x1*x2*0.5)
        /prop13
        +(2.+2.*r1s-4.*r1c+2.*r2s-4.*r1*r2s-3.*x1+2.*r1*x1
        +r1s*x1+r2s*x1-2.*x2+2.*r1*x2-2.*r1s*x2+x1*x2)
        /(6.*prop12)
        +3.*(4.-4.*r1s+4.*r2s-3.*x1+2.*r1*x1+r1s*x1-r2s*x1-5.*x2
        +2.*r1*x2+r1s*x2-r2s*x2+x1*x2+x2s)
        /x3s
        -3.*(3.+r1-5.*r1s+r1c+3.*r2s-r1*r2s-2.*x1+r1*x1+r1s*x1-4.*x2
        +2.*r1s*x2-r2s*x2+x1*x2+x2s)
        /prop23
        +4.*(2.+2.*r1-6.*r1s+2.*r1c+2.*r2s+2.*r1*r2s-x1+r1s*x1-r2s*x1
        -3.*x2-2.*r1*x2+3.*r1s*x2-r2s*x2+x1*x2+x2s)
        /(3.*prop2s);
        rFO2 = 3.*rFO2/4.;
        isSet2 = true;
      }
      if (combi == 4) {
        rLO4 = ps*(1.+r1s-r2s);
        rFO4 = 8.*x1*(-1.-r1s-r2s+x1)
        /(3.*prop1s)
        +6.*(-1-r1s-r2s+2.*x1+r2s*x1-x1s*0.5+x2+r1s*x2-x1*x2*0.5)
        /prop13
        +(2.+2.*r1s+2.*r2s-3.*x1+r1s*x1+r2s*x1-2.*x2-2.*r1s*x2+x1*x2)
        /(3.*prop12)
        +6.*(4.-4.*r1s+4.*r2s-3.*x1+r1s*x1-r2s*x1-5.*x2+r1s*x2-r2s*x2
        +x1*x2+x2s)
        /x3s
        -6.*(3.-5.*r1s+3.*r2s-2.*x1+r1s*x1-4.*x2+2.*r1s*x2-r2s*x2+x1*x2+x2s)
        /prop23
        +8.*(2.-6.*r1s+2.*r2s-x1+r1s*x1-r2s*x1-3.*x2+3.*r1s*x2-r2s*x2
        +x1*x2+x2s)
        /(3.*prop2s);
        rFO4 = 3.*rFO4/8.;
        isSet4 = true;
      }
      break; 
 
    // ~q -> q ~g.
    case 14:
      if (combi == 1 || combi == 3) {
        rLO1 = ps*(1.-r1s-r2s-2.*r1*r2);
        rFO1 = 64.*(1.+r1s+2.*r1*r2+r2s-x1-x2)*(x1+x2)
        /(9.*x3s)
        -16.*(-1.+r1q-2.*r1*r2-2.*r1c*r2-6.*r1s*r2s-2.*r1*r2c+r2q
        +x1-r1s*x1+2.*r1*r2*x1+3.*r2s*x1+x2+r1s*x2-r2s*x2-x1*x2)
        /prop1s
        -16.*(r1s+r1q-2.*r1c*r2+r2s-6.*r1s*r2s-2.*r1*r2c+r2q-r1s*x1
        +r1*r2*x1+2.*r2s*x1+2.*r1s*x2+r1*r2*x2-r2s*x2-x1*x2)
        /prop12
        -64.*(-1.+r1q-2.*r1*r2-2.*r1c*r2-6.*r1s*r2s-2.*r1*r2c+r2q+x1
        -r1s*x1+r2s*x1+x2+3.*r1s*x2+2.*r1*r2*x2-r2s*x2-x1*x2)
        /(9.*prop2s)
        +8.*(-1.+r1q-2.*r1*r2+2.*r1c*r2-2.*r2s-2.*r1*r2c-r2q-2.*r1s*x1
        +2.*r2s*x1+x1s+x2-3.*r1s*x2-2.*r1*r2*x2+r2s*x2+x1*x2)
        /prop13
        -8.*(-1.-2.*r1s-r1q-2.*r1*r2-2.*r1c*r2+2.*r1*r2c+r2q+x1+r1s*x1
        -2.*r1*r2*x1-3.*r2s*x1+2.*r1s*x2-2.*r2s*x2+x1*x2+x2s)
        /(9.*prop23);
        rFO1 = 9.*rFO1/64.;
        isSet1 = true;
      }
      if (combi == 2 || combi == 3) {
        rLO2 = ps*(1.-r1s-r2s+2.*r1*r2);
        rFO2 = 64.*(1.+r1s-2.*r1*r2+r2s-x1-x2)*(x1+x2)
        /(9.*x3s)
        -16.*(-1.+r1q+2.*r1*r2+2.*r1c*r2-6.*r1s*r2s+2.*r1*r2c+r2q+x1
        -r1s*x1-2.*r1*r2*x1+3.*r2s*x1+x2+r1s*x2-r2s*x2-x1*x2)
        /prop1s
        -64.*(-1.+r1q+2.*r1*r2+2.*r1c*r2-6.*r1s*r2s+2.*r1*r2c+r2q+x1
        -r1s*x1+r2s*x1+x2+3.*r1s*x2-2.*r1*r2*x2-r2s*x2-x1*x2)
        /(9.*prop2s)
        +16.*(-r1s-r1q-2.*r1c*r2-r2s+6.*r1s*r2s-2.*r1*r2c-r2q+r1s*x1
        +r1*r2*x1-2.*r2s*x1-2.*r1s*x2+r1*r2*x2+r2s*x2+x1*x2)
        /prop12
        +8.*(-1.+r1q+2.*r1*r2-2.*r1c*r2-2.*r2s+2.*r1*r2c-r2q-2.*r1s*x1
        +2.*r2s*x1+x1s+x2-3.*r1s*x2+2.*r1*r2*x2+r2s*x2+x1*x2)
        /prop13
        -8.*(-1.-2.*r1s-r1q+2.*r1*r2+2.*r1c*r2-2.*r1*r2c+r2q+x1+r1s*x1+
        2.*r1*r2*x1-3.*r2s*x1+2.*r1s*x2-2.*r2s*x2+x1*x2+x2s)
        /(9.*prop23);
        rFO2 = 9.*rFO2/64.;
        isSet2 = true;
      }
      if (combi == 4) {
        rLO4 = ps*(1.-r1s-r2s);
        rFO4 = 128.*(1.+r1s+r2s-x1-x2)*(x1+x2)
        /(9.*x3s)
        -32*(-1.+r1q-6.*r1s*r2s+r2q+x1-r1s*x1+3.*r2s*x1+x2
        +r1s*x2-r2s*x2-x1*x2)
        /prop1s
        -32.*(r1s+r1q+r2s-6.*r1s*r2s+r2q-r1s*x1+2.*r2s*x1+2.*r1s*x2
        -r2s*x2-x1*x2)
        /prop12
        -128.*(-1.+r1q-6.*r1s*r2s+r2q+x1-r1s*x1+r2s*x1+x2+3.*r1s*x2
        -r2s*x2-x1*x2)
        /(9.*prop2s)
        +16.*(-1.+r1q-2.*r2s-r2q-2.*r1s*x1+2.*r2s*x1+x1s
        +x2-3.*r1s*x2+r2s*x2+x1*x2)
        /prop13
        -16.*(-1.-2.*r1s-r1q+r2q+x1+r1s*x1-3.*r2s*x1
        +2.*r1s*x2-2.*r2s*x2+x1*x2+x2s)
        /(9.*prop23);
        rFO4 = 9.*rFO4/128.;
        isSet4 = true;
      }
      break; 
 
    // q -> ~q ~g.
    case 15:
      if (combi == 1 || combi == 3) {
        rLO1 = ps*(1.-r1s+r2s+2.*r2);
        rFO1 = 32*(2.*r2+x2)*(-1.-r1s-r2s+x2)
        /(9.*prop2s)
        +8.*(-1.-r1s-2.*r1s*r2-r2s-2.*r2c+x1+r2*x1+r2s*x1
        +3.*x2*0.5-r1s*x2*0.5+r2*x2-r2s*x2*0.5-x1*x2*0.5)
        /prop12
        +8.*(2.+2.*r1s-2.*r2-2.*r1s*r2-6.*r2s-2.*r2c-3.*x1-r1s*x1
        +2.*r2*x1+3.*r2s*x1+x1s-x2-r1s*x2+r2s*x2+x1*x2)
        /prop1s
        +32.*(4.+4.*r1s-4.*r2s-5.*x1-r1s*x1-2.*r2*x1+r2s*x1+x1s
        -3.*x2-r1s*x2-2.*r2*x2+r2s*x2+x1*x2)
        /(9.*x3s)
        -8.*(3.+3.*r1s-r2+r1s*r2-5.*r2s-r2c-4.*x1-r1s*x1
        +2.*r2s*x1+x1s-2.*x2-r2*x2+r2s*x2+x1*x2)
        /prop13
        -8.*(-1.-r1s+r2+r1s*r2-r2s-r2c+x1+r2*x1+r2s*x1+2.*x2+r1s*x2
        -x1*x2*0.5-x2s*0.5)
        /(9.*prop23);
        rFO1 = 9.*rFO1/32.;
        isSet1 = true;
      }
      if (combi == 2 || combi == 3) {
        rLO2 = ps*(1.-r1s+r2s-2.*r2);
        rFO2 = 32*(2.*r2-x2)*(1.+r1s+r2s-x2)
        /(9.*prop2s)
        +8.*(-1.-r1s+2.*r1s*r2-r2s+2.*r2c+x1-r2*x1+r2s*x1
        +3.*x2*0.5-r1s*x2*0.5-r2*x2-r2s*x2*0.5-x1*x2*0.5)
        /prop12
        +8.*(2.+2.*r1s+2.*r2+2.*r1s*r2-6.*r2s+2.*r2c-3.*x1-r1s*x1
        -2.*r2*x1+3.*r2s*x1+x1s-x2-r1s*x2+r2s*x2+x1*x2)
        /prop1s
        -8.*(3.+3.*r1s+r2-r1s*r2-5.*r2s+r2c-4.*x1-r1s*x1+2.*r2s*x1+x1s
        -2.*x2+r2*x2+r2s*x2+x1*x2)
        /prop13
        +32*(4.+4.*r1s-4.*r2s-5.*x1-r1s*x1+2.*r2*x1+r2s*x1
        +x1s-3.*x2-r1s*x2+2.*r2*x2+r2s*x2+x1*x2)
        /(9.*x3s)
        -8.*(-1.-r1s-r2-r1s*r2-r2s+r2c+x1-r2*x1+r2s*x1+2.*x2+r1s*x2
        -x1*x2*0.5-x2s*0.5)
        /(9.*prop23);
        rFO2 = 9.*rFO2/32.;
        isSet2 = true;
      }
      if (combi == 4) {
        rLO4 = ps*(1.-r1s+r2s);
        rFO4 = 64.*x2*(-1.-r1s-r2s+x2)
        /(9.*prop2s)
        +16.*(-1.-r1s-r2s+x1+r2s*x1+3.*x2*0.5-r1s*x2*0.5
        -r2s*x2*0.5-x1*x2*0.5)
        /prop12
        -16.*(3.+3.*r1s-5.*r2s-4.*x1-r1s*x1+2.*r2s*x1+x1s-2.*x2+r2s*x2
        +x1*x2)
        /prop13
        +64.*(4.+4.*r1s-4.*r2s-5.*x1-r1s*x1+r2s*x1+x1s-3.*x2
        -r1s*x2+r2s*x2+x1*x2)
        /(9.*x3s)
        +16.*(2.+2.*r1s-6.*r2s-3.*x1-r1s*x1+3.*r2s*x1+x1s
        -x2-r1s*x2+r2s*x2+x1*x2)
        /prop1s
        -16.*(-1.-r1s-r2s+x1+r2s*x1+2.*x2+r1s*x2-x1*x2*0.5-x2s*0.5)
        /(9.*prop23);
        rFO4 = 9.*rFO4/64.;
        isSet4 = true;
      }
      break; 
 
    // g -> ~g ~g. Use (9/4)*eikonal. May be changed in the future.
    case 16:
      rLO = ps;
      if (combi == 2) offset = x3s;
      else if (combi == 3) offset = mix * x3s;
      else if (combi == 4) offset = 0.5 * x3s;
      rFO = ps * 4.5 * ( (x1+x2-1.+offset-r1s-r2s)/prop12 
      - r1s/prop2s - r2s/prop1s );
      break; 

    // Eikonal expression for kind == 1; also acts as default.
    default:
      rLO = ps;
      if (combi == 2) offset = x3s;
      else if (combi == 3) offset = mix * x3s;
      else if (combi == 4) offset = 0.5 * x3s;
      rFO = ps * 2. * ( (x1+x2-1.+offset-r1s-r2s)/prop12 
      - r1s/prop2s - r2s/prop1s );
      break;

  // End of ME cases. 
  }

  // Find relevant leading and first order expressions.
  if (combi == 1 && isSet1) {rLO = rLO1; rFO = rFO1;}     
  else if (combi == 2 && isSet2) {rLO = rLO2; rFO = rFO2;}     
  else if (combi == 3 && isSet1 && isSet2) {
    rLO = mix * rLO1 + (1.-mix) * rLO2; 
    rFO = mix * rFO1 + (1.-mix) * rFO2; }
  else if (isSet4) {rLO = rLO4; rFO = rFO4;}     
  else if (combi == 4 && isSet1 && isSet2) {
    rLO = 0.5 * (rLO1 + rLO2);
    rFO = 0.5 * (rFO1 + rFO2); }
  else if (isSet1) {rLO = rLO1; rFO = rFO1;} 

  // Return ratio of first to leading order cross section.     
  return (rFO/rLO);
}  
