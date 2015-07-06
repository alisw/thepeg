// Function definitions (not found in the header) for the
// SpaceShower class, which does spacelike showers. 

#include "Basics.h"
#include "Beam.h"
#include "Shower.h"
#include "SpaceShower.h"
#include "ThePEG/Utilities/LoopGuard.h"

using namespace Pythia7;
using namespace Shower;

//**************************************************************************

// Allowed ISR showers for incoming hadron: = 0: none; = 1: QCD ones;
// = 2 : also allow photon emission;
// = 3 : allow spacelike photon branchings, if supported by the PDF used
// (not the case currently) (does not address VMD part of resolved photons).
long SpaceShower::HADRONSHOWER = 2;

// Allowed ISR showers for incoming lepton: = 0: none; = 1: photon emission;
// = 2 : allow spacelike photon branchings, if supported by the PDF used
// (not the case currently) (does not address VMD part of resolved photons).
long SpaceShower::LEPTONSHOWER = 1;

// Number of allowed quark flavours in g -> q qbar branching.
long SpaceShower::NQUARK = 5;

// Running of alpha_strong in evolution:
// = 0: fixed; = 1: scale Q^2; = 2: scale pT^2. 
// But note that PDF's are always evaluated at Q2.
long SpaceShower::ALPHASMODE = 2;

// Maximum virtuality setting the starting point of the evolution:
// = 0: s; = 1: sHat; = 2: average mT^2; = 3: smallest mT^2.
long SpaceShower::MAXVIRTUALITY = 0;

// Use of matrix element corrections:
// = 0: no; = 1: yes.
long SpaceShower::MEMODE = 1;

// Resum the effect of multiple soft gluon emissions: = 0: no; = 1: yes.
long SpaceShower::SOFTGLUONRESUM = 1;

// Restrict first emission within cone given by colour flow in hard process.
// = 0: no; = 1: yes, isotropic phi angle inside cone; 
// = 2: yes, also with anisotropic phi angle inside cone.
long SpaceShower::FINALCONE = 2;

// Q2 ordering is normal, but could be relaxed (to be developed in future).
// = 0: yes; = 1: no, allow maximum kinematically possible (toy scenario).
long SpaceShower::Q2ORDER = 0;

// Angular ordering; = 0: off, = 1: on.
long SpaceShower::ANGULARORDER = 1;

// Azimuthal asymmetry induced by gluon polarization.
// = 0: no; = 1: yes.
long SpaceShower::PHIPOLASYM = 1;

// Azimuthal asymmetry induced by colour coherence.
// = 0: no; = 1: yes.
long SpaceShower::PHICOHERASYM = 1;

// Use the scale variable of original partons to restrict branchings.
// = 0: no; = 1: yes, the Q2 < scale; = 2: yes, the pT2 < scale,
// = 3: yes, the (E*theta)^2 < scale; = 4: yes, the theta^2 < scale.
// Here theta is the full opening angle of a branching, defined in the
// rest frame of the event. (In all cases relations are approximate.) 
long SpaceShower::RESPECTSCALE = 0;

// Parton shower cut-off mass for QCD emissions.
double SpaceShower::Q0 = 1.;

// Parton shower cut-off mass for photon coupling to coloured particle.
double SpaceShower::Q0CHGQ = 1.;

// Parton shower cut-off mass for pure QED branchings. Assumed <= Q0CHGQ.
double SpaceShower::Q0CHGL = 0.001; 

// Fixed alpha_strong value for ALPHASMODE == 0. 
double SpaceShower::ALPHASFIX = 0.2;

// Lambda_QCD(five flavours) in alpha_strong for ALPHASMODE >= 1 .
double SpaceShower::LAMBDA5 = 0.18;

// Fixed alpha_em value. 
double SpaceShower::ALPHAEMFIX = 0.0073;

// Minimum energy of emitted QCD parton in rest frame of subprocess.
double SpaceShower:: EMINEMITTED = 2.;

// Minimum fraction 1 - z of emitted QCD parton, in addition to above.
double SpaceShower:: ZMINEMITTED = 0.001;

// Minimum x fraction of emitted photon - matched to treatment of photon PDF..
double SpaceShower::XMINEMITTEDCHG = 1E-10;

//*********

// Following should not be interfaced for public changing.

// Smallest particle mass for QED evolution (= electron mass).
double SpaceShower::TINYQCHG = 0.00051;

// Vanishingly small parton density.
double SpaceShower::TINYPDF = 1e-10;

// Vanishingly small product of splitting kernels and parton density ratios.
double SpaceShower::TINYKERNELPDF = 1e-5;

// Vanishingly small recoil mass in branching kinematics reconstruction. 
double SpaceShower::TINYKINPREC = 1e-10;

// Safety margin in x that heavy flavour evolution is at all possible.
double SpaceShower::HEAVYEVOL = 0.8;

// Extra preweight in QED shower evolution, to avoid maximum violation.
double SpaceShower::EXTRAQEDPREWT = 0.1;

// Maximum allowed x when reconstructing back to heavy flavour from 
// gluon or photon.
double SpaceShower::HEAVYXMAX = 0.5;

// Mimimum gap in Q2 values to allow iteration when parton density vanishes.
double SpaceShower::Q2STARTFRAC = 0.9;

//*********

// Top-level driver routine to do a single space-like shower.

void SpaceShower::shower(Event& event, BeamParticle& beam1in, 
  BeamParticle& beam2in, long in1in, long in2in, long MEkindin) {

  // Rare loopback if kinematics reconstruction fails.
  bool isOKkinem; 
  ThePEG::LoopGuard<Veto,void> loopguard(10000);
  do {
    loopguard();
    isOKkinem = true;

    // Read in info on system to be treated.
    read(event, beam1in, beam2in, in1in, in2in);
 
    // Prepare to do shower; skip if no particles may shower.
    if (setUpPrimary(MEkindin) ) {

      // Evolve two colliding partons independently of each other.
      if (entry1[0].canBranch) evolveParton(1); 
      if (entry2[0].canBranch) evolveParton(2); 

      // Construct primary kinematics.
      kinemPrimary();

      // Pick side with highest virtuality, so long as any virtuality left.
      long side;
      while (isOKkinem && (side = pickSide()) > 0) {

        // Prepare to let mother of this parton evolve.
        setUpBranching(side);

        // Evolve  the mother parton on this side.
        evolveParton(side);

        // Construct kinematics for branching.
        isOKkinem = kinemBranching(side);
      }  
    } 

  // Rare loopback if kinematics reconstruction fails.
  } while (!isOKkinem);

  // Write back the generated shower.
  write(event);
  
}

//*********

// Read in info on system to be treated.

void SpaceShower::read(Event& event, BeamParticle& beam1in, 
  BeamParticle& beam2in, long in1in, long in2in) {

  // Read in beams, including parton densities.
  beam1 = &beam1in;
  beam2 = &beam2in;

  // Read in or find positions of colliding partons.
  in1 = in2 = -1;
  if (in1in >= 0 && in2in >= 0 ) {
    in1 = in1in; in2 = in2in;
  } else {
    for (long i = 0; i < event.size(); ++i) {   
      if (event[i].status() == -1) {
	if (in1 == -1) in1 = i;
        else {in2 = i; break;}
      }
    }
  } 

  // Find s and sHat.
  Vec4 pSum = beam1->p() + beam2->p();
  s = pSum.m2calc();
  eCM = sqrtpos(s);
  Vec4 pSub = event[in1].p() + event[in2].p(); 
  sHat = pSub.m2calc();
  mHat = sqrtpos(sHat);
  sHatFinalState = sHat;

  // Find x1 and x2 of the incoming partons.
  double x1 = (pSub.e() + pSub.pz()) / eCM;
  double x2 = (pSub.e() - pSub.pz()) / eCM;
   
  // Set up two incoming partons.
  entry1.resize(1);
  entry1[0] = event[in1];
  entry1[0].x = x1;
  entry1[0].copyOfPrimary = in1;
  entry2.resize(1);
  entry2[0] = event[in2];
  entry2[0].x = x2;
  entry2[0].copyOfPrimary = in2;

  // Read current largest colour index.
  maxColIndx = event.colIndx();

  // Identify outgoing particles, for possible Matrix Elements correction.
  nFinal = 0;
  avgMT2 = 0.;
  minMT2 = sHat;
  for (long i = in2 + 1; i < event.size(); ++i) {   
    if (event[i].mother1() == in1 && event[i].mother2() == in2) {
      ++nFinal;
      if (nFinal == 1) finalId1 = event[i].id();
      if (nFinal == 2) finalId2 = event[i].id();
      avgMT2 += event[i].mT2();
      minMT2 = min( minMT2, event[i].mT2() ); 
    }
  }
  if (nFinal > 0) avgMT2 /= nFinal;

  // Attach vectors of colour and anticolour daughter or partner.
  for (long side = 1; side <=2; ++side) {
    SpaceParticle* incoming = (side == 1) ? &entry1[0] : &entry2[0];
    SpaceParticle* opposite = (side == 1) ? &entry2[0] : &entry1[0]; 
    if (incoming->col() > 0) {
      if (opposite->anticol() == incoming->col()) 
	  incoming->setColVec4(opposite->p());
      else {
	for (long i = in2 + 1; i < event.size(); ++i) { 
          if (event[i].mother1() == in1 && event[i].mother2() == in2
	    && event[i].col() == incoming->col() )
            incoming->setColVec4(event[i].p());
	}  
      }
    }
    if (incoming->anticol() > 0) {
      if (opposite->col() == incoming->anticol()) 
	  incoming->setAntiColVec4(opposite->p());
      else {
	for (long i = in2 + 1; i < event.size(); ++i) { 
          if (event[i].mother1() == in1 && event[i].mother2() == in2
	    && event[i].anticol() == incoming->anticol() )
            incoming->setAntiColVec4(event[i].p());
	}  
      }
    }
  }

  // Boost to rest frame of showering system.
  entry1[0].bst(0., 0., -pSub.pz() / pSub.e());
  entry2[0].bst(0., 0., -pSub.pz() / pSub.e());
  
} 
    
//*********

// Write back shower after treatment.

void SpaceShower::write(Event& event) {

  // Mark primary partons as treated.
  event[in1].addstatus(-1);
  event[in2].addstatus(-1);

  // Longitudinal boost by showers.
  double betaZ = (entry1.back().x - entry2.back().x) / 
    (entry1.back().x + entry2.back().x);

  // Initial values for writeback logic.
  long i1 = entry1.size();
  long mother1 = event[in1].mother1();
  long i2 = entry2.size();
  long mother2 = event[in2].mother1();
  Particle temp;

  // Pick side, in interleaved order of increasing mother Q2.
  // (But with the hard interaction partons at the end.)
  while (i1 > 0 || i2 > 0) {
    long side = 1;
    if (i1 > 2 && i2 > 2 && entry2[i2-3].Q2Now < entry1[i1-3].Q2Now) side = 2;
    if (i1 <= 2 && i2 > 2) side = 2;
    if (i1 <= 0) side = 2;
   
    // Write back pair of spacelike and timelike parton.
    // (Except that the hard interaction parton is unmatched.)
    if (side == 1) {
      temp = entry1[--i1];
      temp.mothers(mother1, -1);    
      temp.bst(0., 0., betaZ);
      mother1 = event.append(temp);
      if (i1 > 1) { 
        temp = entry1[--i1];
        temp.mothers(mother1, -1);    
        temp.bst(0., 0., betaZ);
        event.append(temp);
      }
    } else {
      temp = entry2[--i2];
      temp.mothers(mother2, -1);    
      temp.bst(0., 0., betaZ);
      mother2 = event.append(temp); 
      if (i2 > 1) { 
        temp = entry2[--i2];
        temp.mothers(mother2, -1);    
        temp.bst(0., 0., betaZ);
        event.append(temp);
      }
    }
  }

  // Also write current largest colour index in use.
  event.colIndx(maxColIndx);

}
     
//*********
  
// Set up primary partons for evolution.

bool SpaceShower::setUpPrimary(long MEkindin) {
  
  // Find ME correction kind.
  hasME = false;
  if (MEMODE >= 1) {
    findMEkind(MEkindin);
    if(MEkind > 0) hasME = true;
  }

  // Loop over two event sides; identify beam and colliding parton.   
  for (long side = 1; side <= 2; ++side) {
    BeamParticle* beam = (side == 1) ? beam1 : beam2;
    SpaceParticle* current = (side == 1) ? &entry1[0] :&entry2[0];  
    double sign = (side == 1) ? 1. : -1.; 

    // Set simple kinematics and default showering status.
    double x = current->x;
    current->p(0., 0., sign * 0.5 * x * eCM, 0.5 * x * eCM); 
    current->canBranch = false;
    current->hasBranched = false;
    //    long id = current->id();  LEIF? UNUSED

    // Check colliding parton from hadron/lepton, to see if it can shower.
    if (beam->isHadron() && HADRONSHOWER > 0) current->canBranch = true; 
    else if (beam->isLepton() && LEPTONSHOWER > 0) current->canBranch = true;

    // If yes, set maximum virtuality.
    if (current->canBranch) {
      if (MAXVIRTUALITY == 1) current->Q2Now = sHat;
      else if (MAXVIRTUALITY == 2) current->Q2Now = avgMT2;
      else if (MAXVIRTUALITY == 3) current->Q2Now = minMT2;
      else current->Q2Now = s;
      if (RESPECTSCALE == 1 && current->scale() < current->Q2Now) 
        current->Q2Now = current->scale();
      current->sin2thetaMax = 1.; 

      // Decide which of colour and anticolour should set maximum emission cone.
      current->coneSide = 0;
      if (FINALCONE >= 1 && (current->hasColVec4 
        || current->hasAntiColVec4)) {
        current->coneSide = 1;
        if (current->hasAntiColVec4 && (!current->hasColVec4 
          || Rndm::flat() > 0.5) ) current->coneSide = 2;
      }
    }
  }

  return true; // LEIF? NO ALTERNATIVE RETURN?
  
}

//*********

// Check and set kinematics for primary partons after evolution.

bool SpaceShower::kinemPrimary() {

  // Construct kinematics.
  double m1s = entry1[0].m2Now();
  double m2s = entry2[0].m2Now();
  double e1 = 0.5 * (sHat + m1s - m2s) / mHat; 
  double e2 = 0.5 * (sHat + m2s - m1s) / mHat; 
  double pAbs = pAbsLambdaK(sHat, m1s, m2s); 
  entry1[0].p(0., 0.,  pAbs, e1); entry1[0].m(entry1[0].mNow());
  entry2[0].p(0., 0., -pAbs, e2); entry2[0].m(entry2[0].mNow());

  return true; // LEIF? NO ALTERNATIVE RETURN?
  
}

//*********

// Pick one of the sides for further evolution when required.
  
long SpaceShower::pickSide() {

  // The simple choices when none or only one can be picked.
  SpaceParticle& latest1 = entry1.back();  
  SpaceParticle& latest2 = entry2.back();  
  if (!latest1.hasBranched && !latest2.hasBranched) return 0;
  if (!latest2.hasBranched) return 1;
  if (!latest1.hasBranched) return 2;

  // Else pick side with largest virtuality.
  if (latest1.Q2Now > latest2.Q2Now) return 1;
  return 2;

}

//*********

// Set up mother of parton branching for her subsequent evolution,
// and also insert on-shell (timelike) sister.
  
void SpaceShower::setUpBranching(long side) {

  // Add sister and mother entries to bottom of shower.
  SpaceParticle* daughter;
  SpaceParticle* mother;
  SpaceParticle* sister;
  SpaceParticle* opposite; // LEIF? UNUSED?
  if (side == 1) {
    long now = entry1.size();
    entry1.push_back(SpaceParticle());
    entry1.push_back(SpaceParticle());
    daughter = &entry1[now-1];
    sister = &entry1[now];
    mother = &entry1[now+1];
    opposite = &entry2.back();
  } else {
    long now = entry2.size();
    entry2.push_back(SpaceParticle());
    entry2.push_back(SpaceParticle());
    daughter = &entry2[now-1];
    sister = &entry2[now];
    mother = &entry2[now+1];
    opposite = &entry1.back();
  }

  // Define mother and sister flavours, status and properties.
  long idDaughter = daughter->id();
  long idMother = daughter->idMother;
  long idSister = daughter->idSister;
  mother->id(idMother);
  mother->status(-1);
  mother->canBranch = true;
  mother->hasBranched = false;
  mother->copyOfPrimary = (idMother == idDaughter) 
    ? daughter->copyOfPrimary : -1L;
  mother->coneSide = 0;
  sister->id(idSister);
  sister->status(1);
  sister->canBranch = false;
  sister->m(Mass(idSister));
  sister->coneSide = 0;

  // Define colour flow in branching ...
  ++maxColIndx;

  // ... for q -> q g, qbar -> qbar g.
  if (idSister == 21) {
    if (idDaughter > 0 && idDaughter < 9) {
      mother->cols(maxColIndx, 0);
      sister->cols(maxColIndx, daughter->col());
    } else if (idDaughter < 0 && idDaughter > -9) {
      mother->cols(0, maxColIndx);
      sister->cols(daughter->anticol(), maxColIndx);

  // ... for g -> g g.
    } else if (Rndm::flat() > 0.5) {
      mother->cols(daughter->col(), maxColIndx);
      sister->cols(daughter->anticol(), maxColIndx);
    } else {
      mother->cols(maxColIndx, daughter->anticol());
      sister->cols(maxColIndx, daughter->col());
    } 
     
  // ... for q -> g q, qbar -> g qbar.
  } else if (idDaughter == 21) {          
    if (idMother > 0) {
      mother->cols(daughter->col(), 0);
      sister->cols(daughter->anticol(), 0);
    } else {
      mother->cols(0, daughter->anticol());
      sister->cols(0, daughter->col());
    }

  // ... for g -> q qbar, g -> qbar q.
  } else if (idMother == 21) {          
    if (idDaughter > 0) {
      mother->cols(daughter->col(), maxColIndx);
      sister->cols(0, maxColIndx);
    } else {
      mother->cols(maxColIndx, daughter->anticol());
      sister->cols(maxColIndx, 0);
    }

  // ... for f/fbar -> f/fbar gamma. 
  } else if (idSister == 22) {
    mother->cols(daughter->col(), daughter->anticol());
    sister->cols(0, 0);

  // ... for q -> gamma q, qbar -> gamma qbar, l/lbar -> gamma l/lbar.
  } else if (idDaughter == 22) {
    if (idMother > 0 && idMother < 9) {
      mother->cols(maxColIndx, 0);
    } else if (idMother < 0 && idMother > -9) {
      mother->cols(0, maxColIndx);
    } else {
      mother->cols(0, 0);
    }
    daughter->cols(mother->col(), mother->anticol()); 

  // ... for gamma -> f fbar, gamma -> fbar f.
  } else if (idMother == 22) {
    mother->cols(0, 0);
    sister->cols(daughter->anticol(), daughter->col());
  }

  // Set x and maximal virtuality for mother and sister evolution.
  mother->x = daughter->x / daughter->zMother;
  mother->Q2Now = daughter->Q2Now; 
  if (RESPECTSCALE == 1 && mother->scale() < mother->Q2Now) 
    mother->Q2Now = mother->scale();
  mother->sin2thetaMax = daughter->sin2thetaNow;
  sister->Q2Now = daughter->Q2Now;  
  mother->scale(daughter->scale());
  sister->scale(daughter->scale());

  // Set maximum virtuality for mother evolution when no Q2 ordering.
  if (Q2ORDER >=1) {
    double Q21 = daughter->Q2Now;
    double Q22 = opposite->Q2Now;
    double z = daughter->zMother;
    double Q2Max = 0.5 * (1./z + 1.) * Q21 + 0.5 * (1./z - 1.) *
      (Q22 - sHat + sqrt( pow(sHat + Q22 + Q21, 2) +
      8. * Q22 * Q21 * z / (1. -z) ));
    mother->Q2Now = Q2Max;
    mother->Q2Now = Q21 / z;    
  }
 
  // Update status code of now unresolved daughter.
  daughter->addstatus(-1);

}

//*********

// Check and set kinematics for branching after mother evolution;
// especially pick non-isotropic azimuthal angle.

bool SpaceShower::kinemBranching(long side) {

  // Identify mother, daughter, sister, and opposite-side recoiling.
  SpaceParticle* daughter = (side == 1) ? &entry1[entry1.size()-3] 
    : &entry2[entry2.size()-3];  
  SpaceParticle* sister = (side == 1) ? &entry1[entry1.size()-2] 
    : &entry2[entry2.size()-2];  
  SpaceParticle* mother = (side == 1) ? &entry1.back() : &entry2.back();  
  SpaceParticle* opposite = (side == 1) ? &entry2.back() : &entry1.back();  

  // Kinematical quantities, following ZPC 32 (1986) 67 conventions.
  double Q21 = -daughter->m2Now();
  double Q22 = -opposite->m2Now();
  double Q23 = -mother->m2Now();
  double m24 = sister->m2();
  double z = daughter->zMother;
  double s1 = sHat + Q22 + Q21;
  double s3 = sHat/z + Q22 + Q23;
  double r1 = sqrtpos(s1*s1 - 4. * Q22 * Q21);
  double r3 = sqrtpos(s3*s3 - 4. * Q22 * Q23);

  // Construct pT2 of branching; two cases by physics and/or numerics. 
  double pT2;
  if (opposite->hasBranched && s1 - r1 > TINYKINPREC * s1) {
    double m24Max = 0.5 * (s1 * s3 - r1 * r3) / Q22 - Q21 - Q23;
    pT2 = (m24Max - m24) * (0.5 * (s1 * s3 + r1 * r3)  
      - Q22 * (Q21 + Q23 + m24)) / (r1*r1); 
  } else {
    double m24Max = (Q21/z - Q23) * (sHat/(sHat + Q21) 
      - sHat/(sHat/z + Q23));   
    pT2 = (m24Max - m24) * (sHat/z + Q23) / (sHat + Q21); 
  }
  if (pT2 < 0.) return false;

  // Begin branching kinematics construction in 1+2 rest frame.
  double pT = sqrt(pT2);
  double e3 = 0.5 * (sHat/z + Q22 - Q21 - m24) / mHat; 
  double pz3 = (0.5 * s3 - opposite->e() * e3) / daughter->pz();
  mother->m(mother->mNow());

  // Azimuthal asymmetry coefficient from interference with final colours.
  bool haveConeAsym = false;
  double asymCone = 0.;
  Vec4 coneP;
  if (FINALCONE >= 2 && daughter->coneSide >=1 && sister->id() == 21) {
    Vec4 sisterP(pT, 0., pz3 - daughter->pz(), e3 - daughter->e());
    coneP = (daughter->coneSide == 1) ? daughter->colVec4 
      : daughter->antiColVec4; 
    double thetaGD = theta(sisterP, daughter->p());
    double thetaCD = theta(coneP, daughter->p());
    asymCone = min( 0.95, (thetaGD/thetaCD) * (1. - thetaCD/M_PI));
    if (asymCone > 0.0001) haveConeAsym = true;
  } 

  // Reweight relative phi angle by rejection and renewed kinematics.
  double phi, wtPhi;
  do {
    wtPhi = 1.;

    // Random azimuthal angle. 
    phi = 2. * M_PI * Rndm::flat();

    // Finishranching kinematics construction in 1+2 rest frame.
    double pX = pT * cos(phi);
    double pY = pT * sin(phi);
    mother->p(pX, pY, pz3, e3); 
    sister->p(pX, pY, pz3 - daughter->pz(), e3 - daughter->e());

    // Azimuthal asymmetry weight from interference with initial colours.
    if (haveConeAsym) {
      double cosPhiRel = cosphi(sister->p(), coneP, daughter->p());
      wtPhi *= (1. - asymCone) * (1. - asymCone * cosPhiRel) 
        / (1. + asymCone*asymCone - 2. * asymCone * cosPhiRel);      
    }
  
    // Appropriate rejection of phi angles.
  } while (wtPhi < Rndm::flat());   

  // Boost and rotate whole shower to new rest frame of new originators.
  RotBstMatrix rotBstLocal;
  if (side == 1) rotBstLocal.toCMframe(mother->p(), opposite->p());   
  else rotBstLocal.toCMframe(opposite->p(), mother->p());
  for (long i = 0; i < long(entry1.size()); ++i) {
    entry1[i].rotbst(rotBstLocal);
  }
  for (long i = 0; i < long(entry2.size()); ++i) {
    entry2[i].rotbst(rotBstLocal);
  }

  // Update sHat;
  sHat = sHat / z;
  mHat = sqrtpos(sHat);
  return true;

}

//*********

// Evolve a single parton; main physics routine.
// The daughter is evolved backwards to become unresolved into mother.

void SpaceShower::evolveParton(long side) {

  // Read particle from appropriate side. 
  BeamParticle* beam = (side == 1) ? beam1 : beam2;
  SpaceParticle* daughter = (side == 1) ? &entry1.back() : &entry2.back();  
  // SpaceParticle* opposite = (side == 1) ? &entry2.back() : &entry1.back();  
  // LEIF? UNUSED?

  // Beam nature. Daughter flavour. Reset evolution flag.
  long idBeam = beam->id();
  bool isLeptonBeam = beam->isLepton();
  bool isHadronBeam = beam->isHadron();
  long idDaughter = daughter->id();
  daughter->hasBranched = false;

  // Check whether QCD and QED emissions allowed, separately.
  bool isGluon = (idDaughter == 21) ? true : false;
  long intColour = iColour(idDaughter);
  bool isColoured = (intColour > 1) ? true : false;
  bool isPhoton = (idDaughter == 22) ? true : false;
  long intCharge = iCharge(idDaughter);
  double charge = intCharge/3.;  
  bool isCharged = (intCharge != 0 || isPhoton) ? true : false;
  if (isHadronBeam && HADRONSHOWER <= 1) isCharged = false;
  else if (isLeptonBeam && LEPTONSHOWER <=1 && idDaughter != idBeam 
    && !isPhoton ) isColoured = isCharged = false;
  bool allowPhotonBranch = (isHadronBeam && HADRONSHOWER >= 3) 
    || (isLeptonBeam && LEPTONSHOWER >= 2);
  bool allowLeptonBranch = (isHadronBeam && HADRONSHOWER >= 3) 
    || (isLeptonBeam && LEPTONSHOWER >= 1);

  // Set Q0 scale based on the allowed branchings
  double Q0Now;
  if (isColoured && isCharged) Q0Now = min(Q0, Q0CHGQ);
  else if (isColoured) Q0Now = Q0;
  else if (isCharged) Q0Now = Q0CHGL;
  else return; 

  // Check whether to do Matrix Elements corrections.
  bool doMEcorr = (hasME && ( (side == 1 && entry1.size() == 1)
    || (side == 2 && entry2.size() == 1) )) ? true : false;

  // x and Q2 starting values for Q2 evolution. 
  double xDaughter = daughter->x;
  double Q2 = daughter->Q2Now;
 
  // Extrachecks for heavy flavours that evolution at all possible.
  double Q2Start = Q2;
  bool heavyDaughter = ( (abs(idDaughter) == 4 || abs(idDaughter) == 5
    || abs(idDaughter) == 13 || abs(idDaughter) == 15) 
    && !beam->isValence(idDaughter) ) ? true : false;
  double m2Daughter = pow(Mass(idDaughter), 2);
  if (heavyDaughter && (Q2 < m2Daughter || xDaughter > HEAVYEVOL * Q2 
    / (Q2 + m2Daughter) )) return;
 
  // Minimum x step and thereby allowed z range. Finished if no range.
  double xCutCol = max( EMINEMITTED * 2. * mHat / s, 
    xDaughter * (1. / (1. - ZMINEMITTED) - 1.) ); 
  double xCutChg = XMINEMITTEDCHG;
  if (isColoured && xDaughter >= 1. - 2. * xCutCol) return; 
  else if (isCharged && xDaughter >= 1. - 2. * xCutChg) return; 
  double zMinCol = xDaughter / (1. - xCutCol);
  double zMaxCol = xDaughter / (xDaughter + xCutCol);
  double zMinChg = xDaughter / (1. -  xCutChg);
  double zMaxChg = xDaughter / (xDaughter + xCutChg);

  // Integrals of splitting kernels ... 
  double g2gInt=0.0, q2gInt=0.0, q2qInt=0.0, g2qInt=0.0, q2qIntChg=0.0,
    q2gamInt=0.0, l2lInt1=0.0, l2lInt2=0.0;
  double gam2qInt = 0.; 
  double l2gamInt = 0.;
  double l2lInt = 0.;
  double gam2lInt = 0.;

  // ... for gluons: g -> g, q -> g.
  if (isGluon && isColoured) {
    g2gInt = 6. * log(zMaxCol * (1.-zMinCol) / (zMinCol * (1.-zMaxCol)));
    if (doMEcorr) g2gInt *= calcMEmax(21, 21);
    q2gInt = (16./3.) * (1./sqrt(zMinCol) - 1./sqrt(zMaxCol));
    if (doMEcorr) q2gInt *= calcMEmax(1, 21);

  // .. QCD for quarks: q -> q, g -> q.
  } else if (isColoured) {
    q2qInt = (8./3.) * log( (1. - zMinCol) / (1. - zMaxCol) );
    if (doMEcorr) q2qInt *= calcMEmax(1, 1);
    g2qInt = 0.5 * (zMaxCol - zMinCol);
    if (doMEcorr) g2qInt *= calcMEmax(21, 1);

  // ... QED for quarks: q -> q, gamma -> q (not VMD!). 
    if (isCharged) {
      q2qIntChg = 2. * charge*charge * log( (1. - zMinCol) / (1. - zMaxCol) );
      if (doMEcorr) q2qIntChg *= calcMEmax(1, 1);
      if (allowPhotonBranch) gam2qInt = 3. * charge*charge* (zMaxCol - zMinCol);
      if (doMEcorr) gam2qInt *= calcMEmax(22, 1);
    }

  // ... for photons: q -> gamma, l -> gamma. 
  } else if (isPhoton && isCharged) {
    if(allowLeptonBranch) l2gamInt = (1. / xDaughter) 
      * log( (zMaxChg - xDaughter) / (zMinChg - xDaughter) );  
    if (doMEcorr) l2gamInt *= calcMEmax(11, 22);
    q2gamInt =  4. * (1./sqrt(zMinCol) - 1./sqrt(zMaxCol));
    if (doMEcorr) q2gamInt *= calcMEmax(1, 22);

  // ... for leptons: l -> l, gamma -> l. 
  } else if (isCharged) {
    l2lInt1 = log( (1. - zMinChg) / (1. - zMaxChg) ); 
    l2lInt2 = log( (zMaxChg - xDaughter) / (zMinChg - xDaughter) );  
    if (allowLeptonBranch) l2lInt = 2. * (l2lInt1 + l2lInt2);
    if (doMEcorr) l2lInt *= calcMEmax(11, 11);
    if (allowPhotonBranch) gam2lInt = zMaxChg - zMinChg;
    if (doMEcorr) gam2lInt *= calcMEmax(22, 11);
  }

  // Optional: set up maximal cone of emission for primaries.
  double thetaInitMax = 100.;
  if (daughter->coneSide == 1) thetaInitMax = 
    theta(daughter->p(), daughter->colVec4);
  else if (daughter->coneSide == 2) thetaInitMax = 
    theta(daughter->p(), daughter->antiColVec4);

  // Some constants and variables used inside Q2 evolution loop.
  double LAMBDA2 = LAMBDA5 * LAMBDA5; 
  double alphaLogMin = log(Q0Now*Q0Now / LAMBDA2);
  double m2ChgEvol = pow( max( TINYQCHG, Mass(idBeam) ), 2);
  double evolChgExtra = 1.;
  if (isLeptonBeam) evolChgExtra = EXTRAQEDPREWT / log(s/m2ChgEvol);
  long idMother=0, idSister=0;
  double z=0.0, wtZQ=0.0, xfDaughter=0.0, xfMother[21], xfMotherSum=0.0,
    xfGMother=0.0, kernelPDFcol=0.0, evolCol=0.0,
    // xfLMother, //LEIF? UNUSED?
    xfGamMother=0.0, kernelPDFchg=0.0, evolChg=0.0; 

  // Begin loop over evolution in Q2.
  bool isOKbranch=false;
  bool needNewPDF = true;
  bool freeQ2 = true;
  long loop = 0; 
  do {
    isOKbranch = false;
    ++loop;

    if ( loop > 10000 ) {
      //      cerr << "Pythia7::SpaceShower caught in infinite loop.\n";
      // Silently give up here...
      break;
    }

    // First loop, or when loop went to end, new PDF's should be evaluated.
    if (needNewPDF) {
      xfDaughter = max(TINYPDF, beam->xfx(idDaughter, xDaughter, Q2));
      
      // QCD Evolution coefficient by sum of AP kernels * PDF weights.
      if (isColoured) {
        if (isGluon) {
	  xfMotherSum = 0.;
          for (long i = -5; i <= 5; ++i) {
	    if (i == 0) {
	      xfMother[10] = 0.;
            } else {
              xfMother[i+10] = beam->xfx(i, xDaughter, Q2); 
              xfMotherSum += xfMother[i+10]; 
            }
          } 
          kernelPDFcol = g2gInt + q2gInt * xfMotherSum / xfDaughter;
        } else {
          xfGMother = beam->xfx(21, xDaughter, Q2);
          kernelPDFcol = q2qInt + g2qInt * xfGMother / xfDaughter;
        }  

        // Evolution coefficient depends on alpha_s running scheme.
        if (ALPHASMODE <= 0) evolCol = 2. * M_PI 
          / (ALPHASFIX * max(TINYKERNELPDF, kernelPDFcol) );
        else if (ALPHASMODE == 1) evolCol = 23. 
          / (6. * max(TINYKERNELPDF, kernelPDFcol) ); 
        else evolCol = 23. * alphaLogMin 
          / (6. * max(TINYKERNELPDF, kernelPDFcol) );
      }

      // QED Evolution coefficient by sum of AP kernels * PDF weights.
      if (isCharged) {
        if (isColoured) {
          xfGamMother = beam->xfx(22, xDaughter, Q2);
          kernelPDFchg = q2qIntChg + gam2qInt * xfGamMother / xfDaughter;
        } else if (isPhoton) {
          kernelPDFchg = 0.;
          for (long i = -8; i <= 8; ++i) {
	    if (abs(i) >= 6) {
	      idMother = (i > 0) ? 2 * i - 1 : 2 * i + 1;
              xfMother[i+10] = beam->xfx(idMother, xDaughter, Q2);
              kernelPDFchg += l2gamInt * xfMother[i+10];
	    } else if (i != 0) {
	      double qCharge2 = (abs(i)%2 == 1) ? 1./9. : 4./9.;
              xfMother[i+10] = beam->xfx(i, xDaughter, Q2);
              kernelPDFchg += q2gamInt * qCharge2 * xfMother[i+10];
            }
	  }
	  kernelPDFchg /= xfDaughter; 
        } else {
          xfGamMother = beam->xfx(22, xDaughter, Q2);
          kernelPDFchg = l2lInt + gam2lInt * xfGamMother / xfDaughter;
        }  
        evolChg = 2. * M_PI / (ALPHAEMFIX *  max(TINYKERNELPDF, 
          kernelPDFchg) ); 
      }

      // End evaluation new PDF's and evolution coefficients.
      needNewPDF = false;
    }

    // Pick a Q2 for QCD shower, if relevant.
    double Q2col = Q2;
    if (isColoured) {
      // Evolution with fixed alpha_s (or upper estimate for pT^2 scale).
      if (ALPHASMODE != 1) Q2col *= pow(Rndm::flat(), evolCol); 
      // Evolution with running alpha_s(Q2).
      else Q2col = LAMBDA2 * pow(Q2col / LAMBDA2, pow(Rndm::flat(), evolCol));
      if (kernelPDFcol < TINYKERNELPDF) Q2col = 0.;
    }

    // Pick a Q2 for QED shower, if relevant.
    double Q2chg = Q2;
    double wtChgExtra = 1.;
    // QED evolution in hadron.
    if (isCharged && isHadronBeam) {
      //      Q2chg *= pow(Rndm::flat(), evolChg);
      Q2chg *= exp(max(-100.0, evolChg*log(Rndm::flat())));
      if (kernelPDFchg < TINYKERNELPDF) Q2chg = 0.;
    // QED evolution in lepton: preweighting towards small Q2.
    } else if (isCharged) {
      Q2chg = m2ChgEvol * pow( Q2chg / m2ChgEvol, pow(Rndm::flat(), 
        evolChg * evolChgExtra) );
      if (kernelPDFchg < TINYKERNELPDF) Q2chg = 0.;
      wtChgExtra = evolChgExtra * log(Q2chg / m2ChgEvol);
    }

    // Pick between QCD and QED branching where required.
    // May need to reassign if one evolution below cut while other above.
    bool pickedQCD = (isColoured) ? true : false;
    if (isColoured && isCharged) {
      if (Q2chg > Q2col) pickedQCD = false;
      if (pickedQCD && Q0 > Q0CHGQ && Q2col < Q0*Q0) pickedQCD = false;
      if (!pickedQCD && Q0CHGQ > Q0 && Q2chg < Q0CHGQ*Q0CHGQ) 
        pickedQCD = true;
    } 
    
    // Set new Q2 value accordingly; for heavy flavours sometimes forced.
    if (freeQ2) Q2 = (pickedQCD) ? Q2col : Q2chg; 
    if (heavyDaughter && Q2 < m2Daughter) Q2 = exp(0.05 * min(loop,50L)) 
      * m2Daughter; 
    freeQ2 = true;

    // Evolution finished if below cut-off. 
    if (pickedQCD && Q2 < Q0*Q0) return;
    if (!pickedQCD && isColoured && Q2 < Q0CHGQ*Q0CHGQ) return;  
    if (!pickedQCD && !isColoured &&  Q2 < Q0CHGL*Q0CHGL) return;

    // Select z value (and flavour) of branching, and corrective weight ...
    // ... for gluon daughter.
    if (isGluon) {
      // g -> g (+ g). 
      if (g2gInt > Rndm::flat() * kernelPDFcol) {
        idMother = 21;
        idSister = 21;
        z = 1. / ( 1. + ((1. -zMinCol) / zMinCol) * pow( (zMinCol * 
          (1. - zMaxCol)) / (zMaxCol * (1. - zMinCol)), Rndm::flat() ) );
        wtZQ = pow( 1. - z * (1. - z), 2);
      } else {
      // q -> g (+ q): also select flavour. 
        double temp = xfMotherSum * Rndm::flat();
        idMother = -6;
        do { temp -= xfMother[(++idMother) + 10]; } 
        while (temp > 0 && idMother < 5);  
        idSister = idMother;
        z = (zMinCol * zMaxCol) / pow( sqrt(zMinCol) + Rndm::flat() 
          * ( sqrt(zMaxCol)- sqrt(zMinCol) ), 2);
        wtZQ = 0.5 * (1. + pow(1. - z, 2)) * sqrt(z) 
          * xfDaughter / xfMother[idMother + 10];
      } 

    // ... for quark daughter with QCD branching.
    } else if (pickedQCD) {
      // q -> q (+ g). 
      if (q2qInt > Rndm::flat() * kernelPDFcol) {
        idMother = idDaughter;
        idSister = 21;
        z = 1. - (1. - zMinCol) * pow( (1. - zMaxCol) / (1. - zMinCol),
          Rndm::flat() ); 
        wtZQ = 0.5 * (1. + z*z);
      // g -> q (+ qbar). 
      } else {
        idMother = 21;
        idSister = - idDaughter; 
        z = zMinCol + Rndm::flat() * (zMaxCol - zMinCol);
	wtZQ = (z*z + pow(1. - z, 2)) * xfDaughter / xfGMother ;
      }

    // ... for quark daughter with QED branching.
    } else if (isColoured) {
      // q -> q (+ gamma). 
      if (q2qIntChg > Rndm::flat() * kernelPDFchg) {
        idMother = idDaughter;
        idSister = 22;
        z = 1. - (1. - zMinCol) * pow( (1. - zMaxCol) / (1. - zMinCol),
          Rndm::flat() ); 
        wtZQ = 0.5 * (1. + z*z);
      // gamma -> q (+ qbar) (not VMD!).
      } else {
        idMother = 22;
        idSister = -idDaughter;
        z = zMinCol + Rndm::flat() * (zMaxCol - zMinCol);
        wtZQ = (z*z + pow(1. - z, 2)) * xfDaughter / xfGamMother; 
      }
        
    // ... for photon daughter; from quark or lepton, so select flavour.
    } else if (isPhoton) {
      double temp = Rndm::flat() * kernelPDFchg;
      long i = -9;
      do {
	++i;
	if (abs(i) >= 6) {
          idMother = (i > 0) ? 2 * i - 1 : 2 * i + 1;
          temp -= l2gamInt * xfMother[i+10];
        } else if (i != 0) {
	  idMother = i;
	  double qCharge2 = (abs(i)%2 == 1) ? 1./9. : 4./9.;
          temp -= q2gamInt * qCharge2 * xfMother[i+10];
        }
      } while (temp > 0 && i < 8);  
      idSister = idMother;
      // Selection of z different for quark and lepton.
      if (abs(idMother) < 10) {   
        z = (zMinCol * zMaxCol) / pow( sqrt(zMinCol) + Rndm::flat() 
          * ( sqrt(zMaxCol)- sqrt(zMinCol) ), 2);
        wtZQ = 0.5 * (1. + pow(1. - z, 2)) * sqrt(z); 
      } else {
        z = xDaughter + (zMinChg - xDaughter) * pow( (zMaxChg - xDaughter) 
          / (zMinChg - xDaughter), Rndm::flat() );
        wtZQ = 0.5 * (1. + pow(1. - z, 2)) * xDaughter * (z - xDaughter) / z;
      }
      wtZQ *= xfDaughter / xfMother[i + 10];

    // ... for lepton daughter.
    } else if (l2lInt > kernelPDFchg * Rndm::flat()) {
      // l -> l.
      idMother = idDaughter;
      idSister = 22;
      if (l2lInt1 > Rndm::flat() * (l2lInt1 + l2lInt2)) {
        z = 1. - (1. - zMinChg) * pow( (1. - zMaxChg) / (1. - zMinChg),
          Rndm::flat() ); 
      } else {
        z = xDaughter + (zMinChg - xDaughter) * pow( (zMaxChg - xDaughter) 
          / (zMinChg - xDaughter), Rndm::flat() );
      }
      wtZQ = 0.5 * (1. + z*z) * (z - xDaughter) / (1. - xDaughter); 
    } else {
      // gamma -> l.
      idMother = 22;
      idSister = -idDaughter;
      z = zMinChg + Rndm::flat() * (zMaxChg - zMinChg);
      wtZQ = (z*z + pow(1. - z, 2)) * xfDaughter / xfGamMother; 
    } 

    // Extra correction factor for QED evolution in lepton.
    if (isLeptonBeam && !pickedQCD) wtZQ *= wtChgExtra;
     
    // Option with resummation of soft gluon emission as effective z shift.
    if (pickedQCD && SOFTGLUONRESUM > 0) {
      double deltaSoft = (8./3.) * (1. - zMaxCol) * log(Q2Start/Q2) 
	* 23. / (6. * log(sqrt(Q2 * Q2Start) / LAMBDA2));
      if (isGluon) deltaSoft *= 9./4.;
      double zSoft = 1.;
      do zSoft = 1. + deltaSoft * log(Rndm::flat());
      while (z * zSoft < zMinCol);
      z *= zSoft;
    }

    // Remove kinematically impossible branchings (equivalent to uHat cut).
    double m2Sister = pow(Mass(idSister), 2);
    if ( z > sHat * Q2 / ((sHat + Q2) * (Q2 + m2Sister)) ) continue;

    // Extra requrement when reconstructing back to non-valence heavy quark.
    bool heavyMother = ( (abs(idMother) == 4 || abs(idMother) == 5
      || abs(idMother) == 13 || abs(idMother) == 15)
      && !beam->isValence(idMother) ) ? true : false;
    if ( heavyMother && !heavyDaughter && xDaughter / z > HEAVYXMAX) continue; 

    // Correct to alpha_strong of pT^2. (It would be more time-efficient
    // to reject here, but current solution helps ensure wtZQ < 1.)
    if (pickedQCD && ALPHASMODE >= 2) {
      double pT2 = (1. - z) * Q2;
      if (pT2 < Q0*Q0) continue;
      wtZQ *= alphaLogMin / log(pT2 / LAMBDA2);
    }

    // Restrict maximal cone of emission for primaries (optional).
    if (daughter->coneSide >= 1) {
      double thetaMax = atan( sqrt( 4. * z*z * Q2 / ((1. - z) * sHat) )) 
        + atan( sqrt( 4. * z*z * (1. - z) * Q2 / sHat )); 
      if (thetaMax > thetaInitMax) continue;  
    }

    // Angular ordering requirement in shower itself (complementary to above).
    daughter->sin2thetaNow = 4. * z*z * Q2 / (4. * z*z * Q2 
      + (1. - z) * xDaughter*xDaughter * s); 
    if (ANGULARORDER >= 1 && daughter->sin2thetaNow > daughter->sin2thetaMax)
      continue;

    // Save temporary values. 
    daughter->Q2Now = Q2;
    daughter->idMother = idMother;
    daughter->idSister = idSister;
    daughter->zMother = z;

    // Evaluation of ME correction.
    if (doMEcorr) wtZQ *= calcMEcorr(idMother, idDaughter, z, Q2) 
      / calcMEmax(idMother, idDaughter); 

    // Respect user-set maximum scale for emissions - different alternatives.
    if (RESPECTSCALE == 1 && Q2 > daughter->scale()) wtZQ = 0.;
    if (RESPECTSCALE == 2 && (1. - z) * Q2 > daughter->scale()) wtZQ = 0.;
    if (RESPECTSCALE == 3 && Q2 / (z*z * (1.-z)) > daughter->scale()) wtZQ = 0.;
    if (RESPECTSCALE == 4 && 4. * Q2 / (xDaughter*xDaughter * (1.-z) *s) 
      > daughter->scale()) wtZQ = 0.;

    // Evaluation of new daughter and mother PDF's.
    double xfMotherNew = beam->xfx(idMother, xDaughter/z, Q2);
    double xfDaughterNew = beam->xfx(idDaughter, xDaughter, Q2);
    needNewPDF = true;
    if (xfDaughterNew < TINYPDF) {
      if (idMother == idDaughter) {
        Q2 = Q2Start; 
        q2qInt = g2gInt = l2lInt = q2qIntChg = 0.;
        continue;
      } else if (Q2 < Q2STARTFRAC * Q2Start) {
        Q2 = sqrt(Q2 * Q2Start);
        freeQ2 = false;
        continue; 
      } else xfDaughterNew = TINYPDF;  
    }
    wtZQ *= xfMotherNew / xfDaughterNew;
  
    // Rejection of impossible branchings; else end of Q2 evolution.
    isOKbranch = true;
  } while (!isOKbranch || wtZQ < Rndm::flat() );
  daughter->hasBranched = true;
}

//*********

// Find class of ME correction.

void SpaceShower::findMEkind(long MEkindin) {

  MEkind = 0;
  if (MEkindin > 0) MEkind = MEkindin;
  else if (nFinal == 1) {

    // f + fbar -> vector boson. 
    if ( (finalId1 == 23 || abs(finalId1) == 24 || finalId1 == 32 
      || finalId1 == 33 || abs(finalId1) == 34 || abs(finalId1) == 41)
      && abs(entry1[0].id()) < 20 && abs(entry2[0].id()) < 20)
      MEkind = 1;

    // g + g -> Higgs boson.
    if ( (finalId1 == 25 || finalId1 == 35 || finalId1 == 36)
       && entry1[0].id() == 21 && entry2[0].id() == 21 ) MEkind = 2; 
  }

}

//*********

// Provide maximum of expected ME weight; for preweighting of evolution.

double SpaceShower::calcMEmax(long idMother, long idDaughter) {

  // Currently only one non-unity case, so simplify.
  if (MEkind == 1 && idMother > 20 && abs(idDaughter) < 20 ) return 3.;
  return 1.;

}  

//*********

// Provide actual ME weight for current branching.

double SpaceShower::calcMEcorr(long idMother, long idDaughter, double z, 
  double Q2) {

  // Convert to Mandelstam variables.
  double M2 = sHatFinalState;
  double sH = M2 / z;
  double tH = -Q2;
  double uH = Q2 - M2 * (1. - z) / z;

  // Corrections for f + fbar -> s-channel vector boson.
  if (MEkind == 1) {
    if (abs(idMother) < 20) {
      return (tH*tH + uH*uH + 2. * M2 * sH) / (sH*sH + M2*M2); 
    } else {
      return (sH*sH + uH*uH + 2. * M2 * tH) / (pow(sH - M2, 2) + M2*M2); 
    }

  // Corrections for g + g -> Higgs boson.
  } else if (MEkind == 2) {
    if (abs(idMother) < 20) {
      return (sH*sH + uH*uH) / (sH*sH + pow(sH - M2, 2)); 
    } else {
      return 0.5 * (pow(sH, 4) + pow(tH, 4) + pow(uH, 4) + pow(M2, 4)) 
        / pow(sH*sH - M2 * (sH - M2), 2); 
    }    
  }

  return 1.;

}
