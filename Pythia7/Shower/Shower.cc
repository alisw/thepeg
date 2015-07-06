// Function definitions not found in the header.

#include "Basics.h"
#include "Beam.h"
#include "Shower.h"
#include "TimeShower.h"
#include "SpaceShower.h"

using namespace Pythia7::Shower;

//**************************************************************************

// Find if particle is daughter of another particle.
bool Pythia7::Shower::Event::isDescendent(long daughter, long mother) const {
  long search = daughter;    
  while (entry[search].mother1() > mother) search = entry[search].mother1();
  if (entry[search].mother1() == mother) return true;
  else return false;
}

//*********

// Print an event.

ostream& Pythia7::Shower::operator<<(ostream& os, const Event& event) {
  os << "\n         Event Listing  \n \n  i     id st  m1  m2  p  c"
     << " ac    p_x       p_y       p_z        e         m \n";
  Vec4 pSum;
  for (long i = 0; i < long(event.entry.size()); ++i) 
    {os.width(3); os << i << event.entry[i];
    if (event.entry[i].status() == 1) pSum += event.entry[i].p(); 
  }
  os << "                              " << pSum;
  return os;
}
 
//**************************************************************************

// Master switch for initial-state radiation; = 0: off, = 1: on.
long Shower::ISR = 1;

// Master switch for final-state radiation; = 0: off, = 1: on.
long Shower::FSR = 1;

// Master switch for final-state radiation off initial-state cascade; 
// = 0: off, = 1: on.
long Shower::FSRONISR = 1;

TimeShower *  Shower::theTimeShower= 0;
SpaceShower * Shower::theSpaceShower = 0;

//*********

// Driver routine to handle a combination of initial- and final-state showers.

void Shower::shower(Event& event, BeamParticle& beam1, BeamParticle& beam2,
  bool isInCMframe) {

  // Boost event to rest frame with beams along +-z axis.
  RotBstMatrix rotBstCM; 
  if (!isInCMframe) {
    rotBstCM.toCMframe(beam1.p(), beam2.p());
    for (long i = 0; i < event.size(); ++i) {event[i].rotbst(rotBstCM);}  
  }

  // Initial-state shower, including default values if no such shower..
  sizeNew = 0;
  hadISR = false;
  if (ISR > 0 && beam1.id() != 0 && beam2.id() != 0) 
    doSpaceShower(event, beam1, beam2);
 
  // Final-state showers on timelike branches of initial-state shower.
  if (hadISR && FSRONISR > 0) doTimeOnSpaceShower(event);

  // New copy of final state, with boost from initial-state shower.
  if ( hadISR
       || sizeOld < event.size()) //LEIF? WHY DO WE NEED THIS?
    copyFinalState(event);

  // Final-state showers, sequentially for resonance decays.
  if (FSR > 0) doTimeShower(event);

  // Boost back event to original frame.
  if (!isInCMframe) {
    rotBstCM.invert();
    for (long i = 0; i < event.size(); ++i) {event[i].rotbst(rotBstCM);} 
  } 

}

//*********

// Simpler version of above driver routine to handle only final-state showers.

void Shower::shower(Event& event) {

  // Final-state showers, sequentially for resonance decays.
  sizeNew = 0;
  if (FSR > 0) doTimeShower(event);

}

//*********

// Routine to handle the initial-state shower.

void Shower::doSpaceShower(Event& event, BeamParticle& beam1, 
  BeamParticle& beam2) {

  // Find positions of incoming partons.
  sizeOld = event.size(); 
  in1Old = -1;
  in2Old = -1;
  for (long i = 0; i < event.size(); ++i) {   
    if (event[i].status() == -1) {
      if (in1Old == -1) in1Old = i;
      else {in2Old = i; break;}
    }
  }
  if (in2Old == -1) return; 
   
  // Send away system for spacelike shower.
  theSpaceShower->shower(event, beam1, beam2, in1Old, in2Old);
  sizeMid = event.size();
  hadISR = (sizeMid - sizeOld > 2) ? true : false;  

  // Final ISR values in case no timelike showers on side branches.
  sizeNew = sizeMid;
  in1New = sizeNew - 2;
  in2New = sizeNew - 1;
  
}

//*********

// Rotine to handle final-state showers on side branches of the 
// initial-state shower. Final-state radiation of the hard scattering 
// is handled separately below.

void Shower::doTimeOnSpaceShower(Event& event) {

  // Cumulative boosts when side branches acquire a mass.
  RotBstMatrix cumulBst;

  // Momentum of shower initiators, to deduct timelike branches from.
  Vec4 pSide1, pSide2;
  long motherIn1 = event[in1Old].mother1();
  long motherIn2 = event[in2Old].mother1();
  for (long i = sizeOld; i < sizeMid; ++i) {
    if (event[i].mother1() == motherIn1) pSide1 = event[i].p();
    if (event[i].mother1() == motherIn2) pSide2 = event[i].p();
  }
  // What happens when beams not in event record???

  // Detect and copy side branch parton that may shower shower.
  for (long i = sizeOld; i < sizeMid - 2; ++i) {
    if (event[i].status() > 0) {
      long iCopy = event.append(event[i]);   
      event[iCopy].rotbst(cumulBst);
      event[i].addstatus(3);

      // Set maximum timelike virtuality by spacelike virtuality of sister.
      long iMother = i - 1;
      long iSister = iMother;
      for (long j = i + 1; j < sizeMid; ++j) {               
	if (event[j].mother1() == iMother) iSister = j;
      }
      double Q2max = event[iSister].m2();

      // Define recoling pseudoparticle by system inside considered emission.
      Vec4 pRec;
      for (long j = i + 1; j < sizeMid; ++j) {               
	if (j >= sizeMid - 2 || event[j].status() > 0) pRec += event[j].p();
      }
      pRec.rotbst(cumulBst);
      long iRec = event.append(Particle(0, -3, i+1, sizeMid-1, 0, 0, pRec, 
        pRec.mcalc()));

      // Send off system to shower. 
      vector<long> primary;
      primary.push_back(iCopy);
      primary.push_back(iRec);
      theTimeShower->shower(event, primary, Q2max);
      Vec4 pRecMoved = event[iRec+2].p();
      cumulBst.bst(pRec, pRecMoved);

      // Remove unshowered side branch copy, and two pseudoparticle copies.
      event.remove(iRec+2),
      event.remove(iRec),
      event.remove(iCopy),
      
      // Set modified mother pointers.
      event[iCopy].mothers(i, -1);
      for (long j = iCopy + 1; j < event.size(); ++j) {
	long mother = event[j].mother1();
        if (mother == iCopy + 2) event[j].mothers(mother - 2, -1);
	else event[j].mothers(mother - 3, -1); 
      }

      // Update momenta of incoming spacelike partons.
      if (event.isDescendent(iCopy, motherIn1)) pSide1 -= event[iCopy].p(); 
      if (event.isDescendent(iCopy, motherIn2)) pSide2 -= event[iCopy].p(); 
    }
  }

  // Copy incoming spacelike partons to hard scattering.
  // Record their positions and update their mother pointers.
  in1New = event.append(event[sizeMid - 2]); 
  event[in1New].mothers(sizeMid - 2, -1);
  event[in1New].p(pSide1);
  event[in1New].m(pSide1.mcalc());
  in2New = event.append(event[sizeMid - 1]); 
  event[in2New].mothers(sizeMid - 1, -1);
  event[in2New].p(pSide2);
  event[in2New].m(pSide2.mcalc());
  sizeNew = event.size();

}

//*********

// New copy of final state, with boost from initial-state shower.

void Shower::copyFinalState(Event& event) {

  // Copy final state, with new mothers.
  long shift = sizeNew - in2Old - 1;
  for (long iOld = in2Old + 1; iOld < sizeOld; ++iOld) {
    long iNew = event.append(event[iOld]); 
    long mother = event[iOld].mother1();
    if (mother == in1Old) event[iNew].mother1(in1New);
    else if (mother == in2Old) event[iNew].mother1(in2New);
    else if (mother > in2Old) event[iNew].mother1(mother + shift);
    mother = event[iOld].mother2();
    if (mother == in1Old) event[iNew].mother2(in1New);
    else if (mother == in2Old) event[iNew].mother2(in2New);
    else if (mother > in2Old) event[iNew].mother2(mother + shift);
    event[iNew].prev(iOld);
    event[iOld].status(4);
  }
  
  // Boost and rotation induced by shower imposed on final state.
  RotBstMatrix newCM;
  newCM.toCMframe(event[in1Old].p(), event[in2Old].p());
  newCM.fromCMframe(event[in1New].p(), event[in2New].p());
  for (long iNew = sizeNew; iNew < event.size(); ++iNew) {
    event[iNew].rotbst(newCM);
  }  

}

//*********

// Struct to help keep list of final-state partons that should shower.

struct ShoweringParticle {
  long iBefore, iAfter, mother, status; 
  RotBstMatrix showerBst;
  ShoweringParticle(long iBeforein, long motherin, long statusin) {
  iBefore = iBeforein; iAfter = -1; mother = motherin; status = statusin;} 
}; 

//*********

// Routine to handle the final-state showers, 
// in several steps for sequential resonance decays.

void Shower::doTimeShower(Event& event) {

  // Pick up all partons that could/should be allowed to shower.
  vector<ShoweringParticle> parton;  
  for (long i = sizeNew; i < event.size(); ++i) {
    long status = event[i].status();
    if (status >= 1 && status <= 3) parton.push_back(ShoweringParticle(i, 
      event[i].mother1(), status));
  } 
  if (parton.size() <= 1) return; 

  // Loop over all partons to find systems that may branch.
  long last = 0; 
  long first;
  do {
    first = last;
    vector<long> primary;
    while (last < long(parton.size()) &&
	   parton[last].mother == parton[first].mother) {
      primary.push_back(parton[last].iBefore);
      ++last;
    }

    // Check if common mother has been moved, i.e. has showered earlier.
    long iMove = -1;
    if (first > 1) {
      for (long j = 0; j < first; ++j) {
	if (parton[j].iBefore == parton[first].mother) iMove = j;
      }
    }

    // If so update colour information of system to be treated.
    if (iMove >= 0) {
      long oldCol = event[parton[iMove].iBefore].col();
      long newCol = event[parton[iMove].iAfter].col();
      long oldAntiCol = event[parton[iMove].iBefore].anticol();
      long newAntiCol = event[parton[iMove].iAfter].anticol(); 
      for (long i = first; i < last; ++i) { 
	long ii = parton[i].iBefore;
	if (event[ii].col() == oldCol) event[ii].col(newCol);
	if (event[ii].anticol() == oldAntiCol) event[ii].anticol(newAntiCol);
      }
    }

    // Send away completed system for timelike shower.
    if (last - first == 1) continue; 
    long oldSize = event.size(); 
    theTimeShower->shower(event, primary);
    long newSize = event.size();

    // If common mother has been boosted then also boost showered system.
    if (iMove >= 0) {
      for (long i = oldSize; i < newSize; ++i) {
        event[i].rotbst(parton[iMove].showerBst);
      }
    }

    // Loop in new system: locate showered copy of original decayed parton.
    for (long i = first; i < last; ++i) {
      if (parton[i].status > 1) { 
        Particle& before = event[parton[i].iBefore]; 
        for (long j = oldSize; j < newSize; ++j) {
	  if (event[j].id() == before.id() && event[j].status() <= 3
          && event.isDescendent(j, parton[i].iBefore)) parton[i].iAfter = j; 
        }
        Particle& after = event[parton[i].iAfter]; 

        // Calculate boost of this parton by shower.
        parton[i].showerBst.bst(before.p(), after.p());
      }
    }

  // Keep on looking for more partons that may shower.
  } while (last < long(parton.size()));  

}
