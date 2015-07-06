// This file collects top-level classes for shower evolution. 

#ifndef Shower_H
#define Shower_H

#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>
#include "Pythia7/Config/Pythia7.h"

namespace Pythia7 {
namespace Shower {

using namespace std;

class TimeShower;
class SpaceShower;

//**************************************************************************


/**
 * This class hold a partonic event, which is to shower (on input)
 * or has showered (on output).
 * Used by the internal Pythia7 Shower classes.
 */
class Event {
    
public:
  /**
   * Constructor.
   */
  Event(long capacity = 100) {entry.reserve(capacity); maxColIndx = 0;}


  /**
   * Input of new particle at the end of the event record.
   */
  long append(Particle entryin) {    
    entry.push_back(entryin); 
    if (entryin.col() > maxColIndx) maxColIndx = entryin.col();   
    if (entryin.anticol() > maxColIndx) maxColIndx = entryin.anticol();
    return entry.size() - 1;
  }


  /**
   * Removal of particle.
   */
  void remove(long index) {entry.erase(entry.begin() + index);}


  /**
   * Overload index operator to access element of Particle vector.
   */
  Particle& operator[](long i) {return entry[i];}
  /**
   * Overload index operator to access element of Particle vector.
   */
  const Particle & operator[](long i) const {return entry[i];}
  /**
   * Find if given entry is daughter/granddaughter/... of another one.
   */
  bool isDescendent(long, long) const;

  /**
   * Query or reset size.
   */
  long size() const {return entry.size();}
  /**
   * Query or reset size.
   */
  void zero() {entry.resize(0); maxColIndx = 0;}
  /**
   * Manipulate colour index.
   */
  void colIndx(long indx) {maxColIndx = indx;}
  /**
   * Manipulate colour index.
   */
  long colIndx() const {return maxColIndx;}


  /**
   * Print an event.
   */
  friend ostream& operator<<(ostream&, const Event&) ;  


private: 
  /** NOT DOCUMENTED */
  vector<Particle> entry;
  /** NOT DOCUMENTED */
  long maxColIndx;
};
 
//**************************************************************************


/**
 * The Shower class administrates initial- and final-state radiation,
 * by one call to SpaceShower and repeated calls to TimeShower.
 * (Currently only latter implemented.)
 * Used by the internal Pythia7 Shower classes.
 */
class Shower {

public:
  /**
   * Constant (possibly to be changed externally).
   */
  static long ISR;

  /**
   * Constant (possibly to be changed externally).
   */
  static long FSR;

  /**
   * Constant (possibly to be changed externally).
   */
  static long FSRONISR;

  /**
   * Pointer to the object doing the actual work for time-like showers.
   */
  static TimeShower * theTimeShower;

  /**
   * Pointer to the object doing the actual work for space-like showers.
   */
  static SpaceShower * theSpaceShower;

  /**
   * Driver routine to handle a succession of showers.
   */
  void shower(Event&, BeamParticle&, BeamParticle&, bool = false);
  /**
   * Driver routine to handle a succession of showers.
   * Simpler version for final-state radiation only; no beams required.
   */
  void shower(Event&); 


private: 

  /** @cond NEVERTOBEDOCUMENTED */

  /** NOT DOCUMENTED */
  long in1Old, in2Old, in1New, in2New, sizeOld,sizeMid, sizeNew;
  /** NOT DOCUMENTED */
  bool hadISR;
  /** NOT DOCUMENTED */
  void doSpaceShower(Event&, BeamParticle&, BeamParticle&);
  /** NOT DOCUMENTED */
  void doTimeOnSpaceShower(Event&);
  /** NOT DOCUMENTED */
  void copyFinalState(Event&);
  /** NOT DOCUMENTED */
  void doTimeShower(Event&);
  /** @endcond */

};
 
//**************************************************************************

}
}

#endif // Shower_H
