// -*- C++ -*-
#ifndef ThePEG_CPUClock_H
#define ThePEG_CPUClock_H
//
// This is the declaration of the CPUClock class.
//

#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/Utilities/Named.h"
#include "ThePEG/Utilities/DebugItem.h"
#include <ctime>

namespace ThePEG {

/**
 * The CPUClock class can be used to measure the time spent in a
 * certain function. Somewhere a long lasting object of the CPUClock
 * class should be created giving a name as argument. Typically this
 * is done my a static variable in the function of interest. Then in
 * the beginning of the function, a CPUTimer object is created with
 * the CPUClock as argument, and the time during which the CPUTimer
 * object is alive will be accumulated. Note that the inclusive time
 * spent in the function is measured - also calls to other functions -
 * is counted. If time spent in called functions is not wanted, the
 * timer can be stopped and stated by creating a PauseCPUTimer
 * object. If the function calls itself recursively, the time will not
 * be double-counted. At any time the current values of all CPUClock
 * objects can be dumped to a file with the static dump() function.
 *
 * Note that the call to access the CPU clock takes some time, so the
 * CPUClock class should not be used in very short functions. This
 * time affect the system CPU usage, but this is incluced in the time
 * reported by the CPU clock.
 *
 * For this reason, the cals to the CPU clock is only switched on if
 * the DebugItem ThePEG::CPUTimer is switched on (debug level 2). If
 * not, only the number of calls to the function is recorded.
 */
class CPUClock: public Named {

public:

  /**
   * Everybody needs a friend.
   */
  friend class CPUTimer;

  /**
   * Or two.
   */
  friend class PauseCPUTimer;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The only relevant constructor. The string should typically be on
   * the form "namspace::class::function".
   */
  CPUClock(string clockname)
    : Named(clockname), theNCalls(0), startClock(0), totTime(0), loop(0) {
    clocks().push_back(this);
  }
  //@}

private:

  /**
   * Call this when entring a function, starting the clock.
   */
  void enter() {
    ++theNCalls;
    start();
  }

  /**
   * Start the clock.
   */
  void start() {
    if ( !loop++ ) startClock = clock();
  }

  /**
   * Stop the clock.
   */
  void stop() {
    if ( !--loop ) totTime += clock() - startClock;
  }

public:

  static clock_t clock() {
    static DebugItem on("ThePEG::CPUTimer", 2);
    return on? std::clock(): 0;
  }

  int nCalls() const {
    return theNCalls;
  }

  double timeTot() const {
    return double(totTime)/double(CLOCKS_PER_SEC);
  }

  double timePer() const {
    return timeTot()/double(max(nCalls(), 1));
  }

  static void dump(ostream & os) {
    os << "CPUTimer results:" << endl;
    for ( int i = 0, N = clocks().size(); i < N; ++i )
      os << "  " << clocks()[i]->name() << ": " << clocks()[i]->nCalls()
	   << " calls totalling " << clocks()[i]->timeTot()
	   << "s (" << clocks()[i]->timePer() << "s/call)" << endl;
  }

private:

  static vector<const CPUClock *> & clocks() {
    static vector<const CPUClock *> theClocks;
    return theClocks;
  }

  unsigned long theNCalls;
  // The number of calls made to the corresponding timer.

  clock_t startClock;
  // The clock at the last call to start.

  unsigned long totTime;
  // The total time spent in the corresponding timer.

  unsigned int loop;
  // The loop level

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  CPUClock & operator=(const CPUClock &);

  /**
   * The copy constructor is private and must never be called.
   * In fact, it should not even be implemented.
   */
  CPUClock(const CPUClock &);

  /**
   * The default constructor is private and must never be called.
   * In fact, it should not even be implemented.
   */
  CPUClock();

};

class CPUTimer {

public:

  /**
   * The only relevant constructor starts the given CPUClock.
   */
  CPUTimer(CPUClock & c): theClock(&c) {
    c.enter();
  }

  /**
   * The destructor stops the clock.
   */
  ~CPUTimer() {
    theClock->stop();
  }

private:

  /**
   * The clock under controll.
   */
  CPUClock * theClock;

  /**
   * No default constructor
   */
  CPUTimer();

  /**
   * No copy constructor;
   */
  CPUTimer(const CPUTimer &);

  /**
   * No assignment
   */
  CPUTimer & operator=(const CPUTimer &);

};

class PauseCPUTimer {

public:

  /**
   * The only relevant constructor starts the given CPUClock.
   */
  PauseCPUTimer(CPUClock & c): theClock(&c) {
    c.stop();
  }

  /**
   * The destructor stops the clock.
   */
  ~PauseCPUTimer() {
    theClock->start();
  }

private:

  /**
   * The clock under controll.
   */
  CPUClock * theClock;

  /**
   * No default constructor
   */
  PauseCPUTimer();

  /**
   * No copy constructor;
   */
  PauseCPUTimer(const PauseCPUTimer &);

  /**
   * No assignment
   */
  PauseCPUTimer & operator=(const PauseCPUTimer &);

};

}

#endif /* ThePEG_CPUClock_H */
