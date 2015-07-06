// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ProgressLog class.
//

#include "ProgressLog.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include <sys/times.h>
#include <unistd.h>

using namespace ThePEG;

ProgressLog::ProgressLog(): secstep(0) {}

ProgressLog::~ProgressLog() {}

void ProgressLog::analyze(tEventPtr event, long ieve, int loop, int state) {
  AnalysisHandler::analyze(event, ieve, loop, state);
  long n = generator()->N();
  long i = ieve;
  if ( !statusTime(i, n) ) return;

  double fcpui = fclock();
  time_t timei = time(0);
  double ftime0 = time0;
  double ftime1 = time1;
  double ftimei = timei;
  double eff = 1.0;
  if ( ftimei > ftime1 && fcpui > fcpu1 )
    eff = (fcpui-fcpu1)/(ftimei-ftime1);
  if ( eff >= 1.0 ) eff = 0.999999;
  int ieff = 100*eff;
  double eff0 = 1.0;
  if ( ftimei > ftime0 && fcpui > fcpu0 )
    eff0 = (fcpui-fcpu0)/(ftimei-ftime0);
  if ( eff0 >= 1.0 ) eff0 = 0.999999;
  int ieff0 = 100*eff0;
  double fcpun = fcpu0 + (n*(fcpui-fcpu0))/i;
  time_t timen = (time_t)(ftimei + (fcpun-fcpui)/eff + 30.0);
  time_t timen0 = (time_t)(ftimei + (fcpun-fcpui)/eff0 + 30.0);
  char date[1024];
  char daten[1024];
  char daten0[1024];
  strftime(date,1024,"%y.%m.%d %H:%M",localtime(&timei));
  strftime(daten,1024,"%H:%M",localtime(&timen));
  strftime(daten0,1024,"%H:%M",localtime(&timen0));
  long ii = i;
  if ( n - i < n/10 ) ii = i - n;
  time_t dayn = (timen - timei)/86400;
  time_t dayn0 = (timen0 - timei)/86400;

  ostream & os = generator()->log();

  if ( dayn <= 0 && dayn0 <= 0 ) {
    os << date << " " << setw(8) << ii << "/" << setw(9);
    os.setf(ios::left, ios::adjustfield);
    os << n << " etc:   " << daten << "[";
    os.setf(ios::right, ios::adjustfield);
    os << setw(2) << ieff << "%]   " << daten0 << "[" << ieff0 << "%] "
       << host << ":" << pid << endl << flush;
  } else {
    os << date << " " << setw(8) << ii << "/" << setw(9);
    os.setf(ios::left, ios::adjustfield);
    os << n << " etc: " << dayn << "+" << daten << "[";
    os.setf(ios::right, ios::adjustfield);
    os << setw(2) << ieff << "%] "
       << dayn0 << "+" << daten0 << "[" << ieff0 << "%] "
       << host << ":" << pid << endl << flush;
  }
    
  fcpu1 = fcpui;
  time1 = timei;
  
}

double ProgressLog::fclock() {
  struct tms tmsbuf;
  times(&tmsbuf);
  double d =
    tmsbuf.tms_utime+tmsbuf.tms_stime+tmsbuf.tms_cutime+tmsbuf.tms_cstime;
  d /= sysconf(_SC_CLK_TCK);
  return d;
}

bool ProgressLog::statusTime(long i, long n) const {
  if ( i <= 0 ) return false;
  if ( i == n ) return true;
  if ( i > n/2 ) i = n-i;
  while ( i >= 10 && !(i%10) ) i /= 10;
  if ( i == 1 || i == 2 || i == 5 ) return true;
  if ( secstep > 0 && time(0) > time1 + secstep ) return true;
  return false;
}

IBPtr ProgressLog::clone() const {
  return new_ptr(*this);
}

IBPtr ProgressLog::fullclone() const {
  return new_ptr(*this);
}

void ProgressLog::doinitrun() {
  AnalysisHandler::doinitrun();
  fcpu0 = fcpu1 = fclock();
  time0 = time1 = time(0);
  char name[1024];
  gethostname(name,1024);
  host = name;
  if ( host.find(".") != string::npos ) host = host.substr(0, host.find("."));
  pid = getpid();
  char date[1024];
  strftime(date,1024,"%y.%m.%d %H:%M",localtime(&time0));
  ostream & os = generator()->log();
  os << date << "        0/" << setw(9);
  os.setf(ios::left, ios::adjustfield);
  os << generator()->N();
  os.setf(ios::right, ios::adjustfield);
  os << " Initializing...                "
     << host << ":" << pid << endl << flush;
}


void ProgressLog::persistentOutput(PersistentOStream & os) const {
  os << secstep;
}

void ProgressLog::persistentInput(PersistentIStream & is, int) {
  is >> secstep;
}

ClassDescription<ProgressLog> ProgressLog::initProgressLog;
// Definition of the static class description member.

void ProgressLog::Init() {

  static ClassDocumentation<ProgressLog> documentation
    ("The ProgressLog class will not perform an actual analysis. Instead"
     " it will write out a progress status on the standard log file. By"
     " default it will write on event 1, 2, 5, 10, 20, 50, ... etc. But"
     " optionally it can in addition also write out every given number of"
     " seconds.\n\n"
     "The status line which is written out contains the current date "
     "and time, the number of events processed so far and the total number"
     "of events to be generated, two estimates of the time of completion "
     "(one based on the current cpu usage and one based on the average "
     "cpu usage [the usage is given in brackets]), and the host on which "
     "the program is running, together with its process number.");


  static Parameter<ProgressLog,int> interfaceInterval
    ("Interval",
     "Besides the standard intervals, also write a status line every "
     "given number of seconds.",
     &ProgressLog::secstep, 0, 0, 0,
     true, false, Interface::lowerlim);
  interfaceInterval.setHasDefault(false);


}

