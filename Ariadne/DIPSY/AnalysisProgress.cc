// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AnalysisProgress class.
//

#include "AnalysisProgress.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/Current.h"
#include "DipoleXSec.h"
#include "DipoleEventHandler.h"


#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

#include <sys/times.h>
#include <unistd.h>

using namespace DIPSY;

AnalysisProgress::AnalysisProgress(): secstep(0) {}

AnalysisProgress::~AnalysisProgress() {}

void AnalysisProgress::initialize() {
  ieve = 0;
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

void AnalysisProgress::
analyze(const vector<DipoleStatePtr> & vr, const vector<DipoleStatePtr> & vl,
	const vector<ImpactParameters> & vb, const DipoleXSec & xsec,
	const Vec3D & probs, double jac) {
  long n = Current<DipoleEventHandler>()->preSamples();
  long i = ++ieve;
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

double AnalysisProgress::fclock() {
  struct tms tmsbuf;
  times(&tmsbuf);
  double d =
    tmsbuf.tms_utime+tmsbuf.tms_stime+tmsbuf.tms_cutime+tmsbuf.tms_cstime;
  d /= sysconf(_SC_CLK_TCK);
  return d;
}

bool AnalysisProgress::statusTime(long i, long n) const {
  if ( i <= 0 ) return false;
  if ( i == n ) return true;
  if ( i > n/2 ) i = n-i;
  while ( i >= 10 && !(i%10) ) i /= 10;
  if ( i == 1 || i == 2 || i == 5 ) return true;
  if ( secstep > 0 && time(0) > time1 + secstep ) return true;
  return false;
}

IBPtr AnalysisProgress::clone() const {
  return new_ptr(*this);
}

IBPtr AnalysisProgress::fullclone() const {
  return new_ptr(*this);
}


void AnalysisProgress::persistentOutput(PersistentOStream & os) const {
  os << secstep;
}

void AnalysisProgress::persistentInput(PersistentIStream & is, int) {
  is >> secstep;
}


// Static variable needed for the type description system in ThePEG.
#include "ThePEG/Utilities/DescribeClass.h"
DescribeClass<AnalysisProgress,DIPSY::DipoleAnalysisHandler>
  describeDIPSYAnalysisProgress("DIPSY::AnalysisProgress", "AnalysisProgress.so");


void AnalysisProgress::Init() {

  static ClassDocumentation<AnalysisProgress> documentation
    ("AnalysisProgress writes out the number of events used in the log file.");

  static Parameter<AnalysisProgress,int> interfaceInterval
    ("Interval",
     "Besides the standard intervals, also write a status line every "
     "given number of seconds.",
     &AnalysisProgress::secstep, 0, 0, 0,
     true, false, Interface::lowerlim);

  interfaceInterval.setHasDefault(false);

}

