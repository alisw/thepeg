// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the TestMisc class.
//

#include "TestMisc.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/PDT/MatcherBase.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/EventRecord/Event.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Handlers/FlavourGenerator.h"
#include "ThePEG/Cuts/Cuts.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "TestMisc.tcc"
#endif


using namespace ThePEG;

TestMisc::~TestMisc() {}

NoPIOClassDescription<TestMisc> TestMisc::initTestMisc;
// Definition of the static class description member.

#define checkfail(cond)         \
  if ( cond ) { \
    clog << "passed." << endl;  \
  } else {                      \
    ++fail;                     \
    clog << "failed." << endl;  \
  }                             \
  ++tests

void TestMisc::Init() {
  using std::clog;
  breakThePEG();
  if ( CurrentGenerator::isVoid() ) {
    clog << "No EventGenerator found. Could not perform tests." << endl;
    exit(1);
  }
  clog << endl << "Performing miscellaneous tests ..." << endl << endl;

  int fail = 0;
  int tests = 0;
  EventGenerator & eg = CurrentGenerator::current();

  clog << "Checking initialization of MatcherBase objects ...                ";
  PMPtr match = eg.getObject<MatcherBase>("/Defaults/Matchers/MatchLightQuark");
  set<PDPtr> expected;
  set<PDPtr> obtained;
  expected.insert(eg.getParticleData(ParticleID::u));
  expected.insert(eg.getParticleData(ParticleID::d));
  expected.insert(eg.getParticleData(ParticleID::s));
  for ( ParticleMap::const_iterator it = eg.particles().begin();
	it != eg.particles().end(); ++it )
    if ( match && match->matches(*(it->second)) )
      obtained.insert(it->second);
  checkfail( match && expected == obtained );

  clog << "Checking that Cuts::maxS() returns something sensible ...         ";
  CutsPtr cut = eg.getObject<Cuts>("/Pythia7/Handlers/Cuts/EECuts");
  checkfail(cut && cut->maxS(tcPDVector()) > 0.0*GeV2 );

  clog << "Checking that SimpleFlavour handles meson codes correctly ...     ";
  Ptr<FlavourGenerator>::const_pointer flg =
    eg.getObject<FlavourGenerator>
    ("/Defaults/Handlers/Hadronization/SimpleFlavour");
  long sumch = 0;
  if ( flg ) {
    for ( int iq = 1; iq <= 5; ++iq ) {
      tcPDPtr q = eg.getParticleData(iq);
      tcPDPtr qb = eg.getParticleData(-iq);
      for ( int n = 0; n < 1000; ++n ) {
	tcPDPair ret = flg->alwaysGenerateHadron(q);
	sumch +=
	  abs(ret.first->iCharge() + ret.second->iCharge() - q->iCharge());
	ret = flg->alwaysGenerateHadron(qb);
	sumch +=
	  abs(ret.first->iCharge() + ret.second->iCharge() - qb->iCharge());
      }
    }
  }
  checkfail(flg && sumch == 0);

  clog << "Checking generation of Poissonian distribution ...                ";
  double sum1 = 0.0;
  double sum12 = 0.0;
  double sum2 = 0.0;
  double sum22 = 0.0;
  long N = 10000;
  for ( int i = 0; i < N; ++i ) {
    long l = UseRandom::rndPoisson(3.14);
    sum1 += l;
    sum12 += l*l;
    l = UseRandom::rndPoisson(1000.0);
    sum2 += l;
    sum22 += l*l;
  }
  ios::fmtflags saveflags = clog.setf(ios::fixed, ios::floatfield);
  clog << endl << setprecision(3)
       << "         average/3.14:" << setw(6) << sum1/N/3.14
       << ", var/av:" << setw(6) << sum12/sum1 - sum1/N << endl
       << "         average/1000:" << setw(6) << sum2/N/1000.0
       << ", var/av:" << setw(6) << sum22/sum2 - sum2/N << endl;
  clog.flags(saveflags);

  if ( fail ) {
    clog << endl << "Error: " << fail << " out of " << tests
	 << " test failed." << endl << endl;
    exit(1);
  } else
    clog << endl << "All tests passed." << endl << endl;
}

