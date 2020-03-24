// -*- C++ -*-
//
// utilitiesTestSmearing.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad, 2015 Marco A. Harrendorf
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_Repository_Test_RandomGenerator_H
#define ThePEG_Repository_Test_RandomGenerator_H

#include <boost/test/unit_test.hpp>
#include <boost/iterator/iterator_concepts.hpp>

#include "ThePEG/Repository/RandomGenerator.h"
#include <StandardRandom.h>



/*
 * Helper class to test generated random numbers
 * 
 * The class checks how many random numbers are inside a pre-defined interval and how many outside
 * This class is e.g. useful to check random numbers produced by a gaussian distribution
 */
template <typename Unit>
class BinningCheck {
  private:
    // Do not accept missing interval limits
    BinningCheck() {}
  public:
    // Define interval limits
    BinningCheck(Unit lo, Unit hi) 
      : m_lo(lo), m_hi(hi), m_in(0), m_out(0) {}
    ~BinningCheck() {}
    
    int in() {return m_in;}
    int out() {return m_out;}
    int tot() {return m_out + m_in;}
    // Check if random number is inside / outside interval and increase corresponding counter
    void add(Unit randomNumber) {
      if(randomNumber >= m_lo && randomNumber <= m_hi) {
	m_in++;
      } else {
	m_out++;
      }
    }
    // Reset counters
    void reset() {m_in = 0; m_out = 0;}
  
  private:
    Unit m_lo, m_hi;
    int m_in, m_out; 
};
typedef BinningCheck<double> DoubleBinCheck;
typedef BinningCheck<long> LongBinCheck;

/*
 * Start of BOOST unit tests for Helper class
 * 
 */
BOOST_AUTO_TEST_SUITE(HelperRandomNumberBinning)

BOOST_AUTO_TEST_CASE(HelperDoubleBinning)
{
  DoubleBinCheck a(0, 1);
  a.add(-1.1);
  a.add(-0.1);
  a.add(0.1);
  a.add(0.5);
  a.add(0.8);
  a.add(1.1);
  a.add(100);
  BOOST_CHECK_EQUAL(a.tot(), 7);
  BOOST_CHECK_EQUAL(a.in(), 3);
  BOOST_CHECK_EQUAL(a.out(), 4);
  
  a.reset();
  BOOST_CHECK_EQUAL(a.tot(), 0);
  BOOST_CHECK_EQUAL(a.in(), 0);
  BOOST_CHECK_EQUAL(a.out(), 0);
  
  
  DoubleBinCheck b(-1.5, 0.5);
  b.add(-1.1);
  b.add(-0.1);
  b.add(0.1);
  b.add(0.5);
  b.add(0.8);
  b.add(1.1);
  b.add(100);
  BOOST_CHECK_EQUAL(b.tot(), 7);
  BOOST_CHECK_EQUAL(b.in(), 4);
  BOOST_CHECK_EQUAL(b.out(), 3);
}

/* 
 * End of BOOST unit tests for Helper class
 * 
 */
BOOST_AUTO_TEST_SUITE_END()



/*
 * Local fix to provide randomGenerator object
 * 
 */
struct FixLocal1 {
  FixLocal1() : rng() {
    BOOST_TEST_MESSAGE( "setup local fixture for repositoryTestRandomGenerator" ); 
  }
  
  ~FixLocal1()  { 
    BOOST_TEST_MESSAGE( "teardown local fixture for repositoryTestRandomGenerator" ); 
  }
  
  ThePEG::StandardRandom rng;
};

/*
 * Start of boost unit tests for RandomGenerator.h
 * 
 */
BOOST_FIXTURE_TEST_SUITE(repositoryRandomGenerator, FixLocal1)

/*
 * Boost unit tests 
 * 
 */
BOOST_AUTO_TEST_CASE(rndZeroToOne)
{
  int N = 1000;
  // Check for whole interval 
  DoubleBinCheck posRange(0, 1);
  // Check for flat distribution
  DoubleBinCheck quarter1(0, 0.25);
  DoubleBinCheck quarter3(0.5, 0.75);
  
  // rnd function
  for(int i = 0; i < N; ++i) {
    posRange.add(rng.rnd());
    quarter1.add(rng.rnd());
    quarter3.add(rng.rnd());
  }
  // Prob laying inside of interval should be 1
  BOOST_CHECK_EQUAL(posRange.tot(), N);
  BOOST_CHECK_EQUAL(posRange.in(), N);
  BOOST_CHECK_EQUAL(posRange.out(), 0);
  
  // Prob laying inside of interval should be 0.25
  BOOST_CHECK_EQUAL(quarter1.tot(), N);
  BOOST_CHECK_CLOSE(quarter1.in(), 0.25 * N, 10);
  BOOST_CHECK_CLOSE(quarter1.out(), 0.75 * N, 5);
  
  // Prob laying inside of interval should be 0.25
  BOOST_CHECK_EQUAL(quarter3.tot(), N);
  BOOST_CHECK_CLOSE(quarter3.in(), 0.25 * N, 10);
  BOOST_CHECK_CLOSE(quarter3.out(), 0.75 * N, 5);
  
  
  // repeat for operator()
  posRange.reset();
  quarter1.reset();
  quarter3.reset();
  for(int i = 0; i < N; ++i) {
    posRange.add(rng());
    quarter1.add(rng());
    quarter3.add(rng());
  }
  BOOST_CHECK_EQUAL(posRange.tot(), N);
  BOOST_CHECK_EQUAL(posRange.in(), N);
  BOOST_CHECK_EQUAL(posRange.out(), 0);
  
    // Prob laying inside of interval should be 0.25
  BOOST_CHECK_EQUAL(quarter1.tot(), N);
  BOOST_CHECK_CLOSE(quarter1.in(), 0.25 * N, 10);
  BOOST_CHECK_CLOSE(quarter1.out(), 0.75 * N, 5);
  
  // Prob laying inside of interval should be 0.25
  BOOST_CHECK_EQUAL(quarter3.tot(), N);
  BOOST_CHECK_CLOSE(quarter3.in(), 0.25 * N, 10);
  BOOST_CHECK_CLOSE(quarter3.out(), 0.75 * N, 5);
}

BOOST_AUTO_TEST_CASE(rndZeroTohi)
{
  int N = 1000;
  double hi = 2.5;
    // Check for whole interval 
  DoubleBinCheck zero2hi(0, 2.5);
  // Check for flat distribution
  DoubleBinCheck quarter1(0, 0.5);
  DoubleBinCheck quarter3(1.5, 2.0);
  
  // rnd function
  for(int i = 0; i < N; ++i) {
    zero2hi.add(rng.rnd(hi));
    quarter1.add(rng.rnd(hi));
    quarter3.add(rng.rnd(hi));
  }
  // Prob laying inside of interval should be 1
  BOOST_CHECK_EQUAL(zero2hi.tot(), N);
  BOOST_CHECK_EQUAL(zero2hi.in(), N);
  BOOST_CHECK_EQUAL(zero2hi.out(), 0);
  
  // Prob laying inside of interval should be 0.2
  BOOST_CHECK_EQUAL(quarter1.tot(), N);
  BOOST_CHECK_CLOSE(quarter1.in(), 0.2 * N, 10);
  BOOST_CHECK_CLOSE(quarter1.out(), 0.8 * N, 5);
  
  // Prob laying inside of interval should be 0.2
  BOOST_CHECK_EQUAL(quarter3.tot(), N);
  BOOST_CHECK_CLOSE(quarter3.in(), 0.2 * N, 10);
  BOOST_CHECK_CLOSE(quarter3.out(), 0.8 * N, 5);
  
  
  // repeat for operator(), note it is generally requiring a long!
  zero2hi.reset();
  quarter1.reset();
  quarter3.reset();
  for(int i = 0; i < N; ++i) {
    zero2hi.add(rng(hi));
    quarter1.add(rng(hi));
    quarter3.add(rng(hi));
  }
  BOOST_CHECK_EQUAL(zero2hi.tot(), N);
  BOOST_CHECK_EQUAL(zero2hi.in(), N);
  BOOST_CHECK_EQUAL(zero2hi.out(), 0);
  
    // Prob laying inside of interval should be 0.2
  BOOST_CHECK_EQUAL(quarter1.tot(), N);
  BOOST_CHECK_CLOSE(quarter1.in(), 0.2 * N, 10);
  BOOST_CHECK_CLOSE(quarter1.out(), 0.8 * N, 5);
  
  // Prob laying inside of interval should be 0.2
  BOOST_CHECK_EQUAL(quarter3.tot(), N);
  BOOST_CHECK_CLOSE(quarter3.in(), 0.2 * N, 10);
  BOOST_CHECK_CLOSE(quarter3.out(), 0.8 * N, 5);
  
  // repeat for operator(), note it is requiring a long!
  long hiTwo = 3;
  // Check for whole interval 
  LongBinCheck zero2hiTwo(0, hiTwo);
  // Check for flat distribution
  LongBinCheck quarter1Two(1, 1);
  LongBinCheck quarter3Two(2, 2);
  
  // operator()(long N)
  for(int i = 0; i < N; ++i) {
    zero2hiTwo.add(rng(hiTwo));
    quarter1Two.add(rng(hiTwo));
    quarter3Two.add(rng(hiTwo));
  }
  // Prob laying inside of interval should be 1
  BOOST_CHECK_EQUAL(zero2hiTwo.tot(), N);
  BOOST_CHECK_EQUAL(zero2hiTwo.in(), N);
  BOOST_CHECK_EQUAL(zero2hiTwo.out(), 0);
  
  // Prob laying inside of interval should be 0.333
  BOOST_CHECK_EQUAL(quarter1Two.tot(), N);
  BOOST_CHECK_CLOSE(quarter1Two.in(), 0.333 * N, 10);
  BOOST_CHECK_CLOSE(quarter1Two.out(), 0.666 * N, 5);
  
  // Prob laying inside of interval should be 0.333
  BOOST_CHECK_EQUAL(quarter3Two.tot(), N);
  BOOST_CHECK_CLOSE(quarter3Two.in(), 0.333 * N, 10);
  BOOST_CHECK_CLOSE(quarter3Two.out(), 0.666 * N, 5);  
}

BOOST_AUTO_TEST_CASE(rndIntervallLowerLimitTohi)
{
  int N = 1000;
  double lo = -1.5;
  double hi = 2.5;
  // Check for whole interval 
  DoubleBinCheck zero2hi(lo, hi);
  // Check for flat distribution
  DoubleBinCheck quarter1(-0.5, 0.5);
  DoubleBinCheck quarter3(1.5, 2.5);
  
  // rnd function
  for(int i = 0; i < N; ++i) {
    zero2hi.add(rng.rnd(lo, hi));
    quarter1.add(rng.rnd(lo, hi));
    quarter3.add(rng.rnd(lo, hi));
  }
  // Prob laying inside of interval should be 1
  BOOST_CHECK_EQUAL(zero2hi.tot(), N);
  BOOST_CHECK_EQUAL(zero2hi.in(), N);
  BOOST_CHECK_EQUAL(zero2hi.out(), 0);
  
  // Prob laying inside of interval should be 0.25
  BOOST_CHECK_EQUAL(quarter1.tot(), N);
  BOOST_CHECK_CLOSE(quarter1.in(), 0.25 * N, 10);
  BOOST_CHECK_CLOSE(quarter1.out(), 0.75 * N, 5);
  
  // Prob laying inside of interval should be 0.25
  BOOST_CHECK_EQUAL(quarter3.tot(), N);
  BOOST_CHECK_CLOSE(quarter3.in(), 0.25 * N, 10);
  BOOST_CHECK_CLOSE(quarter3.out(), 0.75 * N, 5);
}

BOOST_AUTO_TEST_CASE(rndZeroToOneVector)
{
  int N = 10;
  int L = 10;
    // Check for whole interval 
  DoubleBinCheck posRange(0, 1);
  // Check for flat distribution
  DoubleBinCheck quarter1(0, 0.25);
  DoubleBinCheck quarter3(0.5, 0.75);
  
  for(int i = 0; i < N; ++i) {
    posRange.reset();
    quarter1.reset();
    quarter3.reset();
    
    std::vector<double> rndvec = rng.rndvec(L);
    BOOST_CHECK_EQUAL(static_cast<int>(rndvec.size()), L);
    for(int j = 0; j < static_cast<int>(rndvec.size()); ++j) {
      posRange.add(rndvec[j]);
      quarter1.add(rndvec[j]);
      quarter3.add(rndvec[j]);
    }
    // Prob laying inside of interval should be 1
    BOOST_CHECK_EQUAL(posRange.tot(), N);
    BOOST_CHECK_EQUAL(posRange.in(), N);
    BOOST_CHECK_EQUAL(posRange.out(), 0);
    
    // Prob laying inside of interval should be 0.25
    BOOST_CHECK_EQUAL(quarter1.tot(), N);
    BOOST_CHECK_CLOSE(quarter1.in(), 0.25 * L, 300);
    BOOST_CHECK_CLOSE(quarter1.out(), 0.75 * L, 300);
    
    // Prob laying inside of interval should be 0.25
    BOOST_CHECK_EQUAL(quarter3.tot(), N);
    BOOST_CHECK_CLOSE(quarter3.in(), 0.25 * L, 300);
    BOOST_CHECK_CLOSE(quarter3.out(), 0.75 * L, 300);
  }
}

BOOST_AUTO_TEST_CASE(rndBoolSingleProbability)
{
  int N = 1000;
  // Check for whole interval 
  LongBinCheck posRange(0, 1);
  // Check for flat distribution
  LongBinCheck zero(0, 0);
  LongBinCheck one(1, 1);
  
  // rndbool function, prob should be 0.5
  for(int i = 0; i < N; ++i) {
    posRange.add(rng.rndbool());
    zero.add(rng.rndbool());
    one.add(rng.rndbool());
  }
  // Prob laying inside of interval should be 1
  BOOST_CHECK_EQUAL(posRange.tot(), N);
  BOOST_CHECK_EQUAL(posRange.in(), N);
  BOOST_CHECK_EQUAL(posRange.out(), 0);
  
  // Prob laying inside of interval should be 0.5
  BOOST_CHECK_EQUAL(zero.tot(), N);
  BOOST_CHECK_CLOSE(zero.in(), 0.5 * N, 10);
  BOOST_CHECK_CLOSE(zero.out(), 0.5 * N, 10);
  
  // Prob laying inside of interval should be 0.5
  BOOST_CHECK_EQUAL(one.tot(), N);
  BOOST_CHECK_CLOSE(one.in(), 0.5 * N, 10);
  BOOST_CHECK_CLOSE(one.out(), 0.5 * N, 10);
  
  
  // Check for whole interval 
  LongBinCheck posRangeTwo(0, 1);
  // Check for flat distribution
  LongBinCheck midPoint(0, 0);
  LongBinCheck hiEdge(1, 1);
  
  // rndbool function, prob should be now 0.4
  for(int i = 0; i < N; ++i) {
    posRangeTwo.add(rng.rndbool(0.4));
    midPoint.add(rng.rndbool(0.4));
    hiEdge.add(rng.rndbool(0.4));
  }
  // Prob laying inside of interval should be 1
  BOOST_CHECK_EQUAL(posRangeTwo.tot(), N);
  BOOST_CHECK_EQUAL(posRangeTwo.in(), N);
  BOOST_CHECK_EQUAL(posRangeTwo.out(), 0);
  
  // Prob laying inside of interval should be 0.6
  BOOST_CHECK_EQUAL(midPoint.tot(), N);
  BOOST_CHECK_CLOSE(midPoint.in(), 0.6 * N, 10);
  BOOST_CHECK_CLOSE(midPoint.out(), 0.4 * N, 10);
  
  // Prob laying inside of interval should be 0.4
  BOOST_CHECK_EQUAL(hiEdge.tot(), N);
  BOOST_CHECK_CLOSE(hiEdge.in(), 0.4 * N, 10);
  BOOST_CHECK_CLOSE(hiEdge.out(), 0.6 * N, 10);
}

BOOST_AUTO_TEST_CASE(prndBoolSingleProbability)
{
  int N = 1000;
  // Check for whole interval 
  LongBinCheck posRange(0, 1);
  // Check for flat distribution
  LongBinCheck zero(0, 0);
  LongBinCheck one(1, 1);
  
  // rndbool function, prob should be 0.5
  for(int i = 0; i < N; ++i) {
    posRange.add(rng.prndbool());
    zero.add(rng.prndbool());
    one.add(rng.prndbool());
  }
  // Prob laying inside of interval should be 1
  BOOST_CHECK_EQUAL(posRange.tot(), N);
  BOOST_CHECK_EQUAL(posRange.in(), N);
  BOOST_CHECK_EQUAL(posRange.out(), 0);
  
  // Prob laying inside of interval should be 0.5
  BOOST_CHECK_EQUAL(zero.tot(), N);
  BOOST_CHECK_CLOSE(zero.in(), 0.5 * N, 15);
  BOOST_CHECK_CLOSE(zero.out(), 0.5 * N, 15);
  
  // Prob laying inside of interval should be 0.5
  BOOST_CHECK_EQUAL(one.tot(), N);
  BOOST_CHECK_CLOSE(one.in(), 0.5 * N, 15);
  BOOST_CHECK_CLOSE(one.out(), 0.5 * N, 15);
  
  
  // Check for whole interval 
  LongBinCheck posRangeTwo(0, 1);
  // Check for flat distribution
  LongBinCheck midPoint(0, 0);
  LongBinCheck hiEdge(1, 1);
  
  // rndbool function, prob should be now 0.4
  for(int i = 0; i < N; ++i) {
    posRangeTwo.add(rng.prndbool(0.4));
    midPoint.add(rng.prndbool(0.4));
    hiEdge.add(rng.prndbool(0.4));
  }
  // Prob laying inside of interval should be 1
  BOOST_CHECK_EQUAL(posRangeTwo.tot(), N);
  BOOST_CHECK_EQUAL(posRangeTwo.in(), N);
  BOOST_CHECK_EQUAL(posRangeTwo.out(), 0);
  
  // Prob laying inside of interval should be 0.6
  BOOST_CHECK_EQUAL(midPoint.tot(), N);
  BOOST_CHECK_CLOSE(midPoint.in(), 0.6 * N, 10);
  BOOST_CHECK_CLOSE(midPoint.out(), 0.4 * N, 10);
  
  // Prob laying inside of interval should be 0.4
  BOOST_CHECK_EQUAL(hiEdge.tot(), N);
  BOOST_CHECK_CLOSE(hiEdge.in(), 0.4 * N, 10);
  BOOST_CHECK_CLOSE(hiEdge.out(), 0.6 * N, 10);
}

BOOST_AUTO_TEST_CASE(rndBoolTwoProbabilities)
{
  int N = 1000;
  double p1 = 0.2;
  double p2 = 0.3;
  // Check for whole interval 
  LongBinCheck posRange(0, 1);
  // Check for flat distribution
  LongBinCheck zero(0, 0);
  LongBinCheck one(1, 1);
  
  for(int i = 0; i < N; ++i) {
    posRange.add(rng.rndbool(p1, p2));
    zero.add(rng.rndbool(p1, p2));
    one.add(rng.rndbool(p1, p2));
  }
  // Prob laying inside of interval should be 1
  BOOST_CHECK_EQUAL(posRange.tot(), N);
  BOOST_CHECK_EQUAL(posRange.in(), N);
  BOOST_CHECK_EQUAL(posRange.out(), 0);
  
  // Prob laying inside of interval should be 0.6
  BOOST_CHECK_EQUAL(zero.tot(), N);
  BOOST_CHECK_CLOSE(zero.in(), 0.6 * N, 10);
  BOOST_CHECK_CLOSE(zero.out(), 0.4 * N, 10);
  
  // Prob laying inside of interval should be 0.2/0.5 = 0.4
  BOOST_CHECK_EQUAL(one.tot(), N);
  BOOST_CHECK_CLOSE(one.in(), 0.4 * N, 10);
  BOOST_CHECK_CLOSE(one.out(), 0.6 * N, 10);
}

BOOST_AUTO_TEST_CASE(prndBoolTwoProbabilities)
{
  int N = 1000;
  double p1 = 0.2;
  double p2 = 0.3;
  // Check for whole interval 
  LongBinCheck posRange(0, 1);
  // Check for flat distribution
  LongBinCheck zero(0, 0);
  LongBinCheck one(1, 1);
  
  for(int i = 0; i < N; ++i) {
    posRange.add(rng.prndbool(p1, p2));
    zero.add(rng.prndbool(p1, p2));
    one.add(rng.prndbool(p1, p2));
  }
  // Prob laying inside of interval should be 1
  BOOST_CHECK_EQUAL(posRange.tot(), N);
  BOOST_CHECK_EQUAL(posRange.in(), N);
  BOOST_CHECK_EQUAL(posRange.out(), 0);
  
  // Prob laying inside of interval should be 0.6
  BOOST_CHECK_EQUAL(zero.tot(), N);
  BOOST_CHECK_CLOSE(zero.in(), 0.6 * N, 10);
  BOOST_CHECK_CLOSE(zero.out(), 0.4 * N, 10);
  
  // Prob laying inside of interval should be 0.2/0.5 = 0.4
  BOOST_CHECK_EQUAL(one.tot(), N);
  BOOST_CHECK_CLOSE(one.in(), 0.4 * N, 10);
  BOOST_CHECK_CLOSE(one.out(), 0.6 * N, 10);
}

BOOST_AUTO_TEST_CASE(rndSign)
{
  int N = 1000;
  double p1 = 0.4;
  double p2 = 0.2;
  double p3 = 0.3;
  // Check for whole interval 
  LongBinCheck range(-1, 1);
  // Check for flat distribution
  LongBinCheck lowEdge(-1, -1);
  LongBinCheck midPoint(0, 0);
  LongBinCheck hiEdge(1, 1);
  
  for(int i = 0; i < N; ++i) {
    range.add(rng.rndsign(p1, p2, p3));
    lowEdge.add(rng.rndsign(p1, p2, p3));
    midPoint.add(rng.rndsign(p1, p2, p3));
    hiEdge.add(rng.rndsign(p1, p2, p3));
  }
  // Prob laying inside of interval should be 1
  BOOST_CHECK_EQUAL(range.tot(), N);
  BOOST_CHECK_EQUAL(range.in(), N);
  BOOST_CHECK_EQUAL(range.out(), 0);
  
  // Prob laying inside of interval should be 0.4
  BOOST_CHECK_EQUAL(lowEdge.tot(), N);
  BOOST_CHECK_CLOSE(lowEdge.in(), 0.4/0.9 * N, 10);
  BOOST_CHECK_CLOSE(lowEdge.out(), 0.5/0.9 * N, 10);
  
  // Prob laying inside of interval should be 0.2
  BOOST_CHECK_EQUAL(midPoint.tot(), N);
  BOOST_CHECK_CLOSE(midPoint.in(), 0.2/0.9 * N, 10);
  BOOST_CHECK_CLOSE(midPoint.out(), 0.7/0.9 * N, 10);
  
  // Prob laying inside of interval should be 0.2
  BOOST_CHECK_EQUAL(hiEdge.tot(), N);
  BOOST_CHECK_CLOSE(hiEdge.in(), 0.3/0.9 * N, 10);
  BOOST_CHECK_CLOSE(hiEdge.out(), 0.6/0.9 * N, 10);
}

BOOST_AUTO_TEST_CASE(prndSign)
{
  int N = 1000;
  double p1 = 0.4;
  double p2 = 0.2;
  double p3 = 0.3;
  // Check for whole interval 
  LongBinCheck range(-1, 1);
  // Check for flat distribution
  LongBinCheck lowEdge(-1, -1);
  LongBinCheck midPoint(0, 0);
  LongBinCheck hiEdge(1, 1);
  
  for(int i = 0; i < N; ++i) {
    range.add(rng.prndsign(p1, p2, p3));
    lowEdge.add(rng.prndsign(p1, p2, p3));
    midPoint.add(rng.prndsign(p1, p2, p3));
    hiEdge.add(rng.prndsign(p1, p2, p3));
  }
  // Prob laying inside of interval should be 1
  BOOST_CHECK_EQUAL(range.tot(), N);
  BOOST_CHECK_EQUAL(range.in(), N);
  BOOST_CHECK_EQUAL(range.out(), 0);
  
  // Prob laying inside of interval should be 0.4
  BOOST_CHECK_EQUAL(lowEdge.tot(), N);
  BOOST_CHECK_CLOSE(lowEdge.in(), 0.4/0.9 * N, 10);
  BOOST_CHECK_CLOSE(lowEdge.out(), 0.5/0.9 * N, 10);
  
  // Prob laying inside of interval should be 0.2
  BOOST_CHECK_EQUAL(midPoint.tot(), N);
  BOOST_CHECK_CLOSE(midPoint.in(), 0.2/0.9 * N, 10);
  BOOST_CHECK_CLOSE(midPoint.out(), 0.7/0.9 * N, 10);
  
  // Prob laying inside of interval should be 0.2
  BOOST_CHECK_EQUAL(hiEdge.tot(), N);
  BOOST_CHECK_CLOSE(hiEdge.in(), 0.3/0.9 * N, 10);
  BOOST_CHECK_CLOSE(hiEdge.out(), 0.6/0.9 * N, 10);
}

BOOST_AUTO_TEST_CASE(rnd2)
{
  int N = 1000;
  double p1 = 0.4;
  double p2 = 0.2;
  // Check for whole interval 
  LongBinCheck posRange(0, 1);
  // Check for flat distribution
  LongBinCheck midPoint(0, 0);
  LongBinCheck hiEdge(1, 1);
  
  for(int i = 0; i < N; ++i) {
    posRange.add(rng.rnd2(p1, p2));
    midPoint.add(rng.rnd2(p1, p2));
    hiEdge.add(rng.rnd2(p1, p2));
  }
  // Prob laying inside of interval should be 1
  BOOST_CHECK_EQUAL(posRange.tot(), N);
  BOOST_CHECK_EQUAL(posRange.in(), N);
  BOOST_CHECK_EQUAL(posRange.out(), 0);
  
  // Prob laying inside of interval should be 0.666
  BOOST_CHECK_EQUAL(midPoint.tot(), N);
  BOOST_CHECK_CLOSE(midPoint.in(), 0.666 * N, 10);
  BOOST_CHECK_CLOSE(midPoint.out(), 0.333 * N, 10);
  
  // Prob laying inside of interval should be 0.333
  BOOST_CHECK_EQUAL(hiEdge.tot(), N);
  BOOST_CHECK_CLOSE(hiEdge.in(), 0.333 * N, 10);
  BOOST_CHECK_CLOSE(hiEdge.out(), 0.666 * N, 10);
}

BOOST_AUTO_TEST_CASE(rnd3)
{
  int N = 1000;
  double p1 = 0.4;
  double p2 = 0.2;
  double p3 = 0.3;
  // Check for whole interval 
  LongBinCheck zero2two(0, 2);
  // Check for flat distribution
  LongBinCheck midPoint(0, 0);
  LongBinCheck hiEdge(1, 1);
  LongBinCheck two(2, 2);
  
  for(int i = 0; i < N; ++i) {
    zero2two.add(rng.rnd3(p1, p2, p3));
    midPoint.add(rng.rnd3(p1, p2, p3));
    hiEdge.add(rng.rnd3(p1, p2, p3));
    two.add(rng.rnd3(p1, p2, p3));
  }
  // Prob laying inside of interval should be 1
  BOOST_CHECK_EQUAL(zero2two.tot(), N);
  BOOST_CHECK_EQUAL(zero2two.in(), N);
  BOOST_CHECK_EQUAL(zero2two.out(), 0);
  
  // Prob laying inside of interval should be 0.4 / 0.9
  BOOST_CHECK_EQUAL(midPoint.tot(), N);
  BOOST_CHECK_CLOSE(midPoint.in(), 0.4/0.9 * N, 10);
  BOOST_CHECK_CLOSE(midPoint.out(), 0.5/0.9 * N, 10);
  
  // Prob laying inside of interval should be 0.2 / 0.9
  BOOST_CHECK_EQUAL(hiEdge.tot(), N);
  BOOST_CHECK_CLOSE(hiEdge.in(), 0.2/0.9 * N, 10);
  BOOST_CHECK_CLOSE(hiEdge.out(), 0.7/0.9 * N, 10);
  
  // Prob laying inside of interval should be 0.3 / 0.9
  BOOST_CHECK_EQUAL(two.tot(), N);
  BOOST_CHECK_CLOSE(two.in(), 0.3/0.9 * N, 10);
  BOOST_CHECK_CLOSE(two.out(), 0.6/0.9 * N, 10);
}

BOOST_AUTO_TEST_CASE(rnd4)
{
  int N = 1000;
  double p1 = 0.4;
  double p2 = 0.2;
  double p3 = 0.3;
  double p4 = 0.5;
  // Check for whole interval 
  LongBinCheck bin_03(0, 3);
  // Check for flat distribution
  LongBinCheck bin0(0, 0);
  LongBinCheck bin1(1, 1);
  LongBinCheck bin2(2, 2);
  LongBinCheck bin3(3, 3);
  
  for(int i = 0; i < N; ++i) {
    bin_03.add(rng.rnd4(p1, p2, p3, p4));
    bin0.add(rng.rnd4(p1, p2, p3, p4));
    bin1.add(rng.rnd4(p1, p2, p3, p4));
    bin2.add(rng.rnd4(p1, p2, p3, p4));
    bin3.add(rng.rnd4(p1, p2, p3, p4));
  }
  // Prob laying inside of interval should be 1
  BOOST_CHECK_EQUAL(bin_03.tot(), N);
  BOOST_CHECK_EQUAL(bin_03.in(), N);
  BOOST_CHECK_EQUAL(bin_03.out(), 0);
  
  // Prob laying inside of interval should be 0.4 / 1.4
  BOOST_CHECK_EQUAL(bin0.tot(), N);
  BOOST_CHECK_CLOSE(bin0.in(),  0.4/1.4 * N, 10);
  BOOST_CHECK_CLOSE(bin0.out(), 1.0/1.4 * N, 10);
  
  // Prob laying inside of interval should be 0.2 / 1.4
  BOOST_CHECK_EQUAL(bin1.tot(), N);
  BOOST_CHECK_CLOSE(bin1.in(),  0.2/1.4 * N, 10);
  BOOST_CHECK_CLOSE(bin1.out(), 1.2/1.4 * N, 10);
  
  // Prob laying inside of interval should be 0.3 / 1.4
  BOOST_CHECK_EQUAL(bin2.tot(), N);
  BOOST_CHECK_CLOSE(bin2.in(),  0.3/1.4 * N, 10);
  BOOST_CHECK_CLOSE(bin2.out(), 1.1/1.4 * N, 10);
  
  // Prob laying inside of interval should be 0.5 / 1.4
  BOOST_CHECK_EQUAL(bin3.tot(), N);
  BOOST_CHECK_CLOSE(bin3.in(),  0.5/1.4 * N, 10);
  BOOST_CHECK_CLOSE(bin3.out(), 0.9/1.4 * N, 10);
}

BOOST_AUTO_TEST_CASE(rndExp)
{
  int N = 1000;
  // Check for whole interval 
  DoubleBinCheck allnums(-DBL_MAX, DBL_MAX);
  // Check for median
  DoubleBinCheck firstHalf(-DBL_MAX, log(2));
  DoubleBinCheck secondHalf(log(2), DBL_MAX);
  // Check for increasing probability
  DoubleBinCheck neg1000(-1000, 0);
  DoubleBinCheck pos1000(0, 1000);
  
  for(int i = 0; i < N; ++i) {
    allnums.add(rng.rndExp());
    firstHalf.add(rng.rndExp());
    secondHalf.add(rng.rndExp());
    neg1000.add(rng.rndExp());
    pos1000.add(rng.rndExp());
  }
  // Prob laying inside of interval should be 1
  BOOST_CHECK_EQUAL(allnums.tot(), N);
  BOOST_CHECK_EQUAL(allnums.in(), N);
  BOOST_CHECK_EQUAL(allnums.out(), 0);
  
  // Prob laying inside of interval should be 0.5
  BOOST_CHECK_EQUAL(firstHalf.tot(), N);
  BOOST_CHECK_CLOSE(firstHalf.in(), 0.5 * N, 15);
  BOOST_CHECK_CLOSE(firstHalf.out(), 0.5 * N, 15);
  
    // Prob laying inside of interval should be 0.5
  BOOST_CHECK_EQUAL(secondHalf.tot(), N);
  BOOST_CHECK_CLOSE(secondHalf.in(), 0.5 * N, 15);
  BOOST_CHECK_CLOSE(secondHalf.out(), 0.5 * N, 15);
  
  // Increasing probability
  BOOST_CHECK(neg1000.in() <= pos1000.in());  
}

BOOST_AUTO_TEST_CASE(rndExpMean)
{
  int N = 1000;
  double meanValue = 5;
  // Check for whole interval 
  DoubleBinCheck allnums(-DBL_MAX, DBL_MAX);
  // Check for median
  DoubleBinCheck firstHalf(-DBL_MAX, log(2) * meanValue);
  DoubleBinCheck secondHalf(log(2) * meanValue, DBL_MAX);
  // Check for increasing probability
  DoubleBinCheck neg1000(-1000, 0);
  DoubleBinCheck pos1000(0, 1000);
  
  for(int i = 0; i < N; ++i) {
    allnums.add(rng.rndExp(meanValue));
    firstHalf.add(rng.rndExp(meanValue));
    secondHalf.add(rng.rndExp(meanValue));
    neg1000.add(rng.rndExp(meanValue));
    pos1000.add(rng.rndExp(meanValue));
  }
  // Prob laying inside of interval should be 1
  BOOST_CHECK_EQUAL(allnums.tot(), N);
  BOOST_CHECK_EQUAL(allnums.in(), N);
  BOOST_CHECK_EQUAL(allnums.out(), 0);
  
  // Prob laying inside of interval should be 0.5
  BOOST_CHECK_EQUAL(firstHalf.tot(), N);
  BOOST_CHECK_CLOSE(firstHalf.in(), 0.5 * N, 15);
  BOOST_CHECK_CLOSE(firstHalf.out(), 0.5 * N, 15);
  
    // Prob laying inside of interval should be 0.5
  BOOST_CHECK_EQUAL(secondHalf.tot(), N);
  BOOST_CHECK_CLOSE(secondHalf.in(), 0.5 * N, 15);
  BOOST_CHECK_CLOSE(secondHalf.out(), 0.5 * N, 15);
  
  // Increasing probability
  BOOST_CHECK(neg1000.in() <= pos1000.in());  
}

BOOST_AUTO_TEST_CASE(rndGauss)
{
  int N = 1000;
  
  // Check for whole interval 
  DoubleBinCheck allnums(-DBL_MAX, DBL_MAX);
  // Check for median
  DoubleBinCheck firstHalf(-DBL_MAX, 0);
  DoubleBinCheck secondHalf(0, DBL_MAX);
  // Check for increasing probability
  DoubleBinCheck oneSigma(-1, 1);
  DoubleBinCheck threeSigma(-3, 3);
  DoubleBinCheck eightSigma(-8, 8);
      
  for(int i = 0; i < N; ++i) {
    allnums.add(rng.rndGauss());
    firstHalf.add(rng.rndGauss());
    secondHalf.add(rng.rndGauss());
    oneSigma.add(rng.rndGauss());
    threeSigma.add(rng.rndGauss());
    eightSigma.add(rng.rndGauss());
  }
  
  // Prob laying inside of interval should be 1
  BOOST_CHECK_EQUAL(allnums.tot(), N);
  BOOST_CHECK_EQUAL(allnums.in(), N);
  BOOST_CHECK_EQUAL(allnums.out(), 0);
  
  // Prob laying inside of interval should be 0.5
  BOOST_CHECK_EQUAL(firstHalf.tot(), N);
  BOOST_CHECK_CLOSE(firstHalf.in(), 0.5 * N, 15);
  BOOST_CHECK_CLOSE(firstHalf.out(), 0.5 * N, 15);
  
    // Prob laying inside of interval should be 0.5
  BOOST_CHECK_EQUAL(secondHalf.tot(), N);
  BOOST_CHECK_CLOSE(secondHalf.in(), 0.5 * N, 15);
  BOOST_CHECK_CLOSE(secondHalf.out(), 0.5 * N, 15);
  
    // Prob laying inside of interval should be 0.68
  BOOST_CHECK_EQUAL(oneSigma.tot(), N);
  BOOST_CHECK_CLOSE(oneSigma.in(), 0.68 * N, 15);
  BOOST_CHECK_CLOSE(oneSigma.out(), 0.32 * N, 15);

  // Prob laying inside of interval should be 99,730 0204%
  BOOST_CHECK_EQUAL(threeSigma.tot(), N);
  BOOST_CHECK_CLOSE(threeSigma.in(), 0.99730 * N, 15);

  // Prob laying inside of interval should be  99,999 999 999 %
  BOOST_CHECK_EQUAL(eightSigma.tot(), N);
  BOOST_CHECK_CLOSE(eightSigma.in(), 0.9999 * N, 15);  
}

BOOST_AUTO_TEST_CASE(rndGaussMeanSigma)
{
  int N = 1000;
  
  for(double mean = 1; mean < 5; mean += 1) {
    for(double sigma = 0.25; sigma < 0.65; sigma += 0.1) {
      // Check for whole interval 
      DoubleBinCheck allnums(-DBL_MAX, DBL_MAX);
      // Check for median
      DoubleBinCheck firstHalf(-DBL_MAX, mean);
      DoubleBinCheck secondHalf(mean, DBL_MAX);
      // Check for increasing probability
      DoubleBinCheck oneSigma(mean - sigma, mean + sigma);
      DoubleBinCheck threeSigma(mean -  3 * sigma, mean + 3 * sigma);
      DoubleBinCheck eightSigma(mean - 8 * sigma, mean + 8 * sigma);
	  
      for(int i = 0; i < N; ++i) {
	allnums.add(rng.rndGauss(sigma, mean));
	firstHalf.add(rng.rndGauss(sigma, mean));
	secondHalf.add(rng.rndGauss(sigma, mean));
	oneSigma.add(rng.rndGauss(sigma, mean));
	threeSigma.add(rng.rndGauss(sigma, mean));
	eightSigma.add(rng.rndGauss(sigma, mean));
      }
      
      // Prob laying inside of interval should be 1
      BOOST_CHECK_EQUAL(allnums.tot(), N);
      BOOST_CHECK_EQUAL(allnums.in(), N);
      BOOST_CHECK_EQUAL(allnums.out(), 0);
      
      // Prob laying inside of interval should be 0.5
      BOOST_CHECK_EQUAL(firstHalf.tot(), N);
      BOOST_CHECK_CLOSE(firstHalf.in(), 0.5 * N, 15);
      BOOST_CHECK_CLOSE(firstHalf.out(), 0.5 * N, 15);
      
	// Prob laying inside of interval should be 0.5
      BOOST_CHECK_EQUAL(secondHalf.tot(), N);
      BOOST_CHECK_CLOSE(secondHalf.in(), 0.5 * N, 15);
      BOOST_CHECK_CLOSE(secondHalf.out(), 0.5 * N, 15);
      
      	// Prob laying inside of interval should be 0.68
      BOOST_CHECK_EQUAL(oneSigma.tot(), N);
      BOOST_CHECK_CLOSE(oneSigma.in(), 0.68 * N, 15);
      BOOST_CHECK_CLOSE(oneSigma.out(), 0.32 * N, 15);
  
      // Prob laying inside of interval should be 99,730 0204%
      BOOST_CHECK_EQUAL(threeSigma.tot(), N);
      BOOST_CHECK_CLOSE(threeSigma.in(), 0.99730 * N, 15);
  
      // Prob laying inside of interval should be  99,999 999 999 %
      BOOST_CHECK_EQUAL(eightSigma.tot(), N);
      BOOST_CHECK_CLOSE(eightSigma.in(), 0.9999 * N, 15);  
    }
  }  
}


BOOST_AUTO_TEST_CASE(rndGaussTwoNumbers)
{
  int N = 1000;
  
  // Check for whole interval 
  DoubleBinCheck allnums(-DBL_MAX, DBL_MAX);
  // Check for median
  DoubleBinCheck firstHalf(-DBL_MAX, 0);
  DoubleBinCheck secondHalf(0, DBL_MAX);
  // Check for increasing probability
  DoubleBinCheck oneSigma(-1, 1);
  DoubleBinCheck threeSigma(-3, 3);
  DoubleBinCheck eightSigma(-8, 8);
  
  double r1, r2;
  for(int i = 0; i < N; ++i) {
    rng.rndGaussTwoNumbers(r1, r2);
    
    allnums.add(r1);
    firstHalf.add(r1);
    secondHalf.add(r1);
    oneSigma.add(r1);
    threeSigma.add(r1);
    eightSigma.add(r1);
    
    allnums.add(r2);
    firstHalf.add(r2);
    secondHalf.add(r2);
    oneSigma.add(r2);
    threeSigma.add(r2);
    eightSigma.add(r2);
  }
  
  // Two random numbers were added in each round
  N *= 2;
  // Prob laying inside of interval should be 1
  BOOST_CHECK_EQUAL(allnums.tot(), N);
  BOOST_CHECK_EQUAL(allnums.in(), N);
  BOOST_CHECK_EQUAL(allnums.out(), 0);
  
  // Prob laying inside of interval should be 0.5
  BOOST_CHECK_EQUAL(firstHalf.tot(), N);
  BOOST_CHECK_CLOSE(firstHalf.in(), 0.5 * N, 15);
  BOOST_CHECK_CLOSE(firstHalf.out(), 0.5 * N, 15);
  
    // Prob laying inside of interval should be 0.5
  BOOST_CHECK_EQUAL(secondHalf.tot(), N);
  BOOST_CHECK_CLOSE(secondHalf.in(), 0.5 * N, 15);
  BOOST_CHECK_CLOSE(secondHalf.out(), 0.5 * N, 15);
  
    // Prob laying inside of interval should be 0.68
  BOOST_CHECK_EQUAL(oneSigma.tot(), N);
  BOOST_CHECK_CLOSE(oneSigma.in(), 0.68 * N, 15);
  BOOST_CHECK_CLOSE(oneSigma.out(), 0.32 * N, 15);

  // Prob laying inside of interval should be 99,730 0204%
  BOOST_CHECK_EQUAL(threeSigma.tot(), N);
  BOOST_CHECK_CLOSE(threeSigma.in(), 0.99730 * N, 15);

  // Prob laying inside of interval should be  99,999 999 999 %
  BOOST_CHECK_EQUAL(eightSigma.tot(), N);
  BOOST_CHECK_CLOSE(eightSigma.in(), 0.9999 * N, 15);  
}

BOOST_AUTO_TEST_CASE(rndGaussTwoNumbersMeanSigma)
{
  int N = 1000;
  // Two random numbers will be added in each round
  int NTwo = 2 * N;
  double r1, r2;
  
  for(double mean = 1; mean < 5; mean += 1) {
    for(double sigma = 0.25; sigma < 0.65; sigma += 0.1) {
      // Check for whole interval 
      DoubleBinCheck allnums(-DBL_MAX, DBL_MAX);
      // Check for median
      DoubleBinCheck firstHalf(-DBL_MAX, mean);
      DoubleBinCheck secondHalf(mean, DBL_MAX);
      // Check for increasing probability
      DoubleBinCheck oneSigma(mean - sigma, mean + sigma);
      DoubleBinCheck threeSigma(mean -  3 * sigma, mean + 3 * sigma);
      DoubleBinCheck eightSigma(mean - 8 * sigma, mean + 8 * sigma);
	  
      for(int i = 0; i < N; ++i) {
	rng.rndGaussTwoNumbers(r1, r2, sigma, mean);
	
	allnums.add(r1);
	firstHalf.add(r1);
	secondHalf.add(r1);
	oneSigma.add(r1);
	threeSigma.add(r1);
	eightSigma.add(r1);
	
	allnums.add(r2);
	firstHalf.add(r2);
	secondHalf.add(r2);
	oneSigma.add(r2);
	threeSigma.add(r2);
	eightSigma.add(r2);
      }
      
      // Prob laying inside of interval should be 1
      BOOST_CHECK_EQUAL(allnums.tot(), NTwo);
      BOOST_CHECK_EQUAL(allnums.in(), NTwo);
      BOOST_CHECK_EQUAL(allnums.out(), 0);
      
      // Prob laying inside of interval should be 0.5
      BOOST_CHECK_EQUAL(firstHalf.tot(), NTwo);
      BOOST_CHECK_CLOSE(firstHalf.in(), 0.5 * NTwo, 15);
      BOOST_CHECK_CLOSE(firstHalf.out(), 0.5 * NTwo, 15);
      
	// Prob laying inside of interval should be 0.5
      BOOST_CHECK_EQUAL(secondHalf.tot(), NTwo);
      BOOST_CHECK_CLOSE(secondHalf.in(), 0.5 * NTwo, 15);
      BOOST_CHECK_CLOSE(secondHalf.out(), 0.5 * NTwo, 15);
      
      	// Prob laying inside of interval should be 0.68
      BOOST_CHECK_EQUAL(oneSigma.tot(), NTwo);
      BOOST_CHECK_CLOSE(oneSigma.in(), 0.68 * NTwo, 15);
      BOOST_CHECK_CLOSE(oneSigma.out(), 0.32 * NTwo, 15);
  
      // Prob laying inside of interval should be 99,730 0204%
      BOOST_CHECK_EQUAL(threeSigma.tot(), NTwo);
      BOOST_CHECK_CLOSE(threeSigma.in(), 0.99730 * NTwo, 15);
  
      // Prob laying inside of interval should be  99,999 999 999 %
      BOOST_CHECK_EQUAL(eightSigma.tot(), NTwo);
      BOOST_CHECK_CLOSE(eightSigma.in(), 0.9999 * NTwo, 15);  
    }
  }  
}

BOOST_AUTO_TEST_CASE(rndPoisson)
{
  int N = 1000;
  
  for(double mean = 5; mean < 21; mean += 5) {
    // Check for whole interval 
    LongBinCheck zero2max(0, LONG_MAX);
    // Check for median
    LongBinCheck zero(0, mean);
    LongBinCheck secondHalf(mean, LONG_MAX);
    // Check for decreasing probability away from mean
    LongBinCheck nearMean(mean - std::sqrt(mean), mean + std::sqrt(mean));
    LongBinCheck awayMean(mean + std::sqrt(mean), mean + 3 * std::sqrt(mean));
	
    for(int i = 0; i < N; ++i) {
      zero2max.add(rng.rndPoisson(mean));
      zero.add(rng.rndPoisson(mean));
      secondHalf.add(rng.rndPoisson(mean));
      nearMean.add(rng.rndPoisson(mean));
      awayMean.add(rng.rndPoisson(mean));
    }
    
    // Prob laying inside of interval should be 1
    BOOST_CHECK_EQUAL(zero2max.tot(), N);
    BOOST_CHECK_EQUAL(zero2max.in(), N);
    BOOST_CHECK_EQUAL(zero2max.out(), 0);
    
    // Prob laying inside of interval should be 0.5
    BOOST_CHECK_EQUAL(zero.tot(), N);
    BOOST_CHECK_CLOSE(zero.in(), 0.5 * N, 30);
    BOOST_CHECK_CLOSE(zero.out(), 0.5 * N, 30);
    
      // Prob laying inside of interval should be 0.5
    BOOST_CHECK_EQUAL(secondHalf.tot(), N);
    BOOST_CHECK_CLOSE(secondHalf.in(), 0.5 * N, 30);
    BOOST_CHECK_CLOSE(secondHalf.out(), 0.5 * N, 30);
    
    // Prob should decrease away from mean
    BOOST_CHECK(nearMean.in() > awayMean.in()); 
  }  
}

BOOST_AUTO_TEST_CASE(rndBW)
{
  int N = 1000;
  
  for(double mean = 5; mean < 21; mean += 5) {
    for (double gamma = 1; gamma < 5; gamma += 1) {
      // Check for whole interval 
      DoubleBinCheck allnums(-DBL_MAX, DBL_MAX);
      // Check for median
      DoubleBinCheck firstHalf(-DBL_MAX, mean);
      DoubleBinCheck secondHalf(mean, DBL_MAX);
      // Check for decreasing probability away from mean
      DoubleBinCheck nearMean(mean - std::sqrt(gamma), mean + std::sqrt(gamma));
      DoubleBinCheck awayMeanOne(mean + std::sqrt(gamma), mean + 3 * std::sqrt(gamma));
      // Additional check for symmetry
      DoubleBinCheck awayMeanTwo(mean - 3 * std::sqrt(gamma), mean - std::sqrt(gamma));
	  
      for(int i = 0; i < N; ++i) {
	allnums.add(rng.rndBW(mean, gamma));
	firstHalf.add(rng.rndBW(mean, gamma));
	secondHalf.add(rng.rndBW(mean, gamma));
	nearMean.add(rng.rndBW(mean, gamma));
	awayMeanOne.add(rng.rndBW(mean, gamma));
	awayMeanTwo.add(rng.rndBW(mean, gamma));
      }
      
      // Prob laying inside of interval should be 1
      BOOST_CHECK_EQUAL(allnums.tot(), N);
      BOOST_CHECK_EQUAL(allnums.in(), N);
      BOOST_CHECK_EQUAL(allnums.out(), 0);
      
      // Prob laying inside of interval should be 0.5
      BOOST_CHECK_EQUAL(firstHalf.tot(), N);
      BOOST_CHECK_CLOSE(firstHalf.in(), 0.5 * N, 30);
      BOOST_CHECK_CLOSE(firstHalf.out(), 0.5 * N, 30);
      
	// Prob laying inside of interval should be 0.5
      BOOST_CHECK_EQUAL(secondHalf.tot(), N);
      BOOST_CHECK_CLOSE(secondHalf.in(), 0.5 * N, 30);
      BOOST_CHECK_CLOSE(secondHalf.out(), 0.5 * N, 30);
      
      // Prob should decrease away from mean
      BOOST_CHECK(nearMean.in() > awayMeanOne.in());
      
      // Prob should be roughly symmetric with a shift to the left or right
      BOOST_CHECK_CLOSE(1.0 * awayMeanTwo.in(), 1.0 * awayMeanOne.in(), 60);
    }
  }  
}

BOOST_AUTO_TEST_CASE(rndBWCut)
{
  int N = 1000;
  
  double cutValue = 2;
  for(double mean = 5; mean < 10; mean += 2) {
    for (double gamma = 1; gamma < 5; gamma += 1) {
      // Check for whole interval 
      DoubleBinCheck allnums(-DBL_MAX, DBL_MAX);
      // Check for cut interval 
      DoubleBinCheck inCut(mean - cutValue, mean + cutValue);
      // Check for median
      DoubleBinCheck firstHalf(-DBL_MAX, mean);
      DoubleBinCheck secondHalf(mean, DBL_MAX);
      // Check for decreasing probability away from mean
      DoubleBinCheck nearMean(mean - std::sqrt(gamma), mean + std::sqrt(gamma));
      DoubleBinCheck awayMeanOne(mean + std::sqrt(gamma), mean + 3 * std::sqrt(gamma));
      // Additional check for symmetry
      DoubleBinCheck awayMeanTwo(mean - 3 * std::sqrt(gamma), mean - std::sqrt(gamma));
	  
      for(int i = 0; i < N; ++i) {
	allnums.add(rng.rndBW(mean, gamma, cutValue));
	inCut.add(rng.rndBW(mean, gamma, cutValue));
	firstHalf.add(rng.rndBW(mean, gamma, cutValue));
	secondHalf.add(rng.rndBW(mean, gamma, cutValue));
	nearMean.add(rng.rndBW(mean, gamma, cutValue));
	awayMeanOne.add(rng.rndBW(mean, gamma, cutValue));
	awayMeanTwo.add(rng.rndBW(mean, gamma, cutValue));
      }
      
      // Prob laying inside of interval should be 1
      BOOST_CHECK_EQUAL(allnums.tot(), N);
      BOOST_CHECK_EQUAL(allnums.in(), N);
      BOOST_CHECK_EQUAL(allnums.out(), 0);
      
      // Prob laying inside of interval should be 1
      BOOST_CHECK_EQUAL(inCut.tot(), N);
      BOOST_CHECK_EQUAL(inCut.in(), N);
      BOOST_CHECK_EQUAL(inCut.out(), 0);
      
      // Prob laying inside of interval should be 0.5
      BOOST_CHECK_EQUAL(firstHalf.tot(), N);
      BOOST_CHECK_CLOSE(firstHalf.in(), 0.5 * N, 30);
      BOOST_CHECK_CLOSE(firstHalf.out(), 0.5 * N, 30);
      
	// Prob laying inside of interval should be 0.5
      BOOST_CHECK_EQUAL(secondHalf.tot(), N);
      BOOST_CHECK_CLOSE(secondHalf.in(), 0.5 * N, 30);
      BOOST_CHECK_CLOSE(secondHalf.out(), 0.5 * N, 30);
      
      // Prob should decrease away from mean
      BOOST_CHECK(nearMean.in() > awayMeanOne.in());
      
      // Prob should be roughly symmetric with a shift to the left or right
      BOOST_CHECK_CLOSE(1.0 * awayMeanTwo.in(), 1.0 * awayMeanOne.in(), 60);
    }
  }  
}

BOOST_AUTO_TEST_CASE(rndRelBW)
{
  int N = 1000;
  
  for(double mean = 5; mean < 21; mean += 5) {
    for (double gamma = 1; gamma < 5; gamma += 1) {
      // Check for whole interval 
      DoubleBinCheck allnums(-DBL_MAX, DBL_MAX);
      // Check for median
      DoubleBinCheck firstHalf(-DBL_MAX, mean);
      DoubleBinCheck secondHalf(mean, DBL_MAX);
      // Check for decreasing probability away from mean
      DoubleBinCheck nearMean(mean - std::sqrt(gamma), mean + std::sqrt(gamma));
      DoubleBinCheck awayMeanOne(mean + std::sqrt(gamma), mean + 3 * std::sqrt(gamma));
      // Additional check for symmetry
      DoubleBinCheck awayMeanTwo(mean - 3 * std::sqrt(gamma), mean - std::sqrt(gamma));
	  
      for(int i = 0; i < N; ++i) {
	allnums.add(rng.rndRelBW<double>(mean, gamma));
	firstHalf.add(rng.rndRelBW<double>(mean, gamma));
	secondHalf.add(rng.rndRelBW<double>(mean, gamma));
	nearMean.add(rng.rndRelBW<double>(mean, gamma));
	awayMeanOne.add(rng.rndRelBW<double>(mean, gamma));
	awayMeanTwo.add(rng.rndRelBW<double>(mean, gamma));
      }
      
      // Prob laying inside of interval should be 1
      BOOST_CHECK_EQUAL(allnums.tot(), N);
      BOOST_CHECK_EQUAL(allnums.in(), N);
      BOOST_CHECK_EQUAL(allnums.out(), 0);
      
      // Prob laying inside of interval should be 0.5
      BOOST_CHECK_EQUAL(firstHalf.tot(), N);
      BOOST_CHECK_CLOSE(firstHalf.in(), 0.5 * N, 50);
      BOOST_CHECK_CLOSE(firstHalf.out(), 0.5 * N, 50);
      
	// Prob laying inside of interval should be 0.5
      BOOST_CHECK_EQUAL(secondHalf.tot(), N);
      BOOST_CHECK_CLOSE(secondHalf.in(), 0.5 * N, 50);
      BOOST_CHECK_CLOSE(secondHalf.out(), 0.5 * N, 50);
      
      // Prob should decrease away from mean
      BOOST_CHECK(nearMean.in() > awayMeanOne.in());
      
      // Prob should be roughly symmetric with a shift to the left or right
      BOOST_CHECK_CLOSE(1.0 * awayMeanTwo.in(), 1.0 * awayMeanOne.in(), 200);
    }
  }  
}

BOOST_AUTO_TEST_CASE(rndRelBWCut)
{
  int N = 1000;
  
  double cutValue = 2;
  for(double mean = 5; mean < 10; mean += 2) {
    for (double gamma = 1; gamma < 5; gamma += 1) {
      // Check for whole interval 
      DoubleBinCheck allnums(-DBL_MAX, DBL_MAX);
      // Check for cut interval 
      DoubleBinCheck inCut(mean - cutValue, mean + cutValue);
      // Check for median
      DoubleBinCheck firstHalf(-DBL_MAX, mean);
      DoubleBinCheck secondHalf(mean, DBL_MAX);
      // Check for decreasing probability away from mean
      DoubleBinCheck nearMean(mean - std::sqrt(gamma), mean + std::sqrt(gamma));
      DoubleBinCheck awayMeanOne(mean + std::sqrt(gamma), mean + 3 * std::sqrt(gamma));
      // Additional check for symmetry
      DoubleBinCheck awayMeanTwo(mean - 3 * std::sqrt(gamma), mean - std::sqrt(gamma));
	  
      for(int i = 0; i < N; ++i) {
	allnums.add(rng.rndRelBW<double>(mean, gamma, cutValue));
	inCut.add(rng.rndRelBW<double>(mean, gamma, cutValue));
	firstHalf.add(rng.rndRelBW<double>(mean, gamma, cutValue));
	secondHalf.add(rng.rndRelBW<double>(mean, gamma, cutValue));
	nearMean.add(rng.rndRelBW<double>(mean, gamma, cutValue));
	awayMeanOne.add(rng.rndRelBW<double>(mean, gamma, cutValue));
	awayMeanTwo.add(rng.rndRelBW<double>(mean, gamma, cutValue));
      }
      
      // Prob laying inside of interval should be 1
      BOOST_CHECK_EQUAL(allnums.tot(), N);
      BOOST_CHECK_EQUAL(allnums.in(), N);
      BOOST_CHECK_EQUAL(allnums.out(), 0);
      
      // Prob laying inside of interval should be 1
      BOOST_CHECK_EQUAL(inCut.tot(), N);
      BOOST_CHECK_EQUAL(inCut.in(), N);
      BOOST_CHECK_EQUAL(inCut.out(), 0);
      
      // Prob laying inside of interval should be 0.5
      BOOST_CHECK_EQUAL(firstHalf.tot(), N);
      BOOST_CHECK_CLOSE(firstHalf.in(), 0.5 * N, 30);
      BOOST_CHECK_CLOSE(firstHalf.out(), 0.5 * N, 30);
      
	// Prob laying inside of interval should be 0.5
      BOOST_CHECK_EQUAL(secondHalf.tot(), N);
      BOOST_CHECK_CLOSE(secondHalf.in(), 0.5 * N, 30);
      BOOST_CHECK_CLOSE(secondHalf.out(), 0.5 * N, 30);
      
      // Prob should decrease away from mean
      BOOST_CHECK(nearMean.in() > awayMeanOne.in());
      
      // Prob should be roughly symmetric with a shift to the left or right
      BOOST_CHECK_CLOSE(1.0 * awayMeanTwo.in(), 1.0 * awayMeanOne.in(), 120);
    }
  }  
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* ThePEG_Repository_Test_RandomGenerator_H */
