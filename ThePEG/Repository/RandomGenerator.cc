// -*- C++ -*-
//
// RandomGenerator.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the RandomGenerator class.
//

#include "RandomGenerator.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "gsl/gsl_randist.h"

extern "C" {
  typedef struct
  {
    ThePEG::RandomGenerator * r;
  }
  thepeg_random_state_t;

  void thepeg_random_set(void *, unsigned long int) {}
  double thepeg_random_get_double(void * s) {
  return static_cast<thepeg_random_state_t *>(s)->r->rnd();
  }
  unsigned long int thepeg_random_get(void * s) {
    return static_cast<unsigned long int>(std::numeric_limits<unsigned long int>::max()*thepeg_random_get_double(s));
  }

  static const gsl_rng_type thepeg_random_type =
  {"thepeg_random",
   (unsigned long int)std::numeric_limits<unsigned>::max(),
   0,
   sizeof(thepeg_random_state_t),
   &thepeg_random_set,
   &thepeg_random_get,
   &thepeg_random_get_double};
   
  const gsl_rng_type *gsl_rng_thepeg_random = &thepeg_random_type;

}
using namespace ThePEG;

RandomGenerator::RandomGenerator()
  : theNumbers(1000), theSize(1000), theSeed(-1),
    savedGauss(0.0), gaussSaved(false) {
  nextNumber = theNumbers.end();
  gsl = gsl_rng_alloc(gsl_rng_thepeg_random);
  static_cast<thepeg_random_state_t *>(gsl->state)->r = this;
}

RandomGenerator::RandomGenerator(const RandomGenerator & rg)
  : Interfaced(rg), theNumbers(rg.theNumbers), theSize(rg.theSize),
    theSeed(rg.theSeed), savedGauss(rg.savedGauss),
    gaussSaved(rg.gaussSaved)  {
  nextNumber = theNumbers.begin() +
    ( RndVector::const_iterator(rg.nextNumber) - rg.theNumbers.begin() );
  gsl = gsl_rng_alloc(gsl_rng_thepeg_random);
  static_cast<thepeg_random_state_t *>(gsl->state)->r = this;
}

RandomGenerator::~RandomGenerator() {
  gsl_rng_free(gsl);
}

void RandomGenerator::doinit() {
  if ( theSeed != 0 ) setSeed(theSeed);
  flush();
}

void RandomGenerator::setSize(size_type newSize) {
  RndVector newNumbers(newSize);
  RndVector::iterator nextNew = newNumbers.end() -
    min( int(theNumbers.end() - nextNumber), int(newSize) );
  for ( RndVector::iterator i = nextNew; i != newNumbers.end(); ++i )
    *i = *nextNumber++;
  RndVector::difference_type pos = nextNew - newNumbers.begin();
  theNumbers.swap(newNumbers);
  nextNumber = theNumbers.begin() + pos;
}

bool RandomGenerator::prndbool(double p) {
  if ( p >= 1.0 ) return true;
  if ( p <= 0.0 ) return false;
  double r = rnd();
  if ( r < p ) {
    push_back(r/p);
    return true;
  } else {
    push_back((r - p)/(1.0 - p));
    return false;
  }
}

int RandomGenerator::rndsign(double p1, double p2, double p3) {
  double sum = p1 + p2 + p3;
  double r = rnd()*sum;
  if ( r < p1 )  return -1;
  else if ( r < p1 + p2 ) return 0;
  else return 1;
}

int RandomGenerator::prndsign(double p1, double p2, double p3) {
  double sum = p1 + p2 + p3;
  double r = rnd()*sum;
  if ( r < p1 ) {
    push_back(r/p1);
    return -1;
  } else if ( r < p1 + p2 ) {
    push_back((r - p1)/p2);
    return 0;
  } else {
    push_back((r - p1 - p2)/p3);
    return 1;
  }
}

int RandomGenerator::rnd4(double p0, double p1, double p2, double p3) {
  double sum = p0 + p1 + p2 + p3;
  double r = rnd()*sum;
  if ( r < p0 )  return 0;
  else if ( r < p0 + p1 ) return 1;
  else if ( r < p0 + p1 + p2 ) return 2;
  else return 3;
}

long RandomGenerator::rndPoisson(double mean) {
  return gsl_ran_poisson(gsl, mean);
}

void RandomGenerator::persistentOutput(PersistentOStream & os) const {
  os << theNumbers
     << RndVector::const_iterator(nextNumber) - theNumbers.begin() << theSize
     << theSeed << savedGauss << gaussSaved;
}

void RandomGenerator::persistentInput(PersistentIStream & is, int) {
  RndVector::difference_type pos;
  is >> theNumbers >> pos >> theSize >> theSeed >> savedGauss >> gaussSaved;
  nextNumber = theNumbers.begin() + pos;
}

ClassDescription<RandomGenerator> RandomGenerator::initRandomGenerator;

void RandomGenerator::Init() {

  static ClassDocumentation<RandomGenerator> documentation
    ("There is no documentation for the ThePEG::RandomGenerator class");

  static Parameter<RandomGenerator,size_type> interfaceSize
    ("CacheSize",
     "The Random numbers are generated in chunks of this size.",
     &RandomGenerator::theSize, 1000, 10, 100000, true, false, true,
     &RandomGenerator::setSize);

  static Parameter<RandomGenerator,long> interfaceSeed
    ("Seed",
     "The seed with which this random generator is initialized. "
     "If set to -1, the default build-in seed will be used. If set to zero, no seed will "
     "be set.",
     &RandomGenerator::theSeed, -1, -1, 100000000, true, false, false);
  interfaceSeed.setHasDefault(false);

  interfaceSize.rank(10);
  interfaceSeed.rank(9);

}

