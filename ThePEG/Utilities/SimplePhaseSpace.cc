// -*- C++ -*-
//
// SimplePhaseSpace.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//

#include "SimplePhaseSpace.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
#include "SimplePhaseSpace.tcc"
#endif

using namespace ThePEG;

Energy SimplePhaseSpace::getMagnitude(Energy2 s, Energy m1, Energy m2)
{
  const Energy2 eps = 10.0*s*Constants::epsilon;
  if ( m1 < ZERO && sqr(m1) < eps ) m1 = ZERO;
  if ( m2 < ZERO && sqr(m2) < eps ) m2 = ZERO;
  if ( m1 >= ZERO && m2 >= ZERO ) {
    Energy2 aa = s - sqr(m1+m2);
    if ( aa < ZERO && aa > -eps ) return ZERO;
    if ( aa < ZERO ) throw ImpossibleKinematics();
    return 0.5*sqrt(aa*(s-sqr(m1-m2))/s);
  }
  Energy2 m12 = m1 < ZERO? -sqr(m1): sqr(m1);
  Energy2 m22 = m2 < ZERO? -sqr(m2): sqr(m2);
  Energy2 r2 = 0.25*(sqr(m12) + sqr(m22 - s) -2.0*m12*(m22 + s))/s;
  if ( r2 < ZERO || r2 + m12 < ZERO || r2 + m22 < ZERO )
    throw ImpossibleKinematics();
  return sqrt(r2);
}

vector<LorentzMomentum> SimplePhaseSpace::
CMSn(Energy m0, const vector<Energy> & m)
{
  using Constants::pi;

  // Setup constants.
  int Np = m.size();
  vector<LorentzMomentum> ret(Np);
  Energy summ = std::accumulate(m.begin(), m.end(), Energy());
  if ( summ >= m0 ) throw ImpossibleKinematics();

  while ( true ) {
    // First get an ordered list of random numbers.
    vector<double> rndv(Np);
    rndv[0] = 1.0;
    rndv.back() = 0.0;
    for ( int i = 1; i < Np - 1; ++i ) rndv[i] = UseRandom::rnd();
    std::sort(rndv.begin() + 1, rndv.end() - 1, std::greater<double>());

    // Now setup masses of subsystems.
    vector<Energy> sm(Np);
    Energy tmass = m0 - summ;
    Energy tmp = summ;
    for ( int i = 0; i < Np; ++i ) {
      sm[i] = rndv[i]*tmass + tmp;
      tmp -= m[i];
    }

    // Now the magnitude of all the momenta can be calculated. This
    // gives the weight.
    double weight = 1.0;
    vector<Energy> p(Np);
    p[Np - 1] = getMagnitude(sqr(sm[Np - 2]), m[Np -2], sm[Np - 1]);
    for ( int i = Np - 2; i >= 0; --i )
      weight *= (p[i] = getMagnitude(sqr(sm[i]), m[i], sm[i + 1]))/sm[i];
    if ( weight > UseRandom::rnd() ) continue;

    // Now we just have to generate the angles.
    ret[Np - 1] = LorentzMomentum(ZERO, ZERO, ZERO, m[Np - 1]);
    for ( int i = Np - 2; i >= 0; --i ) {
      Momentum3 p3 = polar3Vector(p[i], 2.0*UseRandom::rnd() - 1.0,
				  2.0*pi*UseRandom::rnd());
      ret[i] = LorentzMomentum(-p3, sqrt(sqr(p[i]) + sqr(m[i])));
      if ( i == Np -2 ) {
	ret[Np - 1] = LorentzMomentum(p3, sqrt(sqr(m[Np - 1]) + p3.mag2()));
      } else {
	Boost bv = p3*(1.0/sqrt(sqr(p[i]) + sqr(sm[i + 1])));
	if ( bv.mag2() >= 1.0 ) throw ImpossibleKinematics();
	LorentzRotation r(bv);
	for ( int j = i + 1; j < Np; ++j ) ret[j]*=r.one();
      }
    }
    return ret;
  }
}

