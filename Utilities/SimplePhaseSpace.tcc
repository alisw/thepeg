// -*- C++ -*-
//
// SimplePhaseSpace.tcc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined templated member
// functions of the SimplePhaseSpace class.
//

namespace ThePEG {

template <typename PType>
void SimplePhaseSpace::CMS(PType & p1, PType & p2, Energy2 s)
{
  typedef ParticleTraits<PType> Traits;
  Energy m1 = Traits::mass(p1); Energy m2 = Traits::mass(p2);
  Energy z = getMagnitude(s, m1, m2);
  Energy2 m12 = m1 >= ZERO ? sqr(m1) : -sqr(m1);
  Energy2 m22 = m2 >= ZERO ? sqr(m2) : -sqr(m2);
  Energy2 c1 = (s+m12-m22);
  Energy2 c2 = (s-m12+m22);
  Traits::set5Momentum(p1, Lorentz5Momentum(ZERO, ZERO, z, (c1 > ZERO ? 1. : -1.) * 0.5*sqrt(sqr(c1)/s), m1));
  Traits::set5Momentum(p2, Lorentz5Momentum(ZERO, ZERO, -z, (c2 > ZERO ? 1. : -1.) * 0.5*sqrt(sqr(c2)/s), m2));
}

template <typename PType>
void SimplePhaseSpace::CMS(Energy2 s, PType & p1, PType & p2)
{
  CMS(p1, p2, s, 2.0*UseRandom::rnd() - 1.0, Constants::twopi*UseRandom::rnd());
}

template <typename PType>
void SimplePhaseSpace::CMS(PType & p1, PType & p2, Energy2 s,
			   double cthe, double phi)
{
  typedef ParticleTraits<PType> Traits;
  Energy r = getMagnitude(s, Traits::mass(p1), Traits::mass(p2));
  double sthe = sqrt(1.0-sqr(cthe));
  Momentum3 p(r*sthe*cos(phi), r*sthe*sin(phi), r*cthe);
  Traits::set3Momentum(p1, p);
  Traits::set3Momentum(p2, -p);
}

template <typename PType>
void SimplePhaseSpace::
CMS(PType & p1, PType & p2, PType & p3,
    Energy2 s, double x1, double x3)
{
  CMS(p1, p2, p3, s, x1, x3,
      Constants::twopi*UseRandom::rnd(), acos(2.0*UseRandom::rnd() - 1.0),
      Constants::twopi*UseRandom::rnd());
}

template <typename PType>
void SimplePhaseSpace::
CMS(PType & p1, PType & p2, PType & p3, Energy2 s,
    double x1, double x3, double phii, double the, double phi)
{
  typedef ParticleTraits<PType> Traits;
  Energy Etot = sqrt(s);
  Energy m1 = Traits::mass(p1);
  Energy m2 = Traits::mass(p2);
  Energy m3 = Traits::mass(p3);
  Energy e1 = 0.5*x1*Etot;
  Energy e3 = 0.5*x3*Etot;
  Energy e2 = Etot - e1 - e3;
  if ( e1 < m1 || e2 < m2 || e3 < m3 ) throw ImpossibleKinematics();
  Energy r1 = sqrt(sqr(e1)-sqr(m1));
  Energy r2 = sqrt(sqr(e2)-sqr(m2));
  Energy r3 = sqrt(sqr(e3)-sqr(m3));
  Traits::set3Momentum(p1, Momentum3(ZERO, ZERO, r1));
  double cthe2 = (sqr(r3)-sqr(r2)-sqr(r1))/(2.0*r2*r1);
  double cthe3 = (sqr(r2)-sqr(r3)-sqr(r1))/(2.0*r3*r1);
  if ( abs(cthe2) > 1.0 || abs(cthe3) > 1.0 ) throw ImpossibleKinematics();
  double sthe2 = sqrt(1.0-sqr(cthe2));
  Energy px = r2*sthe2*cos(phii);
  Energy py = r2*sthe2*sin(phii);
  Traits::set3Momentum(p2, Momentum3(px, py, r2*cthe2));
  Traits::set3Momentum(p3, Momentum3(-px, -py, r3*cthe3));
  if ( the == 0.0 && phi == 0.0 ) return;
  LorentzRotation r;
  r.rotateZ(phi);
  r.rotateX(the);
  Traits::transform(p1, r);
  Traits::transform(p2, r);
  Traits::transform(p3, r);
}

template <typename PType>
void SimplePhaseSpace::
CMS(PType & p1, PType & p2, Energy2 s, Energy2 t, double phi,
    const PType & p0) {
  typedef ParticleTraits<PType> Traits;
  Energy r = getMagnitude(s, Traits::mass(p1), Traits::mass(p2));
  Energy e = sqrt(sqr(r) + sqr(Traits::mass(p1)));
  Energy r0 = Traits::momentum(p0).rho();
  Energy e0 = Traits::momentum(p0).e();
  double cthe = (t + sqr(e - e0) + sqr(r) + sqr(r0))/(2.0*r*r0);
  if ( abs(cthe) > 1.0 ) throw ImpossibleKinematics();
  double sthe = sqrt(1.0-sqr(cthe));
  Momentum3 p(r*sthe*cos(phi), r*sthe*sin(phi), r*cthe);
  Traits::set3Momentum(p1, p);
  Traits::set3Momentum(p2, -p);
  if ( Traits::momentum(p0).perp2() > ZERO ) {
    LorentzRotation r;
    r.rotateX(Traits::momentum(p0).theta());
    r.rotateZ(Traits::momentum(p0).phi());
    Traits::transform(p1, r);
    Traits::transform(p2, r);
  }

}

template <typename Container>
void SimplePhaseSpace::CMSn(Container & particles, Energy m0)
{
  typedef typename Container::value_type PType;
  typedef typename Container::iterator Iterator;
  if ( particles.size() == 2 ) {
    Iterator it = particles.begin();
    PType & p1 = *it++;
    PType & p2 = *it;
    CMS(sqr(m0), p1, p2);
    return;
  }
  typedef ParticleTraits<PType> Traits;
  vector<Energy> masses(particles.size());
  int j = 0;
  for ( Iterator i = particles.begin();i != particles.end(); ++i, ++j )
    masses[j] = Traits::mass(*i);
  vector<LorentzMomentum> p = CMSn(m0, masses);
  j = 0;
  for ( Iterator i = particles.begin();i != particles.end(); ++i, ++j )
    Traits::set5Momentum(*i, p[j]);
}

}
