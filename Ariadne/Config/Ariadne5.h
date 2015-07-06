// -*- C++ -*-
#ifndef Ariadne5_H
#define Ariadne5_H

// This is the main config header file for Ariande.

#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/Vectors/Lorentz5Vector.h"

namespace Ariadne5 {
using namespace ThePEG;

/** Create a LorentzVector given mass, pt, rapidity, and azimuth
 * angle. */
inline LorentzMomentum
ptRapidity(Energy mass, Energy pt, double y = 0.0, double phi = 0.0) {
  Energy mt = sqrt( sqr(mass) + sqr(pt) );
  return lightCone(mt*exp(y), mt*exp(-y), pt*cos(phi), pt*sin(phi));
}

/**
 * This is a fail-safe version of LorentzMomentum::rapidity().
 */
inline double rapidity(const LorentzMomentum & p,
		double ymax = Constants::MaxRapidity) {
  if ( p.e() > abs(p.z()) ) return p.rapidity();
  return p.z() > ZERO? ymax: ( p.z() < ZERO? -ymax: 0.0 );
}

/** Create a Lorentz5Vector giving its light-cone and transverse
 *  components. */
template <typename Value>
inline Lorentz5Vector<Value>
lightCone5(Value plus, Value minus, Value x, Value y) {
  Lorentz5Vector<Value> r(x, y, 0.5*(plus-minus), 0.5*(plus+minus));
  return r;
}

/** Create a Lorentz5Vector giving its light-cone and transverse
 *  components and a possibly inconsistent mass. */
template <typename Value>
inline Lorentz5Vector<Value>
lightCone5(Value plus, Value minus, Value x, Value y, Value m) {
  Lorentz5Vector<Value> r(x, y, 0.5*(plus-minus), 0.5*(plus+minus), m);
  return r;
}

/** Create a Lorentz5Vector giving its light-cone components. */
template <typename Value>
inline Lorentz5Vector<Value>
lightCone5(Value plus, Value minus) {
  Lorentz5Vector<Value> r(ZERO, ZERO, 0.5*(plus-minus), 0.5*(plus+minus));
  return r;
}

/** Create a Lorentz5Vector giving its light-cone and transverse
 *  components. */
template <typename Value>
inline Lorentz5Vector<Value>
lightCone5(Value plus, Value minus, Transverse<Value> pt) {
  Lorentz5Vector<Value> r(pt.x(), pt.y(), 0.5*(plus-minus), 0.5*(plus+minus));
  return r;
}

/** Create a Lorentz5Vector giving its light-cone and transverse
 *  components and a possibly inconsistent mass. */
template <typename Value>
inline Lorentz5Vector<Value>
lightCone5(Value plus, Value minus, Transverse<Value> pt, Value m) {
  Lorentz5Vector<Value> r(pt.x(), pt.y(), 0.5*(plus-minus), 0.5*(plus+minus), m);
  return r;
}

/** Create a Lorentz5Vector giving its light-cone and transverse
 *  components. If the current Direction<0> is reversed, so is the
 *  z-component. */
template <typename Value>
inline Lorentz5Vector<Value>
lightCone5Dir(Value plus, Value minus,
	     Value x = Value(), Value y = Value()) {
  Lorentz5Vector<Value> r(x, y, Direction<0>::dir()*0.5*(plus - minus),
			  0.5*(plus + minus));
  return r;
}

/** Create a Lorentz5Vector giving its light-cone and transverse
 *  components. If the current Direction<0> is reversed, so is the
 *  z-component. */
template <typename Value>
inline Lorentz5Vector<Value>
lightCone5Dir(Value plus, Value minus, Transverse<Value> pt) {
  Lorentz5Vector<Value> r(pt.x(), pt.y(), Direction<0>::dir()*0.5*(plus - minus),
			  0.5*(plus + minus));
  return r;

}

}

#endif /* Ariadne5_H */

