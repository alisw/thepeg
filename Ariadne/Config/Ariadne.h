// -*- C++ -*-
#ifndef ARIADNE_H
#define ARIADNE_H

// This is the main config header file for Ariande.

#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/Vectors/Lorentz5Vector.h"

namespace Ariadne {
using namespace ThePEG;

/** Create a LorentzVector given mass, pt, rapidity, and azimuth
 * angle. */
inline LorentzMomentum
ptRapidity(Energy mass, Energy pt, double y = 0.0, double phi = 0.0) {
  Energy mt = sqrt( sqr(mass) + sqr(pt) );
  return lightCone(mt*exp(y), mt*exp(-y), pt*cos(phi), pt*sin(phi));
}

}

#endif /* ARIADNE_H */

