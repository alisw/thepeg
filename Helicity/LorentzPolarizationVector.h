// -*- C++ -*-
//
// LorentzPolarizationVector.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2011 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_LorentzPolarizationVector_H
#define ThePEG_LorentzPolarizationVector_H
// This is the declaration of the LorentzPolarizationVector class.

#include "ThePEG/Config/Unitsystem.h"
#include "ThePEG/Config/Complex.h"
#include "ThePEG/Vectors/LorentzVector.h"

namespace ThePEG {
  namespace Helicity {
    /// Convenience typedef.
    typedef LorentzVector<complex<double> > LorentzPolarizationVector;
    /// Convenience typedef.
    typedef LorentzVector<complex<Energy> > LorentzPolarizationVectorE;
  }
}
#endif
