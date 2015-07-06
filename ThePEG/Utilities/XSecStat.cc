// -*- C++ -*-
//
// XSecStat.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the XSecStat class.
//

#include "XSecStat.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

namespace ThePEG {

PersistentOStream & operator<<(PersistentOStream & os, const XSecStat & x) {
  x.output(os);
  return os;
}

PersistentIStream & operator>>(PersistentIStream & is, XSecStat & x) {
  x.input(is);
  return is;
}

}

using namespace ThePEG;

void XSecStat::output(PersistentOStream & os) const {
  os << ounit(theMaxXSec,picobarn) << theAttempts << theAccepted
     << theSumWeights << theSumWeights2 << theLastWeight;
}

void XSecStat::input(PersistentIStream & is) {
  is >> iunit(theMaxXSec,picobarn) >> theAttempts >> theAccepted
     >> theSumWeights >> theSumWeights2 >> theLastWeight;
}

