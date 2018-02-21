// -*- C++ -*-
//
// SubProcess.tcc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined templated member
// functions of the SubProcess class.
//

namespace ThePEG {

template <class InputIterator>
void SubProcess::setIntermediates(InputIterator first, InputIterator last) {
  theIntermediates = ParticleVector(first, last);
}

template <class InputIterator>
void SubProcess::setOutgoing(InputIterator first, InputIterator last) {
  theOutgoing = ParticleVector(first, last);
}

}
