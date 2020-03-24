// -*- C++ -*-
//
// Matcher.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Matcher class.
//

#include "Matcher.h"
#include "StandardMatchers.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
#include "Matcher.tcc"
#endif



#define THEPEG_MATCH_DESC(T)                           \
/**                                                    \
 * This template specialization registers the Matcher  \
 */                                                    \
template <>                                            \
NoPIOClassDescription<T> T::initMatcher                \
 = NoPIOClassDescription<T>();                         \



namespace ThePEG {
  THEPEG_MATCH_DESC(MatchAny)
  THEPEG_MATCH_DESC(MatchStandardQCDParton)
  THEPEG_MATCH_DESC(MatchLightAntiQuark)
  THEPEG_MATCH_DESC(MatchLightQuark)
  THEPEG_MATCH_DESC(MatchLepton)
  THEPEG_MATCH_DESC(MatchDiquark)
  THEPEG_MATCH_DESC(MatchMeson)
  THEPEG_MATCH_DESC(MatchBaryon)
  THEPEG_MATCH_DESC(MatchNegative)
  THEPEG_MATCH_DESC(MatchNeutral)
  THEPEG_MATCH_DESC(MatchPositive)
  THEPEG_MATCH_DESC(MatchCharged)
  THEPEG_MATCH_DESC(MatchNeutrino)
}

namespace {

  using namespace ThePEG;
  static MatchAny m00;
  static MatchStandardQCDParton m01;
  static MatchLightAntiQuark m02;
  static MatchLightQuark m03;
  static MatchLepton m04;
  static MatchDiquark m05;
  static MatchMeson m06;
  static MatchBaryon m07;
  static MatchNegative m08;
  static MatchNeutral m09;
  static MatchPositive m11;
  static MatchCharged m12;
  static MatchNeutrino m13;

}
