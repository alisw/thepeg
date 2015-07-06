// -*- C++ -*-
//
// StandardMatchers.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_StandardMatchers_H
#define ThePEG_StandardMatchers_H
// This is the declaration of the AnyMatcher,


#include "Matcher.h"
#include "ThePEG/PDT/EnumParticles.h"

namespace ThePEG {

/** \file StandardMatchers.h
 *
 * This file declare a set of standard matcher classes. The
 * ChargedMatcher, NegativeMatcher, PositiveMatcher, NeutralMatcher,
 * BaryonMatcher, MesonMatcher, DiquarkMatcher, LeptonMatcher,
 * LightAntiQuarkMatcher, LightQuarkMatcher and
 * StandardQCDPartonMatcher classes can be used by themselves (with
 * their static functions) or together with the Matcher class to
 * define Interfaced objects of the MatcherBase type to be used in the
 * Repository. Suitable typedefs are declared for the latter.
 *
 * @see Matcher
 * @see MatcherBase
 */

/**
 * A Matcher class which matches any particle.
 */
struct AnyMatcher: public MatcherType {
  /** Typedef the class matching the complex conjugate particles. */
  typedef AnyMatcher CC;
  /** The main static function to check if a given particle type \a pd
      matches. */
  static bool Check(const ParticleData & pd) {
    return bool(pd.id());
  }
  /** A simplified but unique class name. */
  static string className() { return "Any"; }
};
/** Gives a MatcherBase class based on AnyMatcher. */
typedef Matcher<AnyMatcher> MatchAny;


/**
 * A Matcher class which matches any charged particle.
 */
struct ChargedMatcher: public MatcherType {
  /** Typedef the class matching the complex conjugate particles. */
  typedef ChargedMatcher CC;
  /** The main static function to check if a given particle type \a pd
      matches. */
  static bool Check(const ParticleData & pd) {
    return PDT::charged(pd.iCharge());
  }
  /** A simplified but unique class name. */
  static string className() { return "Charged"; }
};
/** Gives a MatcherBase class based on ChargedMatcher. */
typedef Matcher<ChargedMatcher> MatchCharged;

struct NegativeMatcher;

/**
 * A Matcher class which matches any positively charged particle.
 */
struct PositiveMatcher: public MatcherType {
  /** Typedef the class matching the complex conjugate particles. */
  typedef NegativeMatcher CC;
  /** The main static function to check if a given particle type \a pd
      matches. */
  static bool Check(const ParticleData & pd) {
    return PDT::positive(pd.iCharge());
  }
  /** A simplified but unique class name. */
  static string className() { return "Positive"; }
};
/** Gives a MatcherBase class based on PositiveMatcher. */
typedef Matcher<PositiveMatcher> MatchPositive;


/**
 * A Matcher class which matches any uncharged particle.
 */
struct NeutralMatcher: public MatcherType {
  /** Typedef the class matching the complex conjugate particles. */
  typedef NeutralMatcher CC;
  /** The main static function to check if a given particle type \a pd
      matches. */
  static bool Check(const ParticleData & pd) {
    return pd.iCharge() == PDT::Charge0;
  }
  /** A simplified but unique class name. */
  static string className() { return "Neutral"; }
};
/** Gives a MatcherBase class based on NeutralMatcher. */
typedef Matcher<NeutralMatcher> MatchNeutral;

/**
 * A Matcher class which matches any negatively charged particle.
 */
struct NegativeMatcher: public MatcherType {
  /** Typedef the class matching the complex conjugate particles. */
  typedef PositiveMatcher CC;
  /** The main static function to check if a given particle type \a pd
      matches. */
  static bool Check(const ParticleData & pd) {
    return PDT::negative(pd.iCharge());
  }
  /** A simplified but unique class name. */
  static string className() { return "Negative"; }
};
/** Gives a MatcherBase class based on NegativeMatcher. */
typedef Matcher<NegativeMatcher> MatchNegative;

/**
 * A Matcher class which matches any baryon.
 */
struct BaryonMatcher: public MatcherType {
  /** Typedef the class matching the complex conjugate particles. */
  typedef BaryonMatcher CC;
  /** The main static function to check if a given particle type \a pd
      matches. */
  static bool Check(const ParticleData & pd) { return Check(pd.id()); }
  /** The main static function to check if a given particle with type
      \a id matches. */
  static bool Check(long id) {
    return (id/10)%10 && (id/100)%10 && (id/1000)%10;
  }
  /** A simplified but unique class name. */
  static string className() { return "Baryon"; }
};
/** Gives a MatcherBase class based on BaryonMatcher. */
typedef Matcher<BaryonMatcher> MatchBaryon;


/**
 * A Matcher class which matches any meson.
 */
struct MesonMatcher: public MatcherType {
  /** Typedef the class matching the complex conjugate particles. */
  typedef MesonMatcher CC;
  /** The main static function to check if a given particle type \a pd
      matches. */
  static bool Check(const ParticleData & pd) { return Check(pd.id()); }
  /** The main static function to check if a given particle with type
      \a id matches. */
  static bool Check(long id) {
    return (id/10)%10 && (id/100)%10 && (id/1000)%10 == 0;
  }
  /** A simplified but unique class name. */
  static string className() { return "Meson"; }
};
/** Gives a MatcherBase class based on MesonMatcher. */
typedef Matcher<MesonMatcher> MatchMeson;


/**
 * A Matcher class which matches any (anti-)diquark.
 */
struct DiquarkMatcher: public MatcherType {
  /** Typedef the class matching the complex conjugate particles. */
  typedef DiquarkMatcher CC;
  /** The main static function to check if a given particle type \a pd
      matches. */
  static bool Check(const ParticleData & pd) { return Check(pd.id()); }
  /** The main static function to check if a given particle with type
      \a id matches. */
  static bool Check(long id) {
    return id/10 && (id/10)%10 == 0 && (id/100)%10 && (id/1000)%10;
  }
  /** A simplified but unique class name. */
  static string className() { return "Diquark"; }
};
/** Gives a MatcherBase class based on DiquarkMatcher. */
typedef Matcher<DiquarkMatcher> MatchDiquark;

/**
 * A Matcher class which matches any (anti-)quark.
 */
struct QuarkMatcher: public MatcherType {
  /** Typedef the class matching the complex conjugate particles. */
  typedef QuarkMatcher CC;
  /** The main static function to check if a given particle type \a pd
      matches. */
  static bool Check(const ParticleData & pd) { return Check(pd.id()); }
  /** The main static function to check if a given particle with type
      \a id matches. */
  static bool Check(long id) {
    return id && abs(id) < 10;
  }
  /** A simplified but unique class name. */
  static string className() { return "Quark"; }
};
/** Gives a MatcherBase class based on QuarkMatcher. */
typedef Matcher<QuarkMatcher> MatchQuark;

/**
 * A Matcher class which matches any lepton.
 */
struct LeptonMatcher: public MatcherType {
  /** Typedef the class matching the complex conjugate particles. */
  typedef LeptonMatcher CC;
  /** The main static function to check if a given particle type \a pd
      matches. */
  static bool Check(const ParticleData & pd) { return Check(pd.id()); }
  /** The main static function to check if a given particle with type
      \a id matches. */
  static bool Check(long id) {
    return abs(id) > 10 && abs(id) <= 20;
  }
  /** A simplified but unique class name. */
  static string className() { return "Lepton"; }
};
/** Gives a MatcherBase class based on LeptonMatcher. */
typedef Matcher<LeptonMatcher> MatchLepton;

struct LightAntiQuarkMatcher;

/**
 * A Matcher class which matches any light quark (d,u or s).
 */
struct LightQuarkMatcher: public MatcherType {
  /** Typedef the class matching the complex conjugate particles. */
  typedef LightAntiQuarkMatcher CC;
  /** The main static function to check if a given particle type \a pd
      matches. */
  static bool Check(const ParticleData & pd) { return Check(pd.id()); }
  /** The main static function to check if a given particle with type
      \a id matches. */
  static bool Check(long id) {
    return id > 0 && id < 4 ;
  }
  /** A simplified but unique class name. */
  static string className() { return "LightQuark"; }
};
/** Gives a MatcherBase class based on LightQuarkMatcher. */
typedef Matcher<LightQuarkMatcher> MatchLightQuark;


/**
 * A Matcher class which matches any light anti-quark
 * (\f$\bar{\mbox{d}}\f$,\f$\bar{\mbox{u}}\f$ or
 * \f$\bar{\mbox{s}}\f$).
 */
struct LightAntiQuarkMatcher: public MatcherType {
  /** Typedef the class matching the complex conjugate particles. */
  typedef LightQuarkMatcher CC;
  /** The main static function to check if a given particle type \a pd
      matches. */
  static bool Check(const ParticleData & pd) { return Check(pd.id()); }
  /** The main static function to check if a given particle with type
      \a id matches. */
  static bool Check(long id) {
    return id < 0 && id > -4 ;
  }
  /** A simplified but unique class name. */
  static string className() { return "LightAntiQuark"; }
};
/** Gives a MatcherBase class based on LightAntiQuarkMatcher. */
typedef Matcher<LightAntiQuarkMatcher> MatchLightAntiQuark;

/**
 * A Matcher class which matches any standard QCD parton, ie. gluons
 * and quarks up to bottom.
 */
struct StandardQCDPartonMatcher: public MatcherType {
  /** Typedef the class matching the complex conjugate particles. */
  typedef StandardQCDPartonMatcher CC;
  /** The main static function to check if a given particle type \a pd
      matches. */
  static bool Check(const ParticleData & pd) { return Check(pd.id()); }
  /** The main static function to check if a given particle with type
      \a id matches. */
  static bool Check(long id) {
    return id && ( abs(id) <= 5 || id == ParticleID::g );
  }
  /** A simplified but unique class name. */
  static string className() { return "StandardQCDParton"; }
};

/** Gives a MatcherBase class based on StandardQCDPartonMatcher. */
typedef Matcher<StandardQCDPartonMatcher> MatchStandardQCDParton;

/**
 * A Matcher class which matches any pseudo scalar meson.
 */
struct PseudoScalarMesonMatcher: public MatcherType {
  /** Typedef the class matching the complex conjugate particles. */
  typedef PseudoScalarMesonMatcher CC;
  /** The main static function to check if a given particle type \a pd
      matches. */
  static bool Check(const ParticleData & pd) { return Check(pd.id()); }
  /** The main static function to check if a given particle with type
      \a id matches. */
  static bool Check(long id) {
    return ( (abs(id)/1000)%1 == 0 && abs(id) > 100 && abs(id)%10 == 1 ) ||
      ( id == ParticleID::K_L0 || id == ParticleID::K_S0 );
  }
  /** A simplified but unique class name. */
  static string className() { return "PseudoScalarMeson"; }
};

/** Gives a MatcherBase class based on PseudoScalarMesonMatcher. */
typedef Matcher<PseudoScalarMesonMatcher> MatchPseudoScalarMeson;


/**
 * A Matcher class which matches any vector meson.
 */
struct VectorMesonMatcher: public MatcherType {
  /** Typedef the class matching the complex conjugate particles. */
  typedef VectorMesonMatcher CC;
  /** The main static function to check if a given particle type \a pd
      matches. */
  static bool Check(const ParticleData & pd) { return Check(pd.id()); }
  /** The main static function to check if a given particle with type
      \a id matches. */
  static bool Check(long id) {
    return (abs(id)/1000)%1 == 0 && abs(id) > 100 && abs(id)%10 == 3;
  }
  /** A simplified but unique class name. */
  static string className() { return "VectorMeson"; }
};

/** Gives a MatcherBase class based on VectorMesonMatcher. */
typedef Matcher<VectorMesonMatcher> MatchVectorMeson;

}

#endif /* ThePEG_StandardMatchers_H */
