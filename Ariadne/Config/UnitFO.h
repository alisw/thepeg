#ifndef ThePEG_UnitFO_H
#define ThePEG_UnitFO_H

#include <iostream>
#include "ThePEG/Vectors/Lorentz5Vector.h"

namespace ThePEG {

using namespace std;

/**
 * The FOUnit< class is used to facilitate output of unitful numbers
 * to an output stream. An Energy can hence be written like this:<BR>
 * <code>os << founit(x, GeV, width, precision);</code><BR> Also
 * containers of unitful numbers may be written like this, as well as
 * LorentzVector and ThreeVector.
 *
 * @see PersistentOStream
 * @see PersistentIStream
 * 
 */
template <typename T, typename UT>
struct FOUnit {

  /** Constructor given an object to be written assuming the given
   *  unit. */
  FOUnit(const T & t, const UT & u, int w, int p, ios::fmtflags f)
    : theX(t), theUnit(u), theWidth(w), thePrecision(p), theFFlags(f) {}

  /** Reference to the object to be written. */
  const T & theX;

  /** The unit assumed when writing the object. */
  const UT & theUnit;

  /**
   * The width to use for the printout.
   */
  int theWidth;

  /**
   * The precision to use for the printout.
   */
  int thePrecision;

  /**
   * The floating point format to use for the printout.
   */
  ios::fmtflags theFFlags;

};

/** Helper function writing out an object with a given unit to an
 *  output stream. */
template <typename T, typename UT>
void founitstream(ostream & os, const T & t, const UT & u,
		  int w, int p, ios::fmtflags f);

/** Helper function creating a FOUnit object given an object and a
 *  unit. */
template <typename T, typename UT>
inline FOUnit<T,UT>
founit(const T & t, const UT & ut,
       int w = -1, int p = -1, ios::fmtflags f = ios::fixed) {
  return FOUnit<T,UT>(t, ut, w, p, f);
}

template <typename T, typename UT>
ostream & operator<<(ostream & os, const FOUnit<T,UT> & u) {
  founitstream(os, u.theX, u.theUnit, u.theWidth, u.thePrecision, u.theFFlags);
  return os;
}

/** Output a ThreeVector with units to a stream. */
template <typename UnitT, typename Value>
void founitstream(ostream & os, const ThreeVector<Value> & v, const UnitT & u ,
		  int w, int p, ios::fmtflags f ) {
  os << founit(v.x(), u, w, p, f) << founit(v.y(), u, w, p, f)
     << founit(v.z(), u, w, p, f);
}

/** Output a LorentzVector with units to a stream. */
template <typename UnitT, typename Value>
void founitstream(ostream & os, const LorentzVector<Value> & v, const UnitT & u ,
		  int w, int p, ios::fmtflags f ) {
  os << founit(v.x(), u, w, p, f) << founit(v.y(), u, w, p, f)
     << founit(v.z(), u, w, p, f) << founit(v.e(), u, w, p, f);
}

/** Output a Lorentz5Vector with units to a stream. */
template <typename UnitT, typename Value>
void founitstream(ostream & os, const Lorentz5Vector<Value> & v, const UnitT & u ,
		  int w, int p, ios::fmtflags f ) {
  os << founit(v.x(), u, w, p, f) << founit(v.y(), u, w, p, f)
     << founit(v.z(), u, w, p, f) << founit(v.e(), u, w, p, f)
     << founit(v.mass(), u, w, p, f);
}

/** Helper function writing out an object with a given unit to an
 *  output stream. */
template <typename T, typename UT>
void founitstream(ostream & os, const T & t, const UT & u,
		  int w, int p, ios::fmtflags f) {
  ios::fmtflags saveflags = os.setf(f, ios::floatfield);
  int saveprec = os.precision();
  if ( p >= 0 )  os.precision(p);
  double x = t/u;
  if ( w >= 0 ) os.width(w);
  os << x;
  os.setf(saveflags);
  os.precision(saveprec);
}

}

#endif
