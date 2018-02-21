// -*- C++ -*-
//
// RhoDMatrix.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2017 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_RhoDMatrix_H
#define ThePEG_RhoDMatrix_H
// This is the declaration of the RhoDMatrix class.

#include "ThePEG/PDT/PDT.h"
#include "ThePEG/Helicity/HelicityDefinitions.h"
#include <cassert>
#include <array>

namespace ThePEG {

/**
 *  The RhoDMatrix class is designed to implement the storage of the
 *  rho and D matrices which are required for the spin correlation
 *  algorithm.  The matrix stores the spin as 2s+1.
 *
 * @author Peter Richardson
 *
 */
class RhoDMatrix {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor with undefined spin.
   */
  RhoDMatrix() = default;

  /**
   * Standard constructor giving the spin as 2s+1. The matrix starts out averaged, 
   * unless the second argument is false, when it is zeroed.
   */
  RhoDMatrix(PDT::Spin inspin, bool average = true) 
  : _spin(inspin), _ispin(abs(int(inspin))) {
    assert(_ispin <= MAXSPIN);
    // initialize to average
    if ( average )
        for(size_t ix=0; ix<_ispin; ++ix)
	    _matrix[ix][ix] = 1./_ispin;
  }
  //@}

public:

  /** @name Access matrix elements. */
  //@{
  /**
   * Return an element of the matrix.
   */
  Complex operator() (size_t ix, size_t iy) const {
    assert(ix < _ispin);
    assert(iy < _ispin);
    return _matrix[ix][iy];
  }

  /**
   * Set an element of the matrix.
   */
  Complex & operator() (size_t ix, size_t iy) {
    assert(ix < _ispin);
    assert(iy < _ispin);
    return _matrix[ix][iy];
  }

  /**
   * renormalise the matrix so it has unit trace
   */
  void normalize() {
#ifndef NDEBUG
    static const double epsa=1e-40, epsb=1e-10;
#endif
    Complex norm = 0.;
    for(size_t ix=0; ix<_ispin; ++ix) 
      norm += _matrix[ix][ix];
    assert(norm.real() > epsa);
    assert(norm.imag()/norm.real() < epsb);
    double invnorm = 1./norm.real();
    for(size_t ix=0; ix<_ispin; ++ix)
      for(size_t iy=0; iy<_ispin; ++iy) 
	_matrix[ix][iy]*=invnorm;
  }
  //@}

  /** @name Access the spin. */
  //@{

  /**
   * Get the spin. The spin is returned as 2J+1 in units of hbar/2.
   */
  PDT::Spin iSpin() const { return _spin; }
  //@}

  /**
   * Output the spin density matrix for debugging purposes.
   */
  friend ostream & operator<<(ostream & os, const RhoDMatrix & rd);

private:

  /**
   * 2s+1 for the particle.
   */
  PDT::Spin _spin;

  /**
   *  Storage of 2s+1 for speed.
   */
  size_t _ispin;

  /**
   * Spin matrix size
   */
  enum { MAXSPIN = 5 };

  /**
   * Storage for the matrix allowing up to spin 2 particles.
   */
  // Deliberately not using vector<> to avoid calls to 'new' 
  // from this commonly used class.
  std::array<std::array<Complex,MAXSPIN>,MAXSPIN> _matrix;

};

/** Output operator */
inline ostream & operator<<(ostream & os, const RhoDMatrix & rd) {
  for (size_t ix = 0; ix < rd._ispin; ++ix) {
    for (size_t iy = 0; iy < rd._ispin; ++iy)
      os << rd._matrix[ix][iy] << "  ";
    os << '\n';
  }
  return os;
}

}

#endif /* ThePEG_RhoDMatrix_H */
