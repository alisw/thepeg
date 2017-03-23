// -*- C++ -*-
//
// epsilon.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2011 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_epsilon_H
#define ThePEG_epsilon_H
//
// This is the declaration of the epsilon class.

#include "ThePEG/Vectors/LorentzVector.h"

namespace ThePEG {
namespace Helicity {

/** \ingroup Helicity
 *  \author Peter Richardson
 *
 *  This class is designed to combine 5-momenta and polarization 
 *  vectors together with the result being the product with the 
 *  eps function. The class is purely static and contains no data.
 *
 *  @see LorentzPolarizationVector
 *  @see Lorentz5Vector
 */

  /**
   *  Return the product 
   *  \f$\epsilon^{\mu\alpha\beta\gamma}v_{1\alpha}v_{2\beta}v_{3\gamma}\f$.
   * @param a The first  vector \f$v_{1\alpha}\f$.
   * @param b The second vector \f$v_{2\alpha}\f$.
   * @param c The third  vector \f$v_{3\alpha}\f$.
   * @return The product 
   * \f$\epsilon^{\mu\alpha\beta\gamma}v_{1\alpha}v_{2\beta}v_{3\gamma}\f$.
   */
  template <typename A, typename B, typename C>
  inline LorentzVector<typename BinaryOpTraits<
			 typename BinaryOpTraits<A,B>::MulT, C>::MulT>
  epsilon(const LorentzVector<A> & a,
	  const LorentzVector<B> & b,
	  const LorentzVector<C> & c) {
    typedef typename BinaryOpTraits<A,B>::MulT ABType;
    typedef typename BinaryOpTraits<ABType, C>::MulT ABCType;
    typedef LorentzVector<ABCType> ResultType;
    
    ABType diffxy = a.x() * b.y()  -  a.y() * b.x();
    ABType diffxz = a.x() * b.z()  -  a.z() * b.x();
    ABType diffxt = a.x() * b.t()  -  a.t() * b.x();
    ABType diffyz = a.y() * b.z()  -  a.z() * b.y();
    ABType diffyt = a.y() * b.t()  -  a.t() * b.y();
    ABType diffzt = a.z() * b.t()  -  a.t() * b.z();
    
    ResultType result;
    result.setX( c.z() * diffyt  - c.t() * diffyz  - c.y() * diffzt); 
    result.setY( c.t() * diffxz  - c.z() * diffxt  + c.x() * diffzt);
    result.setZ(-c.t() * diffxy  + c.y() * diffxt  - c.x() * diffyt);
    result.setT(-c.z() * diffxy  + c.y() * diffxz  - c.x() * diffyz);
    
    return result;
  }


}
}

#endif /* ThePEG_epsilon_H */
