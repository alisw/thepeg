// -*- C++ -*-
//
// epsilon.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2019 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_epsilon_H
#define ThePEG_epsilon_H
//
// This is the declaration of the epsilon class.

#include "ThePEG/Vectors/LorentzVector.h"
#include "LorentzTensor.h"

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
   *  \f$\epsilon^{\alpha\beta\gamma\delta}v_{1\alpha}v_{2\beta}v_{3\gamma}v_{4\delta}\f$.
   * @param a The first  vector \f$v_{1\alpha}\f$.
   * @param b The second vector \f$v_{2\beta}\f$.
   * @param c The third  vector \f$v_{3\gamma}\f$.
   * @param d The fourth vector \f$v_{4\delta}\f$.
   * @return The product 
   * \f$\epsilon^{\alpha\beta\gamma\delta}v_{1\alpha}v_{2\beta}v_{3\gamma}v_{4\delta}\f$
   */
  template <typename A, typename B, typename C, typename D>
  auto epsilon(const LorentzVector<A> & a,
               const LorentzVector<B> & b,
               const LorentzVector<C> & c,
               const LorentzVector<D> & d) 
    -> decltype(a.x()*b.y()*c.z()*d.t())
  {
    auto diffxy = a.x() * b.y()  -  a.y() * b.x();
    auto diffxz = a.x() * b.z()  -  a.z() * b.x();
    auto diffxt = a.x() * b.t()  -  a.t() * b.x();
    auto diffyz = a.y() * b.z()  -  a.z() * b.y();
    auto diffyt = a.y() * b.t()  -  a.t() * b.y();
    auto diffzt = a.z() * b.t()  -  a.t() * b.z();
    
    auto diff2xy = c.x() * d.y()  -  c.y() * d.x();
    auto diff2xz = c.x() * d.z()  -  c.z() * d.x();
    auto diff2xt = c.x() * d.t()  -  c.t() * d.x();
    auto diff2yz = c.y() * d.z()  -  c.z() * d.y();
    auto diff2yt = c.y() * d.t()  -  c.t() * d.y();
    auto diff2zt = c.z() * d.t()  -  c.t() * d.z();
 
    return
      diff2yz*diffxt + diff2zt*diffxy - diff2yt*diffxz - 
      diff2xz*diffyt + diff2xt*diffyz + diff2xy*diffzt;
  }
  
  /**
   *  Return the product 
   *  \f$\epsilon^{\mu\alpha\beta\gamma}v_{1\alpha}v_{2\beta}v_{3\gamma}\f$.
   * @param a The first  vector \f$v_{1\alpha}\f$.
   * @param b The second vector \f$v_{2\beta}\f$.
   * @param c The third  vector \f$v_{3\gamma}\f$.
   * @return The product 
   * \f$\epsilon^{\mu\alpha\beta\gamma}v_{1\alpha}v_{2\beta}v_{3\gamma}\f$.
   */
  template <typename A, typename B, typename C>
  auto epsilon(const LorentzVector<A> & a,
               const LorentzVector<B> & b,
               const LorentzVector<C> & c) 
  -> LorentzVector<decltype(a.x()*b.y()*c.z())>
  {
    auto diffxy = a.x() * b.y()  -  a.y() * b.x();
    auto diffxz = a.x() * b.z()  -  a.z() * b.x();
    auto diffxt = a.x() * b.t()  -  a.t() * b.x();
    auto diffyz = a.y() * b.z()  -  a.z() * b.y();
    auto diffyt = a.y() * b.t()  -  a.t() * b.y();
    auto diffzt = a.z() * b.t()  -  a.t() * b.z();

    using ResultType = LorentzVector<decltype(a.x()*b.x()*c.x())>;    
    ResultType result;
    result.setX( c.z() * diffyt  - c.t() * diffyz  - c.y() * diffzt); 
    result.setY( c.t() * diffxz  - c.z() * diffxt  + c.x() * diffzt);
    result.setZ(-c.t() * diffxy  + c.y() * diffxt  - c.x() * diffyt);
    result.setT(-c.z() * diffxy  + c.y() * diffxz  - c.x() * diffyz);
    
    return result;
  }
  
  /**
   *  Return the product 
   *  \f$\epsilon^{\mu\nu\alpha\beta}v_{1\alpha}v_{2\beta}\f$.
   * @param a The first  vector \f$v_{1\alpha}\f$.
   * @param b The second vector \f$v_{2\beta}\f$.
   * @return The product
   * \f$\epsilon^{\mu\nu\alpha\beta\gamma}v_{1\alpha}v_{2\beta}\f$.
   */
  template <typename A, typename B>
  auto epsilon(const LorentzVector<complex<A> > & a,
               const LorentzVector<complex<B> > & b)
    -> LorentzTensor<decltype(a.x().real()*b.y().real())>
  {
    auto diffxy = a.x() * b.y()  -  a.y() * b.x();
    auto diffxz = a.x() * b.z()  -  a.z() * b.x();
    auto diffxt = a.x() * b.t()  -  a.t() * b.x();
    auto diffyz = a.y() * b.z()  -  a.z() * b.y();
    auto diffyt = a.y() * b.t()  -  a.t() * b.y();
    auto diffzt = a.z() * b.t()  -  a.t() * b.z();
    complex<decltype(a.x().real()*b.x().real())> zero(ZERO);

    using ResultType = LorentzTensor<decltype(a.x().real()*b.x().real())>;
    ResultType result;
    result.setTT(  zero ); result.setTX( diffyz); result.setTY(-diffxz); result.setTZ( diffxy);
    result.setXT(-diffyz); result.setXX(  zero ); result.setXY( diffzt); result.setXZ(-diffyt);
    result.setYT( diffxz); result.setYX(-diffzt); result.setYY(  zero ); result.setYZ( diffxt);
    result.setZT(-diffxy); result.setZX( diffyt); result.setZY(-diffxt); result.setZZ(  zero );

    return result;
  }
  
  /**
   *  Return the product 
   *  \f$\epsilon^{\mu\nu\alpha\beta}v_{1\alpha}v_{2\beta}\f$.
   * @param a The first  vector \f$v_{1\alpha}\f$.
   * @param b The second vector \f$v_{2\beta}\f$.
   * @return The product
   * \f$\epsilon^{\mu\nu\alpha\beta\gamma}v_{1\alpha}v_{2\beta}\f$.
   */
  template <typename A, typename B>
  auto epsilon(const LorentzVector<A>           & a,
               const LorentzVector<complex<B> > & b)
    -> LorentzTensor<decltype(a.x()*b.y().real())>
  {
    auto diffxy = a.x() * b.y()  -  a.y() * b.x();
    auto diffxz = a.x() * b.z()  -  a.z() * b.x();
    auto diffxt = a.x() * b.t()  -  a.t() * b.x();
    auto diffyz = a.y() * b.z()  -  a.z() * b.y();
    auto diffyt = a.y() * b.t()  -  a.t() * b.y();
    auto diffzt = a.z() * b.t()  -  a.t() * b.z();
    complex<decltype(a.x()*b.x().real())> zero(ZERO);

    using ResultType = LorentzTensor<decltype(a.x()*b.x().real())>;
    ResultType result;
    result.setTT(  zero ); result.setTX( diffyz); result.setTY(-diffxz); result.setTZ( diffxy);
    result.setXT(-diffyz); result.setXX(  zero ); result.setXY( diffzt); result.setXZ(-diffyt);
    result.setYT( diffxz); result.setYX(-diffzt); result.setYY(  zero ); result.setYZ( diffxt);
    result.setZT(-diffxy); result.setZX( diffyt); result.setZY(-diffxt); result.setZZ(  zero );

    return result;
  }

  /**
   *  Return the product 
   *  \f$\epsilon^{\mu\nu\alpha\beta}v_{1\alpha}v_{2\beta}\f$.
   * @param a The first  vector \f$v_{1\alpha}\f$.
   * @param b The second vector \f$v_{2\beta}\f$.
   * @return The product
   * \f$\epsilon^{\mu\nu\alpha\beta\gamma}v_{1\alpha}v_{2\beta}\f$.
   */
  template <typename A, typename B>
  auto epsilon(const LorentzVector<complex<A> > & a,
               const LorentzVector<B>           & b)
    -> LorentzTensor<decltype(a.x().real()*b.y())>
  {
    auto diffxy = a.x() * b.y()  -  a.y() * b.x();
    auto diffxz = a.x() * b.z()  -  a.z() * b.x();
    auto diffxt = a.x() * b.t()  -  a.t() * b.x();
    auto diffyz = a.y() * b.z()  -  a.z() * b.y();
    auto diffyt = a.y() * b.t()  -  a.t() * b.y();
    auto diffzt = a.z() * b.t()  -  a.t() * b.z();
    complex<decltype(a.x().real()*b.x())> zero(ZERO);

    using ResultType = LorentzTensor<decltype(a.x().real()*b.x())>;
    ResultType result;
    result.setTT(  zero ); result.setTX( diffyz); result.setTY(-diffxz); result.setTZ( diffxy);
    result.setXT(-diffyz); result.setXX(  zero ); result.setXY( diffzt); result.setXZ(-diffyt);
    result.setYT( diffxz); result.setYX(-diffzt); result.setYY(  zero ); result.setYZ( diffxt);
    result.setZT(-diffxy); result.setZX( diffyt); result.setZY(-diffxt); result.setZZ(  zero );

    return result;
  }

  /**
   *  Return the product
   *  \f$\epsilon^{\mu\nu\alpha\beta}v_{1\alpha}v_{2\beta}\f$.
   * @param a The first  vector \f$v_{1\alpha}\f$.
   * @param b The second vector \f$v_{2\beta}\f$.
   * @return The product
   * \f$\epsilon^{\mu\nu\alpha\beta\gamma}v_{1\alpha}v_{2\beta}\f$.
   */
  template <typename A, typename B>
  auto epsilon(const LorentzVector<A> & a,
               const LorentzVector<B> & b)
    -> LorentzTensor<decltype(a.x()*b.y())>
  {
    auto diffxy = a.x() * b.y()  -  a.y() * b.x();
    auto diffxz = a.x() * b.z()  -  a.z() * b.x();
    auto diffxt = a.x() * b.t()  -  a.t() * b.x();
    auto diffyz = a.y() * b.z()  -  a.z() * b.y();
    auto diffyt = a.y() * b.t()  -  a.t() * b.y();
    auto diffzt = a.z() * b.t()  -  a.t() * b.z();
    complex<decltype(a.x()*b.x())> zero(ZERO);

    using ResultType = LorentzTensor<decltype(a.x()*b.x())>;
    ResultType result;
    result.setTT(  zero ); result.setTX( diffyz); result.setTY(-diffxz); result.setTZ( diffxy);
    result.setXT(-diffyz); result.setXX(  zero ); result.setXY( diffzt); result.setXZ(-diffyt);
    result.setYT( diffxz); result.setYX(-diffzt); result.setYY(  zero ); result.setYZ( diffxt);
    result.setZT(-diffxy); result.setZX( diffyt); result.setZY(-diffxt); result.setZZ(  zero );

    return result;
  }


}
}

#endif /* ThePEG_epsilon_H */
