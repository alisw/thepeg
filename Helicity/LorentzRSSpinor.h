// -*- C++ -*-
//
// LorentzRSSpinor.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2019 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_LorentzRSSpinor_H
#define ThePEG_LorentzRSSpinor_H

#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/Vectors/ThreeVector.h"
#include "HelicityDefinitions.h"
#include "LorentzRSSpinor.fh"
#include "LorentzRSSpinorBar.h"
#include "LorentzSpinorBar.h"
#include "LorentzSpinor.h"
#include "LorentzPolarizationVector.h"

namespace ThePEG{
namespace Helicity{

/**
 *  The LorentzRSSpinor class is designed to store a Rarita-Schwinger
 *  spinor for a spin-3/2 particle. In addition to storing the
 *  components of the spinor information is stored on the type of
 *  spinor, for example u or v type.
 *
 *  At the moment only one choice of the Dirac matrix representation
 *  is supported. For high-energy calculations the choice made by the
 *  HELAS collaboration is more efficient for numerical
 *  calculations. In this representation
 *
 *  \f[
 * \gamma_{i=1,2,3}=\left(\begin{array}{cc}
 *                          0 & \sigma_i \\
 *                          -\sigma_i & 0
 *                        \end{array}\right)
 *          \quad
 * \gamma_0=\left(\begin{array}{cc}
 *                  0 & 1 \\
 *                  1 & 0
 *                \end{array}\right)
 *          \quad
 * \gamma_5=\left(\begin{array}{cc}
 *                  -1 & 0 \\
 *                  0 & 1
 *                \end{array}\right)
 * \f]
 *
 *  The type of the spinor is also stored using the SpinorType
 *  enumeration.  There are three types supported SpinorType::u,
 *  SpinorType::v, SpinorType::unknown.  This information is intended
 *  mainly for use in the case of Majorana particles where matrix
 *  elements can be calculated with either u or v type spinors and
 *  knowledge of which was used will be needed in order to give the
 *  correct correlations. The SpinorType::unknowne is intended for
 *  cases where either the spinor for an off-shell line in a matrix
 *  element calculation or the information is genuinely unknown.
 *
 *  The LorentzRSSpinorBar class is also provided to store the barred
 *  spinor.
 *
 *
 * @see HelicityDefinitions
 * @see LorentzRSSpinorBar
 *
 * \author Peter Richardson
 *
 */
template<typename Value>
class LorentzRSSpinor {

public:

  /** @name Standard constructors. */
  //@{
  /**
   * Default zero constructor, optionally specifying \a t, the type.
   */
  LorentzRSSpinor(SpinorType t = SpinorType::unknown) : _type(t), _spin() {}

  /**
   * Constructor with complex numbers specifying the components,
   * optionally specifying \a t, the type.
   */
  LorentzRSSpinor(complex<Value> a1, complex<Value> b1,
		  complex<Value> c1, complex<Value> d1,
		  complex<Value> a2, complex<Value> b2,
		  complex<Value> c2, complex<Value> d2,
		  complex<Value> a3, complex<Value> b3,
		  complex<Value> c3, complex<Value> d3,
		  complex<Value> a4, complex<Value> b4,
		  complex<Value> c4, complex<Value> d4,
		  SpinorType t=SpinorType::unknown) 
    : _type(t), _spin{{ {{a1,b1,c1,d1}},
                        {{a2,b2,c2,d2}},
                        {{a3,b3,c3,d3}},
                        {{a4,b4,c4,d4}}
                      }} {}
  
  template <typename U>
  LorentzRSSpinor(const LorentzRSSpinor<U> & other)
    : _type(other._type), _spin(other._spin) {}
  //@}

  /** @name Access the components. */
  //@{
  /**
   * Subscript operator to return spinor components
   */
  complex<Value> operator()(int i, int j) const {
    assert( i >= 0 && i <= 3 && j>=0 && j<=3);
    return _spin[i][j];
  }

  /**
   * Set components by index
   */
  complex<Value> & operator () (int i, int j) {
    assert( i >= 0 && i <= 3 && j>=0 && j<=3);
    return _spin[i][j];
  }

  /**
   * Get first spinor component for the x vector
   */
  complex<Value> xs1() const {return _spin[0][0];}

  /**
   * Get second spinor component for the x vector
   */
  complex<Value> xs2() const {return _spin[0][1];}

  /**
   * Get third  spinor component for the x vector
   */
  complex<Value> xs3() const {return _spin[0][2];}

  /**
   * Get fourth  spinor component for the x vector
   */
  complex<Value> xs4() const {return _spin[0][3];}

  /**
   * Get first spinor component for the y vector
   */
  complex<Value> ys1() const {return _spin[1][0];}

  /**
   * Get second spinor component for the y vector
   */
  complex<Value> ys2() const {return _spin[1][1];}
  
  /**
   * Get third spinor component for the y vector
   */
  complex<Value> ys3() const {return _spin[1][2];}
  
  /**
   * Get fourth spinor component for the y vector
   */
  complex<Value> ys4() const {return _spin[1][3];}
  
  /**
   * Get first spinor component for the z vector
   */
  complex<Value> zs1() const {return _spin[2][0];}
  
  /**
   * Get second spinor component for the z vector
   */
  complex<Value> zs2() const {return _spin[2][1];}
  
  /**
   * Get third spinor component for the z vector
   */
  complex<Value> zs3() const {return _spin[2][2];}
  
  /**
   * Get fourth spinor component for the z vector
   */
  complex<Value> zs4() const {return _spin[2][3];}
  
  /**
   * Get first spinor component for the t vector
   */
  complex<Value> ts1() const {return _spin[3][0];}
  
  /**
   * Get second spinor component for the t vector
   */
  complex<Value> ts2() const {return _spin[3][1];}
  
  /**
   * Get third spinor component for the t vector
   */
  complex<Value> ts3() const {return _spin[3][2];}
  
  /**
   * Get fourth spinor component for the t vector
   */
  complex<Value> ts4() const {return _spin[3][3];}
  
  /**
   * Set first spinor component for the x vector
   */
  void setXS1(complex<Value> in) {_spin[0][0]=in;}
  
  /**
   * Set second spinor component for the x vector
   */
  void setXS2(complex<Value> in) {_spin[0][1]=in;}
  
  /**
   * Set third spinor component for the x vector
   */
  void setXS3(complex<Value> in) {_spin[0][2]=in;}
  
  /**
   * Set fourth spinor component for the x vector
   */
  void setXS4(complex<Value> in) {_spin[0][3]=in;}
  
  /**
   * Set first spinor component for the y vector
   */
  void setYS1(complex<Value> in) {_spin[1][0]=in;}
  
  /**
   * Set second spinor component for the y vector
   */
  void setYS2(complex<Value> in) {_spin[1][1]=in;}
  
  /**
   * Set third spinor component for the y vector
   */
  void setYS3(complex<Value> in) {_spin[1][2]=in;}
  
  /**
   * Set fourth spinor component for the y vector
   */
  void setYS4(complex<Value> in) {_spin[1][3]=in;}
  
  /**
   * Set first spinor component for the z vector
   */
  void setZS1(complex<Value> in) {_spin[2][0]=in;}
  
  /**
   * Set second spinor component for the z vector
   */
  void setZS2(complex<Value> in) {_spin[2][1]=in;}
  
  /**
   * Set third spinor component for the z vector
   */
  void setZS3(complex<Value> in) {_spin[2][2]=in;}
  
  /**
   * Set fourth spinor component for the z vector
   */
  void setZS4(complex<Value> in) {_spin[2][3]=in;}
  
  /**
   * Set first spinor component for the t vector
   */
  void setTS1(complex<Value> in) {_spin[3][0]=in;}
  
  /**
   * Set second spinor component for the t vector
   */
  void setTS2(complex<Value> in) {_spin[3][1]=in;}
  
  /**
   * Set third spinor component for the t vector
   */
  void setTS3(complex<Value> in) {_spin[3][2]=in;}
  
  /**
   * Set fourth spinor component for the t vector
   */
  void setTS4(complex<Value> in) {_spin[3][3]=in;}
  //@}

  /// @name Mathematical assignment operators.
  //@{
  template <typename ValueB>
  LorentzRSSpinor<Value> & operator+=(const LorentzRSSpinor<ValueB> & a) {
    for(unsigned int ix=0;ix<4;++ix)
      for(unsigned int iy=0;iy<4;++iy)
	_spin[ix][iy] += a._spin[ix][iy];
    return *this;
  }
  
  template <typename ValueB>
  LorentzRSSpinor<Value> & operator-=(const LorentzRSSpinor<ValueB> & a) {
    for(unsigned int ix=0;ix<4;++ix)
      for(unsigned int iy=0;iy<4;++iy)
	_spin[ix][iy] -= a._spin[ix][iy];
    return *this;
  }
  
  LorentzRSSpinor<Value> & operator*=(double a) {
    for(unsigned int ix=0;ix<4;++ix)
      for(unsigned int iy=0;iy<4;++iy)
	_spin[ix][iy] *=a;
    return *this;
  }
  
  LorentzRSSpinor<Value> & operator/=(double a) {
    for(unsigned int ix=0;ix<4;++ix)
      for(unsigned int iy=0;iy<4;++iy)
	_spin[ix][iy] /=a;
    return *this;
  }
  //@}
  
  /** @name Arithmetic operators. */
  //@{
  /**
   * dot product with a polarization vector
   */
  LorentzSpinor<Value> dot(const LorentzPolarizationVector & vec) const {
    LorentzSpinor<Value> output(_type);
    complex<Value> temp;
    unsigned int ix;
    for(ix=0;ix<4;++ix) {
      temp  = _spin[3][ix]*vec.t();
      temp -= _spin[0][ix]*vec.x();
      temp -= _spin[1][ix]*vec.y();
      temp -= _spin[2][ix]*vec.z();
      output[ix]=temp;
    }
    return output;
  }

  /**
   * dot product with a 4-vector
   */
  LorentzSpinor<Value> dot(const LorentzMomentum & invec) const {
    LorentzSpinor<Value> output(_type);
    complex<Value> temp;
    LorentzVector<double> vec = UnitRemoval::InvE * invec;
    unsigned int ix;
    for(ix=0;ix<4;++ix) {
      temp  = - ( _spin[0][ix]*vec.x() + _spin[1][ix]*vec.y()+
		  _spin[2][ix]*vec.z() ) +  _spin[3][ix]*vec.t();
      output[ix]=temp;
    }
    return output;
  }
  //@}

  /** @name Transformations. */
  //@{
  /**
   * return the barred spinor
   */
  LorentzRSSpinorBar<Value> bar() const;

  /**
   * Standard Lorentz boost specifying the components of the beta vector.
   */
  LorentzRSSpinor & boost(double,double,double);

  /**
   * Standard Lorentz boost specifying the beta vector.
   */
  LorentzRSSpinor & boost(const Boost &);

  /**
   * General transform
   */
  LorentzRSSpinor & transform(const LorentzRotation &);
  //@}

  /** @name Functions related to type. */
  //@{

  /**
   * Return the type of the spinor.
   */
  SpinorType Type() const {return _type;}
  //@}

  /**
   * Scalar product \f$\bar{f}^\alpha(c_LP_L+c_RP_R)f_\alpha\f$for general couplings
   * @param fbar The barred spinor
   * @param left The left-handed coupling, \f$c_L\f$.
   * @param right The right-handed coupling, \f$c_R\f$.
   */
  template <typename ValueB>
  auto generalScalar(LorentzRSSpinorBar<ValueB>& fbar, 
                     Complex left, Complex right) 
  -> decltype(left*fbar(3,0)*this->ts1())
  {
    decltype(left*fbar(3,0)*ts1()) output; 
    unsigned int iz;
    output = 
      left*(fbar(3,0)*_spin[3][0]+fbar(3,1)*_spin[3][1])
      +right*(fbar(3,2)*_spin[3][2]+fbar(3,3)*_spin[3][3]);
    for(iz=0;iz<3;++iz) {
      output -=
	left*(fbar(iz,0)*_spin[iz][0]+fbar(iz,1)*_spin[iz][1])
	+right*(fbar(iz,2)*_spin[iz][2]+fbar(iz,3)*_spin[iz][3]);
    }
    return output;
  }
  
  /**
   *  Current \f$\bar{f}(c_LP_L+c_RP_R)f^\alpha\f$ for general couplings.
   * @param fbar The barred spinor
   * @param left The left-handed coupling, \f$c_L\f$.
   * @param right The right-handed coupling, \f$c_R\f$.
   */
  template <typename ValueB>
  auto generalCurrent(LorentzSpinorBar<ValueB>& fbar, 
                      Complex left, Complex right) 
  -> LorentzVector<decltype(left*fbar.s1()*this->ts1())>
  {
    typedef decltype(left*fbar.s1()*ts1()) ResultT;
    ResultT output[4];
    for(size_t iz=0;iz<4;++iz)
      output[iz]= left*(fbar.s1()*_spin[iz][0]+fbar.s2()*_spin[iz][1])
	+right*(fbar.s3()*_spin[iz][2]+fbar.s4()*_spin[iz][3]);
    return LorentzVector<ResultT>(output[0],output[1],output[2],output[3]);
  }
  
private:

  /**
   * Type of spinor
   */
  SpinorType _type;

  /**
   * Storage of the components.
   */
  std::array<std::array<complex<Value>,4>,4> _spin;
};

/// @name Basic mathematical operations
//@{
template <typename Value>
inline LorentzRSSpinor<double>
operator/(const LorentzRSSpinor<Value> & v, Value a) {
  return LorentzRSSpinor<double>(v.xs1()/a, v.xs2()/a, v.xs3()/a, v.xs4()/a,
			       v.ys1()/a, v.ys2()/a, v.ys3()/a, v.ys4()/a,
			       v.zs1()/a, v.zs2()/a, v.zs3()/a, v.zs4()/a,
			       v.ts1()/a, v.ts2()/a, v.ts3()/a, v.ts4()/a,
			       v.Type());
}

inline LorentzRSSpinor<double>
operator/(const LorentzRSSpinor<double> & v, Complex a) {
  return LorentzRSSpinor<double>(v.xs1()/a, v.xs2()/a, v.xs3()/a, v.xs4()/a,
				 v.ys1()/a, v.ys2()/a, v.ys3()/a, v.ys4()/a,
				 v.zs1()/a, v.zs2()/a, v.zs3()/a, v.zs4()/a,
				 v.ts1()/a, v.ts2()/a, v.ts3()/a, v.ts4()/a,
				 v.Type());
}

template <typename Value>
inline LorentzRSSpinor<Value> operator-(const LorentzRSSpinor<Value> & v) {
  return LorentzRSSpinor<Value>(-v.xs1(),-v.xs2(),-v.xs3(),-v.xs4(),
				-v.ys1(),-v.ys2(),-v.ys3(),-v.ys4(),
				-v.zs1(),-v.zs2(),-v.zs3(),-v.zs4(),
				-v.ts1(),-v.ts2(),-v.ts3(),-v.ts4(),
				v.Type());
}

template <typename ValueA, typename ValueB>
inline LorentzRSSpinor<ValueA>
operator+(LorentzRSSpinor<ValueA> a, const LorentzRSSpinor<ValueB> & b) {
  return a += b;
}

template <typename ValueA, typename ValueB>
inline LorentzRSSpinor<ValueA>
operator-(LorentzRSSpinor<ValueA> a, const LorentzRSSpinor<ValueB> & b) {
  return a -= b;
}

template <typename Value>
inline LorentzRSSpinor<Value>
operator*(const LorentzRSSpinor<Value> & a, double b) {
  return LorentzRSSpinor<Value>(a.xs1()*b, a.xs2()*b, a.xs3()*b, a.xs4()*b,
				a.ys1()*b, a.ys2()*b, a.ys3()*b, a.ys4()*b,
				a.zs1()*b, a.zs2()*b, a.zs3()*b, a.zs4()*b,
				a.ts1()*b, a.ts2()*b, a.ts3()*b, a.ts4()*b,a.Type());
}

template <typename Value>
inline LorentzRSSpinor<Value>
operator*(double b, LorentzRSSpinor<Value> a) {
  return a *= b;
}
  
template <typename Value>
inline LorentzRSSpinor<Value>
operator*(const LorentzRSSpinor<Value> & a, Complex b) {
  return LorentzRSSpinor<Value>(a.xs1()*b, a.xs2()*b, a.xs3()*b, a.xs4()*b,
				a.ys1()*b, a.ys2()*b, a.ys3()*b, a.ys4()*b,
				a.zs1()*b, a.zs2()*b, a.zs3()*b, a.zs4()*b,
				a.ts1()*b, a.ts2()*b, a.ts3()*b, a.ts4()*b,a.Type());
}

template <typename ValueA, typename ValueB>
inline auto operator*(complex<ValueB> a, const LorentzRSSpinor<ValueA> & v) 
  -> LorentzRSSpinor<decltype(a.real()*v.xs1().real())>
{
  return {a*v.xs1(), a*v.xs2(), a*v.xs3(), a*v.xs4(),
          a*v.ys1(), a*v.ys2(), a*v.ys3(), a*v.ys4(),
          a*v.zs1(), a*v.zs2(), a*v.zs3(), a*v.zs4(),
          a*v.ts1(), a*v.ts2(), a*v.ts3(), a*v.ts4(),v.Type()};
}

template <typename ValueA, typename ValueB>
inline auto operator*(const LorentzRSSpinor<ValueA> & v, complex<ValueB> b) 
  -> LorentzRSSpinor<decltype(b.real()*v.xs1().real())>
{
  return b*v;
}

template <typename ValueA, typename ValueB>
inline auto operator/(const LorentzRSSpinor<ValueA> & v, complex<ValueB> b) 
  -> LorentzRSSpinor<decltype(v.xs1().real()/b.real())>
{
  return {v.xs1()/b, v.xs2()/b, v.xs3()/b, v.xs4()/b,
          v.ys1()/b, v.ys2()/b, v.ys3()/b, v.ys4()/b,
          v.zs1()/b, v.zs2()/b, v.zs3()/b, v.zs4()/b,
          v.ts1()/b, v.ts2()/b, v.ts3()/b, v.ts4()/b,v.Type()};
}
  
template <typename ValueA, typename ValueB>
inline auto operator*(ValueB a, const LorentzRSSpinor<ValueA> & v) 
  -> LorentzRSSpinor<decltype(a*v.xs1().real())>
{
  return {a*v.xs1(), a*v.xs2(), a*v.xs3(), a*v.xs4(),
          a*v.ys1(), a*v.ys2(), a*v.ys3(), a*v.ys4(),
          a*v.zs1(), a*v.zs2(), a*v.zs3(), a*v.zs4(),
          a*v.ts1(), a*v.ts2(), a*v.ts3(), a*v.ts4(),v.Type()};
}

template <typename ValueA, typename ValueB>
inline auto operator*(const LorentzRSSpinor<ValueA> & v, ValueB b) 
  -> LorentzRSSpinor<decltype(b*v.xs1().real())>
{
  return b*v;
}

template <typename ValueA, typename ValueB>
inline auto operator/(const LorentzRSSpinor<ValueA> & v, ValueB b) 
  -> LorentzRSSpinor<decltype(v.xs1().real()/b)>
{
  return {v.xs1()/b, v.xs2()/b, v.xs3()/b, v.xs4()/b,
          v.ys1()/b, v.ys2()/b, v.ys3()/b, v.ys4()/b,
          v.zs1()/b, v.zs2()/b, v.zs3()/b, v.zs4()/b,
          v.ts1()/b, v.ts2()/b, v.ts3()/b, v.ts4()/b,v.Type()};
}
//@}

}
}
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
#include "LorentzRSSpinor.tcc"
#endif 

#endif
