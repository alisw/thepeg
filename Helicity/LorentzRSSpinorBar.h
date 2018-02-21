// -*- C++ -*-
//
// LorentzRSSpinorBar.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2017 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_LorentzRSSpinorBar_H
#define ThePEG_LorentzRSSpinorBar_H
// This is the declaration of the LorentzRSSpinorBar class.

#include "ThePEG/Config/ThePEG.h"
#include "ThePEG/Vectors/ThreeVector.h"
#include "HelicityDefinitions.h"
#include "LorentzRSSpinor.fh"
#include "LorentzRSSpinorBar.fh"
#include "LorentzSpinorBar.h"
#include "LorentzSpinor.h"
#include "LorentzPolarizationVector.h"

namespace ThePEG {
namespace Helicity {

/**
 *  The <code>LorentzRSSpinorBar</code> class implements the storage of a
 *  barred Lorentz Rarita-Schwinger Spinor for a spin-3/2 particle.
 *  The design is based on that of the
 *  <code>LorentzRSSpinor</code> class and the details of the implemented
 *  are discussed in more detail in the header file for that class.
 *
 * @see LorentzRSSpinor
 *
 * \author Peter Richardson
 *
 */
template <typename Value>
class LorentzRSSpinorBar {

public:

  /** @name Standard constructors. */
  //@{
  /**
   * Default zero constructor, optionally specifying \a t, the type.
   */
  LorentzRSSpinorBar(SpinorType t = SpinorType::unknown) : _type(t), _spin() {}

  /**
   * Constructor with complex numbers specifying the components,
   * optionally specifying \a t, the type.
   */
  LorentzRSSpinorBar(complex<Value> a1,complex<Value> b1,
		     complex<Value> c1,complex<Value> d1,
		     complex<Value> a2,complex<Value> b2,
		     complex<Value> c2,complex<Value> d2,
		     complex<Value> a3,complex<Value> b3,
		     complex<Value> c3,complex<Value> d3,
		     complex<Value> a4,complex<Value> b4,
		     complex<Value> c4,complex<Value> d4,
		     SpinorType t=SpinorType::unknown)
    : _type(t), _spin{{ {{a1,b1,c1,d1}},
                        {{a2,b2,c2,d2}},
                        {{a3,b3,c3,d3}},
                        {{a4,b4,c4,d4}}
                      }} {}
  //@}

  /** @name Access the components. */
  //@{
  /**
   * Subscript operator to return spinor components
   */
  complex<Value> operator()(int i, int j) const {
    assert( i >= 0 && i <= 3 && j>=0 && j<=3 );
    return _spin[i][j];
  }

  /**
   * Set components by index
   */
  complex<Value> & operator () (int i, int j) {
    assert( i >= 0 && i <= 3 && j>=0 && j<=3 );
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
  void setTS4(complex<Value> in ) {_spin[3][3]=in;}
  //@}

  /** @name Arithmetic operators. */
  //@{
  /**
   * dot product with a polarization vector
   */
  LorentzSpinorBar<Value> dot(const LorentzPolarizationVector & vec) const {
    LorentzSpinorBar<Value> output(_type);
    for(unsigned int ix=0;ix<4;++ix) {
      output[ix]=_spin[3][ix]*vec.t()-_spin[0][ix]*vec.x()
	-_spin[1][ix]*vec.y()-_spin[2][ix]*vec.z();
    }
    return output;
  }

  /**
   * dot product with a 4-momentum
   */
  LorentzSpinorBar<Value> dot(const LorentzMomentum & invec) const {
    LorentzSpinorBar<Value> output(_type);
    LorentzVector<double> vec = UnitRemoval::InvE * invec;
    unsigned int ix;
    for(ix=0;ix<4;++ix) {
      output[ix]=_spin[3][ix]*vec.t()-_spin[0][ix]*vec.x()
	-_spin[1][ix]*vec.y()-_spin[2][ix]*vec.z();
    }
    return output;
  }
  //@}
  
  /** @name Transformations. */
  //@{
  /**
   * return the barred spinor
   */
  LorentzRSSpinor<Value> bar() const;

  /**
   * Standard Lorentz boost specifying the components of the beta vector.
   */
  LorentzRSSpinorBar & boost(double,double,double);

  /**
   * Standard Lorentz boost specifying the beta vector.
   */
  LorentzRSSpinorBar & boost(const Boost &);

  /**
   * General transform
   */
  LorentzRSSpinorBar & transform(const LorentzRotation &);
  //@}

  /** @name Functions related to type. */
  //@{
  /**
   * Return the type of the spinor.
   */
  SpinorType Type() const {return _type;}
  //@}

  /**
   * Current \f$\bar{f}^\alpha(c_LP_L+c_RP_R)f\f$ for general couplings.
   * @param f The unbarred spinor
   * @param left The left-handed coupling, \f$c_L\f$.
   * @param right The right-handed coupling, \f$c_R\f$.
   */
  template <typename ValueB>
  LorentzVector<complex<
    typename BinaryOpTraits<Value,ValueB>::MulT> >
    generalCurrent(LorentzSpinor<ValueB>& f, Complex left, Complex right) {
      typedef complex<typename BinaryOpTraits<Value,ValueB>::MulT> ResultT;
      ResultT output[4];
      unsigned int iz;
      for(iz=0;iz<4;++iz){
	output[iz]=  left*(_spin[iz][0]*f.s1()+_spin[iz][1]*f.s2())
	  +right*(_spin[iz][2]*f.s3()+_spin[iz][3]*f.s4());
      }
      return LorentzVector<ResultT>(output[0],output[1],
				    output[2],output[3]);
    }
  
private:
  /**
   * Type of spinor.
   */
  SpinorType _type;

  /**
   * Storage of the components.
   */
  std::array<std::array<complex<Value>,4>,4> _spin;
};
}
}
#ifndef ThePEG_TEMPLATES_IN_CC_FILE
#include "LorentzRSSpinorBar.tcc"
#endif 

#endif
