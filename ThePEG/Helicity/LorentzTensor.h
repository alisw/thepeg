// -*- C++ -*-
//
// LorentzTensor.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2011 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_LorentzTensor_H
#define ThePEG_LorentzTensor_H
// This is the declaration of the LorentzTensor class.


#include "ThePEG/Config/ThePEG.h"
#include "LorentzPolarizationVector.h"

namespace ThePEG {
namespace Helicity {

// compiler magic needs these pre-declarations to make friend templates work
template<typename Value> class LorentzTensor;

/**
 *  The LorentzTensor class is designed to implement the storage of a
 *  complex tensor to be used to representation the wavefunction of a
 *  spin-2 particle.
 *
 *  At the moment it only implements the storage of the tensor
 *  components but it is envisaged that it will be extended to include
 *  boost methods etc.
 *
 * @author Peter Richardson
 *
 */

template<typename Value> 
class LorentzTensor {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default zero constructor.
   */
  LorentzTensor()  {
    for(unsigned int ix=0;ix<4;++ix)
      for(unsigned int iy=0;iy<4;++iy)
	_tensor[ix][iy]=Value();
  }

  /**
   * Constructor specifyign all components.
   */
  LorentzTensor(complex<Value> xx, complex<Value> xy, 
		complex<Value> xz, complex<Value> xt,
		complex<Value> yx, complex<Value> yy,
		complex<Value> yz, complex<Value> yt,
		complex<Value> zx, complex<Value> zy,
		complex<Value> zz, complex<Value> zt,
		complex<Value> tx, complex<Value> ty,
		complex<Value> tz, complex<Value> tt){
    _tensor[0][0]=xx;_tensor[0][1]=xy;_tensor[0][2]=xz;_tensor[0][3]=xt;
    _tensor[1][0]=yx;_tensor[1][1]=yy;_tensor[1][2]=yz;_tensor[1][3]=yt;
    _tensor[2][0]=zx;_tensor[2][1]=zy;_tensor[2][2]=zz;_tensor[2][3]=zt;
    _tensor[3][0]=tx;_tensor[3][1]=ty;_tensor[3][2]=tz;_tensor[3][3]=tt;
  }

  /**
   * Constructor in terms of two polarization vectors.
   */
  LorentzTensor(const LorentzPolarizationVector & p,
		const LorentzPolarizationVector & q) {
    setXX(p.x() * q.x()); setYX(p.y() * q.x());
    setZX(p.z() * q.x()); setTX(p.t() * q.x());
    setXY(p.x() * q.y()); setYY(p.y() * q.y());
    setZY(p.z() * q.y()); setTY(p.t() * q.y());
    setXZ(p.x() * q.z()); setYZ(p.y() * q.z());
    setZZ(p.z() * q.z()); setTZ(p.t() * q.z());
    setXT(p.x() * q.t()); setYT(p.y() * q.t());
    setZT(p.z() * q.t()); setTT(p.t() * q.t());
  }
  //@}

  /** @name Access individual components. */
  //@{
  /**
   * Get x,x component.
   */
  complex<Value> xx() const {return _tensor[0][0];}

  /**
   * Get y,x component.
   */
  complex<Value> yx() const {return _tensor[1][0];}
  /**
   * Get z,x component.
   */
  complex<Value> zx() const {return _tensor[2][0];}

  /**
   * Get t,x component.
   */
  complex<Value> tx() const {return _tensor[3][0];}

  /**
   * Get x,y component.
   */
  complex<Value> xy() const {return _tensor[0][1];}

  /**
   * Get y,y component.
   */
  complex<Value> yy() const {return _tensor[1][1];}

  /**
   * Get z,y component.
   */
  complex<Value> zy() const {return _tensor[2][1];}

  /**
   * Get t,y component.
   */
  complex<Value> ty() const {return _tensor[3][1];}

  /**
   * Get x,z component.
   */
  complex<Value> xz() const {return _tensor[0][2];}

  /**
   * Get y,z component.
   */
  complex<Value> yz() const {return _tensor[1][2];}

  /**
   * Get z,z component.
   */
  complex<Value> zz() const {return _tensor[2][2];}

  /**
   * Get t,z component.
   */
  complex<Value> tz() const {return _tensor[3][2];}

  /**
   * Get x,t component.
   */
  complex<Value> xt() const {return _tensor[0][3];}

  /**
   * Get y,t component.
   */
  complex<Value> yt() const {return _tensor[1][3];}

  /**
   * Get z,t component.
   */
  complex<Value> zt() const {return _tensor[2][3];}

  /**
   * Get t,t component.
   */
  complex<Value> tt() const {return _tensor[3][3];}

  /**
   * Set x,x component.
   */
  void setXX(complex<Value> a) {_tensor[0][0]=a;}

  /**
   * Set y,x component.
   */
  void setYX(complex<Value> a) {_tensor[1][0]=a;}

  /**
   * Set z,x component.
   */
  void setZX(complex<Value> a) {_tensor[2][0]=a;}

  /**
   * Set t,x component.
   */
  void setTX(complex<Value> a) {_tensor[3][0]=a;}

  /**
   * Set x,y component.
   */
  void setXY(complex<Value> a) {_tensor[0][1]=a;}

  /**
   * Set y,y component.
   */
  void setYY(complex<Value> a) {_tensor[1][1]=a;}

  /**
   * Set z,y component.
   */
  void setZY(complex<Value> a) {_tensor[2][1]=a;}

  /**
   * Set t,y component.
   */
  void setTY(complex<Value> a) {_tensor[3][1]=a;}

  /**
   * Set x,z component.
   */
  void setXZ(complex<Value> a) {_tensor[0][2]=a;}

  /**
   * Set y,z component.
   */
  void setYZ(complex<Value> a) {_tensor[1][2]=a;}

  /**
   * Set z,z component.
   */
  void setZZ(complex<Value> a) {_tensor[2][2]=a;}

  /**
   * Set t,z component.
   */
  void setTZ(complex<Value> a) {_tensor[3][2]=a;}

  /**
   * Set x,t component.
   */
  void setXT(complex<Value> a) {_tensor[0][3]=a;}

  /**
   * Set y,t component.
   */
  void setYT(complex<Value> a) {_tensor[1][3]=a;}

  /**
   * Set z,t component.
   */
  void setZT(complex<Value> a) {_tensor[2][3]=a;}

  /**
   * Set t,t component.
   */
  void setTT(complex<Value> a) {_tensor[3][3]=a;}

  /**
   * Get components by indices.
   */
  complex<Value> operator () (int i, int j) const {
    assert( i>=0 && i<=3 && j>=0 && j<=3);
    return _tensor[i][j];
  }

  /**
   * Set components by indices.
   */
  complex<Value> & operator () (int i, int j) {
    assert( i>=0 && i<=3 && j>=0 && j<=3);
    return _tensor[i][j];
  }
  //@}

  /** @name Transformations. */
  //@{
  /**
   * Standard Lorentz boost specifying the components of the beta vector.
   */
  LorentzTensor & boost(double,double,double);

  /**
   * Standard Lorentz boost specifying the beta vector.
   */
  LorentzTensor<Value> & boost(const Boost & b) {
    return boost(b.x(), b.y(), b.z());
  }

  /**
   * General Lorentz transformation
   */
  LorentzTensor & transform(const SpinOneLorentzRotation & r){
    unsigned int ix,iy,ixa,iya;
    LorentzTensor<Value> output;
    complex<Value> temp;
    for(ix=0;ix<4;++ix) {
      for(iy=0;iy<4;++iy) {
	temp=complex<Value>();
	for(ixa=0;ixa<4;++ixa) {
	  for(iya=0;iya<4;++iya)
	    temp+=r(ix,ixa)*r(iy,iya)*(*this)(ixa,iya);
	}
	output(ix,iy)=temp;
      }
    }
    *this=output;
    return *this;
  }
  
  /**
   * Return the complex conjugate.
   */
  LorentzTensor<Value> conjugate() {
    return LorentzTensor<Value>(conj(xx()), conj(xy()), conj(xz()), conj(xt()),
				conj(yx()), conj(yy()), conj(yz()), conj(yt()),
				conj(zx()), conj(zy()), conj(zz()), conj(zt()),
				conj(tx()), conj(ty()), conj(tz()), conj(tt()));
  }

  //@}

  /** @name Arithmetic operators. */
  //@{
  /**
   * Scaling with a complex number
   */
  LorentzTensor<Value> operator*=(Complex a) {
    for(int ix=0;ix<4;++ix)
      for(int iy=0;iy<4;++iy) _tensor[ix][iy]*=a;
    return *this;
  }

  /**
   * Scalar product with other tensor
   */
  template <typename T, typename U>
  friend complex<typename BinaryOpTraits<T,U>::MulT> 
  operator*(const LorentzTensor<T> & t, const LorentzTensor<U> & u);
    
  /**
   * Addition.
   */
  LorentzTensor<Value> operator+(const LorentzTensor<Value> & in) const {
    return LorentzTensor<Value>(xx()+in.xx(),xy()+in.xy(),xz()+in.xz(),xt()+in.xt(),
				yx()+in.yx(),yy()+in.yy(),yz()+in.yz(),yt()+in.yt(),
				zx()+in.zx(),zy()+in.zy(),zz()+in.zz(),zt()+in.zt(),
				tx()+in.tx(),ty()+in.ty(),tz()+in.tz(),tt()+in.tt());
  }
  
  /**
   * Subtraction.
   */
  LorentzTensor<Value> operator-(const LorentzTensor<Value> & in) const {
    return LorentzTensor<Value>(xx()-in.xx(),xy()-in.xy(),xz()-in.xz(),xt()-in.xt(),
				yx()-in.yx(),yy()-in.yy(),yz()-in.yz(),yt()-in.yt(),
				zx()-in.zx(),zy()-in.zy(),zz()-in.zz(),zt()-in.zt(),
				tx()-in.tx(),ty()-in.ty(),tz()-in.tz(),tt()-in.tt());
  }

  /**
   * Trace
   */
  complex<Value> trace() const {
    return _tensor[3][3]-_tensor[0][0]-_tensor[1][1]-_tensor[2][2];
  }
  //@}

  /**
   *  Various dot products
   */
  //@{
  /**
   *  First index dot product with polarization vector
   */
  LorentzVector<complex<Value> > preDot (const LorentzPolarizationVector & vec) const {
    LorentzVector<complex<Value> > output;
    output.setX(vec.t()*_tensor[3][0]-vec.x()*_tensor[0][0]-
		vec.y()*_tensor[1][0]-vec.z()*_tensor[2][0]);
    output.setY(vec.t()*_tensor[3][1]-vec.x()*_tensor[0][1]-
		vec.y()*_tensor[1][1]-vec.z()*_tensor[2][1]);
    output.setZ(vec.t()*_tensor[3][2]-vec.x()*_tensor[0][2]-
		vec.y()*_tensor[1][2]-vec.z()*_tensor[2][2]);
    output.setT(vec.t()*_tensor[3][3]-vec.x()*_tensor[0][3]-
		vec.y()*_tensor[1][3]-vec.z()*_tensor[2][3]);
    return output;
  }

  /**
   *  Second index dot product with polarization vector
   */
  LorentzVector<complex<Value> > postDot(const LorentzPolarizationVector & vec) const {
    LorentzVector<complex<Value> > output;
    output.setX(vec.t()*_tensor[0][3]-vec.x()*_tensor[0][0]-
		vec.y()*_tensor[0][1]-vec.z()*_tensor[0][2]);
    output.setY(vec.t()*_tensor[1][3]-vec.x()*_tensor[1][0]-
		vec.y()*_tensor[1][1]-vec.z()*_tensor[1][2]);
    output.setZ(vec.t()*_tensor[2][3]-vec.x()*_tensor[2][0]-
		vec.y()*_tensor[2][1]-vec.z()*_tensor[2][2]);
    output.setT(vec.t()*_tensor[3][3]-vec.x()*_tensor[3][0]-
		vec.y()*_tensor[3][1]-vec.z()*_tensor[3][2]);
    return output;
  }

  /**
   *  First index dot product with momentum
   */
  LorentzVector<complex<typename BinaryOpTraits<Value,Energy>::MulT> > 
  preDot (const Lorentz5Momentum & vec) const {
    LorentzVector<complex<typename BinaryOpTraits<Value,Energy>::MulT> > output;
    output.setX(vec.t()*_tensor[3][0]-vec.x()*_tensor[0][0]-
		vec.y()*_tensor[1][0]-vec.z()*_tensor[2][0]);
    output.setY(vec.t()*_tensor[3][1]-vec.x()*_tensor[0][1]-
		vec.y()*_tensor[1][1]-vec.z()*_tensor[2][1]);
    output.setZ(vec.t()*_tensor[3][2]-vec.x()*_tensor[0][2]-
		vec.y()*_tensor[1][2]-vec.z()*_tensor[2][2]);
    output.setT(vec.t()*_tensor[3][3]-vec.x()*_tensor[0][3]-
		vec.y()*_tensor[1][3]-vec.z()*_tensor[2][3]);
    return output;
  }

  /**
   *  Second index dot product with momentum
   */
  LorentzVector<complex<typename BinaryOpTraits<Value,Energy>::MulT> > 
  postDot(const Lorentz5Momentum & vec) const {
    LorentzVector<complex<typename BinaryOpTraits<Value,Energy>::MulT> > output;
    output.setX(vec.t()*_tensor[0][3]-vec.x()*_tensor[0][0]-
		vec.y()*_tensor[0][1]-vec.z()*_tensor[0][2]);
    output.setY(vec.t()*_tensor[1][3]-vec.x()*_tensor[1][0]-
		vec.y()*_tensor[1][1]-vec.z()*_tensor[1][2]);
    output.setZ(vec.t()*_tensor[2][3]-vec.x()*_tensor[2][0]-
		vec.y()*_tensor[2][1]-vec.z()*_tensor[2][2]);
    output.setT(vec.t()*_tensor[3][3]-vec.x()*_tensor[3][0]-
		vec.y()*_tensor[3][1]-vec.z()*_tensor[3][2]);
    return output;
  }
  //@}
private:

  /**
   * The components.
   */
  complex<Value> _tensor[4][4];

};

/**
 * Multiplication by a complex number.
 */
template<typename T, typename U> 
inline LorentzTensor<typename BinaryOpTraits<T,U>::MulT> 
operator*(complex<U> a, const LorentzTensor<T> & t) {
  return LorentzTensor<typename BinaryOpTraits<T,U>::MulT>
    (a*t.xx(), a*t.xy(), a*t.xz(), a*t.xt(),
     a*t.yx(), a*t.yy(), a*t.yz(), a*t.yt(),
     a*t.zx(), a*t.zy(), a*t.zz(), a*t.zt(),
     a*t.tx(), a*t.ty(), a*t.tz(), a*t.tt());
}

/**
 * Multiply a LorentzVector by a LorentzTensor.
 */
template<typename T, typename U> 
inline LorentzVector<typename BinaryOpTraits<complex<T>,U>::MulT>
operator*(const LorentzVector<U> & invec, 
	  const LorentzTensor<T> & inten) {
  LorentzVector<typename BinaryOpTraits<complex<T>,U>::MulT> outvec;
  outvec.setX(invec.t()*inten(3,0)-invec.x()*inten(0,0)
	      -invec.y()*inten(1,0)-invec.z()*inten(2,0));
  outvec.setY(invec.t()*inten(3,1)-invec.x()*inten(0,1)
	      -invec.y()*inten(1,1)-invec.z()*inten(2,1));
  outvec.setZ(invec.t()*inten(3,2)-invec.x()*inten(0,2)
	      -invec.y()*inten(1,2)-invec.z()*inten(2,2));
  outvec.setT(invec.t()*inten(3,3)-invec.x()*inten(0,3)
	      -invec.y()*inten(1,3)-invec.z()*inten(2,3));
  return outvec;
}

/**
 * Multiply a LorentzTensor by a LorentzVector.
 */
template<typename T, typename U> 
inline LorentzVector<typename BinaryOpTraits<complex<T>,U>::MulT>
operator*(const LorentzTensor<T> & inten, const LorentzVector<U> & invec){
  LorentzVector<typename BinaryOpTraits<complex<T>,U>::MulT> outvec;
  outvec.setX(invec.t()*inten(0,3)-invec.x()*inten(0,0)
	      -invec.y()*inten(0,1)-invec.z()*inten(0,2));
  outvec.setY(invec.t()*inten(1,3)-invec.x()*inten(1,0)
	      -invec.y()*inten(1,1)-invec.z()*inten(1,2));
  outvec.setZ(invec.t()*inten(2,3)-invec.x()*inten(2,0)
	      -invec.y()*inten(2,1)-invec.z()*inten(2,2));
  outvec.setT(invec.t()*inten(3,3)-invec.x()*inten(3,0)
	      -invec.y()*inten(3,1)-invec.z()*inten(3,2));
  return outvec;
}

/**
 * Multiply a LorentzTensor by a LorentzTensor
 */
template <typename T, typename U>
inline complex<typename BinaryOpTraits<T,U>::MulT>
operator*(const LorentzTensor<T> & t, const LorentzTensor<U> & u) {
  typedef complex<typename BinaryOpTraits<T,U>::MulT> RetT;
  RetT output=RetT(),temp;
  for(unsigned int ix=0;ix<4;++ix) {
    temp = t._tensor[ix][3]*u._tensor[ix][3];
    for(unsigned int iy=0;iy<3;++iy) {
      temp+= t._tensor[ix][iy]*u._tensor[ix][iy];
    }
    if(ix<3) output-=temp;
    else     output+=temp;
  }
  return output;
}

}
}

#ifndef ThePEG_TEMPLATES_IN_CC_FILE
#include "LorentzTensor.tcc"
#endif 

#endif
