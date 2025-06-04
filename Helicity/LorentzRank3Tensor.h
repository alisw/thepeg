// -*- C++ -*-
//
// LorentzRank3Tensor.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2019 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_LorentzRank3Tensor_H
#define ThePEG_LorentzRank3Tensor_H
// This is the declaration of the LorentzRank3Tensor class.

#include "ThePEG/Config/PhysicalQtyComplex.h"
#include "ThePEG/Config/ThePEG.h"
#include "LorentzTensor.h"

namespace ThePEG {
namespace Helicity {

// compiler magic needs these pre-declarations to make friend templates work
template<typename Value> class LorentzRank3Tensor;

/**
 *  The LorentzRank3Tensor class is designed to implement the storage of a
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
class LorentzRank3Tensor {

public:

  /**
   * Default zero constructor.
   */
  LorentzRank3Tensor() = default;

  /**
   * Get components by indices.
   */
  complex<Value> operator () (int i, int j, int k) const {
    assert( i>=0 && i<=3 && j>=0 && j<=3 && k>=0 && k<=3);
    return _tensor[i][j][k];
  }

  /**
   * Set components by indices.
   */
  complex<Value> & operator () (int i, int j, int k) {
    assert( i>=0 && i<=3 && j>=0 && j<=3 && k>=0 && k<=3);
    return _tensor[i][j][k];
  }
  //@}

  /** @name Transformations. */
  //@{
  /**
   * Standard Lorentz boost specifying the components of the beta vector.
   */
  LorentzRank3Tensor & boost(double,double,double);

  /**
   * Standard Lorentz boost specifying the beta vector.
   */
  LorentzRank3Tensor<Value> & boost(const Boost & b) {
    return boost(b.x(), b.y(), b.z());
  }

  /**
   * General Lorentz transformation
   */
  LorentzRank3Tensor & transform(const SpinOneLorentzRotation & r){
    unsigned int ix,iy,iz,ixa,iya,iza;
    LorentzRank3Tensor<Value> output;
    complex<Value> temp;
    for(ix=0;ix<4;++ix) {
      for(iy=0;iy<4;++iy) {
	for(iz=0;iz<4;++iz) {
	  output(ix,iy,iz) = complex<Value>();
	  for(ixa=0;ixa<4;++ixa) {
	    for(iya=0;iya<4;++iya) {
	      for(iza=0;iza<4;++iza)
		 output(ix,iy,iz) += r(ix,ixa)*r(iy,iya)*r(iz,iza)*(*this)(ixa,iya,iza);
	    }
	  }
	}
      }
    }
    *this=output;
    return *this;
  }
  
  /**
   * Return the complex conjugate.
   */
  LorentzRank3Tensor<Value> conjugate() {
    LorentzRank3Tensor<Value> output;
    for(unsigned int ix=0;ix<4;++ix) {
      for(unsigned int iy=0;iy<4;++iy) {
	for(unsigned int iz=0;iz<4;++iz) {
	  output(ix,iy,iz) = conj(_tensor[ix][iy][iz]);
	}
      }
    }
    return output;
  }

  //@}

  /** @name Arithmetic operators. */
  //@{
  /**
   * Scaling with a complex number
   */
  LorentzRank3Tensor<Value> operator*=(Complex a) {
    for(int ix=0;ix<4;++ix)
      for(int iy=0;iy<4;++iy)
	for(int iz=0;iz<4;++iz) _tensor[ix][iy][iz]*=a;
    return *this;
  }

  /**
   * Scalar product with other tensor
   */
  template <typename T, typename U>
  friend auto
  operator*(const LorentzRank3Tensor<T> & t, const LorentzRank3Tensor<U> & u) -> decltype(t.xx()*u.xx());
    
  /**
   * Addition.
   */
  LorentzRank3Tensor<Value> operator+(const LorentzRank3Tensor<Value> & in) const {
    LorentzRank3Tensor<Value> output;
    for(int ix=0;ix<4;++ix)
      for(int iy=0;iy<4;++iy)
	for(int iz=0;iz<4;++iz) output(ix,iy,iz) = _tensor[ix][iy][iz] + in(ix,iy,iz);
  }
  
  /**
   * Subtraction.
   */
  LorentzRank3Tensor<Value> operator-(const LorentzRank3Tensor<Value> & in) const {
    LorentzRank3Tensor<Value> output;
    for(int ix=0;ix<4;++ix)
      for(int iy=0;iy<4;++iy)
	for(int iz=0;iz<4;++iz) output(ix,iy,iz) = _tensor[ix][iy][iz] - in(ix,iy,iz);
  }
  
  /**
   * Dot product with the ith index
   */
  template<typename ValueB>
  auto dot(const LorentzVector<complex<ValueB> > & vec, unsigned int iloc) const 
    -> LorentzTensor<decltype(ValueB()*Value())> {
    LorentzTensor<decltype(ValueB()*Value())> output;
    if(iloc==0) {
      for(unsigned int iy=0;iy<4;++iy) {
	for(unsigned int iz=0;iz<4;++iz) {
	  output(iy,iz) =
	    vec.t()*_tensor[3][iy][iz] - vec.x()*_tensor[0][iy][iz] -
	    vec.y()*_tensor[1][iy][iz] - vec.z()*_tensor[2][iy][iz];
	}
      }
    }
    else if(iloc==1) {
      for(unsigned int iy=0;iy<4;++iy) {
	for(unsigned int iz=0;iz<4;++iz) {
	  output(iy,iz) =
	    vec.t()*_tensor[iy][3][iz] - vec.x()*_tensor[iy][0][iz] -
	    vec.y()*_tensor[iy][1][iz] - vec.z()*_tensor[iy][2][iz];
	}
      }
    }
    else if(iloc==2) {
      for(unsigned int iy=0;iy<4;++iy) {
	for(unsigned int iz=0;iz<4;++iz) {
	  output(iy,iz) =
	    vec.t()*_tensor[iy][iz][3] - vec.x()*_tensor[iy][iz][0] -
	    vec.y()*_tensor[iy][iz][1] - vec.z()*_tensor[iy][iz][2];
	}
      }
    }
    else
      assert(false);
    return output;
  }
  
  /**
   *  dot product with momentum
   */ 
  auto dot (const Lorentz5Momentum & vec,unsigned int iloc) const 
  -> LorentzTensor<decltype(vec.x()*Value())>
  {
    LorentzTensor<decltype(vec.x()*Value())> output;
    if(iloc==0) {
      for(unsigned int iy=0;iy<4;++iy) {
	for(unsigned int iz=0;iz<4;++iz) {
	  output(iy,iz) =
	    vec.t()*_tensor[3][iy][iz] - vec.x()*_tensor[0][iy][iz] -
	    vec.y()*_tensor[1][iy][iz] - vec.z()*_tensor[2][iy][iz];
	}
      }
    }
    else if(iloc==1) {
      for(unsigned int iy=0;iy<4;++iy) {
	for(unsigned int iz=0;iz<4;++iz) {
	  output(iy,iz) =
	    vec.t()*_tensor[iy][3][iz] - vec.x()*_tensor[iy][0][iz] -
	    vec.y()*_tensor[iy][1][iz] - vec.z()*_tensor[iy][2][iz];
	}
      }
    }
    else if(iloc==2) {
      for(unsigned int iy=0;iy<4;++iy) {
	for(unsigned int iz=0;iz<4;++iz) {
	  output(iy,iz) =
	    vec.t()*_tensor[iy][iz][3] - vec.x()*_tensor[iy][iz][0] -
	    vec.y()*_tensor[iy][iz][1] - vec.z()*_tensor[iy][iz][2];
	}
      }
    }
    else
      assert(false);
    return output;
  } 
  //@}

private:

  /**
   * The components.
   */
  std::array<std::array<std::array<complex<Value>,4>,4>,4> _tensor;

};

/**
 * Multiplication by a complex number.
 */
template<typename T, typename U> 
inline auto
operator*(complex<U> a, const LorentzRank3Tensor<T> & t) -> LorentzRank3Tensor<decltype(a.real()*t.xx().real())> {
  LorentzRank3Tensor<decltype(a.real()*t.xx().real())> output;
  for(int ix=0;ix<4;++ix)
    for(int iy=0;iy<4;++iy)
      for(int iz=0;iz<4;++iz) output(ix,iy,iz) = a*t(ix,iy,iz);
  return output;
}

/**
 * Multiplication by a complex number.
 */
template<typename T, typename U> 
inline auto
operator*(const LorentzRank3Tensor<T> & t,complex<U> a) -> LorentzRank3Tensor<decltype(a.real()*t.xx().real())> {
  LorentzRank3Tensor<decltype(a.real()*t.xx().real())> output;
  for(int ix=0;ix<4;++ix)
    for(int iy=0;iy<4;++iy)
      for(int iz=0;iz<4;++iz) output(ix,iy,iz) = a*t(ix,iy,iz);
  return output;
}

}
}

#ifndef ThePEG_TEMPLATES_IN_CC_FILE
#include "LorentzRank3Tensor.tcc"
#endif 

#endif
