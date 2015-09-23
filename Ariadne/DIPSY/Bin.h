// -*- C++ -*-
#ifndef ThePEG_Bin_H
#define ThePEG_Bin_H
//
// This is the declaration of the Bin class.
//

#include "ThePEG/Config/ThePEG.h"

namespace ThePEG {

/**
 * The Bin class is used to valculate the average of several measurements.
 */
template <typename T, typename T2 = T>
class Bin {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  inline Bin(): sum(T()), sum2(T2()), sumw(0.0), n(0) {}

  /**
   * The constructor fo an initial value.
   */
  inline Bin(const T & t): sum(t), sum2(t*t), sumw(0.0), n(0) {}

  /**
   * The copy constructor fo an initial value.
   */
  inline Bin(const Bin<T,T2> & b): sum(b.sum), sum2(b.sum2),
				     sumw(b.sumw), n(b.n) {}

  /**
   * The assignment operator.
   */
  Bin<T,T2> & operator=(const Bin<T,T2> & b) {
    sum = b.sum;
    sum2 = b.sum2;
    sumw = b.sumw;
    n = b.n;
    return *this;
  }
  //@}

  /**
   * Add a value
   */
  Bin<T,T2> & operator+=(const T & t) {
    return fill(t, 1.0);
  }

  Bin<T,T2> & fill(const T & t, double w) {
    if ( n >= 0 ) {
      sum += t*w;
      sum2 += t*t*w;
      sumw += w;
      ++n;
    }
    return *this;
  }

  T var() const {
    return n > 0? sum2/sumw - average()*average(): sum2;
  }

  T2 err2() const {
    return n > 0? var()/double(abs(n)): sum2;
  }

  T err() const {
    return sqrt(err2());
  }

  T average() const {
    return n > 0? sum/sumw: sum;
  }

  T operator()() const {
    return average();
  }

  Bin<T,T2> & operator+=(const Bin<T,T2> & b) {
    T av = average() + b.average();
    T2 e2 = err2() + b.err2();
    sum = av;
    sum2 = e2;
    sumw = 1.0;
    n = -1;
    return *this;
  }

  Bin<T,T2> & operator-=(const Bin<T,T2> & b) {
    T av = average() - b.average();
    T2 e2 = err2() + b.err2();
    sum = av;
    sum2 = e2;
    sumw = 1.0;
    n = -1;
    return *this;
  }

  Bin<T,T2> & operator*=(const Bin<T,T2> & b) {
    T av = average()*b.average();
    T2 e2 = av*av*(err2()/(average()*average) +
		   b.err2()/(b.average()*b.average()));
    sum = av;
    sum2 = e2;
    sumw = 1.0;
    n = -1;
    return *this;
  }

  Bin<T,T2>  sqr() const {
    Bin<T,T2> ret = *this;
    ret.sum = average()*average();
    ret.sum2 = 4.0*err2()*ret.sum;
    ret.sumw = 1.0;
    ret.n = -1;
    return ret;
  }


  Bin<T,T2> & operator/=(const Bin<T,T2> & b) {
    T av = average()/b.average();
    T2 e2 = av*av*(err2()/(average()*average) +
		   b.err2()/(b.average()*b.average()));
    sum = av;
    sum2 = e2;
    sumw = 1.0;
    n = -1;
    return *this;
  }

  Bin<T,T2> operator+(const Bin<T,T2> & b) const {
    Bin<T,T2> ret = *this;
    return ret += b;
  }

  Bin<T,T2> operator-(const Bin<T,T2> & b) const {
    Bin<T,T2> ret = *this;
    return ret -= b;
  }

  Bin<T,T2> operator*(const Bin<T,T2> & b) const {
    Bin<T,T2> ret = *this;
    return ret *= b;
  }

  Bin<T,T2> operator/(const Bin<T,T2> & b) const {
    Bin<T,T2> ret = *this;
    return ret /= b;
  }

private:

  T sum;
  T2 sum2;
  double sumw;
  long n;

};

}

#endif /* ThePEG_Bin_H */
