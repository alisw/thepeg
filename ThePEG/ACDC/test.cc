#include "ACDCGen.h"
#include "DRand48Traits.h"
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "CLHEP/Random/JamesRandom.h"
#include "fpudebug.h"

using namespace ACDCGenerator;

typedef HepJamesRandom RAND;
// typedef DRAND48 RAND;

namespace Const {
const double eps = 0.00001;
}

struct TestFunc {
  TestFunc() : count(0) {}
  static double sqr(double x) { return x*x;}
  double operator()(const DVector & x);
  static double primGauss(double y);
  static double binGauss(double yl, double yu);

  mutable int count;
  static double a;
};

struct HistBin {

  mutable long sumN;
  double vol;
  DVector lo, up;
  static const long NFac = 100;

  HistBin() : sumN(0), vol(1.0) {}
  HistBin(const HistBin & hb);
  HistBin(const ACDCGenCellInfo & cell);

  void setup(const DVector & x, double dx);

  inline void operator() (const DVector & x) const;

  template <typename ACDCType>
  double getChi(ACDCType &, long N, double integral);

  template <typename ACDCType>
  static vector<HistBin> getBins(ACDCType & gen, int dim, int N, double delta);
  template <typename ACDCType>
  static vector<HistBin> getBins(ACDCType & gen);

};

int main(int argc, char ** argv) {
  fpudebug();
  srand48(940801029);
  long N = 10000;
  if ( argc > 1 ) N = std::atol(argv[1]);
  long n = 10;
  if ( argc > 2 ) n = std::atol(argv[2]);
  int ndim = 6;
  if ( argc > 3 ) ndim = std::atol(argv[3]);

  double sumi = 0.0;
  double sumi2 = 0.0;
  std::ofstream outi("integral.dat");

  vector<HistBin> bins;

  std::ofstream bout("bins.dat");
  std::ofstream boutn("binsn.dat");

  for ( int i = 0; i < n; ++i ) {

    RAND r;
    TestFunc f;
    ACDCGen<RAND,TestFunc*> gen(r);
    gen.clear();
    gen.nTry(1000);
    gen.margin(1.1);
    gen.addFunction(std::abs(ndim), &f);

    int icount = f.count;

    std::cout << std::setw(5) << gen.nBins();
    std::cout << std::setw(5) << gen.depth();
    for ( int I =0 ; I < N; ++I ) {
      gen.generate();
      DVector x = gen.lastPoint();
      if ( ndim == -1 ) x[0] = std::pow(Const::eps, r.flat());
      for ( int ib = 0, Nb = bins.size(); ib < Nb; ++ib ) bins[ib](x);
    }
    outi << std::setw(10) << gen.integral() << std::endl;
    std::cout << std::setw(10) << gen.integral();
    std::cout << std::setw(5) << gen.nBins();
    std::cout << std::setw(5) << gen.depth();
    std::cout << std::setw(12) << double(N)/double(f.count);
    std::cout << std::setw(12) << double(N)/double(f.count - icount);
    std::cout << std::setw(10) << gen.compensating();
    std::cout << std::setw(10) << gen.maxsize << std::endl;

    sumi += gen.integral();
    sumi2 += gen.integral()*gen.integral();

    for ( int ib = 0, Nb = bins.size(); ib < Nb; ++ib ) {
      double chi = bins[ib].getChi(gen, N, gen.integral());
      bout << std::setw(10) << chi << std::endl;
      boutn << std::setw(10) << chi
	    << std::setw(10) << bins[ib].sumN
	    << std::setw(15) << bins[ib].vol
	    << std::setw(15) << bins[ib].lo[0]
	    << std::setw(15) << bins[ib].up[0]
	    << std::endl;
    }

    bins = HistBin::getBins(gen);

  }

  sumi /= double(n);
  sumi2 = sqrt(max(sumi2/double(n) - sumi*sumi, 0.0));
  std::cout << std::setw(12) << sumi
	    << std::setw(12) << sumi2
	    << std::setw(12) << sumi2/sqrt(double(n)) << std::endl;

}


double TestFunc::operator()(const DVector & x) {
  ++count;
  if ( x.size() == 1 ) {
    if ( x[0] < Const::eps ) return 0.0;
    return -1.0/(std::log(Const::eps)*x[0]);
  }
  if ( x.size() < 3 ) return 0.0;
  double r2 = 0.0;
  for ( int i = 0; i < 3; ++i ) r2 += sqr(x[i]);
  double r12 = 0.0;
  for ( int i = x.size() - 3; i < int(x.size()); ++i )  r12 += sqr(1.0-x[i]);
  double rg = sqr(x[x.size()/2]-0.5) + sqr(x[x.size()/2-1]-0.5);
  double rg2 = sqr(x[x.size()/2] - x[x.size()/2 + 1]);
  double f1 = r2 > 0.00000001 && r2 < 1.0? 2.0/(M_PI*r2): 0.0;
  double f2 = r12 > 0.00000001 && r12 < 1.0? 2.0/(M_PI*r12): 0.0;
  double f3 = rg > 700.0/a? 0.0:
    a*std::exp(-a*rg)/(M_PI*sqr(erf(std::sqrt(a)/2)));
  double f4 = rg2 > 700.0/a? 0.0:
    a*std::exp(-a*rg2)/(expm1(-a) + std::sqrt(M_PI*a)*erf(std::sqrt(a)));
  double f5 = 1000.0;
  double d5 = 1.0/std::pow(f5,1.0/double(x.size()))/2.0;
  double xi = d5;
  int D = x.size() + 1;
  double dx = (1.0 - 2.0*d5)/double(D);
  for ( int i = 0; i < D-1; ++i ) {
    xi += dx;
    if ( std::abs(x[i] - xi) > d5 ) {
      f5 = 0.0;
      break;
    }
  }
  return (f1 + f2 + f3 + f4 + f5)/5.0;
}

double TestFunc::a = 100.0;

void HistBin::setup(const DVector & x, double dx) {
  sumN = 0;
  lo = DVector(x.size());
  up = DVector(x.size());
  vol = 1.0;
  for ( DVector::size_type i = 0, imax = x.size(); i < imax; ++i ) {
    lo[i] = max(x[i] - dx, 0.0);
    up[i] = min(x[i] + dx, 1.0);
    vol *= (up[i] - lo[i]);
  }
}

HistBin::HistBin(const HistBin & hb)
  : sumN(hb.sumN), vol(hb.vol), lo(hb.lo), up(hb.up) {}

inline void HistBin::operator() (const DVector & x) const {
  for ( DVector::size_type i = 0, imax = x.size(); i < imax; ++i )
    if ( lo[i] > x[i] || x[i] >= up[i] ) return;
  ++sumN;
}

template <typename ACDCType>
double HistBin::getChi(ACDCType & gen, long N, double integral) {
  double sumF = 0.0;
  double sumF2 = 0.0;
  long Nf = max(sumN, 10L)*NFac;
  if ( up.size() == 1 ) Nf = 0;

  double f0;
  double df0;

  if ( up.size() == 1 ) {
    f0 = std::log(std::max(Const::eps, lo[0])/up[0])/std::log(Const::eps);
    df0 = 0.0;
  } else {
    DVector x(up.size());
    for ( long i = 0; i < Nf; ++i ) {
      gen.rnd(lo, up, x);
      double f = (*gen.lastFunction())(x);
      sumF += f;
      sumF2 += f*f;
    }
    f0 = vol*sumF/double(Nf);
    df0 = vol*sqrt(sumF2 - sumF*sumF/double(Nf))/double(Nf);
  }

  double f = double(sumN)*integral/double(N);
  double df = sqrt(double(sumN + 1))*integral/double(N);
  //  double df = sqrt(double(max(sumN,1L)))*integral/double(N);
  return (f - f0)/sqrt(df*df + df0*df0);
}

  
template <typename ACDCType>
vector<HistBin> HistBin::
getBins(ACDCType & gen, int dim, int N, double delta) {
  vector<HistBin> ret(N);
  for ( int i = 0; i < N; ++i ) {
    DVector x(dim);
    double dx = 0.0;
    while ( dx <= 0.0 ) {
      gen.rnd(dim, x);
      double f = (*gen.lastFunction())(x);
      if ( f <= 0.0 ) continue;
      double dv = delta/f;
      if ( dv > 0.1 ) continue;
      dx = std::pow(dv, 1.0/dim);
    }
    ret[i].setup(x, dx);
  }
  return ret;
}

template <typename ACDCType>
vector<HistBin> HistBin::
getBins(ACDCType & gen) {
  const vector<ACDCGenCellInfo> & cells = gen.extractCellInfo();
  vector<HistBin> ret;
  for ( int i = 0, N = cells.size(); i < N; ++i ) {
    if ( cells[i].iup != 0 || cells[i].ilo != 0 ) continue;
    ret.push_back(HistBin(cells[i]));
  }
  return ret;
}

HistBin::HistBin(const ACDCGenCellInfo & cell)
  : sumN(0), vol(1.0), lo(cell.lo), up(cell.up) {
  for ( DVector::size_type i = 0, imax = up.size(); i < imax; ++i )
    vol *= (up[i] - lo[i]);
}
