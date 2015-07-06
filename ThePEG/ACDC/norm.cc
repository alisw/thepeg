#include <cmath>
#include <iostream>
#include <cstdlib>

int main(int argc, char ** argv) {

  long N = 10000;
  if ( argc > 1 ) N = std::atol(argv[1]);
  double x0 = 5.0;
  if ( argc > 2 ) x0 = std::atof(argv[2]);

  for ( long i = 1; i < N; ++i ) {
    double x = 2.0*x0*double(i - 0.5)/double(N) - x0;
    std::cout << x << '\t' << std::exp(-x*x/2.0)/sqrt(2.0*M_PI) << std::endl;
  }
}

