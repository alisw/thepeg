


#ifdef ThePEG_HAS_UNITS_CHECKING
#include "PhysicalQty.h"
#include "PhysicalQtyOps.h"
#include "PhysicalQtyComplex.h"

namespace ThePEG {

#else

#include <cmath>
namespace ThePEG {

const double ZERO = 0.0;

/// Fractional powers of a double.
template<int P, int R>
double pow(double q) {
  return std::pow(q,double(P)/double(R));
}

#endif

/// Helper class to define unitful quantities.
template <int L, int E, int Q, int DL = 1, int DE = 1, int DQ = 1>
struct QTY {
#ifdef ThePEG_HAS_UNITS_CHECKING
  /// The QTY type is dimensioned.
  typedef Qty<L,E,Q,DL,DE,DQ> Type;
#else
  /// The QTY type is double.
  typedef double	      Type;
#endif
};

}

#endif
