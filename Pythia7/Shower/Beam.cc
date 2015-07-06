// Function definitions not found in the header.

#include "Basics.h"
#include "Beam.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/PDT/ParticleData.h"

using namespace Pythia7;
using namespace Shower;

//**************************************************************************

// Base class for parton distribution functions.

double PDF::xfx(long id, double x, double Q2) { 
  
  if (x != xSav || Q2 != Q2Sav) {xfUpdate(x, Q2); xSav = x; Q2Sav = Q2;}

  // Baryon beam.
  if (abs(idBeam) > 100) { 
    long idNow = (idBeam > 0) ? id : -id;
    long idAbs = abs(id);
    if (idNow == 0 || idAbs == 21) return max(0., xg);  
    if (idNow == 1) return max(0., xd);
    if (idNow == -1) return max(0., xdbar);
    if (idNow == 2) return max(0., xu);
    if (idNow == -2) return max(0., xubar);
    if (idAbs == 3) return max(0., xs);
    if (idAbs == 4) return max(0., xc);
    if (idAbs == 5) return max(0., xb);
    //    if (idAbs == 5) return max(0., xb); //LEIF? DO WE WANT TOP?
    return 0.;

  // Lepton beam.
  } else {
    if (id == idBeam ) return max(0., xlepton);
    if (abs(id) == 22) return max(0., xgamma);
    return 0.;
  }
   
}

//*********

bool PDF::isParton(long id) const {

  // Normal parton content for hadron beams.
  if (abs(idBeam) > 100 && id != 0 && (abs(id) <=5 || id == 21) )
    return true;

  // Normal parton content for lepton beams.
  if ( (abs(idBeam) == 11 || abs(idBeam) == 13 || abs(idBeam) == 15) 
    && (id == idBeam || id == 22) ) return true; 

  // No parton content in unknown beams.
  return false;    

} 

//*********

bool PDF::isValence(long id) const { 

  // Valence content for p/pbar beams.
  if (idBeam == 2212 && (id == 1 || id == 2)) return true;
  if (idBeam == -2212 && (id == -1 || id == -2)) return true;

  // Valence content for lepton beams.
  if ((abs(idBeam) == 11 || abs(idBeam) == 13 || abs(idBeam) == 15)
    && id == idBeam) return true;

  // Default.
  return false;
    
} 

//*********

// Gives the GRV 94 L (leading order) parton distribution function set
// in parametrized form. Authors: M. Glueck, E. Reya and A. Vogt.

void GRV94L::xfUpdate(double x, double Q2) {
 
  // Common expressions.
  double mu2  = 0.23;
  double lam2 = 0.2322 * 0.2322;
  double s  = log (log(Q2/lam2) / log(mu2/lam2));
  double ds = sqrt(s);
  double s2 = s * s;
  double s3 = s2 * s;
 
  // uv :
  double nu  =  2.284 + 0.802 * s + 0.055 * s2;
  double aku =  0.590 - 0.024 * s;
  double bku =  0.131 + 0.063 * s;
  double au  = -0.449 - 0.138 * s - 0.076 * s2;
  double bu  =  0.213 + 2.669 * s - 0.728 * s2;
  double cu  =  8.854 - 9.135 * s + 1.979 * s2;
  double du  =  2.997 + 0.753 * s - 0.076 * s2;
  double uv  = grvv (x, nu, aku, bku, au, bu, cu, du);

  // dv :
  double nd  =  0.371 + 0.083 * s + 0.039 * s2;
  double akd =  0.376;
  double bkd =  0.486 + 0.062 * s;
  double ad  = -0.509 + 3.310 * s - 1.248 * s2;
  double bd  =  12.41 - 10.52 * s + 2.267 * s2;
  double cd  =  6.373 - 6.208 * s + 1.418 * s2;
  double dd  =  3.691 + 0.799 * s - 0.071 * s2;
  double dv  = grvv (x, nd, akd, bkd, ad, bd, cd, dd);
  
  // udb :
  double alx =  1.451;
  double bex =  0.271;
  double akx =  0.410 - 0.232 * s;
  double bkx =  0.534 - 0.457 * s;
  double agx =  0.890 - 0.140 * s;
  double bgx = -0.981;
  double cx  =  0.320 + 0.683 * s;
  double dx  =  4.752 + 1.164 * s + 0.286 * s2;
  double ex  =  4.119 + 1.713 * s;
  double esx =  0.682 + 2.978 * s;
  double udb = grvw (x, s, alx, bex, akx, bkx, agx, bgx, cx,
    dx, ex, esx);

  // del :
  double ne  =  0.082 + 0.014 * s + 0.008 * s2;
  double ake =  0.409 - 0.005 * s;
  double bke =  0.799 + 0.071 * s;
  double ae  = -38.07 + 36.13 * s - 0.656 * s2;
  double be  =  90.31 - 74.15 * s + 7.645 * s2;
  double ce  =  0.;
  double de  =  7.486 + 1.217 * s - 0.159 * s2;
  double del = grvv (x, ne, ake, bke, ae, be, ce, de);
 
  // sb :
  double sts =  0.;
  double als =  0.914;
  double bes =  0.577;
  double aks =  1.798 - 0.596 * s;
  double as  = -5.548 + 3.669 * ds - 0.616 * s;
  double bs  =  18.92 - 16.73 * ds + 5.168 * s;
  double dst =  6.379 - 0.350 * s  + 0.142 * s2;
  double est =  3.981 + 1.638 * s;
  double ess =  6.402;
  double sb  = grvs (x, s, sts, als, bes, aks, as, bs, dst, est, ess);
 
  // cb :
  double stc =  0.888;
  double alc =  1.01;
  double bec =  0.37;
  double akc =  0.;
  double ac  =  0.;
  double bc  =  4.24  - 0.804 * s;
  double dct =  3.46  - 1.076 * s;
  double ect =  4.61  + 1.49  * s;
  double esc =  2.555 + 1.961 * s;
  double chm = grvs (x, s, stc, alc, bec, akc, ac, bc, dct, ect, esc);
 
  // bb :
  double stb =  1.351;
  double alb =  1.00;
  double beb =  0.51;
  double akb =  0.;
  double ab  =  0.;
  double bb  =  1.848;
  double dbt =  2.929 + 1.396 * s;
  double ebt =  4.71  + 1.514 * s;
  double esb =  4.02  + 1.239 * s;
  double bot = grvs (x, s, stb, alb, beb, akb, ab, bb, dbt, ebt, esb);
 
  // gl :
  double alg =  0.524;
  double beg =  1.088;
  double akg =  1.742 - 0.930 * s;
  double bkg =                     - 0.399 * s2;
  double ag  =  7.486 - 2.185 * s;
  double bg  =  16.69 - 22.74 * s  + 5.779 * s2;
  double cg  = -25.59 + 29.71 * s  - 7.296 * s2;
  double dg  =  2.792 + 2.215 * s  + 0.422 * s2 - 0.104 * s3;
  double eg  =  0.807 + 2.005 * s;
  double esg =  3.841 + 0.316 * s;
  double gl  = grvw (x, s, alg, beg, akg, bkg, ag, bg, cg,
    dg, eg, esg);

  // update values
  xg = gl;
  xu = uv + 0.5*(udb - del);
  xd = dv + 0.5*(udb + del); 
  xubar = 0.5*(udb - del); 
  xdbar = 0.5*(udb + del);
  xs = sb;
  xc = chm;
  xb = bot;

} 

//*********

double GRV94L::grvv (double x, double n, double ak, double bk, double a, 
   double b, double c, double d) {

  double dx = sqrt(x);
  return n * pow(x, ak) * (1. + a * pow(x, bk) + x * (b + c * dx)) *
    pow(1. - x, d);

} 

//*********

double GRV94L::grvw (double x, double s, double al, double be, double ak, 
  double bk, double a, double b, double c, double d, double e, double es) {
 
  double lx = log(1./x);
  return (pow(x, ak) * (a + x * (b + x * c)) * pow(lx, bk) + pow(s, al)
    * exp(-e + sqrt(es * pow(s, be) * lx))) * pow(1. - x, d);

}

//*********
  
double GRV94L::grvs (double x, double s, double sth, double al, double be, 
  double ak, double ag, double b, double d, double e, double es) {
 
  if(s <= sth) {
    return 0.;
  } else {
    double dx = sqrt(x);
    double lx = log(1./x);
    return pow(s - sth, al) / pow(lx, ak) * (1. + ag * dx + b * x) *
      pow(1. - x, d) * exp(-e + sqrt(es * pow(s, be) * lx));
  }
 
}

//*********

// Gives electron (or muon, or tau) parton distribution.
 
void Lepton::xfUpdate(double x, double Q2) {
 
  // Common constants.
  double alphaEM = 0.00729735;
  double m2 = pow(Mass(idBeam), 2);

  // Electron inside electron, see R. Kleiss et al., in Z physics at
  // LEP 1, CERN 89-08, p. 34
  double xLog = log(max(1e-10,x));
  double xMinusLog = log( max(1e-10, 1. - x) );
  double Q2Log = log( max(3., Q2/m2) );
  double beta = (alphaEM / M_PI) * (Q2Log - 1.);
  double delta = 1. + (alphaEM / M_PI) * (1.5 * Q2Log + 1.289868) 
    + pow(alphaEM / M_PI, 2) * (-2.164868 * Q2Log*Q2Log 
    + 9.840808 * Q2Log - 10.130464);
  double fPrel =  beta * pow(1. - x, beta - 1.) * sqrtpos(delta)
     - 0.5 * beta * (1. + x) + 0.125 * beta*beta * ( (1. + x) 
     * (-4. * xMinusLog + 3. * xLog) - 4. * xLog / (1. - x) - 5. - x); 

  // Zero distribution for very large x and rescale it for intermediate.
  if (x > 1. - 1e-10) fPrel = 0.;
  else if (x > 1. - 1e-7) fPrel *= pow(1000., beta) / (pow(1000., beta) - 1.); 
  xlepton = x * fPrel; 

  // Photon inside electron (one possible scheme - primitive).
  xgamma = (0.5 * alphaEM / M_PI) * Q2Log * (1. + pow(1. - x, 2));

}

ThePEGPDF::ThePEGPDF(tcPDFPtr pdf, tcPDPtr parent)
    : PDF(parent->id()), thePDF(pdf), theParent(parent), theLeptonID(0) {
  cPDVector partons = thePDF->partons(parent);
  for ( int i = 0, N = partons.size(); i < N; ++i )
    thePartons[partons[i]->id()] = partons[i];
  if ( LeptonMatcher::Check(idBeam) &&
       thePartons.find(idBeam) != thePartons.end() ) theLeptonID = idBeam;
}

void ThePEGPDF::set(long id, double xf) {
  using namespace ParticleID;
  if ( idBeam < 0 && ( abs(id) == u || abs(id) == d ) ) id = -id;
  if ( id == theLeptonID ) xlepton = xf;
  else switch ( id ) {
  case u:
    xu = xf;
    break;
  case ubar:
    xubar = xf;
    break;
  case d:
    xd = xf;
    break;
  case dbar:
    xdbar = xf;
    break;
  case ParticleID::s:
  case sbar:
    xs = xf;
    break;
  case c:
  case cbar:
    xc = xf;
    break;
  case b:
  case bbar:
    xb = xf;
    break;
  case ParticleID::g:
    xg = xf;
    break;
  case ParticleID::gamma:
    xgamma = xf;
  }
}

void ThePEGPDF::xfUpdate(double x, double Q2) {
  for ( map<long, tcPDPtr>::const_iterator it = thePartons.begin();
	it != thePartons.end(); ++it )
    set(it->first, thePDF->xfx(theParent, it->second, Q2*GeV2, x));
}

