// Header file for information on incoming beams, including PDF's.

#ifndef Beam_H
#define Beam_H

#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>
#include "Pythia7/Config/Pythia7.h"
#include "ThePEG/PDF/PDFBase.h"

namespace Pythia7 {
namespace Shower {

using namespace std;

//**************************************************************************

/**
 * Base class for parton distribution functions.
 * Used by the internal Pythia7 Shower classes.
 */
class PDF {

public:
  /** NOT DOCUMENTED */
  PDF(long idBeamin = 2212)
    : idBeam(idBeamin), xSav(-1.0), Q2Sav(-1.0), xu(0.0), xd(0.0), xubar(0.0),
      xdbar(0.0), xs(0.0), xc(0.0), xb(0.0), xg(0.0), xlepton(0.0),
      xgamma(0.0) {}
  /** NOT DOCUMENTED */
  virtual ~PDF() {}
  /** NOT DOCUMENTED */
  double xfx(long, double, double);  
  /** NOT DOCUMENTED */
  virtual bool isParton(long) const; 
  /** NOT DOCUMENTED */
  virtual bool isValence(long) const; 
  
protected:
  /** @cond NEVERTOBEDOCUMENTED */

  /** NOT DOCUMENTED */
  long idBeam;
  /** NOT DOCUMENTED */
  double xSav, Q2Sav;
  /** NOT DOCUMENTED */
  double xu, xd, xubar, xdbar, xs, xc, xb, xg, xlepton, xgamma;
  /** NOT DOCUMENTED */
  virtual void xfUpdate(double x, double Q2) = 0; 

  /** @endcond */

};

//*********


/**
 * Gives the GRV 94 L (leading order) parton distribution function set
 * in parametrized form. Authors: M. Glueck, E. Reya and A. Vogt.
 * Used by the internal Pythia7 Shower classes.
 */
class GRV94L : public PDF {

public:
  /** NOT DOCUMENTED */
  GRV94L(long idBeamin = 2212) : PDF(idBeamin) {;}

private:
  /** NOT DOCUMENTED */
  void xfUpdate(double x, double Q2);
  /** NOT DOCUMENTED */
  double grvv (double x, double n, double ak, double bk, double a, 
    double b, double c, double d);
  /** NOT DOCUMENTED */
  double grvw (double x, double s, double al, double be, double ak, 
    double bk, double a, double b, double c, double d, double e, double es);
  /** NOT DOCUMENTED */
  double grvs (double x, double s, double sth, double al, double be, 
    double ak, double ag, double b, double d, double e, double es);
};

//*********

 
/**
 * Gives electron (or muon, or tau) parton distribution.
 * Used by the internal Pythia7 Shower classes.
 */
class Lepton : public PDF {

public:
  /** NOT DOCUMENTED */
  Lepton(long idBeamin = 11) : PDF(idBeamin) {;}

private:
  /** NOT DOCUMENTED */
  void xfUpdate(double x, double Q2);
};

//**************************************************************************


/**
 * This class holds info on a beam particle in the shower evolution.
 * It derives from the particle class.
 */
class BeamParticle : public Particle {

public:
  /**
   * Constructors.
   */
  BeamParticle() : Particle() {pdfbm = 0;}  
  /** NOT DOCUMENTED */
  BeamParticle(long idin, Vec4 pin = Vec4(0.,0.,0.,0.), double min = 0., 
    PDF* pdfin = 0) :  Particle(idin, -9, -1, -1, 0, 0, pin, min, 0.) 
    {pdf(pdfin);} 
  /** NOT DOCUMENTED */
  BeamParticle(const BeamParticle& pt) : Particle(pt) {pdfbm = pt.pdfbm;}  
  /** NOT DOCUMENTED */
  BeamParticle& operator=(const BeamParticle& pt) {if (this != &pt) {
    pdfbm = pt.pdfbm; return *this;} } 
  /**
   * Construct a BeamParticle from a Particle.
   */
  BeamParticle(const Particle& pt) : Particle(pt) {pdfbm = 0;}  
  /** NOT DOCUMENTED */
  BeamParticle& operator=(const Particle& pt) {Particle::operator=(pt);   
    pdfbm = 0; return *this; }  


  /**
   * Member functions.
   */
  double xfx(long id, double x, double Q2) {return pdfbm->xfx(id, x, Q2);}
  /** NOT DOCUMENTED */
  void pdf(PDF* pdfin) {pdfbm = pdfin;}
  /** NOT DOCUMENTED */
  PDF* pdf() const {return pdfbm;}
  /** NOT DOCUMENTED */
  bool isParton(long id) const {return pdfbm->isParton(id);} 
  /** NOT DOCUMENTED */
  bool isValence(long id) const {return pdfbm->isValence(id);} 
  /** NOT DOCUMENTED */
  bool isHadron() const {if (abs(id()) > 100) return true; return false;}
  /** NOT DOCUMENTED */
  bool isLepton() const {if (abs(id()) == 11 || abs(id()) == 13
    || abs(id()) == 15) return true; return false;}


private: 
  /** NOT DOCUMENTED */
  PDF* pdfbm;
};

/**
 * Interface class to make the ThePEG pdf classes available to the
 * internal Pythia7 Shower classes.
 */
class ThePEGPDF: public PDF {

public:

  /**
   * Constructor taking a ThePEG::PDFBase object and a
   * ThePEG::ParticleData as arguments.
   */
  ThePEGPDF(tcPDFPtr pdf, tcPDPtr parent);

  /**
   * Calculate new values of the parton densities to be cached.
   */
  virtual void xfUpdate(double x, double Q2);

  /** Set the cached value of the density \a xf for parton type \a id. */
  void set(long id, double xf);

private:

  /** The ThePEG pdf object. */
  tcPDFPtr thePDF;
  /** The ThePEG particle type. */
  tcPDPtr theParent;
  /** Map of possible partons indexed by their id. */
  map<long, tcPDPtr> thePartons;
  /** The type of lepton if this is a lepton pdf. */
  long theLeptonID;
  

};



//**************************************************************************
}
}

#endif // Beam_H
