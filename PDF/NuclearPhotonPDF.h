// -*- C++ -*-
#ifndef ThePEG_NuclearPhotonPDF_H
#define ThePEG_NuclearPhotonPDF_H
//
// This is the declaration of the NuclearPhotonPDF class.
//

#include "ThePEG/PDF/PDFBase.h"
#include "ThePEG/Utilities/Maths.h"

namespace ThePEG {

/**
 * The NuclearPhotonPDF class implements the PDF for the radiation of a photon from a heavy ion.
 *
 * @see \ref NuclearPhotonPDFInterfaces "The interfaces"
 * defined for NuclearPhotonPDF.
 */
class NuclearPhotonPDF: public PDFBase {

public:

  /**
   * The default constructor.
   */
  NuclearPhotonPDF();

public:

  /** @name Virtual functions to be overridden by sub-classes. */
  //@{
  /**
   * Return true if this PDF can handle the extraction of partons from
   * the given \a particle.
   */
  virtual bool canHandleParticle(tcPDPtr particle) const;

  /**
   * Return the partons which this PDF may extract from the given
   * \a particle.
   */
  virtual cPDVector partons(tcPDPtr particle) const;

  /**
   * The density. Return the pdf for the given \a parton inside the
   * given \a particle for the virtuality \a partonScale and
   * logarithmic momentum fraction \a l \f$(l=\log(1/x)\f$. The \a
   * particle is assumed to have a virtuality \a particleScale.
   */
  virtual double xfl(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
		     double l, Energy2 particleScale = 0.0*GeV2) const;

  /**
   * The valence density. Return the pdf for the given cvalence \a
   * parton inside the given \a particle for the virtuality \a
   * partonScale and logarithmic momentum fraction \a l
   * \f$(l=\log(1/x)\f$. The \a particle is assumed to have a
   * virtuality \a particleScale. If not overidden by a sub class this
   * will return zero.
   */
  virtual double xfvl(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
		     double l, Energy2 particleScale = 0.0*GeV2) const;

  /**
   * Generate scale (as a fraction of the maximum scale). If the PDF
   * contains strange peaks which can be difficult to handle, this
   * function may be overwritten to return an appropriate scale
   * \f$Q^2/Q^2_{\max}\f$ for a \a z uniformly distributed in
   * ]0,1[. Also the jacobobian of the \f$Q^2/Q^2_{\max}\rightarrow
   * z\f$ variable transformation must multiply the \a jacobian
   * argument. The default version will simply use the function
   * \f$Q^2/Q^2_{\max} = (Q^2_{\max}/Q^2_{\min})^(z-1)\f$ or, if
   * \f$Q^2_{\min}\f$ is zero, \f$Q^2/Q^2_{\max} = z\f$ (where the
   * limits are set by \a cut).
   */
  
  virtual double flattenScale(tcPDPtr particle, tcPDPtr parton,
			       const PDFCuts & cut, double l, double z,
			       double & jacobian) const;

  
  /**
   * Generate a momentum fraction. If the PDF contains strange peaks
   * which can be difficult to handle, this function may be
   * overwritten to return an appropriate \f$l=\log(1/x)\f$ for a \a z
   * uniformly distributed in ]0,1[. Also the jacobobian of the
   * \f$l\rightarrow z\f$ variable transformation must in the function
   * multiply the \a jacobian argument. The default version will
   * simply use the function \f$l(z) = l_{\min} +
   * z*(l_{\max}-l_{\min})\f$ (where the limits are set by \a cut).
   */
  virtual double flattenL(tcPDPtr particle, tcPDPtr parton, const PDFCuts &cut,
			  double z, double & jacobian) const;
  //@}


public:

  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void persistentInput(PersistentIStream & is, int version);
  //@}

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

protected:

  /**
   *  Extract the mass and atomic number for the nuclear PDG code
   */
  pair<int,int> massAndZ(long pid) const {
    pair<int,int> output;
    output.first  = (abs(pid)%10000)/10;
    output.second = (abs(pid)%10000000)/10000;
    return output;
  }

  /**
   *  Dipole form factor
   */
  double dipoleFormFactor(Energy2 q2) const {
    return Math::powi((1 + q2/q02_),-2);
  }

  /**
   *  Convolution of hard sphere and Yukaw for heavy nuclei
   */
  double heavyFormFactor(Energy2 q2, Length R) const {
    double Rq = sqrt(q2)*R/hbarc;
    return 3./pow(Rq,3)*(sin(Rq)-Rq*cos(Rq))/(1.+sqr(yukawaRange_/hbarc)*q2);
  }
  
protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  NuclearPhotonPDF & operator=(const NuclearPhotonPDF &) = delete;
  
private:

  /**
   *  Minimum \f$Q^2\f$ for the photon
   */
  Energy2 q2Min_;

  /**
   *  Maximum \f$Q^2\f$ for the photon
   */
  Energy2 q2Max_;

  /**
   *  Range of the Yukawa potential
   */
  Length yukawaRange_;

  /**
   *  Prefactor for the nucleur radius \f$R_A=aA^{\frac13}\f$
   */
  Length aFact_;

  /**
   *  Fitted scale \f$Q{_0}{^2}=0.71GeV^2\f$ for the dipole form factor
   */
  Energy2 q02_;
  
};

}

#endif /* ThePEG_NuclearPhotonPDF_H */
