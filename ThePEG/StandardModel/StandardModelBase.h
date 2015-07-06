// -*- C++ -*-
//
// StandardModelBase.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_StandardModelBase_H
#define ThePEG_StandardModelBase_H
// This is the declaration of the StandardModelBase class.

#include "ThePEG/Config/ThePEG.h"
#include "AlphaEMBase.h"
#include "CKMBase.h"
#include "AlphaSBase.h"
// #include "StandardModelBase.fh"
// #include "StandardModelBase.xh"

namespace ThePEG {

/**
 * StandardModelBase is used to handle standard model parameters in an
 * EventGenerator. It uses AlphaEMBase, AlphaSBase and CKMBase to help
 * with the implementation of the electro-magnetic and QCD couplings
 * and the flavour mixing matrix. This means that StandardModelBase
 * does not need to be inherited from when it comes to standard model
 * parameters. Beyond the standard model parameters should be
 * implemented as sub-classes.
 *
 * @see \ref StandardModelBaseInterfaces "The interfaces"
 * defined for StandardModelBase.
 * @see EventGenerator
 * @see AlphaEMBase
 * @see AlphaSBase
 * @see CKMBase
 */
class StandardModelBase: public Interfaced {

  /** Declare a pointer to an AlphaEMBase object. */
  typedef Ptr<AlphaEMBase>::pointer AEMPtr;
  /** Declare a pointer to an AlphaSBase object. */
  typedef Ptr<AlphaSBase>::pointer ASPtr;
  /** Declare a pointer to n CKMBase object. */
  typedef Ptr<CKMBase>::pointer CKMPtr;
  /** Declare a transient pointer to an AlphaEMBase object. */
  typedef Ptr<AlphaEMBase>::transient_pointer tAEMPtr;
  /** Declare a transient pointer to an AlphaSBase object. */
  typedef Ptr<AlphaSBase>::transient_pointer tASPtr;
  /** Declare a transient pointer to a CKMBase object. */
  typedef Ptr<CKMBase>::transient_pointer tCKMPtr;

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * Default constructor.
   */
  StandardModelBase();

  /**
   * Destructor.
   */
  virtual ~StandardModelBase();
  //@}

public:

  /**
   * Return the number of families assumed in the standard model.
   */
  unsigned int families() const { return theFamilies; }

public:


  /** @name Functions accessing electro-weak parameters. */
  /**
   *  Return the electroweak scheme used
   */
  unsigned int ewScheme() const { return theElectroWeakScheme; }

  /**
   *  Set the electroweak scheme used
   */
  void ewScheme(unsigned int s) { theElectroWeakScheme = s; }

  //@{
  /**
   * Constant \f$\alpha_{EM}(q^2=0)\f$.
   */
  double alphaEM() const { return theAlphaEM; }

  /**
   * Constant \f$\alpha_{EM}(q^2=m_Z^2)\f$.
   */
  double alphaEMMZ() const { return theAlphaEMMZ; }

  /**
   *  The electromagnetic coupling for vertex classes
   *  in a well defined self-consistent EW scheme if requested
   */
  double alphaEMME(Energy2 scale) const {
    if(theElectroWeakScheme==0)
      return alphaEM(scale);
    else if(scale>1e-6*GeV2)
      return theAlphaEMMZ;
    else
      return theAlphaEM;
  }

  /**
   * Running \f$\alpha_{EM}\f$.
   */
  double alphaEM(Energy2 scale) const {
    return theRunningAlphaEM->value(scale, *this);
  }

  /**
   * Return a pointer to the object handling \f$\alpha_{EM}\f$.
   */
  tAEMPtr alphaEMPtr() const { return theRunningAlphaEM; }

  /**
   * Return \f$\sin^2(\theta_W)\f$.
   */
  double sin2ThetaW() const { return theSin2ThetaW; }

  /**
   *  The Fermi constant
   */
  InvEnergy2 fermiConstant() const {return theGF;}

  /**
   * The neutrino-photon coupling.
   */
  double enu() const { return theEnu; }

  /**
   * The charged lepton-photon coupling.
   */
  double ee() const { return theEe; }

  /**
   * The up-type-photon coupling.
   */
  double eu() const { return theEu; }

  /**
   * The down-type-photon coupling.
   */
  double ed() const { return theEd; }

  /**
   * The vector neutrino-\f$Z^0\f$ coupling.
   */
  double vnu() const { return theVnu; }

  /**
   * The vector charged lepton-\f$Z^0\f$ coupling.
   */
  double ve() const { return theVe; }

  /**
   * The vector up-type-\f$Z^0\f$ coupling.
   */
  double vu() const { return theVu; }

  /**
   * The vector down-type-\f$Z^0\f$ coupling.
   */
  double vd() const { return theVd; }

  /**
   * The axial neutrino-\f$Z^0\f$ coupling.
   */
  double anu() const { return theAnu; }

  /**
   * The axial charged lepton-\f$Z^0\f$ coupling.
   */
  double ae() const { return theAe; }

  /**
   * The axial up-type-\f$Z^0\f$ coupling.
   */
  double au() const { return theAu; }

  /**
   * The axial down-type-\f$Z^0\f$ coupling.
   */
  double ad() const { return theAd; }

  /**
   * Return a pointer to the CKMBase object used.
   */
  tCKMPtr CKM() const { return theCKM; }

  /**
   * Return a square of the element of the Cabibbo-Kobayashi-Maskawa
   * Matrix. The fatrix element for the \a uf up-type family and \a df
   * down-type family.
   */
  double CKM(unsigned int uf, unsigned int df) const;

  /**
   * Return the square of the elements of the Cabibbo-Kobayashi-Maskawa
   * Matrix.
   */
  double CKM(const ParticleData & uType,
		    const ParticleData & dType) const;
  //@}

public:

  /** @name Functions accessing QCD parameters. */
  //@{
  /**
   * Return the number of colours.
   */
  unsigned int Nc() const { return theNc; }

  /**
   * Return the number of avtive quark flavours for a given \a scale.
   */
  unsigned int Nf(Energy2 scale) const {
    return theRunningAlphaS->Nf(scale);
  }

  /**
   * Return the constant strong coupling constant.
   */
  double alphaS() const { return theAlphaS; }

  /**
   * Return the running strong coupling for a given \a scale
   */
  double alphaS(Energy2 scale) const {
    return theRunningAlphaS->value(scale, *this);
  }

  /**
   * Return a pointer to the object handling \f$\alpha_S\f$.
   */
  tASPtr alphaSPtr() const {
    return theRunningAlphaS;
  }

  /**
   * Return the \f$\Lambda_{QCD}\f$ for \a nflav active flavours.
   */
  Energy LambdaQCD(unsigned int nflav) const {
    return theRunningAlphaS->LambdaQCD(nflav);
  }

  /**
   * Return the \f$\Lambda_{QCD}\f$ for the given \a scale.
   */
  Energy LambdaQCD(Energy2 scale) const { return LambdaQCD(Nf(scale)); }
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
   * Standard Init function used to initialize the interface.
   */
  static void Init();

  /**
   * Overloaded function from Interfaced
   */
  virtual bool preInitialize() const {
    return true;
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


protected:

  /** @name Standard Interfaced functions. */
  //@{

  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();
  //@}

private:

  /**
   * The number of families.
   */
  unsigned int theFamilies;

  /**
   * The constant \f$\alpha_{EM}\f$.
   */
  double theAlphaEM;

  /**
   * The constant \f$\alpha_{EM}\f$.
   */
  double theAlphaEMMZ;

  /**
   * Pointer to an object capable of calculating the running
   * \f$\alpha_{EM}\f$.
   */
  AEMPtr theRunningAlphaEM;

  /**
   * The \f$\sin^2(\theta_W)\f$
   */
  double theSin2ThetaW;

  /**
   * The Fermi contants \f$G_F\f$
   */
  InvEnergy2 theGF;

  /**
   * Coupling between a fundamental fermion and the photon.
   */
  double theEnu;

  /**
   * Coupling between a fundamental fermion and the photon.
   */
  double theEe;

  /**
   * Coupling between a fundamental fermion and the photon.
   */
  double theEu;

  /**
   * Coupling between a fundamental fermion and the photon.
   */
  double theEd;

  /**
   * Vector coupling between a fundamental fermion and Z^0.
   */
  double theVnu;

  /**
   * Vector coupling between a fundamental fermion and Z^0.
   */
  double theVe;

  /**
   * Vector coupling between a fundamental fermion and Z^0.
   */
  double theVu;

  /**
   * Vector coupling between a fundamental fermion and Z^0.
   */
  double theVd;

  /**
   * Axial coupling between a fundamental fermions and Z^0.
   */
  double theAnu;

  /**
   * Axial coupling between a fundamental fermions and Z^0.
   */
  double theAe;

  /**
   * Axial coupling between a fundamental fermions and Z^0.
   */
  double theAu;

  /**
   * Axial coupling between a fundamental fermions and Z^0.
   */
  double theAd;

  /**
   * If true, the electro-weak couplings are derived from
   * \f$\theta_W\f$ in the initialization.
   */
  long recalculateEW;

  /**
   * A pointer to an object representing the Cabibbo-Kobayashi-Maskawa
   * matrix.
   */
  CKMPtr theCKM;

  /**
   * The matrix of squared CKM elements set from theCKM at initialization.
   */
  mutable vector< vector<double> > theCKM2Matrix;

  /**
   * The number of colours;
   */
  unsigned int theNc;

  /**
   * The fixed strong coupling.
   */
  double theAlphaS;

  /**
   * Pointer to an object capable of calculating the running
   * \f$\alpha_{S}\f$.
   */
  ASPtr theRunningAlphaS;

  /**
   *  Electroweak scheme
   */
  unsigned int theElectroWeakScheme;

  /**
   *  Option for the calculation of the W/Z widths
   */
  unsigned int theBosonWidthOption;

private:

  /**
   * Describe a concrete class with persistent data.
   */
  static ClassDescription<StandardModelBase> initStandardModelBase;

  /**
   *  Private and non-existent assignment operator.
   */
  StandardModelBase & operator=(const StandardModelBase &);

};

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the base classes
 *  of StandardModelBase. */
template <>
struct BaseClassTrait<StandardModelBase,1>: public ClassTraitsType {
  /** Typedef of the first base class of StandardModelBase. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of the
 *  StandardModelBase class. */
template <>
struct ClassTraits<StandardModelBase>:
    public ClassTraitsBase<StandardModelBase> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::StandardModelBase"; }
};

/** @endcond */

}

#endif /* ThePEG_StandardModelBase_H */
