// -*- C++ -*-
#ifndef THEPEG_CMSlrnsac_H
#define THEPEG_CMSlrnsac_H
//
// This is the declaration of the CMSlrnsac class.
//

#include "ThePEG/Handlers/AnalysisHandler.h"
#include "ThePEG/Utilities/CFileLineReader.h"

namespace ThePEG {

/**
 * Here is the documentation of the CMSlrnsac class.
 *
 * @see \ref CMSlrnsacInterfaces "The interfaces"
 * defined for CMSlrnsac.
 */
class CMSlrnsac: public AnalysisHandler {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  CMSlrnsac();

  /**
   * The copy constructor.
   */
  CMSlrnsac(const CMSlrnsac &);

  /**
   * The destructor.
   */
  virtual ~CMSlrnsac();
  //@}

public:

  /** @name Virtual functions required by the AnalysisHandler class. */
  //@{
  /**
   * Analyze the given vector of particles. The default version calls
   * analyze(tPPtr) for each of the particles.
   * @param particles the vector of pointers to particles to be analyzed
   * @param weight the event weight
   */
  virtual void analyze(const tPVector & particles, double weight);
  //@}

  /**
   * Find the correct pt-bin
   */
  static int ptbin(Energy pt) {
    if ( pt < 1.0*GeV ) return 0;
    if ( pt < 2.0*GeV ) return 1;
    if ( pt < 3.0*GeV ) return 2;
    if ( pt < 4.0*GeV ) return 3;
    return 4;
  }

  /**
   * Find the correct multiplicity bin.
   */
  static int nchbin(int n) {
    if ( n < 35 ) return 0;
    if ( n < 90 ) return 1;
    if ( n < 110 ) return 2;
    return 3;
  }

  static double nozero(double x) {
    return x == 0.0? 1.0: x;
  }

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
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  virtual void dofinish();
  //@}

private:

  vector< vector<tH1DPtr> > Sphi;
  vector< vector<tH1DPtr> > Bphi;
  vector< vector<tH2DPtr> > S2phi;
  vector< vector<tH2DPtr> > B2phi;
  vector< vector<tH1DPtr> > S3phi;
  vector< vector<tH1DPtr> > B3phi;
  vector< vector<tH1DPtr> > multbin;
  vector< vector<tH1DPtr> > multbin2;
  tH2DPtr SMB01;
  tH2DPtr BMB01;
  tH2DPtr SMB13;
  tH2DPtr BMB13;
  tH2DPtr SHN01;
  tH2DPtr BHN01;
  tH2DPtr SHN13;
  tH2DPtr BHN13;
  tH2DPtr S2MB01;
  tH2DPtr B2MB01;
  tH2DPtr S2MB13;
  tH2DPtr B2MB13;
  tH2DPtr S2HN01;
  tH2DPtr B2HN01;
  tH2DPtr S2HN13;
  tH2DPtr B2HN13;
  tH1DPtr multMB01;
  tH1DPtr multMB13;
  tH1DPtr multHN01;
  tH1DPtr multHN13;

  vector< vector< vector<double> > > refphi;
  vector< vector< vector<double> > > refeta;
  vector<int> refNch04;
  vector<double> refweight;
  vector<double> sumrefweight;
  vector<double> sumweight;
  tH1DPtr mult;
  vector<tH1DPtr> bdist;
  double etamax;

  int doRweight;
  string dir;

  string externalFile;
  CFileLineReader file;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  CMSlrnsac & operator=(const CMSlrnsac &);

};

}

#endif /* THEPEG_CMSlrnsac_H */
