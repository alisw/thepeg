// -*- C++ -*-
#ifndef DIPSY_FSAnalysis_H
#define DIPSY_FSAnalysis_H
//
// This is the declaration of the FSAnalysis class.
//

#include "ThePEG/Handlers/AnalysisHandler.h"
#include <iostream>
#include <fstream>

namespace DIPSY {

using namespace ThePEG;

/**
 * Here is the documentation of the FSAnalysis class.
 *
 * @see \ref FSAnalysisInterfaces "The interfaces"
 * defined for FSAnalysis.
 */
class FSAnalysis: public AnalysisHandler {

public:

  /**
   * Internal class to define the rapidity ordering of particles
   */
  struct POrder {
    /** The actual ordering function. */
    bool operator()(tPPtr a, tPPtr b) const {
    double eta1 = a->eta();
    if ( abs(a->eta()) > 1000.0 ) {
      if ( a->momentum().z() > ZERO ) eta1 = 10.0;
      else eta1 = -10.0;
    }
    double eta2 = b->eta();
    if ( abs(b->eta()) > 1000.0 ) {
      if ( b->momentum().z() > ZERO ) eta2 = 10.0;
      else eta2 = -10.0;
    }
      return eta1 < eta2;
    }
  };

  /**
   * A rapidity-ordered set of particles.
   */
  typedef set<PPtr,POrder> PSet;

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  FSAnalysis();

  /**
   * The destructor.
   */
  virtual ~FSAnalysis();
  //@}

public:

  /** @name Virtual functions required by the AnalysisHandler class. */
  //@{
  /**
   * Analyze a given Event. Note that a fully generated event
   * may be presented several times, if it has been manipulated in
   * between. The default version of this function will call transform
   * to make a lorentz transformation of the whole event, then extract
   * all final state particles and call analyze(tPVector) of this
   * analysis object and those of all associated analysis objects. The
   * default version will not, however, do anything on events which
   * have not been fully generated, or have been manipulated in any
   * way.
   * @param event pointer to the Event to be analyzed.
   * @param ieve the event number.
   * @param loop the number of times this event has been presented.
   * If negative the event is now fully generated.
   * @param state a number different from zero if the event has been
   * manipulated in some way since it was last presented.
   */
  virtual void analyze(tEventPtr event, long ieve, int loop, int state);

  /**
   * Transform the event to the desired Lorentz frame and return the
   * corresponding LorentzRotation.
   * @param event a pointer to the Event to be transformed.
   * @return the LorentzRotation used in the transformation.
   */
  virtual LorentzRotation transform(tEventPtr event) const;

  /**
   * Analyze the given vector of particles. The default version calls
   * analyze(tPPtr) for each of the particles.
   * @param particles the vector of pointers to particles to be analyzed
   */
  virtual void analyze(const tPVector & particles);

  /**
   * Analyze the given particle.
   * @param particle pointer to the particle to be analyzed.
   */
  virtual void analyze(tPPtr particle);
  //@}

protected:

  /**
   * The definition of the reaction-plane eccentricity parameter.
   * PLB641 (2006) 260, Eq. 6.
   */
  double reactionPlaneEccentricity(const vector<cPPtr> partons);

/**
 * The definition of the participant eccentricity parameter.
 * PLB641 (2006) 260, Eq. 2.
 */
  double participantEccentricity(const vector<cPPtr> partons);

/**
 * The definition of the average overlap area.
 * PRC81, 024901 (2010), Eq. 11.
 */
  Area averageOverlapArea(const vector<cPPtr> partons);
  Area averageOverlapArea(vector < vector<double> > gluons);

/**
 * Analyses the transverse position of the t = 0 gluons. (v_2 pp paper)
 */
  void v2GluonAnalyze(tEventPtr event);

/**
 * Finds the eccentricities. (v_n AA paper)
 */
  void eccentricities(tEventPtr event);

/**
 * Reads the event from a file with the provided filename (full path to file).
 * Returns a vector of vectors, each subvector being the gluon information on the
 * format
 * particle_number, transverse_position_x(fm), transverse_position_y(fm),
 * transverse_momentum_x(GeV), transverse_momentum_y(GeV),
 * rapidity, colour_neighbour_number, anticolour_neighbour_number
 * 
 */
  vector < vector<double> > extractEvent(string filename);

/**
 * extracts information from an input file of the above format.
 */
  pair<double, double> extractB(string filename);
  double extractWeight(string filename);
  int extractNSpectators(string filename);

/**
 * changes the weight in the file by a factor.
 */
  void changeWeight(string filename, double factor);

/**
 * extracts information from a line of the format of the file above.
 */
  pair<double, double> getB(char line[256]);
  double getWeight(char line[256]);
  int getNSpectators(char line[256]);

/**
 * takes a char line from a file as above, and converts it to a vector of doubles.
 */
  vector<double> getValues(char line[256]);

/**
 * Finds and returns the centre of gravity.
 */
  pair<double, double> findCoG(vector < vector<double> > gluons);

/**
 * Finds and returns the eccentricity for moment n for the provided
 * participants and centre of gravity.
 * the first version use <r^2 sin(3phi)>, the second <r^N sin(3phi)>.
 * The third weights gluons according to max(1,ln(p_T)+1) in the averages.
 * the last takes the two options as argument.
 */
  double eccentricity(int n, vector < vector<double> > gluons, pair<double, double> CoG);
  double eccentricityN(int n, vector < vector<double> > gluons, pair<double, double> CoG);
  double eccentricityW(int n, vector < vector<double> > gluons, pair<double, double> CoG);
  double eccentricity(int n, vector < vector<double> > gluons, pair<double, double> CoG, bool always2, bool weighted);

/**
 * Finds and returns the participant plane for moment n for the provided
 * participants and centre of gravity.
 * the first version use <r^2 sin(3phi)>, the second <r^N sin(3phi)>.
 * The third weights gluons according to max(1,ln(p_T)+1) in the averages.
 * the last takes the two options as arguments.
 */
  double participantAngle(int n, vector < vector<double> > gluons, pair<double, double> CoG);
  double participantAngleN(int n, vector < vector<double> > gluons, pair<double, double> CoG);
  double participantAngleW(int n, vector < vector<double> > gluons, pair<double, double> CoG);
  double participantAngle(int n, vector < vector<double> > gluons, pair<double, double> CoG, bool always2, bool weighted);

/**
 * Saves the time=0 gluons to file.
 */
  void saveGluonsToFile(tEventPtr event);

/**
 * Saves the final state in bins of x,y,eta using the SPheRIO format.
 */
  void exportToSPheRIO(tEventPtr event);


/**
 * Counts and returns the number of spectating nucleons in a HI collision.
 */
  int numberOfSpectators(tEventPtr event);

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

  /** The multiplicity distribution for 1-2 interactions. */
  tH1DPtr mult12;

  /** The reaction plane eccentricity. */
  tH1DPtr epsRP;

  /** The participant eccentricity. */
  tH1DPtr epsPart;

  /** The participant eccentricity squared. */
  tH1DPtr epsPartOver20;

  /** The participant eccentricity to the fourth power. */
  tH1DPtr epsPartOver40;

  /** The distribution of the weights. */
  tH1DPtr weights;

  /** The number of gluons before final state shower. */
  tH1DPtr ngbefore;

  /** The number of gluons after final state shower. */
  tH1DPtr ngafter;

  /** The PT distribution in rapidity for DIPSY. */
  tH1DPtr DIPSYPT;

  /** The PT distribution in rapidity for final state. */
  tH1DPtr finalPT;

  /** The PT distribution in pseudorapidity for final state. */
  tH1DPtr finalPTeta;

  /** The PT distribution in pseudorapidity for final state. */
  tH1DPtr finalChargedPTeta;

  /** The ET distribution in pseudorapidity for final state. */
  tH1DPtr finalETeta;

  /** The Nch distribution in pseudorapidity for final state. */
  tH1DPtr finalNcheta;

  /** The distribution in excited mass dN/dln(M_X^2) for single diffraction. */
  tH1DPtr lnMX2dist;
  tH1DPtr lnMX2distFrac;
  tH1DPtr lnMX2dist0;
  tH1DPtr lnMX2dist1;
  tH1DPtr lnMX2dist2;
  tH1DPtr lnMX2distWeights;
  tH1DPtr lnMX2dist0Weights;
  tH1DPtr lnMX2dist1Weights;
  tH1DPtr lnMX2dist2Weights;

  /** The Nch distribution in pseudorapidity for single diffraction. */
  tH1DPtr SDNchdeta;

  /** The rapidity range versus MX and MX^2. */
  tH1DPtr MXmaxY;
  tH1DPtr MXmaxYWeights;
  tH1DPtr MX2maxY;
  tH1DPtr MX2maxYWeights;

  /** UA4 SD analysis */
  tH1DPtr UA4dndetaMX20;
  tH1DPtr UA4dndetaMX20Weights;
  tH1DPtr UA4dndetaMX60;
  tH1DPtr UA4dndetaMX60Weights;
  tH1DPtr UA4dndetaMX80;
  tH1DPtr UA4dndetaMX80Weights;
  tH1DPtr UA4dndetaMX100;
  tH1DPtr UA4dndetaMX100Weights;
  tH1DPtr UA4dndetaMX115;
  tH1DPtr UA4dndetaMX115Weights;
  tH1DPtr UA4dndetaMX140;
  tH1DPtr UA4dndetaMX140Weights;


  /** The Nch distribution in rapidity (in excited CoM frame) for single diffraction.
   * The MX is the average MX in GeV in the bins 3-8, 8-15 and 15-30 respectively.
   */
  tH1DPtr SDNchdyMX38;
  tH1DPtr SDNchdyMX815;
  tH1DPtr SDNchdyMX1530;
  tH1DPtr SDNchdyMX38Weights;
  tH1DPtr SDNchdyMX815Weights;
  tH1DPtr SDNchdyMX1530Weights;

  /** The energy distribution in rapidity (in excited CoM frame) for single diffraction.
   * The MX is the average MX in GeV in the bins 3-8, 8-18 and 18-30 respectively.
   */
  tH1DPtr SDdedetaMX38;
  tH1DPtr SDdedetaMX818;
  tH1DPtr SDdedetaMX1830;
  tH1DPtr SDdedetaMX38Weights;
  tH1DPtr SDdedetaMX818Weights;
  tH1DPtr SDdedetaMX1830Weights;

  /** The weigth distribution in diffractive events. */
  tH1DPtr SDweights;
  tH1DPtr AU4SDweights;

  /** The number of charged particles above eta = -3 in SD. */
  tH1DPtr SDavN;
  tH1DPtr SDavNWeights;
  tH1DPtr SDavN2;
  tH1DPtr SDavN2Weights;

  /** The PT distribution in rapidity for after FSR, but before hadronisation. */
  tH1DPtr secondLastPT;

  double etaDiffAriadne;
  double NGlueAriadne;

  /** The PT distribution for charged particles in rapidity for final state. */
  tH1DPtr chargedPT;

  /** The number of charged particles in rapidity for final state. */
  tH1DPtr chargedN;

  /** The rap dist of charged particles in mid rapidity, following the ALICE paper. */
  tH1DPtr chargedNALICE;

  /** total charge multiplicity, following the ALICE paper. */
  tH1DPtr chargeMultALICE;

  /** counts the total weight of events with at least one charged track in |eta| < 1 */
  double INEL;

  /** Et density in mid rapidity. */
  Energy EtDensPHENIX;

  //keeps track on the total weight in some regions in impact parameter
  double totalWeight;
  double centralWeight;
  double midWeight;

  /** The distribution of maximum PT among the gluons before FSR. */
  tH1DPtr maxPTBefore;

  /** The sum of pT after DIPSY. */
  tH1DPtr DIPSYSumPT;

  /** The final sum of pT. */
  tH1DPtr finalSumET;


  /** distribution of number of gluon participants */
  tH1DPtr nParticipants;

  /** distribution of number of the eccentricity and area of the participating gluons */
  tH1DPtr RPEcc;
  tH1DPtr PartEcc;
  /** area in femtometer^2 */
  tH1DPtr OverlapArea;

  /** the average Part Ecc as function of the area in femtometers */
  tH1DPtr EccArea;

  /** the average area as function of participants */
  tH1DPtr AreaPart;

  /** the average of the square of the participant eccentricity */
  tH1DPtr PE2;

  /** the average of the 4:th power of the participant beccentricity */
  tH1DPtr PE4;

  /** The average of the square of eq 5 (event by event) binned in dNdeta */
  tH1DPtr v2flow;

  /** The average of the 4:th power of eq 5 (event by event) binned in dNdeta */
  tH1DPtr v4flow;

  /** A constant times the N participants, to simulate final state dNdeta */
  tH1DPtr dNchFlow;

  /** v2 squared, using eq 5 for high density, normal final state for low density. */
  tH1DPtr v2mix;

  /** the distribution of dN_ch/deta for mixed events. */
  tH1DPtr dNchMix;

  /** v2 squared, using final state. */
  tH1DPtr v2final;
  tH1DPtr v2finalWeights;

  /** 4 part corr. */
  tH1DPtr v2quad;
  tH1DPtr v2quadWeights;

  /** v2 with a symmetric gap 0.5 to 2. */
  tH1DPtr v2gap;

  /** the distribution of dN_ch/deta for final state. */
  tH1DPtr dNchFinal;

  /** the ratio of events being dense enough for hydro, binned in multiplicity. */
  tH1DPtr denseRatio;
  tH1DPtr denseRatioWeights;

  /** the average density dN/deta/Area, binned in multiplicity. */
  tH1DPtr density;

  /** the density distribution for Nch multiplicities 10 to 20. */
  tH1DPtr density1020;

  /** the density distribution for Nch multiplicities 40 to 60. */
  tH1DPtr density4060;

  /** the density distribution for Nch multiplicities above 60. */
  tH1DPtr density60plus;

  /** a testing 2D energy distribution. */
  tH2DPtr test2D;
  tH2DPtr test2DWeights;

  /** things to measure average number of participants */
  double averageNParticipants;


  //eps
  /** Number of gluons in mid rapidity as function of impact parameter B. */
  tH1DPtr histNGlue;
  tH1DPtr histNGlueSqr;
  tH1DPtr histNGlueWeights;

  /** Number of spectators as function of impact parameter B. */
  tH1DPtr histNSpectator;
  tH1DPtr histNSpectatorSqr;
  tH1DPtr histNSpectatorWeights;

  /** Number of gluons in mid rapidity as function of number of spectators. */
  tH1DPtr histNGlueSpect;
  tH1DPtr histNGlueSpectSqr;
  tH1DPtr histNGlueSpectWeights;

  /** The area of the central gluons as function of B */
  tH1DPtr histArea;
  tH1DPtr histAreaWeights;

  /** e1 and phi_1 as function of number of gluons in the central rapidity interval. */
  tH1DPtr histe1;
  tH1DPtr histe1Weights;
  tH1DPtr histe1N;
  tH1DPtr histe1NWeights;
  tH1DPtr histe1W;
  tH1DPtr histe1WWeights;
  tH1DPtr histphi1;
  tH1DPtr histphi1central;
  tH1DPtr histphi1mid;
  tH1DPtr histphi1N;
  tH1DPtr histphi1W;

  /** e2 and phi_2 as function of number of gluons in the central rapidity interval. */
  tH1DPtr histe2;
  tH1DPtr histe2Weights;
  tH1DPtr histe2N;
  tH1DPtr histe2NWeights;
  tH1DPtr histe2W;
  tH1DPtr histe2WWeights;
  tH1DPtr histphi2;
  tH1DPtr histphi2central;
  tH1DPtr histphi2mid;
  tH1DPtr histphi2N;
  tH1DPtr histphi2W;

  /** e3 and phi_3 as function of number of gluons in the central rapidity interval. */
  tH1DPtr histe3;
  tH1DPtr histe3Weights;
  tH1DPtr histe3N;
  tH1DPtr histe3NWeights;
  tH1DPtr histe3W;
  tH1DPtr histe3WWeights;
  tH1DPtr histphi3;
  tH1DPtr histphi3central;
  tH1DPtr histphi3mid;
  tH1DPtr histphi3N;
  tH1DPtr histphi3W;

  /** e4 and phi_4 as function of number of gluons in the central rapidity interval. */
  tH1DPtr histe4;
  tH1DPtr histe4Weights;
  tH1DPtr histe4N;
  tH1DPtr histe4NWeights;
  tH1DPtr histe4W;
  tH1DPtr histe4WWeights;
  tH1DPtr histphi4;
  tH1DPtr histphi4central;
  tH1DPtr histphi4mid;
  tH1DPtr histphi4N;
  tH1DPtr histphi4W;

  /** e5 and phi_5 as function of number of gluons in the central rapidity interval. */
  tH1DPtr histe5;
  tH1DPtr histe5Weights;
  tH1DPtr histe5N;
  tH1DPtr histe5NWeights;
  tH1DPtr histe5W;
  tH1DPtr histe5WWeights;
  tH1DPtr histphi5;
  tH1DPtr histphi5central;
  tH1DPtr histphi5mid;
  tH1DPtr histphi5N;
  tH1DPtr histphi5W;

  //with <e^2>
  /** e1 as function of number of gluons in the central rapidity interval. */
  tH1DPtr histe1Sqrd;
  tH1DPtr histe1SqrdWeights;
  /** e2 as function of number of gluons in the central rapidity interval. */
  tH1DPtr histe2Sqrd;
  tH1DPtr histe2SqrdWeights;
  /** e3 as function of number of gluons in the central rapidity interval. */
  tH1DPtr histe3Sqrd;
  tH1DPtr histe3SqrdWeights;
  /** e4 as function of number of gluons in the central rapidity interval. */
  tH1DPtr histe4Sqrd;
  tH1DPtr histe4SqrdWeights;
  /** e5 as function of number of gluons in the central rapidity interval. */
  tH1DPtr histe5Sqrd;
  tH1DPtr histe5SqrdWeights;

  //with <e^4>
  /** e1 as function of number of gluons in the central rapidity interval. */
  tH1DPtr histe1Quad;
  tH1DPtr histe1QuadWeights;
  /** e2 as function of number of gluons in the central rapidity interval. */
  tH1DPtr histe2Quad;
  tH1DPtr histe2QuadWeights;
  /** e3 as function of number of gluons in the central rapidity interval. */
  tH1DPtr histe3Quad;
  tH1DPtr histe3QuadWeights;
  /** e4 as function of number of gluons in the central rapidity interval. */
  tH1DPtr histe4Quad;
  tH1DPtr histe4QuadWeights;
  /** e5 as function of number of gluons in the central rapidity interval. */
  tH1DPtr histe5Quad;
  tH1DPtr histe5QuadWeights;

  //touched
  /** e2 as function of number of gluons in the central rapidity interval. */
  tH1DPtr histe2Touched;
  tH1DPtr histe2TouchedWeights;
  /** e3 as function of number of gluons in the central rapidity interval. */
  tH1DPtr histe3Touched;
  tH1DPtr histe3TouchedWeights;
  /** e4 as function of number of gluons in the central rapidity interval. */
  tH1DPtr histe4Touched;
  tH1DPtr histe4TouchedWeights;
  /** e5 as function of number of gluons in the central rapidity interval. */
  tH1DPtr histe5Touched;
  tH1DPtr histe5TouchedWeights;

  //eta in (-3,-1)
  /** e1 as function of number of gluons in the central rapidity interval. */
  tH1DPtr histe1m3tom1;
  tH1DPtr histe1m3tom1Weights;
  /** e2 as function of number of gluons in the central rapidity interval. */
  tH1DPtr histe2m3tom1;
  tH1DPtr histe2m3tom1Weights;
  /** e3 as function of number of gluons in the central rapidity interval. */
  tH1DPtr histe3m3tom1;
  tH1DPtr histe3m3tom1Weights;
  /** e4 as function of number of gluons in the central rapidity interval. */
  tH1DPtr histe4m3tom1;
  tH1DPtr histe4m3tom1Weights;
  /** e5 as function of number of gluons in the central rapidity interval. */
  tH1DPtr histe5m3tom1;
  tH1DPtr histe5m3tom1Weights;

  //eta in (1, 3)
  /** e1 as function of number of gluons in the central rapidity interval. */
  tH1DPtr histe1p1top3;
  tH1DPtr histe1p1top3Weights;
  /** e2 as function of number of gluons in the central rapidity interval. */
  tH1DPtr histe2p1top3;
  tH1DPtr histe2p1top3Weights;
  /** e3 as function of number of gluons in the central rapidity interval. */
  tH1DPtr histe3p1top3;
  tH1DPtr histe3p1top3Weights;
  /** e4 as function of number of gluons in the central rapidity interval. */
  tH1DPtr histe4p1top3;
  tH1DPtr histe4p1top3Weights;
  /** e5 as function of number of gluons in the central rapidity interval. */
  tH1DPtr histe5p1top3;
  tH1DPtr histe5p1top3Weights;

  /** frequency of the distance between the CoG of the gluons in (-3,-1) and (1,3)  */
  tH1DPtr histCoGdistance;

  /** correaltion between en in (-3,-1) and (1,3) for N_part > 100. */
  double over100Weight;
  double over100under250Weight;
  double over350Weight;
  tH2DPtr histeCorr1;
  tH2DPtr histeCorr2;
  tH2DPtr histeCorr3;
  tH2DPtr histeCorr4;
  tH2DPtr histeCorr5;
  tH1DPtr histphiCorr1;
  tH1DPtr histphiCorr2;
  tH1DPtr histphiCorr2mid;
  tH1DPtr histphiCorr2central;
  tH1DPtr histphiCorr3;
  tH1DPtr histphiCorr4;
  tH1DPtr histphiCorr5;

  //<eps(1,3)*eps(-3,-1)>
  double epFepB1, epF1, epB1, epSqrdF1, epSqrdB1, rho1;
  double epFepB2, epF2, epB2, epSqrdF2, epSqrdB2, rho2;
  double epFepB2mid, epF2mid, epB2mid, epSqrdF2mid, epSqrdB2mid, rho2mid;
  double epFepB2central, epF2central, epB2central, epSqrdF2central, epSqrdB2central, rho2central;
  double epFepB3, epF3, epB3, epSqrdF3, epSqrdB3, rho3;
  double epFepB4, epF4, epB4, epSqrdF4, epSqrdB4, rho4;
  tH1DPtr histe1e1p1top3;
  tH1DPtr histe1e1p1top3Weights;
  tH1DPtr histe2e2p1top3;
  tH1DPtr histe2e2p1top3Weights;
  tH1DPtr histe3e3p1top3;
  tH1DPtr histe3e3p1top3Weights;
  tH1DPtr histe4e4p1top3;
  tH1DPtr histe4e4p1top3Weights;
  tH1DPtr histe5e5p1top3;
  tH1DPtr histe5e5p1top3Weights;

  /** dN/deta in bins of number of spectators */
  tH2DPtr histdNdnSpectdeta;
  tH1DPtr histdNdeta0to100spect;
  tH1DPtr histdNdeta100to200spect;
  tH1DPtr histdNdeta200to300spect;
  tH1DPtr histdNdeta300to400spect;
  double spect0to100Weight;
  double spect100to200Weight;
  double spect200to300Weight;
  double spect300to400Weight;

  /** counts the average number of gluons in rapidity [-0.9,0.9]. */
  double centralGlueMult;

  /** flux tube correlation. */
  tH2DPtr histFluxTube;
  tH2DPtr histFluxTubeWeights;


private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  FSAnalysis & operator=(const FSAnalysis &);

};

}

#endif /* DIPSY_FSAnalysis_H */
