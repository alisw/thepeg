// -*- C++ -*-

#include "ThePEG/Config/algorithm.h"
#include "LundFlavourGenerator.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Utilities/VSelector.h"
#include "ThePEG/Interface/ClassDocumentation.h"


#ifdef PYTHIA7_TEMPLATES_IN_CC_FILE
// #include "LundFlavourGenerator.tcc"
#endif

using namespace Pythia7;


LundFlavourGenerator::LundFlavourGenerator()
  :theBaryonMod(2), DQsup(0.10), IS1DQsup(0.05), BaryonDecupletSup(1.0),
  extraEtaSup(1.0), extraEtapSup(0.4), MesonP(0.5), BssBsup(0.5), BMsBsup(0.5),
  sQsup(0.3), sDQvQsup(0.4), S1lightMesonP(0.5), S1sMesonP(0.6),
  S1hqMesonP(0.75), P_S0L1J1(0.0), P_S1L1J0(0.0), P_S1L1J1(0.0), P_S1L1J2(0.0),
  S1DQsup(0.0), POP_sDQvQsup(0.0), POP_S1DQsup(0.0),
  POP_delta(0.0), DQmaxWeight(0.0) {

  setDefaultMixingAnglesVec();
  setSU6weights();
}

LundFlavourGenerator::LundFlavourGenerator(const LundFlavourGenerator & x)
  :FlavourGenerator(x), theBaryonMod(x.theBaryonMod), DQsup(x.DQsup),
  IS1DQsup(x.IS1DQsup), BaryonDecupletSup(x.BaryonDecupletSup),
  extraEtaSup(x.extraEtaSup), extraEtapSup(x.extraEtapSup), MesonP(x.MesonP),
  BssBsup(x.BssBsup), BMsBsup(x.BMsBsup), sQsup(x.sQsup), sDQvQsup(x.sDQvQsup),
  S1lightMesonP(x.S1lightMesonP), S1sMesonP(x.S1sMesonP),
  S1hqMesonP(x.S1hqMesonP), P_S0L1J1(x.P_S0L1J1), P_S1L1J0(x.P_S1L1J0),
  P_S1L1J1(x.P_S1L1J1),  P_S1L1J2(x.P_S1L1J2),
  theMixingAngles(x.theMixingAngles), DQweight(x.DQweight), S1DQsup(x.S1DQsup),
  POP_sDQvQsup(x.POP_sDQvQsup), POP_S1DQsup(x.POP_S1DQsup),
  POP_delta(x.POP_delta), DQmaxWeight(x.DQmaxWeight),
  theSU6WeightsTable(x.theSU6WeightsTable),
  theMixingProbVec(x.theMixingProbVec){}
  //setDefaultMixingAnglesVec();
  //setSU6weights();


void LundFlavourGenerator::initialize(){
   setMesonFlavourMixingProbs();
   initGenPar();
}


tcPDPair LundFlavourGenerator::generateHadron(tcPDPtr quark) const {
  cPDPair ret;
  ret.first = generateHadron(quark, ret.second, 0);
  if ( ret.second ) ret.second = ret.second->CC();
  return ret;
}

PDPtr LundFlavourGenerator::
generateHadron(tcPDPtr inPD, cPDPtr & newPD, long curtainQ) const {
  PDPtr ret;

  if(!inPD) {
    newPD = PDPtr();
    return ret;
  }

  long inFl = inPD->id();

  do{
    extraSup=0;
    thePopRejection=0;

    // Select Meson/Baryon production
    if(isQuark(inFl) ) {
      (BaryonMod() && DQproduction())? Q2BaryonDQ(inFl): Q2MesonQ(inFl);
    }
    else if(isDiquark(inFl)){
      //the curtain Quark  exists only if BaryonMod >=2
      (curtainQ)? PopDQ2MesonDQ(inFl, curtainQ) : DQ2BaryonQ(inFl);
    }
  }while(extraSuppressed() );

  newPD = getParticleData(newFl);
  ret = getParticleData(theHadron);
  return ret;
}


tcPDPtr LundFlavourGenerator::getHadron(tcPDPtr inPD1, tcPDPtr inPD2) const {
  PDPtr ret;

  long inFl1 = inPD1->id();
  long inFl2 = inPD2->id();

  if(!consistentJoin(inFl1, inFl2)) return ret;

  if( isQuark(inFl1) && isQuark(inFl2) ) {
    theHadron = formMeson(inFl1, inFl2) ;
  }
  else {
    // Note:
    // In contrary to the Pythia6 version the suppression of final BMB
    // configuration is handled by standard Fragmentation Handler (LundFragHandler)
    // i.e if a inQ is joined with a inDQ having the  possibility
    // to give a PopCorn Meson + a new DQ then the full generation
    // is then rejected under the responsbility of the LundFragHandler.

    long inDQ = ( isDiquark(inFl1) )? inFl1 : inFl2;
    long inQ =  inFl1 + inFl2 - inDQ;

    theHadron = formBaryon(inQ, inDQ, getSU6MultWeightVec(inQ, inDQ) );
  }

  ret = getParticleData(theHadron);
  return ret;
}



tcPDPtr LundFlavourGenerator::
getBaryon(long iq1, long iq2, long iq3) const {
  if ( abs(iq1) >= 10 || abs(iq2) >= 10 || abs(iq3) ) return tcPDPtr();
  if ( iq1*iq2*iq3 == 0 ) return tcPDPtr();
  int sign = 0;
  if ( iq1 > 0 && iq2 > 0 && iq3 > 0 ) sign = 1;
  if ( iq1 < 0 && iq2 < 0 && iq3 < 0 ) sign = -1;
  if ( !sign ) return tcPDPtr();
  VSelector< pair<long,long> > sel;
  iq1 = abs(iq1);
  iq2 = abs(iq2);
  iq3 = abs(iq3);
  sel.insert(3.0, make_pair(iq1, 1000*max(iq2, iq3) + 100*min(iq2, iq3) + 3));
  if ( iq2 != iq3 )
    sel.insert(1.0, make_pair(iq1, 1000*max(iq2, iq3) + 100*min(iq2, iq3) + 1));
  sel.insert(3.0, make_pair(iq2, 1000*max(iq3, iq1) + 100*min(iq3, iq1) + 3));
  if ( iq3 != iq1 )
    sel.insert(1.0, make_pair(iq2, 1000*max(iq3, iq1) + 100*min(iq3, iq1) + 1));
  sel.insert(3.0, make_pair(iq3, 1000*max(iq1, iq2) + 100*min(iq1, iq2) + 3));
  if ( iq1 != iq2 )
    sel.insert(1.0, make_pair(iq3, 1000*max(iq1, iq2) + 100*min(iq1, iq2) + 1));
  pair<long,long> qdq = sel[rnd()];
  return getHadron(getParticleData(qdq.first), getParticleData(qdq.second));
}


// ================================================================
//                    Meson Production Methods
// ================================================================

long LundFlavourGenerator::selectQuark() const {
  return 1 + long((2.0 + sQsup)*rnd());
}

long LundFlavourGenerator::selectFlavour() const {
  throw std::runtime_error("Called unimplemented firtual function "
			   "'LundFlavourGenerator::selectFlavour() const'");
  return 1 + long((2.0 + sQsup)*rnd());
}

void LundFlavourGenerator::Q2MesonQ(long inQ) const {
  // Generate inQ + newQ -> newMeson
  newFl = sign(selectQuark(), -inQ);
  theHadron = formMeson(inQ, newFl);
}


long LundFlavourGenerator::formMeson(long inQ, long newQ) const {

  long Hfl = max(abs(inQ), abs(newQ) );
  long Lfl = min(abs(inQ), abs(newQ) );

  int multIdx = selectMesonMultiplet(Hfl);
  int Meson2Jp1 = 2*MesonSpin(multIdx)+1;

  long MesonId=0;

  if(Hfl != Lfl) {
    int Sign = (Hfl==abs(inQ))? sign(1,inQ) : sign(1,-inQ);
    Sign *= Hfl%2? -1: 1;
    MesonId = (100*Hfl + 10*Lfl + Meson2Jp1)*Sign;

  }
  else{
    if(Hfl <= 3)
      MesonId = 110*(RandomFlavourMixing(Hfl, multIdx, rnd()) ) + Meson2Jp1;
    else
      MesonId = 110*Hfl + Meson2Jp1;
  }

  //Reconstruct Meson of high multiplet
  if(multIdx==2 || multIdx==3) MesonId += sign(10000, MesonId);
  else if(multIdx==4) MesonId += sign(20000, MesonId);

  // Extra Suppression of Eta, Eta' Mesons
  //  if( (MesonId == 221 && rnd() > extraEtaSup ) ||
  //      (MesonId == 331 && rnd() > extraEtapSup) )  extraSup = 1;
  if( (MesonId == 221 && !rndbool(extraEtaSup)) ||
      (MesonId == 331 && !rndbool(extraEtapSup)) )  extraSup = 1;

   return MesonId;
}

int LundFlavourGenerator::selectMesonMultiplet(long Hfl) const {
  // Select at random the meson multiplet in function of its heavier quark (Hfl)

  //Select Meson Spin
  int S=0;
  if(Hfl<=2)
    //    S=int(S1lightMesonP + rnd());
    S = rndbool(S1lightMesonP)? 1: 0;
  else if(Hfl==3)
    //    S=int(S1sMesonP + rnd());
    S = rndbool(S1sMesonP)? 1: 0;
  else
    S = rndbool(S1hqMesonP)? 1: 0;

  //Select Meson Multiplet
  int multIdx=S;
  if(S==0 && P_S0L1J1>0.0) {
    //    if(rnd()< P_S0L1J1) multIdx=2;
    if( rndbool(P_S0L1J1) ) multIdx=2;
  }
  else if(S==1) {
    double S1L1MesonP = P_S1L1J0 + P_S1L1J1 + P_S1L1J2;
    if(S1L1MesonP>0.0){
      double rndmultIdx = rnd();
      if(rndmultIdx < P_S1L1J0 )
	multIdx=3;
      else if(rndmultIdx< P_S1L1J0+P_S1L1J1 )
	multIdx=4;
      else if(rndmultIdx< S1L1MesonP )
	multIdx=5;
    }
  }
  return multIdx;
}


// ================================================================
//                 Baryon Production Methods
// ================================================================

void LundFlavourGenerator::Q2BaryonDQ(long inQ) const {
  // Generate  inQ + newDQ -> newBaryon
  // create the newDQ and get the corresponding SU6 Weights Vector.
  // compute the Baryon Weight.
  // Loop until the Baryon Weigh is accepted then the inQ and the newDQ
  // and the weights are sent to the formBaryon() Method to get the Baryon.

  long DQq1, DQq2;
  int DQspin;
  WeightsVecPtr theWeights;
  double theBweight;

  do{
    do{
      // Create the new DiQuark
      double sDQsup = sQsup*sDQvQsup;
      DQq1 = 1+ int((2.0 + sDQsup)*rnd() );
      DQq2 = 1+ int((2.0 + sDQsup)*rnd() );
      DQspin = (DQq1>=DQq2)? 1: 0;
    }while(DQspinSup(DQspin));

    newFl = sign(1000*max(DQq1,DQq2)+ 100*min(DQq1,DQq2)+ (2*DQspin+1), inQ);

    //Get the SU6 Weights for the Baryon production
    theWeights = getSU6MultWeightVec(inQ, newFl);
    theBweight = OctetWT(theWeights) + BaryonDecupletSup*DecupletWT(theWeights);

    if (BaryonMod() >=2) {
      //Rescale of the Baryon weight for PopCorn DQ production
      double DQweight = PopDQweight(newFl);
      if(DQspin==0){
	DQweight /=(3.*POP_S1DQsup);
	theBweight *=(1.+DQweight)/(1.+ DQmaxWeight/(3.*POP_S1DQsup));
      }else if(DQspin==1)
	theBweight *=(1.+DQweight)/(1.+ DQmaxWeight);
    }
  }while( !rndbool(theBweight) );
  //  }while(rnd()>theBweight);

  theHadron = formBaryon(inQ, newFl, theWeights);
}



void LundFlavourGenerator::DQ2BaryonQ(long inDQ) const {   // To check
  // Generate inDQ + newQ -> newBaryon
  // Combine the inDQ with a newly created Quark to produce the final Baryon

  WeightsVecPtr theWeights;
  double theBweight;

  do{
    // create the new quark flavour
    newFl = sign(1 + int((2.0 + sQsup)*rnd()), inDQ);

    //prepare Baryon production
    theWeights = getSU6MultWeightVec(newFl, inDQ);
    theBweight = OctetWT(theWeights) + BaryonDecupletSup*DecupletWT(theWeights);
  }while( !rndbool(theBweight) );
  //  }while( theBweight < rnd() );

  theHadron = formBaryon(newFl, inDQ, theWeights);
}


LundFlavourGenerator::WeightsVecPtr
LundFlavourGenerator::getSU6MultWeightVec(long inQ, long inDQ) const{
  // Returns a Pointer to the Vector storing the relative probabilities
  // for the quark(inQ), diquark(inDQ) to join into a Baryon in the octet or the
  // decuplet multiplet

  inQ = abs(inQ);
  inDQ = abs(inDQ);

  int DQhq,DQlq, DQspin;
  getDQcontent(inDQ, DQhq, DQlq, DQspin);

  int Idx = 2*DQspin;
  if(DQhq == DQlq && DQspin==1) Idx=2*(DQspin+1);
  if(inQ!=DQhq && inQ!=DQlq) ++Idx;

  return(theSU6WeightsTable[Idx].begin());
}



long LundFlavourGenerator::
formBaryon(long inQ, long inDQ, WeightsVecPtr weights) const{
  // Get the Baryon Id taking the incoming Quark(inQ) Diquark(inDQ) flavours
  // and their corresponding SU6 Weights to form the Baryon

  long BaryonId=0;
  int absInQ = abs(inQ);
  int DQhq, DQlq, DQspin;
  getDQcontent(inDQ, DQhq, DQlq, DQspin);

  // get flavour content of the Baryon
  int fla = max(absInQ, DQhq);
  int flc = min(absInQ, DQlq);
  int flb = absInQ +  DQhq + DQlq - fla - flc;

  // get the relative probability for (Q-DQ) to join into a Baryon in the Octect
  double OctetProb = OctetWT(weights) /
    (OctetWT(weights) + BaryonDecupletSup*DecupletWT(weights) );
  //  int Baryon2Sp1 = (rnd() > OctetProb)? 4:2;
  int Baryon2Sp1 = ( !rndbool(OctetProb) )? 4:2;

  // Distinguish Lambda-Sigma like Baryons
  int LambdaLike=0;
  if(Baryon2Sp1==2 && fla>flb && flb>flc){
    if(DQspin==0) LambdaLike = (absInQ==fla)? 1 : int(0.25+rnd());
    if(DQspin==1 && absInQ!=fla ) LambdaLike = int(0.75+rnd());
  }

  BaryonId = (LambdaLike)?
    sign(1000*fla + 100*flc + 10*flb + Baryon2Sp1, inQ):
    sign(1000*fla + 100*flb + 10*flc + Baryon2Sp1, inQ);

  return(BaryonId);
}




// ================================================================
//                         PopCorn methods
// ================================================================


void LundFlavourGenerator::PopDQ2MesonDQ(long inDQ, long absCurtainQ ) const {
  // PopCorn Meson Generation : DQ -> Meson + DQ'
  // Join the free Flavour of the in Diquark with a new created Quark
  // to produce the Pop Corn Meson -> if possible
  // reconstruct the final Diquark
  //
  // WARN : abs(curtainQ) as input here

  long absDQhf = (abs(inDQ)/1000)%10;
  long absDQlf = (abs(inDQ)/100)%10;
  long popFl = sign(absDQhf + absDQlf - absCurtainQ, inDQ);

  //create new favour to form the PopCorn Meson
  long newQ =  sign( 1+int((2.0 + sQsup*POP_sDQvQsup*BMsBsup)*rnd()), -popFl) ;


  // Combine final DQ - Rejected if Spin1 suppressed
  long absNewQ = abs(newQ);
  thePopRejection = 0;

  if(PopConsistentDQ(absCurtainQ, absNewQ)){
    // create the new Flavour i.e the new Diquark
    int DQ2sp1 = (absCurtainQ != absNewQ) ?
      2*int(rnd() + 1.0/(1.0+POP_S1DQsup) ) + 1 : 3 ;

    newFl = sign(1000*max(absNewQ, absCurtainQ)
		 + 100*min(absNewQ, absCurtainQ) + DQ2sp1, -inDQ );

    // form the PopCorn Meson
    theHadron = formMeson(popFl, newQ);
  }
  else{
     thePopRejection = 1;
     newFl=0;
     theHadron=0;
  }
}

bool LundFlavourGenerator::PopConsistentDQ(long q1, long q2) const{
  // Return true if it is possible to join (q1, q2) into a Diquark

  if( q1*q2<0 ) return(0);

  q1 = abs(q1);
  q2 = abs(q2);

  double N = max(2.0, 1.+POP_S1DQsup);
//    if( (q1 != q2 &&   rnd() > (1.+POP_S1DQsup)/N) ||
//        (q1 == q2 && rnd() > 2.0/N) ) return(0);
  if( (q1 != q2 &&   rndbool((1.+POP_S1DQsup)/N) ) ||
      (q1 == q2 && rndbool(2.0/N) ) ) return(0);
  else
    return(1);
}


int LundFlavourGenerator::PopMesonN(long inDQ) const{
  // Return the number of PopCorn Mesons to be produced
  // in between of a Baryon antiBaryon pair according to the
  // PopCorn scheme.

  if(BaryonMod()<2) return(0);

  double DQweight= PopDQweight(inDQ);
  long DQ2sp1 = abs(inDQ)%10;

  if( DQ2sp1 == 1 ) DQweight /=(3.*POP_S1DQsup);

  return((1.0+DQweight)*rnd() >1.0 );
}



long LundFlavourGenerator::PopSelectCurtainFlavour(long inDQ) const{
  // Select a curtain quark of the PopCorn generation
  // Called by the standard Lund Fragmentation Handler to start/initialize
  // the PopCorn generation when it realizes that a Diquark is first produced
  // or if the string being fragmented contains a diquark endpoint.

  //inDQ=abs(inDQ);

  int hq = (inDQ/1000)%10;
  int lq = (inDQ/100)%10;

  int freeFl = hq + (lq - hq)*int(0.5 + rnd()) ;
  int curtainQ = hq + lq - freeFl;

  // Strange Quark suppresion in BBbar and BMBbar configuration
//    if( ( abs(freeFl) ==3  && rnd() > POP_delta) ||
//        ( abs(curtainQ) == 3 && rnd()<POP_delta ) ) {
  if( ( abs(freeFl) ==3  && !rndbool(POP_delta) ) ||
      ( abs(curtainQ) == 3 && rndbool(POP_delta) ) ) {

    freeFl = hq + lq - freeFl;
    curtainQ = hq + lq - curtainQ;
  }

  // WARNING :
  // Return the abs(curtainQ) because |curtainQ| is expected in
  // POPDQ2MesonDQ(...) to be changed

  return abs(curtainQ);
}

double  LundFlavourGenerator::PopDQweight(long inDQ) const{
  // Return the Diquark weight according to its flavour contain
  // for the PopCorn scheme.

  long absHq = (abs(inDQ)/1000)%10;
  long absLq = (abs(inDQ)/100)%10;

  if(absHq == 3) return(DQweight[1]);
  if(absLq == 3) return(DQweight[2]);
  return(DQweight[0]);
}


// ================================================================
//                      Joining Methodes
// ================================================================

bool LundFlavourGenerator::consistentJoin(long fl1, long fl2) const{
  // Returns false if the 2 incoming flavours are not consistent

  if( isQuark(fl1) && isQuark(fl2) && (fl1*fl2>0) ) return(0);
  if( isDiquark(fl1) && isDiquark(fl2) ) return(0);
  if( (isDiquark(fl1) || isDiquark(fl2)) && (fl1*fl2 <0) ) return(0);

  return(1);
}


// ================================================================
//              Non Interfaced Parameters initialization
// ================================================================

void LundFlavourGenerator::initGenPar(){
  // Compute parameters used in the flavour generation depending
  // on the baryon production Scheme.
  // To be called after any user changes of the relevant parameters


  S1DQsup = 3.0*IS1DQsup;                                     //PAR4

  if(BaryonMod()>=2){
    //set POPcorn parameters
    POP_sDQvQsup = sqrt(sDQvQsup);                             //PAR3M
    POP_S1DQsup = 1.0/(3.0*sqrt(IS1DQsup));                     //PAR4M
    POP_delta = BMsBsup/(BMsBsup + POP_sDQvQsup*BssBsup);      //PARDM

    setPopDQweights();

    S1DQsup *= (1.0+DQmaxWeight)/(1.0+DQmaxWeight/(3.0*POP_S1DQsup));  //PAR4
  }
}


void LundFlavourGenerator::setPopDQweights(){
  // Initialize the Diquark weights for the PopCorn Scheme
  // Called by initGenPar() after any user changes of the relevant
  // interfaced parameters

  if( !DQweight.size() ) DQweight = vector<double>(3, 0.0);

  DQweight[0] = MesonP*
    (2.0 + (1.0 + sQsup*POP_sDQvQsup*BMsBsup)*(1.0 + POP_S1DQsup) );

  DQweight[1] = BMsBsup*DQweight[0]/(2.0*POP_sDQvQsup) +
    MesonP*BssBsup*((1.0+POP_S1DQsup) + sQsup*POP_sDQvQsup*BMsBsup);

  DQweight[2] = 2.0*MesonP*BssBsup*BMsBsup*
    (sQsup*BMsBsup + (1.+ POP_S1DQsup)/POP_sDQvQsup);

  DQmaxWeight= max(max(DQweight[0],DQweight[1]), DQweight[2]); //PARSM

}


void  LundFlavourGenerator::setDefaultMixingAnglesVec(){
  // Set the Default mixing angles for the 6 meson multiplets
  // considered in the standard Lund Flavour Generator

  int nMultiplets=6;
  theMixingAngles.push_back(Constants::pi/4.);
  for(int i=1; i!=nMultiplets; ++i) theMixingAngles.push_back(0.0);
}



void LundFlavourGenerator::setMesonFlavourMixingProbs(){
  // Given the set of Mixing angles for the different meson multiplets,
  // defaults (or user defined) initializes the Mixing probalities for the
  // production of diagonal flavour Mesons.
  // To be called after any user changes of the mixing angles

  int nMultiplets= theMixingAngles.size();
  int nFlavours=3;

  // Wave Function Coefficients - should not be touched
  vector<double> alpha(nFlavours), beta(nFlavours);
  alpha[0]=0.0;  alpha[1]=0.0;  alpha[2]=1.0;
  beta[0]=0.5;   beta[1]=0.5;   beta[2]=0.0;

  theMixingProbVec.clear();
  theMixingProbVec.resize(nMultiplets,
			  vector<double>((nFlavours-1)*nFlavours, 0.0));

  for(int multIdx=0; multIdx<nMultiplets; ++multIdx){
    for(int Id=0; Id<nFlavours; ++Id){
      int flProbIdx = 2*Id+1;

      theMixingProbVec[multIdx][flProbIdx] =
	alpha[Id]*sqr(cos(theMixingAngles[multIdx]) ) +
	beta[Id]*sqr(sin(theMixingAngles[multIdx]) );

      theMixingProbVec[multIdx][flProbIdx-1] =
	alpha[Id] + beta[Id];
    }
  }
}

void LundFlavourGenerator::setSU6weights(){
  // Initialize non Interfaced SU6 Weights Table
  // Values taken from reference : NPB197(1982)45
  // Should not be touched by the user.

  //ud0 DQ
  WeightsVec ud0cmFl(2), ud0ncmFl(2) ;
  ud0cmFl[0] = 3./4.  ; ud0ncmFl[0] = 1.0/2.0 ;
  ud0cmFl[1] = 0.0    ; ud0ncmFl[1] = 0.0   ;

  //ud1 DQ
  WeightsVec ud1cmFl(2), ud1ncmFl(2) ;
  ud1cmFl[0] = 1.0/12.0   ; ud1ncmFl[0] = 1.0/6.0 ;
  ud1cmFl[1] = 2.0/3.0    ; ud1ncmFl[1] = 1.0/3.0 ;

  //uu1 DQ
  WeightsVec uu1cmFl(2), uu1ncmFl(2) ;
  uu1cmFl[0] = 0.0   ; uu1ncmFl[0]= 1.0/6.0 ;
  uu1cmFl[1] = 1.0   ; uu1ncmFl[1]= 1.0/3.0 ;


  theSU6WeightsTable = WeightsTable(6);

  theSU6WeightsTable[0] = ud0cmFl;
  theSU6WeightsTable[1] = ud0ncmFl;

  theSU6WeightsTable[2] = ud1cmFl;
  theSU6WeightsTable[3]=ud1ncmFl;

  theSU6WeightsTable[4]=uu1cmFl;
  theSU6WeightsTable[5]=uu1ncmFl;
}


// ================================================================
//                    Standard Interface Methods
// ================================================================

ClassDescription<LundFlavourGenerator>
LundFlavourGenerator::initLundFlavourGenerator;

void LundFlavourGenerator::persistentOutput(PersistentOStream & os) const {
  os << theBaryonMod << DQsup << IS1DQsup <<  BaryonDecupletSup
     << extraEtaSup << extraEtapSup << MesonP << BssBsup << BMsBsup
     << sQsup << sDQvQsup << S1lightMesonP << S1sMesonP << S1hqMesonP
     << P_S0L1J1 << P_S1L1J0 << P_S1L1J1 << P_S1L1J2
     << theMixingAngles << DQweight << S1DQsup << POP_sDQvQsup << POP_S1DQsup
     << POP_delta << DQmaxWeight << theSU6WeightsTable << theMixingProbVec;
}

void LundFlavourGenerator::persistentInput(PersistentIStream & is, int ) {
  is >> theBaryonMod >> DQsup >> IS1DQsup >>  BaryonDecupletSup
     >> extraEtaSup >> extraEtapSup >> MesonP >> BssBsup >> BMsBsup
     >> sQsup >> sDQvQsup >> S1lightMesonP >> S1sMesonP >> S1hqMesonP
     >> P_S0L1J1 >> P_S1L1J0 >> P_S1L1J1 >> P_S1L1J2
     >> theMixingAngles >> DQweight >> S1DQsup >> POP_sDQvQsup >> POP_S1DQsup
     >> POP_delta >>  DQmaxWeight >> theSU6WeightsTable >> theMixingProbVec;
}

void LundFlavourGenerator::Init() {

  static ClassDocumentation<LundFlavourGenerator> documentation
    ("The class responsible for the flavour generation according to the "
     "Lund scheme of fragmentation.");

  static Switch<LundFlavourGenerator, int> interfaceBaryonMod
    ("BaryonMode",
     "Choice of the baryon production mode",
     &LundFlavourGenerator::theBaryonMod, 2, false, false);

  static SwitchOption interfaceBaryonMod0
    (interfaceBaryonMod,
     "NoBaryon",
     "no baryon-antibaryon pair production at all", 0);
  static SwitchOption interfaceBaryonMod1
    (interfaceBaryonMod,
     "Baryon-antiBaryon_Mod",
     "diquark-antidiquark pair production allowed", 1);
  static SwitchOption interfaceBaryonMod2
    (interfaceBaryonMod,
     "PopCorn",
     "diquark-antidiquark pair production allowed with possibility for diquark"
     "\nto be split according to the `popcorn' scheme.", 2);

  static Parameter<LundFlavourGenerator, double> interface_DQsup
    ("DQsup",
     "Is the ratio P(qq)/P(q) : the suppression of diquark-antidiquark pair"
     "\nproduction in the colour field, compared with quark--antiquark production"
     "[no Unit]",
     &LundFlavourGenerator::DQsup, 0.10, 0.0, 1.0, false, false, true);

  static Parameter<LundFlavourGenerator, double> interface_IS1DQsup
    ("S1DQsup",
     "Is (1/3)P(ud_1)/P(ud_0) : the suppression of spin 1 diquarks compared with"
     "\n spin 0 ones (excluding the factor 3 coming from spin counting)."
     "[no Unit]",
     &LundFlavourGenerator::IS1DQsup, 0.05, 0.0, 1.0, false, false, true);

  static Parameter<LundFlavourGenerator, double> interface_BaryonDecupletSup
    ("SBaryonDecupletSup",
     "Extra suppression factor multiplying the ordinary SU(6) weight for spin 3/2"
     "\nbaryons and hence a means to break SU(6) in addition to the dynamic breaking"
     "\nimplied by the parameters sQsup, sDQvQsup, S1DQsup, BssBsup, BMsBsup"
     "[no Unit]",
     &LundFlavourGenerator::BaryonDecupletSup, 1.0, 0.0, 5.0, false, false, true);

  static Parameter<LundFlavourGenerator, double> interface_extraEtaSup
    ("extraEtaSup",
     "Extra suppression factor for Eta production in fragmentation. If an Eta is"
     "\n rejected a new flavour pair is generated and a new hadron formed."
     "[no Unit]",
     &LundFlavourGenerator::extraEtaSup, 1.0, 0.0, 5.0, false, false, true);

  static Parameter<LundFlavourGenerator, double> interface_extraEtapSup
    ("extraEtapSup",
     "Extra suppression factor for Eta' production in fragmentation. If an Eta' is"
     "\n rejected a new flavour pair is generated and a new hadron formed."
     "[no Unit]",
     &LundFlavourGenerator::extraEtapSup, 0.4, 0.0, 1.0, false, false, true);

  static Parameter<LundFlavourGenerator, double> interface_MesonP
    ("MesonP",
     "Determine the relative occurence of baryon production by [B M Bbar] and by"
     "\n[B Bbar] configurations in the popcorn baryon production model. [no Unit]",
     &LundFlavourGenerator::MesonP, 0.5, 0.0, 1.0, false, false, true);

  static Parameter<LundFlavourGenerator, double> interface_BssBsup
    ("BssBsup",
     "Extra suppression for having a [ss_bar] pair shared by the B and Bbar"
     "\nof a [B M Bbar] configuration. [no Unit]",
     &LundFlavourGenerator::BssBsup, 0.5, 0.0, 1.0, false, false, true);

  static Parameter<LundFlavourGenerator, double> interface_BMsBsup
    ("BMsBsup",
     "Extra suppression for having a strange meson M in a [B M Bbar] configuration."
     "[no Unit]",
     &LundFlavourGenerator::BMsBsup, 0.5, 0.0, 1.0, false, false, true);

  static Parameter<LundFlavourGenerator, double> interface_sQsup
    ("sQsup",
     "Is the ratio P(s)/P(u), the suppression of s quark pair production"
     "\nin the field compared with  u(d) pair production [noUnit].",
     &LundFlavourGenerator::sQsup, 0.3, 0.0, 1.0, false, false, true);

  static Parameter<LundFlavourGenerator, double> interface_sDQvQsup
    ("sDQvQsup",
     "Is the ratio (P(us)/P(ud))/(P(s)/P(d)), the extra suppression of strange"
     "diquark production  compared with the normal suppression of strange quarks "
     "[noUnit]",
     &LundFlavourGenerator::sDQvQsup, 0.4, 0.0, 1.0, false, false, true);

  static Parameter<LundFlavourGenerator, double> interface_S1lightMesonP
    ("S1lightMesonP",
     "The probability that a light meson (containing u and d quarks only) has"
     "spin 1 \n(with (1-S1lightMeson) the probability for spin 0) when formed"
     "formed in thefragmentation [noUnit].",
     &LundFlavourGenerator::S1lightMesonP, 0.5, 0.0, 1.0, false, false, true);

  static Parameter<LundFlavourGenerator, double> interface_S1sMesonP
     ("S1sMesonP",
      "The probability that a strange meson has spin 1 [noUnit].",
      &LundFlavourGenerator::S1sMesonP, 0.6, 0.0, 1.0, false, false, true);

  static Parameter<LundFlavourGenerator, double> interface_S1hqMesonP
     ("S1hqMesonP",
      "The probability that a charm or heavier meson has spin 1 [noUnit].",
      &LundFlavourGenerator::S1hqMesonP, 0.75, 0.0, 1.0, false, false, true);

  static Parameter<LundFlavourGenerator, double> interface_P_S0L1J1
    ("P_S0L1J1",
     "The probability that a spin = 0 meson is produced with an orbital angular"
     "\nmomentum 1, for a total spin = 1 [noUnit].",
     &LundFlavourGenerator::P_S0L1J1, 0.0, 0.0, 1.0, false, false, true);

  static Parameter<LundFlavourGenerator, double> interface_P_S1L1J0
    ("P_S1L1J0",
     "The probability that a spin = 1 meson is produced with an orbital angular"
     "\n momentum 1, for a total spin = 0 [noUnit].",
     &LundFlavourGenerator::P_S1L1J0, 0.0, 0.0, 1.0, false, false, true);

  static Parameter<LundFlavourGenerator, double> interface_P_S1L1J1
    ("P_S1L1J1",
     "The probability that a spin = 1 meson is produced with an orbital angular"
     "\nmomentum 1, for a total spin = 1  [noUnit].",
     &LundFlavourGenerator::P_S1L1J1, 0.0, 0.0, 1.0, false, false, true);

  static Parameter<LundFlavourGenerator, double> interface_P_S1L1J2
    ("P_S1L1J2",
     "The probability that a spin = 1 meson is produced with an orbital angular"
     "\nmomentum 1, for a total spin = 2  [noUnit].",
     &LundFlavourGenerator::P_S1L1J2, 0.0, 0.0, 1.0, false, false, true);

  interfaceBaryonMod.rank(10);
  interface_DQsup.rank(9);
  interface_IS1DQsup.rank(8);
  interface_MesonP.rank(7);
  interface_S1lightMesonP.rank(6);
  interface_S1sMesonP.rank(5);
  interface_S1hqMesonP.rank(4);

}


// *** DEBUG ***
void LundFlavourGenerator::DBprint(){
  cout<<"theBaryonMod= "<<theBaryonMod
      <<", DQsup="<<DQsup
      <<", IS1DQsup= "<<IS1DQsup
      <<", \nBaryonDecupletSup= "<<BaryonDecupletSup
      <<", extraEtaSup= "<<extraEtaSup
      <<", extraEtapSup= "<<extraEtapSup
      <<", \nMesonP= "<<MesonP
      <<", BssBsup= "<<BssBsup
      <<", BMsBsup= "<<BMsBsup
      <<", \nsQsup= "<<sQsup
      <<", sDQvQsup= "<<sDQvQsup
      <<", S1lightMesonP= "<<S1lightMesonP
      <<", \nS1sMesonP= "<<S1sMesonP
      <<", S1hqMesonP= "<<S1hqMesonP
      <<", P_S0L1J1= "<<P_S0L1J1
      <<", \nP_S1L1J0= "<<P_S1L1J0
      <<", P_S1L1J1= "<<P_S1L1J1
      <<", P_S1L1J2= "<<P_S1L1J2
      <<endl;

  cout<<"Mixing Probabilities :"<<endl;
  for(unsigned multIdx=0; multIdx<theMixingProbVec.size(); ++multIdx){
    for(unsigned flIdx=0; flIdx<theMixingProbVec[0].size(); ++flIdx){
      cout<<"mixingProb("<<multIdx<<","<<flIdx<<") = "<<
	theMixingProbVec[multIdx][flIdx]<<" ";
    }
    cout<<endl;
  }

  cout<<endl;
}


