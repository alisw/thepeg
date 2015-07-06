// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LundZGenerator class.
//

#include <algorithm>
#include "LundZGenerator.h"
#include "ThePEG/PDT/StandardMatchers.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Repository/EventGenerator.h"


#ifdef PYTHIA7_TEMPLATES_IN_CC_FILE
// #include "LundZGenerator.tcc"
#endif

using namespace Pythia7;


LundZGenerator::~LundZGenerator() {}


void LundZGenerator::setShape(cPDPtr lastPD, cPDPtr newPD, Energy2 mT2 ) const {
  bmT2 = bSym()*mT2;

  theAf = aSym();
  theCf = 1.0;
  
  // set Parameter for DQ fragmentation
  theAf = ( DiquarkMatcher::Check(*newPD) )? 
    aSym() + deltaDQ :
    aSym() ;

  theCf = 1.0;
  if( DiquarkMatcher::Check(*lastPD) ) theCf -= deltaDQ;
  if( DiquarkMatcher::Check(*newPD) ) theCf += deltaDQ;

  // set Heavy Flavour shape
  long hf = ( DiquarkMatcher::Check(*lastPD) ) ? 
    ( (lastPD->id()/1000)%10 ) :
    ( lastPD->id() ) ;

  if( hf > 3 ) 
    theCf += rQ(hf)* bSym() * sqr( lastPD->constituentMass() );  

}


double LundZGenerator::generate(cPDPtr lastPD, cPDPtr newPD, Energy2 mT2 ) const{
  //  Timer<32> timer("LundZGenerator::generate(...)");

  // --- setParameters and flags ---
  setShape(lastPD, newPD, mT2);
  // set flags / init generation 
  bool theSelectedZValue = false;
  double Z = 0.0;

  // --- getZmax(); --- 
  double Zmax;
  if(af() < 0.02){ 
    Zmax =1.0;
    if(cf() > bmT2) Zmax=bmT2/cf(); 
  }else{ 
    if( cCloseToA() ) {
      Zmax = bmT2/(bmT2 +cf());
    }
    else{  
      double Delta = sqrt( sqr(bmT2 - cf()) + 4.0*af()*bmT2 );
      Zmax = 0.5*(bmT2 + cf() - Delta)/(cf() - af());
      if(Zmax > 0.9999 && bmT2 > 100.) Zmax = min(Zmax, (1.- af()/bmT2) );
    }
  }


  do{
    // --- pickZvalue(); --- 
    Z = rnd();
    double preValue = 1.0;

    if( smallZmax(Zmax) ) {  
      // --- weightedLowZvalue(); ---
      //Subdivide Z range peaked near Zero
      double ZdivC = 0.0;
      double Fint = 0.0;
      double Zdiv = 2.75 * Zmax;

      if( cCloseToOne() ) {    
	Fint = 1.0 - log(Zdiv);
      }else{
	double deltaC = 1.0 - cf();
	ZdivC = pow(Zdiv, deltaC );
	Fint = (cf()*ZdivC - 1.0)/(ZdivC*deltaC);
      }
      
      //Choice of Z preweighted for peak at low z
      //      if( Fint*rnd() <= 1.0 ){ 
      if( Fint <= 1.0 || rndbool(1.0/Fint) ){ 
	Z *= Zdiv;
      }else { 
      	if( cCloseToOne() ) {  
	  Z = pow(Zdiv, Z);
	  preValue = Zdiv/Z;
	}else {
	  double xpn = 1.0/(1.0 - cf());
	  Z = pow(ZdivC + Z*(1.0 - ZdivC) , xpn);
	  preValue = pow(Zdiv/Z, cf());
	} 
      } 
    }
    
    if( largeZmax(Zmax) ){  
      // --- weightedHightZvalue(); ---
      //Subdivide Z range peaked near one
      double cOb = cf()/bmT2;
      double Fscb = sqrt(4.0 + sqr(cOb));
      double Zdiv = Fscb - 1.0/Zmax - (cOb)*log(0.5*Zmax*(Fscb + cOb));
      
      if(af() >= 0.02) Zdiv += (af()/bmT2)*log(1.0 - Zmax);  
      
      Zdiv = min(Zmax, max(0.0, Zdiv));
      double Fint = 1.0 + bmT2*(1.0 - Zdiv);
      
      //Choice of Z preweighted for peak at high z
      //      if( Fint*rnd() <= 1.0 ) {
      if ( Fint <= 1.0 || rndbool(1.0/Fint) ) {
	Z = Zdiv + log(Z)/bmT2;
	preValue = exp(bmT2*(Z-Zdiv));
      }else {
	Z = Zdiv + Z*(1.0 - Zdiv); 
      }    
    }
    

    // --- Check Z value found -> restart if necessary ---
    
    if(Z <= 0.0 || Z >= 1.0){
      theSelectedZValue = false;
    }
    else {
      double Fexp = cf()*log(Zmax/Z) + bmT2*(1.0/Zmax - 1.0/Z);
      if(af() >= 0.02) 
	Fexp += af()*log((1.0 - Z)/(1.0 - Zmax));
      double Val = exp( max(-50.0, min(50.0, Fexp)));
      //      if(Val < preValue*rnd()) {
      if ( preValue > Val && !rndbool(Val/preValue) ) {
	theSelectedZValue = false;
      }else {
	theSelectedZValue = true;
      } 
    }
  }while(theSelectedZValue == false);
   
  
  return Z;
}




// *******************************************************************
//                Standard Functions for Interface 
// *******************************************************************

void LundZGenerator::persistentOutput(PersistentOStream & os) const {
  os << asym << ounit(bsym, 1.0/GeV2) << deltaDQ << rQc;  
}

void LundZGenerator::persistentInput(PersistentIStream & is, int) { 
  is >> asym >> iunit(bsym, 1.0/GeV2) >> deltaDQ >> rQc;
}


ClassDescription<LundZGenerator> initLundZGenerator;

void LundZGenerator::Init() {

  static ClassDocumentation<LundZGenerator> documentation
    ("Implements the Lund symmetric fragmentation function.");

  static Parameter<LundZGenerator, double> interfaceAsym
    ("a",
     "The a parameter of the Lund symmetric fragmentation function [noUnit].",
     &LundZGenerator::asym, 0.3, 0.0, 1.0, false, false, true);

  static Parameter<LundZGenerator, InvEnergy2> interfaceBsym
    ("b",
     "The b parameter of the Lund symmetric fragmentation function [GeV^-2].",
     &LundZGenerator::bsym, 1.0/GeV2, 0.58/GeV2, 0.0/GeV2, 1.0/GeV2,
     false, false, true);

  static Parameter<LundZGenerator, double> interfaceDeltaDQ
    ("deltaDQ",
     "Describe the amount by which the effective \\f$a\\f$ parameter in "
     "the Lund flavour dependent symmetric fragmentation function is "
     "assumed to be larger than the normal \\f$a\\f$ when diquarks are "
     "produced. [no Unit]",
     &LundZGenerator::deltaDQ, 0.5, 0.0, 1.0, false, false, true);

  static Parameter<LundZGenerator, double> interfaceRQc
    ("rQc",
     "\\f$r_Q\\f$ factor in the Bowler modification of the Lund symmetric "
     "fragmentation function for heavy endpoint quark. It is assumed "
     "\\f$r_Q\\f$ is the same for all"
     "flavour heavier than the strange quark. [no Unit]", 
     &LundZGenerator::rQc, 1.0, 0.0, 1.0, false, false, true);

  interfaceAsym.rank(10);
  interfaceBsym.rank(9);

}
