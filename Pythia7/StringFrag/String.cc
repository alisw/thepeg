// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the HandlerBase class.
//

#include "String.h"
#include "StringRegion.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Utilities/UtilityBase.h"

#ifdef PYTHIA7_TEMPLATES_IN_CC_FILE
// #include "String.tcc"
#endif

using namespace Pythia7;
using namespace Pythia7::ParticleID;


String::String(const String & sr): 
  Ptotrem(sr.Ptotrem), pplus(sr.pplus), pminus(sr.pminus),
  xPlusRem(sr.xPlusRem), xMinusRem(sr.xMinusRem),
  theStringRegionMap(sr.theStringRegionMap), ns(sr.ns),
  Ptot(sr.Ptot), partons(sr.partons) {}
				
String::~String(){  
  if(theStringRegionMap.size() != 0) clearStringRegionMap();
}

String::String(const PVector & partins)
  : ns(0) {
  partons = partins;
  Ptot = Utilities::sumMomentum(partons);
  ns = partons.size() - 1;
  if ( partons[0]->id() == ParticleID::g ) ++ns;
  pplus  = MomentumVector(ns);
  pminus = MomentumVector(ns);

  for(int is = 0; is < ns; is++) {
    int i = is;
    int j = (i + 1)%partons.size();

    LorentzMomentum p1 = (partons[i]->id() == ParticleID::g? 0.5: 1.0)*
      partons[i]->momentum(); 

    LorentzMomentum p2 = (partons[j]->id() == ParticleID::g? 0.5: 1.0)*
      partons[j]->momentum(); 

    Energy2 Dn = sqrt(sqr(p1*p2) - p1.m2()*p2.m2() ); 
    double k1 = 0.5*((p2*p2 + p1*p2)/Dn - 1.0); 
    double k2 = 0.5*((p1*p1 + p1*p2)/Dn - 1.0); 
    
    pplus[is]  = (1.0 + k1)*p1 -  k2*p2 ;
    pminus[is] = (1.0 + k2)*p2 -  k1*p1 ;  
  }

  // Set Xrem 
  InitXrem();
}

// **************************************************
//           String Region Management  
// ***************************************************

cStringRegionPtr String::getStringRegionPtr(int j, int k) {

  if((j<1 || j>ns) || (k<1 || k>ns) || j > k ) throw Veto(); 

  StringRegionIndex idx = make_pair(j,k);

  StringRegionMap::iterator rit =
    theStringRegionMap.find(idx);

  if(rit != theStringRegionMap.end()){
    return rit->second;
  }
  else{
    return ( theStringRegionMap[idx] = new StringRegion(j, k, this) );
  } 
}


void String::clearStringRegionMap(){
  for(StrRgMapIt it = theStringRegionMap.begin(); 
      it!=  theStringRegionMap.end(); ++it){
    delete it->second;
  }
  theStringRegionMap.clear();
  return;
}


cStringRegionPtr String::nextUp(cStringRegionPtr first){
  return (Oriented::Dir() == Oriented::right)?
    getStringRegionPtr(first->j(), first->k()+1 ):
    getStringRegionPtr(first->j()-1, first->k() ); 
}


cStringRegionPtr String::nextDown(cStringRegionPtr first){
  return (Oriented::Dir()==Oriented::right)?
    getStringRegionPtr(first->j()+1, first->k()):
    getStringRegionPtr(first->j(), first->k()-1 );   
}


// **********************************************************
//                           Dgb Only
// **********************************************************

void String::echo() const{
  //Output parameters for test
  cout<<endl<<"The String Regions : "
      <<"nb. of primary string regions = "<<ns<<endl;
  for (int i=0; i<ns; i++){
    cout<<"\t P+("<<i+1<<")= "<<ounit(Pplus(i+1), GeV);
    cout<<"\t P-("<<i+1<<")= "<<ounit(Pminus(i+1), GeV)<<"\n";
  }
  cout<<endl;
  cout<<"P_{tot rem} [Unit?] = "<<ounit(PtotRem(), GeV)<<"\n";    
  cout<<"W^2_Rem [Unit?] = "<< Wrem2()/GeV2 <<"\n"<<endl;
}
















