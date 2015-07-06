// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the StringRegion class.
//


#include "StringRegion.h"
#include <algorithm>

#ifdef PYTHIA7_TEMPLATES_IN_CC_FILE
// #include "StringRegion.tcc"
#endif

using namespace Pythia7;

#include "ThePEG/Repository/CurrentGenerator.h"

void StringRegion::init() {

  //instead of (p1+p2)^2
  w2=2.0*(Pplus()*Pminus());

  if ( w2 <= 0.0*GeV2 ) return;

  LorentzMomentum p1 = Pplus(), p2=Pminus(); 
  ThreeVector<double> a = p1.vect()*(1.0/p1.e()) - p2.vect()*(1.0/p2.e());
  
  LorentzVector<double> p3, p4;

  if( abs(a.x()) > max(abs(a.y()), abs(a.z()) )) {               
    p3.setZ(1.0);
    p4.setY(1.0);
  }
  else {
    if( abs(a.y()) > abs(a.z())){
      p3.setX(1.0);
      p4.setZ(1.0);
    } 
    else {
      p3.setX(1.0);
      p4.setY(1.0);
    } 
  }
  
  InvEnergy cx1 = (p3*p1)/w2;    
  InvEnergy cx2 = (p3*p2)/w2;
  double cxx = 1.0/sqrt(1.0 + 4.0*cx1*cx2*w2);
      
  InvEnergy cy1 = (p4*p1)/w2;    
  InvEnergy cy2 = (p4*p2)/w2; 
  double cyx = 2.0*cxx*(cx1*cy2 + cx2*cy1)*w2;
  double cyy = 1.0/sqrt(1.0 + 4.0*cy1*cy2*w2 - cyx*cyx);
      
  Ex = cxx*(p3 - 2.0*cx2*p1 - 2.0*cx1*p2 );
  Ey = cyy*(p4 - 2.0*cy2*p1 - 2.0*cy1*p2 - cyx*Ex);
  
}
  


// ************  Only for DBG ********************

void StringRegion::echo() const {
  if(j()==0 && k()==0 ) cout<<"Unknown Str Region !!\n"<<endl;
  cout<<"SR("<<j()<<","<<k()<<")"<<endl;
  cout<<"Oriented SR("<<Ifwd()<<","<<Ibwd()<<")"<<endl;
  cout<<"\t P+{"<<j()<<"} = "<<ounit(Pplus(), GeV);
  cout<<"\t P-{"<<k()<<"} = "<<ounit(Pminus(), GeV)<<endl; 
  cout<<"\t ex ="<<ex()<<"\t ey  ="<<ey()<<endl;
  cout<<"\t W2(SR) = "<<W2()/GeV2<<endl; 
  cout<<"\t --> The SR orientation : "<<
  ((Oriented::Dir()==Oriented::right)? "right":"left")<<endl;
}








































