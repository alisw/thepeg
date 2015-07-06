// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the EndPoint class.
//

#include "EndPoint.h"

#ifdef PYTHIA7_TEMPLATES_IN_CC_FILE
// #include "EndPoint.tcc"
#endif

using namespace Pythia7;

EndPoint EndPoint::CC() const {
  EndPoint retEP(*this);
  retEP.setPt(-pTcomp());
  return retEP;
}

void EndPoint::UpdatedFrom (const EndPoint& currEP) {
  if(this != &currEP){
    theParticle = currEP.PData()->CC();
    theStrRg    = currEP.SR();
    thePTcomp       = - currEP.pTcomp();
    theGamma    = currEP.Gamma();
  }
}


// ************************ DB ******************************************

void EndPoint::echo() const {
  cout<<"\nEndPoint :"<<endl; 
  if(theParticle){
  cout<<"EP id() = "<<theParticle->id()<<endl;
  cout<<"\t EP m() = "<<mass()/GeV<<endl;
  }else{
    cout<<"Particle undefined"<<endl;
  }
  cout<<"\t EP Px() = "<<thePTcomp.first/GeV<<endl;
  cout<<"\t EP Py() = "<<thePTcomp.second/GeV<<endl;
  cout<<"\t EP Gamma() = "<<Gamma()/GeV2<<endl;
  if (theStrRg) {
    theStrRg->echo();
    cout<<"\t Pfwd() = "<<theStrRg->Pfwd()/GeV<<endl;
    cout<<"\t Pbwd() = "<<theStrRg->Pbwd()/GeV<<endl;
  }else cout<<"StringRegion :unknown"<<endl;     

}

