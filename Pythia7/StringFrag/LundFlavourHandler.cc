// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the LundFlavourHandler class.
//

#include "LundFlavourHandler.h"
#include "ThePEG/PDT/StandardMatchers.h"

#ifdef PYTHIA7_TEMPLATES_IN_CC_FILE
// #include "LundFlavourHandler.tcc"
#endif

using namespace Pythia7;

PDPtr LundFlavourHandler::generateHadron(tcPDPtr inPD, cPDPtr& newPD){

  //Leading Baryon, Meson  production
  if(DiquarkMatcher::Check(inPD->id()) && PopN() <0 && !curtainQId() ) {
    PopN() = theFlGen->PopMesonN(inPD->id() );
  }

  PDPtr newH;

  if( PopN() > 0 && !curtainQId() ){
    do{
      curtainQId() = theFlGen->PopSelectCurtainFlavour(inPD->id() );
      newH = theFlGen->generateHadron(inPD, newPD, curtainQId() );
    }while( theFlGen->PopGenRejected()  );
  }
  else {
    newH = theFlGen->generateHadron(inPD, newPD, curtainQId() );

    if( DiquarkMatcher::Check(newPD->id())  &&  PopN() <=0 )
      PopN() = theFlGen->PopMesonN(newPD->id() );
  }

  if( curtainQId() && !(--PopN()) ) curtainQId() = 0;
 
  return newH;
} 

