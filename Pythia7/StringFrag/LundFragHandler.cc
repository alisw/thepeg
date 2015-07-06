// -*- C++ -*-
#include "LundFragHandler.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "ThePEG/Handlers/Hint.h"
#include "ThePEG/PDT/ParticleData.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/Step.h"
#include "ThePEG/PDT/PDT.h"   
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Config/algorithm.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Handlers/ClusterCollapser.h"
#include "ThePEG/Utilities/UtilityBase.h"   
#include "ThePEG/Utilities/SimplePhaseSpace.h"   
#include "ThePEG/PDT/RemnantDecayer.h"   


#ifdef PYTHIA7_TEMPLATES_IN_CC_FILE
// #include "LundFragHandler.tcc"
#endif

using namespace Pythia7;
using std::max;
using std::abs;


LundFragHandler::LundFragHandler()
  : pWmin0(0.8*GeV), pK(2.), pDelta(0.2), pM0(1.0*GeV), m2min(0.09*GeV2),
    m2mini(0.09*GeV2), angmin(0.01), angmini(0.01), pd0(2.5/GeV2/GeV2),
    MaxLoop(100), theCurrentString(0), thefinalSR(0), ntry(0){}


LundFragHandler::LundFragHandler(const LundFragHandler& lfh) 
  : HadronizationHandler(lfh), thePtGen(lfh.thePtGen), theFlGen(lfh.theFlGen ), 
    theZGen(lfh.theZGen), theCollapser(lfh.theCollapser), pWmin0(lfh.pWmin0),
    pK(lfh.pK), pDelta(lfh.pDelta), pM0(lfh.pM0), m2min(lfh.m2min),
    m2mini(lfh.m2mini), angmin(lfh.angmin), angmini(lfh.angmini), pd0(lfh.pd0),
    MaxLoop(lfh.MaxLoop), theCurrentString(0),
    thefinalSR(0), ntry(0){}


LundFragHandler::~LundFragHandler() {
  if ( theCurrentString ) delete theCurrentString;
  if ( theBuffer.size() ) clearBuffer();
} 

void LundFragHandler::
handle(EventHandler & eh, const tPVector & tagged, const Hint & hint) {

  // First find stings and collapse the ones that are too small.
  vector<ColourSinglet> singlets =
    theCollapser->collapse(RemnantDecayer::decayRemnants(tagged, *newStep()),
			   newStep());


  for ( int i = 0, N = singlets.size(); i < N; ++i ) {

    // We don't know how to handle junctions yet.
    if ( singlets[i].nPieces() > 1 ) throw LundFragHdlrException(*this)
      << "We do not know how to handle junction strings yet."
      << Exception::runerror;

    // Perform the actual hadronization and add to the step.
    list<PPtr> hadrons = Hadronize(singlets[i].partons());

    if ( hadrons.empty() ) return;

    newStep()->addDecayProduct(singlets[i].partons().begin(),
			       singlets[i].partons().end(),
			       hadrons.begin(), hadrons.end());
  }

}

// ***********************************************************
//                      Hadronization
// ************************************************************

ParticleList LundFragHandler::Hadronize(const tcPVector& inParticles) {
  ParticleList ret;

  ntry=0;
  if ( theCurrentString ) delete theCurrentString;
  theCurrentString=0;

  m2mini = m2min;
  angmini = angmin;

  for(;;) {
    try {
      if (++ntry > 8*maxLoop() )
	throw LundFragLoopException(*this);
      if ( ntry%maxLoop() == 0 ) {
	m2mini *= 4.0;
	angmini *= 2.0;
	if ( theCurrentString ) delete theCurrentString;
	theCurrentString=0;
     }

      initHadronization(inParticles);
      do{
	getHadron();
      }while(enoughE());
      finalTwoHadrons();
      ret = createParticleList();
      Utilities::transform(ret, cmr);
      return ret;
    }
    catch(Veto) {}
  }
}

// ***********************************************************
//                            Initialization  
// ************************************************************

void LundFragHandler::resetHandler(){  
  GammaM2Solution = false;
  theRightEP.Init();     
  theLeftEP.Init();			
  CurrentEP.Init();	
  thefinalSR = 0;
  newCreatedPD = cPDPtr(); 
  newHadron= Hadron();
  lastHadron= Hadron(); 
}



void LundFragHandler::initHadronization(const tcPVector& inPVec){

  currentBufferIt = theBuffer.begin();

  (inPVec.back()->id() == ParticleID::g)? 
    initClosedString(inPVec) : initOpenString(inPVec); 

  Estatus = true;
}

void  LundFragHandler::initOpenString(const tcPVector& inPVec) {
  if(!theCurrentString){
    theCurrentString = new String(copyRecombine(inPVec.back(), inPVec));

    XhatFwdVector= XhatVector(nSR(), make_pair(1.,1.));
    XhatBwdVector= XhatVector(nSR(), make_pair(0.,0.));

  }else{
    theCurrentString->reset();

    for(int idx=0; idx!=nSR(); ++idx){
      XhatFwdVector[idx].first=1.0; XhatFwdVector[idx].second=1.0;
      XhatBwdVector[idx].first=0.0; XhatBwdVector[idx].second=0.0;
    }
  }

  // Init First End Points
  theRightEP.Init(inPVec.back()->dataPtr(), theCurrentString->firstSR() );
  theLeftEP.Init(inPVec.front()->dataPtr(), theCurrentString->lastSR() ); 

}



void LundFragHandler::initClosedString(const tcPVector& inPVec) {
  if(!theCurrentString){
    theCurrentString = new String(copyRecombine(selectBreakup(inPVec), inPVec));

    XhatFwdVector= XhatVector(nSR(), make_pair(1.,1.));
    XhatBwdVector= XhatVector(nSR(), make_pair(0.,0.));
  }
  else {
    theCurrentString->reset();

    for(int idx=0; idx!=nSR(); ++idx){
      XhatFwdVector[idx].first=1.0; XhatFwdVector[idx].second=1.0;
      XhatBwdVector[idx].first=0.0; XhatBwdVector[idx].second=0.0;
    }
  }
  pickFirstEPts();
}

cPPtr LundFragHandler::selectBreakup(const tcPVector& inPVec) {
  cPPtr first = inPVec.back(); 
  cPPtr next = first;
  tcPPtr nnext;

  Energy2 W2sum=0.0*GeV2;
  do {
    nnext = next->antiColourNeighbour(inPVec.begin(), inPVec.end());
    W2sum += 0.5*(next->momentum() * nnext->momentum());
  } while( ( next = nnext ) && next != first );
  
  //Break the closed string at random
  Energy2 W2rnd = rnd()* W2sum;

  breakup = first; 
    
  while( W2rnd > 0.0*GeV2 ) {
    nnext = breakup->antiColourNeighbour(inPVec.begin(), inPVec.end());
    W2rnd -= 0.5*( breakup->momentum() * nnext->momentum());
    if( W2rnd > 0.0*GeV2 ) breakup = nnext;
  }
  return breakup;
}


void LundFragHandler::pickFirstEPts(){
  //pick at random the starting flavour
  //Warrning : Flavour Generator !!!!!!!

  long RndFlavourId = long(1.0 + (2.0 + SqRatio())*rnd())*(rndbool()? -1: 1);
  
  PDPtr rndPData = getParticleData(RndFlavourId);

  // set random right-left EndPoints
  cPDPtr rightPD;
  generateHadron(rndPData, rightPD);    


  theRightEP.Init(rightPD, theCurrentString->firstSR() );

  cPDPtr leftPD = rightPD->CC();
  theLeftEP.Init(leftPD, theCurrentString->lastSR());

  // pT
  TransverseMomentum initPt = PtGen()->generate();

  //mT2
  Energy2 initSR_E =  theCurrentString->firstSR()->W2();
  Energy2 mT2 = min(25.0*GeV2, 0.1*initSR_E);                   //Units
  
  double z, zR;

  do{
    z = ZGen()->generate(rightPD, leftPD, mT2);
    zR = mT2/(z*initSR_E);
  }while(zR >= 1.0);



  // WARN : call of PYMASS for the quarks here!!!

  Energy2 initGamma = mT2*(1.0 - z)/z;

  //Set Initial EndPoints 
  theRightEP.setPt(initPt);
  theRightEP.Gamma(initGamma);
  theRightEP.SR()->setXrem(1.0 - z, 1.0 - zR);
  //

  // Use first/second to set Xhat Vectors

  XhatFwdVector.front() = make_pair(1.0 - z, 1.0);
  XhatBwdVector.front() = make_pair(zR, 0.0);
  //
  theLeftEP.setPt(-initPt);
  theLeftEP.Gamma(initGamma);
  theLeftEP.SR()->setXrem(z,zR);
  //
  XhatFwdVector.back() = make_pair(1.0,  zR);
  XhatBwdVector.back() = make_pair(0.0, 1.0 - z);
  //
}


// ***********************************************************
//                      produce Hadrons   
// ************************************************************

void LundFragHandler::getHadron(){
 
  Oriented::pickSide(rnd());
  CurrentEP = lastEP();

  newHadron.PData( generateHadron(lastEP().PData(), newCreatedPD) );  
  if(!newHadron.Data) throw Veto();

  // --> Test Heavy Flavour Kfl> 10 !!
  CurrentEP.PData(newCreatedPD);
  CurrentEP.setPt(PtGen()->generate());  

  newHadron.mT2(lastEP(), CurrentEP);

  if(theCurrentString->Wrem2() < 0.1*GeV2) throw Veto();   //WARN : Units 
                                                         //should be here?
  if(theCurrentString->Wrem2() < Wmin2() ) {
    Estatus = false;
    return;
  }

  //
  // Continue Generation
  double z = ZGen()->generate(lastEP().PData(), newCreatedPD, newHadron.mT2() ); 
  //If HQ
  //    - compute re-scaling zprime and mTprime
  //    - test zprime, mTprime   -> Joinning Procedure (c.a.d return ) 
  //    - set z=zprime, mewHadron_mT2= mTprime
  //    - test to skip to 2Hprocedure 
  if ( max(abs(lastEP().PData()->id())%10,
	   abs(lastEP().PData()->id()/1000)%10) >= 4 ||
       max(abs(lastOppEP().PData()->id())%10,
	   abs(lastOppEP().PData()->id()/1000)%10) >= 4 ) {
    Energy2 mtr2 = sqr(lastOppEP().mass() + CurrentEP.mass())
      + sqr(lastOppEP().Px() - CurrentEP.Px())
      + sqr(lastOppEP().Py() - CurrentEP.Py());
    Energy2 mtj2 = newHadron.mT2();
    Energy2 wr2 = theCurrentString->Wrem2();
    Energy2 w12 = sqrt(max(sqr(wr2 - mtj2 - mtr2)
			   - 4.0*mtj2*mtr2, 0.0*GeV2*GeV2));
    z = (wr2 + mtj2 - mtr2 + w12*(2.0*z - 1.0))/(wr2*2.0);
    mtr2 = sqr(lastOppEP().mass() + Wmin0())
      + sqr(lastOppEP().Px() - CurrentEP.Px())
      + sqr(lastOppEP().Py() - CurrentEP.Py());
    if ( (1.0 - z)*(wr2 - mtj2/z) < mtr2 ) {
      Estatus = false;
      return;
    }
  }

  CurrentEP.Gamma((1.0 - z)*(lastEP().Gamma() + newHadron.mT2()/z));   

  // test Stepping is needed  

  if (CurrentSR()->aPrimaryStringRegion() && 
      z*CurrentXremf()*CurrentXremb()*CurrentSR()->W2() >= newHadron.mT2() ) 
    {
    // Continue Hadronization in the Current SR == lastEP
    newHadron.Xf = z * CurrentXremf();
    newHadron.Xb = newHadron.mT2()/(newHadron.Xf*CurrentSR()->W2());
    newHadron.P  = lastEP().pT() +  CurrentEP.pT();
    } 

  else {
    Stepping();
  }

  newHadron.P += newHadron.Xf*CurrentSR()->Pfwd() 
               + newHadron.Xb*CurrentSR()->Pbwd();

  if (newHadron.e() < newHadron.mass()) throw Veto();

  // newHadron.storeMomentum();
  // newHadron.createParticle();

  store(newHadron);  
  loopBack(); 

}


// ***********************************************************
//                      Loop Back
// ************************************************************


void LundFragHandler::loopBack() {

  theCurrentString->updatePtotrem(newHadron.P);
  getLastEP().UpdatedFrom(CurrentEP);  

  int currentIf = CurrentSR()->Ifwd();
  int currentIb = CurrentSR()->Ibwd();

  theCurrentString->updateXremf(currentIf, newHadron.Xf);
  theCurrentString->updateXremb(currentIb, newHadron.Xb);

  Xhatfwd(currentIf) -= newHadron.Xf;
  Xhatbwd(currentIb) += newHadron.Xb;
}


Energy2 LundFragHandler::Wmin2() const{
  Energy  Wmin = Wmin0() + theRightEP.mass() + theLeftEP.mass() 
    + k()*CurrentEP.mass();
  Wmin *= (1.0 + (2.0*rnd() - 1.0)*Delta());
  // --> Test of kfl>10 

  return sqr(max(Wmin, m0() + theRightEP.mass() + theLeftEP.mass())); 
}



// ***********************************************************
//                      Final Two Hadrons 
// ************************************************************

void LundFragHandler::finalTwoHadrons(){

  lastHadron.PData( getHadron(lastOppEP().PData(), newCreatedPD->CC()) );

  if(!lastHadron.Data) throw Veto();
  //lastHadron.createParticle();

  lastHadron.mT2(lastOppEP(), CurrentEP.CC());

  setupCommonFinalRegion();
  solveKinematics();
  return;
}



void LundFragHandler::setupCommonFinalRegion(){

  // set finalSR
  if(CurrentSR() != lastOppEP().SR()){    
    
    thefinalSR = (CurrentSR()->remW2() > lastOppEP().SR()->remW2())?
      CurrentSR() : lastOppEP().SR();
    
    // project final Hadron momentum with repect to final StringRegion axes
    Energy newPx = - firstEP().Px()
      - theCurrentString->PtotRem() * finalSR()->ex(); 
    Energy newPy = - firstEP().Py()
      - theCurrentString->PtotRem() * finalSR()->ey();

    secondEP().setPt(newPx, newPy);

    (finalSR() == CurrentSR())? 
      secondH().mT2(secondEP(), CurrentEP.CC()) :
      secondH().mT2(secondEP(), CurrentEP);
  }
  else { 
    thefinalSR =  CurrentSR();
  }

}

void LundFragHandler::solveKinematics(){

  double rb = 2.0*(theCurrentString->PtotRem()*finalSR()->Pbwd()) 
                  /finalSR()->W2();
  double rf = 2.0*(theCurrentString->PtotRem()*finalSR()->Pfwd()) 
                  /finalSR()->W2();

  Energy2 Erem2 = theCurrentString->Wrem2() 
                + sqr(lastEP().Px() + lastOppEP().Px())
                + sqr(lastEP().Py() + lastOppEP().Py());
  
  double fd = (newHadron.mT()+lastHadron.mT())/sqrt(Erem2);
  
  if(fd >= 1.0) throw Veto();
  //Extra Heavy Flavor Test Here

  double d = d0()*sqr((newHadron.mT2() + lastHadron.mT2()) );
  double Preverse = 0.5*exp(max(-50.0, d*log(fd)));

  Energy2 fa = Erem2 + newHadron.mT2() - lastHadron.mT2();
  Energy2 delta = sqrt(max(0.0*GeV2*GeV2, fa*fa - 4.0*Erem2*newHadron.mT2()));

  Energy2 fb = sign(delta, rnd() - Preverse); 

  //Diquark Test Here : Change FB


  //Momenta of the last 2 hadrons
//   newHadron.P = (lastEP().Px() + CurrentEP.Px()) * CurrentSR()->ex() 
//               + (lastEP().Py() + CurrentEP.Py()) * CurrentSR()->ey() 
//               + (0.5/Erem2)*( rb*(fa+fb)*finalSR()->Pfwd() 
// 		             + rf*(fa-fb)*finalSR()->Pbwd()      );

  newHadron.P = (lastEP().Px() + CurrentEP.Px()) * finalSR()->ex() 
              + (lastEP().Py() + CurrentEP.Py()) * finalSR()->ey() 
              + (0.5/Erem2)*( rb*(fa+fb)*finalSR()->Pfwd() 
	                    + rf*(fa-fb)*finalSR()->Pbwd()      );

  lastHadron.P = theCurrentString->PtotRem() - newHadron.P;



  if( newHadron.e() < newHadron.mass() || 
      lastHadron.e() < lastHadron.mass()  ) throw Veto();
  
  //newHadron.storeMomentum();
  //lastHadron.storeMomentum();

  //newHadron.createParticle();
  //lastHadron.createParticle();

  store(newHadron);
  store(lastHadron, Oriented::OppDir());
  return;
}


// ***********************************************************
//                          Stepping Procedure 
// ************************************************************
 

void  LundFragHandler::Stepping(){

  GammaM2Solution = false; 

  Pzero = LorentzMomentum();

  if(CurrentSR()->aPrimaryStringRegion()) Step();  

  do{
    // re-express pt of new created pair in the newCurrentSR() axes 
    //

    Energy pxp = CurrentEP.Px()*( CurrentSR()->ex() *lastEP().SR()->ex())
                 + CurrentEP.Py()*( CurrentSR()->ex() *lastEP().SR()->ey());

    Energy pyp = CurrentEP.Px()*( CurrentSR()->ey() *lastEP().SR()->ex())
                 +CurrentEP.Py()*( CurrentSR()->ey() *lastEP().SR()->ey());
    
    Energy2 delta = abs(sqr(pxp) + sqr(pyp) -
			sqr(CurrentEP.Px()) - sqr(CurrentEP.Py()));


    if( delta < 0.01*GeV2 ) {
      CurrentEP.setPt(-pxp, -pyp);
    } 

    solveGammaM2System(); 

    //check Solutions - Step to new String Region if needed
    if(newHadron.Xb > CurrentSR()->Xremb()){  
      Step();
    } else{ 
      if (newHadron.Xf > CurrentSR()->Xremf()){
	stepDown();
      }else{ 
	GammaM2Solution =true;
      }
    }         
  }while(!aSolution());  

  return; 
}


void LundFragHandler::Step(){
  //stepUp()
  //Update p0  
  int bidx = CurrentSR()->Ibwd();
  Pzero += theCurrentString->Xremb(bidx)*theCurrentString->Pbwd(bidx);
  //

  //Update of newHadron.Xb here ?!
  Xhatbwd(CurrentSR()->Ibwd()) = 1.0;
  CurrentEP.stepUp();

  if(inconsistentBreakupRegions() ) throw Veto();

  if(CurrentSR()->W2() < 2.0e-2*GeV2 ) {  //WARN: Units
    stepDown();
    //Update of newHadron.Xf here ?! 
  } 
} 

void LundFragHandler::stepDown() {
  //Update p0  
  int fidx = CurrentSR()->Ifwd();
  Pzero +=  theCurrentString->Xremf(fidx)*theCurrentString->Pfwd(fidx);
  //
  
  Xhatfwd(CurrentSR()->Ifwd()) = 0.0;
  CurrentEP.stepDown();
  if(inconsistentBreakupRegions() ) throw Veto();
} 


void  LundFragHandler::solveGammaM2System(){
  //
  // compute P (= pT +p0()) + ci mass coefficients

  //  newHadron.P  = CurrentEP.pT() + lastEP().pT() 
  //  + p0(lastEP().SR(), CurrentSR());

  newHadron.P  = CurrentEP.pT() + lastEP().pT() + Pzero;


  Energy2 massC0 = newHadron.P*newHadron.P;
  Energy2 massC1 = 2.0*newHadron.P*CurrentSR()->Pfwd();
  Energy2 massC2 = 2.0*newHadron.P*CurrentSR()->Pbwd();
  Energy2 massC3 = CurrentSR()->W2();

  //
  //  Compute Coefficients for Gamma expression
  //

  Energy2 gamC0 = 0.0*GeV2, gamC1 = 0.0*GeV2,
    gamC2 = 0.0*GeV2, gamC3 = 0.0*GeV2;
  Energy2 StepE = 0.0*GeV2;

  
  for(OIndex bidx=CurrentSR()->Ifwd() ; bidx<=CurrentSR()->Ibwd(); ++bidx){    
    for(OIndex fidx = CurrentSR()->Ifwd(); fidx<=bidx; ++fidx){

      StepE=2.0*theCurrentString->Pfwd(fidx)
	*theCurrentString->Pbwd(bidx);
    
      gamC0 += Xhatbwd(bidx)*Xhatfwd(fidx)*StepE;
    
      if(fidx == CurrentSR()->Ifwd()) 
	gamC1 -= Xhatbwd(bidx)*StepE;
      
      if(bidx == CurrentSR()->Ibwd()) 
	gamC2 += Xhatfwd(fidx)*StepE;
    
      if(fidx == CurrentSR()->Ifwd() && bidx == CurrentSR()->Ibwd() )
	gamC3 -= StepE;
    }
  }
  

  //  solve (m2, gamma) system solution in the CurrentSR
  //
  Energy4 s1 = massC2*gamC3 - massC3*gamC2;

  if(abs(s1) < (1.e-4*GeV2*GeV2) ) throw Veto();       //WARN : Units


  Energy4 s2 =  massC3*(CurrentEP.Gamma() - gamC0) 
             - massC1*gamC2 - gamC3*(newHadron.m2() - massC0)
             + gamC1*massC2; 

  Energy4 s3 = massC1*(CurrentEP.Gamma() - gamC0)
              - gamC1*(newHadron.m2() - massC0);
  
  newHadron.Xb = 0.5*(sqrt(max(0.0, (s2*s2 - 4.0*s1*s3)/sqr(s1))) - s2/s1);

  if( (massC1 + massC3*newHadron.Xb) <= 0.0*GeV2)  throw Veto();   //WARN: Units
 

  newHadron.Xf =  (newHadron.m2()-massC0 - massC2*newHadron.Xb)  
                 / (massC1 + massC3*newHadron.Xb);
  return;
}


LorentzMomentum
LundFragHandler::p0(cStringRegionPtr initSR, cStringRegionPtr finalSR) {
  LorentzMomentum Pzero;
  
  for(OIndex fidx = initSR->Ifwd(); fidx != finalSR->Ifwd(); ++fidx){
    Pzero +=  theCurrentString->Xremf(fidx)*theCurrentString->Pfwd(fidx);
  }

  for(OIndex bidx = initSR->Ibwd(); bidx!= finalSR->Ibwd(); ++bidx){
    Pzero +=  theCurrentString->Xremb(bidx)*theCurrentString->Pbwd(bidx);
  }

  return Pzero;
}



// ****************************************************************
//                   Flavour generating functions 
// ****************************************************************

tcPDPtr LundFragHandler::
generateHadron(tcPDPtr inPDPtr, cPDPtr& newPDPtr, long curtainQid ) {
  return FlavourGen()->generateHadron(inPDPtr, newPDPtr);
}

tcPDPtr LundFragHandler::getHadron(tcPDPtr inPD1, tcPDPtr inPD2) {
  return FlavourGen()->getHadron(inPD1, inPD2);
}



// ****************************************************************
//                  Storage of produced particles 
// ****************************************************************

void LundFragHandler::store(Hadron& currentH, int Dir){  
  currentH.ProductionSide = Dir;
 
  if(currentBufferIt != theBuffer.end()) {
    *currentBufferIt = currentH;
    ++currentBufferIt;
  }else{
    theBuffer.push_back(currentH);
    currentBufferIt=theBuffer.end();
  }
}



ParticleList LundFragHandler::createParticleList(){
  ParticleList theList;
  ParticleListIt StorageIt = theList.end();

  for(BufferIt Hit = theBuffer.begin(); Hit != currentBufferIt; ++Hit){
    if(Hit->ProductionSide == Oriented::right){
      theList.insert(StorageIt, Hit->createParticle());
    }else{
      theList.insert(StorageIt, Hit->createParticle());
      --StorageIt;
    }
  }

  return theList;
}


// **********************************************************
//                            DB 
// **********************************************************

void LundFragHandler::showEP() const {
    cout<<"[rightEP ] :"<<endl; theRightEP.echo();
    cout<<"\n [leftEP ] :"<<endl; theLeftEP.echo();
}


void LundFragHandler::echoXhat(cStringRegionPtr sr) {

  cout<<"Xhatfwd("<<sr->Ifwd()<<")= "
      <<Xhatfwd(sr->Ifwd())<<"\t (SR_Index)"<<"\t dir="<<Oriented::Dir()<<endl;
    
  cout<<"Xhatbwd("<<sr->Ibwd()<<")= "
      <<Xhatbwd(sr->Ibwd())<<"\t (SR_Index)"<<endl;
  
}

PVector LundFragHandler::copyRecombine(tcPPtr first, const tcPVector & inPVec) {
  PVector partons;
  cPPtr next = first;
  bool loop = first->id() == ParticleID::g;
  do {
    partons.push_back(next->data().produceParticle(next->momentum()));
  } while ( ( next = next->antiColourNeighbour(inPVec.begin(), inPVec.end()) )
	    && next != first );

  cmr = Utilities::boostToCM(partons.begin(), partons.end());
  cmr.invert();

  while ( partons.size() > 2 ) {
    Energy2 m2m = 2.0*m2mini;
    int isel = -1;
    for ( int i = 0, j = 1, N = partons.size(); i < N; ++i, j = (i + 1)%N ) {
      if ( j == 0 && !loop ) break;
      Momentum3 p1 = partons[i]->momentum();
      Momentum3 p2 = partons[j]->momentum();
      Energy2 pap = sqrt(p1.mag2()*p2.mag2());
      Energy2 pvp = p1*p2;
      Energy2 m2i = 4.0*sqr(pap - pvp)/
	max(1.0e-6*GeV2, sqr(angmini)*pap + 2.0*(pap - pvp));
      if ( m2i < m2m ) {
	m2m = m2i;
	isel = i;
      }
    }
    if ( isel < 0 ) break;
    int jsel = (isel + 1)%partons.size();
    int ksel = (isel - 1)%partons.size();
    if ( partons[isel]->id() == ParticleID::g &&
	 partons[jsel]->id() != ParticleID::g ) swap(isel, jsel);
    partons[isel]->setMomentum(partons[isel]->momentum() +
			       partons[jsel]->momentum());
    if ( partons[isel]->id() == ParticleID::g ) {
      LorentzRotation r =
	Utilities::getBoostToCM(make_pair(partons[isel], partons[ksel]));
      Energy2 s = (partons[isel]->momentum() + partons[ksel]->momentum()).m2();
      SimplePhaseSpace::CMS(partons[isel], partons[ksel], s);
      r.invert();
      partons[isel]->transform(r);
      partons[ksel]->transform(r);
    } else {
      partons[isel]->rescaleMass();
    }
    partons.erase(partons.begin() + jsel);
  }
  return partons;
}


// ************************************************************
//                    Interfaced functions
// ************************************************************

ClassDescription<LundFragHandler> LundFragHandler::initLundFragHandler;

void LundFragHandler::persistentOutput(PersistentOStream & os) const {
  os << thePtGen << theFlGen << theZGen << theCollapser << ounit(pWmin0, GeV)
     << pK << pDelta << ounit(pM0, GeV) << ounit(m2min, GeV2)
     << ounit(m2mini, GeV2) << angmin << angmini << ounit(pd0, 1.0/(GeV2*GeV2))
     << MaxLoop;
}

void LundFragHandler::persistentInput(PersistentIStream & is, int) {
  is >> thePtGen >> theFlGen >> theZGen >> theCollapser >> iunit(pWmin0, GeV)
     >> pK >> pDelta >> iunit(pM0, GeV) >> iunit(m2min, GeV2)
     >> iunit(m2mini, GeV2) >> angmin >> angmini >> iunit(pd0, 1.0/(GeV2*GeV2))
     >> MaxLoop;
}


void LundFragHandler::Init(){

  static ClassDocumentation<LundFragHandler> interfaceDocumentation
    ("This is the main handler class for Lund string fragmentation.",
     "Hadronization performed according to the "
     "Lund String Fragmentation scheme \\cite{And83}",
     "\\bibitem{And83} B. Andersson, G. Gustafson, "
     "G. Ingelman and T. Sj""\\""\"ostrand"
     "\n Phys.~Rep.~{\\bf 97} (1983) 31.");
  
  static Reference<LundFragHandler,PtGenerator> interfaceLundPtGen
    ("PtGenerator",
     "The PtGenerator object used by the Lund fragmentation handler. "
     "Should derive from Pythia7::LundPtGenerator.",
     &LundFragHandler::thePtGen, false, false, true, false);

  static Reference<LundFragHandler,ZGenerator> interfaceLundZGen
    ("ZGenerator",
     "Pointer to the ZGenerator object used by the fragmentation handler. "
     "Should derive from Pythia7::LundZGenerator.",
     &LundFragHandler::theZGen, false, false, true, false);

  static Reference<LundFragHandler,ClusterCollapser> interfaceCollapser
    ("ClusterCollapser",
     "Pointer to a ThePEG::ClusterCollapser object used by the Lund "
     "fragmentation handler to avoid too small strings in the fragmentation.",
     &LundFragHandler::theCollapser, false, false, true, false);

  static Reference<LundFragHandler,LundFlavourGenerator> interfaceLundFlGen
    ("FlavourGenerator",
     "Pointer to the FlavourGenerator object used by the fragmentation handler"
     ". Should derive from Pythia7::LundFlavourGenerator.",
     &LundFragHandler::theFlGen, false, false, true, false);

  static Parameter<LundFragHandler, Energy> interfaceWmin0
    ("Wmin0",
     "Used to define, together with quark masses, the remaining energy" 
     "\nbelow which the fragmentation of a jet system is stopped and two final"
     "\nhadrons formed. For the Standard Lund symmetric fragmentation scheme "
     "[in GeV]",
     &LundFragHandler::pWmin0, GeV, 0.8*GeV, 0.0*GeV, 10.0*GeV,
     false, false, true);

  static Parameter<LundFragHandler, double> interfaceK
    ("k",
     "Represents the dependence on the mass of the final quark pair for "
     "defining the stopping point of the fragmentation [noUnit]."
     "Is strongly correlated to the choice of \\f$W_{\\min0}\\f$",
     &LundFragHandler::pK, 2.0, 0.0, 5.0, false, false, true);

  static Parameter<LundFragHandler, double> interfaceDelta
    ("delta",
     "Relative width of the smearing of the stopping point energy [noUnit].",
     &LundFragHandler::pDelta, 0.2, 0.0, 1.0, false, false, true);

  static Parameter<LundFragHandler, Energy> interfaceM0
    ("m_0",
     "Is with quark masses added, used to define the minimum allowable"
     "energy of a colour-singlet jet system [GeV].",
     &LundFragHandler::pM0, GeV, 1.0*GeV, 0.0*GeV, 5.0*GeV,
     false, false, true);

  static Parameter<LundFragHandler, InvEnergy4> interfaceD0
    ("d_0",
     "Refers to the probability for reverse rapidity ordering of the final"
     "\ntwo hadrons in the standard Lund symmetric fragmentation scheme "
     "[GeV^-4].",
     &LundFragHandler::pd0, 1.0/GeV2/GeV2, 2.5/GeV2/GeV2, 0.0/GeV2/GeV2,
     5.0/GeV2/GeV2, false, false, true);

  static Parameter<LundFragHandler, long> interfaceMaxLoop
    ("maxLoop",
     "Number of maximun try alowed for a String Fragmentation [noUnit]."
     "The hadronization of a given string will be tried 8*maxLoop times "
     "and the M2Min and MinAngel will be increased every maxLoop attempts.",
     &LundFragHandler::MaxLoop, 100, 0, 1000, false, false, true);
  interfaceMaxLoop.setHasDefault(false);

  static Parameter<LundFragHandler,Energy2> interfaceM2Min
    ("M2Min",
     "The effective cut-off in squared mass, below which partons may be "
     "recombined to simplify (machine precision limited) kinematics of "
     "string fragmentation. (Default chosen to be of the order of a light "
     "quark mass, or half a typical light meson mass.)",
     &LundFragHandler::m2min, GeV2, 0.09*GeV2, 0.0*GeV2, 10.0*GeV2,
     true, false, true);


  static Parameter<LundFragHandler,double> interfaceMinAngle
    ("MinAngle",
     "The effective angular cut-off in radians for recombination of "
     "partons, used in conjunction with M2Min.",
     &LundFragHandler::angmin, 0.01, 0.0, 1.0,
     true, false, true);

  interfaceLundZGen.rank(10);
  interfaceLundPtGen.rank(9);
  interfaceLundFlGen.rank(8);

}

LundFragInitException::LundFragInitException(const  LundFragHandler & lfh) {
  theMessage <<" The Fragmentation Handler "<< lfh.name() << "does not have "
	     <<" pointers to the necessary generators correctly setup ";
  severity(maybeabort);
}

LundFragLoopException::LundFragLoopException(const LundFragHandler & lfh) { 

  theMessage << "The maximum number of attempts ("<< lfh.maxLoop() <<") " 
	     << "to fragment a string exceeded" 
	     <<"\nin Fragmentation process handled by "<< lfh.name() <<"."
	     << "\nThe corresponding event has been rejected." ;

  severity(eventerror);    
}

