#ifndef THEP8I_RopeUserHooks_H
#define THEP8I_RopeUserHooks_H

#include "Pythia.h"
#include "ThePEG/EventRecord/Particle.h"
#include "TheP8I/Hadronization/ParameterHandler.h"
#include "TheP8I/Hadronization/Ropewalk.h"
#include "ThePEG/Utilities/Throw.h"



namespace TheP8I {

class RopeUserHooks : public Pythia8::UserHooks {

// Convenient typedefs
typedef map<string,double> PytPars;
typedef map<Energy2, Ropewalk::Dipole *> DipMass;
typedef map<tcPPtr,DipMass> FlavourEnd;
 public:

 RopeUserHooks() : _ph(NULL), _h(-1.0), _m0(0*GeV), _pTcut(100*GeV), _r0(1*femtometer), _window(false) {
  }

  ~RopeUserHooks() {
  }

  virtual bool canChangeFragPar() {
    return true;
  }

  virtual bool doChangeFragPar(Pythia8::StringFlav* flavPtr, Pythia8::StringZ* zPtr, Pythia8::StringPT * pTPtr, int endFlavour, double m2Had, vector<int> iParton) {
    // We're not using it here
    (void) iParton;
    // Get new parameters
    PytPars newPar = fetchParameters(endFlavour,m2Had);
    if(newPar.find("null") != newPar.end()) {
      Throw<RopeException>()
	<< "Problem fetching parameters in RopeUserHook. "
	<< "Ropes switched off in this string."
	<< Exception::warning;
      _h = 1.0;
      newPar = fetchParameters(endFlavour,m2Had);
      _h = -1.0;
    }
    for(PytPars::iterator itr = newPar.begin(); itr!=newPar.end();++itr) settingsPtr->parm(itr->first,itr->second);

    // Re-initialize all three
    flavPtr->init(*settingsPtr,rndmPtr);
    zPtr->init(*settingsPtr,*particleDataPtr,rndmPtr);
    pTPtr->init(*settingsPtr,*particleDataPtr,rndmPtr);

    return true;
  }

  void setParameterHandler(ParameterHandler * ph){
    _ph = ph;
  }


  void setEnhancement(double h){
    // Will overrule others if set.
    _h = h;
  }

  void setWindow(bool window, Energy pTcut){
    // Can set a window where rope should not apply
    _window = window;
    _pTcut = pTcut;
  }

  bool setDipoles(const pair<ColourSinglet *, vector<Ropewalk::Dipole *> > *  dm, Energy m0, Length r0, bool throwHadr = false, double alpha = 1.0){

    // If we set the dipoles, we should also use them.
    _h = -1.0;
    _m0 = m0;
    _r0 = r0;
    _throwHadr = throwHadr;
    _alpha = alpha;

    vector<Ropewalk::Dipole*> dip = dm->second;
     // Closed gluon loop -- we will return false and let StringFragmentation figure out what to do
    if(dip.front()->pa->id() == 21 && dip.front()->pc->id() == 21 &&
        dip.back()->pa->id() == 21 && dip.back()->pc->id() == 21)
            return false;
        
   
    // Get flavour of quark ends -- get the one that is not a gluon.
    // We always assign left to the front of the 
    // dipole vector and right to the back. 
    DipMass leftC;
    DipMass rightC;
    leftC.clear();
    rightC.clear();
    if(dip.size()==1){
      _leftP = dip[0]->pc;
      _rightP = dip[0]->pa;
    }
    else{
      _leftP = (dip.front()->pa->id() != 21 ? dip.front()->pa : dip.front()->pc); 
      _rightP = (dip.back()->pc->id() != 21 ? dip.back()->pc : dip.back()->pa); 
    }
    if (_leftP->id() == _rightP->id() )
      Throw<RopeException>()
	<< "Flavours in strings corrupted. "
	<< "This is a serious error - please contact the authors."
	<< Exception::abortnow;

    // A dipole gets all of quark mass, half of gluon
    // The end quarks energy m0 

    // Start by adding half of the first quark
    // from both left and right

    LorentzMomentum momL(0*GeV,0*GeV,0*GeV,0*GeV);
    LorentzMomentum momR(0*GeV,0*GeV,0*GeV,0*GeV);

    
    for(size_t i = 0, N = dip.size(); i < N; ++i ){
      Ropewalk::Dipole * dPtrL = dip[i];
      Ropewalk::Dipole * dPtrR = dip[N - i - 1];

      if(!dPtrL || dPtrL->broken || !dPtrR || dPtrR->broken){
      Throw<RopeException>()
	<< "Broken dipole in RopeUserHooks. "
	<< "This is a serious error - please contact the authors."
	<< Exception::abortnow;
        return false;
      }
      momL += 0.5*(dPtrL->pc->momentum());
      momL += 0.5*(dPtrL->pa->momentum());
      momR += 0.5*(dPtrR->pc->momentum());
      momR += 0.5*(dPtrR->pa->momentum());

      if(i == 0){
        momL += 0.5 * _leftP->momentum();
        momR += 0.5 * _rightP->momentum();
      }
      if(i == N - 1){
        momR += 0.5 * _leftP->momentum();
        momL += 0.5 * _rightP->momentum();
      }
      leftC.insert(DipMass::value_type(momL.m2(),dPtrL));
      rightC.insert(DipMass::value_type(momR.m2(),dPtrR));
    }

    flavourDipoles.clear();
    flavourDipoles.insert( FlavourEnd::value_type(_leftP,leftC) );
    flavourDipoles.insert( FlavourEnd::value_type(_rightP,rightC) );

    return true;
  } 

  bool setDipoles(const pair<ColourSinglet * const, vector<Ropewalk::Dipole *> > *  dm,
		  Energy m0, Length r0, bool throwHadr = false, double alpha = 1.0){

    // If we set the dipoles, we should also use them.
    _h = -1.0;
    _m0 = m0;
    _r0 = r0;
    _throwHadr = throwHadr;
    _alpha = alpha;

    vector<Ropewalk::Dipole*> dip = dm->second;
     // Closed gluon loop -- we will return false and let StringFragmentation figure out what to do
    if(dip.front()->pa->id() == 21 && dip.front()->pc->id() == 21 &&
        dip.back()->pa->id() == 21 && dip.back()->pc->id() == 21)
            return false;
        
   
    // Get flavour of quark ends -- get the one that is not a gluon.
    // We always assign left to the front of the 
    // dipole vector and right to the back. 
    DipMass leftC;
    DipMass rightC;
    leftC.clear();
    rightC.clear();
    if(dip.size()==1){
      _leftP = dip[0]->pc;
      _rightP = dip[0]->pa;
    }
    else{
      _leftP = (dip.front()->pa->id() != 21 ? dip.front()->pa : dip.front()->pc); 
      _rightP = (dip.back()->pc->id() != 21 ? dip.back()->pc : dip.back()->pa); 
    }
    if (_leftP->id() == _rightP->id() )
      Throw<RopeException>()
	<< "Flavours in strings corrupted. "
	<< "This is a serious error - please contact the authors."
	<< Exception::abortnow;

    // A dipole gets all of quark mass, half of gluon
    // The end quarks energy m0 

    // Start by adding half of the first quark
    // from both left and right

    LorentzMomentum momL(0*GeV,0*GeV,0*GeV,0*GeV);
    LorentzMomentum momR(0*GeV,0*GeV,0*GeV,0*GeV);

    
    for(size_t i = 0, N = dip.size(); i < N; ++i ){
      Ropewalk::Dipole * dPtrL = dip[i];
      Ropewalk::Dipole * dPtrR = dip[N - i - 1];

      if(!dPtrL || dPtrL->broken || !dPtrR || dPtrR->broken){
      Throw<RopeException>()
	<< "Broken dipole in RopeUserHooks. "
	<< "This is a serious error - please contact the authors."
	<< Exception::abortnow;
        return false;
      }
      momL += 0.5*(dPtrL->pc->momentum());
      momL += 0.5*(dPtrL->pa->momentum());
      momR += 0.5*(dPtrR->pc->momentum());
      momR += 0.5*(dPtrR->pa->momentum());

      if(i == 0){
        momL += 0.5 * _leftP->momentum();
        momR += 0.5 * _rightP->momentum();
      }
      if(i == N - 1){
        momR += 0.5 * _leftP->momentum();
        momL += 0.5 * _rightP->momentum();
      }
      leftC.insert(DipMass::value_type(momL.m2(),dPtrL));
      rightC.insert(DipMass::value_type(momR.m2(),dPtrR));
    }

    flavourDipoles.clear();
    flavourDipoles.insert( FlavourEnd::value_type(_leftP,leftC) );
    flavourDipoles.insert( FlavourEnd::value_type(_rightP,rightC) );

    return true;
  } 

public:

  /**
   * Exception class.
   */
  struct RopeException: public Exception {};

 private:

  PytPars fetchParameters(int endFlavour, double m2Had){
    
    // Did we remember to set the callback pointer? 
    // Do not just construct a new one, as there will be no sensible initial values, resulting
    // in very subtle errors.
    if(!(_ph)){
      Throw<RopeException>()
	<< "Missing parameter handler in RopeUserHooks. "
	<< "This is a serious error - please contact the authors."
	<< Exception::abortnow;
      PytPars p;
      p.insert(make_pair<const string,double>("null",0.));
      return p;
    }
    
    // If we have manually set enhancement, we should just use that, and avoid complicated behavior.
    if(_h > 0) return _ph->GetEffectiveParameters(_h);

    // Test the string ends...
    if( (endFlavour == _leftP->id() && endFlavour == _rightP->id()) ||
         (endFlavour != _leftP->id() && endFlavour != _rightP->id()) ){
      Throw<RopeException>()
	<< "Flavours in strings corrupted. "
	<< "This is a serious error - please contact the authors."
	<< Exception::abortnow;
      PytPars p;
      p.insert(make_pair<const string,double>("null",0.));
      return p;
      }
    
    // Find the right end
    tcPPtr thisEnd = (_leftP->id() == endFlavour ? _leftP : _rightP);
    FlavourEnd::iterator feItr = flavourDipoles.find( thisEnd);
    
    if( feItr == flavourDipoles.end()){
      Throw<RopeException>()
	<< "Could not find dipole end. "
	<< "Ropes switched off in this string."
	<< Exception::warning;
      PytPars p;
      p.insert(make_pair<const string,double>("null",0.));
      return p;
    }
    DipMass dms = feItr->second;

    // Find the first dipole where the total invariant mass sq. from the string
    // (left or right) exceeds  m2Had    
    DipMass::iterator dmsItr = dms.lower_bound(m2Had*GeV2);
    if( dmsItr == dms.end()){
      Throw<RopeException>()
	<< "Could not get correct dipole. mpythia: " << m2Had
	<< " " << (--dmsItr)->first/GeV2 << ". "
	<< "Ropes switched off in this string."
	<< Exception::warning;
      PytPars p;
      p.insert(make_pair<const string,double>("null",0.));
      return p;
    }  

    Ropewalk::Dipole * dPtr = dmsItr->second; 
    if(!dPtr){
      Throw<RopeException>()
	<< "Missing dipole in RopeUserHooks. "
	<< "This is a serious error - please contact the authors."
	<< Exception::abortnow;
      PytPars p;
      p.insert(make_pair<const string,double>("null",0.));
      return p;
    }
      
    // Find out how long in we are on the dipole
    double dipFrac;
    double mBig = dmsItr->first/GeV2;

    if(m2Had == 0) dipFrac = 0;
    else if( dmsItr == dms.begin() ) dipFrac = sqrt(m2Had/mBig);
    else{
      double mSmall = double((--dmsItr)->first/GeV2);        
      dipFrac = (sqrt(m2Had) - sqrt(mSmall)) / (sqrt(mBig) - sqrt(mSmall));
    } 
    // We need dipFrag to be fraction from the pc parton in the dipole.
    // Check if the starting end is pa or pc. If it is already a pc, do noting,
    // otherwise take fraction from other side
    if( thisEnd == dms.lower_bound(0*GeV2)->second->pa) dipFrac = 1 - dipFrac;

    // use the window here
    if(_window && ((dPtr->pc->momentum().perp() > _pTcut && dipFrac < 0.5) ||
            (dPtr->pa->momentum().perp() > _pTcut && dipFrac >=0.5)) ){
          return _ph->GetEffectiveParameters(1.0);
      }

    dPtr->reinit(dipFrac,_r0,_m0);
    if (_throwHadr){
      (*dPtr).hadr = true;
      return _ph->GetEffectiveParameters(dPtr->firstKappa(_alpha));
    }
    //  cout << "p = " << dPtr->p << " q = " << dPtr->q << " kappa " <<  << endl;  
    return _ph->GetEffectiveParameters(dPtr->kappaEnhancement());    
  
}

  tcPPtr _leftP, _rightP;
  FlavourEnd flavourDipoles;
  ParameterHandler * _ph;
  double _h, _alpha;
  Energy _m0, _pTcut;
  Length _r0;
  bool _throwHadr, _window;
};

}
#endif
