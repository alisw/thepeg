// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the StringFragmentation class.
//

#include "StringFragmentation.h"
#include "Ropewalk.h"
#include "RandomAverageHandler.h"
#include "RandomHandler.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ParVector.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/PDT/RemnantDecayer.h"   
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/Step.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DebugItem.h"
#include "ThePEG/Utilities/EnumIO.h"
#include "ThePEG/Utilities/MaxCmp.h"
#include "ThePEG/Utilities/HoldFlag.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include <algorithm>
using namespace TheP8I;


StringFragmentation::StringFragmentation()
  :  pythia(), fScheme(0), stringR0(0.5*femtometer), stringm0(1.0*GeV), junctionDiquark(1.0), alpha(1.0), average(true),
     stringyCut(2.0), stringpTCut(6*GeV), fragmentationMass(0.134),
     baryonSuppression(0.5), window(false), throwaway(false), _eventShapes(0),
     useThrustAxis(0), opHandler(0), maxTries(2), doAnalysis(0), analysisPath(""),
#include "StringFragmentation-init.h"
      {
     }

StringFragmentation::~StringFragmentation() {

  if ( opHandler ) delete opHandler;

}

void StringFragmentation::handle(EventHandler & eh, const tPVector & tagged,
	    const Hint & h) {
  ++nev;
  static unsigned long lastEventId = 0;
  bool secondary = false;
  // Get the event
  tPVector all = RemnantDecayer::decayRemnants(tagged, *newStep());
  HoldFlag<int> oldFS(fScheme, fScheme);
  // Prevent funny behavior if hadronizing decays of eg. Upsilon.
  if ( newStep()->collision()->uniqueId == lastEventId ){
    secondary = true;
    if(fScheme == 4 || fScheme == 5 || fScheme == 1) fScheme = 99;
    else fScheme = 0;
  }
  else 
    lastEventId = newStep()->collision()->uniqueId;

  if(_eventShapes) _eventShapes->reset(all);
  vector<ColourSinglet> singlets;
  singlets.clear();

  if ( theCollapser && !secondary )
    singlets =  theCollapser->collapse(all, newStep());
  else
    singlets = ColourSinglet::getSinglets(all.begin(), all.end());
  // Goto correct hadronization scheme
  switch(fScheme){
    case 0: // Pythia
      {
      hadronizeSystems(*opHandler->GetPythiaPtr(1.0,nev%10000==0), singlets, all);
      break;
      }
    case 1: // For playing around and printing stuff
      {
        if(throwaway){
          Ropewalk ropewalk_init(singlets, stringR0, stringm0, baryonSuppression, throwaway, false);
          vector< pair<ColourSinglet,double> > allStrings = ropewalk_init.getSinglets(stringyCut);
          tPVector allParticles;
          for(vector< pair<ColourSinglet,double> >::iterator sItr = allStrings.begin(); sItr != allStrings.end(); ++sItr)
              for(tcPVector::iterator pItr = sItr->first.partons().begin(); pItr != sItr->first.partons().end(); ++pItr)
                 allParticles.push_back(const_ptr_cast<tPPtr>(*pItr));
          singlets =  theCollapser->collapse(allParticles, newStep());   
        }
        Ropewalk ropewalk(singlets, stringR0, stringm0, baryonSuppression, false, false);


        vector< pair<ColourSinglet,double> > strings = ropewalk.getSinglets(stringyCut);

	static ofstream out(( generator()->runName() + ".txt").c_str());        

	if(!isnan(ropewalk.lambdaSum()+ropewalk.getkW()+ropewalk.getNb()))
	out << generator()->currentEvent()->weight() << " " << ropewalk.lambdaSum() << " " << ropewalk.getkW() << " " << ropewalk.getNb() << endl;
		//vector<double> mspec = ropewalk.mspec;
		//for(int i = 0, N = mspec.size(); i < N; ++i) out <<  generator()->currentEvent()->weight() << " " <<  mspec[i] << endl; 
	break;

    }
    case 2: // Ropewalk/Dipole
     {
     Ropewalk ropewalk(singlets, stringR0, stringm0, junctionDiquark, throwaway, false);
     vector< pair<ColourSinglet,double> > strings = ropewalk.getSinglets(stringyCut);
     vector<ColourSinglet> toHadronization;    
     for(int i = 0, N = strings.size(); i < N; ++i ){
        Pythia8Interface * pytp = opHandler->GetPythiaPtr(strings[i].second,nev%10000==0);
        toHadronization.clear();
        toHadronization.push_back(strings[i].first);  
        hadronizeSystems(*pytp,toHadronization,all);
      }
      break;
      }
     case 3: // Pipes with average
     {
      RandomAverageHandler avg_enhancer(throwaway);
      // TODO: Implement throwaway functionality in this
      avg_enhancer.clear();
      vector<StringPipe> pipes;
      // Make all the pipes
      for(vector<ColourSinglet>::iterator sItr = singlets.begin(); sItr!=singlets.end(); ++sItr) 
          pipes.push_back(StringPipe(&(*sItr),stringR0,_eventShapes));
      
      avg_enhancer.SetEvent(pipes);        
      for(vector<StringPipe>::iterator it = pipes.begin(); it!=pipes.end(); ++it){
          vector<ColourSinglet> toHadronization;    
          double h = avg_enhancer.KappaEnhancement(*it);
          if(h > 0){
            toHadronization.push_back(*(*it).GetSingletPtr());    
           forcerun = false;
           if(!hadronizeSystems(*(opHandler->GetPythiaPtr(h,nev%10000==0)),toHadronization,all)) 
              hadronizeSystems(*(opHandler->GetPythiaPtr(1.0,false)),toHadronization,all);
          }
        }
      break;
     }
     case 4: // Average kappa over whole strings
     {
       Ropewalk ropewalk(singlets, stringR0, stringm0, baryonSuppression,throwaway, false);
        vector< pair<ColourSinglet,double> > strings = ropewalk.getSinglets(stringyCut);
        vector<ColourSinglet> toHadronization;    
        double avge = 0;
        for(int i = 0, N = strings.size(); i < N; ++i ){
          avge += strings[i].second;
          pythia.getRopeUserHooksPtr()->setEnhancement(strings[i].second);
          toHadronization.clear();
          toHadronization.push_back(strings[i].first);  
          hadronizeSystems(pythia,toHadronization,all);
      }
       // ofstream out(( "/scratch/galette/bierlich/work/dipsy/pprange/" + generator()->runName() + "/stringplot.txt").c_str(),ios::app);
       // out << generator()->currentEvent()->weight() << " " << strings.size() << " " << avge << endl;
      break;
     }
      case 5: // Altering kappa per dipole
     {
        if (throwaway) {
          Ropewalk ropewalk_init(singlets, stringR0, stringm0,
				 baryonSuppression, throwaway, false);
	  vector< pair<ColourSinglet,double> > allStrings =
	    ropewalk_init.getSinglets(stringyCut);
	  tPVector allParticles;
	  for ( vector< pair<ColourSinglet,double> >::iterator sItr = allStrings.begin();
		sItr != allStrings.end(); ++sItr )
	    for( tcPVector::iterator pItr = sItr->first.partons().begin();
			pItr != sItr->first.partons().end(); ++pItr )
	      allParticles.push_back(const_ptr_cast<tPPtr>(*pItr));
	  singlets =  theCollapser->collapse(allParticles, newStep());	 
        }
        Ropewalk ropewalk(singlets, stringR0, stringm0, baryonSuppression, false, false);
	typedef multimap<Energy,Ropewalk::DipoleMap::const_iterator> OrderedMap;
	OrderedMap ordered;
        for ( Ropewalk::DipoleMap::iterator itr = ropewalk.begin();
	      itr != ropewalk.end(); ++itr )
	  ordered.insert(make_pair(maxPT(*itr->first, stringyCut), itr));
	for ( OrderedMap::iterator it = ordered.begin(); it != ordered.end(); ++it ) {
	  if ( !pythia.getRopeUserHooksPtr()->setDipoles(&(*it->second), stringm0,
							 stringR0, !average, alpha)) {
	    pythia.getRopeUserHooksPtr()->
	      setEnhancement(ropewalk.averageKappaEnhancement(it->second, stringyCut));
	    for ( int i = 0, N = it->second->second.size(); i < N; ++i )
	      it->second->second[i]->hadr = true;
	  }
	  vector<ColourSinglet> toHadronization(1, *(it->second->first));
    forcerun = false;

	  if(!hadronizeSystems(pythia, toHadronization, all))
	    hadronizeSystems(pythia, toHadronization, all);
	   
	}
	
      break;
     }
      case 99: // New Pythia
      {
      pythia.getRopeUserHooksPtr()->setEnhancement(1.0);
      hadronizeSystems(pythia, singlets, all);
      break;
      }
    default:
      cout << "We should really not be here. This is bad..." << endl;
  }

  //      fScheme = oldFS;
  
    // Do short analysis here:
  if(doAnalysis){
      ofstream out((analysisPath + "analysis.txt").c_str(),ios::app);
        // out << "<event>" << endl;
        // out << PrintStringsToMatlab(singlets) << endl;
        // out << "</event>" << endl;
      out.close();
    }
 
}

void StringFragmentation::dofinish(){
  //if(doAnalysis!=0){
  //  YODA::mkWriter("yoda");
  //  YODA::WriterYODA::write(analysisPath + generator()->runName() + "1D.yoda",_histograms.begin(),_histograms.end());
  //  YODA::WriterYODA::write(analysisPath + generator()->runName() + "2D.yoda",_histograms2D.begin(),_histograms2D.end());

  //}
}

Energy StringFragmentation::maxPT(const ColourSinglet & cs, double deltaY) {
  MaxCmp<Energy2> maxpt2;
  for ( int i = 0, N = cs.partons().size(); i < N; ++i ) {
    if ( deltaY > 0.0 && abs(cs.partons()[i]->eta()) > deltaY ) continue;
    maxpt2(cs.partons()[i]->momentum().perp2());
  }
  if ( !maxpt2 ) return ZERO;
  return sqrt(maxpt2.value());
}

bool StringFragmentation::
hadronizeSystems(Pythia8Interface & pyt, const vector<ColourSinglet> & singlets, const tPVector & all) {

  TheP8EventShapes * _es = NULL;
  Pythia8::Event & event = pyt.event();
  pyt.clearEvent();

  for ( int i = 0, N = singlets.size(); i < N; ++i ) {

    if ( singlets[i].nPieces() == 3 ) {
      // Simple junction.
      // Save place where we will store dummy particle.
      int nsave = event.size();
      event.append(22, -21, 0, 0, 0, 0, 0, 0, 0.0, 0.0, 0.0, 0.0);
      int npart = 0;
      for ( int ip = 1; ip <= 3; ++ip )
        for ( int j = 0, M = singlets[i].piece(ip).size(); j < M; ++j ) {
          pyt.addParticle(singlets[i].piece(ip)[j], 23, nsave, 0);
          ++npart;
        }
      event[nsave].daughters(nsave + 1, nsave + npart);  
    }
    else if ( singlets[i].nPieces() == 5 ) {
      // Double, connected junction
      // Save place where we will store dummy beam particles.
      int nb1 = event.size();
      int nb2 = nb1 + 1;
      event.append(2212, -21, 0, 0, nb1 + 2, nb1 + 4, 0, 0,
                            0.0, 0.0, 0.0, 0.0);
      event.append(-2212, -21, 0, 0, nb2 + 4, nb2 + 6, 0, 0,
                            0.0, 0.0, 0.0, 0.0);

      // Find the string piece connecting the junctions, and the
      // other loose string pieces.
      int connector = 0;
      int q1 = 0;
      int q2 = 0;
      int aq1 = 0;
      int aq2 = 0;
      for ( int ip = 1; ip <= 5; ++ip ) {
        if ( singlets[i].sink(ip).first && singlets[i].source(ip).first )
            connector = ip;
        else if ( singlets[i].source(ip).first ) {
          if ( q1 ) q2 = ip;
          else q1 = ip;
        }
        else if ( singlets[i].sink(ip).first ) {
          if ( aq1 ) aq2 = ip;
          else aq1 = ip;
        }
      }
      if ( !connector || !q1 || !q2 || !aq1 || ! aq2 )
        Throw<StringFragError>()
          << name() << " found complicated junction string. Although Pythia8 can "
          << "hadronize junction strings, this one was too complicated."
          << Exception::runerror;
      
      // Insert the partons of the loose triplet ends.
      int start = event.size();
      for ( int j = 0, M = singlets[i].piece(q1).size(); j < M; ++j )
        pyt.addParticle(singlets[i].piece(q1)[j], 23, nb1, 0);
      for ( int j = 0, M = singlets[i].piece(q2).size(); j < M; ++j )
        pyt.addParticle(singlets[i].piece(q2)[j], 23, nb1, 0);
      // Insert dummy triplet incoming parton with correct colour code.
      int col1 = singlets[i].piece(connector).empty()? event.nextColTag():
        pyt.addColourLine(singlets[i].piece(connector).front()->colourLine());
      int dum1 = event.size();
      event.append(2, -21, nb1, 0, 0, 0, col1, 0, 0.0, 0.0, 0.0, 0.0, 0.0);
      event[nb1].daughters(start, start + singlets[i].piece(q1).size() +
                                    singlets[i].piece(q2).size());

      // Insert the partons of the loose anti-triplet ends.
      start = event.size();
      for ( int j = 0, M = singlets[i].piece(aq1).size(); j < M; ++j )
        pyt.addParticle(singlets[i].piece(aq1)[j], 23, nb2, 0);
      for ( int j = 0, M = singlets[i].piece(aq2).size(); j < M; ++j )
        pyt.addParticle(singlets[i].piece(aq2)[j], 23, nb2, 0);
      // Insert dummy anti-triplet incoming parton with correct colour code.
      int col2 = singlets[i].piece(connector).empty()? col1:
        pyt.addColourLine(singlets[i].piece(connector).back()->antiColourLine());
      int dum2 = event.size();
      event.append(-2, -21, nb2, 0, 0, 0, 0, col2, 0.0, 0.0, 0.0, 0.0, 0.0);
      event[nb2].daughters(start, start + singlets[i].piece(aq1).size() +
                                    singlets[i].piece(aq2).size());

      // Add the partons from the connecting string piece.
      for ( int j = 0, M = singlets[i].piece(connector).size(); j < M; ++j )
        pyt.addParticle(singlets[i].piece(connector)[j], 23, dum1, dum2);
    }
    else if ( singlets[i].nPieces() > 1 ) {
      // We don't know how to handle other junctions yet.
      Throw<StringFragError>()
        << name() << " found complicated junction string. Although Pythia8 can "
        << "hadronize junction strings, that interface is not ready yet."
        << Exception::runerror;
    } else {
      // Normal string
      for ( int j = 0, M = singlets[i].partons().size(); j < M; ++j )
		pyt.addParticle(singlets[i].partons()[j], 23, 0, 0);
    }

  }

  for ( int i = 0, N = all.size(); i < N; ++i )
    if ( !all[i]->coloured() )
      pyt.addParticle(all[i], 23, 0, 0);

  int oldsize = event.size();

  Pythia8::Event saveEvent = event;
  int ntry = maxTries;
  CurrentGenerator::Redirect stdout(cout, false);
  while ( !pyt.go() && --ntry ) event = saveEvent;

  if ( !ntry ) {
    static DebugItem printfailed("TheP8I::PrintFailed", 10);
    if ( printfailed ) {
      double ymax = -1000.0;
      double ymin = 1000.0;
      double sumdy = 0.0;
      for ( int i = 0, N = singlets.size(); i < N; ++i )
        for ( int j = 0, M = singlets[i].partons().size(); j < M; ++j ) {
          const Particle & p = *singlets[i].partons()[j];
          ymax = max(ymax, (_es ? _es->yT(p.momentum()) : p.momentum().rapidity()));
          //ymax = max(ymax, p.momentum().rapidity());
          ymin = min(ymin,(_es ? _es->yT(p.momentum()) : p.momentum().rapidity()));
          //ymin = min(ymin, p.momentum().rapidity());
          cerr << setw(5) << j << setw(14) << p.momentum().rapidity()
               << setw(14) << p.momentum().perp()/GeV
               << setw(14) << p.momentum().m()/GeV;
          if ( j == 0 && p.data().iColour() == PDT::Colour8 ) {
            cerr << setw(14) << (p.momentum() + singlets[i].partons().back()->momentum()).m()/GeV;
            sumdy += abs((_es ? _es->yT(p.momentum()) : p.momentum().rapidity()) ) -(_es ? _es->yT(singlets[i].partons().back()->momentum())
             : singlets[i].partons().back()->momentum().rapidity());
            //sumdy += abs(p.momentum().rapidity() - singlets[i].partons().back()->momentum().rapidity());
          }
          else if ( j > 0 ) {
            cerr << setw(14) << (p.momentum() + singlets[i].partons()[j-1]->momentum()).m()/GeV;
            sumdy += abs((_es ? _es->yT(p.momentum()) : p.momentum().rapidity()) ) -(_es ? _es->yT(singlets[i].partons()[j-i]->momentum())
             : singlets[i].partons()[j-i]->momentum().rapidity());

            //sumdy += abs(p.momentum().rapidity() - singlets[i].partons()[j-1]->momentum().rapidity());
            if ( j + 1 < M )
              cerr << setw(14) << (p.momentum() + singlets[i].partons()[j-1]->momentum()).m()*
                (p.momentum() + singlets[i].partons()[j+1]->momentum()).m()/
                (p.momentum() + singlets[i].partons()[j+1]->momentum() +
                 singlets[i].partons()[j-1]->momentum()).m()/GeV;
          }
          cerr << endl;
        }
      cerr << setw(14) << ymax - ymin << setw(14) << sumdy/(ymax - ymin) << endl;
    }
    CurrentGenerator::Redirect stdout2(cout, true);
    //event.list();
    pyt.errorlist();
    cout << "ThePEG event listing:\n" << *(generator()->currentEvent());
    Throw<StringFragError>()
      << "Pythia8 failed to hadronize partonic state:\n"
      << stdout2.str() << "This event will be discarded!\n"
      << Exception::warning;
    throw Veto();
  }
if(window){
  //	cout << stringyCut << " " << stringpTCut/GeV << endl;
	if(forcerun) forcerun = false;
	else{
	  for( int i = 1, N = event.size(); i < N; ++i ) {
	  	tPPtr p = pyt.getParticle(i);
	//  cout << p->momentum().perp()/GeV << " " << p->rapidity() << " " << p->id() << endl;  	
	    if (p->momentum().perp() > stringpTCut && abs(p->rapidity()) < stringyCut ){  
       /* if(abs(p->id()) == 310 || abs(p->id()) == 3122){
          static ofstream fout(( generator()->runName() + ".txt").c_str(),ios::app);        
          tcPVector pVector = singlets[0].partons();
          fout << p->id() << " " << generator()->currentEvent()->weight() << " ";
          for(size_t j=0;j<pVector.size();++j){
            //if(pVector[j]->id()<20){
              fout << pVector[j]->id() << " " << pVector[j]->momentum().perp()/GeV;
              fout.close();
            //}
          }
          */
        
	  		 forcerun = true;
	  		 return false;
	  		 }
      }
		}
  }

  //  event.list(cerr);
  map<tPPtr, set<tPPtr> > children;
  map<tPPtr, set<tPPtr> > parents;
  for ( int i = 1, N = event.size(); i < N; ++i ) {
    tPPtr p = pyt.getParticle(i);
    int d1 = event[i].daughter1();
    if ( d1 <= 0 ) continue;
    children[p].insert(pyt.getParticle(d1));
    parents[pyt.getParticle(d1)].insert(p);
    int d2 = event[i].daughter2();
    if ( d2 > 0 ) {
      children[p].insert(pyt.getParticle(d2));
      parents[pyt.getParticle(d2)].insert(p);
    }
    for ( int di = d1 + 1; di < d2; ++di ) {
      children[p].insert(pyt.getParticle(di));
      parents[pyt.getParticle(di)].insert(p);
    }
  }

  for ( int i = oldsize, N = event.size(); i < N; ++i ) {
    PPtr p = pyt.getParticle(i);
    set<tPPtr> & pars = parents[p];
    if ( !p ) {
      Throw<StringFragError>()
        << "Failed to reconstruct hadronized state from Pythia8:\n"
        << stdout.str() << "This event will be discarded!\n" << Exception::warning;
      throw Veto();
    }
    if ( isnan(p->momentum().perp()/GeV) ){
            Throw<StringFragError>()
        << "Failed to reconstruct hadronized state from Pythia8:\n"
        << stdout.str() << "This event will be discarded!\n" << Exception::warning;
      throw Veto();
    }
    if ( pars.empty() ) {
      Pythia8::Particle & pyp = event[i];
      if ( pyp.mother1() > 0 ) pars.insert(pyt.getParticle(pyp.mother1()));
      if ( pyp.mother2() > 0 ) pars.insert(pyt.getParticle(pyp.mother2()));
      for ( int im = pyp.mother1() + 1; im < pyp.mother2(); ++im )
        pars.insert(pyt.getParticle(im));
      if ( pars.empty() ) {
        Throw<StringFragError>()
          << "Failed to reconstruct hadronized state from Pythia8:\n"
          << stdout.str() << "This event will be discarded!\n" << Exception::warning;
        throw Veto();
      }
    }
    if ( pars.size() == 1 ) {
      tPPtr par = *pars.begin();
      if ( children[par].size() == 1 &&
           *children[par].begin() == p && par->id() == p->id() )
        newStep()->setCopy(par, p);
      else
        newStep()->addDecayProduct(par, p);
    } else {
      newStep()->addDecayProduct(pars.begin(), pars.end(), p);
    }
  }
  return true;
}

string StringFragmentation::PrintStringsToMatlab(vector<ColourSinglet>& singlets) {
    stringstream ss;
    vector<vector<double> > drawit;
    vector<char> names;
    for(int i=0;i<26;++i) names.push_back(char(i+97));
    vector<char>::iterator nItr = names.begin();

 
    for(vector<ColourSinglet>::iterator sItr = singlets.begin(); sItr!=singlets.end(); ++sItr){
      drawit.clear();
      for (tcPVector::iterator pItr = sItr->partons().begin(); pItr!=sItr->partons().end();++pItr) {
          vector<double> tmp;
          tmp.clear();
          tmp.push_back((*pItr)->rapidity());
          tmp.push_back((*pItr)->vertex().x()/femtometer);
          tmp.push_back((*pItr)->vertex().y()/femtometer);
          drawit.push_back(tmp);
      }
      ss << *nItr << " = [";
        for(int i=0, N=drawit.size();i<N;++i){
          ss << drawit[i].at(0) << ", " << drawit[i].at(1) << ", " << drawit[i].at(2) << ";" << endl;
        }
      ss << "]';\n\n" << endl;
      ++nItr;
    }
    return ss.str();
}


IBPtr StringFragmentation::clone() const {
  return new_ptr(*this);
}

IBPtr StringFragmentation::fullclone() const {
  return new_ptr(*this);
}


void StringFragmentation::doinitrun() {
  HadronizationHandler::doinitrun();

   theAdditionalP8Settings.push_back("ProcessLevel:all = off");
   theAdditionalP8Settings.push_back("HadronLevel:Decay = off");
   theAdditionalP8Settings.push_back("Check:event = off");
   theAdditionalP8Settings.push_back("Next:numberCount = 0");
   theAdditionalP8Settings.push_back("Next:numberShowLHA = 0");
   theAdditionalP8Settings.push_back("Next:numberShowInfo = 0");
   theAdditionalP8Settings.push_back("Next:numberShowProcess = 0");
   theAdditionalP8Settings.push_back("Next:numberShowEvent = 0");
   theAdditionalP8Settings.push_back("Init:showChangedSettings = 0");
   theAdditionalP8Settings.push_back("Init:showAllSettings = 0");
   theAdditionalP8Settings.push_back("Init:showChangedParticleData = 0");
   theAdditionalP8Settings.push_back("Init:showChangedResonanceData = 0");
   theAdditionalP8Settings.push_back("Init:showAllParticleData = 0");
   theAdditionalP8Settings.push_back("Init:showOneParticleData = 0");
   theAdditionalP8Settings.push_back("Init:showProcesses = 0");

  if(fScheme == 4 || fScheme == 5 || fScheme == 1){ // We don't need multiple Pythia objects anymore!
    pythia.enableHooks();
    pythia.init(*this,theAdditionalP8Settings);
    if(pythia.version() > 0 && pythia.version() - 1.234 > 1.0){
       cout << "Pythia version: " <<  pythia.version() << endl;
       cout << "The chosen fragmentation scheeme requires a tweaked "
	    << "Pythia v. 1.234 for option. I will default you to "
	 "option 'pythia'." << endl;
      fScheme = 0;
    }
    else{
      PytPars p;
  
     p.insert(make_pair<const string,double>("StringPT:sigma",theStringPT_sigma));
     p.insert(make_pair<const string,double>("StringZ:aLund",theStringZ_aLund));
     p.insert(make_pair<const string,double>("StringZ:bLund",theStringZ_bLund));
     p.insert(make_pair<const string,double>("StringFlav:probStoUD",theStringFlav_probStoUD));
     p.insert(make_pair<const string,double>("StringFlav:probSQtoQQ",theStringFlav_probSQtoQQ));
     p.insert(make_pair<const string,double>("StringFlav:probQQ1toQQ0",theStringFlav_probQQ1toQQ0));
     p.insert(make_pair<const string,double>("StringFlav:probQQtoQ",theStringFlav_probQQtoQ));

     phandler.init(fragmentationMass*fragmentationMass,baryonSuppression,p);
     pythia.getRopeUserHooksPtr()->setParameterHandler(&phandler);
     //pythia.getRopeUserHooksPtr()->setWindow(window,stringpTCut);
    }
  }
  if( !( fScheme == 4  || fScheme == 5 || fScheme == 1) ){ // Not 'else' as it can change above...

   vector<string> moresettings = theAdditionalP8Settings;
   moresettings.push_back("StringPT:sigma = " + convert(theStringPT_sigma));
   moresettings.push_back("StringZ:aLund = " + convert(theStringZ_aLund));
   moresettings.push_back("StringZ:bLund = " + convert(theStringZ_bLund));
   moresettings.push_back("StringFlav:probStoUD = " +
  			 convert(theStringFlav_probStoUD));
   moresettings.push_back("StringFlav:probSQtoQQ = " +
  			 convert(theStringFlav_probSQtoQQ));
   moresettings.push_back("StringFlav:probQQ1toQQ0 = " +
  			 convert(theStringFlav_probQQ1toQQ0));
   moresettings.push_back("StringFlav:probQQtoQ = " +
  			 convert(theStringFlav_probQQtoQ));
   moresettings.push_back("OverlapStrings:fragMass = " +
  			 convert(fragmentationMass));
   moresettings.push_back("OverlapStrings:baryonSuppression = " + convert(baryonSuppression));

  // Initialize the OverlapPythia Handler
  if ( opHandler ) delete opHandler;
  opHandler = new OverlapPythiaHandler(this,moresettings);
  
  }

  // Should we do event shapes?
  if(_eventShapes) delete _eventShapes;
  _eventShapes = NULL;
  if(useThrustAxis==1){
    _eventShapes = new TheP8EventShapes();
  }

  //if( doAnalysis ) {
  // _histograms.push_back(new YODA::Histo1D(100,1,15,"/h_lowpT","h_lowpT"));
  //_histograms.push_back(new YODA::Histo1D(100,1,15,"/h_midpT","h_midpT"));
  //  _histograms.push_back(new YODA::Histo1D(100,1,15,"/h_highpT","h_highpT"));
  //  _histograms2D.push_back(new YODA::Histo2D(10,0.,15,10,0,15,"/m_vs_pT","m_vs_pT"));
  //}
    nev = 0;
    
}

void StringFragmentation::persistentOutput(PersistentOStream & os) const {
  os
#include "StringFragmentation-output.h"
    << fScheme << ounit(stringR0,femtometer) << ounit(stringm0,GeV) << junctionDiquark << alpha << average << stringyCut
    << ounit(stringpTCut,GeV) << fragmentationMass << baryonSuppression
    << oenum(window) << oenum(throwaway) << useThrustAxis << maxTries
    << theCollapser << doAnalysis << analysisPath;

}

void StringFragmentation::persistentInput(PersistentIStream & is, int) {
  is
#include "StringFragmentation-input.h"
    >> fScheme >> iunit(stringR0,femtometer) >> iunit(stringm0,GeV) >> junctionDiquark >> alpha >> average >> stringyCut
    >> iunit(stringpTCut,GeV) >> fragmentationMass >> baryonSuppression
    >> ienum(window) >> ienum(throwaway) >> useThrustAxis >> maxTries
    >> theCollapser >> doAnalysis >> analysisPath;
}

ClassDescription<StringFragmentation> StringFragmentation::initStringFragmentation;
// Definition of the static class description member.

void StringFragmentation::Init() {
#include "StringFragmentation-interfaces.h"
  static Reference<StringFragmentation,ClusterCollapser> interfaceCollapser
    ("Collapser",
     "A ThePEG::ClusterCollapser object used to collapse colour singlet "
     "clusters which are too small to fragment. If no object is given the "
     "MinistringFragmentetion of Pythia8 is used instead.",
     &StringFragmentation::theCollapser, true, false, true, true, false);
  
  static Switch<StringFragmentation,int> interfaceFragmentationScheme
    ("FragmentationScheme",
     "Different options for how to handle overlapping strings.",
     &StringFragmentation::fScheme, 0, true, false);
  static SwitchOption interfaceFragmentationSchemepythia
    (interfaceFragmentationScheme,
     "pythia",
     "Plain old Pythia fragmentation",
     0);
  static SwitchOption interfaceFragmentationSchemedep1
    (interfaceFragmentationScheme,
     "dep1",
     "Not sure about this.",
     1);
  static SwitchOption interfaceFragmentationSchemedep2
    (interfaceFragmentationScheme,
     "dep2",
     "Not sure about this.",
     2);
  static SwitchOption interfaceFragmentationSchemenone
    (interfaceFragmentationScheme,
     "none",
     "Plain old Pythia fragmentation",
     0);
  static SwitchOption interfaceFragmentationSchemepipe
    (interfaceFragmentationScheme,
     "pipe",
     "Each string is enclosed by a cylinder. Effective parameters are calculated from the overlap of cylinders.",
     3);
  static SwitchOption interfaceFragmentationSchemeaverage
    (interfaceFragmentationScheme,
     "average",
     "The overlap is calculated for each dipole in all strings, and for each string the effective parameters are obtained from the average overlap.",
     4);
  static SwitchOption interfaceFragmentationSchemedipole
    (interfaceFragmentationScheme,
     "dipole",
     "Effective parameters are calculated for each breakup by determining the overlap of the corresponding dipole with other dipoles.",
     5);

  static Parameter<StringFragmentation,int> interfaceUseThrustAxis
    ("UseThrustAxis",
     "Whether or not to use rapidity wrt. Thrust Axis when counting strings "
     "(put 1 or 0).",
     &StringFragmentation::useThrustAxis, 0, 0, 0,
     true, false, Interface::lowerlim);

  static Parameter<StringFragmentation,int> interfaceDoAnalysis
    ("DoAnalysis",
     "Can do a short analysis where histograms of the effective parameters are "
     "plotted, and the strings in (bx,by,y)-space are printed for a couple of events. "
     "(put 1 or 0).",
     &StringFragmentation::doAnalysis, 0, 0, 0,
     true, false, Interface::lowerlim);

  static Parameter<StringFragmentation,string> interfaceAnalysisPath
    ("AnalysisPath",
     "Set the system path where the (optional) analysis output should appear "
     " (directory must exist!)" ,
     &StringFragmentation::analysisPath, "");

  static Switch<StringFragmentation,bool> interfaceThrowAway
    ("ThrowAway",
     "Throw away strings with weight:"
     "Use 1 - (p + q)/(m + n) for (0,-1) configuration." ,
     &StringFragmentation::throwaway, false, true, false);
  static SwitchOption interfaceThrowAwayTrue
    (interfaceThrowAway,"True","enabled.",true);
  static SwitchOption interfaceThrowAwayFalse
    (interfaceThrowAway,"False","disabled.",false);

  static Switch<StringFragmentation,bool> interfaceWindow
    ("StringWindow",
     "Enable the 'window'-cut procedure, off by default."
     "Parameters will only affect if this is switched on. " ,
     &StringFragmentation::window, false, true, false);
  static SwitchOption interfaceWindowTrue
    (interfaceWindow,"True","enabled.",true);
  static SwitchOption interfaceWindowFalse
    (interfaceWindow,"False","disabled.",false);

  static Parameter<StringFragmentation,Energy> interfaceStringpTCut
    ("StringpTCut",
     "No enhancement of strings with a constituent pT higher than StringpTCut (GeV), and within "
     " StringyCut (in GeV).",
     &StringFragmentation::stringpTCut, GeV, 6.0*GeV,
     0.0*GeV, 0*GeV,
     true, false, Interface::lowerlim);

  static Parameter<StringFragmentation,Energy> interfaceStringm0
    ("Stringm0",
     ".",
     &StringFragmentation::stringm0, GeV, 1.0*GeV,
     0.0*GeV, 0*GeV,
     true, false, Interface::lowerlim);

  static Parameter<StringFragmentation,double> interfaceJunctionDiquark
    ("JunctionDiquark",
     "Suppress diquark production in breaking of junctions with this amount",
     &StringFragmentation::junctionDiquark, 0, 1.0,
     0.0, 0,
     true, false, Interface::lowerlim);

  static Parameter<StringFragmentation,double> interfaceAlpha
    ("Alpha",
     "Boost the enhanced string tension with additional factor",
     &StringFragmentation::alpha, 0, 1.0,
     1.0, 0,
     true, false, Interface::lowerlim);

  static Switch<StringFragmentation,bool> interfaceAverage
    ("Average",
     "Disable to try and hadronise strings with individual tensions"
     "instead of average value." ,
     &StringFragmentation::average, false, true, false);
  static SwitchOption interfaceAverageTrue
    (interfaceAverage,"True","enabled.",true);
  static SwitchOption interfaceAverageFalse
    (interfaceAverage,"False","disabled.",false);

  static Parameter<StringFragmentation,double> interfaceStringyCut
    ("StringyCut",
      "No enhancement of strings with a constituent pT higher than StringpTCut (GeV), and within "
     " StringyCut.",
       &StringFragmentation::stringyCut, 0, 2.0,
     0.0, 0,
     true, false, Interface::lowerlim);

  static Parameter<StringFragmentation,Length> interfaceStringR0
    ("StringR0",
     "In the string overlap model the string R0 dictates the minimum radius "
     " of a string (in fm).",
     &StringFragmentation::stringR0, femtometer, 0.5*femtometer,
     0.0*femtometer, 0*femtometer,
     true, false, Interface::lowerlim);

  static Parameter<StringFragmentation,double> interfaceFragmentationMass
    ("FragmentationMass",
     "Set the mass used for the f(z) integral in the overlap string model "
     " default is pion mass (units GeV).",
     &StringFragmentation::fragmentationMass, 0, 0.135,
     0., 1.,
     true, false, Interface::lowerlim);

  static Parameter<StringFragmentation,double> interfaceBaryonSuppression
    ("BaryonSuppression",
     "Set the fudge factor used for baryon suppression in the overlap string model. "
     " One day this should be calculated properly. This day is not today. Probably around 2...",
     &StringFragmentation::baryonSuppression, 0, 0.5,
     0.0, 1.0,
     true, false, Interface::limited);

  static Parameter<StringFragmentation,int> interfaceMaxTries
    ("MaxTries",
     "Sometimes Pythia gives up on an event too easily. We therefore allow "
     "it to re-try a couple of times.",
     &StringFragmentation::maxTries, 2, 1, 0,
     true, false, Interface::lowerlim);

}

string StringFragmentation::convert(double d) {
  ostringstream os;
  os << d;
  return os.str();
}
