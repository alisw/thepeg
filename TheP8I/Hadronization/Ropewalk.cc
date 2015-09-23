// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Ropewalk class.
//

#include "Ropewalk.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/EventRecord/Step.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/PDT/EnumParticles.h"
#include "ThePEG/Utilities/Selector.h"
#include "ThePEG/Utilities/SimplePhaseSpace.h"
#include "ThePEG/Utilities/UtilityBase.h"
#include "ThePEG/Utilities/Throw.h"
#include "ThePEG/Handlers/StepHandler.h"



using namespace TheP8I;

Ropewalk::Ropewalk(const vector<ColourSinglet> & singlets, Length width,
		   Energy scale, double jDP, bool throwaway, bool rap, bool deb) 
  : R0(width), m0(scale), junctionDiquarkProb(jDP), rapidityOverlap(rap), debug(deb),
    strl0(0.0), strl(0.0), avh(0.0), avp(0.0), avq(0.0) {
  for ( int is = 0, Ns = singlets.size(); is < Ns; ++is ){
    stringdipoles[cloneToFinal(singlets[is])];
  }
  getDipoles();

  lambdaBefore = lambdaSum();
  if ( debug ) CurrentGenerator::log() << singlets.size() << "s "
				       << dipoles.size() << "d ";
  calculateOverlaps();
  if(throwaway) doBreakups();
  
}

Ropewalk::~Ropewalk() {
  for ( DipoleMap::iterator it = stringdipoles.begin();
	it != stringdipoles.end(); ++it ) delete it->first;
}

ColourSinglet * Ropewalk::cloneToFinal(const ColourSinglet & cs) {
  tcParticleSet final;
  for ( int i = 0, N = cs.partons().size(); i < N; ++i )
    final.insert(cs.partons()[i]->final());
  tcPPtr p = *final.begin();
  if ( p->colourLine() ) return new ColourSinglet(p->colourLine(), final);
  if ( p->antiColourLine() ) return new ColourSinglet(p->antiColourLine(), final);
  Throw<ColourException>()
    << "Cloning ColourSinglets failed in Ropewalk. "
    << "This is a serious error - please contact the authors."
    << Exception::abortnow;
  return 0;
}

LorentzPoint Ropewalk::
propagate(const LorentzPoint & b, const LorentzMomentum & p) const {
  return b + p/(p.e()*m0/hbarc);
  //  return b;
}

void Ropewalk::getDipoles() {

  // First extract all dipoles from all strings (pieces).
  for ( DipoleMap::iterator it = stringdipoles.begin();
	it != stringdipoles.end(); ++it ) {
    ColourSinglet & string = *it->first;
    for ( int isp = 1, Nsp = string.nPieces(); isp <= Nsp; ++isp ) {
      const ColourSinglet::StringPiece & sp = string.piece(isp);
      if ( sp.size() < 2 ) continue;
      bool forward = sp[0]->colourLine() &&
	sp[0]->colourLine() == sp[1]->antiColourLine();
      for ( int i = 0, N = sp.size(); i < (N - 1); ++i ) {
	if ( forward )
	  dipoles.push_back(Dipole(sp[i], sp[i + 1], string));
	else
	  dipoles.push_back(Dipole(sp[i + 1], sp[i], string));
      }
      if ( !forward && sp.front()->colourLine() &&
	   sp.front()->colourLine() == sp.back()->antiColourLine() )
	dipoles.push_back(Dipole(sp.front(), sp.back(), string));
      else if ( forward && sp.front()->antiColourLine() &&
		sp.front()->antiColourLine() == sp.back()->colourLine() )
	dipoles.push_back(Dipole(sp.back(), sp.front(), string));
    }
  }

  for ( int i = 0, N = dipoles.size(); i < N; ++i ) {
    stringdipoles[dipoles[i].string].push_back(&dipoles[i]);
    strl0 += dipoles[i].yspan(1.0*GeV);
  }

}

double Ropewalk::limitrap(const LorentzMomentum & p) const {
  if ( p.z() == ZERO ) return 0.0;
  Energy mt = sqrt(max(p.perp2() + p.m2(), sqr(m0)));
  double rap = log((p.t() + abs(p.z()))/mt);
  return p.z() > ZERO? rap: -rap;
} 

Ropewalk::OverlappingDipole::
OverlappingDipole(const Dipole & d, const LorentzRotation R, const Ropewalk * rw)
  : dipole(&d), dir(1) {
  // Get the boost to other dipole's rest frame and calculate
  // rapidities.
  bc = R*rw->propagate(d.pc->vertex(), d.pc->momentum());
  ba = R*rw->propagate(d.pa->vertex(), d.pa->momentum());
  yc = rw->limitrap(R*d.pc->momentum());
  ya = rw->limitrap(R*d.pa->momentum());
  if ( yc < ya ) dir = -1;
}

void Ropewalk::calculateOverlaps() {

  // Then calculate all overlaps between pairs of dipoles.
  for ( int i1 = 0, Nd = dipoles.size(); i1 < Nd; ++i1 ) {
    Dipole & d1 = dipoles[i1];
    if ( d1.s() <= sqr(m0) ) continue;
    // Get the boost to dipoles rest frame and calculate rapidities.
    LorentzMomentum pc1 = d1.pc->momentum();
    LorentzMomentum pa1 = d1.pa->momentum();
    LorentzRotation R = Utilities::boostToCM(make_pair(&pc1, &pa1));
    d1.bc = R*propagate(d1.pc->vertex(), d1.pc->momentum());
    d1.ba = R*propagate(d1.pa->vertex(), d1.pa->momentum());
    double yc1 = limitrap(pc1);
    double ya1 = limitrap(pa1);
    double n = 0.0;
    double m = 0.0;
    if ( yc1 <= ya1 ) continue; // don't care about too small strings.
    for ( int i2 = 0; i2 < Nd; ++i2 ) {
      if ( i1 == i2 ) continue;
      // Boost to the rest frame of dipole 1.
      if ( dipoles[i2].s() <= sqr(m0) ) continue;
      OverlappingDipole od(dipoles[i2], R, this);
      // Ignore if not overlapping in rapidity.
      if(rapidityOverlap)
	if ( min(od.yc, od.ya) > yc1 || max(od.yc, od.ya) < ya1 || od.yc == od.ya )
	continue;
      // Calculate rapidity overlap.
      double yomax = min(max(od.yc, od.ya), yc1);
      double yomin = max(min(od.yc, od.ya), ya1);
      // Sample a random point in the rapidity overlap region and get
      // positions in space and check overlap.
      double yfrac = (yomax - yomin)/(yc1 - ya1);
      double y = UseRandom::rnd(yomin, yomax);
      if ( !od.overlap(y, d1.ba + (d1.bc - d1.ba)*(y - ya1)/(yc1 - ya1), R0) )
	yfrac = 0.0;
      if(!rapidityOverlap)
	yfrac = 1.0;
      // Sum overlaps.
      if ( od.dir > 0 ) {
	n += yfrac;
      } else {
	m += yfrac;
      }
      od.yfrac = yfrac*od.dir;
      d1.overlaps.push_back(od);
    }
    d1.n = int(n + UseRandom::rnd());
    d1.m = int(m + UseRandom::rnd());
    d1.initWeightMultiplet();
    avp += d1.p*d1.yspan(1.0*GeV);
    avq += d1.q*d1.yspan(1.0*GeV);
  } 
}

double Ropewalk::Dipole::reinit(double ry, Length R0, Energy m0) {
  // First calculate the rapidity in this restframe.
  double y = yspan(m0)*(ry - 0.5);
  // Then get the position in space of the string at this rapidity.
  LorentzPoint b = ba + (bc - ba)*ry;
  p = 1;
  q = 0;
  // Now go through all overlapping dipoles.
  for ( int i = 0, N = overlaps.size(); i < N; ++i ) {
    if ( overlaps[i].dipole->broken ) continue;
    if( overlaps[i].dipole->hadr ) continue;
    if ( overlaps[i].overlap(y, b, R0) ) {
      if ( overlaps[i].dir > 0 ) ++p;
      else ++q;
    }
  }
  return kappaEnhancement();
} 

void Ropewalk::Dipole::initMultiplet() {
  typedef pair<int,int> Plet;

  int ns = 0;
  int ms = 0;
  while ( ns < n || ms < m ) {
    Plet diff;
    double octetveto = double(ns + ms + 1 - p - q)/double(ns + ms + 1);
    if ( double(n - ns) > UseRandom::rnd()*double(n + m - ns - ms) ) {
      Selector<Plet> sel;
      sel.insert(multiplicity(p + 1, q), Plet(1, 0));
      sel.insert(multiplicity(p, q - 1)*octetveto, Plet(0, -1));
      sel.insert(multiplicity(p - 1, q + 1), Plet(-1, 1));
      diff = sel[UseRandom::rnd()];
      ++ns;
    } else {
      Selector<Plet> sel;
      sel.insert(multiplicity(p, q + 1), Plet(0, 1));
      sel.insert(multiplicity(p - 1, q)*octetveto, Plet(-1, 0));
      sel.insert(multiplicity(p + 1, q - 1), Plet(1, -1));
      diff = sel[UseRandom::rnd()];
      ++ms;
    }
    if ( diff.first == -diff.second ) nb++;
    p += diff.first;
    q += diff.second;
  }
  /* *** ATTENTION *** maybe better if Christian implements this */

}

void Ropewalk::Dipole::initNewMultiplet() {
  typedef pair<int,int> Plet;

  int ns = 0;
  int ms = 0;
  for ( int i = 0, N = overlaps.size(); i < N; ++i ) {
    if ( abs(overlaps[i].yfrac) <  UseRandom::rnd() ) continue;
    Plet diff(0,0);
    double octetveto = double(ns + ms + 1 - p - q)/double(ns + ms + 1);
    if ( overlaps[i].dir > 0 ) {
      Selector<Plet> sel;
      sel.insert(multiplicity(p + 1, q), Plet(1, 0));
      sel.insert(multiplicity(p, q - 1)*octetveto, Plet(0, -1));
      sel.insert(multiplicity(p - 1, q + 1), Plet(-1, 1));
      diff = sel[UseRandom::rnd()];
      ++ns;
    } else {
      Selector<Plet> sel;
      sel.insert(multiplicity(p, q + 1), Plet(0, 1));
      sel.insert(multiplicity(p - 1, q)*octetveto, Plet(-1, 0));
      sel.insert(multiplicity(p + 1, q - 1), Plet(1, -1));
      diff = sel[UseRandom::rnd()];
      ++ms;
    }
    if ( diff.first == -diff.second ) nb++;
    p += diff.first;
    q += diff.second;
  }
  n = ns;
  m = ms;
  /* *** ATTENTION *** maybe better if Christian implements this */

}

void Ropewalk::Dipole::initWeightMultiplet() {
  typedef pair<int,int> Plet;

  int ns = 0;
  int ms = 0;
  double mab = sqrt(s())/GeV;
  for ( int i = 0, N = overlaps.size(); i < N; ++i ) {
    if ( abs(overlaps[i].yfrac) <  UseRandom::rnd() ) continue;
    const Dipole * od = overlaps[i].dipole;
    double mcd = sqrt(od->s())/GeV;
    double w8 = 1.0/sqr(mab*mcd); //The weight for 3 + 3bar -> 8
    double w6 = w8; // The weight for 3 + 3 -> 6
    Plet diff(0,0);
    double octetveto = double(ns + ms + 1 - p - q)/double(ns + ms + 1);
    if ( overlaps[i].dir > 0 ) {
      double mac = (pc->momentum() + od->pc->momentum()).m()/GeV;
      double mbd = (pa->momentum() + od->pa->momentum()).m()/GeV;
      double w1 = octetveto/sqr(mac*mbd); // The weight for 3 + 3bar -> 1
      double w3 = 1.0/(mab*mcd*mac*mbd); //The weight for 3 + 3 -> 3bar
      Selector<Plet> sel;
      sel.insert(multiplicity(p + 1, q) +
		 multiplicity(p, q - 1)*w8/(w8 + w1) +
		 multiplicity(p - 1, q + 1)*w6/(w6 + w3), Plet(1, 0));
      sel.insert(multiplicity(p, q - 1)*w1/(w1 + w8), Plet(0, -1));
      sel.insert(multiplicity(p - 1, q + 1)*w3/(w3 + w6), Plet(-1, 1));
      diff = sel[UseRandom::rnd()];
      ++ns;
    } else {
      double mac = (pa->momentum() + od->pc->momentum()).m()/GeV;
      double mbd = (pc->momentum() + od->pa->momentum()).m()/GeV;
      if ( mac*mbd <= 0.0 ) { 
	q++;
	continue;
      }
      double w1 = octetveto/sqr(mac*mbd); // The weight for 3 + 3bar -> 1
      double w3 = 1.0/(mab*mcd*mac*mbd); //The weight for 3 + 3 -> 3bar
      Selector<Plet> sel;
      sel.insert(multiplicity(p, q + 1) +
		 multiplicity(p - 1, q)*w8/(w8 + w1) +
		 multiplicity(p + 1, q - 1)*w6/(w6 + w3), Plet(0, 1));
      sel.insert(multiplicity(p - 1, q)*w1/(w1 + w8), Plet(-1, 0));
      sel.insert(multiplicity(p + 1, q - 1)*w3/(w3 + w6), Plet(1, -1));
      diff = sel[UseRandom::rnd()];
      ++ms;
    }
    if ( diff.first == -diff.second ) nb++;
    p += diff.first;
    q += diff.second;
  }
  n = ns;
  m = ms;
  /* *** ATTENTION *** maybe better if Christian implements this */

}

double Ropewalk::getNb(){
  double nba = 0;
   for ( int i = 0, N = dipoles.size(); i < N; ++i ){
     nba += double(dipoles[i].nb);
     //cout << dipoles[i].n << " " << dipoles[i].m << " " << dipoles[i].p << " " << dipoles[i].q << endl;
     //cout << double(dipoles[i].n + dipoles[i].m + 1 - dipoles[i].p - dipoles[i].q) << endl;
     //cout << " --- " << endl;  
     //nba += double(dipoles[i].n + dipoles[i].m + 1 - dipoles[i].p - dipoles[i].q);
  }
   // cout << nba << endl;
   return nba;
}

double Ropewalk::getkW(){
    double kw = 0;
     for ( int i = 0, N = dipoles.size(); i < N; ++i )
       kw += dipoles[i].kappaEnhancement()*log(dipoles[i].s()/sqr(m0));
     return kw;
}

pair<double,double> Ropewalk::getMN(){
  pair<double,double> ret = make_pair<double,double>(0,0);
  for (int i = 0, N = dipoles.size(); i < N; ++i){
    ret.first += dipoles[i].m;
    ret.second += dipoles[i].n;
  }
  return ret;
}

double Ropewalk::lambdaSum(){
  double ret = 0;
  for (int i = 0, N = dipoles.size(); i < N; ++i)
    ret += dipoles[i].s() > 0.1*GeV2 ? log(dipoles[i].s()/sqr(m0)) : 0;
  return ret;
}


double Ropewalk::Dipole::breakupProbability() const {
  if ( !breakable() || n + m <= 0.0 || n + m + 1 == p + q ) return 0.0;
  double sum = 0.0;
  for (int j = 0, N = overlaps.size(); j < N; ++j )
    if ( overlaps[j].dipole->breakable() )
      sum += abs(overlaps[j].yfrac);
  if ( sum <= 0.0 ) return 1.0;
  return double(n + m + 1 - p - q)/(sum + 1.0);
}

Step & Ropewalk::step() const {
  return *(CurrentGenerator()->currentStepHandler()->newStep());
}

void Ropewalk::doBreakups() {

  typedef multimap<double, const Dipole *> BMap;
  BMap breakups;

  // First order alldipoles in increasing probability of breaking
  for ( int i = 0, N = dipoles.size(); i < N; ++i )
    breakups.insert(make_pair(dipoles[i].breakupProbability(), &dipoles[i]));

  // Then start breaking in reverse order
  mspec.clear();
  for ( BMap::reverse_iterator it = breakups.rbegin();
	it != breakups.rend(); ++it ) {
    const Dipole & d = *it->second;
    if ( d.breakupProbability() < UseRandom::rnd() ) continue;
    mspec.push_back(d.s()/GeV2);
    double hinv = 1.0 / d.kappaEnhancement();  
    // The default parameters should be set                                                                                                                                        
    // properly from the outside. Now just Pythia8 defaults                                                                                                                        
    double rho = pow(0.19, hinv);
    double x = pow(1.0,hinv);
    double y = pow(0.027,hinv);

    Selector<int> qsel;
    qsel.clear();
    if(UseRandom::rnd() < junctionDiquarkProb){  
    // Take just ordinary quarks
      qsel.insert(1.0, ParticleID::ubar);
      qsel.insert(1.0, ParticleID::dbar);
      qsel.insert(rho, ParticleID::sbar);
    }
    else {
      // Take diquarks
      // TODO (Pythia 8.2) Set status code of these diquarks to 74
      qsel.insert(1.0, ParticleID::ud_0);
      qsel.insert(3.0*y, ParticleID::dd_1);
      qsel.insert(3.0*y, ParticleID::ud_1);
      qsel.insert(3.0*y, ParticleID::uu_1);
      qsel.insert(rho*x, ParticleID::su_0);
      qsel.insert(3*x*rho*y, ParticleID::su_1);
      qsel.insert(3*y*x*x*rho*rho, ParticleID::ss_1);
      qsel.insert(rho*x, ParticleID::sd_0);
      qsel.insert(3*x*rho*y, ParticleID::sd_1);
    }
    int flav = qsel[UseRandom::rnd()];
    PPtr dic = CurrentGenerator()->getParticle(-flav);
    PPtr dia = CurrentGenerator()->getParticle(flav);
    //cout << (dic->momentum() + dia->momentum()).m2() / GeV2 << " " <<  d.s() / GeV2 << endl;
    if((dic->momentum() + dia->momentum()).m2() > d.s()) continue;
    // Get boost to dipole rest frame and place the di-quarks there.
    LorentzRotation R = Utilities::getBoostFromCM(make_pair(d.pc, d.pa));

    SimplePhaseSpace::CMS(dic, dia, d.s(), 1.0, 0.0);
    dia->transform(R);
    dic->transform(R);
    // Add the di-quarks as children to the decayed gluons and split
    // the resulting string.
    if ( !( step().addDecayProduct(d.pc, dic, true) &&
	    step().addDecayProduct(d.pa, dia, true) &&
	    step().addDecayProduct(d.pc, dia, false) &&
	    step().addDecayProduct(d.pa, dic, false) ) )
      Throw<ColourException>()
	<< "Adding decay products failed in Ropewalk. "
	<< "This is a serious error - please contact the authors."
	<< Exception::abortnow;
    tcParticleSet pset(d.string->partons().begin(), d.string->partons().end());
    pset.erase(d.pc);
    pset.erase(d.pa);
    pset.insert(dic);
    pset.insert(dia);
    // Reset the neighboring dipoles to point to the new diquarks.
    d.broken = true;
    ColourSinglet * olds = d.string;
    const vector<Dipole *> & oldd = stringdipoles[olds];
    for ( int id = 0, Nd = oldd.size(); id < Nd; ++id ) {
      if ( oldd[id]->pc == d.pa ) oldd[id]->pc = dia;
      if ( oldd[id]->pa == d.pc ) oldd[id]->pa = dic;
    }
    // remove the old string and insert the new ones fixing the
    // pointers from the dipolmes.
    vector<ColourSinglet> newstrings =
      ColourSinglet::getSinglets(pset.begin(), pset.end());
    for ( int is = 0, Ns = newstrings.size(); is < Ns; ++is ) {
      ColourSinglet * news = new ColourSinglet(newstrings[is]);
      tcParticleSet pset(news->partons().begin(), news->partons().end());
      vector<Dipole *> & newd = stringdipoles[news];
      for ( int id = 0, Nd = oldd.size(); id < Nd; ++id ) {
	if ( !oldd[id]->broken && pset.find(oldd[id]->pc) != pset.end() ) {
	  newd.push_back(oldd[id]);
	  oldd[id]->string = news;
	}
      }
      sort(newd, news->partons()[0]);
    }
    stringdipoles.erase(olds);
    delete olds;
  }
}

vector< pair<ColourSinglet,double> > Ropewalk::getSinglets(double DeltaY) const {
  vector< pair<ColourSinglet,double> > ret;
  ret.reserve(stringdipoles.size());
  double toty = 0.0;
  double totyh = 0.0;
  double toth = 0.0;
  double toth2 = 0.0;
  double minh = 10.0;
  double maxh = 0.0;
  // Calculate the kappa-enhancement of the remaining strings.
  for ( DipoleMap::const_iterator it = stringdipoles.begin();
	it != stringdipoles.end(); ++it ) {
    double sumyh = 0.0;
    double sumy = 0.0;
    for ( int id = 0, Nd = it->second.size(); id < Nd; ++id ) {
      double y = it->second[id]->yspan(m0);
      double h = it->second[id]->kappaEnhancement();
      avh += h*it->second[id]->yspan(1.0*GeV);
      strl += it->second[id]->yspan(1.0*GeV);
      if ( DeltaY > 0.0 &&
	   abs(it->second[id]->pc->eta()) > DeltaY &&
	   abs(it->second[id]->pa->eta()) > DeltaY ) continue;
      sumy += y;
      sumyh += y*h;
    }
    toty += sumy;
    totyh += sumyh;
    
    if ( sumy > 0.0 ) {
      double h = 0.1*int(10.0*sumyh/sumy + UseRandom::rnd());
      ret.push_back(make_pair(*it->first, h));
      toth += sumyh/sumy;
      toth2 += sqr(sumyh/sumy);
      minh = min(minh, sumyh/sumy);
      maxh = max(maxh, sumyh/sumy);
    } else {
      ret.push_back(make_pair(*it->first, 1.0));
      toth += 1.0;
      toth2 += 1.0;
      minh = min(minh, 1.0);
      maxh = max(maxh, 1.0);
    }
  }
  if ( debug )
    CurrentGenerator::log() << ret.size() << "ns "
			    << totyh/max(toty, 0.00001) << "hwa "
			    << minh << "<h "
			    << toth/ret.size() << "ha "
			    << maxh << ">h "
			    << sqrt(toth2/ret.size() - sqr(toth/ret.size())) << "hsd " << endl;

  if ( avh > 0.0 ) avh /= strl;
  if ( avp > 0.0 ) avp /= strl0;
  if ( avq > 0.0 ) avq /= strl0;

  return ret;

}

double Ropewalk::averageKappaEnhancement(const vector<Dipole *> & dipoles,
					 double DeltaY) const {
    double sumyh = 0.0;
    double sumy = 0.0;
    for ( int id = 0, Nd = dipoles.size(); id < Nd; ++id ) {
      dipoles[id]->reinit(UseRandom::rnd(), R0, m0);
      double y = dipoles[id]->yspan(m0);
      double h = dipoles[id]->firstKappa();
      if ( DeltaY > 0.0 &&
	   abs(dipoles[id]->pc->eta()) > DeltaY &&
	   abs(dipoles[id]->pa->eta()) > DeltaY ) continue;
      sumy += y;
      sumyh += y*h;
    }
    if ( sumy <= 0.0 ) return 1.0;
    return sumyh/sumy;
}

void Ropewalk::sort(vector<Dipole *> & dips, tcPPtr p) {
  if ( p->colourLine() ) {
    map<tcPPtr, Dipole *> dmap;
    Dipole * current = 0;
    for ( int i = 0, N = dips.size(); i < N; ++i ) {
      if ( dips[i]->pa == p ) current =  dips[i];
      else dmap[dips[i]->pa] = dips[i];
    }
    dips.clear();
    while ( current ) {
      dips.push_back(current);
      map<tcPPtr, Dipole *>::iterator it = dmap.find(current->pc);
      if ( it == dmap.end() ) current = 0;
      else {
	current = it->second;
	dmap.erase(it);
      }
    }
    if ( !dmap.empty() )
      Throw<ColourException>()
	<< "Failed to sort dipoles in Ropewalk. "
	<< "This is a serious error - please contact the authors."
	<< Exception::abortnow;
  } else {
    map<tcPPtr, Dipole *> dmap;
    Dipole * current = 0;
    for ( int i = 0, N = dips.size(); i < N; ++i ) {
      if ( dips[i]->pc == p ) current =  dips[i];
      else dmap[dips[i]->pc] = dips[i];
    }
    dips.clear();
    while ( current ) {
      dips.push_back(current);
      map<tcPPtr, Dipole *>::iterator it = dmap.find(current->pa);
      if ( it == dmap.end() ) current = 0;
      else {
	current = it->second;
	dmap.erase(it);
      }
    }
    if ( !dmap.empty() )
      Throw<ColourException>()
	<< "Failed to sort dipoles in Ropewalk. "
	<< "This is a serious error - please contact the authors."
	<< Exception::abortnow;
  }
}

