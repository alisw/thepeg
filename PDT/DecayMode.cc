// -*- C++ -*-
//
// DecayMode.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DecayMode class.
//

#include "DecayMode.h"
#include "DecayMode.xh"
#include "ThePEG/Repository/Repository.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/Rebinder.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Reference.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/EnumIO.h"

using namespace ThePEG;

DecayMode::DecayMode()
  : theBrat(0.0), isOn(false) {}

DecayMode::DecayMode(tPDPtr newParticle, double newBrat, bool newOn)
  : theBrat(newBrat), isOn(newOn), theParent(newParticle) {}

DecayMode::DecayMode(const DecayMode & dm)
  : Interfaced(dm), theTag(dm.theTag), theBrat(dm.theBrat), isOn(dm.isOn),
    theParent(dm.theParent), theProducts(dm.theProducts),
    theOrderedProducts(dm.theOrderedProducts),
    theCascadeProducts(dm.theCascadeProducts), theMatchers(dm.theMatchers),
    theWildMatcher(dm.theWildMatcher), theExcluded(dm.theExcluded),
    theOverlap(dm.theOverlap), theDecayer(dm.theDecayer),
    theAntiPartner(dm.theAntiPartner), theLinks(dm.theLinks) {}

DecayMode::~DecayMode() {}

DMPtr DecayMode::Create(tPDPtr newParent, double newBrat, bool newOn) {
  DMPtr dm = new_ptr(DecayMode(newParent, newBrat, newOn));
  Repository::Register(dm, newParent->fullName() + "/NEWMODE");
  if ( !newParent->CC() ) return dm;
  DMPtr adm = new_ptr(DecayMode(newParent->CC(), newBrat, newOn));
  Repository::Register(adm, newParent->CC()->fullName() + "/NEWMODE");
  dm->theAntiPartner = adm;
  adm->theAntiPartner = dm;
  return dm;
}

void DecayMode::readSetup(istream & is) {
  string decnam;
  is >> theBrat >> ienum(isOn) >> decnam;
  if ( decnam.empty() ) return;
  BaseRepository::DirectoryAppend(decnam);
  setDecayer(BaseRepository::GetObject<tDecayerPtr>(decnam));	       
}

IBPtr DecayMode::clone() const {
  return dmclone();
}

DMPtr DecayMode::dmclone() const {
  return new_ptr(*this);
}

IBPtr DecayMode::fullclone() const {
  DMPtr dm = dmclone();
  Repository::Register(dm);
  if ( !CC() ) return dm;
  DMPtr adm = CC()->dmclone();
  Repository::Register(adm);
  dm->theAntiPartner = adm;
  adm->theAntiPartner = dm;
  return dm;
}

void DecayMode::doupdate() {
  Interfaced::doupdate();
  bool redo = touched();
  UpdateChecker::check(decayer(), redo);
  if ( !redo ) return;
  if ( theBrat > 0.0 && isOn && !decayer()->accept(*this) )
    throw DecModNoAccept(tag(), decayer()->name());
}

void DecayMode::rebind(const TranslationMap & trans) {
  try {
    theParent = trans.alwaysTranslate(theParent);
    ParticleMSet newProducts;
    trans.alwaysTranslate(inserter(newProducts),
			  products().begin(), products().end());
    products().swap(newProducts);
    tPDVector newOrdered;
    trans.alwaysTranslate(inserter(newOrdered),
			  orderedProducts().begin(), orderedProducts().end());
    theOrderedProducts.swap(newOrdered);
    ModeMSet newCasc;
    trans.alwaysTranslate(inserter(newCasc),
			  cascadeProducts().begin(), cascadeProducts().end());
    cascadeProducts().swap(newCasc);
    MatcherMSet newMatchers;
    trans.alwaysTranslate(inserter(newMatchers),
			  productMatchers().begin(), productMatchers().end());
    productMatchers().swap(newMatchers);
    wildProductMatcher() = trans.alwaysTranslate(wildProductMatcher());
    ParticleMSet newExclude;
    trans.alwaysTranslate(inserter(newExclude),
			  excluded().begin(), excluded().end());
    excluded().swap(newExclude);
    theAntiPartner = trans.alwaysTranslate(CC());

    for ( int i = 0, N = theLinks.size(); i < N; ++i ) {
      theLinks[i].first = trans.alwaysTranslate(theLinks[i].first);
      theLinks[i].second = trans.alwaysTranslate(theLinks[i].second);
    }
  }
  catch (IBPtr ip) {
    throw DecModRebind(name(), ip->name());
  }
  catch (...) {
    throw DecModRebind(name(), "<unknown>");
  }
}

IVector DecayMode::getReferences() {
  IVector ret;
  ret.push_back(theParent);
  ret.insert(ret.end(), products().begin(), products().end());
  ret.insert(ret.end(), cascadeProducts().begin(), cascadeProducts().end());
  ret.insert(ret.end(), productMatchers().begin(), productMatchers().end());
  if ( wildProductMatcher() ) ret.push_back(wildProductMatcher());
  ret.insert(ret.end(), excluded().begin(), excluded().end());
  if ( CC() ) ret.push_back(CC());
  return ret;
}

bool DecayMode::addOverlap(tcDMPtr d) {
  bool inc = includes(*d);
  if ( !inc ) return false;
  if ( find(theOverlap.begin(), theOverlap.end(), d) != theOverlap.end() )
    return true;
  theOverlap.push_back(d);
  return true;
}

void DecayMode::resetOverlap() {
  theOverlap.clear();
}

bool DecayMode::
compareId(const ParticleMSet & s1, const ParticleMSet & si2) const {
  if ( generator() ) return s1 == si2;
  ParticleMSet s2 = si2;
  for ( ParticleMSet::const_iterator p1 = s1.begin(); p1 != s1.end(); ++p1 ) {
    ParticleMSet::const_iterator p2 = s2.begin();
    while ( p2 != s2.end() && (**p2).id() != (**p1).id() ) ++p2;
    if ( p2 == s2.end() ) return false;
    s2.erase(p2);
  }
  return s2.empty();
}

ParticleMSet::const_iterator
DecayMode::findId(const ParticleMSet & s, const ParticleData & p) const {
  for ( ParticleMSet::const_iterator pit = s.begin(); pit != s.end(); ++pit )
    if ( (**pit).id() == p.id() ) return pit;
  return s.end();
}

struct IdCmp {
  bool operator()(tcPDPtr p1, tcPDPtr p2) { return p1->id() == p2->id(); }
};

bool DecayMode::includes(const DecayMode & d) const {
  // Fast check for ordinary decay modes.
  if (  cascadeProducts().empty() && productMatchers().empty() &&
	excluded().empty() &&  !wildProductMatcher() &&
	d.cascadeProducts().empty() &&  d.productMatchers().empty() &&
	d.excluded().empty() &&  !d.wildProductMatcher() ) {
    if ( links().size() != d.links().size() ) return false;
    if ( !compareId(products(), d.products()) ) return false;
    LinkVector dlinks = d.links();
    for ( int i = 0, N = links().size(); i < N; ++i ) {
      for ( int j = 0, M = dlinks.size(); j < M; ++j ) {
	if ( ( links()[i].first->id() == dlinks[j].first->id() &&
	       links()[i].second->id() == dlinks[j].second->id() ) ||
	     ( links()[i].first->id() == dlinks[j].second->id() &&
	       links()[i].second->id() == dlinks[j].first->id() ) ) {
	  dlinks.erase(dlinks.begin() + j);
	  break;
	}
      }
      return false;
    }
    return dlinks.empty();
  }

  // First check that none of the excluded products in this are
  // present in the other.
  ParticleMSet::const_iterator pit;
  for ( pit = excluded().begin(); pit != excluded().end(); ++pit ) {
    if ( findId(d.products(), **pit ) != d.products().end() ) return false;
    for ( ModeMSet::const_iterator mit = d.cascadeProducts().begin();
	  mit != d.cascadeProducts().end(); ++mit )
      if ( (**pit).id() == (**mit).parent()->id() ) return false;
  }

  // Check that all cascade decays in this overlaps with one in the
  // other. Save the ones that are left
  ModeMSet cascleft = d.cascadeProducts();
  for ( ModeMSet::iterator mit = cascadeProducts().begin();
	mit != cascadeProducts().end(); ++mit ) {
    ModeMSet::iterator mit2 = cascleft.begin();
    while ( mit2 != cascleft.end() && !(**mit).includes(**mit2) ) ++mit2;
    if ( mit2 == cascleft.end() ) return false;
  }


  // Check that all cascade product parents in the other matches
  // something in this. Otherwise expand the cascade product.
  ParticleMSet partleft = d.products();
  MatcherMSet matchleft = d.productMatchers();
  ParticleMSet excludeleft = d.excluded();
  MatcherMSet wildleft;
  if ( d.wildProductMatcher() ) wildleft.insert(wildProductMatcher());
  ParticleMSet part = products();
  MatcherMSet match = productMatchers();
  while ( cascleft.size() ) {
    ModeMSet::iterator cdmit = cascleft.begin();
    cDMPtr cdm = *cdmit;
    ParticleMSet::iterator pit = findId(part, *(cdm->parent()));
    if ( pit != part.end() ) {
      cascleft.erase(cdmit);
      part.erase(pit);
    } else {
      MatcherMSet::iterator mit = match.begin();
      while ( mit != match.end() && !(**mit).matches(*(cdm->parent())) ) ++mit;
      if ( mit != match.end() ) {
	cascleft.erase(cdmit);
	match.erase(mit);
      } else {
	if ( wildProductMatcher() &&
	     wildProductMatcher()->matches(*(cdm->parent())) ) {
	  cascleft.erase(cdmit);
	} else {
	  cascleft.erase(cdmit);
	  partleft.insert(cdm->products().begin(), cdm->products().end());
	  matchleft.insert(cdm->productMatchers().begin(),
			   cdm->productMatchers().end());
	  if ( cdm->wildProductMatcher() )
	    wildleft.insert(cdm->wildProductMatcher());
	  excludeleft.insert(cdm->excluded().begin(), cdm->excluded().end());
	  cascleft.insert(cdm->cascadeProducts().begin(),
			  cdm->cascadeProducts().end());
	}
      }
    }
  }

  // Check that all excluded left in the other are absent in this.
  if ( find_first_of(excludeleft.begin(), excludeleft.end(), part.begin(),
		     part.end(), IdCmp()) != excludeleft.end() ) return false;

  // Now all particles and matches left in this must match something
  // in the other.
  pit = part.begin();
  while ( pit != part.end() ) {
    ParticleMSet::iterator pit2 = findId(partleft, **pit++);
    if ( pit2 == partleft.end() ) return false;
    partleft.erase(pit2);
  }
  MatcherMSet::const_iterator pmit = match.begin();
  while ( pmit != match.end() ) {
    ParticleMSet::iterator pit2 = partleft.begin();
    while ( pit2 != partleft.end() && ! (**pmit).matches(**pit2) ) ++pit2;
    if ( pit2 != partleft.end() ) {
      partleft.erase(pit2);
    } else {
      MatcherMSet::iterator pmit2 = matchleft.begin();
      while ( pmit2 != matchleft.end() && ! (**pmit).matches(**pmit2) ) ++pmit2;
      if ( pmit2 != matchleft.end() ) {
	matchleft.erase(pmit2);
      } else
	return false;
    }
  }

  // Now all particles and matchers left in the other must be matched
  // by the wild match in this.
  if ( wildProductMatcher() ) {
    pit = partleft.begin();
    while ( pit != partleft.end() )
      if ( !(wildProductMatcher()->matches(**pit++)) ) return false;
    pmit = matchleft.begin();
    while (pmit != matchleft.end() )
      if ( !(wildProductMatcher()->matches(**pmit++)) ) return false;
    pmit = wildleft.begin();
    while (pmit != wildleft.end() )
      if ( !(wildProductMatcher()->matches(**pmit++)) ) return false;
  } else
    return partleft.empty() && matchleft.empty() && wildleft.empty();
  return true;
}

DMPtr DecayMode::clone(tPDPtr pd) const {
  DMPtr dm = dmclone();
  dm->theParent = pd;
  Repository::Register(dm, pd->fullName() + "/" + dm->name());
  if ( !theDecayer || !theDecayer->accept(*dm) ) dm->isOn = false;
  if ( pd->CC() ) {
    DMPtr adm = CC()? CC()->dmclone(): dmclone();
    adm->theParent = pd->CC();
    Repository::Register(adm, pd->CC()->fullName() + "/" + adm->name());
    dm->theAntiPartner = adm;
    adm->theAntiPartner = dm;
    if ( !adm->theDecayer->accept(*adm) ) adm->isOn = false;
  } else
    dm->theAntiPartner = DMPtr();
  return dm;
}

void DecayMode::synchronize() {
  if ( !CC() ) return;
  theBrat = CC()->theBrat;
  isOn = CC()->isOn;
  theDecayer = CC()->theDecayer;
}

struct ParticleOrdering {
  bool operator()(tcPDPtr p1, tcPDPtr p2) {
    return abs(p1->id()) > abs(p2->id()) ||
      ( abs(p1->id()) == abs(p2->id()) && p1->id() > p2->id() ) ||
      ( p1->id() == p2->id() && p1->fullName() > p2->fullName() );
  }
};

struct ModeOrdering {
  bool operator()(tcDMPtr d1, tcDMPtr d2) {
    ParticleOrdering ord;
    return ord(d1->parent(), d2->parent()) ||
      ( !ord(d2->parent(), d1->parent()) &&
	( d1->tag() < d2->tag() ||
	  ( d1->tag() == d2->tag() && d1->fullName() < d2->fullName() ) ) );
  }
};

struct MatcherOrdering {
  bool operator()(tcPMPtr m1, tcPMPtr m2) {
    return m1->name() < m2->name() ||
      ( m1->name() == m2->name() && m1->fullName() < m2->fullName() );
  }
};

PVector DecayMode::produceProducts() const {
  PVector ret;
  for ( int i = 0, N = orderedProducts().size(); i < N; ++i )
    ret.push_back(orderedProducts()[i]->produceParticle());
  return ret;
}

string DecayMode::makeTag() const {
  string ret;
  typedef multiset<tcPDPtr,ParticleOrdering> OrderedParticles;
  typedef multiset<tcPMPtr,MatcherOrdering> OrderedMatchers;
  typedef multiset<tcDMPtr,ModeOrdering> OrderedModes;
  LinkVector dlinks = links();

  ret = theParent->PDGName() + "->";
  if ( dlinks.empty() ) {
    OrderedParticles prod(products().begin(), products().end());
    for (OrderedParticles::iterator pit = prod.begin();
	 pit != prod.end(); ++pit )
      ret += (**pit).PDGName() + ",";
  } else {
    unsigned int dl = 0;
    for ( int i = 0, N = orderedProducts().size(); i < N; ++i ) {
      if ( dl < dlinks.size() && orderedProducts()[i] == dlinks[dl].first ) {
	ret +=  orderedProducts()[i]->PDGName() + "=";
	++dl;
      } else 
	ret +=  orderedProducts()[i]->PDGName() + ",";
    }
  }

  OrderedModes casc(cascadeProducts().begin(), cascadeProducts().end());
  for ( OrderedModes::iterator dmit = casc.begin();dmit != casc.end(); ++dmit )
    ret += "[" + (**dmit).tag() + "],";
  OrderedMatchers match(productMatchers().begin(), productMatchers().end());
  for ( OrderedMatchers::iterator mit = match.begin();
	mit != match.end(); ++mit ) ret += "?" +(**mit).name() + ",";
  if ( theWildMatcher )
    ret += "*" + theWildMatcher->name() + ",";
  OrderedParticles ex(excluded().begin(), excluded().end());
  for ( OrderedParticles::iterator pit = ex.begin(); pit != ex.end(); ++pit )
    ret += "!" + (**pit).PDGName() + ",";
  ret[ret.size()-1] = ';';
  return ret;
}

void DecayMode::brat(double newBrat) {
  theBrat = newBrat;
  if ( theBrat <= 0.0 ) switchOff();
  if ( CC() && parent()->synchronized() ) CC()->theBrat = newBrat;
}

double DecayMode::brat() const {
  return isOn? theDecayer->brat(*this, *theParent, theBrat): 0.0;
}

double DecayMode::brat(const Particle & p) const {
  return isOn && p.dataPtr() == parent()?
    theDecayer->brat(*this, p, theBrat): 0.0;
}

void DecayMode::switchOn() {
  isOn = true;
  if ( CC() && parent()->synchronized() ) CC()->isOn = true;
}

void DecayMode::switchOff() {
  isOn = false;
  if ( CC() && parent()->synchronized() ) CC()->isOn = false;
}

void DecayMode::addProduct(tPDPtr pd) {
  products().insert(pd);
  theOrderedProducts.push_back(pd);
  if ( CC() ) {
    CC()->products().insert(pd->CC()? pd->CC(): pd);
    CC()->theOrderedProducts.push_back(pd->CC()? pd->CC(): pd);
  }
  resetTag();
}

void DecayMode::addLink(tPDPtr a, tPDPtr b) {
  theLinks.push_back(make_pair(a, b));
  if ( CC() ) CC()->theLinks.push_back(make_pair(a->CC()? a->CC(): a,
						 b->CC()? b->CC(): b));
  resetTag();
}

void DecayMode::addCascadeProduct(tDMPtr dm) {
  cascadeProducts().insert(dm);
  if ( CC() ) CC()->cascadeProducts().insert(dm->CC()? dm->CC(): dm);
  resetTag();
}

void DecayMode::addProductMatcher( tPMPtr pm) {
  productMatchers().insert(pm);
  if ( CC() ) CC()->productMatchers().insert(pm->CC()? pm->CC(): pm);
  resetTag();
}

void DecayMode::setWildMatcher(tPMPtr pm) {
  wildProductMatcher() = pm;
  if ( CC() ) CC()->wildProductMatcher() = pm->CC()? pm->CC(): pm;
  resetTag();
}

void DecayMode::addExcluded(tPDPtr pd) {
  excluded().insert(pd);
  if ( CC() ) CC()->excluded().insert(pd->CC()? pd->CC(): pd);
  resetTag();
}

void DecayMode::decayer(tDecayerPtr dec) {
  if ( !dec || !dec->accept(*this) )
    throw DecModSetupNoAccept(tag(), dec->name());
  if ( CC() && parent()->synchronized() ) {
    if ( !dec->accept(*CC()) )
      throw DecModSetupNoAccept(CC()->tag(), dec->name());
    CC()->theDecayer = dec;
  }
  theDecayer = dec;
}

DMPtr DecayMode::constructDecayMode(string & tag, vector<DMPtr> * save) {
  DMPtr rdm;
  DMPtr adm;
  int level = 0;
  string::size_type end = 0;
  while ( end < tag.size() && ( tag[end] != ']' || level ) ) {
    switch ( tag[end++] ) {
    case '[':
      ++level;
      break;
    case ']':
      --level;
      break;
    }
  }
  string::size_type next = tag.find("->");
  if ( next == string::npos ) return rdm;
  if ( tag.find(';') == string::npos ) return rdm;
  tPDPtr pd = Repository::findParticle(tag.substr(0,next));
  if ( !pd ) return rdm;

  rdm = ptr_new<DMPtr>();
  rdm->parent(pd);
  if ( pd->CC() ) {
    adm = ptr_new<DMPtr>();
    adm->parent(pd->CC());
    rdm->theAntiPartner = adm;
    adm->theAntiPartner = rdm;
  }
  bool error = false;
  tag = tag.substr(next+2);
  tPDPtr lastprod;
  bool dolink = false;
  do {
    switch ( tag[0] ) {
    case '[':
      {
	tag = tag.substr(1);
	DMPtr cdm = constructDecayMode(tag, save);
	if ( save ) save->push_back(cdm);
	if ( cdm ) rdm->addCascadeProduct(cdm);
	else error = true;
      } break;
    case '=':
      dolink = true;
    case ',':
    case ']':
      tag = tag.substr(1);
      break;
    case '?':
      {
	next = min(tag.find(','), tag.find(';'));
	tPMPtr pm = Repository::findMatcher(tag.substr(1,next-1));
	if ( pm ) rdm->addProductMatcher(pm);
	else error = true;
	tag = tag.substr(next);
      } break;
    case '!':
      {
	next = min(tag.find(','), tag.find(';'));
	tPDPtr pd = Repository::findParticle(tag.substr(1,next-1));
	if ( pd ) rdm->addExcluded(pd);
	else error = true;
	tag = tag.substr(next);
      } break;
    case '*':
      {
	next = min(tag.find(','), tag.find(';'));
	tPMPtr pm = Repository::findMatcher(tag.substr(1,next-1));
	if ( pm ) rdm->setWildMatcher(pm);
	else error = true;
	tag = tag.substr(next);
      } break;
    default:
      {
	next = min(tag.find('='), min(tag.find(','), tag.find(';')));
	tPDPtr pdp = Repository::findParticle(tag.substr(0,next));
	if ( pdp ) rdm->addProduct(pdp);
	else error = true;
	tag = tag.substr(next);
	if ( dolink && lastprod ) {
	  rdm->addLink(lastprod, pdp);
	  dolink = false;
	}
	lastprod = pdp;
      } break;
    }
  } while ( tag[0] != ';' && tag.size() );
  if ( tag[0] != ';' || error ) {
    return DMPtr();
  }

  tag = tag.substr(1);
  
  for ( DecaySet::const_iterator dit = pd->decayModes().begin();
	dit != pd->decayModes().end(); ++dit )
    if ( (**dit).tag() == rdm->tag() ) return *dit;

  if ( save ) {
    save->push_back(rdm);
    save->push_back(adm);
  } else {
    pd->addDecayMode(rdm);
    Repository::Register(rdm, pd->fullName() + "/" + rdm->tag());
    if ( adm )
      Repository::Register(adm, pd->CC()->fullName() + "/" + adm->tag());
  }
  return rdm;
}

void DecayMode::persistentOutput(PersistentOStream & os) const {
  multiset<tcPDPtr,ParticleOrdering> prod(products().begin(), products().end());
  multiset<tcDMPtr,ModeOrdering>
    casc(cascadeProducts().begin(), cascadeProducts().end());
  multiset<tcPMPtr,MatcherOrdering>
    match(productMatchers().begin(), productMatchers().end());
  multiset<tcPDPtr,ParticleOrdering>  ex(excluded().begin(), excluded().end());
  multiset<tcDMPtr,ModeOrdering> ovlap(overlap().begin(), overlap().end());

  os << theTag << theBrat << isOn << theParent << prod << theOrderedProducts
     << casc << match << theWildMatcher << ex << ovlap << theDecayer
     << theAntiPartner << theLinks;
}

void DecayMode::persistentInput(PersistentIStream & is, int) {
  is >> theTag >> theBrat >> isOn >> theParent >> theProducts
     >> theOrderedProducts >> theCascadeProducts >> theMatchers
     >> theWildMatcher >> theExcluded >> theOverlap >> theDecayer
     >> theAntiPartner >> theLinks;
}

ClassDescription<DecayMode> DecayMode::initDecayMode;

void DecayMode::setOn(long i) {
  isOn = i;
}

long DecayMode::getOn() const {
  return isOn;
}

void DecayMode::setDecayer(DecayerPtr dp) {
  decayer(dp);
}

void DecayMode::Init() {

  static ClassDocumentation<DecayMode> documentation
    ("Represents a specific decay channel of a particle.");

  static Parameter<DecayMode,double> interfaceBrat
    ("BranchingRatio",
     "The branching fraction for this decay mode. Note that if the sum of "
     "branching ratios for one particle is always renormalized to 1. Also, "
     "the decaying particle may change this branching ratio if it has a "
     "ThePEG::WidthGenerator object assigned to it. ",
     &DecayMode::theBrat, 0.0, 0.0, 1.0, false, false, true);
  interfaceBrat.setHasDefault(false);

  static Switch<DecayMode> interfaceOn
    ("OnOff",
     "Indicates if the decay mode is switched on or off.",
     0, 0, false, false, &DecayMode::setOn, &DecayMode::getOn);
  static SwitchOption interfaceOnYes
    (interfaceOn, "On", "The decay channel is switched on.", 1);
  static SwitchOption interfaceOnNo
    (interfaceOn, "Off", "The decay channel is switched off.", 0);
  interfaceOn.setHasDefault(false);

  static Switch<DecayMode> interfaceActive
    ("Active",
     "Indicates if the decay mode is switched on or off.",
     0, 0, false, false, &DecayMode::setOn, &DecayMode::getOn);
  static SwitchOption interfaceActiveYes
    (interfaceActive, "Yes", "The decay channel is switched on.", 1);
  static SwitchOption interfaceActiveNo
    (interfaceActive, "No", "The decay channel is switched off.", 0);
  interfaceActive.setHasDefault(false);

  static Reference<DecayMode,Decayer> interfaceDecayer
    ("Decayer",
     "The ThePEG::Decayer object responsible for performing this decay.",
     &DecayMode::theDecayer, false, false, true, false,
     &DecayMode::setDecayer);

  interfaceBrat.rank(10);
  interfaceDecayer.rank(9);
  interfaceOn.rank(8);
  interfaceActive.rank(8);

}

DecModNoAccept::DecModNoAccept(string tag, string dec) {
  theMessage << "The Decayer '" << dec << "' is not capable to "
	     << "perform the decay in the DecayMode '" << tag << "'.";
  severity(warning);
}

DecModSetupNoAccept::DecModSetupNoAccept(string tag, string dec) {
  theMessage << "The Decayer '" << dec << "' is not capable to "
	     << "perform the decay in the DecayMode '" << tag << "'.";
  severity(warning);
}

DecModRebind::DecModRebind(string tag, string obj) {
  theMessage << "'Rebind' of DecayMode '" << tag << "' failed because "
	     << "the object '" << obj << "' refered to lacked a translation.";
  severity(abortnow);
}
