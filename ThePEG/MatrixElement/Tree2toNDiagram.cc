// -*- C++ -*-
//
// Tree2toNDiagram.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the Tree2toNDiagram class.
//

#include "Tree2toNDiagram.h"
#include "ThePEG/EventRecord/SubProcess.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/Rebinder.h"
#include "ThePEG/Utilities/UtilityBase.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace ThePEG;

Tree2toNDiagram::~Tree2toNDiagram() {}

Tree2toNDiagram & Tree2toNDiagram::add(tcPDPtr pd) {
  if ( thePartons.size() < theNSpace ) addSpacelike(pd);
  else addTimelike(pd, nextOrig);
  return *this;
}

void Tree2toNDiagram::addTimelike(tcPDPtr pd, size_type orig) {
  if ( allPartons().size() < theNSpace ||
       orig >= allPartons().size())
    throw Tree2toNDiagramError();
  thePartons.push_back(pd);
  theParents.push_back(orig);
}

tPVector Tree2toNDiagram::
construct(SubProPtr sp, const StandardXComb & xc, const ColourLines & cl) const {
  tPVector out;
  vector<Lorentz5Momentum> pout(xc.meMomenta().begin() + 2,
				xc.meMomenta().end());
//   Utilities::transform(pout.begin(), pout.end(),
// 		       Utilities::getBoostFromCM(xc.lastPartons()));
  tPPair in = xc.lastPartons();
  if ( xc.mirror() ) swap(in.first, in.second);

  tPVector ret;
  if ( in.first->dataPtr() != allPartons()[0] ||
       in.second->dataPtr() != allPartons()[nSpace() - 1] )
    throw Tree2toNDiagramError();

  PVector slike;
  slike.push_back(in.first);
  for ( int i = 1; i < nSpace() - 1; ++i )
    slike.push_back(allPartons()[i]->produceParticle());
  slike.push_back(in.second);
  ret = tPVector(slike.begin(), slike.end());
  for ( size_type i = 1; i < slike.size() - 1; ++i ) {
    slike[i-1]->addChild(slike[i]);
    sp->addIntermediate(slike[xc.mirror()? i: slike.size() - 1 - i], false);
  }
  int io = pout.size();
  PVector tlike(allPartons().size() - nSpace());
  ParticleSet done;
  for ( int i = allPartons().size() - 1; i >=  nSpace(); --i ) {
    int it = i - nSpace();
    pair<int,int> ch = children(i);
    bool iso = ch.first < 0;
    if ( iso ) {
      tlike[it] = allPartons()[i]->produceParticle(pout[--io]);
      done.insert(tlike[it]);
    } else {
      Lorentz5Momentum p = tlike[ch.first - nSpace()]->momentum() +
	tlike[ch.second - nSpace()]->momentum();
      tlike[it] = allPartons()[i]->produceParticle(p);
    }
    if ( parent(i) < nSpace() ) {
      slike[parent(i)]->addChild(tlike[it]);
      if ( parent(i) == nSpace() - 2 )
	slike[parent(i) + 1]->addChild(tlike[it]);
    }
    if ( !iso ) {
      tlike[it]->addChild(tlike[ch.first - nSpace()]);
      tlike[it]->addChild(tlike[ch.second - nSpace()]);
    }
    if ( iso ) out.push_back(tlike[it]);
    else sp->addIntermediate(tlike[it], false);
  }
  ret.insert(ret.end(), tlike.begin(), tlike.end());
  for ( int i = 0, N = out.size(); i < N; ++i )
    sp->addOutgoing(out[xc.mirror()? i: out.size() - i - 1], false);
  for ( PVector::size_type i = 0; i < slike.size() - 2; ++i ) {
    pair<int,int> ch = children(i);
    slike[ch.first]->set5Momentum(slike[i]->momentum() -
				  tlike[ch.second - nSpace()]->momentum());
  }

  cl.connect(ret);

  return out;
}

tcPDVector Tree2toNDiagram::outgoing() const {
  tcPDVector pdv;
  for ( size_type i = nSpace(); i < allPartons().size(); ++i )
    if ( children(i).first < 0 ) pdv.push_back(allPartons()[i]);
  return pdv;
}

tcPDVector Tree2toNDiagram::external() const {
  tcPDVector pdv;
  pdv.push_back(allPartons()[0]);
  pdv.push_back(allPartons()[nSpace() - 1]);
  for ( size_type i = nSpace(); i < allPartons().size(); ++i )
    if ( children(i).first < 0 ) pdv.push_back(allPartons()[i]);
  return pdv;
}

tcPDPair Tree2toNDiagram::incoming() const {
  return tcPDPair(allPartons()[0], allPartons()[nSpace() - 1]);
}

pair<int,int> Tree2toNDiagram::children(int ii) const {
  pair<int,int> ret = make_pair(-1, -1);
  for ( size_type i = 0; i < theParents.size(); ++i ) {
    if ( parent(i) == ii ) {
      if ( ret.first < 0 ) ret.first = i;
      else if ( ret.second < 0 ) ret.second = i;
      else throw Tree2toNDiagramError();
    }
  }
  return ret;
}

void Tree2toNDiagram::check() {
  vector< pair<int,int> > children(allPartons().size(), make_pair(-1, -1));
  theNOutgoing = 0;
  for ( size_type i = nSpace(); i < allPartons().size(); ++i ) {
    if ( children[parent(i)].first < 0 ) children[parent(i)].first = i;
    else if ( children[parent(i)].second < 0 ) children[parent(i)].second = i;
    else throw Tree2toNDiagramError();
  }
  for ( size_type i = nSpace(); i < allPartons().size(); ++i ) {
    if ( children[i].first < 0 && children[i].second < 0 ) ++theNOutgoing;
    else if ( children[i].first < 0 || children[i].second < 0 )
      throw Tree2toNDiagramError();
  }
  cPDVector parts(2);
  parts[0] = incoming().first;
  parts[1] = incoming().second;
  tcPDVector out(outgoing());
  parts.insert(parts.end(), out.begin(), out.end());
  partons(2, parts, nextOrig + 1);
}

bool Tree2toNDiagram::isSame (tcDiagPtr diag) const {
  Ptr<Tree2toNDiagram>::tcptr cmp = 
    dynamic_ptr_cast<Ptr<Tree2toNDiagram>::tcptr>( diag );
  if ( !cmp )
    return false;
  return equals(cmp) && external() == cmp->external();
}

bool Tree2toNDiagram::isSame (tcDiagPtr diag, map<int,int>& remap) const {
  Ptr<Tree2toNDiagram>::tcptr cmp = 
    dynamic_ptr_cast<Ptr<Tree2toNDiagram>::tcptr>( diag );
  if ( !cmp )
    return false;
  remap.clear();
  remap[0] = 0;
  return equals(cmp,remap);
}

bool Tree2toNDiagram::equals(Ptr<Tree2toNDiagram>::tcptr diag, 
			     int start, int startCmp) const {

  if ( start < 0 && startCmp < 0 )
    return true;

  if ( allPartons()[start] != diag->allPartons()[startCmp] )
    return false;

  pair<int,int> ch = children(start);
  pair<int,int> chCmp = diag->children(startCmp);

  bool match =
    equals(diag,ch.first,chCmp.first) &&
    equals(diag,ch.second,chCmp.second);

  // also try swapped outgoing legs on same vertex
  if ( !match && start > nSpace() - 1 &&
       children(ch.first).first < 0 && children(ch.second).first < 0 &&
       diag->children(chCmp.first).first < 0 && diag->children(chCmp.second).first < 0 )
    match = 
      equals(diag,ch.first,chCmp.second) &&
      equals(diag,ch.second,chCmp.first);

  return match;
    
}

bool Tree2toNDiagram::equals(Ptr<Tree2toNDiagram>::tcptr diag, 
			     map<int,int>& remap,
			     int start, int startCmp) const {

  if ( start < 0 && startCmp < 0 )
    return true;

  if ( allPartons()[start] != diag->allPartons()[startCmp] )
    return false;

  pair<int,int> ch = children(start);
  pair<int,int> chCmp = diag->children(startCmp);

  if ( ch.first < 0 && chCmp.first < 0 ) {
    remap[externalId(start)] = diag->externalId(startCmp);
  }

  bool match =
    equals(diag,remap,ch.first,chCmp.first) &&
    equals(diag,remap,ch.second,chCmp.second);

  // also try swapped outgoing legs on same vertex
  if ( !match && start > nSpace() - 1 &&
       children(ch.first).first < 0 && children(ch.second).first < 0 &&
       diag->children(chCmp.first).first < 0 && diag->children(chCmp.second).first < 0 )
    match = 
      equals(diag,remap,ch.first,chCmp.second) &&
      equals(diag,remap,ch.second,chCmp.first);

  return match;

}

int Tree2toNDiagram::externalId(int id) const {
  if ( id < 0 )
    return -1;
  if ( id == 0 )
    return 0;
  if ( id == nSpace() - 1 )
    return 1;
  int k = 1;
  for ( size_type i = nSpace(); i < allPartons().size(); ++i ) {
    if ( children(i).first < 0 ) ++k;
    if ( i == size_type(id) )
      break;
  }
  return k;
}

int Tree2toNDiagram::diagramId(int id) const {
  if ( id < 0 )
    return -1;
  if ( id == 0 ) return 0;
  if ( id == 1 ) return nSpace() - 1;
  int k = 1;
  size_type i = nSpace();
  for ( ; i < allPartons().size(); ++i ) {
    if ( children(i).first < 0 ) ++k;
    if ( k == id )
      break;
  }
  return i;
}

int Tree2toNDiagram::mergeEmission(int emitter, int id, map<int,int>& remap) {

  if ( id < 2 )
    return -1;

  if ( remap.find(emitter) != remap.end() ) {
    remap.erase(emitter);
  }
  if ( remap.find(id) != remap.end() ) {
    remap.erase(id);
  }

  for ( map<int,int>::iterator rm = remap.begin();
	rm != remap.end(); ++rm ) {
    if ( rm->first == 0 || rm->first == 1 ) {
      rm->second = rm->first;
    } else {
      rm->second = diagramId(rm->first);
    }
  }

  // translate to diagram id
  int did = diagramId(id);
  int demitter = diagramId(emitter);

  if ( children(did) != make_pair(-1,-1) )
    return -1;

  // now get the parent
  int p = parent(did);

  int npos = -1;
  if ( p == 0 || p == nSpace() - 2 ) {
    npos = ( p == 0 ? 0 : 1 );
  } else if ( p >= nSpace() ) {
    if ( id > emitter )
      npos = emitter;
    else
      npos = emitter - 1;
  }

  pair<int,int> remove;

  size_type theNSpaceBackup = theNSpace;
  int theNOutgoingBackup = theNOutgoing;
  int nextOrigBackup = nextOrig;
  cPDVector thePartonsBackup = thePartons;
  vector<int> theParentsBackup = theParents;

  int deltaFlow = 0;
  if ( npos == 1 ) {
    if ( thePartons[did]->CC() )
      deltaFlow -= ( thePartons[did]->id() < 0 ? -1 : 1 );
    if ( thePartons[nSpace()-1]->CC() )
      deltaFlow += ( thePartons[nSpace()-1]->id() < 0 ? -1 : 1 );
  }

  // emitted from spacelike
  if ( p == 0 || p == nSpace() - 2 ) {
    if ( p == 0 && p != demitter )
      return -1;
    if ( p == nSpace() - 2 && demitter != nSpace()-1 )
      return -1;
    if ( p == 0 )
      remove = make_pair(p,did);
    else
      remove = make_pair(nSpace()-1,did);
    --theNSpace;
    --theNOutgoing;
  } else if ( p >= nSpace() ) {
    remove = children(p);
    if ( remove.first != demitter )
      swap(remove.first,remove.second);
    if ( remove != make_pair(demitter,did) )
      return -1;
    --theNOutgoing;
  } else {
    return -1;
  }

  if ( remove.first > remove.second )
    swap(remove.first,remove.second);

  for ( map<int,int>::iterator rm = remap.begin();
	rm != remap.end(); ++rm ) {
    if ( rm->first > 1 ) {
      if ( rm->second > remove.first &&
	   rm->second < remove.second )
	rm->second -= 1;
      else if ( rm->second > remove.second )
	rm->second -= 2;
    }
  }

  for ( unsigned int k = remove.first + 1; k < theParents.size(); ++k ) {
    if ( theParents[k] >= remove.first && 
	 theParents[k] < remove.second &&
	 theParents[k] >= 0 )
      theParents[k] -= 1;
    else if ( theParents[k] > remove.second && theParents[k] > 0 )
      theParents[k] -= 2;
  }
  thePartons.erase(thePartons.begin() + remove.second);
  theParents.erase(theParents.begin() + remove.second);
  thePartons.erase(thePartons.begin() + remove.first);
  theParents.erase(theParents.begin() + remove.first);

  if ( npos > 1 )
    if ( npos != externalId(p) ) {
      pair<int,int> swapDiagIds(p,diagramId(npos));
      swap(thePartons[swapDiagIds.first],thePartons[swapDiagIds.second]);
      swap(theParents[swapDiagIds.first],theParents[swapDiagIds.second]);
      for ( map<int,int>::iterator rm = remap.begin();
	    rm != remap.end(); ++rm ) {
	if ( rm->first > 1 ) {
	  if ( rm->second == swapDiagIds.first ) {
	    rm->second = swapDiagIds.second;
	  } else if ( rm->second == swapDiagIds.second ) {
	    rm->second = swapDiagIds.first;
	  }
	}
      }
    }

  for ( map<int,int>::iterator rm = remap.begin();
	rm != remap.end(); ++rm ) {
    if ( rm->first > 1 ) {
      rm->second = externalId(rm->second);
    }
  }

  if ( npos == 1 ) {
    if ( thePartons[nSpace()-1]->CC() )
      deltaFlow -= ( thePartons[nSpace()-1]->id() < 0 ? -1 : 1 );

    if ( deltaFlow != 0 )
      thePartons[nSpace()-1] = thePartons[nSpace()-1]->CC();

  }

  try {
    check();
  } catch (Tree2toNDiagramError&) {
    theNSpace = theNSpaceBackup;
    theNOutgoing = theNOutgoingBackup;
    nextOrig = nextOrigBackup;
    thePartons = thePartonsBackup;
    theParents = theParentsBackup;
    return -1;
  }

  return npos;

}
 
ClassDescription<Tree2toNDiagram> Tree2toNDiagram::initTree2toNDiagram;

void Tree2toNDiagram::persistentInput(PersistentIStream & is, int) {
  is >> theNSpace >> theNOutgoing >> thePartons >> theParents >> nextOrig;
}

void Tree2toNDiagram::persistentOutput(PersistentOStream & os) const {
  os << theNSpace << theNOutgoing << thePartons << theParents << nextOrig;
}

Tree2toNDiagramError::Tree2toNDiagramError() {
  theMessage << "An error occurred while setting up a diagram of class "
	     << "'Tree2toNDiagram'.";
  severity(abortnow);
}

