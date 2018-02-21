// -*- C++ -*-
//
// BlobDiagram.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the BlobDiagram class.
//

#include "BlobDiagram.h"
#include "ThePEG/EventRecord/SubProcess.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Utilities/Rebinder.h"
#include "ThePEG/Utilities/UtilityBase.h"
#include "ThePEG/Utilities/Exception.h"
#include "ThePEG/Handlers/StandardXComb.h"
#include "ThePEG/PDT/EnumParticles.h"

using namespace ThePEG;

BlobDiagram::~BlobDiagram() {}

tPVector BlobDiagram::
construct(SubProPtr sp, const StandardXComb & xc, const ColourLines & cl) const {

  tPVector out;
  vector<Lorentz5Momentum> pout(xc.meMomenta().begin() + 2, xc.meMomenta().end());
  if ( xc.needsReshuffling() )
    xc.reshuffle(pout); 
  tPPair in = xc.lastPartons();
  if ( xc.mirror() ) swap(in.first, in.second);

  tPVector ret;
  if ( in.first->dataPtr() != partons()[0] ||
       in.second->dataPtr() != partons()[1] )
    throw Exception() << "incoming partons in XComb do not match incoming partons in BlobDiagram"
		      << Exception::setuperror;

  PVector slike;
  slike.push_back(in.first);
  slike.push_back(in.second);

  ret = tPVector(slike.begin(), slike.end());
  for ( size_type i = 1; i < slike.size() - 1; ++i ) {
    slike[i-1]->addChild(slike[i]);
    sp->addIntermediate(slike[xc.mirror()? i: slike.size() - 1 - i], false); 
  }

  int io = pout.size();
  PVector tlike(partons().size() - 2);
  ParticleSet done;
  
  for ( int i = partons().size() - 1; i >=  2; --i ) {
    int it = i - 2;
    tlike[it] = partons()[i]->produceParticle(pout[--io]);
    done.insert(tlike[it]);
    // add the time-like parton as the child of both incoming (space-like) partons. 
    slike[0]->addChild(tlike[it]);
    slike[1]->addChild(tlike[it]);
    out.push_back(tlike[it]); 
  }

  ret.insert(ret.end(), tlike.begin(), tlike.end());
  for ( int i = 0, N = out.size(); i < N; ++i )
    sp->addOutgoing(out[xc.mirror()? i: out.size() - i - 1], false);

  cl.connect(ret);

  return out;

}

ClassDescription<BlobDiagram> BlobDiagram::initBlobDiagram;

void BlobDiagram::persistentInput(PersistentIStream &, int) {}

void BlobDiagram::persistentOutput(PersistentOStream &) const {}
