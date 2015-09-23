// -*- C++ -*-
//
// SpinInfo.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2011 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the SpinInfo class.
//
// Author: Peter Richardson
//
#include "SpinInfo.h"
#include "ThePEG/Repository/CurrentGenerator.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "Event.h"

using namespace ThePEG;

const double SpinInfo::_eps=1.0e-8;

SpinInfo::SpinInfo(const SpinInfo & x)
  : EventInfoBase(x), _production(x._production), _decay(x._decay),
    _timelike(x._timelike),
    _prodloc(x._prodloc), _decayloc(x._decayloc),
    _decayed(x._decayed), _developed(x._developed),_rhomatrix(x._rhomatrix),
    _Dmatrix(x._Dmatrix),_spin(x._spin),
    _productionmomentum(x._productionmomentum),
    _decaymomentum(x._decaymomentum),
    _currentmomentum(x._currentmomentum) 
{
  x._production=VertexPtr();
  x._decay=VertexPtr();
  // set the vertex so it now points to the copy
  if(_production) {
    // for timelike
    if(_timelike) _production->resetOutgoing(this,_prodloc); 
    // for spacelike
    else          _production->resetIncoming(this,_prodloc);
  }
}

EIPtr SpinInfo::clone() const {
  tcSpinPtr temp=this;
  return const_ptr_cast<SpinPtr>(temp);
}

void SpinInfo::rebind(const EventTranslationMap & trans) {
  if(_production) _production = trans.translate(_production);
  if(_decay)      _decay      = trans.translate(_decay);
  EventInfoBase::rebind(trans);
}


NoPIOClassDescription<SpinInfo> SpinInfo::initSpinInfo;
// Definition of the static class description member.

void SpinInfo::Init() {}

void SpinInfo::update() const {
  // number of instances fo this object
  int nref=referenceCount();
  if(nref<2||nref>3) return;
  // work out the number of references there should be
  int nmin=0;
  // check the production pointers
  if(_production) {
    if(_timelike) { if(_production->outgoing()[_prodloc]==this) ++nmin;}
    else          { if(_production->incoming()[_prodloc]==this) ++nmin;}
  }
  // check the decay pointers
  if(_decay) {
    if(_decay->incoming()[_decayloc]==this) ++nmin;
  }
  // delete the pointers
  SpinPtr temp;
  if(nmin+1==nref) {
    // delete the production pointers
    if(_production) {
      if(_timelike) {
	if(_production->outgoing()[_prodloc]==this)
	  _production->resetOutgoing(SpinPtr(),_prodloc);
      }
      else {
	if(_production->incoming()[_prodloc]==this)
	  _production->resetIncoming(SpinPtr(),_prodloc);
      }
    }
    // delete the decay pointers
    if(_decay) {
      if(_decay->incoming()[_decayloc]==this)
	_decay->resetIncoming(SpinPtr(),_decayloc);
    }
  }
}

void SpinInfo::decay(bool recursive) const {
  // if the particle has already been decayed do nothing
  if(_decayed) return;
  // otherwise we need to obtain the correct rho (timelike) or D (spacelike) matrix
  assert(_developed!=NeedsUpdate);
  if(_timelike) {
    if(_developed==Developed&&iSpin()!=PDT::Spin0) _developed=NeedsUpdate;
    if(productionVertex()) {
      if(recursive) redecay();
      else _rhomatrix = productionVertex()->getRhoMatrix(_prodloc,true);
    }
  }
  else {
    if(_developed==Developed&&iSpin()!=PDT::Spin0) _developed=NeedsUpdate;
    if(_production) _Dmatrix = _production->getDMatrix(_prodloc);
  }
  _decaymomentum = _currentmomentum;
  _decayed=true;
}

void SpinInfo::redevelop() const {
  assert(developed()==NeedsUpdate);
  // calculate rho/D matrix
  if(_timelike) {
    _Dmatrix   = decayVertex() ? 
      decayVertex()->getDMatrix(decayLocation()) : RhoDMatrix(iSpin());
  }
  else {
    _rhomatrix = decayVertex() ? 
      decayVertex()->getRhoMatrix(decayLocation(),false) :  RhoDMatrix(iSpin());
  }
  // update the D matrix of this spininfo
  _developed = Developed;
  // update the parent if needed
  if(productionVertex() &&
     productionVertex()->incoming().size()==1) {
    tcSpinPtr parent = _timelike ? 
      productionVertex()->incoming()[0] : productionVertex()->outgoing()[0];
    parent->needsUpdate();
    parent->redevelop();
  }
}

void SpinInfo::develop() const {
  // if the particle has already been developed do nothing
  switch(_developed) {
  case Developed:
    return;
  case NeedsUpdate:
    redevelop();
    return;
  case Undeveloped:
    if(_timelike) {
      if(_decay) _Dmatrix = _decay->getDMatrix(_decayloc);
      else       _Dmatrix = RhoDMatrix(iSpin());
    }
    else {
      if(_decay) _rhomatrix = _decay->getRhoMatrix(_decayloc,false);
      else       _rhomatrix = RhoDMatrix(iSpin());
    }
    _developed=Developed;
    return;
  }
}

void SpinInfo::redecay() const {
  if(!productionVertex()) return;
  if(productionVertex()->incoming().size()==1) {
    tcSpinPtr parent;
    if(productionVertex()->incoming()[0]->timelike())
      parent = productionVertex()->incoming()[0];
    else {
      if(productionVertex()->outgoing()[0]!=this)
	parent = productionVertex()->outgoing()[0];
      else
	parent = productionVertex()->outgoing()[1];
    }
    parent->redecay();
  }
  _rhomatrix = productionVertex()->getRhoMatrix(_prodloc,true);
}
