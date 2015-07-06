// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the String class.
//

#include "String.h"
#include "DipoleState.h"
#include "Parton.h"
#include "Dipole.h"
#include "EMDipole.h"

#ifdef ThePEG_TEMPLATES_IN_CC_FILE
// #include "String.tcc"
#endif

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Ariadne;

void String::split(tParPtr g, tParPtr q, tParPtr qbar, bool remg) {

  // Connect dipoles
  tDipPtr odip = g->oDip();
  tDipPtr idip = g->iDip();
  odip->iPart(qbar);
  qbar->oDip(odip);
  idip->oPart(q);
  q->iDip(idip);
  if ( remg ) state()->remove(g);


  // Remove the em dipole from DipoleState and the string
  if(EMDip()){
    state()->removeEmitter(EMDip());
    EMDip(tEMDipPtr());
  }

  // If string was closed, simply break it.
  if ( endpoints().first->isG() ) {
    endpoints(qbar, q);
    qbar->string(this);
    q->string(this);
    return;
  }

  // If open string create a new one starting it qith qbar, ending
  // this with q.
  StrPtr s = state()->create<String>();
  s->endpoints(qbar, endpoints().second);
  endpoints(endpoints().first, q);

  // Connect partons to correct string.
  q->string(this);

  tParPtr p = qbar;
  do {
    p->string(s);
  } while ( p = p->next() );

}

void String::moveEndpoint(tParPtr oldend, tParPtr newend){
  tParPair end = endpoints();

  if(end.first == oldend || end.second == oldend){
    // Remove the em dipole from DipoleState and the string
    if(EMDip()){
      state()->removeEmitter(EMDip());
      EMDip(tEMDipPtr());
    }

    if(end.first == oldend){
      endpoints( newend, end.second );
    }
    else{
      endpoints( end.first, newend );
    }
  }
}

void String::join(tParPtr q, tParPtr g){
  tStrPtr str = q->string();
  if(str == tStrPtr(this)){
    return;
  }
  g->string(this);

  //Move all the partons to the current string.
  tParPair endp = str->endpoints();
  tParPtr p = endp.first;
  do {
    p->string(this);
  } while ( p = p->next() );

  //Connect the dipoles.
  if(endp.first == q){
    g->oDip(q->oDip());
    g->iDip(endpoints().second->iDip());
    g->oDip()->iPart(g);
    g->iDip()->oPart(g);
    endpoints(endpoints().first, endp.second);
  }
  else{
    g->iDip(q->iDip());
    g->oDip(endpoints().first->oDip());
    g->oDip()->iPart(g);
    g->iDip()->oPart(g);
    endpoints(endp.first, endpoints().second);
  }

  // Remove the string from DipoleState
  state()->removeString(str);

  // Remove the em dipole from DipoleState and the string
  if(EMDip()){
    state()->removeEmitter(EMDip());
    EMDip(tEMDipPtr());
  }
  if(str->EMDip()){
    state()->removeEmitter(str->EMDip());
    str->EMDip(tEMDipPtr());
  }
}

ClonePtr String::clone() const {
  return new_ptr(*this);
}

void String::rebind(const TranslationMap & trans) {
  CascadeBase::rebind(trans);
  theEndpoints.first = trans.translate(theEndpoints.first);
  theEndpoints.second = trans.translate(theEndpoints.second);
  theEMDipole = trans.translate(theEMDipole);
}

void String::persistentOutput(PersistentOStream & os) const {
  os << theEndpoints;
}

void String::persistentInput(PersistentIStream & is, int) {
  is >> theEndpoints;
}

ClassDescription<String> String::initString;
// Definition of the static class description member.

void String::Init() {}

void String::debugme() const {
  CascadeBase::debugme();
}

