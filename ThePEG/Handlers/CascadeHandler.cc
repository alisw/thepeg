// -*- C++ -*-
//
// CascadeHandler.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the CascadeHandler class.
//

#include "CascadeHandler.h"
#include "ThePEG/Handlers/EventHandler.h"
#include "ThePEG/Handlers/Hint.h"
#include "ThePEG/Handlers/XComb.h"
#include "ThePEG/PDF/PartonExtractor.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace ThePEG;

CascadeHandler::~CascadeHandler() {}

void CascadeHandler::
handle(EventHandler & eh, const tPVector & tagged,
       const Hint & hint)
{
  static tPVector notags;
  theLastXComb = eh.lastXCombPtr();
  theHint = &hint;
  //  if ( hint.tagged() )
    theTagged = &tagged;
    //  else
    //    theTagged = &notags;
    //
  if ( !lastXComb().empty() )
    thePDFs = make_pair(pdf<PDF>(lastPartons().first),
			pdf<PDF>(lastPartons().second));

  cascade();

}

void CascadeHandler::resetPDFs(const pair<tcPDFPtr,tcPDFPtr> & pdfpair) {
  thePDFs.first = PDF(pdfpair.first, thePDFs.first.particle());
  thePDFs.second = PDF(pdfpair.second, thePDFs.second.particle());
}

void CascadeHandler::resetPDFs(const pair<tcPDFPtr,tcPDFPtr> & pdfpair, PBPair ppair) {
  thePDFs.first = PDF(pdfpair.first, ppair.first->particle());
  thePDFs.second = PDF(pdfpair.second, ppair.second->particle());
}

double CascadeHandler::reweightCKKW(int, int) {
  return 1.0;
}

void CascadeHandler::setXComb(tXCombPtr xc) {
  theLastXComb = xc;
}

AbstractNoPIOClassDescription<CascadeHandler>
CascadeHandler::initCascadeHandler;

void CascadeHandler::Init() {

  static ClassDocumentation<CascadeHandler> documentation
    ("This is the base class to be used by all models of partonic cascades.");

}

