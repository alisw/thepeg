// -*- C++ -*-
//
// ColourSinglet.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the ColourSinglet class.
//

#include "ColourSinglet.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Utilities/UtilityBase.h"
#include "ThePEG/Repository/UseRandom.h"

using namespace ThePEG;

LorentzMomentum ColourSinglet::momentum() const {
  return Utilities::sumMomentum(partons().begin(), partons().end());
}

tcPDVector ColourSinglet::getTripletData() const {
  tcPDVector ret;
  for ( int i = 0, N = partons().size(); i < N; ++i )
    if ( parton(i)->data().iColour() != PDT::Colour8 )
      ret.push_back(parton(i)->dataPtr());
  return ret;
}

vector<ColourSinglet> ColourSinglet::getSinglets(tcParticleSet & left) {
  vector<ColourSinglet> ret;

  while ( !left.empty() ) {
    tcPPtr p = *left.begin();

    // First just remove colour singlets.
    if ( !p->coloured() || !p->hasColourInfo() ) {
      left.erase(left.begin());
      continue;
    }

    // Find a connected coulour line.
    tcColinePtr cl = p->colourLine();
    if ( !cl ) cl = p->antiColourLine();

    // Get the Colour singlet corresponding to this line.
    ret.push_back(ColourSinglet(cl, left));

  }
  return ret;
}

ColourSinglet::ColourSinglet(tcColinePtr cl, tcParticleSet & left) {

  // Follow colour line forward and add coloured partons to the first
  // string piece.
  addPiece();
  if ( !fill(1, true, cl, left) )
    // If needed also follow colourline backward and add
    // anti-coloured partons.
    fill(1, false, cl, left);

  for ( Index i = 1, N = nPieces(); i <= N; ++i )
    partons().insert(partons().end(), piece(i).begin(), piece(i).end());

}

bool ColourSinglet::
fill(Index s0, bool forward, tcColinePtr cl, tcParticleSet & left) {
  tcColinePtr first = cl;
  tcPPtr p;
  while ( (p = cl->getColouredParticle(left.begin(), left.end(), !forward)) ) {
    left.erase(p);
    if ( forward ) piece(s0).push_back(p);
    else piece(s0).push_front(p);

    if ( p->hasColourLine(first, forward) ) return true;
    if ( !( cl = p->colourLine(forward) ) ) return false;
  }

  // If we get here we have ended up in a colour source or sink.
  tColinePair fork = cl->sourceNeighbours(forward);
  if ( !fork.first || !fork.second )
    throw ColourSingletException()
      << "Inconsistent Colour flow." << Exception::eventerror;
  Junction j = addJunction(s0, forward);
  fill(j.first, !forward, fork.first, left);
  fill(j.second, !forward, fork.second, left);
  return false;
}

ColourSinglet::Junction ColourSinglet::addJunction(Index s0, bool forward) {
  // Add two new string pieces.
  Index s1 = addPiece();
  Index s2 = addPiece();

  // Connect the new pieces with the original.
  junction(s0, forward).first = s1;
  junction(s0, forward).second = s2;
  junction(s1, forward).first = s0;
  junction(s1, forward).second = s2;
  junction(s2, forward).first = s0;
  junction(s2, forward).second = s1;

  // Return the indices of the new pieces.
  return junction(s0, forward);

}

ColourSinglet ColourSinglet::splitInternal(Index sp) {
  if ( !sp ) {
    // If no piece is specified, first find all internal string pieces.
    vector<Index> internals;
    for ( Index i = 1, N = nPieces(); i <= N; ++i )
      if ( sink(i).first && sink(i).second &&
	   source(i).first && source(i).second )
	internals.push_back(i);

    // If no internal lines are found, return empty singlet.
    if ( internals.empty() ) return ColourSinglet();

    // Otherwise pick randomly between the internal lines.
    sp = UseRandom::irnd(internals.size());
  }

  vector<bool> assign(piece(sp).size());
  for ( int i = 0, N = piece(sp).size(); i < N; ++i )
    assign[i] = UseRandom::rndbool();
  return splitInternal(sp, source(sp).first,
		       UseRandom::rndbool()? sink(sp).first: sink(sp).second,
		       assign);
}

ColourSinglet ColourSinglet::
splitDiDiQuark(tcPPair qq1, tcPPair qq2, const vector<bool> & assign) {
  ColourSinglet ret;

  // Add teo new string pieces.
  Junction j = addJunction(1, true);
  sink(1) = sink(j.first) = sink(j.second) =
    source(1) = source(j.first) = source(j.second) = Junction();

  // Add the second quarks to one piece each
  piece(j.first).push_back(qq2.first);
  piece(j.second).push_back(qq2.second);

  // Add intermediate colour octets to either of the new pieces.
  for ( unsigned i = piece(1).size() - 2; i > 0; --i )
    ((i < assign.size()? assign[i]: UseRandom::rndbool() )?
     piece(j.first): piece(j.second)).push_back(piece(1)[i]);

  // Add the first quarks to one piece each
  piece(j.first).push_back(qq1.first);
  piece(j.second).push_back(qq1.second);

  // Create two new singlets and let this be one of them and return
  // the other.
  ColourSinglet csa(*this, j.first);
  ColourSinglet csb(*this, j.second);
  ret.swap(csa);
  swap(csb);
  return ret;
    
}

ColourSinglet ColourSinglet::
splitDiQuarkJunction(Index sp, tcPPtr diq, tcPPair qq,
		     const vector<bool> & assign) {
  ColourSinglet ret;
  Junction j;

  if ( source(sp).first && source(sp).second ) {
    // This is a diquark, so we remove the source junction and add
    // string piece sp to the source neighbours according to assign
    // and end with one of the quarks.
    j = source(sp);
    source(j.first) = source(j.second) = Junction();
    for ( unsigned i = 0, N = piece(sp).size(); i < N; ++i ) {
      if ( diq == piece(sp)[i] ) {
	piece(j.first).push_back(qq.first);
	piece(j.second).push_back(qq.second);
      }
      ((i < assign.size()? assign[i]: UseRandom::rndbool() )?
       piece(j.first): piece(j.second)).push_back(piece(sp)[i]);
    }
  }
  else if ( sink(sp).first && sink(sp).second ) {
    // This is a anti-diquark, so we remove the source junction and add
    // string piece sp to the source neighbours according to assign
    // and end with one of the anti quarks.
    j = sink(sp);
    sink(j.first) = sink(j.second) = Junction();
    // 'i' can't be unsigned here. would get infinite loop!
    for ( int i = piece(sp).size() - 1; i >= 0; --i ) {
      if ( diq == piece(sp)[i] ) {
	piece(j.first).push_front(qq.first);
	piece(j.second).push_front(qq.second);
      }
      ((static_cast<size_t>(i) < assign.size()? assign[i]: UseRandom::rndbool() )?
       piece(j.first): piece(j.second)).push_front(piece(sp)[i]);
    }
  }
  else
    return ret;

  // Create two new singlets and let this be one of them and return
  // the other.
  ColourSinglet csa(*this, j.first);
  ColourSinglet csb(*this, j.second);
  ret.swap(csa);
  swap(csb);
  return ret;
    
}

ColourSinglet ColourSinglet::
splitInternal(Index sp, Index sa, Index sc, const vector<bool> & assign) {
  ColourSinglet ret;

  // If the selected string piece was not internal, return nothing.
  if ( !sink(sp).first || !sink(sp).second ||
       !source(sp).first || !source(sp).second ) return ret;

  // If the selected adjacent pieces are not connected, return
  // nothing.
  Index sb = 0;
  if ( sa == source(sp).first ) sb = source(sp).second;
  else if ( sa == source(sp).second ) sb = source(sp).first;
  else return ret;
  Index sd = 0;
  if ( sc == sink(sp).first ) sd = sink(sp).second;
  else if ( sc == sink(sp).second ) sd = sink(sp).first;
  else return ret;

  // Copy the partons from the split string piece randomly to the
  // pieces to be joined.
  for ( unsigned i = 0, N = piece(sp).size(); i < N; ++i ) {
    if ( (i < assign.size()? assign[i]: UseRandom::rndbool() ) )
      piece(sa).push_back(piece(sp)[i]);
    else
      piece(sb).push_back(piece(sp)[i]);
  }

  // Join the to string pieces
  piece(sa).insert(piece(sa).end(), piece(sc).begin(), piece(sc).end());
  piece(sb).insert(piece(sb).end(), piece(sd).begin(), piece(sd).end());
  source(sa) = source(sc);
  source(sb) = source(sd);

  // Create new colour singlets
  ColourSinglet csa(*this, sa);
  ColourSinglet csb(*this, sb);
  ret.swap(csa);
  swap(csb);
  return ret;
}

ColourSinglet::ColourSinglet(const ColourSinglet & cs, Index si) {
  addPiece();
  fill(1, true, cs, si);
  fill(1, false, cs, si);

  for ( Index i = 1, N = nPieces(); i <= N; ++i )
    partons().insert(partons().end(), piece(i).begin(), piece(i).end());

}

void ColourSinglet::
fill(Index i0, bool forward, const ColourSinglet & cs, Index i1) {
  // first copy the string piece from cs.
  piece(i0) = cs.piece(i1);

  // If the piece in cs had a (forward) junction, add copy also the
  // connected string pieces
  if ( cs.junction(i1, forward).first || cs.junction(i1, forward).second ) {
    Junction j = addJunction(i0, forward);
    fill(j.first, !forward, cs, cs.junction(i1, forward).first);
    fill(j.second, !forward, cs, cs.junction(i1, forward).second);
  }
  return;
}
