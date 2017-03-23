// -*- C++ -*-
//
// VVVVVertex.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2011 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VVVVVertex class.
//

#include "VVVVVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace ThePEG;
using namespace Helicity;
    
AbstractNoPIOClassDescription<VVVVVertex> VVVVVertex::initVVVVVertex;
// Definition of the static class description member.
    
void VVVVVertex::Init() {
      
static ClassDocumentation<VVVVVertex> documentation
  ("The VVVVVertex class is the implementation of the 4-vector vertex");
 
}

// calculate the vertex
Complex VVVVVertex::evaluate(Energy2 q2 , int iopt, 
			     const VectorWaveFunction & vec1,
			     const VectorWaveFunction & vec2,
			     const VectorWaveFunction & vec3,
			     const VectorWaveFunction & vec4) {
  // workout the coupling
  setCoupling(q2,vec1.particle(),vec2.particle(),
	      vec3.particle(),vec4.particle());
  Complex vertex,ii(0.,1.);
  // calculate the vertex
  // QCD type vertex
  assert(_itype>=1&&_itype<=2);
  if(_itype==1) {
    // dot products we need
    Complex dotv1v2 = vec1.wave().dot(vec2.wave());
    Complex dotv3v4 = vec3.wave().dot(vec4.wave());
    Complex dotv1v4 = vec1.wave().dot(vec4.wave());
    Complex dotv2v3 = vec3.wave().dot(vec2.wave());
    // first the 4-point part of the vertex
    vertex = dotv1v2*dotv3v4-dotv1v4*dotv2v3;
    // now the virtual gluon exchange if needed
    if(iopt!=0) {
      // dot products
      Complex dotv1v3 = vec1.wave().dot(vec3.wave());
      Complex dotv2v4 = vec4.wave().dot(vec2.wave());
      complex<Energy> dotv1p13 = vec1.wave().dot(2.*vec3.momentum() +
						 vec1.momentum());
      complex<Energy> dotv2p24 = vec2.wave().dot(2.*vec4.momentum() +
						 vec2.momentum());
      complex<Energy> dotv3p13 = vec3.wave().dot(2.*vec1.momentum() +
						 vec3.momentum());
      complex<Energy> dotv4p24 = vec4.wave().dot(2.*vec2.momentum() +
						 vec4.momentum());
      LorentzPolarizationVectorE veca = 
	dotv3p13*vec1.wave() - dotv1p13*vec3.wave() +
	dotv1v3*(vec3.momentum()-vec1.momentum());
      LorentzPolarizationVectorE vecb = 
	dotv4p24*vec2.wave() - dotv2p24*vec4.wave() +
	dotv2v4*(vec4.momentum()-vec2.momentum());
      InvEnergy2 numerator = 1./(vec1.momentum()+vec3.momentum()).m2();
      vertex += numerator*veca.dot(vecb);
    }
  }
  // EW type vertex
  else if(_itype==2) {
    Complex dotv1v2 = vec1.wave().dot(vec2.wave());
    Complex dotv1v3 = vec1.wave().dot(vec3.wave());
    Complex dotv1v4 = vec1.wave().dot(vec4.wave());
    Complex dotv2v3 = vec2.wave().dot(vec3.wave());
    Complex dotv2v4 = vec2.wave().dot(vec4.wave());
    Complex dotv3v4 = vec3.wave().dot(vec4.wave());
    // evaluate the vertex
    // need to sort the order out here
    if(( _iorder[0]==0 && _iorder[1]==1 && _iorder[2]==2 && _iorder[3]==3)||
       ( _iorder[0]==1 && _iorder[1]==0 && _iorder[2]==2 && _iorder[3]==3)||
       ( _iorder[0]==0 && _iorder[1]==1 && _iorder[2]==3 && _iorder[3]==2)||
       ( _iorder[0]==1 && _iorder[1]==0 && _iorder[2]==3 && _iorder[3]==2)||
       ( _iorder[0]==2 && _iorder[1]==3 && _iorder[2]==0 && _iorder[3]==1)||
       ( _iorder[0]==2 && _iorder[1]==3 && _iorder[2]==1 && _iorder[3]==0)||
       ( _iorder[0]==3 && _iorder[1]==2 && _iorder[2]==0 && _iorder[3]==1)||
       ( _iorder[0]==3 && _iorder[1]==2 && _iorder[2]==1 && _iorder[3]==0)) {
      // contact term
      vertex = 2.*dotv1v2*dotv3v4-dotv1v3*dotv2v4-dotv1v4*dotv2v3;
      // now for the u- and t-channel terms if needed
      if(iopt!=0) {
	// dot products of momenta and wavefunction
	complex<Energy> dotv1p13 =
	  +vec1.wave().dot(vec1.momentum() + 2.*vec3.momentum() );
	complex<Energy> dotv1p14 =
	  +vec1.wave().dot(vec1.momentum() + 2.*vec4.momentum() );
	complex<Energy> dotv2p23 =
	  +vec2.wave().dot(vec2.momentum() + 2.*vec3.momentum() );
	complex<Energy> dotv2p24 =
	  +vec2.wave().dot(vec2.momentum() + 2.*vec4.momentum() );
	complex<Energy> dotv3p31 = 
	  +vec3.wave().dot(vec3.momentum() + 2.*vec1.momentum() );
	complex<Energy> dotv3p32 = 
	  +vec3.wave().dot(vec3.momentum() + 2.*vec2.momentum() );
	complex<Energy> dotv4p41 = 
	  +vec4.wave().dot(vec4.momentum() + 2.*vec1.momentum() );
	complex<Energy> dotv4p42 = 
	  +vec4.wave().dot(vec4.momentum() + 2.*vec2.momentum() );
	LorentzPolarizationVectorE ja = 
	  (vec3.momentum()-vec1.momentum())*dotv1v3
	  +dotv3p31*vec1.wave()-dotv1p13*vec3.wave();
	LorentzPolarizationVectorE jb =
	  (vec4.momentum()-vec2.momentum())*dotv2v4
	  +dotv4p42*vec2.wave()-dotv2p24*vec4.wave(); 
	LorentzPolarizationVectorE jc =
	  (vec4.momentum()-vec1.momentum())*dotv1v4
	  +dotv4p41*vec1.wave()-dotv1p14*vec4.wave();
	LorentzPolarizationVectorE jd =
	  (vec3.momentum()-vec2.momentum())*dotv2v3
	  +dotv3p32*vec2.wave()-dotv2p23*vec3.wave();
	// dot products of these vectors
	complex<Energy2> dotjajb = ja.dot(jb);
	complex<Energy2> dotjcjd = jc.dot(jd);
	complex<Energy2> dotjaq  = ja.dot(vec1.momentum()+vec3.momentum());
	complex<Energy2> dotjbq  = jb.dot(vec1.momentum()+vec3.momentum());
	complex<Energy2> dotjck  = jc.dot(vec1.momentum()+vec4.momentum());
	complex<Energy2> dotjdk  = jd.dot(vec1.momentum()+vec4.momentum());
	Energy2 q2 = (vec1.momentum()+vec3.momentum()).m2();
	Energy2 k2 = (vec1.momentum()+vec4.momentum()).m2();
	// compute the term we need
	Energy2 mass2;
	for(int ix=0;ix<2;++ix) {
	  if(_inter[ix]) {
	    mass2 = sqr(_inter[ix]->mass());
	    if(mass2!=Energy2()) {
	      vertex += UnitRemoval::InvE2 *
		_coup[ix]*propagator(iopt,q2,_inter[ix])*
		(dotjajb-dotjaq*dotjbq/mass2);
	      vertex += UnitRemoval::InvE2 *
		_coup[ix]*propagator(iopt,k2,_inter[ix])*
		(dotjcjd-dotjck*dotjdk/mass2);
	    }
	    else {
	      vertex+=UnitRemoval::InvE2 *_coup[ix]*
		propagator(iopt,q2,_inter[ix])*dotjajb;
	      vertex+=UnitRemoval::InvE2 *_coup[ix]*
		propagator(iopt,k2,_inter[ix])*dotjcjd;
	    }
	  }
	}
      }
    }
    else if(( _iorder[0]==0 && _iorder[1]==2 && _iorder[2]==1 && _iorder[3]==3)||
	    ( _iorder[0]==2 && _iorder[1]==0 && _iorder[2]==1 && _iorder[3]==3)||
	    ( _iorder[0]==0 && _iorder[1]==2 && _iorder[2]==3 && _iorder[3]==1)||
	    ( _iorder[0]==2 && _iorder[1]==0 && _iorder[2]==3 && _iorder[3]==1)||
	    ( _iorder[0]==1 && _iorder[1]==3 && _iorder[2]==0 && _iorder[3]==2)||
	    ( _iorder[0]==1 && _iorder[1]==3 && _iorder[2]==2 && _iorder[3]==0)||
	    ( _iorder[0]==3 && _iorder[1]==1 && _iorder[2]==0 && _iorder[3]==2)||
	    ( _iorder[0]==3 && _iorder[1]==1 && _iorder[2]==2 && _iorder[3]==0)) {
      // contact term
      vertex = 2.*dotv1v3*dotv2v4-dotv1v2*dotv3v4-dotv1v4*dotv2v3;
      // now for the u- and t-channel terms if needed
      if(iopt!=0) {
	// dot products of momenta and wavefunction
	complex<Energy> dotv1p12 =
	  vec1.wave().dot(vec1.momentum() + 2.*vec2.momentum() );
	complex<Energy> dotv1p14 =
	  vec1.wave().dot(vec1.momentum() + 2.*vec4.momentum() );
	complex<Energy> dotv3p32 =
	  vec3.wave().dot(vec3.momentum() + 2.*vec2.momentum() );
	complex<Energy> dotv3p34 =
	  vec3.wave().dot(vec3.momentum() + 2.*vec4.momentum() );
	complex<Energy> dotv2p21 = 
	  vec2.wave().dot(vec2.momentum() + 2.*vec1.momentum() );
	complex<Energy> dotv2p23 = 
	  vec2.wave().dot(vec2.momentum() + 2.*vec3.momentum() );
	complex<Energy> dotv4p41 = 
	  vec4.wave().dot(vec4.momentum() + 2.*vec1.momentum() );
	complex<Energy> dotv4p43 = 
	  vec4.wave().dot(vec4.momentum() + 2.*vec3.momentum() );
	LorentzPolarizationVectorE ja = 
	  (vec2.momentum() - vec1.momentum() )*dotv1v2 +
	  dotv2p21*vec1.wave() - dotv1p12*vec2.wave();
	LorentzPolarizationVectorE jb = 
	  (vec4.momentum() - vec3.momentum())*dotv3v4 +
	  dotv4p43*vec3.wave() - dotv3p34*vec4.wave();
	LorentzPolarizationVectorE jc = 
	  (vec4.momentum() - vec1.momentum())*dotv1v4 +
	  dotv4p41*vec1.wave() - dotv1p14*vec4.wave();
	LorentzPolarizationVectorE jd = 
	  (vec2.momentum() - vec3.momentum())*dotv2v3 + 
	  dotv2p23*vec3.wave() - dotv3p32*vec2.wave();
	// dot products of these vectors
	complex<Energy2> dotjajb = ja.dot(jb);
	complex<Energy2> dotjcjd = jc.dot(jd);
	complex<Energy2> dotjaq  = ja.dot(vec1.momentum()+vec2.momentum());
	complex<Energy2> dotjbq  = jb.dot(vec1.momentum()+vec2.momentum());
	complex<Energy2> dotjck  = jc.dot(vec1.momentum()+vec4.momentum());
	complex<Energy2> dotjdk  = jd.dot(vec1.momentum()+vec4.momentum());
	Energy2 q2 = (vec1.momentum()+vec2.momentum()).m2();
	Energy2 k2 = (vec1.momentum()+vec4.momentum()).m2();
	// compute the term we need
	Energy2 mass2;
	for(int ix=0;ix<2;++ix) {
	  if(_inter[ix]) {
	    mass2 = (_inter[ix]->mass())*(_inter[ix]->mass());
	    if(mass2!=Energy2()) {
	      vertex+=UnitRemoval::InvE2 *
		_coup[ix]*propagator(iopt,q2,_inter[ix])*
		(dotjajb-dotjaq*dotjbq/mass2);
	      vertex+=UnitRemoval::InvE2 *
		_coup[ix]*propagator(iopt,k2,_inter[ix])*
		(dotjcjd-dotjck*dotjdk/mass2);
	    }
	    else {
	      vertex+=UnitRemoval::InvE2 *_coup[ix]*
		propagator(iopt,q2,_inter[ix])*dotjajb;
	      vertex+=UnitRemoval::InvE2 *_coup[ix]*
		propagator(iopt,k2,_inter[ix])*dotjcjd;
	    }
	  }
	}
      }
    }
    else if(( _iorder[0]==0 && _iorder[1]==3 && _iorder[2]==1 && _iorder[3]==2)||
	    ( _iorder[0]==3 && _iorder[1]==0 && _iorder[2]==1 && _iorder[3]==2)||
	    ( _iorder[0]==0 && _iorder[1]==3 && _iorder[2]==2 && _iorder[3]==1)||
	    ( _iorder[0]==3 && _iorder[1]==0 && _iorder[2]==2 && _iorder[3]==1)||
	    ( _iorder[0]==1 && _iorder[1]==2 && _iorder[2]==0 && _iorder[3]==3)||
	    ( _iorder[0]==2 && _iorder[1]==1 && _iorder[2]==0 && _iorder[3]==3)||
	    ( _iorder[0]==1 && _iorder[1]==2 && _iorder[2]==3 && _iorder[3]==0)||
	    ( _iorder[0]==2 && _iorder[1]==1 && _iorder[2]==3 && _iorder[3]==0)) {
      // contact term
      vertex = 2.*dotv1v4*dotv2v3-dotv1v3*dotv2v4-dotv1v2*dotv3v4;
      // now for the u- and t-channel terms if needed
      if(iopt!=0) {
	// dot products of momenta and wavefunction
	complex<Energy> dotv1p12 =
	  vec1.wave().dot(vec1.momentum() + 2.*vec2.momentum() );
	complex<Energy> dotv1p13 =
	  vec1.wave().dot(vec1.momentum() + 2.*vec3.momentum() );
	complex<Energy> dotv2p24 = 
	  vec2.wave().dot(vec2.momentum() + 2.*vec4.momentum() );
	complex<Energy> dotv2p21 = 
	  vec2.wave().dot(vec2.momentum() + 2.*vec1.momentum() );
	complex<Energy> dotv3p31 = 
	  vec3.wave().dot(vec3.momentum() + 2.*vec1.momentum() );
	complex<Energy> dotv3p34 =
	  vec3.wave().dot(vec3.momentum() + 2.*vec4.momentum() );
	complex<Energy> dotv4p43 = 
	  vec4.wave().dot(vec4.momentum() + 2.*vec3.momentum() );
	complex<Energy> dotv4p42 =
	  vec4.wave().dot(vec4.momentum() + 2.*vec2.momentum() );
	LorentzPolarizationVectorE ja = 
	  (vec2.momentum()-vec1.momentum())*dotv1v2
	  +dotv2p21*vec1.wave()-dotv1p12*vec2.wave();
	LorentzPolarizationVectorE jb = 
	  (vec3.momentum()-vec4.momentum())*dotv3v4
	  +dotv3p34*vec4.wave()-dotv4p43*vec3.wave();
	LorentzPolarizationVectorE jc = 
	  (vec3.momentum()-vec1.momentum())*dotv1v3
	  +dotv3p31*vec1.wave()-dotv1p13*vec3.wave();
	LorentzPolarizationVectorE jd = 
	  (vec2.momentum()-vec4.momentum())*dotv2v4
	  +dotv2p24*vec4.wave()-dotv4p42*vec2.wave();
	// dot products of these vectors
	complex<Energy2> dotjajb = ja.dot(jb);
	complex<Energy2> dotjcjd = jc.dot(jd);
	complex<Energy2> dotjaq = ja.dot(vec1.momentum()+vec2.momentum());
	complex<Energy2> dotjbq = jb.dot(vec1.momentum()+vec2.momentum());
	complex<Energy2> dotjck = jc.dot(vec1.momentum()+vec3.momentum());
	complex<Energy2> dotjdk = jd.dot(vec1.momentum()+vec3.momentum());
	Energy2 q2 = (vec1.momentum()+vec2.momentum()).m2();
	Energy2 k2 = (vec1.momentum()+vec3.momentum()).m2();
	// compute the term we need
	Energy2 mass2;
	for(int ix=0;ix<2;++ix) {
	  if(_inter[ix]) {
	    mass2 = sqr(_inter[ix]->mass());
	    if(mass2!=Energy2()) {
	      vertex+=UnitRemoval::InvE2 *
		_coup[ix]*propagator(iopt,q2,_inter[ix])*
		(dotjajb-dotjaq*dotjbq/mass2);
	      vertex+=UnitRemoval::InvE2 *
		_coup[ix]*propagator(iopt,k2,_inter[ix])*
		(dotjcjd-dotjck*dotjdk/mass2);
	    }
	    else {
	      vertex+=UnitRemoval::InvE2 *_coup[ix]*
		propagator(iopt,q2,_inter[ix])*dotjajb;
	      vertex+=UnitRemoval::InvE2 *_coup[ix]*
		propagator(iopt,k2,_inter[ix])*dotjcjd;
	    }
	  }
	}
      }
    }
    else
      throw HelicityConsistencyError() << "Unknown order of particles in "
				       << "VVVVVertex::evaluate()"
				       << Exception::runerror;
  }
  // return the answer
  return  -ii*norm()*vertex;
}
