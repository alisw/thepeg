// -*- C++ -*-
//
// HepMCHelper_HepMC.h is a part of Herwig++ - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2011 The Herwig Collaboration
//
// Herwig++ is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is a helper header to implement HepMC conversions
//
#include "ThePEG/Vectors/HepMCConverter.h"
#include "HepMC/GenEvent.h"
#include "HepMC/GenVertex.h"
#include "HepMC/GenParticle.h"
#include "HepMC/Version.h"

namespace ThePEG {
/**
 * Struct for HepMC conversion
 */
template<> 
struct HepMCTraits<HepMC::GenEvent> 
  : public HepMCTraitsBase<HepMC::GenEvent,
			   HepMC::GenParticle,
			   HepMC::GenVertex,
			   HepMC::Polarization,
			   HepMC::PdfInfo>
{
#ifdef HEPMC_VERSION_CODE
  // This is version 3!

  /** Create an event object with number \a evno and \a weight. */
  static EventT * newEvent(long evno, double,
			   const map<string,double>&) {
    EventT * e = new EventT(HepMC::Units::GEV, HepMC::Units::MM);
    e->set_event_number(evno);
//     e->weights().push_back(weight);
//     for ( map<string,double>::const_iterator w = optionalWeights.begin();
// 	  w != optionalWeights.end(); ++w ) {
// #ifdef HEPMC_HAS_NAMED_WEIGHTS
//       e->weights()[w->first] = w->second;
// #else
//       e->weights().push_back(w->second);
// #endif
//     }
    return e;
  }

  /** Set the \a scale, \f$\alpha_S\f$ (\a aS) and \f$\alpha_{EM}\f$
      (\a aEM) for the event \a e. The scale will be scaled with \a
      unit before given to the GenEvent. */
  static void setScaleAndAlphas(EventT &, Energy2,
				double,  double, Energy) {
    // e.set_event_scale(sqrt(scale)/unit);
    // e.set_alphaQCD(aS);
    // e.set_alphaQED(aEM);
  }

  /** Set the colour line (with index \a indx) to \a coline for
      particle \a p. */
  static void setColourLine(ParticleT &, int, int) {
    //    p.set_flow(indx, coline);
  }

  /** Add an incoming particle, \a p, to the vertex, \a v. */
  static void addIncoming(VertexT & v, ParticleT * p) {
    v.add_particle_in(HepMC::GenParticlePtr(p));
  }

  /** Add an outgoing particle, \a p, to the vertex, \a v. */
  static void addOutgoing(VertexT & v, ParticleT * p) {
    v.add_particle_out(HepMC::GenParticlePtr(p));
  }

  /** Set the primary vertex, \a v, for the event \a e. */
  static void setSignalProcessVertex(EventT & e, VertexT * v) {
    e.add_vertex(HepMC::GenVertexPtr(v));
  }

  /** Set a vertex, \a v, for the event \a e. */
  static void addVertex(EventT & e, VertexT * v) {
    e.add_vertex(HepMC::GenVertexPtr(v));
  }

  /** Set the beam particles for the event.*/
  static void setBeamParticles(EventT & e, ParticleT * p1, ParticleT * p2) {
    //    e.set_beam_particles(p1,p2);
    p1->set_status(4);
    p2->set_status(4);
    e.set_beam_particles(HepMC::GenParticlePtr(p1), HepMC::GenParticlePtr(p2));
  }

  /** Create a new particle object with momentum \a p, PDG number \a
      id and status code \a status. The momentum will be scaled with
      \a unit which according to the HepMC documentation should be
      GeV. */
  static ParticleT * newParticle(const Lorentz5Momentum & p,
				 long id, int status, Energy unit) {
    // Note that according to the documentation the momentum is stored in a
    // HepLorentzVector in GeV (event though the CLHEP standard is MeV).
    HepMC::FourVector p_scalar(p.x()/unit, p.y()/unit, p.z()/unit, p.e()/unit);
    ParticleT * genp = new ParticleT(p_scalar, id, status);
    genp->set_generated_mass(p.mass()/unit);
    return genp;
  }

  /** Set the polarization directions, \a the and \a phi, for particle
      \a p. */
  static void setPolarization(ParticleT &, double, double) {
    //    genp.set_polarization(PolarizationT(the, phi));
  }

  /** Set the position \a p for the vertex, \a v. The length will be
      scaled with \a unit which normally should be millimeters. */
  static void setPosition(VertexT & v, const LorentzPoint & p, Length unit) {
    HepMC::FourVector v_scaled(p.x()/unit, p.y()/unit, p.z()/unit, p.e()/unit);
    v.set_position(v_scaled);
  }

#endif

};
}
