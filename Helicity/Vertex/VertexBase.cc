// -*- C++ -*-
//
// VertexBase.cc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 2003-2017 Peter Richardson, Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the VertexBase class.
//

#include "VertexBase.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Utilities/Rebinder.h"
#include <ThePEG/PDT/ParticleData.h>
#include <ThePEG/PDT/WidthGenerator.h>
#include <iterator>

using namespace ThePEG;
using namespace ThePEG::Helicity;

VertexBase::VertexBase(VertexType::T name, bool kine) 
  : _npoint(1), _norm(0), _calckinematics(kine), 
    _kine(), _theName(name), 
    _ordergEM(0), _ordergS(0),
    _coupopt(0), _gs(sqrt(4.*Constants::pi*0.3)), 
    _ee(sqrt(4.*Constants::pi/137.04)),
    _sw(sqrt(0.232)) 
{
  assert ( name != VertexType::UNDEFINED ); 
  // Count number of lines from length of 'name'
  while ( name /= 10 ) ++_npoint;
}

// setup the lists of particles 
// should only be called from child class constructors
void VertexBase::addToList(long ida, long idb, long idc, long idd) {
  if ( idd == 0 ) addToList({ ida, idb, idc });
  else            addToList({ ida, idb, idc, idd });
}

void VertexBase::addToList(const vector<long> & ids) {
  assert( ids.size() == _npoint );
  vector<PDPtr> tmp;
  int chargeSum = 0;
  for ( auto id : ids ) {
    tPDPtr p = getParticleData(id);
    if ( !p ) return; // needed e.g. to deal with chi_5 in MSSM
    tmp.push_back(p);
    chargeSum += p->iCharge();
  }
  assert( tmp.size() == _npoint );
  if ( chargeSum != 0 ) {
    cerr << "Problem with the addToList() calls in "
	 << fullName() << ":\n"
	 << "Vertex particles ";
    copy (ids.begin(), ids.end(),
          std::ostream_iterator<long>(cerr," "));
    cerr << "have non-zero electric charge " << chargeSum << "/3.\n";
    assert( false );
  }
  _particles.push_back(tmp);
}

void VertexBase::doinit() {
  Interfaced::doinit();
  // set up the incoming and outgoing particles
  if ( !_outpart.empty() || !_inpart.empty() )
    return;
  for ( const auto & pvec : _particles ) {
    for ( tPDPtr p : pvec ) {
      assert( p );
      tPDPtr cc = p->CC();
      _inpart.insert( cc ? cc : p );
      _outpart.insert(p);
    }
  }
  // check the couplings
  if(Debug::level>1&&_npoint!=2+_ordergEM+_ordergS)
    generator()->log() << fullName() << " has inconsistent number of "
		       << "external particles and coupling order\nQED = " 
		       << _ordergEM << " QCD = " << _ordergS << " for"
		       << " a perturbative interaction. Either it's an"
		       << " effective vertex or something is wrong.\n";
  assert(_npoint<=2+_ordergEM+_ordergS);
}

    
void VertexBase::persistentOutput(PersistentOStream & os) const {
  os << _npoint << _inpart << _outpart 
     << _particles << _calckinematics
     << _coupopt << _gs << _ee << _sw;
}

void VertexBase::persistentInput(PersistentIStream & is, int) {
  is >> _npoint >> _inpart >> _outpart 
     >> _particles >> _calckinematics
     >> _coupopt >> _gs >> _ee >> _sw;
}

// Static variable needed for the type description system in ThePEG.
DescribeAbstractClass<VertexBase,Interfaced>
describeThePEGVertexBase("ThePEG::VertexBase", "");

void VertexBase::Init() {
  
  static Switch<VertexBase,bool> interfaceCalculateKinematics
    ("CalculateKinematics",
     "Calculate kinematic invariants at the vertices. This is"
     " mainly needed for loop vertices.",
     &VertexBase::_calckinematics, false, false, false);
  static SwitchOption interfaceCalculateKinematicsCalculate
    (interfaceCalculateKinematics,
     "Calculate",
     "Calculate the kinematics",
     true);
  static SwitchOption interfaceCalculateKinematicsNoKinematics
    (interfaceCalculateKinematics,
     "NoKinematics",
     "Do not calculate the kinematics",
     false);

  static ClassDocumentation<VertexBase> documentation
    ("The VertexBase class is designed to be the base class"
     "of all vertices.");

  static Switch<VertexBase,unsigned int> interfaceCoupling
    ("Coupling",
     "Treatment of the running couplings",
     &VertexBase::_coupopt, 0, false, false);
  static SwitchOption interfaceCouplingRunning
    (interfaceCoupling,
     "Running",
     "Use the running couplings from the StandardModel object",
     0);
  static SwitchOption interfaceCouplingFixedSM
    (interfaceCoupling,
     "FixedSM",
     "Use the fixed values from the StandardModel object",
     1);
  static SwitchOption interfaceCouplingFixedLocal
    (interfaceCoupling,
     "FixedLocal",
     "Use the local fixed values",
     2);

  static Parameter<VertexBase,double> interfaceStrongCoupling
    ("StrongCoupling",
     "The fixed value of the strong coupling to use",
     &VertexBase::_gs, sqrt(4.*Constants::pi*0.3), 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<VertexBase,double> interfaceElectroMagneticCoupling
    ("ElectroMagneticCoupling",
     "The fixed value of the electromagnetic coupling to use",
     &VertexBase::_ee, sqrt(4.*Constants::pi/128.91), 0.0, 10.0,
     false, false, Interface::limited);

  static Parameter<VertexBase,double> interfaceSinThetaW
    ("SinThetaW",
     "The fixed value of sin theta_W to use",
     &VertexBase::_sw, sqrt(0.232), 0.0, 10.0,
     false, false, Interface::limited);

}

// find particles with a given id    
vector<long> VertexBase::search(unsigned int iloc,long idd) const {
  assert( iloc < _npoint );
  vector<long> out;
  for(const auto & pvec : _particles ) {
    bool found = pvec[iloc]->id() == idd;
    if(found) {
      for( tcPDPtr p : pvec ) {
	out.push_back(p->id());
      }
    }
  }
  return out;
}

// find particles with a given id    
vector<tPDPtr> VertexBase::search(unsigned int iloc,tcPDPtr idd) const {
  assert( iloc < _npoint );
  vector<tPDPtr> out;
  for(const auto & pvec : _particles) {
    if(pvec[iloc] == idd) {
      out.insert(out.end(),pvec.begin(), pvec.end());
    }
  }
  return out;
}

// check a given combination is allowed for a three point vertex
bool VertexBase::allowed(long ida, long idb, long idc, long idd) const {
  assert( ( _npoint==3 && idd == 0 ) || _npoint == 4 );
  vector<long> out = search(0,ida);
  for ( size_t ix = 0; ix < out.size(); ix += _npoint ) {
    if ( out[ix+1] == idb && out[ix+2] == idc
	 && ( idd == 0 || out[ix+3] == idd ) ) {
      return true; 
    }
  }
  return false;
}
  
// output the information
ostream & ThePEG::Helicity::operator<<(ostream & os, const VertexBase & in) {
  os << "Information on Vertex" << endl;
  os << "This is an " << in._npoint << " vertex\n";
  os << string( in._calckinematics ? 
		"The kinematic invariants are calculated" : 
		"The kinematics invariants are not calculated" ) << "\n";
  os << " Particles allowed for this Vertex\n";
  for(unsigned int ix=0;ix<in._particles.size();++ix) {
    for(unsigned int iy=0;iy<in._particles[ix].size();++iy) {
      os << in._particles[ix][iy]->PDGName() << "   ";
    }
    os << '\n';
  }
  return os;
}

// calculate the propagator for a diagram
Complex VertexBase::propagator(int iopt, Energy2 p2,tcPDPtr part,
			       complex<Energy> mass, complex<Energy> width) {
  if(mass.real() < ZERO) mass = part->mass();
  const complex<Energy2> mass2 = sqr(mass);

  if(width.real() < ZERO) {
    const tcWidthGeneratorPtr widthgen = part->widthGenerator();
    width = widthgen && (iopt==2 || iopt==6 ) ? 
      widthgen->width(*part,sqrt(p2)) : part->width();
  }
  const Complex ii(0.,1.);

  complex<Energy2> masswidth;
  if(iopt==5) {
    return Complex(UnitRemoval::E2/p2);
  }
  else if(iopt==4) {
    return 1.0;
  }
  else if(p2 < ZERO) {
    if(iopt!=7)
      masswidth = ZERO;
    else
      masswidth = ii * mass * width;
  }
  else {
    switch (iopt) {
    case 1: case 2: case 7:
      masswidth = ii * mass * width;         
      break;
    case 3: 
      masswidth = ZERO;                    
      break;
    case 6: 
      masswidth = ii * mass2 * width / sqrt(p2);
      return Complex(UnitRemoval::E2 * (mass2/p2) / (p2-mass2+masswidth));
    default:
      assert( false );
      return -999.999;
    }
  }
  return Complex(UnitRemoval::E2/(p2-mass2+masswidth));
}

void VertexBase::rebind(const TranslationMap & trans) {
  for (auto cit = _particles.begin(); cit != _particles.end(); ++cit) {
    for (auto cjt = cit->begin(); cjt != cit->end(); ++cjt) {
      *cjt = trans.translate(*cjt);
    }
  }
  set<tPDPtr> newinpart;
  for (auto it = _inpart.begin(); it != _inpart.end(); ++it) {
    newinpart.insert(trans.translate(*it));
  }
  _inpart = newinpart;
  set<tPDPtr> newoutpart;
  for (auto it = _outpart.begin(); it != _outpart.end(); ++it) {
    newoutpart.insert(trans.translate(*it));
  }
  _outpart = newoutpart;
  Interfaced::rebind(trans);
}

IVector VertexBase::getReferences() {
  IVector ret = Interfaced::getReferences();
  for (auto cit = _particles.begin(); cit != _particles.end(); ++cit) {
    for (auto cjt = cit->begin(); cjt != cit->end(); ++cjt) {
      ret.push_back(*cjt);
    }
  }
  for (auto it = _inpart.begin(); it != _inpart.end(); ++it) {
      ret.push_back(*it);
  }
  for (auto it = _outpart.begin(); it != _outpart.end(); ++it) {
      ret.push_back(*it);
  }
  return ret;
}


