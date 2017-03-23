// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the AbstractFFVVertex class.
//

#include "AbstractFFVVertex.h"
#include "ThePEG/Interface/ClassDocumentation.h"

using namespace Helicity;

AbstractNoPIOClassDescription<AbstractFFVVertex> 
AbstractFFVVertex::initAbstractFFVVertex;
// Definition of the static class description member.

void AbstractFFVVertex::Init() {

  static ClassDocumentation<AbstractFFVVertex> documentation
    ("The AbstractFFVVertex class provides a base class for all"
     " fermion-fermion-vector interactions");

}
SpinorWaveFunction 
AbstractFFVVertex::evaluateSmall(Energy2 q2,int iopt, tcPDPtr out,
				 const SpinorWaveFunction & sp1,
				 const VectorWaveFunction & vec3,
				 unsigned int , unsigned int ,
				 double , double , double ,
				 bool ,
				 SmallAngleDirection ,
				 Energy mass, Energy width) {
  return evaluate(q2,iopt,out,sp1,vec3,mass,width);
}

SpinorBarWaveFunction 
AbstractFFVVertex::evaluateSmall(Energy2 q2,int iopt, tcPDPtr out,
				 const SpinorBarWaveFunction & sbar2,
				 const VectorWaveFunction & vec3,
				 unsigned int , unsigned int ,
				 double , double , double ,
				 bool ,
				 SmallAngleDirection ,
				 Energy mass, Energy width) {
  return evaluate(q2,iopt,out,sbar2,vec3,mass,width);
}
