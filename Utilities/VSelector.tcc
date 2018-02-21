// -*- C++ -*-
//
// VSelector.tcc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//

namespace ThePEG {

template <typename T, typename WeightType>
WeightType VSelector<T,WeightType>::erase(const T & t) {
  theSum = WeightType();
  int j = 0;
  for ( int i = 0, N = theWeights.size(); i < N; ++i ) {
    if ( theObjects[i] == t ) continue;
    theSums[j] = (theSum += theWeights[i]);
    if ( i != j ) {
      theWeights[j] = theWeights[i];
      theObjects[j] = theObjects[i];
    }
    ++j;
  }
  theSums.erase(theSums.begin() + j, theSums.end());
  theWeights.erase(theWeights.begin() + j, theWeights.end());
  theObjects.erase(theObjects.begin() + j, theObjects.end());
  return theSum;
}

template <typename T, typename WeightType>
WeightType VSelector<T,WeightType>::reweight(WeightType d, const T & t) {
  d = max(d, WeightType());
  theSum = WeightType();
  for ( int i = 0, N = theWeights.size(); i < N; ++i ) {
    if ( theObjects[i] == t ) theWeights[i] = d;
    theSums[i] = ( theSum += theWeights[i] );
  }
  return theSum;
}

template <typename T, typename WeightType>
typename VSelector<T,WeightType>::size_type VSelector<T,WeightType>::
iselect(double rnd, double * remainder) const {
  if ( rnd <= 0 )
    throw range_error("Random number out of range in VSelector::select.");
  WeightType sum = rnd*theSum;
  WIterator it = upper_bound(theSums.begin(), theSums.end(), sum);
  if ( it == theSums.end() )
    throw range_error("Empty Selector, or random number out of range "
		      "in Selector::select");
  size_type i = it - theSums.begin();
  if ( remainder ) *remainder = 1.0 - (theSums[i] - sum)/theWeights[i];
  return i;
}

template <typename T, typename WeightType>
template <typename OStream>
void VSelector<T,WeightType>::output(OStream & os) const {
  os << theSum << theSums << theWeights << theObjects;
}

template <typename T, typename WeightType>
template <typename IStream>
void VSelector<T,WeightType>::input(IStream & is) {
  clear();
  is >> theSum >> theSums >> theWeights >> theObjects;
}

}
