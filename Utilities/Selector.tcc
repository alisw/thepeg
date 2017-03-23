// -*- C++ -*-
//
// Selector.tcc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//

namespace ThePEG {

template <typename T, typename WeightType>
WeightType Selector<T,WeightType>::erase(const T & t) {
  Selector<T,WeightType> newSelector;
  WeightType oldsum = WeightType();
  for ( iterator it = theMap.begin();
	it != theMap.end(); ++it ) {
    WeightType r = it->first - oldsum;
    oldsum = it->first;
    if ( it->second != t ) newSelector.insert(r, it->second);
  }
  theMap.swap(newSelector.theMap);
  return theSum = newSelector.theSum;
}

template <typename T, typename WeightType>
const T & Selector<T,WeightType>::
select(double rnd, double * remainder) const throw(range_error) {
  if ( rnd <= 0 )
    throw range_error("Random number out of range in Selector::select.");
  const_iterator it = theMap.upper_bound(rnd*theSum);
  if ( it == theMap.end() )
    throw range_error("Empty Selector, or random number out of range "
		      "in Selector::select");
  if ( remainder ) {
    if ( it == theMap.begin() )
      *remainder = rnd*theSum/(it->first);
    else {
      const_iterator prit = it;
      --prit;
      *remainder = (rnd*theSum - prit->first)/(it->first - prit->first);
    }
  }
  return it->second;
}

template <typename T, typename WeightType>
T & Selector<T,WeightType>::
select(double rnd, double * remainder) throw(range_error) {
  if ( rnd <= 0 )
    throw range_error("Random number out of range in Selector::select.");
  iterator it = theMap.upper_bound(rnd*theSum);
  if ( it == theMap.end() )
    throw range_error("Empty Selector, or random number out of range "
		      "in Selector::select");
  if ( remainder ) {
    if ( it == theMap.begin() )
      *remainder = rnd*theSum/(it->first);
    else {
      const_iterator prit = it;
      --prit;
      *remainder = (rnd*theSum - prit->first)/(it->first - prit->first);
    }
  }
  return it->second;
}

template <typename T, typename WeightType>
template <typename OStream>
void Selector<T,WeightType>::output(OStream & os) const {
  os << ounit(theSum,TypeTraits<WeightType>::baseunit) << size();
  for ( const_iterator it = theMap.begin(); it != theMap.end(); ++it )
    os << ounit(it->first,TypeTraits<WeightType>::baseunit) << it->second;
}

template <typename T, typename WeightType>
template <typename IStream>
void Selector<T,WeightType>::input(IStream & is) {
  typedef typename MapType::value_type value_type;
  clear();
  T t;
  WeightType weightsum;
  long n;
  is >> iunit(theSum,TypeTraits<WeightType>::baseunit) >> n;
  while ( n-- ) {
    is >> iunit(weightsum,TypeTraits<WeightType>::baseunit) >> t;
    theMap.insert(theMap.end(), value_type(weightsum, t));
  }
}

}
