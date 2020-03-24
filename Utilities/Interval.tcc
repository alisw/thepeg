// -*- C++ -*-
//
// Interval.tcc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2019 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//



namespace ThePEG {

template <typename T, typename CMP>
template <typename Iterator>
bool Interval<T,CMP>::check(Iterator first, Iterator last) {
  while ( first != last ) if ( includes(*first++) ) return true;
  return false;
}

template <typename T, typename CMP>
template <typename Iterator>
bool Interval<T,CMP>::checkAll(Iterator first, Iterator last) {
  while ( first != last ) if ( !includes(*first++) ) return false;
  return true;
}

template <typename T, typename CMP>
std::vector< Interval<T,CMP> > Interval<T,CMP>::split(Interval<T,CMP> i, T x) {
  std::vector< Interval<T,CMP> > intervals;
  Interval<T,CMP> low = i.chopLower(x);
  if ( low.check() ) intervals.push_back(low);
  intervals.push_back(x);
  return intervals;
}
  

template <typename T, typename CMP>
template <typename Iterator>
std::vector< Interval<T,CMP> >
Interval<T,CMP>::split(Interval<T,CMP> i, Iterator first, Iterator last) {
  typedef std::vector< Interval<T,CMP> > IVec;
  IVec intervals;
  if ( first != last ) {
    Interval<T,CMP> low = i.chopLower(*first++);
    IVec loints;
    if ( low.check() ) loints = split(low, first, last);
    IVec hiints = split(i, first, last);
    intervals.insert(intervals.end(), loints.begin(), loints.end());
    intervals.insert(intervals.end(), hiints.begin(), hiints.end());
  } else
    intervals = IVec(1, i);
  return intervals;
}

template <typename T, typename CMP>
CMP Interval<T,CMP>::cmp;



}
