// -*- C++ -*-
//
// BaseRepository.tcc is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined templated member
// functions of the BaseRepository class.
//

namespace ThePEG {

template <typename T>
typename Ptr<T>::pointer BaseRepository::GetPtr(const T & t) {
  typedef typename Ptr<T>::pointer ptr;
  ObjectSet::iterator it = allObjects().find
    (IBPtr(const_cast<InterfacedBase *>(dynamic_cast<const InterfacedBase*>(&t))));
  return it == allObjects().end()? ptr():
    dynamic_ptr_cast<ptr>(*it);
}

template <typename PtrType>
PtrType BaseRepository::GetPtr(string name) {
  return dynamic_ptr_cast<PtrType>(GetPointer(name));
}

template <typename PtrType>
PtrType BaseRepository::GetObject(string name) {
  typedef ClassTraits<typename PtrTraits<PtrType>::value_type> Traits;
  IBPtr ip = GetPointer(name);
  if ( !ip ) throw RepositoryNotFound(name);
  PtrType p = dynamic_ptr_cast<PtrType>(ip);
  if ( !p ) throw RepositoryClassMisMatch(*ip, Traits::className());
  return p;
}

template <typename T>
typename Ptr<T>::pointer BaseRepository::clone(const T & t) {
  typedef typename Ptr<T>::pointer ptr;
  const InterfacedBase & ib = t;
  ptr ret;
  try {
    ret = dynamic_ptr_cast<ptr>(ib.clone());
  }
  catch ( ... ) {
    throw BadClone(t);
  }
  if ( !ret ) throw BadClassClone(t);
  return ret;
}

template <typename T>
typename Ptr<T>::pointer BaseRepository::fullclone(const T & t) {
  typedef typename Ptr<T>::pointer ptr;
  ptr ret;
  try {
    ret = dynamic_ptr_cast<ptr>(t.fullclone());
  }
  catch ( ... ) {
    throw BadClone(t);
  }
  if ( !ret ) throw BadClassClone(t);
  return ret;
}

template <typename Cont>
vector< pair<IBPtr, const InterfaceBase *> >
BaseRepository::getNonDefaultInterfaces(const Cont & c) {
  vector< pair<IBPtr, const InterfaceBase *> > ret;
  for ( typename Cont::const_iterator it = c.begin(); it != c.end(); ++it ) {
    InterfaceMap im = getInterfaces(typeid(**it));
    for ( InterfaceMap::iterator iit = im.begin(); iit != im.end(); ++iit )
      if ( iit->second->notDefault(**it) )
	ret.push_back(make_pair(*it, iit->second));
  }
  return ret;
}

}
