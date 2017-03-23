// -*- C++ -*-
//
// AnyReference.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2007 Leif Lonnblad
// Copyright (C) 2010-2011 Simon Platzer
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_AnyReference_H
#define ThePEG_AnyReference_H

#include "ThePEG/Config/ThePEG.h"
#include "Exception.h"

namespace ThePEG {

  /**
   * AnyReference is inspired by boost::any to hold a reference to an
   * object of arbitrary type.
   */
  class AnyReference {

    struct ReferenceHolderBase {

      /**
       * The destructor.
       */
      virtual ~ReferenceHolderBase() {}

      /**
       * Return the type held.
       */
      virtual const std::type_info& type() const = 0;

      /**
       * Clone this reference holder.
       */
      virtual ReferenceHolderBase* clone() const = 0;

    };

    template<class T>
    struct ReferenceHolder 
      : public ReferenceHolderBase {

      /**
       * The reference held.
       */
      T& value;

      /**
       * Static member to initialize the reference.
       */
      static T& init() {
	static T v; return v;
      }

      /**
       * The default constructor
       */
      ReferenceHolder()
	: value(init()) {}

      /**
       * Construct from given reference.
       */
      explicit ReferenceHolder(T& v)
	: value(v) {}

      /**
       * Return the type held.
       */
      virtual const std::type_info& type() const {
	return typeid(T);
      }

      /**
       * Clone this reference holder.
       */
      virtual ReferenceHolderBase* clone() const {
	return new ReferenceHolder(*this);
      }

    };

    /**
     * The reference holder used.
     */
    ReferenceHolderBase* holder;

  public:

    /**
     * The default constructor
     */
    AnyReference()
      : holder(0) {}

    /**
     * The standard constructor
     */
    template<class T>
    AnyReference(T& v)
      : holder(new ReferenceHolder<T>(v)) {}

    /**
     * The copy constructor
     */
    AnyReference(const AnyReference& other)
      : holder(other.holder ? other.holder->clone() : 0) {}

    /**
     * The destructor.
     */
    ~AnyReference() {
    	delete holder; holder = 0;
    }

  public:

    /**
     * Return the type held.
     */
    const std::type_info& type() const { return holder ? holder->type() : typeid(void); }

    /**
     * Return true, if no reference is held.
     */
    bool empty() const { return !holder; }

  public:

    /**
     * Swap the references held.
     */
    AnyReference& swap(AnyReference& other) { 
      std::swap(holder,other.holder); return *this;
    }

    /**
     * Assign from definite type
     */
    template<class T>
    AnyReference& operator=(T& v) {
      AnyReference(v).swap(*this);
      return *this;
    }

    /**
     * Assign from AnyReference
     */
    AnyReference& operator=(AnyReference other) {
      other.swap(*this);
      return *this;
    }

    /**
     * Reset to not keep track of any reference.
     */
    void reset() {
      delete holder; holder = 0;
    }

  public:

    /**
     * Extract the held reference.
     */
    template<class T>
    T& cast() const {
      if ( !empty() && type() == typeid(T) ) {
	return static_cast<ReferenceHolder<T>*>(holder)->value;
      }
      throw Exception() << "Bad cast in AnyReference" << Exception::abortnow;
      return ReferenceHolder<T>::init();
    }

  };

}

#endif // ThePEG_AnyReference_H
