// -*- C++ -*-
//
// HoldFlag.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2011 Leif Lonnblad
//
// ThePEG is licenced under version 2 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_HoldFlag_H
#define ThePEG_HoldFlag_H
// This is the declaration of the HoldFlag class.

namespace ThePEG {

template <typename T = bool>
/**
 * <code>HoldFlag</code> objects are used to temporarily change the value
 * of an object, restoring the original value when the
 * <code>HoldFlag</code> object is destructed.
 *
 * @see Level
 */
class HoldFlag {

public:

  /**
   * Constructor setting a temporary value for the given object.
   * @param newFlag the object which value is temporarily changed.
   * @param holdFlag the temporary value for the newFlag object.
   */
  HoldFlag(T & newFlag, const T & holdFlag)
    : theFlag(newFlag), oldFlag(holdFlag) { std::swap(theFlag, oldFlag); }

  /**
   * Constructor setting the a temporary value for the given object.
   * @param newFlag the object which value is temporarily changed.
   * @param holdFlag the temporary value for the newFlag object.
   * @param finalFlag the newFlag object will be given the value
   * finalFlag when the HoldFlag object is destroyed.
   */
  HoldFlag(T & newFlag, const T & holdFlag, const T & finalFlag)
    : theFlag(newFlag), oldFlag(holdFlag) 
  {
    std::swap(theFlag, oldFlag);
    oldFlag = finalFlag;
  }

  /**
   * Destructor. Restores the corresponding object to its original
   * value.
   */
  ~HoldFlag() { std::swap(theFlag, oldFlag); }

private:

  /**
   * The object to be changed.
   */
  T & theFlag;

  /**
   * The value which will be restored when this is destroyed.
   */
  T oldFlag;

  /**
   * Default constructor is private and not implemented.
   */
  HoldFlag();

  /**
   * Copy constructor is private and not implemented.
   */
  HoldFlag(const HoldFlag &);

  /**
   * Assignment is private and not implemented.
   */
  HoldFlag & operator=(const HoldFlag &);

};

/**
 * Specialization of HoldFlag for boolean variables.
 */
template <>
class HoldFlag<bool> {

public:

  /**
   * Constructor setting the a temporary value for the bool variable.
   * @param newFlag the boolean variable which value is temporarily changed.
   * @param holdFlag the temporary value for the newFlag variable.
   */
  HoldFlag(bool & newFlag, bool holdFlag = true)
    : theFlag(newFlag), oldFlag(newFlag) { theFlag = holdFlag; }

  /**
   * Constructor setting the a temporary value for the bool variable.
   * @param newFlag the boolean variable which value is temporarily changed.
   * @param holdFlag the temporary value for the newFlag variable.
   * @param finalFlag the newFlag variable will be given the value
   * finalFlag when the HoldFlag object is destroyed.
   */
  HoldFlag(bool & newFlag, bool holdFlag, bool finalFlag)
    : theFlag(newFlag), oldFlag(finalFlag) { theFlag = holdFlag; }

  /**
   * Destructor. Restores the corresponding variable to its original
   * value.
   */
  ~HoldFlag() { theFlag = oldFlag; }

private:

  /**
   * The variable to be changed.
   */
  bool & theFlag;

  /**
   * The value which will be restored when this is destroyed.
   */
  bool oldFlag;

  /**
   * Default constructor is private and not implemented.
   */
  HoldFlag();

  /**
   * Copy constructor is private and not implemented.
   */
  HoldFlag(const HoldFlag &);

  /**
   * Assignment is private and not implemented.
   */
  HoldFlag & operator=(const HoldFlag &);

};

}

#endif /* ThePEG_HoldFlag_H */
