// -*- C++ -*-
//
// LoopGuard.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_LoopGuard_H
#define ThePEG_LoopGuard_H
// This is the declaration of the LoopGuard class.

namespace ThePEG {

/**
 * A LoopGuard object can be used to throw an exception if a loop is
 * iterated too many times. It is used by constructing an object
 * before the loop giving the maximum number of iterations allowed and
 * a message to be used as argument to the constructor of the
 * exception to be thrown. Inside the loop the parenthesis is called
 * without argument, which will increment and check an internal counter.
 * 
 *
 * @see Level
 */
template <typename ExceptionT = Exception,
          typename MessageT = const char *>
class LoopGuard {

public:

  /**
   * Create a loop guard object which will throw an exception of type
   * ExceptionT, constructed with 'mess' as argument, if the maximum
   * number of iterations is exceeded.
   */
  LoopGuard(const MessageT & mess, long maxc = 1000000 )
    : count(0), maxCount(maxc), message(mess) {}

  /**
   * Increase the iteration count and throw an ExceptionT if the
   * maximum number of iterations is exceeded.
   */
  void operator()()
  {
    if ( ++count > maxCount ) throw ExceptionT(message);
  }

private:

  /**
   * The number of counts so far.
   */
  long count;

  /**
   * The maximum number of counts allowed.
   */
  long maxCount;

  /**
   * The message with which the thrown ExceptionT object will be
   * initialized.
   */
  const MessageT & message;

private:

  /**
   * Default constructor is private and not implemented.
   */
  LoopGuard();

  /**
   * Copy constructor is private and not implemented.
   */
  LoopGuard(const LoopGuard &);

};

/**
 * A LoopGuard object can be used to throw an exception if a loop is
 * iterated too many times. It is used by constructing an object
 * before the loop giving the maximum number of iterations allowed and
 * a message to be used as argument to the constructor of the
 * exception to be thrown. Inside the loop the parenthesis is called
 * without argument, which will increment and check an internal
 * counter.  This specialization is for the case where the exception
 * class cannot be created with a message.
 * 
 *
 * @see Level
 */
template <typename ExceptionT>
class LoopGuard<ExceptionT,void> {

public:

  /**
   * Create a loop guard object which will throw an exception of type
   * ExceptionT, constructed with 'mess' as argument, if the maximum
   * number of iterations is exceeded.
   */
  LoopGuard(long maxc = 1000000 )
    : count(0), maxCount(maxc) {}

  /**
   * Increase the iteration count and throw an ExceptionT if the
   * maximum number of iterations is exceeded.
   */
  void operator()()
  {
    if ( ++count > maxCount ) throw ExceptionT();
  }

private:

  /**
   * The number of counts so far.
   */
  long count;

  /**
   * The maximum number of counts allowed.
   */
  long maxCount;

private:

  /**
   * Default constructor is private and not implemented.
   */
  LoopGuard();

  /**
   * Copy constructor is private and not implemented.
   */
  LoopGuard(const LoopGuard &);

};

}

#endif /* ThePEG_LoopGuard_H */
