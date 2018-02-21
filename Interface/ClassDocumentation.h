// -*- C++ -*-
//
// ClassDocumentation.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef ThePEG_ClassDocumentation_H
#define ThePEG_ClassDocumentation_H
// This is the declaration of the ClassDocumentation class.

#include "ThePEG/Config/ThePEG.h"
#include "ClassDocumentation.fh"

namespace ThePEG {

/**
 * The ClassDocumentationBase class is used to communicate
 * documetation about an Interfaced class to the Repository.
 * Similarly to classes inheriting from InterfaceBase, only one static
 * object of the templated ClassDocumentation, which inherits from
 * ClassDocumentationBase, should be created for each Interfaced
 * class. This object will then automatically register itself with the
 * static Repository.
 *
 * The ClassDocumentationBase contains three strings with information
 * which are all specified in the constructor:
 *
 * The <i>documentation</i> contains a brief description of the
 * corresponding class which can be displayed by the Repository (or
 * user interfaces derived from it).
 *
 * The <i>model description</i> contains a very brief description of
 * the model of the process implemented in the step handler, given in
 * the form of a LaTeX \\item. This is written to a file after a run
 * by an EventGenerator.
 *
 * The <i>model references</i> contains possible LaTeX \\bibitems
 * corresponding to \\cite commands in the <i>model
 * description</i>. This is also written to a file after a run by an
 * EventGenerator.
 *
 * @see Interfaced
 * @see Repository
 * 
 */
class ClassDocumentationBase  {

protected:

  /**
   * The standard constructor can only be used from subclasses.
   * @param newDocumentation the <i>documentation</i> for the
   * corresponding class.
   * @param newModelDescription the <i>model description</i> for the
   * corresponding class.
   * @param newModelReferences the <i>model references</i> of the
   * corresponding class..
   * @param newTypeInfo the type_info object of the corresponding
   * class.
   */
  ClassDocumentationBase(string newDocumentation,
			 string newModelDescription,
			 string newModelReferences,
			 const type_info & newTypeInfo);

public:

  /**
   * The destructor.
   */
  virtual ~ClassDocumentationBase() {}

public:

  /**
   * Return the brief <i>documentation</i> of the corresponding class.
   */
  string documentation() const { return theDocumentation; }

  /**
   * Return the <i>model description</i> of the corresponding class.
   */
  string modelDescription() const { return theModelDescription; }

  /**
   * Return the <i>model references</i> of the corresponding class.
   */
  string modelReferences() const { return theModelReferences; }

private:

  /**
   * The brief <i>documentation</i> of the corresponding class.
   */
  string theDocumentation;

  /**
   * The <i>model description</i> of the corresponding class.
   */
  string theModelDescription;

  /**
   * The <i>model references</i> of the corresponding class.
   */
  string theModelReferences;

private:

  /**
   * Private and unimplemented default constructor.
   */
  ClassDocumentationBase();

  /**
   * Private and unimplemented copy constructor.
   */
  ClassDocumentationBase(const ClassDocumentationBase &);

  /**
   * Private and unimplemented assignment operator.
   */
  ClassDocumentationBase & operator=(const ClassDocumentationBase &);

};


/**
 * The <code>ClassDocumentation</code> class is used to communicate
 * documetation about an Interfaced class to the Repository.
 * Similarly to classes inheriting from InterfaceBase, only one static
 * object of the templated ClassDocumentation, which inherits from
 * ClassDocumentationBase, should be created for each Interfaced
 * class. This object will then automatically register itself with
 * the static Repository.
 *
 * The ClassDocumentation should in the constructor specify three
 * strings with information:
 *
 * The <i>documentation</i> contains a brief description of the
 * corresponding class which can be displayed by the Repository (or
 * user interfaces derived from it).
 *
 * The <i>model description</i> contains a very brief description of
 * the model of the process implemented in the step handler, given in
 * the form of a LaTeX \\item. This is written to a file after a run
 * by an EventGenerator.
 *
 * The <i>model references</i> contains possible LaTeX \\bibitems
 * corresponding to \\cite commands in the <i>model
 * description</i>. This is also written to a file after a run by an
 * EventGenerator.
 *
 * @see Interfaced
 * @see Repository
 * 
 */
template <typename T>
class ClassDocumentation: public ClassDocumentationBase {

public:

  /**
   * The standard constructor. All other constructors are private.
   * @param newDocumentation the <i>documentation</i> for the
   * corresponding class.
   * @param newModelDescription the <i>model description</i> for the
   * corresponding class.
   * @param newModelReferences the <i>model references</i> of the
   * corresponding class..
   */
  ClassDocumentation(string newDocumentation,
		     string newModelDescription = "",
		     string newModelReferences = "")
    : ClassDocumentationBase(newDocumentation, newModelDescription,
			     newModelReferences, typeid(T)) {}

private:

  /**
   * Private and unimplemented default constructor.
   */
  ClassDocumentation();

  /**
   * Private and unimplemented copy constructor.
   */
  ClassDocumentation(const ClassDocumentation &);

  /**
   * Private and unimplemented assignment operator.
   */
  ClassDocumentation & operator=(const ClassDocumentation &);

};

}

#endif /* ThePEG_ClassDocumentation_H */
