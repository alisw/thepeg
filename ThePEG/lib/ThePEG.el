(require 'dired)

(defun ThePEG-fheader-file (class)
  "Create a .fh file with a skeleton suitable for a class CLASS."
  (interactive "sClass Name: ")
  (setq namespace (thepeg-get-namespace class))
  (setq class (thepeg-get-class class))
  (thepeg-fheaderfile namespace class))

(defun ThePEG-class-files (class)
  "Create .h and .cc files with skeletons suitable for a class
CLASS. The user will be prompted for the base class and the main include file.
The class may or may not be INTERFACED, PERSISTENT and/or CONCRETE."
  (interactive "sClass Name: ")
  (setq base (read-from-minibuffer "Base class name: " "HandlerBase"))
  (setq namespace (thepeg-get-namespace class))
  (setq class (thepeg-get-class class))
  (setq baseheader
	(cond ((string-equal base "")
	       "ThePEG/Config/ThePEG.h")
	      ((string-equal base "Base")
	       "ThePEG/Config/ThePEG.h")
	      ((string-equal base "PersistentBase") 
	       "ThePEG/Config/ThePEG.h")
	      ((string-equal base "Interfaced")
	       "ThePEG/Interface/Interfaced.h")
	      ((string-equal base "HandlerBase")
	       "ThePEG/Handlers/HandlerBase.h")
	      ((string-equal base "AnalysisHandler")
	       "ThePEG/Handlers/AnalysisHandler.h")
	      ((string-equal base "StepHandler")
	       "ThePEG/Handlers/StepHandler.h")
	      ((string-equal base "CascadeHandler")
	       "ThePEG/Handlers/CascadeHandler.h")
	      ((string-equal base "HadronizationHandler")
	       "ThePEG/Handlers/HadronizationHandler.h")
	      ((string-equal base "DecayHandler")
	       "ThePEG/Handlers/DecayHandler.h")
	      ((string-equal base "MEBase")
	       "ThePEG/MatrixElement/MEBase.h")
	      ((string-equal base "Decayer")
	       "ThePEG/PDT/Decayer.h")
	      ((string-equal base "PDFBase")
	       "ThePEG/PDF/PDFBase.h")
	      (t (concat namespace "/"
			 (car ( cdr (reverse
				     (split-string default-directory "/"))))
			 "/" base ".h"))))
  (setq baseheader (read-from-minibuffer "#include file for the base class: "
					 baseheader))
  (cond ((string-equal base "Interfaced") (setq interfaced t))
	((string-equal base "HandlerBase") (setq interfaced t))
	((string-equal base "StepHandler") (setq interfaced t))
	((string-equal base "CascadeHandler") (setq interfaced t))
	((string-equal base "HadronizationHandler") (setq interfaced t))
	((string-equal base "DecayHandler") (setq interfaced t))
	((string-equal base "MEBase") (setq interfaced t))
	((string-equal base "PDFBase") (setq interfaced t))
	((string-equal base "Decayer") (setq interfaced t))
	(t (setq interfaced (y-or-n-p "Will this class be interfaced "))))

  (setq declarations (cond ((string-equal base "PDFBase")
			    (thepeg-PDF-declare class base))
			   ((string-equal base "MEBase")
			    (thepeg-ME-declare class base))
			   ((string-equal base "StepHandler")
			    (thepeg-StepHandler-declare class base))
			   ((string-equal base "AnalysisHandler")
			    (thepeg-AnalysisHandler-declare class base))
			   ((string-equal base "Decayer")
			    (thepeg-Decayer-declare class base))
			   (t "")))
  (setq implementations (cond ((string-equal base "PDFBase")
			       (thepeg-PDF-implement class base))
			      ((string-equal base "MEBase")
			       (thepeg-ME-implement class base))
			      ((string-equal base "StepHandler")
			       (thepeg-StepHandler-implement class base))
			      ((string-equal base "AnalysisHandler")
			       (thepeg-AnalysisHandler-implement class base))
			      ((string-equal base "Decayer")
			       (thepeg-Decayer-implement class base))
			      (t "")))

  (setq persist (y-or-n-p "Will this class be persistent "))
  (setq concrete (y-or-n-p "Will this class be concrete "))
  (thepeg-sourcefile namespace class base persist interfaced concrete implementations)
;  (thepeg-iheaderfile namespace class base interfaced concrete)
;  (thepeg-fheaderfile namespace class)
  (thepeg-headerfile namespace class base baseheader interfaced
		persist concrete declarations))

(defun ThePEG-ME-class-files ()
  (interactive)
  (thepeg-specific-class-files "MEBase" "ThePEG/MatrixElement/MEBase.h"
			       'thepeg-ME-declare 'thepeg-ME-implement))

(defun ThePEG-PDF-class-files ()
  (interactive)
  (thepeg-specific-class-files "PDFBase" "ThePEG/PDF/PDFBase.h"
			       'thepeg-PDF-declare 'thepeg-PDF-implement))

(defun ThePEG-StepHandler-class-files ()
  (interactive)
  (thepeg-specific-class-files "StepHandler" "ThePEG/Handlers/StepHandler.h"
			       'thepeg-StepHandler-declare
			       'thepeg-StepHandler-implement))

(defun ThePEG-AnalysisHandler-class-files ()
  (interactive)
  (thepeg-specific-class-files "AnalysisHandler" "ThePEG/Handlers/AnalysisHandler.h"
			       'thepeg-AnalysisHandler-declare
			       'thepeg-AnalysisHandler-implement))

(defun ThePEG-Decayer-class-files ()
  (interactive)
  (thepeg-specific-class-files "Decayer" "ThePEG/PDT/Decayer.h"
			       'thepeg-Decayer-declare
			       'thepeg-Decayer-implement))

(defun thepeg-specific-class-files (base baseheader declfn implfn)
  (setq class (read-from-minibuffer "Class name: "))
  (setq base (read-from-minibuffer "Base class name: " base))
  (setq namespace (thepeg-get-namespace class))
  (setq class (thepeg-get-class class))
  (setq baseheader (read-from-minibuffer "#include file for the base class: "
					 baseheader))
  (setq interfaced t)
  (setq persist (y-or-n-p "Will this class be persistent "))
  (setq concrete (y-or-n-p "Will this class be concrete "))
  (setq declare (funcall declfn class base))

  (setq implement (funcall implfn class base))

  (thepeg-sourcefile namespace class base persist interfaced concrete implement)
;  (thepeg-iheaderfile namespace class base interfaced concrete)
;  (thepeg-fheaderfile namespace class)
  (thepeg-headerfile namespace class base baseheader
		interfaced persist concrete declare))

(defun thepeg-get-namespace (class)
  "Return the namespace part in an \"namepace::class\" string"
  (cond ((cdr (split-string class "::"))
	 (car (split-string class "::")))
	(t "ThePEG")))

(defun thepeg-get-class (class)
  "Return the namespace part in an \"namepace::class\" string"
  (cond ((cdr (split-string class "::"))
	 (car (cdr (split-string class "::"))))
	(t class)))

(defun thepeg-sourcefile (namespace class baseclass persistent interfaced concrete specialfn)
  "Create a file suitable for the implementation of a class CLASS.
The class may or may not be PERSISTENT and/or CONCRETE. SPECIALFN may be
used to include special function definitions"
  (find-file (concat class ".cc"))
  (c++-mode)
  (cond ((> (buffer-size) 0))
	(t (insert-string (thepeg-source namespace class baseclass persistent
					 interfaced concrete specialfn))
	   (beginning-of-buffer))))

(defun thepeg-source (namespace class baseclass persistent interfaced concrete specialfn)
  "Return a skeleton suitable for the implementation file of a class CLASS.
The class may or may not be PERSISTENT and/or CONCRETE. SPECIALFN may be
used to include special function definitions"

(setq interface (cond (interfaced (concat (cond (concrete "
IBPtr THECLASS::clone() const {
  return new_ptr(*this);
}

IBPtr THECLASS::fullclone() const {
  return new_ptr(*this);
}
")
						  (t ""))
					    "

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

"))
			(t "")))

  (setq piostring (cond (persistent "
void THECLASS::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
}

void THECLASS::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
}
")
			(t "")))

  (setq pioinclude (cond (persistent "
#include \"ThePEG/Persistency/PersistentOStream.h\"
#include \"ThePEG/Persistency/PersistentIStream.h\"")
			(t "")))

  (setq description (cond (persistent (cond (concrete "
DescribeClass<THECLASS,THEBASE>
  describeTHENAMESPACETHECLASS(\"THENAMESPACE::THECLASS\", \"THECLASS.so\");")
					    (t "
DescribeAbstractClass<THECLASS,THEBASE>
  describeTHENAMESPACETHECLASS(\"THENAMESPACE::THECLASS\", \"THECLASS.so\");")))
			  (t (cond (concrete "
DescribeNoPIOClass<THECLASS,THEBASE>
  describeTHENAMESPACETHECLASS(\"THENAMESPACE::THECLASS\", \"THECLASS.so\");")
				   (t "
DescribeAbstractNoPIOClass<THECLASS,THEBASE>
  describeTHENAMESPACETHECLASS(\"THENAMESPACE::THECLASS\", \"THECLASS.so\");")))))

  (thepeg-replace
   "THENAMESPACE" namespace
   (thepeg-replace
    "THECLASS" class
    (thepeg-replace
     "THEBASE" baseclass (concat "// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the THECLASS class.
//

#include \"THECLASS.h\"
#include \"ThePEG/Interface/ClassDocumentation.h\"
#include \"ThePEG/EventRecord/Particle.h\"
#include \"ThePEG/Repository/UseRandom.h\"
#include \"ThePEG/Repository/EventGenerator.h\"
#include \"ThePEG/Utilities/DescribeClass.h\"

" pioinclude "

using namespace " namespace ";

THECLASS::THECLASS() {}

THECLASS::~THECLASS() {}
" specialfn interface piostring "

// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found)." description "

void THECLASS::Init() {

  static ClassDocumentation<THECLASS> documentation
    (\"There is no documentation for the THECLASS class\");

}

")))))

(defun thepeg-iheaderfile (namespace class base interfaced concrete)
  "Create a icc-header file suitable for the inline function definitions
of a class CLASS which derives from the base class BASE. The class may or
may not be INTERFACED and/or CONCRETE."
  (find-file (concat class ".icc"))
  (c++-mode)
  (cond ((> (buffer-size) 0))
	(t (insert-string (thepeg-iheader namespace class base
					  interfaced concrete))
	   (beginning-of-buffer))))

(defconst thepeg-typetraits-decl "
#include \"ThePEG/Utilities/ClassTraits.h\"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of THECLASS. */
template <>
struct BaseClassTrait<THEFULLCLASS,1> {
  /** Typedef of the first base class of THECLASS. */
  typedef THEBASE NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the THECLASS class and the shared object where it is defined. */
template <>
struct ClassTraits<THEFULLCLASS>
  : public ClassTraitsBase<THEFULLCLASS> {
  /** Return a platform-independent class name */
  static string className() { return \"THEFULLCLASS\"; }
  /**
   * The name of a file containing the dynamic library where the class
   * THECLASS is implemented. It may also include several, space-separated,
   * libraries if the class THECLASS depends on other classes (base classes
   * excepted). In this case the listed libraries will be dynamically
   * linked in the order they are specified.
   */
  static string library() { return \"THECLASS.so\"; }
};

/** @endcond */

}
")

(defconst thepeg-interfaced-impl "
void THECLASS::doupdate() {
  THEBASE::doupdate();
  // First update base class.
  bool redo = touched();
  // redo if touched.
//  UpdateChecker::check(aDependentMember, redo);
  // Update referenced objects on which this depends redo is set to true
  // if the dependent object is touched.
//  for_each(ContainerOfDependencies, UpdateChecker(redo));
  // Update a container of references.
//  for_each(MapOfDependencies, UpdateMapChecker(redo));
  // Update a map of references.
  if ( !redo ) return;
  // return if nothing has been touched. Otherwise do the actual update.
//  touch()
  // Touch if anything has changed.
}

void THECLASS::doinit() {
  THEBASE::doinit();
}

void THECLASS::dofinish() {
  THEBASE::dofinish();
}

void THECLASS::doinitrun() {
  THEBASE::doinitrun();
}

void THECLASS::rebind(const TranslationMap & trans) {
  // dummy = trans.translate(dummy);
  THEBASE::rebind(trans);
}

IVector THECLASS::getReferences() {
  IVector ret = THEBASE::getReferences();
  // ret.push_back(dummy);
  return ret;
}
")

(defun thepeg-iheader (namespace class base interfaced concrete)
  "Return a skeleton suitable for the icc-header file for the inline
function definitions of a class CLASS which derives from the base
class BASE."

  (setq interface (cond (interfaced (concat (cond (concrete "
inline IBPtr THECLASS::clone() const {
  return new_ptr(*this);
}

inline IBPtr THECLASS::fullclone() const {
  return new_ptr(*this);
}
")
						  (t ""))
					    "

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).

"))
			(t "")))
  (thepeg-replace "THECLASS" class
		  (thepeg-replace "THEBASE" base (concat "// -*- C++ -*-
//
// This is the implementation of the inlined member functions of
// the THECLASS class.
//

namespace " namespace " {

inline THECLASS::THECLASS() {}

inline THECLASS::THECLASS(const THECLASS & x)
  : THEBASE(x) {}
" interface "
}
"))))

(defun thepeg-fheaderfile (namespace class)
  "Create a ifh-header file suitable for the forward declarations
of a class CLASS."
  (find-file (concat class ".fh"))
  (c++-mode)
  (cond ((> (buffer-size) 0))
	(t (insert-string (thepeg-fheader namespace class))
	   (beginning-of-buffer))))


(defun thepeg-fheader (namespace class)
  "Return a skeleton suitable for the fh-header file for the forward
declarations of a class CLASS."

  (setq NAMESPACE namespace)
  
  (thepeg-replace "THECLASS" class
		  (concat "// -*- C++ -*-
//
// This is the forward declaration of the THECLASS class.
//
#ifndef " NAMESPACE "_THECLASS_FH
#define " NAMESPACE "_THECLASS_FH

#include \"ThePEG/Config/ThePEG.h\"

namespace " namespace " {

class THECLASS;

ThePEG_DECLARE_POINTERS(" namespace "::THECLASS,THECLASSPtr);

}

#endif
")))

(defun thepeg-headerfile (namespace class base baseheader interfaced
				    persistent concrete specialfn)
  "Create a header file for a class CLASS
derived from a base class BASE which may be found in the include file
BASEHEADER. The class may or may not be INTERFACED, PERSISTENT and/or CONCRETE.
SPECIALFN may be used to include special function definitions"
  (find-file (concat class ".h"))
  (c++-mode)
  (cond ((> (buffer-size) 0))
	(t (insert-string (thepeg-header namespace class base baseheader
					 persistent concrete interfaced
					 specialfn))
	   (beginning-of-buffer))))

(defconst thepeg-persistent-decl "
  /** @name Functions used by the persistent I/O system. */
  //@{
  /**
   * Function used to write out object persistently.
   * @param os the persistent output stream written to.
   */
  void persistentOutput(PersistentOStream & os) const;

  /**
   * Function used to read in object persistently.
   * @param is the persistent input stream read from.
   * @param version the version number of the object when written.
   */
  void persistentInput(PersistentIStream & is, int version);
  //@}
")

(defconst thepeg-clone-decl "
protected:

  /** @name Clone Methods. */
  //@{
  /**
   * Make a simple clone of this object.
   * @return a pointer to the new object.
   */
  virtual IBPtr clone() const;

  /** Make a clone of this object, possibly modifying the cloned object
   * to make it sane.
   * @return a pointer to the new object.
   */
  virtual IBPtr fullclone() const;
  //@}
")

(defconst thepeg-interfaced-decl "
protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * Check sanity of the object during the setup phase.
   */
  virtual void doupdate();

  /**
   * Initialize this object after the setup phase before saving an
   * EventGenerator to disk.
   * @throws InitException if object could not be initialized properly.
   */
  virtual void doinit();

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  virtual void dofinish();

  /**
   * Rebind pointer to other Interfaced objects. Called in the setup phase
   * after all objects used in an EventGenerator has been cloned so that
   * the pointers will refer to the cloned objects afterwards.
   * @param trans a TranslationMap relating the original objects to
   * their respective clones.
   * @throws RebindException if no cloned object was found for a given
   * pointer.
   */
  virtual void rebind(const TranslationMap & trans);

  /**
   * Return a vector of all pointers to Interfaced objects used in this
   * object.
   * @return a vector of pointers.
   */
  virtual IVector getReferences();
  //@}
")

(defun ThePEG-persistent-decl ()
  "Insert declarations of the standard ThePEG persistent I/O functions."
  (interactive)
  (insert-string thepeg-persistent-decl))

(defun ThePEG-clone-decl ()
  "Insert declarations of the standard ThePEG clone functions."
  (interactive)
  (insert-string thepeg-clone-decl))

(defun ThePEG-interfaced-decl ()
  "Insert declarations of the standard ThePEG interfaced functions."
  (interactive)
  (insert-string thepeg-interfaced-decl))

(defun ThePEG-interfaced-impl ()
  "Insert default implementations of the standard ThePEG interfaced functions."
  (interactive)
  (setq class (read-from-minibuffer "Class Name: "
				    (file-name-sans-extension
				     (file-name-nondirectory
				      (buffer-file-name)))))
  (setq base (read-from-minibuffer "Base class Name: "))
  (insert-string (thepeg-replace "THEBASE" base
				 (thepeg-replace "THECLASS" class
						 thepeg-interfaced-impl))))

(defun ThePEG-typetraits-decl ()
  "Insert default declarations of the standard ThePEG type description system."
  (interactive)
  (setq fullclass (read-from-minibuffer "Class Name: "
					(file-name-sans-extension
					 (file-name-nondirectory
					  (buffer-file-name)))))
  (setq base (read-from-minibuffer "Base class Name: "))
  (setq class (thepeg-get-class fullclass))
  
  (insert-string (thepeg-replace
		  "THEFULLCLASS" fullclass
		  (thepeg-replace
		   "THEBASE" base
		   (thepeg-replace
		    "THECLASS" class
		    thepeg-typetraits-decl)))))

(defconst thepeg-type-desc "
// Static variable needed for the type description system in ThePEG.
#include \"ThePEG/Utilities/DescribeClass.h\"
DescribeABSTRACTNOPIOClass<THECLASS,THEBASE>
  describeTHENAMESPACETHECLASS(\"THENAMESPACE::THECLASS\", \"THELIBRARY.so\");
")

(defun ThePEG-type-desc ()
  "Insert default declarations of the standard ThePEG type description system."
  (interactive)
  (setq fullclass (read-from-minibuffer "Class Name: "
					(file-name-sans-extension
					 (file-name-nondirectory
					  (buffer-file-name)))))
  (setq base (read-from-minibuffer "Base class Name: "))
  (setq persist (y-or-n-p "Will this class be persistent "))
  (setq concrete (y-or-n-p "Will this class be concrete "))
  (setq mainlib (y-or-n-p "Will this class be in a main library file "))
  (setq class (thepeg-get-class fullclass))
  (setq namespace (thepeg-get-namespace fullclass))
  (setq abstract (cond (concrete "") (t "Abstract")))
  (setq nopio (cond (persist "") (t "NoPIO")))
  (setq library (cond (mainlib (concat "lib" namespace)) (t class)))

  (insert-string (thepeg-replace
		  "THELIBRARY" library
		  (thepeg-replace
		   "ABSTRACT" abstract
		   (thepeg-replace
		    "NOPIO" nopio
		    (thepeg-replace
		     "THENAMESPACE" namespace
		     (thepeg-replace
		      "THEBASE" base
		      (thepeg-replace
		       "THECLASS" class
		       thepeg-type-desc))))))))
		 
(defun thepeg-header (namespace class base baseheader persistent
				concrete interfaced specialfn)
  "Return a skeleton suitable for the header file of a class CLASS
derived from a base class BASE which may be found in the include file
BASEHEADER. The class may or may not be PERSISTENT, INTERFACED and/or CONCRETE.
SPECIALFN may be used to include special function definitions"

  (setq piostring (cond (persistent thepeg-persistent-decl)
			(t "")))
  (setq cloning (cond ((and concrete interfaced) thepeg-clone-decl)
		      (t "")))
  (setq interface (cond (interfaced "

// If needed, insert declarations of virtual function defined in the
// InterfacedBase class here (using ThePEG-interfaced-decl in Emacs).

")
			(t "")))
  (setq doxygen-see-interfaces (cond (interfaced (concat "
 *
 * @see \\ref " class "Interfaces \"The interfaces\"
 * defined for " class "."))
			(t "")))
  (setq include (cond ((string-equal baseheader "") "ThePEG/Config/ThePEG.h")
		      (t baseheader)))
  (setq special (cond ((string-equal specialfn "") "")
		      (t (concat "
public:
" specialfn))))
  
  (setq using (cond ((string-match namespace "ThePEG") "")
		    (t "
using namespace ThePEG;
")))
  
  (setq NAMESPACE  namespace)
  
  (setq namespacequalifyer (cond ((string-match namespace "ThePEG") "")
				 (t (concat namespace "::"))))
  
  (setq basedeclare (cond ((string-equal base "") "")
			  (t (concat ": public " base))))
  
  (thepeg-replace "THECLASS" class (concat "// -*- C++ -*-
#ifndef " NAMESPACE "_THECLASS_H
#define " NAMESPACE "_THECLASS_H
//
// This is the declaration of the THECLASS class.
//

#include \"" include "\"

namespace " namespace " {
" using "
/**
 * Here is the documentation of the THECLASS class." doxygen-see-interfaces "
 */
class THECLASS" basedeclare " {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  THECLASS();

  /**
   * The destructor.
   */
  virtual ~THECLASS();
  //@}
" specialfn "
public:
" piostring "
  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();
" cloning interface "
private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  THECLASS & operator=(const THECLASS &);

};

}

#endif /* " NAMESPACE "_THECLASS_H */
")))

(defun thepeg-replace (regexp newtext string)
  (dired-string-replace-match regexp string newtext t t))

(defun thepeg-PDF-implement (class base)
(concat "

#include \"ThePEG/PDT/ParticleData.h\"
#include \"ThePEG/PDT/EnumParticles.h\"

bool " class "::canHandleParticle(tcPDPtr particle) const {
  // We assume that all protons and neutrons can be handled
  return ( abs(particle->id()) == abs(long(ParticleID::pplus)) ||
	   abs(particle->id()) == abs(long(ParticleID::n0)) );
}

cPDVector " class "::partons(tcPDPtr particle) const {
  // We assume that all standard partons can be extracted.
  cPDVector ret;
  if ( canHandleParticle(particle) ) {
    ret.push_back(getParticleData(ParticleID::g));
    for ( int i = 1; i <= 5; ++i ) {
      ret.push_back(getParticleData(i));
      ret.push_back(getParticleData(-i));
    }
  }
  return ret;
}

double " class "::xfl(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
                      double l, Energy2 particleScale) const {
  // Here we should return the actual density.
  return 0.0;
}

double " class "::xfvl(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
		     double l, Energy2 particleScale) const {
  // Here we should return the actual valence density.
  return 0.0;
}

"))

(defun thepeg-ME-implement (class base)
(concat "

#include \"ThePEG/PDT/EnumParticles.h\"
#include \"ThePEG/MatrixElement/Tree2toNDiagram.h\"

void " class "::getDiagrams() const {
  // Here is an example on how to specify diagrams.

  tcPDPtr g = getParticleData(ParticleID::g);
  for ( int i = 1; i <= 5; ++i ) {
    tcPDPtr q = getParticleData(i);
    tcPDPtr qb = q->CC();

    // For each flavour we add:
    add(new_ptr((Tree2toNDiagram(3), q, qb, qb, 1, g, 2, g, -1)));
    // t-channel q + qbar -> g + g
    add(new_ptr((Tree2toNDiagram(3), q, qb, qb, 2, g, 1, g, -2)));
    // u-channel q + qbar -> g + g
    add(new_ptr((Tree2toNDiagram(2), q, qb, 1, g , 3, g, 3, g, -3)));
    // s-channel q + qbar -> g + g
  }
}

Energy2 " class "::scale() const {
  return sHat();
}

int " class "::nDim() const {
  return 1;
}

void " class "::setKinematics() {
  " base "::setKinematics(); // Always call the base class method first.
}

bool " class "::generateKinematics(const double * r) {
  // Here you can use nDim() random numbers in the vector provided
  // to generate the internal kinematics. Note that sHat() has
  // already been given from the outside.

  // Save the jacobian dPS/dr for later use.
  jacobian(1.0);

  return false; // Return false if matrix element was zero.
}

double " class "::me2() const {
  return 0.0;
}

CrossSection " class "::dSigHatDR() const {
  return sqr(hbarc)*me2()*jacobian()/sHat(); // Here we can add other prefactors
                                             // coming from the phase space
                                             // integration.
}

unsigned int " class "::orderInAlphaS() const {
  return 0;
}

unsigned int " class "::orderInAlphaEW() const {
  return 0;
}

Selector<MEBase::DiagramIndex>
" class "::diagrams(const DiagramVector & diags) const {
  // This example corresponds to the diagrams specified in the example
  // in the getDiagrams() function.

  Selector<DiagramIndex> sel;
  for ( DiagramIndex i = 0; i < diags.size(); ++i ) 
    if ( diags[i]->id() == -1 ) sel.insert(1.0, i);
    else if ( diags[i]->id() == -2 )  sel.insert(1.0, i);
    else if ( diags[i]->id() == -3 )  sel.insert(1.0, i);
  // You probably do not want equal weights here...
  return sel;

  // If there is only one possible diagram you can override the
  // MEBase::diagram function instead.

}

Selector<const ColourLines *>
" class "::colourGeometries(tcDiagPtr diag) const {
  // This example corresponds to the diagrams specified in the example
  // in the getDiagrams() function.

  static ColourLines ctST(\"1 4, -4 -2 5, -5 -3\");
  static ColourLines ctSU(\"1 5, -5 -2 4, -4 -3\");

  Selector<const ColourLines *> sel;
  if ( diag->id() == -1 || diag->id() == -3 )
    sel.insert(1.0, &ctST);
  else
    sel.insert(1.0, &ctSU);
  return sel;

  // If there is only one possible colour geometry you can override the
  // MEBase::selectColourGeometry function instead.

}

"))

(defun thepeg-PDF-declare (class base)
  (concat "
public:

  /** @name Virtual functions to be overridden by sub-classes. */
  //@{
  /**
   * Return true if this PDF can handle the extraction of partons from
   * the given \\a particle.
   */
  virtual bool canHandleParticle(tcPDPtr particle) const;

  /**
   * Return the partons which this PDF may extract from the given
   * \\a particle.
   */
  virtual cPDVector partons(tcPDPtr particle) const;

  /**
   * The density. Return the pdf for the given \\a parton inside the
   * given \\a particle for the virtuality \\a partonScale and
   * logarithmic momentum fraction \\a l \\f$(l=\\log(1/x)\\$f. The \\a
   * particle is assumed to have a virtuality \\a particleScale.
   */
  virtual double xfl(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
		     double l, Energy2 particleScale = 0.0*GeV2) const;

  /**
   * The valence density. Return the pdf for the given cvalence \\a
   * parton inside the given \\a particle for the virtuality \\a
   * partonScale and logarithmic momentum fraction \\a l
   * \\f$(l=\\log(1/x)\\$f. The \\a particle is assumed to have a
   * virtuality \\a particleScale. If not overidden by a sub class this
   * will return zero.
   */
  virtual double xfvl(tcPDPtr particle, tcPDPtr parton, Energy2 partonScale,
		     double l, Energy2 particleScale = 0.0*GeV2) const;
  //@}

"))

(defun thepeg-ME-declare (class base)
  (concat "
public:

  /** @name Virtual functions required by the MEBase class. */
  //@{
  /**
   * Return the order in \\f$\\alpha_S\\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInAlphaS() const;

  /**
   * Return the order in \\f$\\alpha_{EW}\\f$ in which this matrix
   * element is given.
   */
  virtual unsigned int orderInAlphaEW() const;

  /**
   * The matrix element for the kinematical configuration
   * previously provided by the last call to setKinematics(), suitably
   * scaled by sHat() to give a dimension-less number.
   * @return the matrix element scaled with sHat() to give a
   * dimensionless number.
   */
  virtual double me2() const;

  /**
   * Return the scale associated with the last set phase space point.
   */
  virtual Energy2 scale() const;

  /**
   * Set the typed and momenta of the incoming and outgoing partons to
   * be used in subsequent calls to me() and colourGeometries()
   * according to the associated XComb object. If the function is
   * overridden in a sub class the new function must call the base
   * class one first.
   */
  virtual void setKinematics();

  /**
   * The number of internal degrees of freedom used in the matrix
   * element.
   */
  virtual int nDim() const;

  /**
   * Generate internal degrees of freedom given nDim() uniform
   * random numbers in the interval \\f$ ]0,1[ \\f$. To help the phase space
   * generator, the dSigHatDR should be a smooth function of these
   * numbers, although this is not strictly necessary.
   * @param r a pointer to the first of nDim() consecutive random numbers.
   * @return true if the generation succeeded, otherwise false.
   */
  virtual bool generateKinematics(const double * r);

  /**
   * Return the matrix element squared differential in the variables
   * given by the last call to generateKinematics().
   */
  virtual CrossSection dSigHatDR() const;

  /**
   * Add all possible diagrams with the add() function.
   */
  virtual void getDiagrams() const;

  /**
   * Get diagram selector. With the information previously supplied with the
   * setKinematics method, a derived class may optionally
   * override this method to weight the given diagrams with their
   * (although certainly not physical) relative probabilities.
   * @param dv the diagrams to be weighted.
   * @return a Selector relating the given diagrams to their weights.
   */
  virtual Selector<DiagramIndex> diagrams(const DiagramVector & dv) const;

  /**
   * Return a Selector with possible colour geometries for the selected
   * diagram weighted by their relative probabilities.
   * @param diag the diagram chosen.
   * @return the possible colour geometries weighted by their
   * relative probabilities.
   */
  virtual Selector<const ColourLines *>
  colourGeometries(tcDiagPtr diag) const;
  //@}

"))

(defun thepeg-StepHandler-declare (class base)
  (concat "
public:

  /** @name Virtual functions required by the StepHandler class. */
  //@{
  /**
    * The main function called by the EventHandler class to
    * perform a step. Given the current state of an Event, this function
    * performs the event generation step and includes the result in a new
    * Step object int the Event record.
    * @param eh the EventHandler in charge of the Event generation.
    * @param tagged if not empty these are the only particles which should
    * be considered by the StepHandler.
    * @param hint a Hint object with possible information from previously
    * performed steps.
    * @throws Veto if the StepHandler requires the current step to be discarded.
    * @throws Stop if the generation of the current Event should be stopped
    * after this call.
    * @throws Exception if something goes wrong.
    */
  virtual void handle(EventHandler & eh, const tPVector & tagged,
		      const Hint & hint);
  //@}
"))

(defun thepeg-StepHandler-implement (class base)
  (concat "
void " class "::
handle(EventHandler & eh, const tPVector & tagged,
       const Hint & hint) {
  // Implement the Handle method here.
  // Note that if the method actually does anything to the current event
  // the changes should be inserted in a new step which should be obtained
  // by 'ch.newStep()'.
  // Note also that the general advice is to only consider the particles in
  // the 'tagged' vector.
}
"))

(defun thepeg-AnalysisHandler-declare (class base)
  (concat "
public:

  /** @name Virtual functions required by the AnalysisHandler class. */
  //@{
  /**
   * Analyze a given Event. Note that a fully generated event
   * may be presented several times, if it has been manipulated in
   * between. The default version of this function will call transform
   * to make a lorentz transformation of the whole event, then extract
   * all final state particles and call analyze(tPVector) of this
   * analysis object and those of all associated analysis objects. The
   * default version will not, however, do anything on events which
   * have not been fully generated, or have been manipulated in any
   * way.
   * @param event pointer to the Event to be analyzed.
   * @param ieve the event number.
   * @param loop the number of times this event has been presented.
   * If negative the event is now fully generated.
   * @param state a number different from zero if the event has been
   * manipulated in some way since it was last presented.
   */
  virtual void analyze(tEventPtr event, long ieve, int loop, int state);

  /**
   * Return a LorentzTransform which would put the event in the
   * desired Lorentz frame.
   * @param event a pointer to the Event to be considered.
   * @return the LorentzRotation used in the transformation.
   */
  virtual LorentzRotation transform(tcEventPtr event) const;

  /**
   * Analyze the given vector of particles. The default version calls
   * analyze(tPPtr) for each of the particles.
   * @param particles the vector of pointers to particles to be analyzed
   * @param weight the weight of the current event.
   */
  virtual void analyze(const tPVector & particles, double weight);

  /**
   * Analyze the given particle.
   * @param particle pointer to the particle to be analyzed.
   * @param weight the weight of the current event.
   */
  virtual void analyze(tPPtr particle, double weight);
  //@}

protected:

  /**
   * Initialize this object. Called in the run phase just before
   * a run begins.
   */
  virtual void doinitrun();

  /**
   * Finalize this object. Called in the run phase just after a
   * run has ended. Used eg. to write out statistics.
   */
  virtual void dofinish();

"))

(defun thepeg-AnalysisHandler-implement (class base)
  (concat "


#ifndef LWH_AIAnalysisFactory_H
#ifndef LWH 
#define LWH ThePEGLWH
#endif
#include \"ThePEG/Analysis/LWH/AnalysisFactory.h\"
#endif

void " class "::analyze(tEventPtr event, long ieve, int loop, int state) {
  " base "::analyze(event, ieve, loop, state);
  // Rotate to CMS, extract final state particles and call analyze(particles).
}

LorentzRotation " class "::transform(tcEventPtr event) const {
  return LorentzRotation();
  // Return the Rotation to the frame in which you want to perform the analysis.
}

void " class "::analyze(const tPVector & particles, double weight) {
  " base "::analyze(particles);
  // Calls analyze() for each particle.
}

void " class "::analyze(tPPtr, double weight) {}

void " class "::dofinish() {
  " base "::dofinish();
  // *** ATTENTION *** Normalize and post-process histograms here.
}

void " class "::doinitrun() {
  " base "::doinitrun();
  // *** ATTENTION *** histogramFactory().registerClient(this); // Initialize histograms.
  // *** ATTENTION *** histogramFactory().mkdirs(\"/SomeDir\"); // Put histograms in specal directory.
}

"))

(defun thepeg-Decayer-declare (class base)
  (concat "
public:

  /** @name Virtual functions required by the Decayer class. */
  //@{
  /**
   * Check if this decayer can perfom the decay specified by the
   * given decay mode.
   * @param dm the DecayMode describing the decay.
   * @return true if this decayer can handle the given mode, otherwise false.
   */
  virtual bool accept(const DecayMode & dm) const;

  /**
   * Perform a decay for a given DecayMode and a given Particle instance.
   * @param dm the DecayMode describing the decay.
   * @param p the Particle instance to be decayed.
   * @return a ParticleVector containing the decay products.
   */
  virtual ParticleVector decay(const DecayMode & dm, const Particle & p) const;
  //@}
"))

(defun thepeg-Decayer-implement (class base)
  (concat "
bool " class "::accept(const DecayMode & dm) const {
  return false;
}

ParticleVector " class "::decay(const DecayMode & dm,
				  const Particle & parent) const {
  ParticleVector children = dm.produceProducts();
  return children;
}

"))

(defun ThePEG-see-interfaces ()
  "Create a doxygen comment referring to the interfaces of a class."
  (interactive)
  (setq class (read-from-minibuffer "Class Name: "
				    (file-name-sans-extension
				     (file-name-nondirectory
				      (buffer-file-name)))))
  (insert-string (concat " * @see \\ref " class "Interfaces \"The interfaces\"
 * defined for " class ".
")))
  
(defun ThePEG-command ()
  "Create a Command variable suitable for inclusion in an Init() function."
  (interactive)
  (setq class (read-from-minibuffer "Class Name: "
				    (file-name-sans-extension
				     (file-name-nondirectory
				      (buffer-file-name)))))
  (setq name (read-from-minibuffer "Name: "))
  (setq desc (read-from-minibuffer "Description: "))
  (setq memb (read-from-minibuffer "Member function: "
				   (concat "&" class "::do" name)))
  (setq safe (y-or-n-p "Is this parameter dependency safe? "))
  (insert-string (concat "
  static Command<" class "> interface" name "
    (\"" name "\",
     \"" desc "\",
     " memb ", " (cond (safe "true")(t "false")) ");
")))

(defun ThePEG-parameter ()
  "Create a Parameter variable suitable for inclusion in an Init() function."
  (interactive)
  (setq class (read-from-minibuffer "Class Name: "
				    (file-name-sans-extension
				     (file-name-nondirectory
				      (buffer-file-name)))))
  (setq type (read-from-minibuffer "Type: " "double"))
  (setq isstring (cond ((string-equal type "string") t) (t nil)))
  (setq name (read-from-minibuffer "Name: "))
  (setq desc (read-from-minibuffer "Description: "))
  (setq memb (read-from-minibuffer "Member variable: "
				   (concat "&" class "::the" name)))
  (setq hasu (cond (isstring nil)
		   (t (y-or-n-p "Does this parameter have a unit? "))))
  (cond (hasu (setq unit (read-from-minibuffer "Unit: " "GeV"))
	      (setq stru (concat "*" unit)))
	(t    (setq stru "")))
  (setq defa (read-from-minibuffer "Default value: " (concat "1.0" stru)))
  (setq limlo (cond (isstring nil)
		    (t (y-or-n-p "Is this parameter limited from below? "))))
  (cond (limlo (setq mini (read-from-minibuffer "Minimum value: "
						(concat "0.0" stru))))
	(t    (setq mini (concat "0" stru))))
  (setq limup (cond (isstring nil)
		    (t (y-or-n-p "Is this parameter limited from above? "))))
  (cond (limup (setq maxi (read-from-minibuffer "Maximum value: "
						(concat "10.0" stru))))
	(t    (setq maxi (concat "0" stru))))
  (setq safe (y-or-n-p "Is this parameter dependency safe? "))
  (setq ronl (y-or-n-p "Is this parameter read-only? "))
  (insert-string (concat "
  static Parameter<" class "," type "> interface" name "
    (\"" name "\",
     \"" desc "\",
     " memb ", " (cond (hasu (concat unit ", "))(t "")) defa
     (cond (isstring "")(t (concat ", " mini ", " maxi)))
     ",
     " (cond (safe "true")(t "false")) ", " (cond (ronl "true")(t "false"))
     (cond (isstring "")
	   ((and limup limlo) ", Interface::limited")
	   (limup ", Interface::upperlim")
	   (limlo ", Interface::lowerlim")
	   (t ", Interface::nolimits"))))
  (cond ((y-or-n-p "Are there any set/get functions? ")
	 (insert-string (concat ",
     " (thepeg-ifn "Set-function: " (concat "&" class "::set" name)
		   (concat "(void(" class "::*)(" type "))(0)"))
     ", "
     (thepeg-ifn "Get-function: " (concat "&" class "::get" name)
		   (concat "(" type "(" class "::*)()const)(0)"))
     (cond (isstring "")
	   (t
	    (concat ",
     "
		    (thepeg-ifn "Min-function: "
				(concat "&" class "::min" name)
				(concat "(" type "(" class "::*)()const)(0)"))
		    ", "
		    (thepeg-ifn "Max-function: "
				(concat "&" class "::max" name)
				(concat "(" type "(" class "::*)()const)(0)"))
		    )))
     ",
     "
     (thepeg-ifn "Default-function: "
			   (concat "&" class "::def" name)
			   (concat "(" type "(" class "::*)()const)(0)"))))))
  (insert-string ");
"))

(defun ThePEG-parvector ()
  "Create a ParVector variable suitable for inclusion in an Init() function."
  (interactive)
  (setq class (read-from-minibuffer "Class Name: "
				    (file-name-sans-extension
				     (file-name-nondirectory
				      (buffer-file-name)))))
  (setq type (read-from-minibuffer "Type: " "double"))
  (setq name (read-from-minibuffer "Name: "))
  (setq desc (read-from-minibuffer "Description: "))
  (setq memb (read-from-minibuffer "Member variable: "
				   (concat "&" class "::the" name)))
  (setq hasu (y-or-n-p "Does the parameters in this vector have a unit? "))
  (cond (hasu (setq unit (read-from-minibuffer "Unit: " "GeV"))
	      (setq stru (concat "*" unit)))
	(t    (setq stru "")))
  (setq size (read-from-minibuffer "Size of vector (varying if < 0): " "-1"))
  (setq defa (read-from-minibuffer "Default value: " (concat "1.0" stru)))
  (setq limlo (y-or-n-p "Is this parameter limited from below? "))
  (cond (limlo (setq mini (read-from-minibuffer "Minimum value: "
						(concat "0.0" stru))))
	(t    (setq mini (concat "0" stru))))
  (setq limup (y-or-n-p "Is this parameter limited from above? "))
  (cond (limup (setq maxi (read-from-minibuffer "Maximum value: "
						(concat "10.0" stru))))
	(t    (setq maxi (concat "0" stru))))
  (setq safe (y-or-n-p "Is this parameter vector dependency safe? "))
  (setq ronl (y-or-n-p "Are the parameters read-only? "))
  (insert-string (concat "
  static ParVector<" class "," type "> interface" name "
    (\"" name "\",
     \"" desc "\",
     " memb ", " (cond (hasu (concat unit ", "))(t "")) size ", " defa ", "
     mini ", " maxi ",
     " (cond (safe "true")(t "false")) ", " (cond (ronl "true")(t "false"))
     (cond ((and limup limlo) ", Interface::limited")
	   (limup ", Interface::lowerlim")
	   (limlo ", Interface::upperlim")
	   (t ", Interface::nolimits"))))
  (cond ((y-or-n-p "Are there any set/get functions? ")
	 (insert-string (concat ",
     " (thepeg-ifn "Set-function" (concat "&" class "::set" name)
		   (concat "(void(" class "::*)(" type ",int))(0)"))
     ", "
     (thepeg-ifn "Insert-function: " (concat "&" class "::ins" name)
		 (concat "(void(" class "::*)(" type ",int))(0)"))
     ",
     "
     (thepeg-ifn "Remove-function: " (concat "&" class "::del" name)
		 (concat "(void(" class "::*)(int))(0)"))
     ", "
     (thepeg-ifn "Get-function: " (concat "&" class "::get" name)
		 (concat "(vector<" type ">(" class "::*)()const)(0)"))
     ",
     "
     (thepeg-ifn "Default-function: " (concat "&" class "::def" name)
		 (concat "(" type "(" class "::*)(int)const)(0)"))
     ", "
     (thepeg-ifn "Min-function: " (concat "&" class "::min" name)
		 (concat "(" type "(" class "::*)(int)const)(0)"))
     ",
     "
     (thepeg-ifn "Max-function: " (concat "&" class "::max" name)
		 (concat "(" type "(" class "::*)(int)const)(0)"))
     ", "
     (thepeg-ifn "Strin-get-function: " (concat "&" class "::max" name)
		 (concat "(vector<string>(" class "::*)()const)(0)"))))))
  (insert-string ");
"))

(defun thepeg-ifn (desc def zero)
  (cond ((y-or-n-p (concat "Is there a " desc "?"))
	 (read-from-minibuffer (concat desc ": ") def))
	(t zero)));


(defun ThePEG-reference ()
  "Create a Reference variable suitable for inclusion in an Init() function."
  (interactive)
  (setq class (read-from-minibuffer "Class Name: "
				    (file-name-sans-extension
				     (file-name-nondirectory
				      (buffer-file-name)))))
  (setq type (read-from-minibuffer "Type: " "Interfaced"))
  (setq name (read-from-minibuffer "Name: "))
  (setq desc (read-from-minibuffer "Description: "))
  (setq memb (read-from-minibuffer "Member variable: "
				   (concat "&" class "::the" name)))
  (setq safe (y-or-n-p "Is this reference dependency safe? "))
  (setq ronl (y-or-n-p "Is this reference read-only? "))
  (setq rebi (y-or-n-p "Should this reference be automatically rebound? "))
  (setq null (y-or-n-p "Can this reference be null? "))
  (setq defn (y-or-n-p "Should this reference be given a default if null? "))

  (insert-string (concat "
  static Reference<" class "," type "> interface" name "
    (\"" name "\",
     \"" desc "\",
     " memb ", "
     (cond (safe "true")(t "false")) ", " (cond (ronl "true")(t "false")) ", "
     (cond (rebi "true")(t "false")) ", " (cond (null "true")(t "false")) ", "
     (cond (defn "true")(t "false"))))
  (cond ((y-or-n-p "Are there any set/get functions? ")
	 (insert-string (concat ",
     " (thepeg-ifn "Set-function: " (concat "&" class "::set" name)
		   (concat "(void(" class "::*)(Ptr<" type ">::pointer))(0)"))
     ",
     "
     (thepeg-ifn "Get-function: " (concat "&" class "::get" name)
		 (concat "(Ptr<" type ">::pointer(" class "::*)()const)(0)"))
     ",
     "
     (thepeg-ifn "Check-function: "
		 (concat "&" class "::check" name)
		 (concat "(bool(" class "::*)(Ptr<"
			 type ">::const_pointer)const)(0)"))))))
  (insert-string ");
"))

(defun ThePEG-refvector ()
  "Create a RefVector variable suitable for inclusion in an Init() function."
  (interactive)
  (setq class (read-from-minibuffer "Class Name: "
				    (file-name-sans-extension
				     (file-name-nondirectory
				      (buffer-file-name)))))
  (setq type (read-from-minibuffer "Type: " "Interfaced"))
  (setq name (read-from-minibuffer "Name: "))
  (setq desc (read-from-minibuffer "Description: "))
  (setq memb (read-from-minibuffer "Member variable: "
				   (concat "&" class "::the" name)))
  (setq size (read-from-minibuffer "Size of vector (varying if < 0): " "-1"))
  (setq safe (y-or-n-p "Is this reference vector dependency safe? "))
  (setq ronl (y-or-n-p "Is this reference vector read-only? "))
  (setq rebi (y-or-n-p
	      "Should the references in the vector be automatically rebound? "))
  (setq null (y-or-n-p "Can the references in the vector be null? "))
  (setq defn
	(y-or-n-p
	 "Should the reference in the vector be given a default if null? "))

  (insert-string (concat "
  static RefVector<" class "," type "> interface" name "
    (\"" name "\",
     \"" desc "\",
     " memb ", " size ", "
     (cond (safe "true")(t "false")) ", " (cond (ronl "true")(t "false")) ", "
     (cond (rebi "true")(t "false")) ", " (cond (null "true")(t "false")) ", "
     (cond (defn "true")(t "false"))))
  (cond ((y-or-n-p "Are there any set/get functions? ")
	 (insert-string (concat ",
     " (thepeg-ifn "Set-function: " (concat "&" class "::set" name)
		   (concat "(void(" class "::*)(Ptr<"
			   type ">::pointer,int))(0)"))
     ",
     "
     (thepeg-ifn "Insert-function: " (concat "&" class "::ins" name)
		 (concat "(void(" class "::*)(Ptr<"
			 type ">::pointer,int))(0)"))
     ",
     "
     (thepeg-ifn "Remove-function: " (concat "&" class "::del" name)
		 (concat "(void(" class "::*)(int))(0)"))
     ",
     "
     (thepeg-ifn "Get-function: " (concat "&" class "::get" name)
		 (concat "(vector<Ptr<"
			 type ">::pointer>(" class "::*)()const)(0)"))
     ",
     "
     (thepeg-ifn "Check-function: "
		 (concat "&" class "::check" name)
		 (concat "(bool(" class "::*)(Ptr<" type
			 ">::const_pointer,int)const)(0)"))))))
  (insert-string ");
"))

(defun thepeg-switch-option (switchname)
  "Create a SwitchOption variable suitable for inclusion in an Init() function."
  (cond ((y-or-n-p "Do you want to add an option? ")
	 (setq name (read-from-minibuffer "Name: "))
	 (setq desc (read-from-minibuffer "Description: "))
	 (setq valu (read-from-minibuffer "Value: "))
	 (insert-string (concat "
  static SwitchOption " switchname name "
    (" switchname ",
     \"" name "\",
     \"" desc "\",
     " valu ");"))
	 (thepeg-switch-option switchname))
	(t (insert-string "
"))))

(defun ThePEG-switch ()
  "Create a Switch variable suitable for inclusion in an Init() function."
  (interactive)
  (setq class (read-from-minibuffer "Class Name: "
				    (file-name-sans-extension
				     (file-name-nondirectory
				      (buffer-file-name)))))
  (setq type (read-from-minibuffer "Integer type: " "int"))
  (setq name (read-from-minibuffer "Name: "))
  (setq desc (read-from-minibuffer "Description: "))
  (setq memb (read-from-minibuffer "Member variable: "
				   (concat "&" class "::the" name)))
  (setq defa (read-from-minibuffer "Default value: " "0"))
  (setq safe (y-or-n-p "Is this reference dependency safe? "))
  (setq ronl (y-or-n-p "Is this reference read-only? "))
  (insert-string (concat "
  static Switch<" class "," type "> interface" name "
    (\"" name "\",
     \"" desc "\",
     " memb ", " defa ", "
     (cond (safe "true")(t "false")) ", " (cond (ronl "true")(t "false"))))
  (cond ((y-or-n-p "Are there any set/get functions? ")
	 (insert-string (concat ",
     " (thepeg-ifn "Set-function: " (concat "&" class "::set" name)
		   (concat "(void(" class "::*)(" type "))(0)"))
     ",
     "
     (thepeg-ifn "Get-function: " (concat "&" class "::get" name)
		   (concat "(" type "(" class "::*)()const)(0)"))
     ", "
     (thepeg-ifn "Default-function: "
			   (concat "&" class "::def" name)
		   (concat "(" type "(" class "::*)()const)(0)"))))))
  (insert-string ");")
  (thepeg-switch-option (concat "interface" name)))

(defun ThePEG-comment-group ()
  "Group a number of functions (or anything within the currently
   marked region) together to a common Doxygen section."
  (interactive)
  (setq desc (read-from-minibuffer "Description: "))
  (save-excursion
    (goto-char (region-beginning))
    (beginning-of-line)
    (cond ((string-equal desc ""))
	  (t (insert-string (concat "  /** @name " desc " */
"))))
    (insert-string "  //@{
"))
  (save-excursion
    (goto-char (region-end))
    (beginning-of-line)
    (insert-string "  //@}
")))

(defun ThePEG-cond-section ()
  "Make a region of code a conditional section for doxygen."
  (interactive)
  (setq sect (read-from-minibuffer "Section: "))
  (cond ((string-equal sect "") (setq sect "DEVELOPER")))
  (save-excursion
    (goto-char (region-beginning))
    (beginning-of-line)
    (insert-string (concat "  /** @cond " sect " */
")))
  (save-excursion
    (goto-char (region-end))
    (insert-string "  /** @endcond */
")))

(defun ThePEG-variable-cutnpaste ()
  "Create code fragments for a variable declaration together with
   corresponding set and get functions, I/O functions and constructors."
  (interactive)
  (setq class (read-from-minibuffer "Class Name: "
				    (file-name-sans-extension
				     (file-name-nondirectory
				      (buffer-file-name)))))
  (setq type (read-from-minibuffer "Type: " "double"))
  (setq name (read-from-minibuffer "Name (of access function): "))
  (setq theName (concat (cond ((string-equal type "bool") "use")
			      (t "the")) (upcase-initials name)))
  (setq theName (read-from-minibuffer "Name (of variable): " theName))
  (setq desc (read-from-minibuffer "Description: "))
  (setq unit (read-from-minibuffer "Unit: "))
  (setq outname (cond ((string-equal unit "") theName)
		      (t (concat "ounit(" theName ", " unit ")"))))
  (setq inname (cond ((string-equal unit "") theName)
		     (t (concat "iunit(" theName ", " unit ")"))))
  (setq def (cond ((string-equal unit "") "0")
		  (t (concat "0*" unit))))
  (setq buf (get-buffer "*ThePEG-scratch*"))
  (setq win (cond (buf (get-buffer-window buf 'visible))
		  (t nil)))
  (setq getdesc (replace-regexp-in-string "^The" "Get the" desc t))
  (setq setdesc (replace-regexp-in-string "^The" "Set the" desc t))
  (cond (win (select-window win))
	(t (switch-to-buffer-other-window "*ThePEG-scratch*")))
  (c++-mode)
  (beginning-of-buffer)
  (insert-string (concat "  /**
   * " desc "
   */
  " type " " theName ";

  /**
   * " getdesc "
   */
  inline " type " " name "() const {
    return " theName ";
  }

  /**
   * " setdesc "
   */
  inline void " name "(" type " x) {
    " theName " = x;
  }

" class "::" class "() : Base(), " theName "(" def ") {}

" class "::" class "(const " class " & x) : Base(x), " theName "(x." theName ") {}

void " class "::persistentOutput(PersistentOStream & os) const {
  os << " outname ";
}

void " class "::persistentInput(PersistentIStream & is, int) {
  is >> " inname ";
}

"
))
(beginning-of-buffer))

(defun ThePEG-comment-function ()
  "Insert a Doxygen comment block to preceed a function declaration."
  (interactive)
  (insert-string "
  /**
   * 
   */
  "))

(global-set-key "\C-c\C-f" 'ThePEG-comment-function)

(define-generic-mode 'ThePEG-repository-mode
;; Comment characters
      '("#")
;; Language Keywords 
;; codegrep 'verb == \"[[:alnum:]]*\"' -o \
;; | cut -d: -f3 | cut -d\  -f3 | sort | uniq | column
      '(
	"appendpath"		"library"		"popd"
	"baseclasses"		"load"			"prependpath"
	"cd"			"ls"			"pushd"
	"check"			"lsclass"		"pwd"
	"cp"			"lspaths"		"rcp"
	"create"		"lsruns"		"read"
	"decaymode"		"makeanti"		"removerun"
	"def"			"makerun"		"rm"
	"defaultparticle"	"max"			"rmdir"
	"describe"		"mdef"			"rmgloballibrary"
	"describeclass"		"mdo"			"rmrun"
	"dirs"			"merase"		"rrmdir"
	"DISABLEREADONLY"	"mget"			"run"
	"do"			"min"			"save"
	"doxygendump"		"minsert"		"saverun"
	"ENABLEREADONLY"	"mkdir"			"saverunfile"
	"erase"			"mmax"			"send"
	"EXITONERROR"		"mmin"			"set"
	"fulldescribe"		"msend"			"setdef"
	"get"			"mset"			"setup"
	"globallibrary"		"msetdef"               "newdef"
	"insert"		"mv"                    "clear"
	)
;; Further Regex highlighting
      '(
	;; anchored on "read" highlight the second pattern
	("read" "[[:alnum:]]+\\.in" nil nil (0 'font-lock-string-face) )
	;; simplest pattern highlight
	;;("[[:alnum:]]+\\.so" . 'font-lock-function-name-face)
	("[[:alnum:]]+\\.so" . 'font-lock-type-face)
	("<[[:alnum:]]+>" . 'font-lock-function-name-face)
	("([[:alnum:]]+)" . 'font-lock-function-name-face)
	("\\([[:alnum:]]+::\\)+[[:alnum:]]+" . 'font-lock-type-face)


	(":\\([]_/[:alnum:]\\[]+\\)" . (1 'font-lock-constant-face) )


	;; highlight everything that isn't painted yet: 'keep'
	("create.*"  . (0 'font-lock-variable-name-face keep)	)

;;;;;;;;;; Directories in the repo
;;	("/[/[:alnum:]]+/" . 'font-lock-builtin-face )
;;	("\\(rrm\\|rm\\|mk\\)dir.*\\|cd.*"  
;;	 . (0 'font-lock-builtin-face keep)	)
;;;;;;;;;;;;;;;;;;;;;
	)
      ;; file pattern to apply this mode to
      nil ;;'(".in\\'")
      nil
      "Major mode for editing ThePEG repository input files.")
