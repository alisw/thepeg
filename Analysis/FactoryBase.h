// -*- C++ -*-
//
// FactoryBase.h is a part of ThePEG - Toolkit for HEP Event Generation
// Copyright (C) 1999-2017 Leif Lonnblad
//
// ThePEG is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
#ifndef THEPEG_FactoryBase_H
#define THEPEG_FactoryBase_H
//
// This is the declaration of the FactoryBase class.
//

#include "ThePEG/Interface/Interfaced.h"
#include "FactoryBase.fh"

namespace AIDA {
  class IHistogram1D;
  class IHistogram2D;
  class IDataPointSet;
  class IAnalysisFactory;
  class ITree;
  class IHistogramFactory;
  class IDataPointSetFactory;
}

namespace ThePEG {

/**
 * Here is the documentation of the FactoryBase class. This
 * abstract class is used to wrap the interface to a particular
 * AIDA-compliant histogram package to be used in
 * <code>AnalysisHandler</code>s in ThePEG. Concrete subclasses must
 * implement the doinitrun() function to create an object of a class
 * inheriting from AIDA::IAnalysisFactory and assign it with the
 * analysisFactory(AIDA::IAnalysisFactory*) function before calling
 * doinitrun() for the FactoryBase base class.
 *
 * A FactoryBase object should be assigned to the EventGenerator
 * object controlling a run, and <code>AnalysisHandler</code>s should
 * access it via the Generator with the generator() function.
 *
 * @see \ref FactoryBaseInterfaces "The interfaces"
 * defined for FactoryBase.
 */
class FactoryBase: public Interfaced {

public:

  /**
   * Convenient typedef for pointer to AIDA::IHistogram1D.
   */
  typedef AIDA::IHistogram1D * tH1DPtr;

  /**
   * Convenient typedef for pointer to const AIDA::IHistogram1D.
   */
  typedef const AIDA::IHistogram1D * tcH1DPtr;

  /**
   * Convenient typedef for pointer to AIDA::IHistogram2D.
   */
  typedef AIDA::IHistogram2D * tH2DPtr;

  /**
   * Convenient typedef for pointer to const AIDA::IHistogram2D.
   */
  typedef const AIDA::IHistogram2D * tcH2DPtr;

  /**
   * Convenient typedef for pointer to AIDA::IHistogram1D.
   */
  typedef AIDA::IDataPointSet * tDSetPtr;

  /**
   * Convenient typedef for pointer to const AIDA::IHistogram1D.
   */
  typedef const AIDA::IDataPointSet * tcDSetPtr;

public:

  /**
   * DataFiller is a helper class to facilitate adding data to a
   * DataPointSet. For a D-dimensional DataPointSet N*3*D numbers should
   * be added with the standard &lt;&lt; operator ordered as x-value,
   * x-upper-error, x-lower-error, y-value, y-upper-error, etc.. Only
   * when the DataFIller object is detleted will the points be added.
   */
  class DataFiller {

  public:

    /**
     * The standard constructor needs a IDataPointSet as argument.
     */
    DataFiller(AIDA::IDataPointSet * dps) : dset(dps) {}

    /**
     * Copy constructor.
     */
    DataFiller(const DataFiller & df) : dset(df.dset) {}

    /**
     * Destructor. Will commit the pints filled to the underlying
     * IDataPointSet.
     */
    ~DataFiller();

    /**
     * Add a number to measurement currently being read.
     */
    DataFiller & operator<<(double x) {
      v.push_back(x);
      return *this;
    }

    /**
     * Automatic conversion to the underlying IDataPointSet.
     */
    operator AIDA::IDataPointSet * () { return dset; }

  private:

    /**
     * The underlying IDataPointSet.
     */
    AIDA::IDataPointSet * dset;

    /**
     * The collected numbers to be committed to the IDataPointSet.
     */
    deque<double> v;

  };

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  FactoryBase();

  /**
   * The copy constructor.
   */
  FactoryBase(const FactoryBase &);

  /**
   * The destructor.
   */
  virtual ~FactoryBase();
  //@}

public:

  /** @name Simple access functions. */
  //@{
  /**
   * Together with suffix(), the name of the file where the resulting
   * histograms will be stored. If empty, generator()->filename() will
   * be used instead.
   */
  const string & filename() const { return theFilename; }

  /**
   * Together with filename(), the name of the file where the
   * resulting histograms will be stored.
   */
  const string & suffix() const { return theSuffix; }

  /**
   * The format in which the histograms are stored in the output file.
   */
  const string & storeType() const { return theStoreType; }
  //@}

  /** @name Manipulate histograms */
  //@{
  /**
   * Rescale the given \a histogram so that the integral over the bins
   * will give the correct integrated cross section for the observable
   * in the given \a unit.
   */
  virtual void
  normalizeToXSec(tH1DPtr histogram, CrossSection unit = picobarn) const = 0;

  /**
   * Rescale the given \a histogram so that the integral over the bins
   * will give the correct integrated cross section for the observable
   * in the given \a unit.
   */
  virtual void
  normalizeToXSec(tH2DPtr histogram, CrossSection unit = picobarn) const = 0;

  /**
   * Rescale the given \a histogram so that the integral over the bins
   * gives the fraction of the total cross section generated which is
   * contained in the bins.
   */
  virtual void normalizeToXSecFraction(tH1DPtr histogram) const = 0;

  /**
   * Rescale the given \a histogram so that the integral over the bins
   * gives the fraction of the total cross section generated which is
   * contained in the bins.
   */
  virtual void normalizeToXSecFraction(tH2DPtr histogram) const = 0;

  /**
   * Rescale the given \a histogram so that the integral over the bins
   * gives one.
   */
  virtual void normalizeToUnity(tH1DPtr histogram) const = 0;

  /**
   * Rescale the given \a histogram so that the integral over the bins
   * gives one.
   */
  virtual void normalizeToUnity(tH2DPtr histogram) const = 0;
  //@}

  /** @name Access the underlying AIDA objects. */
  //@{
  /**
   * Access the underlying AIDA::IAnalysisFactory object.
   */
  AIDA::IAnalysisFactory & analysisFactory() const {
    return *theAnalysisFactory;
  }

  /**
   * Access the underlying AIDA::ITree object.
   */
  AIDA::ITree & tree() const;

  /**
   * A pointer to the underlying AIDA::IHistogramFactory object.
   */
  AIDA::IHistogramFactory & histogramFactory() const;

  /**
   * A pointer to the underlying AIDA::IDataPointSetFactory object.
   */
  AIDA::IDataPointSetFactory & dataSetFactory() const;

  /**
   * Create a new directory in the underlying AIDA tree.
   */
  void mkdir(const string & path);

  /**
   * Create a new directory in the underlying AIDA tree.
   */
  void mkdirs(const string & path);

  /**
   * Set the default working directory for the underlying AIDA tree.
   */
  void cd(const string & path);

  /**
   * Create and return a AIDA::IHistogram1D object in the underlying
   * AIDA histogram factory. Note that the histogram factory is
   * responsible for deleting this histogram.
   * @param path the full path of where the histogram should be placed
   * in the underlying AIDA tree (on the form
   * "/dir/subdir/histogramname"). Not that the directory part of the
   * path typically must already exist in the tree. The directories
   * can be created with mkdir(string) or mkdirs(string).
   * The title of the histogram will be set to the name part of the path.
   * @param nb the number of bins in the histogram.
   * @param lo the lower edge of the histogram.
   * @param up the upper edge of the histogram.
   * @return a pointer to the created AIDA::IHistogram1D object.
   */
  tH1DPtr createHistogram1D(const string & path, int nb, double lo, double up);

  /**
   * Create and return a AIDA::IHistogram1D object in the underlying
   * AIDA histogram factory. Note that the histogram factory is
   * responsible for deleting this histogram.
   * @param path the full path of where the histogram should be placed
   * in the underlying AIDA tree (on the form
   * "/dir/subdir/histogramname"). Not that the directory part of the
   * path typically must already exist in the tree. The directories
   * can be created with mkdir(string) or mkdirs(string).
   * @param title the title of the histogram.
   * @param nb the number of bins in the histogram.
   * @param lo the lower edge of the histogram.
   * @param up the upper edge of the histogram.
   * @return a pointer to the created AIDA::IHistogram1D object.
   */
  tH1DPtr createHistogram1D(const string & path, const string & title, int nb,
			    double lo, double up);

  /**
   * Create and return a AIDA::IHistogram1D object in the underlying
   * AIDA histogram factory. Note that the histogram factory is
   * responsible for deleting this histogram.
   * @param path the full path of where the histogram should be placed
   * in the underlying AIDA tree (on the form
   * "/dir/subdir/histogramname"). Not that the directory part of the
   * path typically must already exist in the tree. The directories
   * can be created with mkdir(string) or mkdirs(string).
   * @param title the title of the histogram.
   * @param edges A vector of bin edges specifying th bins.
   * @return a pointer to the created AIDA::IHistogram1D object.
   */
  tH1DPtr createHistogram1D(const string & path, const string & title,
			    const std::vector<double> & edges);

  /**
   * Create and return a AIDA::IHistogram2D object in the underlying
   * AIDA histogram factory. Note that the histogram factory is
   * responsible for deleting this histogram.
   * @param path the full path of where the histogram should be placed
   * in the underlying AIDA tree (on the form
   * "/dir/subdir/histogramname"). Not that the directory part of the
   * path typically must already exist in the tree. The directories
   * can be created with mkdir(string) or mkdirs(string).
   * The title of the histogram will be set to the name part of the path.
   * @param nbx the number of x-bins in the histogram.
   * @param xlo the lower x-edge of the histogram.
   * @param xup the upper x-edge of the histogram.
   * @param nbx the number of y-bins in the histogram.
   * @param xlo the lower y-edge of the histogram.
   * @param xup the upper y-edge of the histogram.
   * @return a pointer to the created AIDA::IHistogram1D object.
   */
  tH2DPtr createHistogram2D(const string & path,
			    int nbx, double xlo, double xup,
			    int nby, double ylo, double yup);

  /**
   * Create and return a AIDA::IHistogram2D object in the underlying
   * AIDA histogram factory. Note that the histogram factory is
   * responsible for deleting this histogram.
   * @param path the full path of where the histogram should be placed
   * in the underlying AIDA tree (on the form
   * "/dir/subdir/histogramname"). Not that the directory part of the
   * path typically must already exist in the tree. The directories
   * can be created with mkdir(string) or mkdirs(string).
   * @param title the title of the histogram.
   * @param nbx the number of x-bins in the histogram.
   * @param xlo the lower x-edge of the histogram.
   * @param xup the upper x-edge of the histogram.
   * @param nby the number of y-bins in the histogram.
   * @param ylo the lower y-edge of the histogram.
   * @param yup the upper y-edge of the histogram.
   * @return a pointer to the created AIDA::IHistogram1D object.
   */
  tH2DPtr createHistogram2D(const string & path, const string & title,
			    int nbx, double xlo, double xup,
			    int nby, double ylo, double yup);

  /**
   * Create and return a AIDA::IHistogram2D object in the underlying
   * AIDA histogram factory. Note that the histogram factory is
   * responsible for deleting this histogram.
   * @param path the full path of where the histogram should be placed
   * in the underlying AIDA tree (on the form
   * "/dir/subdir/histogramname"). Not that the directory part of the
   * path typically must already exist in the tree. The directories
   * can be created with mkdir(string) or mkdirs(string).
   * @param title the title of the histogram.
   * @param xedges A vector of bin edges specifying the x-bins.
   * @param yedges A vector of bin edges specifying the y-bins.
   * @return a pointer to the created AIDA::IHistogram1D object.
   */
  tH2DPtr createHistogram2D(const string & path, const string & title,
			    const std::vector<double> & xedges,
			    const std::vector<double> & yedges);

  /**
   * Create a IDataPointSet with the given \a path and \a title and
   * containing points with the given number of dimensions, \a
   * dim. The returned object is a DataFiller which can be used to
   * facilitate the addition of data points to the set or be converted
   * to a pointer to the created IDataPointSet.
   */
  DataFiller createDataSet(const string & path, const string & title, int dim);

  /**
   * Used by a \a client object to indicate that he has required
   * histograms from this factory. It is guaranteed that the clients
   * finish() function is called before the underlying AIDA::ITree is
   * committed and the AIDA::IHistogramFactory is deleted together
   * with all histograms.
   */
  void registerClient(tIPtr client);
  //@}

protected:

  /**
   * Set the underlying AIDA::IAnalysisFactory object. Note that this
   * surrenders the controll of the factory to the FactoryBase
   * object which will delete it in the finish() function. Typically
   * this function should be called by a concrete subclass in the
   * doinitrun() function before the doinitrun() function of this
   * class is called.
   */
  void analysisFactory(AIDA::IAnalysisFactory * x) {
    theAnalysisFactory = x;
  }

  /**
   * Delete all associated AIDA objects. Note that the tree is not
   * explicitly committed.
   */
  void clear();


public:

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

  /**
   * The standard Init function used to initialize the interfaces.
   * Called exactly once for each class by the class description system
   * before the main function starts or
   * when this class is dynamically loaded.
   */
  static void Init();

protected:

  /** @name Standard Interfaced functions. */
  //@{
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
  //@}

private:

  /**
   * Together with theSuffix, the name of the file where the resulting
   * histograms will be stored. If empty, generator()->filename() will
   * be used instead.
   */
  string theFilename;

  /**
   * Together with theFilename, the name of the file where the
   * resulting histograms will be stored.
   */
  string theSuffix;

  /**
   * The format in which the histograms are stored in the output file.
   */
  string theStoreType;

  /**
   * A pointer to the underlying AIDA::IAnalysisFactory object.
   */
  AIDA::IAnalysisFactory * theAnalysisFactory;

  /**
   * A pointer to the underlying AIDA::ITree object.
   */
  AIDA::ITree * theTree;

  /**
   * A pointer to the underlying AIDA::IHistogramFactory object.
   */
  AIDA::IHistogramFactory * theHistogramFactory;

  /**
   * A pointer to the underlying AIDA::IDataPointSetFactory object.
   */
  AIDA::IDataPointSetFactory * theDataSetFactory;

  /**
   * A set of client objects which have required histograms from this
   * factory.
   */
  set<IPtr> clients;

private:

  /**
   * The static object used to initialize the description of this class.
   * Indicates that this is an abstract class with persistent data.
   */
  static AbstractClassDescription<FactoryBase> initFactoryBase;

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  FactoryBase & operator=(const FactoryBase &);

};

}

#include "ThePEG/Utilities/ClassTraits.h"

namespace ThePEG {

/** @cond TRAITSPECIALIZATIONS */

/** This template specialization informs ThePEG about the
 *  base classes of FactoryBase. */
template <>
struct BaseClassTrait<FactoryBase,1> {
  /** Typedef of the first base class of FactoryBase. */
  typedef Interfaced NthBase;
};

/** This template specialization informs ThePEG about the name of
 *  the FactoryBase class and the shared object where it is defined. */
template <>
struct ClassTraits<FactoryBase>
  : public ClassTraitsBase<FactoryBase> {
  /** Return a platform-independent class name */
  static string className() { return "ThePEG::FactoryBase"; }
};

/** @endcond */

}

#endif /* THEPEG_FactoryBase_H */
