// -*- C++ -*-
#ifndef ThePEG_Settings_H
#define ThePEG_Settings_H
//
// This is the declaration of the Settings class.
//

#include "ThePEG/Handlers/HandlerBase.h"

namespace ThePEG {

/**
 * Here is the documentation of the Settings class.
 *
 * @see \ref SettingsInterfaces "The interfaces"
 * defined for Settings.
 */
class Settings: public HandlerBase {

public:

  /** @name Standard constructors and destructors. */
  //@{
  /**
   * The default constructor.
   */
  Settings();

  /**
   * The destructor.
   */
  virtual ~Settings();
  //@}

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

public:

  /**
   * The function used for the set command.
   */
  string setFunction(string);

protected:

  /** @name Standard Interfaced functions. */
  //@{
  /**
   * This object needs to be initialized before all
   * other objects.
   */
  virtual bool preInitialize() const;

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

private:

  /**
   * The objects that needs setting.
   */
  vector<string> setObjects;

  /**
   * The interfaces that needs setting.
   */
  vector<string> setInterfaces;

  /**
   * The values to be set for the interfaces of the objects.
   */
  vector<string> setValues;

private:

  /**
   * The assignment operator is private and must never be called.
   * In fact, it should not even be implemented.
   */
  Settings & operator=(const Settings &);

};

}

#endif /* ThePEG_Settings_H */
