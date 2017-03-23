package ThePEG;

import java.util.*;

public class ObjectRef {

  SetupThePEG thepeg;
  String name;
  String fullName;
  boolean dir = false;
  String pclass;

  public ObjectRef(SetupThePEG tp, String nm, String fn, boolean isd) {
    thepeg = tp;
    name = nm;
    fullName = fn;
    dir = isd;
  }

  public ObjectRef(SetupThePEG tp, String f) {
    thepeg = tp;
    fullName = f;
    if ( f.charAt(f.length() - 1) == '/' ) {
      f = f.substring(0, f.length() - 1);
      dir = true;
    }
    name = f.substring(f.lastIndexOf('/') + 1);
  }

  public String toString() {
    return name;
  }

  public String getClassName() {
    if ( pclass != null ) return pclass;
    if ( isDir() ) return "";
    LinkedList ret = thepeg.exec("fulldescribe " + getFullName());
    if ( ret.size() > 0 ) ret.remove(0);
    if ( ret.size() > 0 ) pclass = (String)ret.remove(0);  
    return pclass;
  }

  public String getName() {
    return name;
  }

  public String getFullName() {
    return fullName;
  }

  public boolean isDir() {
    return dir;
  }

  public void open() {
    if ( isDir() ) return;
    thepeg.openObject(getFullName());
  }

  public static void classcheck() {}

}

