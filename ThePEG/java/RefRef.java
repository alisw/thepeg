package ThePEG;

import javax.swing.*;
import java.util.*;
import java.awt.*;
import java.awt.event.*;

public class RefRef extends JLabel {

  String name;
  String fullname;
  int index;

  public RefRef(int indx, String nm, String fnm) {
    setup(indx, nm, fnm);
  }

  public RefRef(int indx, String fnm) {
    if ( fnm.equals("NULL") ) setup(indx, fnm, fnm);
    else setup(indx, fnm.substring(fnm.lastIndexOf('/')+1), fnm);
  }

  public RefRef(String nm, String fnm) {
    setup(-1, nm, fnm);
  }

  public RefRef(String fnm) {
    setup(-1, fnm);
  }

  public RefRef() {
    setup(-1, "");
  }

  public void setup(int indx, String nm, String fnm) {
    index = indx;
    name = nm;
    fullname = fnm;
    String text = name;
    if ( index >= 0 && !name.equals(" - end - ") )
      text = Long.toString(index) + ": " + name;
    setText(text);
    setToolTipText(fullname);
  }

  public String toString() {
    return getText();
  }    

  public void setup(int indx, String fnm) {
    setup(indx, fnm.substring(fnm.lastIndexOf('/')+1), fnm);
  }

  public void setup(String nm, String fnm) {
    setup(-1, nm, fnm);
  }

  public void setup(String fnm) {
    setup(-1, fnm);
  }

  public int getIndex() {
    return index;
  }

  public String getName() {
    return name;
  }

  public String getFullName() {
    return fullname;
  }

  public void set(String fnm) {
    setup(index, fnm);
  }

  public void set(String nm, String fnm) {
    setup(index, nm, fnm);
  }

  public static void classcheck() {}

}
