package ThePEG;

import javax.swing.*;
import java.util.*;
import java.awt.*;
import java.awt.event.*;

public class SwitchOption extends JLabel {

  String name;
  String comment;
  long index;

  public SwitchOption(long indx, String nm, String cm) {
    index = indx;
    name = nm;
    comment = cm;
    setText(Long.toString(index) + ": " + name);
    setToolTipText(SetupThePEG.ttt + comment);
    setBorder(BorderFactory.createEmptyBorder(2,2,2,2));
  }

  public long getIndex() {
    return index;
  }

  public String getName() {
    return name;
  }

  public String getComment() {
    return comment;
  }

  public String toString() {
    return getText();
  }

  public static void classcheck() {}

}
