package ThePEG;

import javax.swing.*;
import java.util.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.text.html.HTMLDocument;

public class Interface extends JDialog {

  SetupThePEG owner;
  String object;
  String objectName;
  String interfaceName;
  String description = "";
  boolean readonly;

  //  JTextArea info = new JTextArea(new HTMLDocument());
  JEditorPane info = new JEditorPane();

  public Interface(SetupThePEG own, ObjectFrame obj, LinkedList input) {
    super(obj);
    owner = own;
    object = obj.getFullName();
    objectName = obj.getName();
    info.addHyperlinkListener(obj);
  }

  protected boolean setup(LinkedList input) {
    String s = (String)input.remove(0);
    if ( input.size() > 0 ) interfaceName = (String)input.remove(0);
    else return false;
    while ( input.size() > 0 && ( s = (String)input.remove(0) ) != null &&
	    !s.equals("-*-mutable-*-") && !s.equals("-*-readonly-*-") )
      description += s + " ";
    if ( s.equals("-*-mutable-*-") ) readonly = false;
    else if ( s.equals("-*-readonly-*-") ) readonly = true;
    else return false;
    return true;
  }

  protected boolean reset() {
    return setup(owner.exec("fulldescribe " + object + ":" + interfaceName));
  }

  protected void setupFrame(int width, int height) {
    setDefaultCloseOperation(DISPOSE_ON_CLOSE);
    owner.setLocation(this);
    setSize(width, height);
    setVisible(true);
  }

  public void setTitle(String type) {
    super.setTitle("ThePEG " + type + " " + getObjectName() + ":" + getName());
  }
    
  protected JScrollPane getDescriptionArea() {
    info.setContentType("text/html");
    info.setText("<font size=-1>" + ObjectFrame.htmlFormat(description) +
		 "</font>");
    info.setEditable(false);
    info.setCaretPosition(0);
    return new JScrollPane(info);
  }

  protected boolean isReadonly() {
    return readonly;
  }

  protected void exec(String cmd) {
    owner.action(cmd);
  }

  protected void setValue(String val) {
    exec("set " + getObject() + ":" + getName() + " " + val);
    reset();
  }

  protected void setValue(int indx, String val) {
    exec("set " + getObject() + ":" + getName() + " [" +
	 Integer.toString(indx) + "] " + val);
    reset();
  }

  protected void insert(int indx, String val) {
    exec("insert " + getObject() + ":" + getName() + "[" +
	 Long.toString(indx) + "] " + val);
    reset();    
  }

  protected void erase(int indx) {
    exec("erase " + getObject() + ":" + getName() + "[" +
	 Long.toString(indx) + "] ");
    reset();    
  }

  public String getName() {
    return interfaceName;
  }

  public String getObject() {
    return object;
  }

  public String getObjectName() {
    return objectName;
  }

  public void dispose() {
    owner.removeLocation(this);
    super.dispose();
  }

  public static void classcheck() {}

}
