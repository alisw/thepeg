package ThePEG;

import javax.swing.*;
import java.util.*;
import java.awt.*;
import java.awt.event.*;

public class Reference extends Interface implements ActionListener {

  boolean nullable;
  boolean defnull;
  String current;
  RefRef selected = new RefRef();
  String dir;
  String pclass;

  JButton ok = new JButton("Ok");
  JButton cancel = new JButton("Cancel");
  JButton apply = new JButton("Apply");
  JButton reset = new JButton("Reset");
  JButton setnull = new JButton("Null");
  JButton open = new JButton("Open");
  JButton sel = new JButton("Select");

  public Reference(SetupThePEG own, ObjectFrame obj, LinkedList input) {
    super(own, obj, input);
    if ( !setup(input) ) {
      JOptionPane.showMessageDialog(own, "Could not create Reference view",
				    "Error", JOptionPane.ERROR_MESSAGE);
      return;
    }

    getContentPane().setLayout(new BorderLayout());
    getContentPane().add(getDescriptionArea(), BorderLayout.CENTER);
    JPanel buttons = new JPanel();
    buttons.add(cancel);
    if ( !isReadonly() ) buttons.add(reset);
    if ( !isReadonly() && nullable ) buttons.add(setnull);
    if ( !isReadonly() ) buttons.add(sel);
    buttons.add(open);
    if ( !isReadonly() ) buttons.add(apply);
    buttons.add(ok);
    cancel.addActionListener(this);
    if ( !isReadonly() ) reset.addActionListener(this);
    if ( !isReadonly() ) sel.addActionListener(this);
    if ( !isReadonly() && nullable ) setnull.addActionListener(this);
    open.addActionListener(this);
    if ( !isReadonly() ) apply.addActionListener(this);
    ok.addActionListener(this);
    fixButtons();
    getContentPane().add(buttons, BorderLayout.SOUTH);
    Font f = selected.getFont();
    selected.setFont(new Font(f.getName(), f.getStyle(), f.getSize() + 2));
    selected.addMouseListener(new MouseAdapter() {
	public void mousePressed(MouseEvent e) {
	  if ( e.getSource() == selected && e.getClickCount() >= 2 &&
	       ! selected.getName().equals("NULL") ) owner.openObject(current);
	}
      });
    //    getContentPane().add(selected, BorderLayout.NORTH);
    JPanel top = new JPanel();
    top.add(new JLabel("Current Value: "));
    top.add(selected);
    getContentPane().add(top, BorderLayout.NORTH);
    setTitle("Reference");
    setupFrame(500,150);
  }

  protected boolean setup(LinkedList input) {
    String s = (String)input.getFirst();
    pclass = s.substring(2, s.length() - 1);
    if ( !super.setup(input) ) return false;
    s = (String)input.remove(0);
    if ( s.equals("nullable") ) nullable = true;
    else if ( s.equals("nevernull") ) nullable = false;
    else return false;
    s = (String)input.remove(0);
    if ( s.equals("defnull") ) defnull = true;
    else if ( s.equals("nodefnull") ) defnull = false;
    else return false;
    current = (String)input.remove(0);
    selected.setup(current);
    dir = "";
    if ( current.lastIndexOf("/") >= 0 )
      dir = current.substring(0, current.lastIndexOf("/") + 1);
    if ( dir.equals("") ) dir = object.substring(0,object.lastIndexOf("/") + 1);
    if ( dir.equals("") ) dir = "/";
    return true;
  }

  private void setValue() {
    setValue(selected.getFullName());
  }

  public void actionPerformed(ActionEvent e) {
    if ( e.getSource() == cancel ) {
      dispose();
    }
    else if ( e.getSource() == ok ) {
      if ( !isReadonly() ) setValue();
      dispose();
    }
    else if ( e.getSource() == apply ) {
      setValue();
    }
    else if ( e.getSource() == setnull ) {
      selected.setup("NULL", "NULL");
    }
    else if ( e.getSource() == sel ) {
      ObjectSelector ch = new ObjectSelector(owner, dir, pclass, this);
      if ( ch.selected() != null ) selected.setup(ch.selected().getFullName());
    }
    else if ( e.getSource() == open ) {
      if ( ! selected.getName().equals("NULL") ) owner.openObject(current);
    }
    else if ( e.getSource() == reset ) {
      selected.setup(current);
    }
    fixButtons();
  }

  public void fixButtons() {
    setnull.setEnabled(!selected.getName().equals("NULL"));
    open.setEnabled(!selected.getName().equals("NULL"));
    apply.setEnabled(!selected.getFullName().equals(current));
    reset.setEnabled(!selected.getFullName().equals(current));
  }

  public static void classcheck() {}

}


