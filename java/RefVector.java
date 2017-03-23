package ThePEG;

import javax.swing.*;
import javax.swing.event.*;
import java.util.*;
import java.awt.*;
import java.awt.event.*;

public class RefVector extends Interface
  implements ActionListener, ListSelectionListener {

  boolean nullable;
  boolean defnull;
  Vector current;
  Vector values;
  String pclass;
  String dir;
  int size;

  JButton ok = new JButton("Ok");
  JButton cancel = new JButton("Cancel");
  JButton apply = new JButton("Apply");
  JButton reset = new JButton("Reset");
  JButton setnull = new JButton("Null");
  JButton open = new JButton("Open");
  JButton sel = new JButton("Select");
  JButton insert = new JButton("Insert");
  JButton erase = new JButton("Erase");
  JList selector = new JList() {
      public String getToolTipText(MouseEvent e) {
	int i = locationToIndex(e.getPoint());
	if ( i < 0 || i > current.size() - 2 ) return null;
	RefRef r = (RefRef)current.get(i);
	if ( r == null ) return null;
	return r.getFullName();
      }
    };

  public RefVector(SetupThePEG own, ObjectFrame obj, LinkedList input) {
    super(own, obj, input);
    if ( !setup(input) ) {
      JOptionPane.showMessageDialog(own,
				    "Could not create Reference Vector view",
				    "Error", JOptionPane.ERROR_MESSAGE);
      return;
    }

    getContentPane().setLayout(new BorderLayout());
    getContentPane().add(getDescriptionArea(), BorderLayout.CENTER);
    selector.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
    selector.addListSelectionListener(this);
    selector.addMouseListener(new MouseAdapter() {
	public void mousePressed(MouseEvent e) {
	  int i = selector.locationToIndex(e.getPoint());
	  if ( i < 0 || i > current.size() - 2 ) return;
	  RefRef r = (RefRef)current.get(i);
	  if ( r == null ) return;
	  if ( e.getClickCount() >= 2) owner.openObject(r.getFullName());
	}
      });
    getContentPane().add(new JScrollPane(selector), BorderLayout.NORTH);
    JPanel buttons = new JPanel();
    buttons.add(cancel);
    if ( !isReadonly() ) buttons.add(reset);
    if ( !isReadonly() && nullable ) buttons.add(setnull);
    if ( !isReadonly() ) buttons.add(sel);
    if ( !isReadonly() && size <= 0 ) buttons.add(insert);
    if ( !isReadonly() && size <= 0 ) buttons.add(erase);
    buttons.add(open);
    if ( !isReadonly() ) buttons.add(apply);
    buttons.add(ok);
    cancel.addActionListener(this);
    if ( !isReadonly() ) reset.addActionListener(this);
    if ( !isReadonly() ) sel.addActionListener(this);
    if ( !isReadonly() && size <= 0 ) insert.addActionListener(this);
    if ( !isReadonly() && size <= 0 ) erase.addActionListener(this);
    if ( !isReadonly() && nullable ) setnull.addActionListener(this);
    open.addActionListener(this);
    if ( !isReadonly() ) apply.addActionListener(this);
    ok.addActionListener(this);
    getContentPane().add(buttons, BorderLayout.SOUTH);
    fixButtons();

    setTitle("Reference Vector");
    setupFrame(600,300);
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
    size = Integer.parseInt((String)input.remove(0));
    int nref = Integer.parseInt((String)input.remove(0));
    if ( current == null ) current = new Vector();
    else current.clear();
    if ( values == null ) values = new Vector();
    else values.clear();
    for ( int iref = 0; iref < nref; ++iref ) {
      values.add(new RefRef(iref, (String)input.getFirst()));
      current.add(new RefRef(iref, (String)input.remove(0)));
    }
    current.add(new RefRef(nref, " - end - ", "NULL"));
    values.add(new RefRef(nref, " - end - ", "NULL"));
    int sel = selector.getSelectedIndex();
    selector.setListData(current);
    if ( sel >= 0 && sel < current.size() ) selector.setSelectedIndex(sel);

    return true;
  }

  protected void setValue() {

    for ( int iref = 0; iref < current.size() - 1; ++iref ) {
      String cmd = "set " + getObject() + ":" + getName() + "[" +
	Long.toString(iref) + "] " + ((RefRef)current.get(iref)).getFullName();
      exec(cmd);
    }
    reset();
  }

  protected void insert() {
    if ( selector.getSelectedIndex() < 0 )
      selector.setSelectedIndex(current.size() - 1);
    String newname = "NULL";
    ObjectSelector ch = new ObjectSelector(owner, dir, pclass, this);
    if ( ch.selected() == null ) {
      if ( !nullable )
	newname = ((RefRef)selector.getSelectedValue()).getFullName();
    } else {
      newname = ch.selected().getFullName();
    }
    if ( newname != "NULL" || nullable )
      insert(selector.getSelectedIndex(), newname);
  }

  protected void erase() {
    if ( selector.getSelectedIndex() >= 0 ) erase(selector.getSelectedIndex());
  }

  public String selectedValue() {
    RefRef sel = (RefRef)selector.getSelectedValue();
    if ( sel == null ) return null;
    return sel.getFullName();
  }

  public void fixButtons() {
    dir = "";
    String val = selectedValue();
    if ( val != null && val.lastIndexOf("/") >= 0 )
	dir = val.substring(0, val.lastIndexOf("/") + 1);
    if ( dir.equals("") ) dir = object.substring(0,object.lastIndexOf("/") + 1);
    if ( dir.equals("") ) dir = "/";
    if ( isReadonly() ) return;
    boolean changed = false;
    for ( int iref = 0; iref < current.size(); ++iref )
      if ( !(((RefRef)current.get(iref)).getFullName().
	     equals(((RefRef)current.get(iref)).getFullName())) )
	changed = true;
    apply.setEnabled(changed);
    reset.setEnabled(changed);
    sel.setEnabled(false);
    open.setEnabled(false);
    setnull.setEnabled(false);
    erase.setEnabled(false);
    insert.setEnabled(size <= 0 && !changed);
    if ( selector.getSelectedIndex() >= 0 &&
	 selector.getSelectedIndex() < current.size() - 1 ) {
      sel.setEnabled(true);
      open.setEnabled(selectedValue() != null && selectedValue() != "NULL");
      setnull.setEnabled(nullable);
      erase.setEnabled(size <= 0 && !changed);
    }
      
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
      ((RefRef)selector.getSelectedValue()).set("NULL", "NULL");
    }
    else if ( e.getSource() == sel ) {
      ObjectSelector ch = new ObjectSelector(owner, dir, pclass, this);
      if ( ch.selected() != null )
	((RefRef)selector.getSelectedValue()).set(ch.selected().getFullName());
    }
    else if ( e.getSource() == open ) {
      owner.openObject(((RefRef)selector.getSelectedValue()).getFullName());
    }
    else if ( e.getSource() == reset ) {
      reset();
    }
    else if ( e.getSource() == insert ) {
      insert();
    }
    else if ( e.getSource() == erase ) {
      erase();
    }
    fixButtons();
  }

  public void valueChanged(ListSelectionEvent e) {
    fixButtons();
  }    

}


