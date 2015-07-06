package ThePEG;

import javax.swing.*;
import javax.swing.event.*;
import java.util.*;
import java.math.*;
import java.awt.*;
import java.awt.event.*;

public class ParVector extends Interface
  implements ActionListener, ListSelectionListener, ChangeListener {

  Vector def;
  Vector min;
  Vector max;
  Vector current;
  Vector val;
  int sel = -1;
  boolean integer;
  String defdef;
  int size;

  JButton ok = new JButton("Ok");
  JButton cancel = new JButton("Cancel");
  JButton apply = new JButton("Apply");
  JButton reset = new JButton("Reset");
  JButton setdef = new JButton("Default");
  JButton insert = new JButton("Insert");
  JButton erase = new JButton("Erase");
  JList selector = new JList();
  FullSlider slider;
  JLabel rolabel = new JLabel("No value selected");

  class Entry {
    public String s;
    public Entry(String in) {
      s = in;
    }
    public String toString() {
      return s;
    }
    public void set(String in) {
      s = in;
    }
  }

  public ParVector(SetupThePEG own, ObjectFrame obj,
		   LinkedList input, boolean isint) {
    super(own, obj, input);
    if ( !setup(input) ) {
      JOptionPane.showMessageDialog(own,
				    "Could not create Parameter Vector view",
				    "Error", JOptionPane.ERROR_MESSAGE);
      return;
    }

    getContentPane().setLayout(new BorderLayout());
    JSplitPane split = new JSplitPane();
    split.setRightComponent(getDescriptionArea());
    selector.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
    selector.addListSelectionListener(this);
    split.setLeftComponent(selector);
    getContentPane().add(split, BorderLayout.CENTER);
    rolabel.setBorder(BorderFactory.createEmptyBorder(2,2,2,2));
    slider = new FullSlider(-1.0, 0.0, 1.0, 0.0, false);
    slider.setEnabled(false);
    slider.setOwner(this);
    getContentPane().add(slider, BorderLayout.NORTH);
    if ( isReadonly() ) {
      ok.addActionListener(this);
      getContentPane().add(ok, BorderLayout.SOUTH);
    } else {
      JPanel buttons = new JPanel();
      buttons.add(cancel);
      buttons.add(reset);
      buttons.add(setdef);
      if ( size <= 0 ) buttons.add(insert);
      if ( size <= 0 ) buttons.add(erase);
      buttons.add(apply);
      buttons.add(ok);
      cancel.addActionListener(this);
      reset.addActionListener(this);
      setdef.addActionListener(this);
      if ( size <= 0 ) insert.addActionListener(this);
      if ( size <= 0 ) erase.addActionListener(this);
      apply.addActionListener(this);
      ok.addActionListener(this);
      getContentPane().add(buttons, BorderLayout.SOUTH);
      fixButtons();
    }
    setTitle("Parameter Vector");
    setupFrame(600,300);   
  }

  protected boolean setup(LinkedList input) {
    if ( !super.setup(input) ) return false;
    try {
      size = Integer.parseInt((String)input.remove(0));
      int nvals = Integer.parseInt((String)input.remove(0));
      if ( current == null ) {
	current = new Vector();
	def = new Vector();
	min = new Vector();
	max = new Vector();
	val = new Vector();
      } else {
	def.clear();
	min.clear();
	max.clear();
	val.clear();
	current.clear();
      }
      for ( int i = 0; i < nvals; ++i ) {
	val.add(new Entry((String)input.getFirst()));
	current.add((String)input.remove(0));
	min.add((String)input.remove(0));
	def.add((String)input.remove(0));
	max.add((String)input.remove(0));
      }
      defdef = (String)input.remove(0);
    }
    catch ( NumberFormatException e ) {
      return false;
    }
    val.add(new Entry(" - end - "));
    sel = -1;
    int ssel = selector.getSelectedIndex();
    selector.setListData(val);
    if ( ssel >= 0 && ssel < current.size() ) selector.setSelectedIndex(ssel);
    return true;
  }

  public void pushValue() {
    if ( sel < 0 || sel >= current.size() ) return;
    ((Entry)val.get(sel)).set(integer? Long.toString(slider.getInt()):
			      Double.toString(slider.getDouble()));
    fixButtons();
    selector.repaint();
  }
    

  protected void setValue() {
    for ( int i = 0; i < current.size(); ++i ) {
      String cmd = "set " + getObject() + ":" + getName() + "[" +
	Long.toString(i) + "] " + val.get(i);
      exec(cmd);
    }
    reset();
  }

  protected void insert() {
    int ssel = selector.getSelectedIndex();
    if ( ssel < 0 || ssel > current.size() ) {
      ssel = current.size();
      selector.setSelectedIndex(ssel);
    }
    insert(ssel, defdef);
  }

  protected void erase() {
    int ssel = selector.getSelectedIndex();
    if ( ssel >= 0 || ssel < current.size() ) erase(ssel);
  }

  public void fixButtons() {
    int ssel = selector.getSelectedIndex();
    if ( isReadonly() ) return;
    boolean changed = false;
    for ( int i = 0; i < current.size(); ++i )
      if ( !((String)current.get(i)).equals(((Entry)val.get(i)).toString()) )
	changed = true;
    apply.setEnabled(changed);
    reset.setEnabled(changed);
    insert.setEnabled(ssel >= 0 && ssel < val.size() );
    erase.setEnabled(ssel >= 0 && ssel < current.size() );
    setdef.setEnabled(false);
    if ( ssel >= 0 && ssel < current.size() &&
	 !((String)def.get(ssel)).equals(((Entry)val.get(ssel)).toString()) )
      setdef.setEnabled(true);
  }


  public void actionPerformed(ActionEvent e) {
    if ( e.getSource() == cancel ) {
      dispose();
    }
    else if ( e.getSource() == ok ) {
      if ( !readonly ) setValue();
      dispose();
    }
    else if ( e.getSource() == apply ) {
      setValue();
    }
    else if ( e.getSource() == setdef ) {
      setDefault(selector.getSelectedIndex());
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

  public void setDefault(int ssel) {
    if ( ssel < 0 || ssel >= current.size() ) return;
    ((Entry)val.get(ssel)).set((String)def.get(ssel));
    selectSlider(ssel);
  }    

  public void selectSlider(int ssel) {
    if ( ssel >= 0 && ssel < current.size() ) {
      try {
	double dval = Double.parseDouble(((Entry)val.get(ssel)).toString());
	double ddef = Double.parseDouble((String)def.get(ssel));
	String cmin = (String)min.get(ssel);
	String cmax = (String)max.get(ssel);
	sel = ssel;
	slider.set(cmin, dval, cmax, ddef);
	slider.setEnabled(!isReadonly());
      }
      catch ( NumberFormatException ex ) {
	return;
      }
    } else {
	slider.setEnabled(false);
	sel = -1;
    }
    fixButtons();
  }



  public void valueChanged(ListSelectionEvent e) {
    if ( e.getSource() == selector )
      selectSlider(selector.getSelectedIndex());
  }

  public void stateChanged(ChangeEvent e) {
    fixButtons();
  }

  public static void classcheck() {}

}
