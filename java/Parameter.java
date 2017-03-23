package ThePEG;

import javax.swing.*;
import javax.swing.event.*;
import java.util.*;
import java.awt.*;
import java.awt.event.*;

public class Parameter extends Interface
  implements ActionListener, ChangeListener {

  double def;
  double min;
  double max;
  double current;
  boolean integer;

  JButton ok = new JButton("Ok");
  JButton cancel = new JButton("Cancel");
  JButton apply = new JButton("Apply");
  JButton reset = new JButton("Reset");
  JButton setdef = new JButton("Default");
  FullSlider slider;

  public Parameter(SetupThePEG own, ObjectFrame obj, LinkedList input, boolean in) {
    super(own, obj, input);
    integer = in;
    if ( !setup(input) ) {
      JOptionPane.showMessageDialog(own, "Could not create Parameter view",
				    "Error", JOptionPane.ERROR_MESSAGE);
      return;
    }

    getContentPane().setLayout(new BorderLayout());
    getContentPane().add(getDescriptionArea(), BorderLayout.CENTER);
    if ( isReadonly() ) {
      JLabel lab = new JLabel("");
      if ( integer )
	lab.setText("Readonly value is " + slider.getInt());
      else
	lab.setText("Readonly value is " + slider.getDouble());
      lab.setBorder(BorderFactory.createEmptyBorder(2,2,2,2));
      getContentPane().add(lab, BorderLayout.NORTH);
      ok.addActionListener(this);
      getContentPane().add(ok, BorderLayout.SOUTH);
    } else {
      JPanel buttons = new JPanel();
      buttons.add(cancel);
      buttons.add(reset);
      buttons.add(setdef);
      buttons.add(apply);
      buttons.add(ok);
      cancel.addActionListener(this);
      reset.addActionListener(this);
      setdef.addActionListener(this);
      apply.addActionListener(this);
      ok.addActionListener(this);
      getContentPane().add(buttons, BorderLayout.SOUTH);
      getContentPane().add(slider, BorderLayout.NORTH);
      slider.addChangeListener(this);
      fixButtons();
    }

    setTitle("Parameter");
    setupFrame(500,150);
    
  }

  protected boolean setup(LinkedList input) {
    if ( !super.setup(input) ) return false;
    try {
      current = Double.parseDouble((String)input.remove(0));
      String cmin = (String)input.remove(0);
      def = Double.parseDouble((String)input.remove(0));
      String cmax = (String)input.remove(0);

      if ( slider == null )
	slider = new FullSlider(cmin, current, cmax, def, integer);
      else
	slider.set(cmin, current, cmax, def);
    }
    catch ( NumberFormatException e ) {
      return false;
    }
    return true;
  }

  private void setValue() {
    setValue(integer? Long.toString(slider.getInt()):
      Double.toString(slider.getDouble()));
    reset();
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
      slider.set(min, def, max, def);
    }
    else if ( e.getSource() == reset ) {
      slider.set(min, current, max, def);
    }
    fixButtons();
  }

  public void stateChanged(ChangeEvent e) {
    fixButtons();
  }

  public void fixButtons() {
    double val = slider.getDouble();
    setdef.setEnabled(val != def);
    reset.setEnabled(val != current);
    apply.setEnabled(val != current);
  }

  public static void classcheck() {}

}

