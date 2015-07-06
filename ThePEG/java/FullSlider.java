package ThePEG;

import javax.swing.*;
import javax.swing.event.*;
import java.util.*;
import java.text.*;
import java.awt.*;
import java.awt.event.*;

public class FullSlider extends JPanel
                        implements ActionListener, ChangeListener {

  double val;
  double max;
  double min;
  double cmax;
  double cmin;
  double def;
  boolean integer;
  JTextField value = new JTextField(8);
  JLabel minl = new JLabel("");
  JLabel maxl = new JLabel("");
  JSlider slider = new JSlider(0, 1000);
  DecimalFormat nf = new DecimalFormat();
  ParVector owner = null;

  public FullSlider(double mi, double va, double ma, double de, boolean in) {
    integer = in;
    set(mi, va, ma, de);
    slider.addChangeListener(this);
    value.addActionListener(this);
    add(minl);
    add(slider);
    add(maxl);
    add(value);
    slider.addMouseListener(new MouseAdapter() {
	public void mouseReleased(MouseEvent e) {
	  setValue(cmin + slider.getValue()*(cmax - cmin)*0.001);
	}
      });
  }

  public FullSlider(String mi, double va, String ma, double de, boolean in) {
    integer = in;
    set(mi, va, ma, de);
    slider.addChangeListener(this);
    value.addActionListener(this);
    add(minl);
    add(slider);
    add(maxl);
    add(value);
    slider.addMouseListener(new MouseAdapter() {
	public void mouseReleased(MouseEvent e) {
	  setValue(cmin + slider.getValue()*(cmax - cmin)*0.001);
	}
      });
  }

  public void set(String cmi, double va, String cma, double de) {
    double mi = 0.0;
    double ma = 0.0;
    try {
      if ( cmi.equals("-inf") ) mi = de + Math.abs(de) + 1.0;
      else mi = Double.parseDouble(cmi);
      if ( cma.equals("inf") ) ma = de - Math.abs(de) - 1.0;
      else ma = Double.parseDouble(cma);
    }
    catch ( NumberFormatException e ) {}
    set(mi, va, ma, de);
  }

  public void set(double mi, double va, double ma, double de) {
    min = mi;
    max = ma;
    val = va;
    def = de;
    setValue(va);
  }

  public void setEnabled(boolean en) {
    slider.setEnabled(en);
    value.setEditable(en);
  }

  private void setValue(double v) {
    if ( ( min > def || val >= min ) && ( max < def || val <= max ) ) val = v;
    double scale = Math.max(Math.abs(def), Math.abs(val));
    if ( scale == 0.0 && min <= def ) scale = Math.abs(min)/10.0;
    if ( scale == 0.0 && max >= def ) scale = Math.abs(max)/10.0;
    if ( scale == 0.0 ) scale = 1.0;
    cmax = max;
    if ( max < def ) cmax = scale*10.0;
    else if ( val != 0.0 ) cmax = Math.min(scale*10.0, max);
    if ( min > def ) cmin = -scale*10.0;
    else if ( val != 0.0 ) cmin = Math.max(-scale*10.0, min);
    slider.setValue((int)(1000.0*(val - cmin)/(cmax - cmin)));
    if ( integer ) {
      if ( min > def )
	minl.setText("-inf [" + Math.round(cmin) + "]");
      else
	minl.setText("" + Math.round(min) + " [" + Math.round(cmin) + "]");
      if ( max < def )
	maxl.setText("[" + Math.round(cmax) + "] inf");
      else
	maxl.setText("[" + Math.round(cmax) + "] " + Math.round(max));
    } else {
      if ( min > def )
	minl.setText("-inf [" + format(cmin) + "]");
      else
	minl.setText("" + format(min) + " [" + format(cmin) + "]");
      if ( max < def )
	maxl.setText("[" + format(cmax) + "] inf");
      else
	maxl.setText("[" + format(cmax) + "] " + nf.format(max));
    }
    if ( integer )
      value.setText("" + getInt());
    else
      value.setText("" + format(val));
    value.setCaretPosition(0);
    if ( owner != null ) owner.pushValue();
  }

  public void setOwner(ParVector pv) {
    owner = pv;
  }

  public String format(double d) {
    if ( d == 0 ) return "0";
    if ( Math.abs(d) < 0.1 ) {
      nf.applyPattern("0.####E0");
    }
    else if ( Math.abs(d) < 1 ) {
      nf.applyPattern("0.####");
    }
    else if ( Math.abs(d) < 10 ) {
      nf.applyPattern("0.###");
    }
    else if ( Math.abs(d) < 100 ) {
      nf.applyPattern("00.##");
    }
    else if ( Math.abs(d) < 1000 ) {
      nf.applyPattern("000.#");
    }
    else if ( Math.abs(d) < 10000 ) {
      nf.applyPattern("0000");
    }
    else {
      nf.applyPattern("0.####E0");
    }
    return nf.format(d);
  }

  public long getInt() {
    return Math.round(val);
  }

  public double getDouble() {
    return val;
  }

  public void changeText() {
    double v = cmin + slider.getValue()*(cmax - cmin)*0.001;
    if ( integer )
      value.setText("" + Math.round(v));
    else
      value.setText("" + format(v));
    value.setCaretPosition(0);
  }


  public void addActionListener(ActionListener l) {
    value.addActionListener(l);
  }

  public void addChangeListener(ChangeListener l) {
    slider.addChangeListener(l);
  }

  public void actionPerformed(ActionEvent e) {
    if ( e.getSource() == value ) {
      try {
	setValue(Double.parseDouble(value.getText()));
      }
      catch ( NumberFormatException ex ) {
	setValue(val);
      }
    }
  }

  public void stateChanged(ChangeEvent e) {
    if ( e.getSource() == slider ) changeText();
  }

  public static void classcheck() {}

}
