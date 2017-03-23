package ThePEG;

import javax.swing.*;
import java.util.*;
import java.awt.*;
import java.awt.event.*;

public class Switch extends Interface implements ActionListener {

  long def;
  long current;

  JButton ok = new JButton("Ok");
  JButton cancel = new JButton("Cancel");
  JButton apply = new JButton("Apply");
  JButton reset = new JButton("Reset");
  JButton setdef = new JButton("Default");
  JComboBox selector = new JComboBox();
  HashMap options = new HashMap();

  public Switch(SetupThePEG own, ObjectFrame obj, LinkedList input) {
    super(own, obj, input);
    if ( !setup(input) ) {
      JOptionPane.showMessageDialog(own, "Could not create Switch view",
				    "Error", JOptionPane.ERROR_MESSAGE);
      return;
    }

    getContentPane().setLayout(new BorderLayout());
    getContentPane().add(getDescriptionArea(), BorderLayout.CENTER);
    if ( isReadonly() ) {
      SwitchOption sel = (SwitchOption)selector.getSelectedItem();
      getContentPane().add(sel, BorderLayout.NORTH);
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
      selector.addActionListener(this);
      getContentPane().add(buttons, BorderLayout.SOUTH);
      getContentPane().add(selector, BorderLayout.NORTH);
      fixButtons();
    }
    setTitle("Switch");
    setupFrame(500,150);
    
  }

  protected boolean setup(LinkedList input) {
    if ( !super.setup(input) ) return false;
    try {
      current = Long.parseLong((String)input.remove(0));
      def = Long.parseLong((String)input.remove(0));
      long nopt = Long.parseLong((String)input.remove(0));
      options.clear();
      selector.removeAllItems();
      for ( long iopt = 0; iopt < nopt; ++iopt ) {
	long val = Long.parseLong((String)input.remove(0));
	SwitchOption swopt = new SwitchOption(val, (String)input.remove(0),
					      (String)input.remove(0));
	selector.addItem(swopt);
	options.put(new Long(val), swopt);
      }
      SwitchOption sel = (SwitchOption)options.get(new Long(current));
      selector.setSelectedItem(sel);
      selector.setToolTipText(SetupThePEG.ttt + sel.getComment());
    }
    catch ( NumberFormatException e ) {
      return false;
    }
    return true;
  }

  public long getSelected() {
    SwitchOption opt = (SwitchOption)(selector.getSelectedItem());
    return opt.getIndex();
  }

  private void setValue() {
    setValue(Long.toString(getSelected()));
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
      selector.setSelectedItem(options.get(new Long(def)));
    }
    else if ( e.getSource() == reset ) {
      selector.setSelectedItem(options.get(new Long(current)));
    }
    else if ( e.getSource() == selector ) {
      SwitchOption sel = (SwitchOption)selector.getSelectedItem();
      if ( sel != null ) selector.setToolTipText(sel.getComment());
    }
    fixButtons();
  }

  public void fixButtons() {
    if ( selector.getSelectedIndex() < 0 ) return;
    apply.setEnabled(getSelected() != current);
    reset.setEnabled(getSelected() != current);
    setdef.setEnabled(getSelected() != def);
  }

  public static void classcheck() {}

}

