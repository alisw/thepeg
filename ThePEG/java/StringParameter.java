package ThePEG;

import javax.swing.*;
import javax.swing.event.*;
import java.util.*;
import java.awt.*;
import java.awt.event.*;

public class StringParameter extends Interface
  implements ActionListener {

  String current;
  String def;

  JButton ok = new JButton("Ok");
  JButton cancel = new JButton("Cancel");
  JButton apply = new JButton("Apply");
  JButton reset = new JButton("Reset");
  JButton setdef = new JButton("Default");
  JButton browse = new JButton("Browse");
  JTextField text;

  static final int NOFILE = 0;
  static final int FILE = 1;
  static final int DIRECTORY = 2;

  int filetype = 0;

  public StringParameter(SetupThePEG own, ObjectFrame obj, LinkedList input,
			 int filetype) {
    super(own, obj, input);
    this.filetype = filetype;
    if ( !setup(input) ) {
      JOptionPane.showMessageDialog(own, "Could not create Parameter view",
				    "Error", JOptionPane.ERROR_MESSAGE);
      return;
    }

    getContentPane().setLayout(new BorderLayout());
    getContentPane().add(getDescriptionArea(), BorderLayout.CENTER);
    if ( isReadonly() ) {
      JLabel lab = new JLabel("");
      lab.setText("Readonly value is '" + current + "'");
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
      JPanel textpanel = new JPanel();
      textpanel.add(new JLabel("Value: "));
      textpanel.add(text);
      if ( filetype != NOFILE ) {
	browse.addActionListener(this);
	textpanel.add(browse);
      }
      getContentPane().add(textpanel, BorderLayout.NORTH);
      text.addActionListener(this);
      fixButtons();
    }

    setTitle("Parameter");
    setupFrame(500,150);
    
  }

  protected boolean setup(LinkedList input) {
    if ( !super.setup(input) ) return false;
    current = (String)input.remove(0);
    String next = (String)input.remove(0);
    while ( !next.equals("-inf") && input.size() > 0 ) {
      current += "\n" + next;
      next = (String)input.remove(0);
    }
    def = (String)input.remove(0);

    if ( text == null )
      text = new JTextField(current, 30);
    else
      text.setText(current);
    text.setCaretPosition(0);
    return true;
  }

  private void setValue() {
    setValue(text.getText());
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
      text.setText(def);
      text.setCaretPosition(0);
    }
    else if ( e.getSource() == reset ) {
      text.setText(current);
      text.setCaretPosition(0);
    }
    else if ( e.getSource() == browse ) {
      if ( filetype == NOFILE ) return;
      JFileChooser fc = new JFileChooser(".");
      fc.setDialogType(JFileChooser.OPEN_DIALOG);
      if ( filetype == FILE ) {
	fc.setFileSelectionMode(JFileChooser.FILES_ONLY);
	fc.setDialogTitle("Select file");
      } else if ( filetype == DIRECTORY ) {
	fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
	fc.setDialogTitle("Select directory");
      }
      if ( fc.showDialog(this, "Select") != JFileChooser.APPROVE_OPTION )
	return;
      text.setText(fc.getSelectedFile().getAbsolutePath());
    }
    fixButtons();
  }

  public void fixButtons() {
    String val = text.getText();
    setdef.setEnabled(!val.equals(def));
    reset.setEnabled(!val.equals(current));
    apply.setEnabled(!val.equals(current));
  }

  public static void classcheck() {}

}

