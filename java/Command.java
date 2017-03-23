package ThePEG;

import javax.swing.*;
import javax.swing.event.*;
import java.util.*;
import java.awt.*;
import java.awt.event.*;

public class Command extends Interface
  implements ActionListener {

  JButton ok = new JButton("Ok");
  JButton cancel = new JButton("Cancel");
  JButton apply = new JButton("Apply");
  JTextField text;

  public Command(SetupThePEG own, ObjectFrame obj, LinkedList input) {
    super(own, obj, input);
    if ( !setup(input) ) {
      JOptionPane.showMessageDialog(own, "Could not create Command view",
				    "Error", JOptionPane.ERROR_MESSAGE);
      return;
    }

    getContentPane().setLayout(new BorderLayout());
    getContentPane().add(getDescriptionArea(), BorderLayout.CENTER);
    JPanel buttons = new JPanel();
    buttons.add(cancel);
    buttons.add(apply);
    buttons.add(ok);
    cancel.addActionListener(this);
    apply.addActionListener(this);
    ok.addActionListener(this);
    getContentPane().add(buttons, BorderLayout.SOUTH);
    JPanel textpanel = new JPanel();
    textpanel.add(new JLabel("Enter arguments: "));
    textpanel.add(text);
    getContentPane().add(textpanel, BorderLayout.NORTH);
    text.addActionListener(this);

    setTitle("Command");
    setupFrame(500,150);
    
  }

  protected boolean setup(LinkedList input) {
    if ( !super.setup(input) ) return false;
    if ( text == null ) text = new JTextField("", 30);
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
    else if ( e.getSource() == ok || e.getSource() == text) {
      if ( !readonly ) setValue();
      dispose();
    }
    else if ( e.getSource() == apply ) {
      setValue();
    }
  }

  public static void classcheck() {}

}


