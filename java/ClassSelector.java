package ThePEG;

import javax.swing.*;
import java.util.*;
import java.awt.*;
import java.awt.event.*;

public class ClassSelector extends JDialog
                           implements MouseListener,
				      ActionListener {
  String selectedClass;
  boolean done = false;
  JList classList = new JList();
  JButton cancelButton = new JButton("Cancel");
  JButton okButton = new JButton("Select");
  SetupThePEG thepeg;

  public ClassSelector(SetupThePEG thepeg) {
    super(thepeg, "Choose a class", true);
    init(thepeg, "ThePEG::Interfaced", null);
  }

  public ClassSelector(SetupThePEG thepeg, String cl) {
    super(thepeg, "Choose a class", true);
    init(thepeg, cl, null);
  }

  public ClassSelector(SetupThePEG thepeg, JDialog owner, String cl) {
    super(owner, "Choose a class", true);
    init(thepeg, cl, owner);
  }

  public void init(SetupThePEG thepeg, String cl, JDialog owner) {
    this.thepeg = thepeg;
    setDefaultCloseOperation(DISPOSE_ON_CLOSE);
    LinkedList ret = thepeg.exec("lsclass " + cl);
    if ( ret.size() <= 0 ) {
      JOptionPane.showMessageDialog(this, "No classes found.", "Error",
				    JOptionPane.ERROR_MESSAGE);
      dispose();
      return;
    }

    classList.setListData(ret.toArray());
    classList.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
    classList.addMouseListener(this);

    getContentPane().setLayout(new BorderLayout());
    getContentPane().add(new JScrollPane(classList), BorderLayout.CENTER);

    JPanel panel = new JPanel();
    cancelButton.addActionListener(this);
    okButton.addActionListener(this);
    panel.add(cancelButton);
    panel.add(okButton);
    getContentPane().add(panel, BorderLayout.SOUTH);

    setSize(300,300);
    thepeg.setLocation(this);

    setVisible(true);

  }

  public String selected() {
    if ( selectedClass == null || !done ) return null;
    return selectedClass;
  }

  public void actionPerformed(ActionEvent e) {
    if ( e.getSource() == cancelButton ) {
      selectedClass = null;
      dispose();
    }
    if ( e.getSource() == okButton ) {
      if ( selectedClass != null ) {
	done = true;
	dispose();
      } else if ( classList.getSelectedValue() != null ) {
	selectedClass = (String)classList.getSelectedValue();
	done = true;
	dispose();
      }
    }
  }

  public void mouseClicked(MouseEvent e) {
    if ( e.getSource() == classList && e.getClickCount() >= 2 ) {
      selectedClass = (String)classList.getSelectedValue();
      done = true;
      dispose();
    }
  }

  public void mouseEntered(MouseEvent e) {}
  public void mouseExited(MouseEvent e) {}
  public void mousePressed(MouseEvent e) {}
  public void mouseReleased(MouseEvent e) {}

  public void dispose() {
    thepeg.removeLocation(this);
    super.dispose();
  }

  public static void classcheck() {}

}
