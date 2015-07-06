package ThePEG;

import javax.swing.*;
import java.util.*;
import java.awt.*;
import java.awt.event.*;

public class ObjectSelector extends JDialog
                           implements MouseListener,
				      ActionListener {
  ObjectRef selectedObject;
  boolean done = false;
  BrowserTree tree;
  JButton cancelButton = new JButton("Cancel");
  JButton okButton = new JButton("Select");
  SetupThePEG thepeg;

  public ObjectSelector(SetupThePEG thepeg, String dir,
		       String selclass, JDialog owner) {
    super(owner, "Choose an object", true);
    init(thepeg, dir, selclass);
  }

  public ObjectSelector(SetupThePEG thepeg, String dir, String selclass) {
    super(thepeg, "Choose an object", true);
    init(thepeg, dir, selclass);
  }

  public void init(SetupThePEG thepeg, String dir, String selclass) {
    this.thepeg = thepeg;
    setDefaultCloseOperation(DISPOSE_ON_CLOSE);
    tree = new BrowserTree(thepeg, selclass, true, dir);
    tree.addMouseListener(this);

    getContentPane().setLayout(new BorderLayout());
    getContentPane().add(tree, BorderLayout.CENTER);

    JPanel panel = new JPanel();
    cancelButton.addActionListener(this);
    okButton.addActionListener(this);
    panel.add(cancelButton);
    panel.add(okButton);
    getContentPane().add(panel, BorderLayout.SOUTH);

    setSize(300,400);

    thepeg.setLocation(this);

    setVisible(true);

  }

  public ObjectRef selected() {
    if ( !done ) return null;
    return selectedObject;
  }

  public void actionPerformed(ActionEvent e) {
    if ( e.getSource() == cancelButton ) {
      selectedObject = null;
      dispose();
    }
    if ( e.getSource() == okButton ) {
      if ( selectedObject != null ) {
	done = true;
	dispose();
      } else if ( tree.getSelectedObject() != null ) {
	selectedObject = tree.getSelectedObject();
	done = true;
	dispose();
      }
    }
  }

  public void mouseClicked(MouseEvent e) {
    if ( e.getClickCount() >= 2 ) {
      selectedObject = tree.getSelectedObject();
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
