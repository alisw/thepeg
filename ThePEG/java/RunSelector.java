package ThePEG;

import javax.swing.*;
import java.util.*;
import java.awt.*;
import java.awt.event.*;

public class RunSelector extends JDialog
                           implements MouseListener,
				      ActionListener {
  Object selectedRun;
  boolean done = false;
  JList runList = new JList();
  JButton cancelButton = new JButton("Cancel");
  JButton okButton = new JButton("Select");
  SetupThePEG thepeg;

  public RunSelector(SetupThePEG thepeg) {
    super(thepeg, "Choose a run", true);
    this.thepeg = thepeg;
    setDefaultCloseOperation(DISPOSE_ON_CLOSE);
    LinkedList ret = thepeg.exec("lsruns");

    Vector runs = new Vector();
    while ( ret.size() > 0 ) runs.add(thepeg.getRun((String)ret.removeFirst()));
    runs.add("Create a new run");
    
    runList.setListData(runs);
    runList.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
    runList.addMouseListener(this);

    getContentPane().setLayout(new BorderLayout());
    getContentPane().add(new JScrollPane(runList), BorderLayout.CENTER);

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

  public Object selected() {
    if ( selectedRun == null || !done ) return null;
    return selectedRun;
  }

  public void actionPerformed(ActionEvent e) {
    if ( e.getSource() == cancelButton ) {
      selectedRun = null;
      dispose();
    }
    if ( e.getSource() == okButton ) {
      if ( selectedRun != null ) {
	done = true;
	dispose();
      } else if ( runList.getSelectedValue() != null ) {
	selectedRun = runList.getSelectedValue();
	done = true;
	dispose();
      }
    }
  }

  public void mouseClicked(MouseEvent e) {
    if ( e.getSource() == runList && e.getClickCount() >= 2 ) {
      selectedRun = runList.getSelectedValue();
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
