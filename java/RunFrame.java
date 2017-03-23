package ThePEG;

import java.io.*;
import java.util.*;
import java.awt.*;
import java.awt.event.*;
import javax.swing.*;

public class RunFrame extends JFrame
  implements ActionListener, Runnable {

  SetupThePEG thepeg;

  String name;

  JTextField dirfield = new JTextField(".", 25);
  JTextField hostfield = new JTextField("localhost", 25);
  JTextField nevefield = new JTextField("default", 25);
  JTextField argfield = new JTextField("", 25);
  JProgressBar progbar = new JProgressBar(0, 1000);
  JLabel progfield = new JLabel("Not running");
  JLabel etcfield = new JLabel("Not running");
  JPanel progpanel = new JPanel();
  CardLayout progcards = new CardLayout();

  JButton closebutton = new JButton("Close");
  JButton removebutton = new JButton("Remove");
  JButton runbutton = new JButton("Run");

  File runfile = null;
  boolean running = false;
  boolean done = false;
  boolean stop = false;
  long progpercent;

  Thread runthread = null;
  Process runthepeg;
  String errors;

  GridBagConstraints gc = new GridBagConstraints();

  public RunFrame(SetupThePEG thepeg, String name) {
    super("ThePEG run " + name);
    this.thepeg = thepeg;
    this.name = name;
    setIconImage(thepeg.icon);

    getContentPane().setLayout(new GridBagLayout());
    gc.fill = GridBagConstraints.HORIZONTAL;
    gc.insets = new Insets(2, 10, 2, 10);

    progbar.setStringPainted(true);
    progpanel.setLayout(progcards);
    progpanel.add(progbar, "running");
    progpanel.add(progfield, "stopped");
    progcards.show(progpanel, "stopped");

    add(new JLabel("Run name: ", JLabel.RIGHT), 0, 0, 1, 1.2);
    add(new JLabel("Run in directory: ", JLabel.RIGHT), 0, 1, 1, 1.2);
    add(new JLabel("Run on host: ", JLabel.RIGHT), 0, 2, 1, 1.2);
    add(new JLabel("Number of events: ", JLabel.RIGHT), 0, 3, 1, 1.2);
    add(new JLabel("Additional arguments: ", JLabel.RIGHT), 0, 4, 1, 1.2);
    add(new JLabel("Progress: ", JLabel.RIGHT), 0, 5, 1, 1.2);
    add(new JLabel("Estimated time left: ", JLabel.RIGHT), 0, 6, 1, 1.2);
	
    add(new JLabel(name, JLabel.LEFT), 1, 0, 1.99, 1);
    add(dirfield, 1, 1, 1.99, 1);
    add(hostfield, 1, 2, 1.99, 1);
    add(nevefield, 1, 3, 1.99, 1);
    add(argfield, 1, 4, 1.99, 1);
    add(progpanel, 1, 5, 1.99, 1);
    add(etcfield, 1, 6, 1.99, 1);

    dirfield.setToolTipText("The (existing) directory where the " +
			    "run should be done.");
    hostfield.setToolTipText(SetupThePEG.ttt +
			     "The host where the run should be done. " +
			     "Must be reachable with <i>ssh</i> without " + 
			     "password and must have ThePEG installed.");
    nevefield.setToolTipText(SetupThePEG.ttt +
			     "If not a number, the default number of events " +
			     "of the EventGenerator object will be used.");
    argfield.setToolTipText("Eg. \"-d 2\" for debug level.");
    
    JPanel buttons = new JPanel();
    buttons.add(closebutton);
    buttons.add(removebutton);
    buttons.add(runbutton);
    add(buttons, 0, 7, 2.99, 1.2);

    closebutton.addActionListener(this);
    removebutton.addActionListener(this);
    runbutton.addActionListener(this);

    pack();
  }

  private void startRun() {
    runthread = new Thread(this);
    runthread.start();

  }

  private void abortRun() {
    stop = true;
    if ( runthepeg != null ) runthepeg.destroy();
    try {
      runthread.join();
    } catch (InterruptedException e) {}
    setAborted();
  }

  private void setAborted() {
    progcards.show(progpanel, "stopped");
    etcfield.setText("Aborted");
    progfield.setText("Aborted");
  }    

  public void showError(String mess) {
    JOptionPane.showMessageDialog(this, mess, "Error",
				  JOptionPane.ERROR_MESSAGE);

  }

  public void run() {
    runbutton.setText("Abort");
    dirfield.setEnabled(false);
    hostfield.setEnabled(false);
    argfield.setEnabled(false);
    nevefield.setEnabled(false);
    running = true;
    done = false;
    progcards.show(progpanel, "running");
    etcfield.setText("Not known");

    Vector args = new Vector();
    args.add("runThePEG");
    args.add("-cat");
    try {
      runfile = File.createTempFile("thepeg-", ".run");
    }
    catch (IOException e) {
      JOptionPane.showMessageDialog(this, "Could not save run file.", "Error",
				    JOptionPane.ERROR_MESSAGE);
      return;
    }
    thepeg.exec("saverunfile " + getName() + " " + runfile.getPath());
    args.add(runfile.getPath());
    String host = hostfield.getText();
    if ( host.matches("\\w+") && !host.equals("localhost") ) {
      args.add("-host");
      args.add(host);
    }
    args.add("-dir");
    args.add(dirfield.getText());
    try {
      long i = Long.parseLong(nevefield.getText());
      if ( i > 0 ) {
	args.add("-N");
	args.add("" + i);
      }
    } catch ( NumberFormatException e) {}
    args.add("--tics");
    StringTokenizer st = new StringTokenizer(argfield.getText());
    while ( st.hasMoreTokens() ) args.add(st.nextToken());

    if ( thepeg.debug ) System.err.print("> Starting: ");
    String [] arg = new String[args.size()];
    for ( int i = 0; i < args.size(); ++i ) {
      arg[i] = args.get(i).toString();
      if ( thepeg.debug ) System.err.print( " " + arg[i]);
    }
    if ( thepeg.debug ) System.err.println();
    etcfield.setText("Unknown");
    try {
      runthepeg = Runtime.getRuntime().exec(arg);

      BufferedReader outStream =
	new BufferedReader(new InputStreamReader(runthepeg.getInputStream()));
      BufferedReader errStream =
	new BufferedReader(new InputStreamReader(runthepeg.getErrorStream()));
      errors = "";

      String s;
      Date t0 = null;
      Date tp = null;
      long ip = 0;
      long ieve = 0;
      long neve = 0;
      while ( ( s = errStream.readLine() ) != null ) {
	st = new StringTokenizer(s);
	if ( st.countTokens() >= 3 && st.nextToken().equals("tic>") ) {
	  try {
	    ieve = Long.parseLong(st.nextToken());
	    neve = Long.parseLong(st.nextToken());
	    progbar.setValue((int)((ieve*1000)/neve));
	    progpercent = (ieve*100)/neve;
	    if ( t0 == null ) {
	      t0 = new Date();
	      tp = t0;
	      ip = ieve;
	    } else if ( ieve != neve ) {
	      Date t = new Date();
	      long tn = (t.getTime() - t0.getTime())*(neve - ieve)/ieve;
	      long tnp = (t.getTime() - tp.getTime())*(neve - ieve)/(ieve - ip);
	      tn = ((neve - ieve)*tn + ieve*tnp)/neve;
	      setETC(tn);
	      tp = t;
	      ip = ieve;
	    }
	  } catch ( NumberFormatException e) {
	    showError(s);
	  }
	} else
	  showError(s);
      }
      if ( ieve == neve ) {
	progfield.setText("Done");
	progcards.show(progpanel, "stopped");
	etcfield.setText("Done");
      } else
	setAborted();
    } catch ( IOException ex ) {
      setAborted();
    }
    progbar.setValue(0);
    runbutton.setText("Rerun");
    dirfield.setEnabled(true);
    hostfield.setEnabled(true);
    argfield.setEnabled(true);
    nevefield.setEnabled(true);
    running = false;
    if ( !stop ) done = true;
    stop = false;
    runfile.delete();
    runfile = null;
    runthepeg = null;
  }

  public void setETC(long t) {
    String str;
    if ( t >= 3600000 ) {
      long h = t/3600000;
      str = "" + h + " hour";
      if ( h > 1 ) str += "s";
      long m = (t%3600000)/60000;
      if ( h < 10 && m > 0 ) {
	str += " and " + m + " minute";
	if ( m > 1 ) str += "s";
      }
    } else if ( t >= 60000 ) {
      long m = t/60000;
      str = "" + m + " minute";
      if ( m > 1 ) str += "s";
      long s = (t%60000)/1000;
      if ( m < 10 && s > 0 ) {
	str += " and " + s + " second";
	if ( s > 1 ) str += "s";
      }
    } else if ( t >= 1000 ) {
      long s = t/1000;
      str = "" + s + " second";
      if ( s > 1 ) str += "s";
    } else
      str = "0 seconds";
    etcfield.setText(str);
  }

  public void actionPerformed(ActionEvent e) {
    if ( e.getSource() == closebutton ) {
      dispose();
    }
    else if ( e.getSource() == removebutton ) {
      if ( running ) {
	if ( JOptionPane.showConfirmDialog
	     (this, "This run is currently executing. Removing it will " +
	      "abort the run. Really remove?", "Really Remove?",
	      JOptionPane.YES_NO_CANCEL_OPTION, JOptionPane.QUESTION_MESSAGE )
	     != JOptionPane.YES_OPTION ) return;
	abortRun();
      }
      thepeg.remove(this);
      dispose();
    }
    else if ( e.getSource() == runbutton ) {
      if ( running )
	abortRun();
      else
	startRun();
    }
  }

  public String toString() {
    if ( done ) return name + " (done)";
    else if ( running ) return name + " (" + progpercent + "%)";
    else return name;
  }

  public String getName() {
    return name;
  }

  public void add(Component item, int x, int y, double dx, double dy) {
    int ox = gc.gridx;
    gc.gridx = x;
    int oy = gc.gridy;
    gc.gridy = y;
    int odx = gc.gridwidth;
    gc.gridwidth = (int)dx;
    int ody = gc.gridheight;
    gc.gridheight = (int)dy;
    double owx = gc.weightx;
    gc.weightx = dx - gc.gridwidth;
    double owy = gc.weighty;
    gc.weighty = dy - gc.gridheight;
    getContentPane().add(item, gc);
    gc.gridx = ox;
    gc.gridy = oy;
    gc.gridwidth = odx;
    gc.gridheight = ody;
    gc.weightx = owx;
    gc.weighty = owy;
    

  }

  public static void classcheck() {}

}
