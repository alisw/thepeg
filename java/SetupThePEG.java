package ThePEG;

import javax.swing.filechooser.FileFilter;
import javax.swing.*;
import java.io.*;
import java.util.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.image.*;

public class SetupThePEG extends JFrame implements ActionListener {

  
  Process theProcess;
  PrintWriter cmdStream;
  BufferedReader inStream;
  BufferedReader errStream;
  BufferedReader stdin;

  BrowserTree tree;
  HashSet allTrees = new HashSet();

  JEditorPane helparea = new JEditorPane();
  String inithelp;

  TreeMap objects = new TreeMap();
  ButtonGroup showgroup = new ButtonGroup();
  JMenu showmenu = new JMenu("Show");

  Vector actions = new Vector();
  File actionFile;
  boolean actionsPerformed = false;

  JMenuItem loadrepo = new JMenuItem("Load New Repository ...");
  JMenuItem saverepo = new JMenuItem("Save Repository ...");

  JMenuItem running = new JMenuItem("Run ...");

  JMenuItem showonly = new JMenuItem("Show only class ...");

  JMenuItem readActions = new JMenuItem("Read Actions ...");
  JMenuItem saveActionsAs = new JMenuItem("Save Actions As ...");
  JMenuItem saveActions = new JMenuItem("Save Actions");

  JMenuItem appendpath = new JMenuItem("Append Library Path ...");
  JMenuItem prependpath = new JMenuItem("Prepend Library Path ...");
  JMenuItem loadlib = new JMenuItem("Load Dynamic Library ...");

  JMenuItem quit = new JMenuItem("Quit");

  boolean debug;

  static public final String ttt = "<html><body width=250>";

  ObjectRef clip = null;
  boolean clipIsClone = false;

  ObjectRef dragged = null;
  
  Vector locations = new Vector();

  JFileChooser fc = new JFileChooser(".");

  TreeMap runs = new TreeMap();

  FileFilter sofilter = new FileFilter() {
      public boolean accept(File file) {
	if ( file.isDirectory() ) return true;
	String name = file.getName();
	if ( name.lastIndexOf('.') < 0 ) return false;
	name = name.substring(name.lastIndexOf('.'), name.length());
	return name.equals(".so") || name.equals(".dylib");
      }
      public String getDescription() {
	return "Dynamic libraries";
      }
    };

  FileFilter infilter = new FileFilter() {
      public boolean accept(File file) {
	if ( file.isDirectory() ) return true;
	String name = file.getName();
	if ( name.lastIndexOf('.') < 0 ) return false;
	name = name.substring(name.lastIndexOf('.'), name.length());
	return name.equals(".in");
      }
      public String getDescription() {
	return "Files containing ThePEG setup commands";
      }
    };

  FileFilter rpofilter = new FileFilter() {
      public boolean accept(File file) {
	if ( file.isDirectory() ) return true;
	String name = file.getName();
	if ( name.lastIndexOf('.') < 0 ) return false;
	name = name.substring(name.lastIndexOf('.'), name.length());
	return name.equals(".rpo");
      }
      public String getDescription() {
	return "Files containing ThePEG repositories";
      }
    };

    protected static Image createImage() {
        //Create a 16x16 pixel image.
        BufferedImage bi = new BufferedImage(16, 16, BufferedImage.TYPE_INT_RGB);

        //Draw into it.
        Graphics g = bi.getGraphics();
        g.setColor(Color.DARK_GRAY);
        g.fillRect(0, 0, 16, 16);

        g.setColor(Color.LIGHT_GRAY);
	// T
        g.drawLine(1, 1, 3, 1);
        g.drawLine(2, 1, 2, 7);
	// h
        g.drawLine(5, 1, 5, 7);
        g.drawLine(5, 4, 7, 4);
        g.drawLine(7, 4, 7, 7);

	// e
        g.drawLine(9, 5, 9, 5);
        g.drawLine(10, 4, 10, 4);
        g.drawLine(10, 6, 11, 6);
        g.drawLine(11, 5, 11, 5);
        g.drawLine(10, 7, 11, 8);

        g.setColor(Color.RED);
	// P
	g.drawLine(1, 6, 1, 14);
	g.drawLine(1, 6, 3, 6);
	g.drawLine(1, 10, 3, 10);
	g.drawLine(4, 7, 4, 9);

	// E
	g.drawLine(6, 6, 6, 14);
	g.drawLine(6, 6, 9, 6);
	g.drawLine(6, 14, 9, 14);
	g.drawLine(6, 10, 8, 10);

	// G
	g.drawLine(11, 7, 11, 13);
	g.drawLine(12, 6, 13, 6);
	g.drawLine(14, 7, 14, 7);
	g.drawLine(12, 14, 13, 14);
	g.drawLine(13, 10, 14, 10);
	g.drawLine(14, 10, 14, 13);

        //Clean up.
        g.dispose();

        //Return it.
        return bi;
    }

  public Image icon = createImage();

  SetupThePEG(String [] args) {
    super("Setup ThePEG");

    if ( args.length <= 0 ) System.exit(1);
    String [] cmdarray = new String[args.length + 1];
    cmdarray[0] = args[0];
    cmdarray[1] = "--java";

    debug = false;

    for ( int i = 1; i < args.length; ++i ) {
      cmdarray[i + 1] = args[i];
      if ( args[i].length() >= 2 && args[i].substring(0,2).equals("-d") )
	debug = true;
    }
    setIconImage(icon);

    try {
      theProcess = Runtime.getRuntime().exec(cmdarray);
      cmdStream =
	new PrintWriter(new OutputStreamWriter(theProcess.getOutputStream()));
      inStream =
	new BufferedReader(new InputStreamReader(theProcess.getInputStream()));
      errStream =
	new BufferedReader(new InputStreamReader(theProcess.getErrorStream()));
      stdin =
	new BufferedReader(new InputStreamReader(System.in));
      
      JMenuBar menubar = new JMenuBar();
      JMenu file = new JMenu("File");
      loadrepo.addActionListener(this);
      file.add(loadrepo);
      saverepo.addActionListener(this);
      file.add(saverepo);
      file.add(new JSeparator());
      running.addActionListener(this);
      file.add(running);
      file.add(new JSeparator());
      saveActions.addActionListener(this);
      saveActions.setEnabled(false);
      file.add(saveActions);
      saveActionsAs.addActionListener(this);
      file.add(saveActionsAs);
      readActions.addActionListener(this);
      file.add(readActions);
      file.add(new JSeparator());
      loadlib.addActionListener(this);
      file.add(loadlib);
      appendpath.addActionListener(this);
      file.add(appendpath);
      prependpath.addActionListener(this);
      file.add(prependpath);
      file.add(new JSeparator());
      quit.addActionListener(this);
      file.add(quit);
      menubar.add(file);
      setJMenuBar(menubar);

      addWindowListener(new WindowAdapter() {
	  public void windowClosing(WindowEvent e) {
	    safeQuit();
	  }
	});

      inStream.readLine();
      getContentPane().setLayout(new BorderLayout());
      tree = new BrowserTree(this);

      JMenu edit = tree.getExternalMenu();
      showonly.addActionListener(this);
      edit.add(new JSeparator(), 0);
      showmenu.add(showonly);
      JRadioButtonMenuItem cr = new JRadioButtonMenuItem("All classes", true);
      showgroup.add(cr);
      showmenu.add(cr);
      cr.addActionListener(this);
      cr = new JRadioButtonMenuItem("Event generators", true);
      showgroup.add(cr);
      showmenu.add(cr);
      cr.addActionListener(this);

      menubar.add(edit);
      menubar.add(showmenu);

      createHelpText();
      helparea.setContentType("text/html");
      helparea.setText(inithelp);
      helparea.setEditable(false);
      helparea.setCaretPosition(0);
      JScrollPane s = new JScrollPane(helparea);
      s.setMinimumSize(new Dimension(300,250));
      tree.setMinimumSize(new Dimension(200,250));

      JSplitPane split = new JSplitPane();
      split.setLeftComponent(tree);
      split.setRightComponent(s);
      split.setOneTouchExpandable(true);

      getContentPane().add(split, BorderLayout.CENTER);


      setDefaultCloseOperation(DO_NOTHING_ON_CLOSE);
      setSize(600,450);
      setLocation(100,100);
      setVisible(true);
    }
    catch (IOException ex) {
      System.err.println("Could not start sub process. Exiting.");
      System.exit(1);
    }
  }

  public void lostProcess() {
    System.err.println("Lost contact with sub process. Exiting.");
    System.exit(2);
  }

  public synchronized LinkedList action(String cmd) {
    LinkedList ret = exec(cmd);
    if ( ret == null ) return ret;
    actions.add(cmd);
    actionsPerformed = true;
    saveActions.setEnabled(true);
    return ret;
  }

  public synchronized LinkedList exec(String cmd) {
    LinkedList ret = new LinkedList();
    boolean error = false;
    try {
      if ( debug ) System.err.println("> " + cmd);
      cmdStream.println(cmd);
      cmdStream.flush();
      String s;
      while ( ( s = inStream.readLine() ) != null && !s.equals("-*-ready-*-")) {
	if ( debug ) System.err.println("< " + s);
	if ( s.length() >= 6 && s.substring(0,6).equals("Error:") )
	  error = true;
	if ( ret.size() == 0 && s.equals(cmd) ) continue;
	ret.add(s);
      }
      if ( s == null ) lostProcess();
    } catch ( IOException ex ) {
      lostProcess();
    }

    if ( error ) {
      String mess = "";
      while ( ret.size() > 0 ) mess += (String)ret.remove(0) + "\n";
      JOptionPane.showMessageDialog(this, mess, "Error",
				    JOptionPane.ERROR_MESSAGE);
      return null;
    }
      
    return ret;
  }

  public void createHelpText() {
    inithelp =
      "<h1>T<font size=+1>HE</font>PEG</h1>This is the java-based setup " +
      "program for T<font size=-1>HE</font>PEG, the Toolkit for High Energy " +
      "Physics Event Generation. You can use it manipulate the objects in " +
      "the <i>repository</i> organized in a directory structure represented " +
      "by the tree view on the left.<p>If you double-click on an object you " +
      "will get a window representing it where you can modify its behavior " +
      "by manipulating its <i>interfaces</i>. These interfaces can be " +
      "parameters, switches or references to other objects.<p>In the end you " +
      "should end up with an <i>event generator</i> object which you can " +
      "initialize and save to a file from which it can be read into another " +
      "program and be used to generate events. The repository already " +
      "contains a number of ready-built event generators, and you are " +
      "probably better off modifying one of them, rather than trying to " +
      "build your own from scratch.";
  }

  public void setLocation(Component c) {
    for ( int i = 0; i < locations.size(); ++i )
      if ( locations.get(i) == null ) {
	locations.set(i, c);
	setLocation(c, i);
	return;
      }
    locations.add(c);
    setLocation(c, locations.size() - 1);
  }

  public void setLocation(Component c, int i) {
    Point pos = getLocationOnScreen();
    pos.translate(200 + i*50, 100 + i*50);
    c.setLocation(pos);
  }

  public void removeLocation(Component c) {
    for ( int i = 0; i < locations.size(); ++i )
      if ( locations.get(i) == c ) locations.set(i, null);
  }

  public void openObject(String name) {
    if ( name == null ) return;
    Object o = objects.get(name);
    if ( o == null ) {
      ObjectFrame obj = new ObjectFrame(this, name);
      objects.put(name, obj);
    } else {
      ObjectFrame of = (ObjectFrame)o;
      if ( !of.isVisible() ) {
	setLocation(of);
	of.setState(JFrame.NORMAL);
	of.setVisible(true);
      }
      of.toFront();
    }
  }

  public void run() {
    RunSelector r = new RunSelector(this);
    if ( r.selected() == null ) return;
    if ( r.selected().equals("Create a new run") ) {
      ObjectSelector os =
	new ObjectSelector(this, "/", "ThePEG::EventGenerator");
      if ( os.selected() == null ) return;
      openObject(os.selected().getFullName());
      return;
    }
    RunFrame run = (RunFrame)r.selected();
    setLocation(run);
    run.setVisible(true);
  }

  public void remove(RunFrame r) {
    runs.remove(r.getName());
    action("rmrun " + r.getName());
  }

  public RunFrame getRun(String name) {
    RunFrame r = (RunFrame)runs.get(name);
    if ( r == null ) {
      r = new RunFrame(this, name);
      runs.put(name, r);
      return r;
    }
    return r;
  }

  public void remove(String obj) {
    objects.remove(obj);
  }

  public int confirm(String question, String title) {
    return JOptionPane.showConfirmDialog(this, question, title,
					 JOptionPane.YES_NO_CANCEL_OPTION,
					 JOptionPane.QUESTION_MESSAGE );
  }

  public void maybeSaveActions(String question) {
    if ( actionsPerformed &&
	 confirm(question, "Save Actions?") == JOptionPane.CANCEL_OPTION )
      saveChanges();
  }


  public void safeQuit() {
    if ( actionsPerformed ) {
      switch ( confirm("Save actions before quitting?", "Really Quit?") ) {
      case JOptionPane.CANCEL_OPTION:
	return;
      case JOptionPane.YES_OPTION:
	saveChanges();
      case JOptionPane.NO_OPTION:
      }				     
    }
    System.exit(0);
  }

  public void readFile() {
    fc.setFileSelectionMode(JFileChooser.FILES_ONLY);
    fc.setDialogTitle("Read file of commands");
    fc.setDialogType(JFileChooser.OPEN_DIALOG);
    fc.setFileFilter(infilter);
    fc.setSelectedFile(new File(""));
    fc.setApproveButtonToolTipText("Read the selected file " +
				   "containing ThePEG setup commands.");
    if ( fc.showDialog(this, "Load") != JFileChooser.APPROVE_OPTION ) return;
    action("read " + fc.getSelectedFile().getAbsolutePath());
    updateTrees();
  }
    
  public void saveRepository() {
    fc.setFileSelectionMode(JFileChooser.FILES_ONLY);
    fc.setDialogTitle("Save the current repository");
    fc.setDialogType(JFileChooser.SAVE_DIALOG);
    fc.setFileFilter(rpofilter);
    fc.setSelectedFile(new File(""));
    fc.setApproveButtonToolTipText("Save the current repository to " +
				   "the selected file");
    if ( fc.showDialog(this, "Save") != JFileChooser.APPROVE_OPTION ) return;
    exec("save " + fc.getSelectedFile().getAbsolutePath());
  }
    
  public void loadRepository() {
    maybeSaveActions("Save actions before loading new repository?");
    fc.setFileSelectionMode(JFileChooser.FILES_ONLY);
    fc.setDialogTitle("Load a new repository");
    fc.setDialogType(JFileChooser.OPEN_DIALOG);
    fc.setFileFilter(rpofilter);
    fc.setSelectedFile(new File(""));
    fc.setApproveButtonToolTipText("Load a new repository from the " +
				   "selected file.");
    if ( fc.showDialog(this, "Load") != JFileChooser.APPROVE_OPTION ) return;
    exec("load " + fc.getSelectedFile().getAbsolutePath());
    updateTrees();
    actions.clear();
    actionsPerformed = false;
    saveActions.setEnabled(false);
  }

  public void saveChanges() {
    if ( actionFile == null ) saveChangesAs();
    else {
      try {
	PrintWriter out = new PrintWriter(new FileWriter(actionFile));
	for ( Iterator it = actions.iterator(); it.hasNext(); ) {
	  String cmd = (String)it.next();
	  out.println(cmd);
	}
	out.close();
	actionsPerformed = false;
      }
      catch ( IOException ex ) {
	JOptionPane.showMessageDialog(this, "Save failed!", ex.getMessage(),
				      JOptionPane.ERROR_MESSAGE);
      }
    }
  }

  public void loadLibrary() {
    fc.setFileSelectionMode(JFileChooser.FILES_ONLY);
    fc.setDialogTitle("Load Dynamic Library");
    fc.setDialogType(JFileChooser.OPEN_DIALOG);
    fc.setSelectedFile(new File(""));
    fc.setFileFilter(sofilter);
    fc.setApproveButtonToolTipText(ttt + "Load the selected library. " +
				   "Note that the directory must " +
				   "be in the search path.");
    if ( fc.showDialog(this, "Load") != JFileChooser.APPROVE_OPTION ) return;
    action("library " + fc.getSelectedFile().getName());
  }

  public void addLibraryPath(boolean prepend) {
    fc.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
    fc.setDialogType(JFileChooser.OPEN_DIALOG);
    String action = "Append";
    if ( prepend ) action = "Prepend";
    fc.setDialogTitle(action + " directory to load path");
    fc.setDialogType(JFileChooser.OPEN_DIALOG);
    fc.setApproveButtonToolTipText(ttt + action +
				   " the delected directory to the list" +
				   "of search paths for loading dynaic " +
				   "libraries.");
    fc.setSelectedFile(new File(""));
    if ( fc.showDialog(this, action) != JFileChooser.APPROVE_OPTION ) return;
    if ( prepend )
      action("prependpath " + fc.getSelectedFile().getAbsolutePath());
    else
      action("appendpath " + fc.getSelectedFile().getAbsolutePath());
  }    

  public void saveChangesAs() {
    fc.setFileSelectionMode(JFileChooser.FILES_ONLY);
    fc.setDialogTitle("Save actions to file");
    fc.setDialogType(JFileChooser.SAVE_DIALOG);
    fc.setFileFilter(infilter);
    fc.setSelectedFile(new File(""));
    fc.setApproveButtonToolTipText("Save the actions performed on the " +
				   "repository to the selected file.");
    if ( fc.showSaveDialog(this) != JFileChooser.APPROVE_OPTION ) return;
    actionFile = fc.getSelectedFile();
    if ( actionFile != null ) saveChanges();
  }

  public void showOnly() {
    ClassSelector sel = new ClassSelector(this);
    String cl = sel.selected();
    if ( cl == null ) return;
    Enumeration it = showgroup.getElements();
    while ( it.hasMoreElements() ) {
      JRadioButtonMenuItem b = (JRadioButtonMenuItem)it.nextElement();
      if ( b.getText().equals(cl) ) {
	showgroup.remove(b);
	showmenu.remove(b);
      }
    }
    JRadioButtonMenuItem b = new JRadioButtonMenuItem(cl, true);
    b.addActionListener(this);
    showgroup.add(b);
    showmenu.add(b);
    b.doClick();
  }

  public void showOnly(String cl) {
    if ( cl.equals("All classes") )
      tree.setClass("");
    else if ( cl.equals("Event generators") )
      tree.setClass("ThePEG::EventGenerator");
    else 
      tree.setClass(cl);
    tree.update();
  }

  public void actionPerformed(ActionEvent e) {
    if ( e.getSource() == quit ) safeQuit();
    else if ( e.getSource() == saveActions ) saveChanges();
    else if ( e.getSource() == saveActionsAs ) saveChangesAs();
    else if ( e.getSource() == running ) run();
    else if ( e.getSource() == readActions ) readFile();
    else if ( e.getSource() == loadlib ) loadLibrary();
    else if ( e.getSource() == loadrepo ) loadRepository();
    else if ( e.getSource() == saverepo ) saveRepository();
    else if ( e.getSource() == showonly ) showOnly();
    else if ( e.getSource() == prependpath ||
	      e.getSource() == appendpath )
      addLibraryPath(e.getSource() == prependpath);
    else if ( e.getSource() instanceof JRadioButtonMenuItem ) {
      JRadioButtonMenuItem b = (JRadioButtonMenuItem)e.getSource();
      showOnly(b.getText());
    }
  }

  public ObjectRef getClip() {
    return clip;
  }

  public ObjectRef getDragged() {
    return dragged;
  }

  public void setDragged(ObjectRef o) {
    dragged = o;
  }

  public boolean isClipClone() {
    return clipIsClone;
  }

  public void copy(ObjectRef obj) {
    clip = obj;
    clipIsClone = false;
  }

  public void clone(ObjectRef obj) {
    clip = obj;
    clipIsClone = true;
  }

  public void delete(ObjectRef obj) {
    if ( obj.isDir() )
      action("rmdir " + obj.getFullName());
    else
      action("rm " + obj.getFullName());
    updateTrees();
  }

  public void rename(ObjectRef obj, String newName) {
    if ( obj.isDir() ) return;
    action("mv " + obj.getFullName() + " " + newName);
    updateTrees();
  }

  public void addTree(BrowserTree t) {
    allTrees.add(t);
  }

  public void removeTree(BrowserTree t) {
    allTrees.remove(t);
  }

  public void updateTrees() {
    for ( Iterator it = allTrees.iterator(); it.hasNext(); )
      ((BrowserTree)it.next()).update();
  }

  public static void main(String [] args) {
    for ( int iarg = 0; iarg < args.length; ++iarg )
      if ( args[iarg].equals("--classcheck") ) {
	BrowserTree.classcheck();
	ClassSelector.classcheck();
	Command.classcheck();
	FullSlider.classcheck();
	Interface.classcheck();
	ObjectFrame.classcheck();
	ObjectNode.classcheck();
	ObjectRef.classcheck();
	ObjectSelector.classcheck();
	ParVector.classcheck();
	Parameter.classcheck();
	PopCardPanel.classcheck();
	RefRef.classcheck();
	RefVector.classcheck();
	Reference.classcheck();
	RunSelector.classcheck();
	RunFrame.classcheck();
	StringParameter.classcheck();
	Switch.classcheck();
	SwitchOption.classcheck();
	return;
      }
    SetupThePEG setup = new SetupThePEG(args);
  }

}
