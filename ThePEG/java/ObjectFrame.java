package ThePEG;

import javax.swing.*;
import javax.swing.event.*;
import java.util.*;
import java.awt.*;
import java.awt.event.*;

public class ObjectFrame extends JFrame
                         implements ActionListener, MouseListener,
				    HyperlinkListener {

  SetupThePEG thepeg;
  String object = "";
  String name = "";
  String pclass = "";
  JList parameterList;
  JList switchList;
  JList referenceList;
  JList commandList;
  JList primaryList;
  JList secondaryList;
  int ctype = 0;
  static final int EVENT_GENERATOR = 1;
  PopCardPanel interfaces = new PopCardPanel();
  JButton close = new JButton("Close");
  JButton makerun = new JButton("Make run");
  JButton saverun = new JButton("Save run");
  JEditorPane info = new JEditorPane();

  HashMap imap = new HashMap();

  public ObjectFrame(SetupThePEG frame, String obj) {
    setIconImage(frame.icon);
    name = obj.substring(obj.lastIndexOf("/") + 1);
    object = obj;
    thepeg = frame;
    setup();
    setTitle(name + " [" + pclass + "]");
    if ( primaryList != null ) {
      primaryList.addMouseListener(this);
      addCard("Primary interfaces:", primaryList);
    }
    if ( parameterList != null ) {
      parameterList.addMouseListener(this);
      addCard("All Parameters:", parameterList);
    }
    if ( switchList != null ) {
      switchList.addMouseListener(this);
      addCard("All Switches:", switchList);
    }
    if ( referenceList != null ) {
      referenceList.addMouseListener(this);
      addCard("All References:", referenceList);
    }
    if ( commandList != null ) {
      commandList.addMouseListener(this);
      addCard("All Commands:", commandList);
    }
    if ( secondaryList != null ) {
      secondaryList.addMouseListener(this);
      addCard("Secondary interfaces:", secondaryList);
    }

    JSplitPane split = new JSplitPane(JSplitPane.VERTICAL_SPLIT);
    JScrollPane s = new JScrollPane(info);
    s.setMinimumSize(new Dimension(300,250));
      
    split.setLeftComponent(s);

    if ( interfaces.getItemCount() > 0 ) {
      split.setRightComponent(interfaces);
    } else {
      JPanel noinf = new JPanel();
      noinf.add(new JLabel("There are no interfaces defined for " + pclass));
      split.setRightComponent(noinf);
    }
    getContentPane().setLayout(new BorderLayout());
    getContentPane().add(split, BorderLayout.CENTER);

    close.addActionListener(this);
    JPanel buttons = new JPanel();
    buttons.add(close);

    switch ( ctype ) {
    case EVENT_GENERATOR:
      makerun.addActionListener(this);
      buttons.add(makerun);
      saverun.addActionListener(this);
      buttons.add(saverun);
      break;
    }

    getContentPane().add(buttons, BorderLayout.SOUTH);

    thepeg.setLocation(this);

    setSize(400, 500);

    setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
    setVisible(true);
    
  }

  public String getFullName() {
    return object;
  }

  public String getName() {
    return name;
  }

  private void addCard(String name, JList list) {
    if ( list == null ) return;
    interfaces.addCard(name, new JScrollPane(list));
  }

  public static String htmlFormat(String text) {
    return text.replaceAll("<interface>(\\w+::)?(\\w+)</interface>",
			   "<a href=\"interface:$2\">$2</a>");
  }

  private void fixDescription(String desc) {
    if ( pclass.equals("ThePEG::EventGenerator") ) ctype = EVENT_GENERATOR;
    String text = "<h2>Object \"" + getName() + "\" of class <tt>"
      + pclass + "</tt></h2>";
    if ( desc.length() < 25 ||
	 !desc.substring(0, 25).equals("There is no documentation") )
      text += htmlFormat(desc);
    LinkedList ret = thepeg.exec("baseclasses " + pclass);
    if ( ret.size() > 0 ) ret.remove(0);
    while ( ret.size() > 0 ) {
      String bclass = (String)ret.remove(0);
      if ( bclass.equals("ThePEG::HandlerBase") ) break;
      text += "<h3>Inherits from class <tt>" + bclass + "</tt></h3>";
      LinkedList bret = thepeg.exec("describeclass " + bclass);
      String inf = "";
      while ( bret.size() > 0 ) inf += (String)bret.remove(0);
      if ( inf.length() < 25 ||
	   !inf.substring(0, 25).equals("There is no documentation") )
      text += htmlFormat(inf);
    }
    text += "<br>You can modify the behavior of this object by " +
      "manipulating the interfaces below.";
    info.setContentType("text/html");
    info.setText("<font size=-3>" + text + "</font>");
    info.setEditable(false);
    info.setCaretPosition(0);
    info.addHyperlinkListener(this);
  }

  private void setup() {
    LinkedList ret = thepeg.exec("fulldescribe " + object);
    if ( ret.size() > 0 ) ret.remove(0);
    if ( ret.size() > 0 ) pclass = (String)ret.remove(0);
    String line = "";
    String description = "";
    while ( ret.size() > 0 &&
	    !(line = (String)ret.remove(0)).equals("Interfaces:") )
      description += line + " ";
    Vector parameters = new Vector();
    Vector switches = new Vector();
    Vector references = new Vector();
    Vector commands = new Vector();
    Vector secondary = new Vector();
    Vector primary = new Vector();
    boolean prim = true;
    while ( ret.size() > 0 ) {
      line = (String)ret.remove(0);
      if ( line.length() == 0 ) continue;
      String name = "";
      if ( line.length() > 3 ) name = line.substring(3);
      switch ( line.charAt(0) ) {
      case '0':
	prim = false;
	break;
      case 'P':
	parameters.add(name);
	break;
      case 'S':
	switches.add(name);
	break;
      case 'C':
	commands.add(name);
	break;
      case 'D':
	continue;
      case 'R':
	name = line.substring(line.lastIndexOf('>') + 2);
	references.add(name);
	break;
      case 'V':
	if ( line.charAt(1) == '<' ) {
	  name = line.substring(line.lastIndexOf('>') + 2);
	  references.add(name);
	} else
	  parameters.add(name);
      }
      if ( prim ) primary.add(name);
      else secondary.add(name);
    }
    primaryList = new JList(primary);
    if ( parameters.size() > 0 ) parameterList = new JList(parameters);
    if ( switches.size() > 0 ) switchList = new JList(switches);
    if ( references.size() > 0 ) referenceList = new JList(references);
    if ( commands.size() > 0 ) commandList = new JList(commands);
    if ( secondary.size() > 0 ) secondaryList = new JList(secondary);
    fixDescription(description);
  }

  public void dispose() {
    super.dispose();
    thepeg.removeLocation(this);
  }

  private void showInterface(String name) {
    if ( imap.get(name) != null ) {
      Interface i = (Interface)imap.get(name);
      if ( !i.isVisible() ) {
	thepeg.setLocation(i);
	i.setVisible(true);
      }
      i.toFront();
      return;
    }

    if ( name == null ) return;
    LinkedList ret = thepeg.exec("fulldescribe " + object + ":" + name);
    if ( ret.size() == 0 ) return;
    String par = (String)ret.getFirst();
    Interface ifc = null;
    switch ( par.charAt(0) ) {
    case 'P':
      if ( par.charAt(1) == 'i' )
	ifc = new Parameter(thepeg, this, ret, true);
      else if ( par.charAt(1) == 's' )
	ifc = new StringParameter(thepeg, this, ret, StringParameter.NOFILE);
      else if ( par.charAt(1) == 'D' )
	ifc = new StringParameter(thepeg, this, ret, StringParameter.DIRECTORY);
      else if ( par.charAt(1) == 'F' )
	ifc = new StringParameter(thepeg, this, ret, StringParameter.FILE);
      else
	ifc = new Parameter(thepeg, this, ret, false);
      break;
    case 'S':
      ifc = new Switch(thepeg, this, ret);
      break;
    case 'C':
      ifc = new Command(thepeg, this, ret);
      break;
    case 'R':
      ifc = new Reference(thepeg, this, ret);
      break;
    case 'V':
      if ( par.charAt(1) == '<' )
	ifc = new RefVector(thepeg, this, ret);
      else if ( par.charAt(1) == 'i' )
	ifc = new ParVector(thepeg, this, ret, true);
      else if ( par.charAt(1) != 's' )
	ifc = new ParVector(thepeg, this, ret, false);
      break;
    }
    imap.put(name, ifc);
  }

  public void actionPerformed(ActionEvent e) {
    if ( e.getSource() == close )
      dispose();
    else if ( e.getSource() == makerun )
      showInterface("MakeRun");
    else if ( e.getSource() == saverun )
      showInterface("SaveRun");
  }

  public void mouseClicked(MouseEvent e) {
    if (e.getClickCount() >= 2 ) {
      if ( e.getSource() == parameterList ) {
	showInterface((String)parameterList.getSelectedValue());
      }
      else if ( e.getSource() == switchList ) {
	showInterface((String)switchList.getSelectedValue());
      }
      else if ( e.getSource() == referenceList ) {
	showInterface((String)referenceList.getSelectedValue());
      }
      else if ( e.getSource() == commandList ) {
	showInterface((String)commandList.getSelectedValue());
      }
      else if ( e.getSource() == primaryList ) {
	showInterface((String)primaryList.getSelectedValue());
      }
      else if ( e.getSource() == secondaryList ) {
	showInterface((String)secondaryList.getSelectedValue());
      }
    }
  }

  public void mouseEntered(MouseEvent e) {}
  public void mouseExited(MouseEvent e) {}
  public void mousePressed(MouseEvent e) {}
  public void mouseReleased(MouseEvent e) {}

  public void hyperlinkUpdate(HyperlinkEvent e) {
    if (e.getEventType() == HyperlinkEvent.EventType.ACTIVATED) {
      String target = e.getDescription();
      String type = target.substring(0, target.indexOf(':'));
      target = target.substring(target.indexOf(':') + 1, target.length());
      showInterface(target);
    }
  }

  public static void classcheck() {}

}

