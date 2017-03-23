package ThePEG;

import javax.swing.JEditorPane;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSplitPane;
import javax.swing.UIManager;

import javax.swing.*;
import javax.swing.tree.*;
import javax.swing.event.*;

import java.net.URL;
import java.io.IOException;
import java.awt.*;
import java.awt.event.*;
import java.awt.dnd.*;
import java.awt.datatransfer.*;
import java.util.*;

public class BrowserTree extends JPanel
                      implements TreeSelectionListener, TreeWillExpandListener,
				 ActionListener, MouseListener,
				 TreeModelListener {

  private SetupThePEG thepeg;
  private JTree tree;
  private DefaultTreeModel treeModel;
  private ObjectNode top;
  private static boolean DEBUG = false;
  JPopupMenu edit = new JPopupMenu("Edit");
  JMenuItem open;
  JMenuItem copy;
  JMenuItem clone;
  JMenuItem remove;
  JMenuItem create;
  JMenuItem newdir;
  JMenuItem paste;
  JMenu exedit = new JMenu("Edit");
  JMenuItem exopen;
  JMenuItem excopy;
  JMenuItem exclone;
  JMenuItem exremove;
  JMenuItem excreate;
  JMenuItem exnewdir;
  JMenuItem expaste;
  DropTarget dropper;
  String pclass;
  boolean modal = false;

  public BrowserTree(SetupThePEG owner) {
    this(owner, "");
  }

  public BrowserTree(SetupThePEG owner, String pc) {
    this(owner, pc, false);
  }

  public BrowserTree(SetupThePEG owner, String pc, boolean mod, String dir) {
    this(owner, pc, mod);
    expandToDir(top, dir);
  }

  public BrowserTree(SetupThePEG owner, String pc, boolean mod) {
    super(new GridLayout(1,0));

    modal = mod;
    pclass = pc;
    thepeg = owner;
    thepeg.addTree(this);

    // Create the menus.
    open = createMenuItem("Open", edit);
    edit.addSeparator();
    create = createMenuItem("New object", edit);
    newdir = createMenuItem("New directory", edit);
    remove = createMenuItem("Delete", edit);
    edit.addSeparator();
    copy = createMenuItem("Copy", edit);
    clone = createMenuItem("Clone", edit);
    paste = createMenuItem("Paste", edit);

    exopen = createMenuItem("Open", exedit);
    exedit.addSeparator();
    excreate = createMenuItem("New object", exedit);
    exnewdir = createMenuItem("New directory", exedit);
    exremove = createMenuItem("Delete", exedit);
    exedit.addSeparator();
    excopy = createMenuItem("Copy", exedit);
    exclone = createMenuItem("Clone", exedit);
    expaste = createMenuItem("Paste", exedit);

    // Create the nodes.
    top = new ObjectNode(thepeg, "<Root>", "/", true);
    top.setExpanded(true);

    // Create a tree that allows one selection at a time.
    tree = new JTree(top) {
	public String getToolTipText(MouseEvent e) {
	  TreePath selPath =
	    tree.getPathForLocation(e.getX(), e.getY());
	  if ( selPath == null ) return "";
	  ObjectRef o = getLastObject(selPath);
	  if ( o.isDir() ) return "Directory";
	  else if ( o.getFullName().equals("") ) return "Null object";
	  else return "Object of class " + o.getClassName();
	}
      };
    ToolTipManager.sharedInstance().registerComponent(tree);
    tree.getSelectionModel().setSelectionMode
      (TreeSelectionModel.SINGLE_TREE_SELECTION);
    tree.setEditable(true);
    tree.setDragEnabled(true);

    tree.getModel().addTreeModelListener(this);

    // Listen for when the selection changes.
    tree.addTreeSelectionListener(this);
    tree.addTreeWillExpandListener(this);
    tree.addMouseListener(this);
    treeModel = (DefaultTreeModel)tree.getModel();
    updateNode(top);

    dropper = new DropTarget(tree, new DropTargetAdapter() {
	public void drop(DropTargetDropEvent e) {
	  TreePath selPath =
	    tree.getPathForLocation(e.getLocation().x,
				    e.getLocation().y);
	  ObjectRef dragged = thepeg.getDragged();
	  thepeg.setDragged(null);
	  if ( selPath == null || dragged == null ) {
	    e.rejectDrop();
	    return;
	  }
	  ObjectRef o = getLastObject(selPath);
	  if ( !o.isDir() || dragged.isDir() ) {
	    e.rejectDrop();
	    return;
	  }
	  tree.expandPath(selPath);
	  thepeg.rename(dragged, o.getFullName() + dragged.getName());
	}
      });

    //Create the scroll pane and add the tree to it. 
    JScrollPane treeView = new JScrollPane(tree);

    add(treeView);
    tree.expandRow(0);

  }

  public JMenuItem createMenuItem(String name, JMenu menu) {
    JMenuItem it = new JMenuItem(name);
    it.addActionListener(this);
    menu.add(it);
    it.setEnabled(false);
    return it;
  }

  public JMenuItem createMenuItem(String name, JPopupMenu menu) {
    JMenuItem it = new JMenuItem(name);
    it.addActionListener(this);
    menu.add(it);
    it.setEnabled(false);
    return it;
  }

  public JMenu getExternalMenu() {
    return exedit;
  }

  public void dispose() {
    thepeg.removeTree(this);
  }

  public void addMouseListener(MouseListener l) {
    tree.addMouseListener(l);
  }

  public void valueChanged(TreeSelectionEvent e) {
    updateNode((ObjectNode)tree.getLastSelectedPathComponent());
  }

  public ObjectRef getObjectRef(ObjectNode node) {
    return node.getObject();
  }
    
  public void update() {
    updateNode(top);
  }

  public void addNode(ObjectNode parent, ObjectNode child) {
    insertNode(parent, child, parent.getChildCount());
  }

  public void insertNode(ObjectNode parent, ObjectNode child, int indx) {
    treeModel.insertNodeInto(child, parent, indx);
  }

  public void removeNode(ObjectNode child) {
    treeModel.removeNodeFromParent(child);
  }

  public void removeChildren(ObjectNode parent) {
    while ( parent.getChildCount() > 0 )
      removeNode((ObjectNode)parent.getFirstChild());
  }

  public void updateNode(ObjectNode node) {
    if ( node == null || !node.isDir() ) return;

    if ( !node.isExpanded() ) {
      if ( node.getChildCount() == 0 )
	addNode(node, new ObjectNode(thepeg, "<Empty>", "", false));
      return;
    }

    if ( node.getChildCount() == 1 &&
	 ((ObjectNode)node.getFirstChild()).getFullName().equals("") )
      removeNode((ObjectNode)node.getFirstChild());
      
    LinkedList ret = thepeg.exec("ls " + node.getFullName() + " " + pclass);
    if ( ret.size() == 0 ) {
      removeChildren(node);
      addNode(node, new ObjectNode(thepeg, "<Empty>", "", false));
      return;
    }
    
    int nobj = ret.size();
    HashSet names = new HashSet(ret);
    for ( int i = node.getChildCount() - 1; i >= 0; --i ) {
      ObjectNode n = (ObjectNode)node.getChildAt(i);
      if ( !names.contains(n.getFullName()) ) removeNode(n);
    }
    for ( int i = 0; i < nobj; ++i ) {
      String name = (String)ret.remove(0);
      if ( i >= node.getChildCount() ||
	   !((ObjectNode)node.getChildAt(i)).getFullName().equals(name) )
	insertNode(node, new ObjectNode(thepeg, name), i);
    }

    for ( int i = node.getChildCount() - 1; i >= 0; --i )
      updateNode((ObjectNode)node.getChildAt(i));

  }

  public void setClass(String cl) {
    pclass = cl;
  }

  public boolean expandToDir(ObjectNode node, String name) {
    int l = node.getFullName().length();
    if ( l > name.length() ) return false;
    if ( !name.substring(0, l).equals(node.getFullName()) ) return false;
    if ( !node.isDir() ) return true;
    expandNode(node);
    for ( int i = node.getChildCount() - 1; i >= 0; --i )
      if ( expandToDir((ObjectNode)node.getChildAt(i), name) ) return true;
    return false;
  }

  public ObjectNode getLastNode(TreePath path) {
    if ( path == null ) return null;
    return (ObjectNode)(path.getLastPathComponent());
  }

  public ObjectRef getLastObject(TreePath path) {
    return getLastNode(path).getObject();
  }

  public ObjectRef getSelectedObject() {
    ObjectNode node = getSelectedNode();
    if ( node == null ) return null;
    return node.getObject();
  }

  public ObjectNode getSelectedNode() {
    return getLastNode(tree.getSelectionPath());
  }

  public void openObject(TreePath path) {
    getLastObject(path).open();
  }

  public TreePath getNodePath(ObjectNode node) {
    return new TreePath(node.getPath());
  }

  public void expandNode(ObjectNode node) {
    tree.expandPath(getNodePath(node));
  }

  public void collapseNode(ObjectNode node) {
    tree.collapsePath(getNodePath(node));
  }

  public void checkSelection(ObjectRef node) {
    if ( node.isDir() ) {
      newdir.setEnabled(true);
      open.setEnabled(!modal);
      copy.setEnabled(false);
      remove.setEnabled(true);
      clone.setEnabled(false);
      create.setEnabled(true);
      paste.setEnabled(thepeg.getClip() != null);
      exnewdir.setEnabled(true);
      exopen.setEnabled(true);
      excopy.setEnabled(false);
      exremove.setEnabled(true);
      exclone.setEnabled(false);
      excreate.setEnabled(true);
      expaste.setEnabled(thepeg.getClip() != null);
    } else {
      newdir.setEnabled(false);
      open.setEnabled(!modal);
      copy.setEnabled(true);
      remove.setEnabled(true);
      clone.setEnabled(true);
      create.setEnabled(false);
      paste.setEnabled(false);
      exnewdir.setEnabled(false);
      exopen.setEnabled(true);
      excopy.setEnabled(true);
      exremove.setEnabled(true);
      exclone.setEnabled(true);
      excreate.setEnabled(false);
      expaste.setEnabled(false);
    }
  }

  public ObjectNode getVisibleNode(ObjectNode node, String name) {
    for ( int i = 0; i < node.getChildCount(); ++i ) {
      ObjectNode n = (ObjectNode)node.getChildAt(i);
      if ( name.equals(n.getFullName()) ) return n;
      if ( n.isDir() && n.isExpanded() ) {
	ObjectNode ret = getVisibleNode(n, name);
	if ( ret != null ) return ret;
      }
    }
    return null;
  }

  public void treeWillCollapse(TreeExpansionEvent event) {
    getLastNode(event.getPath()).setExpanded(false);
  }

  public void treeWillExpand(TreeExpansionEvent event) {
    ObjectNode node = getLastNode(event.getPath());
    node.setExpanded(true);
    updateNode(node);
  }

  public void actionPerformed(ActionEvent e) {
    if ( e.getSource() == open || e.getSource() == exopen ) {
      ObjectRef ref = getSelectedObject();
      if ( ref.isDir() ) tree.expandPath(tree.getSelectionPath());
      else ref.open();
    }
    else if ( e.getSource() == copy || e.getSource() == excopy ) {
      thepeg.copy(getSelectedObject());
    }
    else if ( e.getSource() == remove || e.getSource() == exremove ) {
      thepeg.delete(getSelectedObject());
    }
    else if ( e.getSource() == clone || e.getSource() == exclone ) {
      thepeg.clone(getSelectedObject());
    }
    else if ( e.getSource() == create || e.getSource() == excreate ) {
      ObjectNode dir = getSelectedNode();
      if ( !dir.isDir() ) return;
      ClassSelector sel = new ClassSelector(thepeg);
      String cl = sel.selected();
      if ( cl == null ) return;
      String newName = dir.getFullName() + "NewObject";
      thepeg.action("create " + cl + " " + newName);
      expandNode(dir);
      update();
      ObjectNode n = getVisibleNode(dir, newName);
      if ( n != null ) tree.startEditingAtPath(getNodePath(n));
    }
    else if ( e.getSource() == newdir || e.getSource() == exnewdir ) {
      ObjectNode dir = getSelectedNode();
      String newName = dir.getFullName() + "NewDir/";
      thepeg.action("mkdir " + newName);
      expandNode(dir);
      update();
      ObjectNode n = getVisibleNode(dir, newName);
      if ( n != null ) tree.startEditingAtPath(getNodePath(n));
    }
    else if ( e.getSource() == paste || e.getSource() == expaste ) {
      ObjectNode dir = getSelectedNode();
      if ( !dir.isDir() ) return;
      if ( thepeg.isClipClone() ) {
	String newName = dir.getFullName();
	thepeg.action("rcp " + thepeg.getClip().getFullName() + " " + newName);
      } else {
	String newName = dir.getFullName() + thepeg.getClip().getName();
	thepeg.action("cp " + thepeg.getClip().getFullName() + " " + newName);
      }
      expandNode(dir);
      thepeg.updateTrees();
    }
	
  }

  public void mousePressed(MouseEvent e) {
    TreePath selPath = tree.getPathForLocation(e.getX(), e.getY());
    if( selPath != null ) {
      tree.setSelectionPath(selPath);
      checkSelection(getLastObject(selPath));
      if ( e.isPopupTrigger() ) {
	edit.show(e.getComponent(), e.getX(), e.getY());
      }
      else if(e.getClickCount() == 2) {
	if ( !modal ) openObject(selPath);
      }
      else {
	ObjectRef dragged = getLastObject(selPath);
	if ( dragged.isDir() ) dragged = null;
	thepeg.setDragged(dragged);
      }
    }
  }

  public void mouseClicked(MouseEvent e) {}

  public void mouseEntered(MouseEvent e) {}

  public void mouseExited(MouseEvent e) {}

  public void mouseReleased(MouseEvent e) {}

  public void treeNodesChanged(TreeModelEvent e) {
    ObjectNode node = getSelectedNode();
    if ( node.isRoot() ) return;
    if ( node.toString().equals(node.getObject().getName()) ) return;
    if ( node.isDir() ) {
      if ( node.getChildCount() > 1 ) node.setUserObject(node.getObject());
      else if ( node.getChildCount() == 1 &&
		!((ObjectNode)node.getFirstChild()).getFullName().equals("") )
	  node.setUserObject(node.getObject());
      else {
	thepeg.action("rmdir " + node.getFullName());
	String newName = node.getFullName();
	newName = newName.substring(0, newName.length() - 1);
	newName = newName.substring(0, newName.lastIndexOf('/') + 1);
	newName += node.toString();
	thepeg.action("mkdir " + newName);
	node.setUserObject(new ObjectRef(thepeg, newName));
	thepeg.updateTrees();
      }
      return;
    }
    ObjectNode par = (ObjectNode)node.getParent();
    thepeg.rename(node.getObject(),
		  par.getObject().getFullName() + node.toString());
  }

  public void treeNodesInserted(TreeModelEvent e) {}

  public void treeNodesRemoved(TreeModelEvent e) {}

  public void treeStructureChanged(TreeModelEvent e) {}

  public static void classcheck() {}

}
