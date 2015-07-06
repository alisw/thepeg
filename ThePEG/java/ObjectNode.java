package ThePEG;

import javax.swing.tree.*;

public class ObjectNode extends DefaultMutableTreeNode {

  ObjectRef ref;
  boolean expanded = false;

  public ObjectNode(SetupThePEG tp, String nm, String fn, boolean isd) {
    super(new ObjectRef(tp, nm, fn, isd));
    ref = (ObjectRef)getUserObject();
    setAllowsChildren(ref.isDir());
  }

  public ObjectNode(SetupThePEG tp, String f) {
    super(new ObjectRef(tp, f));
    ref = (ObjectRef)getUserObject();
    setAllowsChildren(ref.isDir());
  }

  public ObjectRef getObject() {
    return ref;
  }

  public boolean isDir() {
    return ref.isDir();
  }

  public String getFullName() {
    return ref.getFullName();
  }

  public void setExpanded(boolean ex) {
    expanded = ex;
  }

  public boolean isExpanded() {
    return expanded;
  }

  public static void classcheck() {}

}

