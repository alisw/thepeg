package ThePEG;

import javax.swing.*;
import java.util.*;
import java.awt.*;
import java.awt.event.*;

public class PopCardPanel extends JPanel implements ActionListener {

  CardLayout cards = new CardLayout();
  JComboBox selector = new JComboBox();
  JPanel panel = new JPanel();

  public PopCardPanel() {
    setLayout(new BorderLayout());
    panel.setLayout(cards);
    selector.addActionListener(this);
    add(panel, BorderLayout.CENTER);
    add(selector, BorderLayout.NORTH);
  }

  public int getItemCount() {
    return selector.getItemCount();
  }

  public void addCard(String name, JComponent obj) {
    panel.add(obj, name);
    selector.addItem(name);
  }

  public void actionPerformed(ActionEvent e) {
    if ( e.getSource() == selector )
      cards.show(panel, (String)selector.getSelectedItem());
  }

  public static void classcheck() {}

}
