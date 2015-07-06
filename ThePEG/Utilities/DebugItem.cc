// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the DebugItem class.
//

#include "DebugItem.h"
#include "ThePEG/Utilities/Debug.h"

using namespace ThePEG;

DebugItem::DebugItem(string itemname, int level): debug(false) {
  if ( level <= Debug::level ) debug = true;
  items().insert(make_pair(itemname, this));
  map<string,long>::iterator it = nametics().find(itemname);
  if ( it != nametics().end() ) {
    if ( ticker() >= it->second ) debug = true;
    else itemtics().insert(make_pair(it->second, this));
  } else {
    while ( itemname.rfind("::") != string::npos ) {
      itemname = itemname.substr(0, itemname.rfind("::")) + "::all";
      it = nametics().find(itemname);
      if ( it != nametics().end() ) {
	if ( ticker() >= it->second ) debug = true;
	else itemtics().insert(make_pair(it->second, this));
      }
      itemname = itemname.substr(0, itemname.rfind("::"));
    }
  }
}

void DebugItem::tic() {
  ticker()++;
  multimap<long,DebugItem*>::iterator it = itemtics().begin();
  while ( it != itemtics().end() &&
	  ticker() >= it->first ) (it++)->second->debug = true;
  itemtics().erase(itemtics().begin(), it);
}

void DebugItem::setDebugItem(string itemname, long after) {
  typedef multimap<string,DebugItem*>::iterator ItemIt;
  if ( itemname.rfind('=') != string::npos ) {
    after =  atoi(itemname.substr(itemname.rfind('=') + 1).c_str());
    itemname = itemname.substr(0, itemname.rfind('='));
  }
  nametics()[itemname] = after;
  pair<ItemIt,ItemIt> range = items().equal_range(itemname);
  while ( range.first != range.second ) {
    if ( ticker() >= after ) (range.first++)->second->debug = true;
    else itemtics().insert(make_pair(after, (range.first++)->second));
  }
  if ( itemname.substr(itemname.length() - 5) == "::all" ) {
    itemname = itemname.substr(itemname.length() - 3);
    for ( ItemIt it = items().begin(); it != items().end(); ++it )
      if ( it->first.substr(0, itemname.length()) == itemname ) {
	if ( ticker() >= after ) it->second->debug = true;
	else itemtics().insert(make_pair(after, it->second));
      }	
  }
}

long & DebugItem::ticker() {
  static long tics = 0;
  return tics;
}

multimap<string,DebugItem*> & DebugItem::items() {
  static multimap<string,DebugItem*> itemmap;
  return itemmap;
}

multimap<long,DebugItem*> & DebugItem::itemtics() {
  static multimap<long,DebugItem*> itemmap;
  return itemmap;
}

map<string,long> & DebugItem::nametics() {
  static map<string,long> namemap;
  return namemap;
}
