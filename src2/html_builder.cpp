#include "html_builder.h"

#include <iostream>

using namespace std;
using namespace farn;

HtmlBuilder::HtmlBuilder(string flnm) : dest(flnm.c_str()) {
  if (!dest) {
    cerr << "can't open " << flnm << " in htmlBuilder" << endl;
    return;
  }

  dest << "<html>" << endl;
  dest << "<head>" << endl;
  dest << "  <title>Different frequencies</title>" << endl;
  dest << "</head>" << endl;
  dest << "<body>" << endl;
  dest << "  <table>" << endl;
}

HtmlBuilder::~HtmlBuilder() {
  dest << "  </table>" << endl;
  dest << "</body>" << endl;
  dest << "</html>" << endl;
}

void
HtmlBuilder::beginRow() {
  dest << "<tr>" << endl;
}

void
HtmlBuilder::endRow() {
  dest << "</tr>" << endl;
}

void
HtmlBuilder::setPicture(string name) {
  dest << "  <td><img src=\"" << name << "\"></td>" << endl;
}

void
HtmlBuilder::setHeaderItem(string itm) {
  dest << "  <td align=\"middle\">" << itm << "</td>" << endl;
}
