#pragma once

#include <fstream>

namespace farn {
  //simple builder of table with images
  class HtmlBuilder {
    std::ofstream dest;
  public:
    HtmlBuilder(std::string);
    ~HtmlBuilder();

    void beginRow();
    void endRow();
    void setPicture(std::string name);
    void setHeaderItem(std::string itm);
  };
}
