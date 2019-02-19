#include <QApplication>

#include "main_dialog.h"

int main(int argc, char* argv[]){

  QApplication app(argc, argv);
  MainDialog dial;
  dial.show();
  return app.exec();

  return 0;
}
