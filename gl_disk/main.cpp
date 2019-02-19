#include <QApplication>

#include "main_dialog.h"

int main(int argc, char * argv[]){
  QApplication app(argc, argv);
  MainDialog md;
  md.show();
  return app.exec();
}
