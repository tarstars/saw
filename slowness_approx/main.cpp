#include <QApplication>
#include "main_dialog.h"

int main(int argc, char* argv[]){
  QApplication app(argc, argv);
  MainDialog sam;
  sam.show();
  return app.exec();
}
