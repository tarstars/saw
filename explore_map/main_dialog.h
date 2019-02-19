#ifndef MAIN_DIALOG
#define MAIN_DIALOG

#include <QDialog>
#include <QImage>

#include "ui_main_dialog.h"

class QPaintEvent;
class QMouseEvent;

class MainDialog : public QDialog, public Ui::MainDialog {
  Q_OBJECT

  QImage img;
  int shx, shy;

 public:

  MainDialog(QWidget* = 0, Qt::WindowFlags = 0);

  public slots:
    void getText();
    void openFile();

 protected:
    void paintEvent(QPaintEvent*);
    void mouseMoveEvent(QMouseEvent*);
    void mousePressEvent(QMouseEvent*);
};

#endif
