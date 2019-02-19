#ifndef MAIN_DIALOG
#define MAIN_DIALOG

#include <QDialog>
#include "ui_main_dialog.h"

class SlowGlView;

class MainDialog : public QDialog, public Ui::mainDialog {
  Q_OBJECT

  SlowGlView *pglwidget;
 public:
  MainDialog(QWidget * = 0, Qt::WindowFlags = 0);
};
#endif
