#ifndef MAIN_DIALOG
#define MAIN_DIALOG

#include <QDialog>
#include "ui_main_dialog.h"

class MainDialog : public QDialog, public Ui::MainDialog{
  Q_OBJECT

    public:

  MainDialog(QWidget* = 0, Qt::WindowFlags = 0);
};

#endif