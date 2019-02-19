#include "main_dialog.h"
#include "gl_view.h"

MainDialog::MainDialog(QWidget* parent, Qt::WindowFlags f):
  QDialog(parent, f){

  setupUi(this);

  glLayout -> addWidget(new GlView(this));

}
