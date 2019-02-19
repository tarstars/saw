#include "main_dialog.h"
#include "gl_view.h"

MainDialog::MainDialog(QWidget* parent, Qt::WindowFlags f):
  QDialog(parent, f){
  setupUi(this);

  pGlView = new GlView(this);
  glLayout -> addWidget(pGlView);

  connect(dialTheta, SIGNAL(valueChanged(int)), pGlView, SLOT(setTheta(int)));
  connect(dialPhi,   SIGNAL(valueChanged(int)), pGlView, SLOT(setPhi(int)));

  connect(dialRotx, SIGNAL(valueChanged(int)), pGlView, SLOT(setRotx(int)));
  connect(dialRoty, SIGNAL(valueChanged(int)), pGlView, SLOT(setRoty(int)));
  connect(dialRotz, SIGNAL(valueChanged(int)), pGlView, SLOT(setRotz(int)));

  dialTheta -> setValue(100);
  dialPhi -> setValue(100);

  dialRotx -> setValue(125);
  dialRoty -> setValue(180);
  dialRotz -> setValue(70);
}
