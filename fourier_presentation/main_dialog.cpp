#include "main_dialog.h"
#include "gl_view.h"
#include "gl_test_arrow.h"

MainDialog::MainDialog(QWidget* parent, Qt::WindowFlags f):QDialog(parent, f){
  setupUi(this);

  pSurfaceView = new GlView(this);
  pGlTestArrow = new GlTestArrow(this);
    
  glSurfaceLayout -> addWidget(pSurfaceView);
  glTestArrowLayout -> addWidget(pGlTestArrow);

  connect(rotx, SIGNAL(valueChanged(int)), pSurfaceView, SLOT(setRotx(int)));
  connect(roty, SIGNAL(valueChanged(int)), pSurfaceView, SLOT(setRoty(int)));
  connect(rotz, SIGNAL(valueChanged(int)), pSurfaceView, SLOT(setRotz(int)));
}
